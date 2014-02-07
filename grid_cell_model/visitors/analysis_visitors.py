'''
Visitors that perform data analysis on data.
'''
import numpy as np

import analysis.signal as asignal
import analysis.spikes as aspikes
import data_storage.sim_models.ei as simei
from interface        import DictDSVisitor 
from otherpkg.log     import log_info, log_warn
from analysis.signal  import localExtrema, butterBandPass, autoCorrelation


import logging
logger = logging.getLogger(__name__)


__all__ = ['AutoCorrelationVisitor', 'CrossCorrelationVisitor',
        'FiringRateVisitor', 'SpikeTrainXCVisitor', 'SpikeStatsVisitor']


def findFreq(ac, dt, ext_idx, ext_t):
    '''
    Find the first local maximum in an autocorrelation function and extract the
    frequency and the value of the autocorrelation at the detected frequency.

    Parameters
    ----------
    ac : numpy vector
        A vector containing the autocorrelation function
    dt : float
        Sampling rate
    ext_idx : numpy array
        An array containig the indexes of local extrema (both maxima and
        minima) in 'ac'
    ext_t : numpy array
        The array of the same size as 'ext_idx', which contains the types of
        the local extrema ( >0 ... max, <0 ... min)
    output
        A tuple containig ('freq', 'corr'), where 'freq' is the frequency at
        the first local maximum and 'corr' is the value of the autocorrelation
        function at the same point.
    '''
    max_idx = np.nonzero(ext_t > 0)[0]
    if (len(max_idx) == 0):
        return (np.nan, np.nan)
    
    # First local maximum ac[0] excluded
    max1_idx = ext_idx[max_idx[0]]
    max1_t   = max1_idx * dt
    max1     = ac[max1_idx]

    return (1./max1_t, max1)




class AutoCorrelationVisitor(DictDSVisitor):
    '''
    A visitor to compute autocorrelations of state monitor data and extract
    information from them.

    The autocorrelation visitor takes as an input a state monitor, computes
    autocorrelations (with specified lag) of the specified synaptic currents
    and detects the major frequency, and power at that frequency, for all the
    monitored neurons. The results will be stored to the dictionary where the
    data came from.
    '''
    def __init__(self, monName, stateList, dtMult=1e-3, tStart=None, tEnd=None,
            norm=True, bandStart=20, bandEnd=200, forceUpdate=False):
        '''
        Initialise the visitor.

        Parameters
        ----------
        monName : string
            Name of the monitor; key in the data set dictionary
        stateList : list of strings
            A list of strings naming the state variables to extract (and sum)
        dtMult : float, optional
            dt Multiplier to transform dt into seconds
        tStart : float, optional
            Start time of the analysis. If None, the signal will not be
            cropped. The first value of the signal array is treated as time
            zero.
        tEnd : float, optional
            End time of the analysis. If None, the signal will not be cropped.
        norm : bool, optional
            Whether the autocorrelation function should be normalized
        bandStart : float, optional
            Bandpass start frequency
        bandEnd   : float, optional
            Bandpass end frequency
        forceUpdate : bool
            Whether to compute and store all the data even if they already
            exist in the data set.
        '''
        self.monName     = monName
        self.stateList   = stateList
        self.maxLag      = None
        self.dtMult      = dtMult
        self.tStart      = tStart
        self.tEnd        = tEnd
        self.norm        = norm
        self.bandStart   = bandStart
        self.bandEnd     = bandEnd
        self.forceUpdate = forceUpdate



    def extractACStat(self, mon):
        '''
        Extract autocorrelation statistics from a monitor.

        For each monitored neuron, extract the (highest) frequency, value of
        the autocorrelation at the frequency and the autocorrelation function
        itself.
    
        Parameters
        ----------
        mon : list of dicts
            A list of (NEST) state monitors' status dictionaries
        output : tuple
            A tuple (freq, acval, acVec), containing the arrays of frequencies
            for the monitored neurons, autocorrelation values at the
            corresponding frequencies, and autocorrelation functions of all the
            neurons.
        '''
        freq   = [] # Frequency of input signal
        acval  = [] # Auto-correlation at the corresponding frequency
        acVec  = []
        for n_id in range(len(mon)):
        #for n_id in range(5):
            #print "n_id: ", n_id
            sig, dt = simei.sumAllVariables(mon, n_id, self.stateList)
            startIdx = 0
            endIdx   = len(sig)
            if (self.tStart is not None):
                startIdx = int(self.tStart / dt)
            if (self.tEnd is not None):
                endIdx = int(self.tEnd / dt)
            sig = sig[startIdx:endIdx]
            sig = butterBandPass(sig, dt*self.dtMult, self.bandStart,
                    self.bandEnd)
            ac = autoCorrelation(sig - np.mean(sig), max_lag=self.maxLag/dt,
                    norm=self.norm)
            ext_idx, ext_t = localExtrema(ac)
            acVec.append(ac)
    
            f, a = findFreq(ac, dt*self.dtMult, ext_idx, ext_t)
            freq.append(f)
            acval.append(a)
    
        return freq, acval, acVec, dt


    def visitDictDataSet(self, ds, **kw):
        '''
        Visit the dictionary data set and extract frequency, autocorrelation
        for the detected frequency, and autocorrelation functions, for all the
        monitored neurons.  The parameters are defined by the constructor of
        the object.

        If the analysed data is already present, the analysis and storage of
        the data will be skipped.

        Parameters
        ----------
        ds : a dict-like object
            A data set to perform analysis on.
        '''
        data = ds.data
        if (not self.folderExists(data, ['analysis'])):
            data['analysis'] = {}
        a = data['analysis']

        if (('freq' not in a.keys()) or self.forceUpdate):
            log_info("AutoCorrelationVisitor", "Analysing a dataset")
            o = data['options']
            self.maxLag = 1. / (o['theta_freq'] * 1e-3)
            freq, acVal, acVec, dt = self.extractACStat(data[self.monName])
            print freq
            a['freq']  = np.array(freq)
            a['acVal'] = np.array(acVal)
            a['acVec'] = np.array(acVec)
            a['ac_dt'] = dt
        else:
            log_info("AutoCorrelationVisitor", "Data present. Skipping analysis.")



class CrossCorrelationVisitor(DictDSVisitor):
    '''
    A visitor that computes a cross-correlation function between a set of state
    monitors.

    The crosscorrelation visitor takes a list of state monitors, each monitor
    collecting data from one neuron (Vm, Isyn, etc.) and computes the cross
    correlation function between all the pairs of these signals.

    The results will be stored to the input data dictionary.
    '''
    def __init__(self, monName, stateList, maxLag=None, tStart=None, tEnd=None,
            norm=True, forceUpdate=False):
        '''
        Initialise the visitor.

        Parameters
        ----------
        monName : string
            Name of the monitor; key in the data set dictionary
        stateList : list of strings
            A list of strings naming the state variables to extract (and sum)
        maxLag : int
            Maximal lag (positive and negative) of the cross-correlation
            function.
        tStart : float, optional
            Start time of the analysis. If None, the signal will not be
            cropped. The first value of the signal array is treated as time
            zero.
        tEnd : float, optional
            End time of the analysis. If None, the signal will not be cropped.
        norm : bool, optional
            Whether the cross-correlation function should be normalized
        forceUpdate : bool
            Whether to compute and store all the data even if they already
            exist in the data set.
        '''
        self.monName     = monName
        self.stateList   = stateList
        self.maxLag      = maxLag
        self.tStart      = tStart
        self.tEnd        = tEnd
        self.norm        = norm
        self.forceUpdate = forceUpdate


    def extractCCStat(self, mon, out):
        '''
        Extract x-correlation statistics from a monitor.

        For each pair of monitored neurons    
        Parameters
        ----------
        mon : list of dicts
            A list of (NEST) state monitors' status dictionaries
        out : dictionary
            Output data dictionary.
        '''
        out['x-corr'] = dict(
                correlations=[])
        xcOut = out['x-corr']['correlations']
        for n_id1 in range(len(mon)):
            sig1, dt1 = simei.sumAllVariables(mon, n_id1, self.stateList)
            xcOut.append([])
            xcOut2 = xcOut[n_id1]
            for n_id2 in range(len(mon)):
                print('n_id1, n_id2 = {0}, {1}'.format(n_id1, n_id2))
                sig2, dt2 = simei.sumAllVariables(mon, n_id2, self.stateList)
                if (dt1 != dt2):
                    raise ValueError('dt1 != dt2')

                dt        = dt1
                startIdx  = 0
                lag_start = -int(self.maxLag/dt)
                lag_end   = -lag_start
                endIdx1   = len(sig1)
                endIdx2   = len(sig2)

                if (self.tStart is not None):
                    startIdx = int(self.tStart / dt)
                if (self.tEnd is not None):
                    endIdx1 = int(self.tEnd / dt)
                    endIdx2 = endIdx1
                sig1 = sig1[startIdx:endIdx1]
                sig2 = sig2[startIdx:endIdx2]
                C = asignal.corr(sig1, sig2, mode='range',
                        lag_start=lag_start, lag_end=lag_end)
                if (self.norm):
                    C /= np.max(C)
                xcOut2.append(C)
        out['x-corr']['lags'] = np.arange(lag_start, lag_end+1) * dt
        
        


    def visitDictDataSet(self, ds, **kw):
        '''
        Visit the dictionary data set and extract the cross-correlation
        functions, for all pairs of the monitored neurons monitored neurons.
        The parameters are defined by the constructor of the object.

        If the analysed data is already present, the analysis and storage of
        the data will be skipped.

        Parameters
        ----------
        ds : a dict-like object
            A data set to perform analysis on.
        '''
        data = ds.data
        if (not self.folderExists(data, ['analysis'])):
            data['analysis'] = {}
        a = data['analysis']

        if (('x-corr' not in a.keys()) or self.forceUpdate):
            log_info("CrossCorrelationVisitor", "Analysing a dataset")
            o = data['options']
            if (self.maxLag is None):
                self.maxLag = 1. / (o['theta_freq'] * 1e-3)
            self.extractCCStat(data[self.monName], a)
        else:
            log_info("CrossCorrelationVisitor", "Data present. Skipping analysis.")



class FiringRateVisitor(DictDSVisitor):
    '''
    Determine various firing rate statistics of a population of neurons on the
    spiking data dataset:
        * Average firing rate in the middle of the theta cycle

    Save the data to the original data set.
    '''

    def __init__(self, winLen=0.5, thetaStartT=None, thetaFreq=None, tEnd=None, 
            forceUpdate=False):
        '''
        Initialize the visitor.

        Parameters
        ----------
        winLen : float (ms)
            Length of the firing rate window as a fraction of the theta cycle
            time (0, 1>.
        thetaStartT : float (ms)
            Start time of the theta signal. The center of the firing rate
            window will be in the middle of the theta signal. Therefore it is
            up to the user to ensure that the peak of the theta signal is in
            the middle. If None, extract from the data when performing analysis
        thetaFreq : float (Hz)
            Theta signal frequency. If None, extract from the data
        tEnd : float (ms)
            Analysis end time. If None, extract from the data
        forceUpdate : boolean, optional
            Whether to do the data analysis even if the data already exists.
        '''
        self.thetaStartT = thetaStartT
        self.thetaFreq   = thetaFreq
        self.tEnd        = tEnd
        self.winLen      = winLen
        self.forceUpdate = forceUpdate


    def _getSpikeTrain(self, data, monName, dimList):
        senders, times, N = DictDSVisitor._getSpikeTrain(self, data, monName,
                dimList)
        return aspikes.ThetaSpikeAnalysis(N, senders, times)

    def visitDictDataSet(self, ds, **kw):
        data = ds.data
        if (not self.folderExists(data, ['analysis'])):
            data['analysis'] = {}
        a = data['analysis']

        thetaStartT = self._checkAttrIsNone(self.thetaStartT,
                'theta_start_t', data)
        thetaFreq = self._checkAttrIsNone(self.thetaFreq, 'theta_freq',
                data)
        tEnd = self._checkAttrIsNone(self.tEnd, 'time', data)
        if (not self.folderExists(a, ['FR_e']) or self.forceUpdate):
            log_info('FiringRateVisitor', "Analysing...")
            eSp = self._getSpikeTrain(data, 'spikeMon_e', ['Ne_x', 'Ne_y'])
            eFR = eSp.avgFiringRate(thetaStartT, tEnd)
            a['FR_e'] = {
                    'all' : eFR,
                    'avg' : np.mean(eFR)
            }

            iSp = self._getSpikeTrain(data, 'spikeMon_i', ['Ni_x', 'Ni_y'])
            iFR = iSp.avgFiringRate(thetaStartT, tEnd)
            a['FR_i'] = {
                    'all' : iFR,
                    'avg' : np.mean(iFR)
            }
        else:
            log_info("FiringRateVisitor", "Data present. Skipping analysis.")


##############################################################################

class SpikeTrainXCVisitor(DictDSVisitor):
    '''
    Compute spike train crosscorrelations between spikes of neuron of a
    population.
    '''

    def __init__(self, monitorName, bins, lagRange=None, neuronIdx=None,
            forceUpdate=False):
        '''
        Parameters:

        monitorName : string
            Name of the monitor in the data hierarchy.
        bins : int
            Number of bins for the cross-correlation histogram. Bin centers
        lagRange : (tStart, tEnd)
            Start/end time of the lags in histograms. This values will define
            the left and right edges of the histogram. All the values outside
            this range will be ignored.
        neuronIdx : list of ints or None
            List of neurons for which to compute the pairwise
            cross-correlation. If None, use all the neurons.
        '''
        self.allowedMonitors = ['spikeMon_e', 'spikeMon_i']
        if (not monitorName in self.allowedMonitors):
            msg = "monitorName must be one of {0}".format(allowedMonitors)
            raise ValueError(msg)

        self.monitorName = monitorName
        self.lagRange    = lagRange
        self.bins        = bins
        self.forceUpdate = forceUpdate
        self.neuronIdx   = neuronIdx

        if (self.monitorName == "spikeMon_e"):
            self.NName = "net_Ne"
            self.outputName = "XCorrelation_e"
        elif (self.monitorName == "spikeMon_i"):
            self.NName = "net_Ni"
            self.outputName = "XCorrelation_i"


    def visitDictDataSet(self, ds, **kw):
        data = ds.data
        
        if (not self.folderExists(data, ['analysis'])):
            data['analysis'] = {}
        a = data['analysis']

        if (self.outputName in a.keys() and not self.forceUpdate):
            log_info("SpikeTrainXCorrelation", "Data present. Skipping analysis.")
            return

        spikes = simei.MonitoredSpikes(data, self.monitorName, self.NName)
        if (self.neuronIdx is None):
            idx1 = range(spikes.N)
        else:
            idx1 = self.neuronIdx
        correlations, bin_centers, bin_edges = spikes.spikeTrainXCorrelation(idx1, None,
                self.lagRange, self.bins)

        a[self.outputName] = dict(
                neuronIdx    = idx1,
                correlations = correlations,
                bin_edges    = bin_edges,
                bin_centers  = bin_centers)


class SpikeStatsVisitor(DictDSVisitor):
    def __init__(self, monitorName, forceUpdate=False):
        '''
        Parameters:

        monitorName : string
            Name of the monitor in the data hierarchy.
        '''
        self.allowedMonitors = ['spikeMon_e', 'spikeMon_i']
        if (not monitorName in self.allowedMonitors):
            msg = "monitorName must be one of {0}".format(allowedMonitors)
            raise ValueError(msg)

        self.monitorName = monitorName
        self.forceUpdate = forceUpdate

        if (self.monitorName == "spikeMon_e"):
            self.NName = "net_Ne"
            self.outputName = "CV_e"
        elif (self.monitorName == "spikeMon_i"):
            self.NName = "net_Ni"
            self.outputName = "CV_i"


    def visitDictDataSet(self, ds, **kw):
        data = ds.data
        
        if (not self.folderExists(data, ['analysis'])):
            data['analysis'] = {}
        a = data['analysis']

        if (self.outputName in a.keys() and not self.forceUpdate):
            log_info("SpikeStatsVisitor", "Data present. Skipping analysis.")
            return

        spikes = simei.MonitoredSpikes(data, self.monitorName, self.NName)
        a[self.outputName] = np.array(spikes.ISICV())
