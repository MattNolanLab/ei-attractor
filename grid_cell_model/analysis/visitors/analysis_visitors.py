#
#   analysis_visitors.py
#
#   Visitors that perform data analysis on data.
#
#       Copyright (C) 2012  Lukas Solanka <l.solanka@sms.ed.ac.uk>
#       
#       This program is free software: you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation, either version 3 of the License, or
#       (at your option) any later version.
#       
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#       
#       You should have received a copy of the GNU General Public License
#       along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
import numpy as np
from interface        import DictDSVisitor, extractStateVariable, \
        extractSpikes, sumAllVariables
from otherpkg.log     import log_info
from analysis.signal  import localExtrema, butterBandPass, autoCorrelation
from analysis.image   import Position2D, fitGaussianBumpTT
from analysis.spikes  import slidingFiringRateTuple, ThetaSpikeAnalysis


__all__ = ['AutoCorrelationVisitor', 'BumpFittingVisitor', 'FiringRateVisitor']


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
    def __init__(self, monName, stateList, dtMult=1e-3, norm=True,
            bandStart=20, bandEnd=200, forceUpdate=False):
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
            sig, dt = sumAllVariables(mon, n_id, self.stateList)
            sig = butterBandPass(sig, dt*self.dtMult, self.bandStart,
                    self.bandEnd)
            ac = autoCorrelation(sig - np.mean(sig), max_lag=self.maxLag/dt,
                    norm=self.norm)
            ext_idx, ext_t = localExtrema(ac)
            acVec.append(ac)
    
            f, a = findFreq(ac, dt*1e-3, ext_idx, ext_t)
            freq.append(f)
            acval.append(a)
    
        return freq, acval, acVec


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
            freq, acVal, acVec = self.extractACStat(data[self.monName])
            a['freq']  = np.array(freq),
            a['acVal'] = np.array(acVal),
            a['acVec'] = np.array(acVec)
        else:
            log_info("AutoCorrelationVisitor", "Data present. Skipping analysis.")



class BumpFittingVisitor(DictDSVisitor):
    '''
    The bump fitting visitor takes spikes of the neural population (in the
    dictionary data set provided) and tries to
    fit a Gaussian shaped function to the population firing rate map (which is
    assumed to be a twisted torus). It saves the data under the 'analysis' key
    into the dataset.

    If the corresponding data is already present the visitor skips the data
    analysis and saving.
    '''
    
    def __init__(self, forceUpdate=False):
        self.forceUpdate = forceUpdate

        # Population FR parameters
        # All time units in msec
        self.dt     = 20.0
        self.winLen = 250.0

    def fitGaussianToMon(self, mon, Nx, Ny, tstart, tend):
        '''
        Fit a Gaussian function to the monitor mon, and return the results.
        '''
        N = Nx * Ny

        senders, times = extractSpikes(mon)
        F, Ft = slidingFiringRateTuple((senders, times), N, tstart, tend,
                self.dt, self.winLen)
            
        bumpT = tend - 2*self.winLen
        bumpI = bumpT / self.dt
        bump = np.reshape(F[:, bumpI], (Ny, Nx))
        dim = Position2D()
        dim.x = Nx
        dim.y = Ny
        return fitGaussianBumpTT(bump, dim)

    def visitDictDataSet(self, ds, **kw):
        '''
        Apply the bump fitting procedure onto the dataset 'ds' if necessary,
        and save the data into the dataset.
        '''
        data = ds.data
        if (('analysis' not in data.keys())):
            data['analysis'] = {}

        a = data['analysis']
        tstart = 0.0
        tend   = self.getOption(data, 'time')
        if (('bump_e' not in a.keys()) or self.forceUpdate):
            log_info("BumpFittingVisitor", "Analysing a dataset")
            # Fit the Gaussian onto E neurons
            Nx  = self.getNetParam(data, 'Ne_x')
            Ny  = self.getNetParam(data, 'Ne_y')
            mon = data['spikeMon_e']
            (A, mu_x, mu_y, sigma), err2 = self.fitGaussianToMon(mon, Nx, Ny,
                    tstart, tend)
            a['bump_e'] = {
                    'A' : A,
                    'mu_x' : mu_x,
                    'mu_y' : mu_y,
                    'sigma' : np.abs(sigma),
                    'err2'  : np.sum(err2)
            }

            # Fit the Gaussian onto I neurons
            # TODO
        else:
            log_info("BumpFittingVisitor", "Data present. Skipping analysis.")


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
        return ThetaSpikeAnalysis(N, senders, times)

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


