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
from scipy.optimize import leastsq
from os.path        import splitext

import analysis.signal as asignal
import analysis.spikes as aspikes
import data_storage.sim_models.ei as simei
from interface        import DictDSVisitor 
from otherpkg.log     import log_info, log_warn
from analysis.signal  import localExtrema, butterBandPass, autoCorrelation
from analysis.image   import Position2D, fitGaussianBumpTT


__all__ = ['AutoCorrelationVisitor', 'CrossCorrelationVisitor',
        'BumpFittingVisitor', 'FiringRateVisitor', 'BumpVelocityVisitor',
        'SpikeTrainXCVisitor', 'SpikeStatsVisitor']


def findFreq(ac, dt, ext_idx, ext_t):
    '''
    Find the first local maximum in an autocorrelation function and extract the
    frequency and the value of the autocorrelation at the detected frequency.

    Parameter, 'BumpVelocityVisitor']
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

        senders, times = simei.extractSpikes(mon)
        F, Ft = aspikes.slidingFiringRateTuple((senders, times), N, tstart, tend,
                self.dt, self.winLen)
            
        bumpT = tend - 2*self.winLen
        bumpI = bumpT / self.dt
        bump = np.reshape(F[:, bumpI], (Ny, Nx))
        dim = Position2D()
        dim.x = Nx
        dim.y = Ny
        return fitGaussianBumpTT(bump, dim), bump

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
            log_info("BumpFittingVisitor", "Analysing an E dataset")
            # Fit the Gaussian onto E neurons
            Nx  = self.getNetParam(data, 'Ne_x')
            Ny  = self.getNetParam(data, 'Ne_y')
            mon = data['spikeMon_e']
            ((A, mu_x, mu_y, sigma), err2), bump = self.fitGaussianToMon(mon,
                    Nx, Ny, tstart, tend)
            a['bump_e'] = {
                    'A' : A,
                    'mu_x' : mu_x,
                    'mu_y' : mu_y,
                    'sigma' : np.abs(sigma),
                    'err2'  : np.sum(err2),
                    'bump_e_rateMap' : bump
            }
        else:
            log_info("BumpFittingVisitor", "E data present. Skipping analysis.")

        # Only export population firing rates of I neurons
        if ('bump_i' not in a.keys() or self.forceUpdate):
            log_info("BumpFittingVisitor", "Analysing an I dataset")
            Nx  = simei.getNetParam(data, 'Ni_x')
            Ny  = simei.getNetParam(data, 'Ni_y')
            senders, times = simei.extractSpikes(data['spikeMon_i'])
            sheetSize = (Nx, Ny)
            spikes = aspikes.TorusPopulationSpikes(senders, times, sheetSize)

            F, Ft = spikes.slidingFiringRate(tstart, tend, self.dt,
                    self.winLen)
            bumpT = tend - 2*self.winLen
            bumpI = bumpT / self.dt
            bump = F[:, :, bumpI]

            a['bump_i'] = dict(
                    bump_i_rateMap=bump)
        else:
            log_info("BumpFittingVisitor", "I data present. Skipping analysis.")


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


###############################################################################


def fitCircularSlope(bumpPos, times, normFac):
    '''
    Fit a (circular) slope line onto velocity response of the bump and extract
    the slope (velocity), in neurons/s

    Parameters
    ----------
    bumpPos : numpy array
        An array of bump positions
    times : numpy array
        A corresponding vector of times, of the same size as bumpPos
    normFac : float
        Normalizing factor for the positional vector. In fact this is the size
        of the toroidal sheet (X direction).
    output : float
        Estimated bump speed (neurons/time unit).
    '''
    t = np.array(times) - times[0]
    bumpPos_norm = np.unwrap(1.0 * bumpPos / normFac * 2 * np.pi) # normalise to 2*Pi
    func = lambda X: X[0]*t - bumpPos_norm
    x0 = np.array([0.0])  # slope
    x = leastsq(func, x0)
    return x[0][0] / 2. / np.pi * normFac
    

def getLineFit(X, Y):
    '''
    Fit a line to data
    '''
    func = lambda P: P[0]*X  - Y
    P0 = np.array([0.0]) # slope
    P = leastsq(func, P0)
    return P[0][0]*X, P[0][0]


def fitBumpSpeed(IvelVec, bumpSpeed, fitRange):
    '''
    Fit a line to bumpSpeed vs IvelVec, useing IvelVec[0:fitRange+1]
    '''
    fitSpeed    = np.array(bumpSpeed)[:, 0:fitRange+1]
    fitIvelVec  = np.repeat([IvelVec[0:fitRange+1]], fitSpeed.shape[0],
            axis=0)
    fitSpeed    = fitSpeed.flatten()
    fitIvelVec   = fitIvelVec.flatten()
    line, slope = getLineFit(fitIvelVec, fitSpeed)
    lineFitErr  = np.abs(line - fitSpeed) / len(fitIvelVec)

    # Compose return values; as a function of IvelVec[0:fitRange+1]
    retLine = np.array(IvelVec[0:fitRange+1]) * slope
    retIvelVec = IvelVec[0:fitRange+1]
    return retLine, slope, lineFitErr, retIvelVec


class BumpVelocityVisitor(DictDSVisitor):
    '''
    A visitor that estimates the relationship between injected velocity current
    and bump speed.
    '''

    def __init__(self, bumpSpeedMax, win_dt=20.0, winLen=250.0,
            forceUpdate=False, printSlope=False):
        self.bumpSpeedMax = bumpSpeedMax # cm/s
        self.win_dt = win_dt
        self.winLen = winLen
        self.forceUpdate = forceUpdate
        self.printSlope = printSlope

        assert(self.bumpSpeedMax is not None)

    def _getSpikeTrain(self, data, monName, dimList):
        N_x = self.getNetParam(data, dimList[0])
        N_y = self.getNetParam(data, dimList[1])
        senders, times = simei.extractSpikes(data[monName])
        return senders, times, (N_x, N_y)


    def visitDictDataSet(self, ds, **kw):
        '''
        Visit a data set that contains all the trials and Ivel simulations.
        '''
        data = ds.data
        trials = data['trials']

        slopes = []
        for trialNum in xrange(len(trials)):
            log_info('BumpVelocityVisitor', "Trial no. {0}.".format(trialNum))
            slopes.append([])
            IvelVec = trials[trialNum]['IvelVec']
            for IvelIdx in xrange(len(IvelVec)):
                iData = trials[trialNum]['IvelData'][IvelIdx]
                if 'analysis' in iData.keys() and not self.forceUpdate:
                    log_info('BumpVelocityVisitor', "Data present. Skipping analysis.")
                    slopes[trialNum].append(iData['analysis']['slope'])
                    continue

                senders, times, sheetSize =  self._getSpikeTrain(iData,
                        'spikeMon_e', ['Ne_x', 'Ne_y'])
                pop = aspikes.TorusPopulationSpikes(senders, times, sheetSize)
                tStart = self.getOption(iData, 'theta_start_t')
                tEnd   = self.getOption(iData, 'time')
                bumpPos, bumpPos_t = pop.populationVector(tStart, tEnd,
                        self.win_dt, self.winLen)
                # NOTE: the bump moves in an opposite direction; we have to
                # negate the speed
                slope = -fitCircularSlope(bumpPos[:, 0], bumpPos_t,
                        sheetSize[0]/2.0)*1e3
                slopes[trialNum].append(slope)
                iData['analysis'] = {
                        'bumpPos'   : bumpPos,
                        'bumpPos_t' : bumpPos_t,
                        'slope'     : slope
                }
        slopes = np.array(slopes)

        analysisTop = {'bumpVelAll' : slopes}

        if (self.printSlope and 'fileName' not in kw.keys()):
            msg = 'printSlope requested, but did not receive the fileName ' + \
                    'as a keyword argument.'
            log_warn('BumpVelocityVisitor', msg)
            return
        elif (self.printSlope):
            if (len(trials) == 0):
                msg = 'Something wrong: len(trials) == 0'
                log_warn('BumpVelocityVisitor', msg)
                return


            # Plot the estimated bump velocities (nrns/s)
            from matplotlib.pyplot import figure, errorbar, xlabel, ylabel, \
                    plot, title, savefig, legend
            figure()
            IvelVec = trials[0]['IvelVec'] # All the same
            avgSlope = np.mean(slopes, axis=0)
            stdErrSlope = np.std(slopes, axis=0) / np.sqrt(len(trials))
            errorbar(IvelVec, avgSlope, stdErrSlope, fmt='o-')
            xlabel('Velocity current (pA)')
            ylabel('Bump velocity (neurons/s)')

            # Fit a line (nrns/s/pA)
            # results
            line       = None
            lineFitErr = None
            slope      = None
            fitIvelVec = None
            errSum     = None
            # All, in case we need the one with max range
            line_all       = []
            lineFitErr_all = []
            slope_all      = []
            fitIvelVec_all = []
            maxRange_all   = []
            for fitRange in xrange(1, len(IvelVec)):
                log_info('BumpVelocityVisitor', "fitRange: {0}".format(fitRange))
                newLine, newSlope, newErr, newIvelVecRange = \
                        fitBumpSpeed(IvelVec, slopes, fitRange)
                line_all.append(newLine)
                lineFitErr_all.append(newErr)
                slope_all.append(newSlope)
                fitIvelVec_all.append(newIvelVecRange)
                maxRange_all.append(newLine[-1])
                # Keep only lines that covers the desired range and have a
                # minimal error of fit
                if ((newLine[-1] >= self.bumpSpeedMax)):
                    errSumNew = np.sum(newErr)
                    if ((line is None) or (errSumNew <= errSum)):
                        line       = newLine
                        lineFitErr = newErr
                        slope      = newSlope
                        fitIvelVec = newIvelVecRange
                        errSum     = errSumNew

            if (line is None):
                msg = 'No suitable fits that cover <0, {0:.2f}> neurons/s.' +\
                ' Using the fit with the max. bump speed range.'
                log_info('BumpVelocityVisitor', msg.format(self.bumpSpeedMax))
                msg = "Bump speed maxima: {0}".format(maxRange_all)
                log_info('BumpVelocityVisitor', msg.format(self.bumpSpeedMax))
                maxRangeIdx = np.argmax(maxRange_all)
                line        = line_all[maxRangeIdx]
                lineFitErr  = lineFitErr_all[maxRangeIdx]
                slope       = slope_all[maxRangeIdx]
                fitIvelVec  = fitIvelVec_all[maxRangeIdx]

            plot(fitIvelVec, line, 'o-')
            plot(IvelVec, slopes.T, 'o', color='blue', alpha=0.4)
            t = "Line fit slope: {0:.3f} nrns/s/pA, error: {1:.3f} " + \
                    "neurons/s (norm)"
            title(t.format(slope, np.sum(lineFitErr)))
            legend(['Average bump speed', 'Line fit', 'Estimated bump speed'], loc='best')
            
            fileName = splitext(kw['fileName'])[0] + '.pdf'
            savefig(fileName)

            analysisTop.update({
                'lineFitLine'  : line,
                'lineFitSlope' : slope,
                'lineFitErr'   : lineFitErr,
                'fitIvelVec'   : fitIvelVec
            })

        data['analysis'] = analysisTop


        
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
