#
#   spikes.py
#
#   Functions for analysis of spike trains
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
'''
Classes and functions for spike train analysis.

.. currentmodule:: analysis.spikes

Classes
-------
.. autosummary::

    PopulationSpikes
    TorusPopulationSpikes
    ThetaSpikeAnalysis


Functions
---------
These functions are deprecated and should not be used unless really needed.
Please use the classes defined above.

.. autosummary::
    slidingFiringRateTuple
    torusPopulationVector


'''
import numpy as np
import scipy
import collections
from scipy import weave

import _spikes
from otherpkg.log import log_warn

__all__ = [ 'slidingFiringRateTuple', 'torusPopulationVector',
        'torusPopulationVectorFromRates', 'SpikeTrain', 'PopulationSpikes',
        'ThetaSpikeAnalysis', 'TorusPopulationSpikes']


#def firingRate(spikeTrain, tstart, tend):
#    '''
#    Compute an average firing rate of a spike train, between tstart and tend
#    spikeTrain      A list or a 1D array of spike times
#    tstart, tend    Must be the same units as spikeTrain
#    '''
#    return np.sum(np.logical_and(spikeTrain >= tstart, spikeTrain <= tend)) / \
#            (tend - tstart)
#
#
#def multipleFiringRate(spikeTrainArray, tstart, tend):
#    '''
#    Compute an average firing rate of an array of spike trains. For each spike
#    train, return the firing rate between tstart and tend. Thus, the result will
#    be of len(spikeTrainArray)
#    '''
#    aLen = len(spikeTrainArray)
#    result = np.ndarray((aLen, ))
#    for it in xrange(aLen):
#        result[it] = firingRate(spikeTrainArray[it].flatten(), tstart, tend)
#
#    return result
#
#
#
#def firingRateSlidingWindow(spikeTrain, tstart, tend, dt, winLen):
#    '''
#    Compute sliding window firing rate for the neuron
#    dt      resolution of firing rate
#    winLen  Length of the sliding window
#
#    The spikes are computed from tstart to tend, so that the resulting array
#    length is int((tend-tstart)/dt)+1 long.
#    dt does not have to be relevant to simulation dt at all
#    '''
#    #lg.debug('Start firing rate processing')
#    szRate = int((tend-tstart)/dt)+1
#    r = np.ndarray((len(spikeTrain), szRate))
#    times = np.ndarray(szRate)
#    for n_i in xrange(len(spikeTrain)):
#        tmp = np.array(spikeTrain[n_i])
#        for t_i in xrange(szRate):
#            t = tstart + t_i*dt
#            r[n_i][t_i] = np.sum(np.logical_and(tmp > t-winLen/2, tmp <
#                t+winLen/2))
#            times[t_i] = t
#
#    #lg.debug('End firing rate processing')
#    return (r/winLen, times)


def slidingFiringRateTuple(spikes, N, tstart, tend, dt, winLen):
    '''
    Compute a firing rate with a sliding window from a tuple of spike data:
    spikes is a tuple(n_id, times), in which n_id is a list/array of neuron id
    and times is a list/array of spike times

    Parameters
    ----------
    spikes  A pair (n_id, spikes)
    N       Total number of neurons
    tstart  When the firing rate will start (ms)
    tend    End time of firing rate (ms)
    dt      Sliding window dt - not related to simulation time (ms)
    winLen  Length of the sliding window (ms). Must be >= dt.

    return  An array of shape (N, int((tend-tstart)/dt)+1
    '''
    #print "Start sliding firing rate.."
    
    szRate      = int((tend-tstart)/dt)+1
    n_ids       = np.array(spikes[0])
    spikeTimes  = np.array(spikes[1])
    lenSpikes   = len(spikeTimes)
    bitSpikes   = np.zeros((N, szRate))
    fr          = np.zeros((N, szRate))
    dtWlen      = int(winLen/dt)
    times       = np.arange(tstart, tend+dt, dt)
    N           = int(N)

    #print "max(n_ids): ", np.max(n_ids)
    #print 'szRate: ', szRate
    #print 'N: ', N
    #print 'dtWlen: ', dtWlen

    code = """
        for (int i = 0; i < lenSpikes; i++)
        {
            int spikeSteps = (spikeTimes(i) - tstart) / dt;
            if (spikeSteps >= 0 && spikeSteps < szRate)
            {
                int n_id = n_ids(i);
                bitSpikes(n_id, spikeSteps) += 1;
            }
        }

        for (int n_id = 0; n_id < N; n_id++)
            for (int t = 0; t < szRate; t++)
            {
                fr(n_id, t) = .0;
                for (int s = 0; s < dtWlen; s++)
                    if ((t+s) < szRate)
                        fr(n_id, t) += bitSpikes(n_id, t+s);
            }
        """

    err = weave.inline(code,
            ['N', 'szRate', 'dtWlen', 'lenSpikes', 'n_ids', 'spikeTimes',
                'tstart', 'dt', 'bitSpikes', 'fr'],
            type_converters=weave.converters.blitz,
            compiler='gcc',
            extra_compile_args=['-O3'],
            verbose=2)

    #print "End sliding firing rate"

    return fr/(winLen*1e-3), times

        



def torusPopulationVector(spikes, sheetSize, tstart=0, tend=-1, dt=0.02, winLen=1.0):
    '''
    This function is deprecated. Use the OO version instead
    '''
    log_warn('This function is deprecated')
    N = sheetSize[0]*sheetSize[1]
    F, tsteps = slidingFiringRateTuple(spikes, N, tstart, tend, dt, winLen)

    return torusPopulationVectorFromRates((F, tsteps), sheetSize)



class SpikeTrain(object):
    '''
    A base class for handling spike trains.
    '''
    def __init__(self):
        raise NotImplementedError()



class PopulationSpikes(SpikeTrain, collections.Sequence):
    '''
    Class to handle a population of spikes and a set of methods to do analysis
    on the *whole* population.
    '''
    def __init__(self, N, senders, times):
        '''
        N : int
            Number of neurons in the population
        senders : 1D array
            Neuron numbers corresponding to the spikes
        times : 1D array
            Spike times. The shape of this array must be the same as for
            `senders`.
        '''
        self._N        = N
        if (N < 0):
            msg = "Number of neurons in the spike train must be " +\
                    "non-negative! Got {0}."
            raise ValueError(msg.format(N))

        # We are expecting senders and times as numpy arrays, if they are not,
        # convert them. Moreover, senders.dtype must be int, for indexing.
        self._senders  = np.asarray(senders, dtype=int)
        self._times    = np.asarray(times)
        self._unpacked = [None] * self._N # unpacked version of spikes


    @property
    def N(self):
        '''
        Number of neurons in the population
        '''
        return self._N


    def avgFiringRate(self, tStart, tEnd):
        '''
        Compute and average firing rate for all the neurons between 'tstart'
        and 'tend'. Return an array of firing rates, one item for each neuron
        in the population.

        Parameters
        ----------
        tStart : float (ms)
            Start time.
        tEnd   : float (ms)
            End time.
        output : numpy array
            Firing rate in Hz for each neuron in the population.
        '''
        result  = np.zeros((self._N, ))
        times   = self._times
        senders = self._senders
        N       = int(self._N)
        ts      = float(tStart)
        te      = float(tEnd)
        code = '''
            for (int i = 0; i < senders.size(); i++)
            {
                int t = times(i);
                int s = senders(i);
                if (s >= 0 && s < N && t >= ts && t <= te)
                    result(s)++; 
                else if (s < 0 || s >= N)
                    std::cout << "senders is outside range <0, N)" << 
                            std::endl;
            }
        '''

        err = weave.inline(code,
                ['N', 'times', 'senders', 'ts', 'te', 'result'],
                type_converters=weave.converters.blitz,
                compiler='gcc',
                extra_compile_args=['-O3'],
                verbose=2)
        return 1e3 * result / (tEnd - tStart)


    def slidingFiringRate(self, tStart, tEnd, dt, winLen):
        '''
        Compute a sliding firing rate over the population of spikes, by taking
        a rectangular window of specified length.

        Parameters
        ----------
        tStart : float
            Start time of the firing rate analysis.
        tEnd : float
            End time of the analysis
        dt : float
            Firing rate window time step
        winLen : float
            Lengths of the windowing function (rectangle)
        output : a tuple
            A pair (F, t), specifying the vector of firing rates and
            corresponding times. F is a 2D array of the shape (N, Ntimes), in
            which N is the number of neurons and Ntimes is the number of time
            steps. 't' is a vector of times corresponding to the time windows
            taken.
        '''
        spikes = (self._senders, self._times)
        return slidingFiringRateTuple(spikes, self._N, tStart, tEnd, dt,
                winLen)


    def windowed(self, tLimits):
        '''
        Return population spikes restricted to tLimits.

        Parameters
        ----------
        tLimits : a pair
            A tuple (tStart, tEnd). The spikes in the population must satisfy
            tStart >= t <= tEnd.
        output : PopulationSpikes instance
            A copy of self with only a subset of spikes, limited by the time
            window.
        '''
        tStart = tLimits[0]
        tEnd   = tLimits[1]
        tIdx = np.logical_and(self._times >= tStart, self._times <= tEnd)
        return PopulationSpikes(self._N, self._senders[tIdx],
                self._times[tIdx])


    def rasterData(self, neuronList=None):
        '''
        Extract the senders and corresponding spike times for a raster plot.
        TODO: implement neuronList

        Parameters
        ==========
        neuronList : list, optional
            Extract only neurons given in this list
        output : a tuple
            A pair containing (senders, times).
        '''
        if (neuronList is not None):
            raise NotImplementedError()

        return self._senders, self._times




    def spikeTrainDifference(self, idx1, idx2=None, full=True, reduceFun=None):
        '''
        Compute time differences between pairs of spikes of two neurons or a
        list of neurons.
        
        Parameters
        ----------
        idx1 : int, or a sequence of ints
            Index of the first neuron or a list of neurons for which to compute
            the correlation histogram.
        idx2 : int, or a sequence of ints, or None
            Index of the second neuron or a list of indexes for the second set
            of spike trains.
        full : bool, optional
            Not fully implemented yet. Must be set to True.
        reduceFun : callable, optional
            Any callable object that computes a function over an array of each
            spike train difference. The function must take one input argument,
            which will be the array of spike time differences for a pair of
            neurons. The output of this function will be stored instead of the
            default output.
        output : A 2D or 1D array of spike train autocorrelation histograms for all
            the pairs of neurons.
        
        The computation takes the following steps:
        
         * If ``idx1`` or ``idx2`` are integers, they will be converted to a list
           of size 1.
         * If ``idx2`` is None, then the result will be a list of lists of pairs
           of cross-correlations between the neurons. Even if there is only one
           neuron. If ``full == True``, the output will be an upper triangular
           matrix of all the pairs, i.e. it will exclude the duplicated.
           Otherwise there will be cross correlation histograms between all the
           pairs.
         * if ``idx2`` is not None, then ``idx1`` and ``idx2`` must be arrays
           of the same length, specifying the pairs to compute autocorrelation
           for
        '''
        if (full == False):
            raise NotImplementedError()

        if (reduceFun is None):
            reduceFun = lambda x: x

        if (not isinstance(idx1, collections.Iterable)):
            idx1 = [idx1]

        if (idx2 is None):
            idx2 = idx1
            res = [[] for x in idx1]
            for n1 in idx1:
                for n2 in idx2:
                    print n1, n2, len(self[n1]), len(self[n2])
                    res[n1].append(reduceFun(_spikes.spike_time_diff(self[n1],
                        self[n2])))
            return res
        elif (not isinstance(idx2, collections.Iterable)):
            idx2 = [idx2]

        # Two arrays of pairs
        if (len(idx1) != len(idx2)):
            raise TypeError('Length of neuron indexes do not match!')

        res = [None] * len(idx1)
        for n in xrange(len(idx1)):
            res[n] = reduceFun(_spikes.spike_time_diff(self[idx1[n]],
                self[idx2[n]]))

        return res


    class CallableHistogram(object):
        def __init__(self, **kw):
            self.kw = kw

        def __call__(self, x):
            '''
            Perform the histogram on x and return the result of
            numpy.histogram, without bin_edges
            '''
            res, _ = np.histogram(x, **self.kw)
            return res

        def get_bin_edges(self):
            _, bin_edges = np.histogram([], **self.kw)
            return bin_edges


    def spikeTrainXCorrelation(self, idx1, idx2, range, bins=50,
            **kw):
        '''
        Compute the spike train crosscorrelation function for all pairs of
        spike trains in the population.

        Parameters (for explanation of how ``idx1`` and ``idx2`` are treated,
        see :meth:`~PopulationSpikes.spikeTrainDifference`):

        idx1 : int, or a sequence of ints
            Index of the first neuron or a list of neurons for which to compute
            the correlation histogram.
        idx2 : int, or a sequence of ints, or None
            Index of the second neuron or a list of indexes for the second set
            of spike trains.
        range : (lag_start, lag_end)
            Limits of the cross-correlation function. The bins will always be
            **centered** on the values.
        bins : int, optional
            Number of bins
        kw : dict
            Keyword arguments passed on to the numpy.histogram function

        output : a 2D or 1D list
            see :meth:`~PopulationSpikes.spikeTrainDifference`.
        '''
        lag_start = range[0]
        lag_end   = range[1]
        binWidth = (lag_end - lag_start) / (bins - 1)
        bin_edges = np.linspace(lag_start - binWidth/2.0, lag_end +
                binWidth/2.0, bins+1) 
        h = self.CallableHistogram(bins=bin_edges, **kw)
        XC = self.spikeTrainDifference(idx1, idx2, full=True, reduceFun=h)
        bin_edges = h.get_bin_edges()
        bin_centers = (bin_edges[0:-1] + bin_edges[1:])/2.0
        return XC, bin_centers, bin_edges
        
        
        
        

    #######################################################################
    # Functions implementing collections.Sequence
    def __getitem__(self, key):
        '''Retrieve a spike train for one neuron.'''
        if self._unpacked[key] is not None:
            return self._unpacked[key]
        ret = self._times[self._senders == key]
        self._unpacked[key] = ret
        return ret

    def __len__(self):
        return self._N






class TorusPopulationSpikes(PopulationSpikes):
    '''
    Spikes of a population of neurons on a twisted torus.
    '''

    def __init__(self, senders, times, sheetSize):
        self._sheetSize = sheetSize
        N = sheetSize[0]*sheetSize[1]
        PopulationSpikes.__init__(self, N, senders, times)


    def getXSize(self):
        return self._sheetSize[0]
    def getYSize(self):
        return self._sheetSize[1]

    
    def populationVector(self, tStart, tEnd, dt, winLen):
        '''
        Compute the population vector on a torus, from the spikes present. Note
        that this method will have a limited functionality on a twisted torus,
        but can be used if the population activity translates in the X
        dimension only.

        Parameters
        ----------
        tStart : float
            Start time of analysis
        tEnd : float
            End time of analysis
        dt : float
            Time step of the (rectangular) windowing function
        winLen : float
            Length of the windowing function
        output : tuple
            A pair (r, t) in which r is a 2D vector of shape
            (int((tEnd-tStart)/dt)+1), 2), corresponding to the population
            vector for each time step of the windowing function, and t is a
            vector of times, of length the first dimension of r.
        '''
        sheetSize_x = self.getXSize()
        sheetSize_y = self.getYSize()
        N = sheetSize_x*sheetSize_y
        
        F, tsteps = PopulationSpikes.slidingFiringRate(self, tStart, tEnd, dt,
                winLen)
        P = np.ndarray((len(tsteps), 2), dtype=complex)
        X, Y = np.meshgrid(np.arange(sheetSize_x), np.arange(sheetSize_y))
        X = np.exp(1j*(X - sheetSize_x/2)/sheetSize_x*2*np.pi).ravel()
        Y = np.exp(1j*(Y - sheetSize_y/2)/sheetSize_y*2*np.pi).ravel()
        for t_it in xrange(len(tsteps)):
            P[t_it, 0] = np.dot(F[:, t_it], X)
            P[t_it, 1] = np.dot(F[:, t_it], Y)

        return (np.angle(P)/2/np.pi*self._sheetSize, tsteps)


    def slidingFiringRate(self, tStart, tEnd, dt, winLen):
        '''
        Compute a sliding firing rate over the population of spikes, by taking
        a rectangular window of specified length. However, unlike the ancestor
        method (PopulationSpikes.slidingFiringRate), return a 3D array, a
        succession of 2D population firing rates in time.

        Parameters
        ----------
        tStart : float
            Start time of the firing rate analysis.
        tEnd : float
            End time of the analysis
        dt : float
            Firing rate window time step
        winLen : float
            Lengths of the windowing function (rectangle)
        output : a tuple
            A pair (F, t), specifying the vector of firing rates and
            corresponding times. F is a 3D array of the shape (Nx, Ny, Ntimes),
            in which Nx/Ny are the number of neurons in X and Y dimensions,
            respectively, and Ntimes is the number of time steps. 't' is a
            vector of times corresponding to the time windows taken.
        '''
        spikes = (self._senders, self._times)
        F, Ft = slidingFiringRateTuple(spikes, self._N, tStart, tEnd, dt,
                winLen)
        Nx = self.getXSize()
        Ny = self.getYSize()
        return np.reshape(F, (Ny, Nx, len(Ft))), Ft





class ThetaSpikeAnalysis(PopulationSpikes):
    '''
    Analyse population spike trains for theta oscillation-related information.
    '''
    def __init__(self, N, senders, times):
        PopulationSpikes.__init__(self, N, senders, times)


    def firingRateMiddleTheta(self, thetaStartT, thetaFreq, tEnd, winLen):
        '''
        Compute firing rate for every neuron in the population. For each
        neuron, return an array of firing rates for every theta cycle.

        Parameters
        ----------
        thetaStartT : float (ms)
            Start time of the theta signal. The center of the firing rate
            window will be in the middle of the theta signal. Therefore it is
            up to the user to ensure that the peak of the theta signal is in
            the middle.
        thetaFreq : float (Hz)
            Theta signal frequency
        tEnd : float (ms)
            Analysis end time
        winLen : float
            Length of the firing rate window as a fraction of the theta cycle
            time. Will be derived from thetaFreq, and should be in the range
            (0, 1>.
        '''
        thetaT = 1e3 / thetaFreq # ms
        spikes = (self._senders, self._times)
        return slidingFiringRateTuple(spikes, self._N, thetaStartT + .5*thetaT,
                tEnd, thetaT, winLen * thetaT)


    def avgFiringRateMiddleTheta(self, thetaStartT, thetaFreq, tEnd, winLen):
        '''
        Do the same thing as getFiringRateMiddleTheta, but for each neuron also
        compute its average firing rate from theta_start_t to tend
        '''
        fr, times = self.firingRateMiddleTheta(theta_start_t, theta_freq, tend, winlen)
        return np.mean(fr, 1)



