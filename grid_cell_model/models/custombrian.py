from brian.monitor import SpikeMonitor, StateMonitor
from scipy.io import savemat

import numpy as np
import operator as op
import logging as lg

from datetime import datetime

from analysis.spikes import firingRateSlidingWindow, torusPopulationVectorFromRates



class ExtendedSpikeMonitor(SpikeMonitor):
    '''
    SpikeMonitor class that stores spikes implicitly as an array of lists,
    instead of a list of pairs (n, t) as does the original SpikeMonitor class
    '''
    def __init__(self, source, record=True, delay=0, function=None):
        SpikeMonitor.__init__(self, source, record, delay, function)

        self.aspikes = np.ndarray(len(source), dtype=list)
        self.reinit()

    def reinit(self):
        self._newspikes = True
        self.nspikes = 0
        for i in xrange(len(self.source)):
            self.aspikes[i] = []

    def propagate(self, spikes):
        if len(spikes):
            self._newspikes = True
        self.nspikes += len(spikes)
        if self.record:
            tmp = np.ndarray(len(spikes), dtype=list)
            tmp[:] = [[self.source.clock.t]]*len(spikes)
            self.aspikes[spikes] += tmp


    def __getitem__(self, i):
        return np.array(self.aspikes[i])

    def getFiringRate(self, tstart, tend, dt, winLen):
        '''
        Compute sliding window firing rate for the neuron
        dt      resolution of firing rate
        winLen  Length of the sliding window

        The spikes are computed from tstart to tend, so that the resulting array
        length is int((tend-tstart)/dt)+1 long.
        dt does not have to be relevant to simulation dt at all
        '''
        return firingRateSlidingWindow(self.aspikes, tstart, tend, dt, winLen)



    def getFiringRateMiddleTheta(self, theta_start_t, theta_freq, tend, winlen):
        '''
        Compute firing rate for every neuron with window of length winlen, in
        the middle of each theta cycle.
        This will produce an array of average firing rates for each theta cycle,
        until tend.
        '''
        theta_T = 1. / theta_freq
        return self.getFiringRate(theta_start_t + .5*theta_T, tend, theta_T, winlen)


    def getAvgFiringRateMiddleTheta(self, theta_start_t, theta_freq, tend, winlen):
        '''
        Do the same thing as getFiringRateMiddleTheta, but for each neuron also
        compute its average firing rate from theta_start_t to tend
        '''
        fr, times = self.getFiringRateMiddleTheta(theta_start_t, theta_freq, tend, winlen)
        return np.mean(fr, 1)


    def getNSpikes(self):
        total = np.ndarray(len(self.aspikes))
        for n_it in xrange(len(self.aspikes)):
            total[n_it] = len(self.aspikes[n_it])
        return total

    def getTimedFiringRate(self, startT, endT):
        total = np.ndarray(len(self.aspikes))
        for n_it in xrange(len(self.aspikes)):
            n_spikes = self.aspikes[n_it]
            n_spikes_filt = len(n_spikes[np.logical_and(n_spikes >= startT, n_spikes <= startT)])
            total[n_it] = n_spikes_filt / (endT - startT)
        return total

    def torusPopulationVector(self, sheetSize, tstart=0, tend=-1, dt=0.02, winLen=1.0):
        (F, tsteps) = self.getFiringRate(tstart, tend, dt, winLen)
        return torusPopulationVectorFromRates((F, tsteps), sheetSize)



def matFormatSaver(fileName, items, do_compression=False):
    '''Saves either instances of StateMonitor or Spikemonitor, or objects
    saveable by scipy.io.savemat. items is a list of pairs (name, obj)'''

    outData = {}

    for name, obj in items:
        if isinstance(obj, SpikeMonitor):
            outData[name + '_spikeCell'] = obj.aspikes
        elif isinstance(obj, StateMonitor):
            outData[name + '_times'] = obj.times
            outData[name + '_values'] = obj.values
        else:
            outData[name] = obj

    outData['timeSnapshot'] = datetime.now().strftime("%Y-%m-%dT%H-%M-%S")

    savemat(fileName, outData, do_compression=do_compression)

