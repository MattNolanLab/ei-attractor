from brian.monitor import SpikeMonitor

import numpy as np
import operator as op
import logging as lg


class NeuronSpikes:
    def __init__(self, aspikes):
        self.aspikes = aspikes

    def getFiringRate(self, tstart, tend, dt, winLen):
        '''
        Compute sliding window firing rate for the neuron
        dt      resolution of firing rate
        winLen  Length of the sliding window

        The spikes are computed from tstart to tend, so that the resulting array
        length is int((tend-tstart)/dt)+1 long.
        dt does not have to be relevant to simulation dt at all
        '''
        szRate = int((tend-tstart)/dt)+1
        r = np.ndarray((len(self.aspikes), szRate))
        times = np.ndarray(szRate)
        for n_i in xrange(len(self.aspikes)):
            tmp = np.array(self.aspikes[n_i])
            t = tstart
            for t_i in xrange(szRate):
                t = t_i*dt
                r[n_i][t_i] = np.sum(np.logical_and(tmp > t-winLen/2, tmp <
                    t+winLen/2))
                times[t_i] = t

        return (r/winLen, times)



class ExtendedSpikeMonitor(SpikeMonitor):
    '''
    SpikeMonitor class that stores spikes implicitly as an array of lists,
    instead of a list of pairs (n, t) as does the original SpikeMonitor class
    '''
    def __init__(self, source, record=True, delay=0, function=None):
        SpikeMonitor.__init__(self, source, record, delay, function)

        self.aspikes = np.ndarray(len(source), dtype=list)
        self._spiketimes = {}
        self._spikes = None
        self.reinit()

    def reinit(self):
        self._newspikes = True
        self.nspikes = 0
        for i in xrange(len(self.source)):
            self.aspikes[i] = []
            self.spiketimes[i] = self.aspikes[i]

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

    def getspiketimes(self):
        return self._spiketimes

    def getSpikes(self):
        ''' Just to make sure everyone gets what they want'''
        if self._newspikes:
            self._newspikes = False
            self._spikes = []
            for i in xrange(len(self.source)):
                self._spikes += zip([i]*len(self.source), self.aspikes[i])
            lg.debug('Spike sort start...')
            self._spikes.sort(key=op.itemgetter(1))
            lg.debug('Spike sort end.')
        return self._spikes

    def _setSpikes(self, val):
        self._spikes = val

    def getFiringRate(self, tstart, tend, dt, winLen):
        return NeuronSpikes(self.aspikes).getFiringRate(tstart, tend, dt, winLen)

    spiketimes = property(fget=getspiketimes)
    spikes = property(fget=getSpikes, fset=_setSpikes)

