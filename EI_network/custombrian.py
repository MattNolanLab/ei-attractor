from brian.monitor import SpikeMonitor

import numpy as np
import operator as op
import logging as lg


class ExtendedSpikeMonitor(SpikeMonitor):
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

    spiketimes = property(fget=getspiketimes)
    spikes = property(fget=getSpikes, fset=_setSpikes)

