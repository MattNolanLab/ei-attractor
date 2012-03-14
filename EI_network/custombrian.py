from brian.monitor import SpikeMonitor, StateMonitor
from scipy.io import savemat

import numpy as np
import operator as op
import logging as lg

from datetime import datetime



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
        lg.debug('Start firing rate processing')
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

        lg.debug('End firing rate processing')
        return (r/winLen, times)

    def getNSpikes(self):
        total = np.ndarray(len(self.aspikes))
        for n_it in xrange(len(self.aspikes)):
            total[n_it] = len(self.aspikes[n_it])
        return total

    def torusPopulationVector(self, sheetSize, tstart=0, tend=-1, dt=0.02, winLen=1.0):
        (F, tsteps) = self.getFiringRate(tstart, tend, dt, winLen)
        
        P = np.ndarray((len(tsteps), 2), dtype=complex)
        X, Y = np.meshgrid(np.arange(sheetSize), np.arange(sheetSize))
        X = np.exp(1j*(X - sheetSize/2)/sheetSize*2*np.pi).ravel()
        Y = np.exp(1j*(Y - sheetSize/2)/sheetSize*2*np.pi).ravel()
        for t_it in xrange(len(tsteps)):
            P[t_it, 0] = np.dot(F[:, t_it], X)
            P[t_it, 1] = np.dot(F[:, t_it], Y)

        return (np.angle(P)/2/np.pi*sheetSize, tsteps)
        

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
