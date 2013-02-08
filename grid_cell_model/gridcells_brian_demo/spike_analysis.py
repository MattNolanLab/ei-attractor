#
#   spike_analysis.py
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
import numpy as np


def firingRate(spikeTrain, tstart, tend):
    '''
    Compute an average firing rate of a spike train, between tstart and tend
    spikeTrain      A list or a 1D array of spike times
    tstart, tend    Must be the same units as spikeTrain
    '''
    return np.sum(np.logical_and(spikeTrain >= tstart, spikeTrain <= tend)) / \
            (tend - tstart)


def multipleFiringRate(spikeTrainArray, tstart, tend):
    '''
    Compute an average firing rate of an array of spike trains. For each spike
    train, return the firing rate between tstart and tend. Thus, the result will
    be of len(spikeTrainArray)
    '''
    aLen = len(spikeTrainArray)
    result = np.ndarray((aLen, ))
    for it in xrange(aLen):
        result[it] = firingRate(spikeTrainArray[it].flatten(), tstart, tend)

    return result


def spikeMonitorFiringRate(spikeMonitor, tstart, tend):
    N = len(spikeMonitor.source)
    firingRate = np.ndarray((N, ))
    for n_it in xrange(N):
        spikes = spikeMonitor[n_it]
        firingRate[n_it] = np.sum(np.logical_and(spikes >= tstart, spikes <= tend)) / (tend - tstart)

    return firingRate




def firingRateSlidingWindow(spikeTrain, tstart, tend, dt, winLen):
    '''
    Compute sliding window firing rate for the neuron
    dt      resolution of firing rate
    winLen  Length of the sliding window

    The spikes are computed from tstart to tend, so that the resulting array
    length is int((tend-tstart)/dt)+1 long.
    dt does not have to be relevant to simulation dt at all
    '''
    #lg.debug('Start firing rate processing')
    szRate = int((tend-tstart)/dt)+1
    r = np.ndarray((len(spikeTrain), szRate))
    times = np.ndarray(szRate)
    for n_i in xrange(len(spikeTrain)):
        tmp = np.array(spikeTrain[n_i])
        for t_i in xrange(szRate):
            t = tstart + t_i*dt
            r[n_i][t_i] = np.sum(np.logical_and(tmp > t-winLen/2, tmp <
                t+winLen/2))
            times[t_i] = t

    #lg.debug('End firing rate processing')
    return (r/winLen, times)


