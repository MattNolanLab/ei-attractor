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
import numpy as np
import scipy

from scipy import weave

__all__ = ['firingRate', 'multipleFiringRate', 'firingRateFromPairs',
        'firingRateSlidingWindow', 'slidingFiringRateTuple',
        'torusPopulationVector', 'torusPopulationVectorFromRates']


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

def firingRateFromPairs(N, N_ids, times, tstart, tend):
    '''
    Compute average firing rate for all the neurons in the range <0, N), given
    two arrays: N_ids (neuron ids) and times of spikes corresponding to N_ids).
    '''
    result = np.ndarray((N, ))
    for n_it in xrange(N):
        result[n_it] = firingRate(times[N_ids == n_it], tstart, tend)
    return result



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
    winLen  Length of the sliding window (ms)

    return  a n array of shape (N, int((tend-tstart)/dt)+1
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
    N = sheetSize[0]*sheetSize[1]
    F, tsteps = slidingFiringRateTuple(spikes, N, tstart, tend, dt, winLen)

    return torusPopulationVectorFromRates((F, tsteps), sheetSize)


def torusPopulationVectorFromRates(FR, sheetSize):
        F = FR[0]
        tsteps = FR[1]

        sheetSize_x = sheetSize[0]
        sheetSize_y = sheetSize[1]
        
        P = np.ndarray((len(tsteps), 2), dtype=complex)
        X, Y = np.meshgrid(np.arange(sheetSize_x), np.arange(sheetSize_y))
        X = np.exp(1j*(X - sheetSize_x/2)/sheetSize_x*2*np.pi).ravel()
        Y = np.exp(1j*(Y - sheetSize_y/2)/sheetSize_y*2*np.pi).ravel()
        for t_it in xrange(len(tsteps)):
            P[t_it, 0] = np.dot(F[:, t_it], X)
            P[t_it, 1] = np.dot(F[:, t_it], Y)

        return (np.angle(P)/2/np.pi*sheetSize, tsteps)

