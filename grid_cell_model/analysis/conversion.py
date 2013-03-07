#
#   conversion.py
#
#   Conversion module. Convert between different formats of data, types, etc.
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



## Convert a pair (n_idx, spikeTimes) into an (object) array in which every row
# contains spike times for each neuron in n_idx.
#
# For BIG arrays this will be incredibly slow!!
#
# @param n_idx      An array of neuron indexes.
# @param spikeTimes An array of spike times for each neuron no. in n_idx
# @return An object array, in which every row contains spike times for the
#   corresponding neuron
#
def spikePairsToArray(n_idx, spikeTimes):
    if (n_idx.size == 0):
        return np.ndarray((0,), dtype=object)

    n_sz = np.max(n_idx) + 1
    res = np.ndarray((n_sz, ), dtype=object)

    for n_i in xrange(n_sz):
        res[n_i] = spikeTimes[n_idx == n_i]

    return res



## Simply extract all the spikes of a neuron with index n_idx
#
# @param n_idx      An index of a neuron to extract spikes for.
# @param all_n_idx  Indexes of neurons for the spikeTimes array.
# @param spikeTimes The array of spike times, each value corresponds to a
# @param sorted     Sort the result before returning
#   neuron with an id in all_n_idx
# @return An array of spike times for the corresponding neuron. Sorted if
#   needed.
#
def neuronSpikes(n_idx, all_n_idx, spikeTimes, sorted=False):
    if (sorted):
        return np.sort(spikeTimes[all_n_idx == n_idx])
    else:
        return spikeTimes[all_n_idx == n_idx]



###############################################################################
#                               Tests
###############################################################################
if (__name__ == "__main__"):
    
    N            = np.array([1,   20,   31,   20])
    N_spikeTimes = np.array([1.1, 20.1, 31.1, 20.05])

    spikes = spikePairsToArray(N, N_spikeTimes)

    
    # Pick one neuron only
    # no. 20 - 2 values
    spikes_full = neuronSpikes(20, N, N_spikeTimes)
    print spikes_full
    # empty
    spikes_em   = neuronSpikes(0, N, N_spikeTimes)
    print spikes_em

    # sorted
    spikes_sort = neuronSpikes(20, N, N_spikeTimes, sorted=True)
    print spikes_sort

    ## Big array
    #N      = 4000
    #spPerN = 50
    #T      = 1e3
    #N_id   = np.random.randint(0, N, N*spPerN)
    #spiket = T * np.random.rand(N*spPerN)

    #spikes = spikePairsToArray(N_id, spiket)
    #print spikes


