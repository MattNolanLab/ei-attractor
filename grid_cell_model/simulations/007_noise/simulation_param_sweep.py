#
#   simulation_param_sweep.py
#
#   Main simulation run: parameter sweep runs (noise)
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

from matplotlib.backends.backend_pdf import PdfPages

from scipy      import linspace
from scipy.io   import loadmat
from scipy.io   import savemat
from optparse   import OptionParser

from parameters              import *
from grid_cell_network_brian import *
from custombrian             import *
from tools                   import *
from plotting                import *

import time
import math
import sys
import numpy as np
import logging as lg
import random


lg.basicConfig(level=lg.DEBUG)

parser          = getOptParser()
parser.add_option("--theta_start_mon_t",  type="float",    help="theta start monitoring time")

(options, args) = parser.parse_args()
options         = setOptionDictionary(parser, options)

# Other
figSize = (12,8)
histFigSize = (6, 4)

def correlogram_local(T1, T2, width=20 * ms, bin=1 * ms, T=None):
    '''
    Returns a cross-correlogram with lag in [-width,width] and given bin size.
    T is the total duration (optional) and should be greater than the duration of T1 and T2.
    The result is in Hz (rate of coincidences in each bin).

    N.B.: units are discarded.
    TODO: optimise?
    N.B.1: Adapted from the Brian simulator code
    '''
    if (len(T1) == 0) or (len(T2) == 0): # empty spike train
        return NaN, NaN
    # Remove units
    width = float(width)
    T1 = array(T1)
    T2 = array(T2)
    i = 0
    j = 0
    n = int(ceil(width / bin)) # Histogram length
    l = []
    for t in T1:
        while i < len(T2) and T2[i] < t - width: # other possibility use searchsorted
            i += 1
        while j < len(T2) and T2[j] < t + width:
            j += 1
        l.extend(T2[i:j] - t)
    H, bins = histogram(l, bins=arange(2 * n + 1) * bin - n * bin) #, new = True)

    # Divide by time to get rate
    if T is None:
        T = max(T1[-1], T2[-1]) - min(T1[0], T2[0])
    # Windowing function (triangle)
    W = zeros(2 * n)
    W[:n] = T - bin * arange(n - 1, -1, -1)
    W[n:] = T - bin * arange(n)

    return H / W, bins


def plotHistBars(binEdges, hist):
    binCenters = 0.5*(binEdges[1:] + binEdges[:-1])
    bar(binCenters, hist, width=binCenters[1]-binCenters[0])

################################################################################
#                              Network setup
################################################################################
print "Starting network and connections initialization..."
start_time=time.time()
total_start_t = time.time()

options.ndim = 'twisted_torus'
ei_net = BrianGridCellNetwork(options, simulationOpts=None)
ei_net.uniformInhibition()
ei_net.uniformExcitation()
ei_net.setConstantCurrent()
ei_net.setStartCurrent()
ei_net.setThetaCurrentStimulation()

duration=time.time()-start_time
print "Network setup time:",duration,"seconds"
#                            End Network setup
################################################################################

simulationClock = ei_net._getSimulationClock()

rec_all_spikes = True
if rec_all_spikes:
    nrecSpike_e = ei_net.Ne_x*ei_net.Ne_y
    nrecSpike_i = ei_net.Ni_x*ei_net.Ni_y
else:
    nrecSpike_e = 200
    nrecSpike_i = 50

state_record_e = [ei_net.Ne_x/2 - 1, ei_net.Ne_y/2*ei_net.Ne_x + ei_net.Ne_x/2 - 1]
state_record_i = [ei_net.Ni_x/2 - 1, ei_net.Ni_y/2*ei_net.Ni_x + ei_net.Ni_x/2 - 1]

spikeMon_e          = ExtendedSpikeMonitor(ei_net.E_pop[0:nrecSpike_e])
spikeMon_i          = ExtendedSpikeMonitor(ei_net.I_pop[0:nrecSpike_i])

stateMon_e          = StateMonitor(ei_net.E_pop, 'vm',     record = state_record_e, clock=simulationClock)
stateMon_i          = StateMonitor(ei_net.I_pop, 'vm',     record = state_record_i, clock=simulationClock)
stateMon_ge_e       = StateMonitor(ei_net.E_pop, 'ge',     record = state_record_e, clock=simulationClock)
stateMon_Iclamp_e   = StateMonitor(ei_net.E_pop, 'Iclamp', record = state_record_e, clock=simulationClock)
stateMon_Iclamp_i   = StateMonitor(ei_net.I_pop, 'Iclamp', record = state_record_i, clock=simulationClock)
stateMon_Iext_e     = StateMonitor(ei_net.E_pop, 'Iext',   record = state_record_e, clock=simulationClock)
stateMon_Iext_i     = StateMonitor(ei_net.I_pop, 'Iext',   record = state_record_i, clock=simulationClock)


ei_net.net.add(spikeMon_e, spikeMon_i)
ei_net.net.add(stateMon_e, stateMon_i, stateMon_Iclamp_e, stateMon_Iclamp_i)
ei_net.net.add(stateMon_ge_e)
ei_net.net.add(stateMon_Iext_e, stateMon_Iext_i)


#x_lim = [options.time-0.5, options.time]
x_lim = [options.time/1e3 - 1, options.time/1e3]

################################################################################
#                              Main cycle
################################################################################
for trial_it in range(ei_net.no.ntrials):
    print "Starting trial no. " + str(trial_it) + "..."
    print "Simulation running..."
    start_time=time.time()


    print "  Network initialisation..."
    ei_net.net.run(options.theta_start_mon_t*msecond, report='stdout')

    spikeMon_e.reinit()
    spikeMon_i.reinit()
    stateMon_e.reinit()
    stateMon_i.reinit()
    stateMon_Iclamp_e.reinit()
    stateMon_Iclamp_i.reinit()
    stateMon_ge_e.reinit()
    stateMon_Iext_e.reinit()
    stateMon_Iext_i.reinit()

    print "  Theta stimulation..."
    ei_net.net.run((options.time - options.theta_start_mon_t)*msecond, report='stdout')
    duration=time.time()-start_time
    print "Simulation time:",duration,"seconds"
    
    
    output_fname = "{0}/{1}job{2:04}_trial{3:04}".format(options.output_dir,
            options.fileNamePrefix, options.job_num, trial_it)


    F_tstart = options.theta_start_mon_t*1e-3
    F_tend = options.time*1e-3
    F_dt = 0.05
    F_winLen = 0.25
    Fe, Fe_t = spikeMon_e.getFiringRate(F_tstart, F_tend, F_dt, F_winLen) 
    Fi, Fi_t = spikeMon_i.getFiringRate(F_tstart, F_tend, F_dt, F_winLen)



    figure()
    ax = subplot(211)
    plot(stateMon_e.times, stateMon_e.values[0:2, :].T/mV)
    ylabel('E membrane potential (mV)')
    subplot(212, sharex=ax)
    plot(stateMon_i.times, stateMon_i.values[0:2, :].T/mV)
    xlabel('Time (s)')
    ylabel('I membrane potential (mV)')
    xlim(x_lim)
    savefig(output_fname + '_Vm.pdf')
    

    figure()
    plot(stateMon_ge_e.times, stateMon_ge_e.values[0:2, :].T/nS)
    ylabel('E cell ge (nS)')
    xlabel('Time (s)')
    xlim(x_lim)
    savefig(output_fname + '_ge.pdf')
    

    
    figure()
    ax = subplot(211)
    plot(stateMon_Iclamp_e.times, stateMon_Iclamp_e.values[0:2, :].T/pA)
    ylabel('E synaptic current (pA)')
    #ylim([0, 3000])
    subplot(212, sharex=ax)
    plot(stateMon_Iclamp_i.times, stateMon_Iclamp_i.values[0:2, :].T/pA)
    xlabel('Time (s)')
    ylabel('I synaptic current (pA)')
    xlim(x_lim)
    savefig(output_fname + '_Isyn.pdf')
    
    figure()
    ax = subplot(211)
    plot(stateMon_Iext_e.times, -stateMon_Iext_e.values[1, :].T/pA)
    ylabel('E external current (pA)')
    subplot(212, sharex=ax)
    plot(stateMon_Iext_i.times, -stateMon_Iext_i.values[0, :].T/pA)
    xlabel('Time (s)')
    ylabel('I external current (pA)')
    xlim(x_lim)
    savefig(output_fname + '_Iext.png')
    
    figure()
    pcolormesh(np.reshape(Fe[:, len(Fe_t)/2], (ei_net.Ne_y, ei_net.Ne_x)))
    xlabel('E neuron no.')
    ylabel('E neuron no.')
    colorbar()
    axis('equal')
    savefig(output_fname + '_firing_snapshot_e.png')

    figure()
    pcolormesh(np.reshape(Fi[:, len(Fi_t)/2], (ei_net.Ni_y, ei_net.Ni_x)))
    xlabel('I neuron no.')
    ylabel('I neuron no.')
    colorbar()
    axis('equal')
    savefig(output_fname + '_firing_snapshot_i.png')



    ###################################################################### 
    #                  
    ###################################################################### 



    ei_net.reinit()
#                            End main cycle
################################################################################

total_time = time.time()-total_start_t
print "Overall time: ", total_time, " seconds"

