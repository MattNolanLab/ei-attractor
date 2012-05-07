#
#   simulation.py
#
#   Main simulation run.
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

import time
import math
import sys
import numpy as np
import logging as lg


lg.basicConfig(level=lg.DEBUG)

parser = getOptParser()

parser.add_option("--Ivel",         type="float", help="Velocity input (pA)")
parser.add_option("--pAMPA_mu",     type="float", help="AMPA profile center (normalised)")
parser.add_option("--pAMPA_sigma",  type="float", help="AMPA profile spread (normalised)")
parser.add_option("--pGABA_sigma",  type="float", help="GABA A profile spread (normalised)")
parser.add_option("--NMDA_amount",  type="float", help="NMDA portion relative to AMPA (%)")
parser.add_option("--Iext_e_min",   type=float,
        help="Minimal external current onto E cells (theta stim.) (A)")
parser.add_option("--Iext_i_min", type=float,
        help="Minimal external current onto I cells (theta stim.) (I)")
parser.add_option("--g_extraGABA_total", type=float,
        help="Uniform inhibition (E-->I only) total conductance (S)")
parser.add_option("--extraGABA_density", type=float,
        help="Uniform inhibition (E-->I only) connection density")
parser.add_option("--prefDirC", type=float, help="Preferred directtion multiplier")

(options, args) = parser.parse_args()
options = setOptionDictionary(parser, options)

# Other
figSize = (12,8)


################################################################################
#                              Network setup
################################################################################
print "Starting network and connections initialization..."
start_time=time.time()
total_start_t = time.time()

options.ndim = 'twisted_torus'
ei_net = BrianGridCellNetwork(options, simulationOpts=None)
ei_net.setConstantCurrent()
ei_net.setStartCurrent()


duration=time.time()-start_time
print "Network setup time:",duration,"seconds"
#                            End Network setup
################################################################################


state_record_e = [ei_net.Ne_x/2 -1 , ei_net.Ne_y/2*ei_net.Ne_x + ei_net.Ne_x/2 - 1]
state_record_i = [ei_net.Ni_x/2 - 1, ei_net.Ni_y/2*ei_net.Ni_x + ei_net.Ni_x/2 - 1]

spikeMon_e = ExtendedSpikeMonitor(ei_net.E_pop)
spikeMon_i = ExtendedSpikeMonitor(ei_net.I_pop)

stateMon_e = StateMonitor(ei_net.E_pop, 'vm', record = state_record_e, clock=simulationClock)
stateMon_i = StateMonitor(ei_net.I_pop, 'vm', record = state_record_i, clock=simulationClock)
stateMon_Iclamp_e = StateMonitor(ei_net.E_pop, 'Iclamp', record = state_record_e, clock=simulationClock)
stateMon_Iclamp_i = StateMonitor(ei_net.I_pop, 'Iclamp', record = state_record_i, clock=simulationClock)
stateMon_Iext_e = StateMonitor(ei_net.E_pop, 'Iext', record = state_record_e, clock=simulationClock)
stateMon_Iext_i = StateMonitor(ei_net.I_pop, 'Iext', record = state_record_i, clock=simulationClock)


ei_net.net.add(spikeMon_e, spikeMon_i)
ei_net.net.add(stateMon_e, stateMon_i, stateMon_Iclamp_e, stateMon_Iclamp_i)
ei_net.net.add(stateMon_Iext_e, stateMon_Iext_i)


#x_lim = [options.time-0.5, options.time]
x_lim = [0, options.time]

################################################################################
#                              Main cycle
################################################################################
for trial_it in range(ei_net.o.ntrials):
    print "Starting trial no. " + str(trial_it) + "..."
    print "Simulation running..."
    start_time=time.time()
    
    print "  Network initialisation..."
    ei_net.net.run(theta_start_mon_t, report='stdout',
            report_period=options.update_interval*second)

    theta_spikeMon_e.reinit()
    theta_spikeMon_i.reinit()
    theta_stateMon_Iclamp_e.reinit()
    theta_stateMon_Iclamp_i.reinit()

    print "  Theta stimulation..."
    ei_net.net.run(options.time*second - theta_start_mon_t, report='stdout',
            report_period=options.update_interval*second)
    duration=time.time()-start_time
    print "Simulation time:",duration,"seconds"
    
    
    output_fname = "{0}/{1}job{2:04}_trial{3:04}".format(options.output_dir,
            options.fileNamePrefix, options.job_num, trial_it)
    
    
    F_tstart = 0
    F_tend = options.time
    F_dt = 0.02
    F_winLen = 0.25
    Fe, Fe_t = spikeMon_e.getFiringRate(F_tstart, F_tend, F_dt, F_winLen) 
    Fi, Fi_t = spikeMon_i.getFiringRate(F_tstart, F_tend, F_dt, F_winLen)

    # plot firing rates
    figure(figsize=figSize)
    subplot(211)
    T, FR = np.meshgrid(Fe_t, np.arange(ei_net.net_Ne))
    pcolormesh(T, FR, Fe)
    ylabel('E Neuron no.')
    colorbar()
    subplot(212)
    T, FR = np.meshgrid(Fi_t, np.arange(ei_net.net_Ni))
    pcolormesh(T, FR, Fi)
    xlabel('Time (s)')
    ylabel('I Neuron no.')
    colorbar()
    savefig(output_fname + '_firing_rate_e.png')

    figure()
    ax = subplot(211)
    plot(stateMon_e.times, stateMon_e.values[0:2].T/mV)
    ylabel('E membrane potential (mV)')
    subplot(212, sharex=ax)
    plot(stateMon_i.times, stateMon_i.values[0:2].T/mV)
    xlabel('Time (s)')
    ylabel('I membrane potential (mV)')
    xlim(x_lim)
    savefig(output_fname + '_Vm.pdf')
    
    
    figure()
    ax = subplot(211)
    plot(stateMon_Iclamp_e.times, stateMon_Iclamp_e.values[0:2].T/pA)
    ylabel('E synaptic current (pA)')
    subplot(212, sharex=ax)
    plot(stateMon_Iclamp_i.times, stateMon_Iclamp_i.values[0:2].T/pA)
    xlabel('Time (s)')
    ylabel('I synaptic current (pA)')
    xlim(x_lim)
    savefig(output_fname + '_Isyn.pdf')
    
    figure()
    ax = subplot(211)
    plot(stateMon_Iext_e.times, -stateMon_Iext_e.values[1].T/pA)
    ylabel('E external current (pA)')
    subplot(212, sharex=ax)
    plot(stateMon_Iext_i.times, -stateMon_Iext_i.values[0].T/pA)
    xlabel('Time (s)')
    ylabel('I external current (pA)')
    xlim(x_lim)
    savefig(output_fname + '_Iext.pdf')
    
    # High pass filter these signals
    figure()
    ax = subplot(211)
    plot(stateMon_Iclamp_e.times, butterHighPass(stateMon_Iclamp_e.values[1].T/pA, options.sim_dt, 40))
    plot(stateMon_Iext_e.times, -(stateMon_Iext_e.values[0]/pA - stim_e_DC/pA))
    ylabel('E current (pA)')
    ylim([-500, 500])
    subplot(212, sharex=ax)
    plot(stateMon_Iclamp_i.times, butterHighPass(stateMon_Iclamp_i.values[0].T/pA, options.sim_dt, 40))
    #plot(stateMon_Iclamp_i.times, stateMon_Iext_i.values[0]/pA)
    xlabel('Time (s)')
    ylabel('I current (pA)')
    xlim(x_lim)
    ylim([-500, 500])
    savefig(output_fname + '_Isyn_filt.pdf')
    
    
    
    Ne = options.Ne
    figure()
    pcolormesh(np.reshape(ei_net.AMPA_conn.W.todense()[57*68 + 33, :], (ei_net.Ni_y,
        ei_net.Ni_x)));
    xlabel('I neuron no.')
    ylabel('I neuron no.')
    colorbar()
    axis('equal')
    savefig(output_fname + '_E2I_conn.png')

    Ni = options.Ni
    figure()
    pcolormesh(np.reshape(ei_net.GABA_conn1.W.todense()[0, :], (ei_net.Ne_y,
        ei_net.Ne_x)));
    xlabel('E neuron no.')
    ylabel('E neuron no.')
    colorbar()
    axis('equal')
    savefig(output_fname + '_I2E_conn.png')

    #figure()
    #pcolormesh(np.reshape(np.dot(ei_net.AMPA_conn.W.todense(),
    #    ei_net.GABA_conn1.W.todense())[15, :], (ei_net.Ne_y, ei_net.Ne_x)));
    #xlabel('E neuron no.')
    #ylabel('E neuron no.')
    #colorbar()
    #savefig(output_fname + '_E2E_conn.png')


    figure()
    pcolormesh(np.reshape(Fe[:, len(Fe_t)/2], (ei_net.Ne_y, ei_net.Ne_x)))
    xlabel('E neuron no.')
    ylabel('E neuron no.')
    colorbar()
    axis('equal')
    savefig(output_fname + '_firing_snapshot_e.png')


    
    # Print a plot of bump position
    (pos, bumpPos_times) = theta_spikeMon_e.torusPopulationVector([ei_net.Ne_x,
        ei_net.Ne_y], theta_start_t, options.time, F_dt, F_winLen)
    figure(figsize=figSize)
    plot(bumpPos_times, pos)
    xlabel('Time (s)')
    ylabel('Bump position (neurons)')
    ylim([-ei_net.Ne_x/2 -5, ei_net.Ne_y/2 + 5])
    
    savefig(output_fname + '_bump_position.pdf')



    print "End of trial no. " + str(trial_it) + "..."
    print 

    #ei_net.reinit()
#                            End main cycle
################################################################################

total_time = time.time()-total_start_t
print "Overall time: ", total_time, " seconds"

