#
#   simulation_GABA_rev.py
#
#   Main simulation run: Change of GABA_A reversal potential
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
from models.gc_net_brian     import *
from custombrian             import *

import time
import math
import sys
import numpy as np
import logging as lg


lg.basicConfig(level=lg.DEBUG)

parser          = getOptParser()
parser.add_option("--ndumps",  type="int",    help="Number of data output dumps during the simulation")

(options, args) = parser.parse_args()
options         = setOptionDictionary(parser, options)

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
ei_net.setThetaCurrentStimulation()
ei_net.setPlaceCurrentInput()

#const_v = [0.0, 1.0]
#ei_net.setConstantVelocityCurrent_e(const_v)
ei_net.setVelocityCurrentInput_e()

duration=time.time()-start_time
print "Network setup time:",duration,"seconds"
#                            End Network setup
################################################################################

simulationClock = ei_net._getSimulationClock()

rec_all_spikes = False
if rec_all_spikes:
    nrecSpike_e = ei_net.Ne_x*ei_net.Ne_y
    nrecSpike_i = ei_net.Ni_x*ei_net.Ni_y
else:
    nrecSpike_e = 200
    nrecSpike_i = 50


state_record_e = [ei_net.Ne_x/2 -1 , ei_net.Ne_y/2*ei_net.Ne_x + ei_net.Ne_x/2 - 1]
state_record_i = [ei_net.Ni_x/2 - 1, ei_net.Ni_y/2*ei_net.Ni_x + ei_net.Ni_x/2 - 1]

spikeMon_e          = ExtendedSpikeMonitor(ei_net.E_pop[0:nrecSpike_e])
spikeMon_i          = ExtendedSpikeMonitor(ei_net.I_pop[0:nrecSpike_i])

stateMon_e          = RecentStateMonitor(ei_net.E_pop, 'vm',     duration=options.stateMonDur*ms,   record = state_record_e, clock=simulationClock)
stateMon_i          = RecentStateMonitor(ei_net.I_pop, 'vm',     duration=options.stateMonDur*ms,   record = state_record_i, clock=simulationClock)
stateMon_Iclamp_e   = RecentStateMonitor(ei_net.E_pop, 'Iclamp', duration=options.stateMonDur*ms,   record = state_record_e, clock=simulationClock)
stateMon_Iclamp_i   = RecentStateMonitor(ei_net.I_pop, 'Iclamp', duration=options.stateMonDur*ms,   record = state_record_i, clock=simulationClock)
stateMon_Iext_e     = RecentStateMonitor(ei_net.E_pop, 'Iext',   duration=options.stateMonDur*ms,   record = state_record_e, clock=simulationClock)
stateMon_Iext_i     = RecentStateMonitor(ei_net.I_pop, 'Iext',   duration=options.stateMonDur*ms,   record = state_record_i, clock=simulationClock)

ei_net.net.add(spikeMon_e, spikeMon_i)
ei_net.net.add(stateMon_e, stateMon_i, stateMon_Iclamp_e, stateMon_Iclamp_i)
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
    
    for dump_it in range(options.ndumps):
        ei_net.net.run(options.time/options.ndumps*msecond, report='stdout')
        duration=time.time()-start_time
        print "Simulation time:",duration,"seconds"
        
        
        output_fname = "{0}/{1}job{2:04}_trial{3:04}_dump{4:03}".format(options.output_dir,
                options.fileNamePrefix, options.job_num, trial_it, dump_it)


        #F_tstart = 0
        #F_tend = options.time*1e-3
        #F_dt = 0.02
        #F_winLen = 0.25
        #Fe, Fe_t = spikeMon_e.getFiringRate(F_tstart, F_tend, F_dt, F_winLen) 
        #Fi, Fi_t = spikeMon_i.getFiringRate(F_tstart, F_tend, F_dt, F_winLen)


        ## plot firing rates
        #figure(figsize=figSize)
        #subplot(211)
        #T, FR = np.meshgrid(Fe_t, np.arange(ei_net.net_Ne))
        #pcolormesh(T, FR, Fe)
        #ylabel('E Neuron no.')
        #colorbar()
        #subplot(212)
        #T, FR = np.meshgrid(Fi_t, np.arange(ei_net.net_Ni))
        #pcolormesh(T, FR, Fi)
        #xlabel('Time (s)')
        #ylabel('I Neuron no.')
        #colorbar()
        #savefig(output_fname + '_firing_rate_e.png')

        #figure()
        #ax = subplot(211)
        #plot(stateMon_e.times, stateMon_e.values[:, 0:2]/mV)
        #ylabel('E membrane potential (mV)')
        #subplot(212, sharex=ax)
        #plot(stateMon_i.times, stateMon_i.values[:, 0:2]/mV)
        #xlabel('Time (s)')
        #ylabel('I membrane potential (mV)')
        #xlim(x_lim)
        #savefig(output_fname + '_Vm.pdf')
        #
        #
        #figure()
        #ax = subplot(211)
        #plot(stateMon_Iclamp_e.times, stateMon_Iclamp_e.values[:, 0:2]/pA)
        #ylabel('E synaptic current (pA)')
        #subplot(212, sharex=ax)
        #plot(stateMon_Iclamp_i.times, stateMon_Iclamp_i.values[:, 0:2]/pA)
        #ylim([-1000, 0])
        #xlabel('Time (s)')
        #ylabel('I synaptic current (pA)')
        #xlim(x_lim)
        #savefig(output_fname + '_Isyn.pdf')
        #
        #figure()
        #ax = subplot(211)
        #plot(stateMon_Iext_e.times, -stateMon_Iext_e.values[:, 1]/pA)
        #ylabel('E external current (pA)')
        #subplot(212, sharex=ax)
        #plot(stateMon_Iext_i.times, -stateMon_Iext_i.values[:, 0]/pA)
        #xlabel('Time (s)')
        #ylabel('I external current (pA)')
        #xlim(x_lim)
        #savefig(output_fname + '_Iext.pdf')
    
        #figure()
        #pcolormesh(np.reshape(Fe[:, len(Fe_t)/2], (ei_net.Ne_y, ei_net.Ne_x)))
        #xlabel('E neuron no.')
        #ylabel('E neuron no.')
        #colorbar()
        #axis('equal')
        #savefig(output_fname + '_firing_snapshot_e.png')


        ## Print a plot of bump position
        #F_dt = 0.02
        #F_winLen = 0.25
        #(pos, bumpPos_times) = spikeMon_e.torusPopulationVector([ei_net.Ne_x,
        #    ei_net.Ne_y], options.theta_start_t*1e-3, options.time*1e-3, F_dt, F_winLen)
        #figure(figsize=figSize)
        #plot(bumpPos_times, pos)
        #xlabel('Time (s)')
        #ylabel('Bump position (neurons)')
        #ylim([-ei_net.Ne_x/2 -5, ei_net.Ne_y/2 + 5])
        #
        #savefig(output_fname + '_bump_position.pdf')

        
        outData = ei_net.getRatData()
        #outData['timeSnapshot'] = datetime.now().strftime("%Y-%m-%dT%H-%M-%S")

        #outData['bumpPos'] = pos
        #outData['bumpPos_times'] = bumpPos_times

        #outData['Fe'] = Fe
        #outData['Fe_t'] = Fe_t

        outData['spikeCell_e']              = spikeMon_e.aspikes[0:200]
        outData['spikeCell_i']              = spikeMon_i.aspikes[0:50]
        outData['options']                  = options._einet_optdict
        outData['velocityStart']            = options.theta_start_t

        outData['stateMon_times']           = stateMon_Iclamp_e.times
        outData['stateMon_Iclamp_e_values'] = stateMon_Iclamp_e.values
        outData['stateMon_Iclamp_i_values'] = stateMon_Iclamp_i.values
        outData['stateMon_e_values']        = stateMon_e.values
        outData['stateMon_i_values']        = stateMon_i.values
        
        savemat(output_fname + '_output.mat', outData, do_compression=False)

        print "Dump after " + str(simulationClock.t)


    print "End of trial no. " + str(trial_it) + "..."
    print 

    ei_net.reinit()
#                            End main cycle
################################################################################

total_time = time.time()-total_start_t
print "Overall time: ", total_time, " seconds"

