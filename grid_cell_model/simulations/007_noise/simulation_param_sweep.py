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
ei_net.uniformExcitation()
ei_net.setConstantCurrent()
ei_net.setStartCurrent()
ei_net.setThetaCurrentStimulation()

duration=time.time()-start_time
print "Network setup time:",duration,"seconds"
#                            End Network setup
################################################################################

simulationClock = ei_net._getSimulationClock()
stateRecordClock = Clock(dt=2*ms)

rec_all_spikes = True
if rec_all_spikes:
    nrecSpike_e = ei_net.Ne_x*ei_net.Ne_y
    nrecSpike_i = ei_net.Ni_x*ei_net.Ni_y
else:
    nrecSpike_e = 200
    nrecSpike_i = 50

NStateRecord_e = 100
NStateRecord_i = 50

state_record_e = random.sample(xrange(ei_net.Ne_x*ei_net.Ne_y), NStateRecord_e)
state_record_i = random.sample(xrange(ei_net.Ni_x*ei_net.Ni_y), NStateRecord_i)

spikeMon_e          = ExtendedSpikeMonitor(ei_net.E_pop[0:nrecSpike_e])
spikeMon_i          = ExtendedSpikeMonitor(ei_net.I_pop[0:nrecSpike_i])

stateMon_e          = StateMonitor(ei_net.E_pop, 'vm',     record = state_record_e, clock=stateRecordClock)
stateMon_i          = StateMonitor(ei_net.I_pop, 'vm',     record = state_record_i, clock=stateRecordClock)
stateMon_Iclamp_e   = StateMonitor(ei_net.E_pop, 'Iclamp', record = state_record_e, clock=stateRecordClock)
stateMon_Iclamp_i   = StateMonitor(ei_net.I_pop, 'Iclamp', record = state_record_i, clock=stateRecordClock)
stateMon_Iext_e     = StateMonitor(ei_net.E_pop, 'Iext',   record = state_record_e, clock=stateRecordClock)
stateMon_Iext_i     = StateMonitor(ei_net.I_pop, 'Iext',   record = state_record_i, clock=stateRecordClock)


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


    print "  Network initialisation..."
    ei_net.net.run(options.theta_start_mon_t*msecond, report='stdout')

    #spikeMon_e.reinit()
    #spikeMon_i.reinit()
    #stateMon_e.reinit()
    #stateMon_i.reinit()
    #stateMon_Iclamp_e.reinit()
    #stateMon_Iclamp_i.reinit()
    #stateMon_ge_e.reinit()
    #stateMon_Iext_e.reinit()
    #stateMon_Iext_i.reinit()

    print "  Theta stimulation..."
    ei_net.net.run((options.time - options.theta_start_mon_t)*msecond, report='stdout')
    duration=time.time()-start_time
    print "Simulation time:",duration,"seconds"
    
    
    output_fname = "{0}/{1}job{2:05}_trial{3:04}".format(options.output_dir,
            options.fileNamePrefix, options.job_num, trial_it)


    #F_tstart = options.theta_start_mon_t*1e-3
    #F_tend = options.time*1e-3
    #F_dt = 0.05
    #F_winLen = 0.25
    #Fe, Fe_t = spikeMon_e.getFiringRate(F_tstart, F_tend, F_dt, F_winLen) 
    #Fi, Fi_t = spikeMon_i.getFiringRate(F_tstart, F_tend, F_dt, F_winLen)



    #figure()
    #ax = subplot(211)
    #plot(stateMon_e.times, stateMon_e.values[0:2, :].T/mV)
    #ylabel('E membrane potential (mV)')
    #subplot(212, sharex=ax)
    #plot(stateMon_i.times, stateMon_i.values[0:2, :].T/mV)
    #xlabel('Time (s)')
    #ylabel('I membrane potential (mV)')
    #xlim(x_lim)
    #savefig(output_fname + '_Vm.pdf')
    #

    #figure()
    #ax = subplot(211)
    #plot(stateMon_Iclamp_e.times, stateMon_Iclamp_e.values[0:2, :].T/pA)
    #ylabel('E synaptic current (pA)')
    ##ylim([0, 3000])
    #subplot(212, sharex=ax)
    #plot(stateMon_Iclamp_i.times, stateMon_Iclamp_i.values[0:2, :].T/pA)
    #xlabel('Time (s)')
    #ylabel('I synaptic current (pA)')
    #xlim(x_lim)
    #savefig(output_fname + '_Isyn.pdf')
    #
    #figure()
    #ax = subplot(211)
    #plot(stateMon_Iext_e.times, -stateMon_Iext_e.values[1, :].T/pA)
    #ylabel('E external current (pA)')
    #subplot(212, sharex=ax)
    #plot(stateMon_Iext_i.times, -stateMon_Iext_i.values[0, :].T/pA)
    #xlabel('Time (s)')
    #ylabel('I external current (pA)')
    #xlim(x_lim)
    #savefig(output_fname + '_Iext.png')
    
    #figure()
    #pcolormesh(np.reshape(Fe[:, len(Fe_t)/2], (ei_net.Ne_y, ei_net.Ne_x)))
    #xlabel('E neuron no.')
    #ylabel('E neuron no.')
    #colorbar()
    #axis('equal')
    #savefig(output_fname + '_firing_snapshot_e.png')

    #figure()
    #pcolormesh(np.reshape(Fi[:, len(Fi_t)/2], (ei_net.Ni_y, ei_net.Ni_x)))
    #xlabel('I neuron no.')
    #ylabel('I neuron no.')
    #colorbar()
    #axis('equal')
    #savefig(output_fname + '_firing_snapshot_i.png')



    ###################################################################### 
    #                           Export data
    ###################################################################### 
    outData = dict()
    #outData['timeSnapshot'] = datetime.now().strftime("%Y-%m-%dT%H-%M-%S")

    outData['spikeCell_e']              = spikeMon_e.aspikes
    outData['spikeCell_i']              = spikeMon_i.aspikes
    outData['options']                  = options._einet_optdict
    outData['theta_start_mon_t']        = options.theta_start_mon_t

    outData['stateMon_times']           = stateMon_Iclamp_e.times
    outData['stateMon_Iclamp_e_values'] = stateMon_Iclamp_e.values
    outData['stateMon_Iclamp_i_values'] = stateMon_Iclamp_i.values
    outData['stateMon_e_values']        = stateMon_e.values
    outData['stateMon_i_values']        = stateMon_i.values
    outData['state_record_e']           = state_record_e
    outData['state_record_i']           = state_record_i
    
    savemat(output_fname + '_output.mat', outData, do_compression=False)


    ei_net.reinit()
#                            End main cycle
################################################################################

total_time = time.time()-total_start_t
print "Overall time: ", total_time, " seconds"

