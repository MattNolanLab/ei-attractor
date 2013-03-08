#
#   simulation_step_current.py
#
#   Main simulation run: step current simulations
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

from models.parameters       import *
from models.gc_net_brian     import *
from custombrian             import *

import time
import math
import sys
import numpy as np
import logging as lg


lg.basicConfig(level=lg.DEBUG)

parser          = getOptParser()
parser.add_option("--ngenerations",      type="int",    help="Number of network generation cycles")
parser.add_option("--stateMonDuration",  type="float",  help="Time windown for the state monitor (ms)")
parser.add_option("--velModulationType", type="string", help="Type of velocity modulation (excitatory, inhibitory, weights, conjunctive)")
parser.add_option("--stepStart_t",       type="float",  help="Step current start time (ms)")
parser.add_option("--stepEnd_t",         type="float",  help="Step current end time (ms)")

(options, args) = parser.parse_args()
options         = setOptionDictionary(parser, options)

# Other
figSize = (12,8)


for gen_it in range(options.ngenerations):
    print "Network generation no. " + str(gen_it)
    ################################################################################
    #                              Network setup
    ################################################################################
    print "Starting network and connections initialization..."
    start_time=time.time()
    total_start_t = time.time()

    ei_net = BrianGridCellNetwork(options, simulationOpts=None)
    ei_net.setConstantCurrent()
    ei_net.setStartCurrent()
    ei_net.setThetaCurrentStimulation()
    const_v = [1.0, 0.0]
    if options.velModulationType == "excitatory":
        #Modulation of excitatory neurons (weight shifts are E-->I)
        ei_net.setConstantVelocityCurrent_e(const_v, start_t=options.stepStart_t, end_t=options.stepEnd_t)
    elif options.velModulationType == "inhibitory":
        #Modulation of inhibitory neurons (weight shifts are I-->E)
        ei_net.setConstantVelocityCurrent_i(const_v, start_t=options.stepStart_t, end_t=options.stepEnd_t)
    else:
        raise Exception("Unknown velocity modulation type")
    
    
    duration=time.time()-start_time
    print "Network setup time:",duration,"seconds"
    #                            End Network setup
    ################################################################################
    
    simulationClock = ei_net._getSimulationClock()
    
    
    state_record_e = [ei_net.Ne_x/2 -1 , ei_net.Ne_y/2*ei_net.Ne_x + ei_net.Ne_x/2 - 1]
    state_record_i = [ei_net.Ni_x/2 - 1, ei_net.Ni_y/2*ei_net.Ni_x + ei_net.Ni_x/2 - 1]
    
    spikeMon_e          = ExtendedSpikeMonitor(ei_net.E_pop)
    spikeMon_i          = ExtendedSpikeMonitor(ei_net.I_pop)
    
    stateMon_e          = StateMonitor(ei_net.E_pop, 'vm',       record = state_record_e, clock=simulationClock)
    stateMon_i          = StateMonitor(ei_net.I_pop, 'vm',       record = state_record_i, clock=simulationClock)
    stateMon_Iclamp_e   = StateMonitor(ei_net.E_pop, 'Iclamp',   record = state_record_e, clock=simulationClock)
    stateMon_Iclamp_i   = StateMonitor(ei_net.I_pop, 'Iclamp',   record = state_record_i, clock=simulationClock)
    stateMon_Iext_e     = StateMonitor(ei_net.E_pop, 'Iext',     record = state_record_e, clock=simulationClock)
    stateMon_Iext_vel_e = StateMonitor(ei_net.E_pop, 'Iext_vel', record = state_record_i, clock=simulationClock)
    
    ei_net.net.add(spikeMon_e, spikeMon_i)
    ei_net.net.add(stateMon_e, stateMon_i, stateMon_Iclamp_e, stateMon_Iclamp_i)
    ei_net.net.add(stateMon_Iext_e, stateMon_Iext_vel_e)
    
    
    #x_lim = [options.time-0.5, options.time]
    x_lim = [options.time/1e3 - 1, options.time/1e3]
    
    ################################################################################
    #                              Main cycle
    ################################################################################
    print "Simulation running..."
    start_time=time.time()
    
    ei_net.net.run(options.time*msecond, report='stdout')
    duration=time.time()-start_time
    print "Simulation time:",duration,"seconds"
    
    
    output_fname = "{0}/{1}job{2:04}_gen{3:04}".format(options.output_dir,
            options.fileNamePrefix, options.job_num, gen_it)
    

    # Print a plot of bump position
    F_dt = 0.02
    F_winLen = 0.25
    (pos, bumpPos_times) = spikeMon_e.torusPopulationVector([ei_net.Ne_x,
        ei_net.Ne_y], options.theta_start_t*1e-3, options.time*1e-3, F_dt, F_winLen)
    figure(figsize=figSize)
    plot(bumpPos_times, pos)
    xlabel('Time (s)')
    ylabel('Bump position (neurons)')
    ylim([-ei_net.Ne_x/2 -5, ei_net.Ne_y/2 + 5])
    
    savefig(output_fname + '_bump_position.pdf')

    
    outData = {}
    outData['bumpPos']                    = pos
    outData['bumpPos_times']              = bumpPos_times
    outData['options']                    = options._einet_optdict
    outData['stepStart_t']                = options.stepStart_t
    outData['stepEnd_t']                  = options.stepEnd_t
    outData['Ivel']                       = options.Ivel
    outData['stateMon_Iext_vel_e_times']  = stateMon_Iext_vel_e.times
    outData['stateMon_Iext_vel_e_values'] = stateMon_Iext_vel_e.values
    savemat(output_fname + '_output.mat', outData, do_compression=True)


    print "End of generation no. " + str(gen_it) + "..."
    print 
#                            End main cycle
################################################################################

total_time = time.time()-total_start_t
print "Overall time: ", total_time, " seconds"

