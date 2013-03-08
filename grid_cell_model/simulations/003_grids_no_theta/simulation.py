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

#const_v = [0.0, 1.0]
#ei_net.setConstantVelocityCurrent_e(const_v)
ei_net.setVelocityCurrentInput_e()

duration=time.time()-start_time
print "Network setup time:",duration,"seconds"
#                            End Network setup
################################################################################

simulationClock = ei_net._getSimulationClock()


state_record_e = [ei_net.Ne_x/2 -1 , ei_net.Ne_y/2*ei_net.Ne_x + ei_net.Ne_x/2 - 1]
state_record_i = [ei_net.Ni_x/2 - 1, ei_net.Ni_y/2*ei_net.Ni_x + ei_net.Ni_x/2 - 1]

spikeMon_e          = ExtendedSpikeMonitor(ei_net.E_pop)
spikeMon_i          = ExtendedSpikeMonitor(ei_net.I_pop)

stateMon_e          = StateMonitor(ei_net.E_pop, 'vm', record = state_record_e, clock=simulationClock)
stateMon_i          = StateMonitor(ei_net.I_pop, 'vm', record = state_record_i, clock=simulationClock)
stateMon_Iclamp_e   = StateMonitor(ei_net.E_pop, 'Iclamp', record = state_record_e, clock=simulationClock)
stateMon_Iclamp_i   = StateMonitor(ei_net.I_pop, 'Iclamp', record = state_record_i, clock=simulationClock)
stateMon_Iext_e     = StateMonitor(ei_net.E_pop, 'Iext', record = state_record_e, clock=simulationClock)
stateMon_Iext_i     = StateMonitor(ei_net.I_pop, 'Iext', record = state_record_i, clock=simulationClock)

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

        outData['stateMon_Iclamp_e_times']  = stateMon_Iclamp_e.times
        outData['stateMon_Iclamp_e_values'] = stateMon_Iclamp_e.values
        outData['stateMon_Iclamp_i_times']  = stateMon_Iclamp_i.times
        outData['stateMon_Iclamp_i_values'] = stateMon_Iclamp_i.values
        
        savemat(output_fname + '_output.mat', outData, do_compression=True)

        print "Dump after " + str(simulationClock.t)


    print "End of trial no. " + str(trial_it) + "..."
    print 

    ei_net.reinit()
#                            End main cycle
################################################################################

total_time = time.time()-total_start_t
print "Overall time: ", total_time, " seconds"

