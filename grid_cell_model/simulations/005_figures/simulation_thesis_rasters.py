'''Run the simulation as in simulation_fig_model, but only export spikes.'''
from scipy      import linspace
from optparse   import OptionParser

from parameters              import *
from grid_cell_network_brian import *
from custombrian             import *
from tools                   import *
from plotting                import *

import time
import numpy as np
import logging as lg
from simtools.storage import DataStorage


lg.basicConfig(level=lg.DEBUG)

parser          = getOptParser()
parser.add_option("--theta_start_mon_t",  type="float",    help="theta start monitoring time")

(options, args) = parser.parse_args()
options         = setOptionDictionary(parser, options)

################################################################################
#                              Network setup
################################################################################
print "Starting network and connections initialization..."
start_time=time.time()
total_start_t = time.time()

options.ndim = 'twisted_torus'
ei_net = BrianGridCellNetwork(options, simulationOpts=None)
ei_net.uniformInhibition()
ei_net.setConstantCurrent()
ei_net.setStartCurrent()
ei_net.setThetaCurrentStimulation()

duration=time.time()-start_time
print "Network setup time:",duration,"seconds"
#                            End Network setup
################################################################################

simulationClock = ei_net._getSimulationClock()

spikeMon_e = SpikeMonitor(ei_net.E_pop)
spikeMon_i = SpikeMonitor(ei_net.I_pop)

ei_net.net.add(spikeMon_e, spikeMon_i)

################################################################################
#                              Main cycle
################################################################################

output_fname = "{0}/thesis_rasters.h5".format(options.output_dir)
ds = DataStorage.open(output_fname, 'w')
ds['trials'] = []

for trial_it in range(ei_net.no.ntrials):
    print "Starting trial no. " + str(trial_it) + "..."
    print "Simulation running..."
    start_time=time.time()

    print "  Network initialisation..."
    ei_net.net.run(options.theta_start_mon_t*msecond, report='stdout')

    print "  Theta stimulation..."
    ei_net.net.run((options.time - options.theta_start_mon_t)*msecond, report='stdout')
    duration=time.time()-start_time
    print "Simulation time:",duration,"seconds"

    senders_e, times_e = zip(*spikeMon_e.spikes)
    senders_i, times_i = zip(*spikeMon_i.spikes)
    ds['trials'].append(
        {
            'spikeMon_e' : {
                'events' : {
                    'senders' : np.asanyarray(senders_e),
                    'times': np.asanyarray(times_e),
                }
            },
            'spikeMon_i': {
                'events' : {
                    'senders' : np.asanyarray(senders_i),
                    'times': np.asanyarray(times_i),
                }
            },

            'net_attr': {
                'net_Ne': ei_net.net_Ne,
                'net_Ni': ei_net.net_Ni,
            },
        }
    )
#                            End main cycle
################################################################################

ds.close()

total_time = time.time()-total_start_t
print "Overall time: ", total_time, " seconds"

