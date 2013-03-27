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
import numpy as np
import logging  as lg

from matplotlib.pyplot  import *

from numpy.random       import choice   # numpy >= 1.7.0
from optparse           import OptionParser

from models.parameters  import *
from models.gc_net_nest import BasicGridCellNetwork
from data_storage       import DataStorage

import time
import nest


lg.basicConfig(level=lg.DEBUG)

parser          = getOptParser()
parser.add_option("--gammaNSample",   type="float",   help="Fraction of neurons in the network to sample from, for the frequency analysis.")

(options, args) = parser.parse_args()
options         = setOptionDictionary(parser, options)


################################################################################
#                              Network setup
################################################################################
print "Starting network and connections initialization..."
start_time=time.time()
total_start_t = time.time()

ei_net = BasicGridCellNetwork(options, simulationOpts=None)

const_v = [00.0, 0.0]
ei_net.setConstantVelocityCurrent_e(const_v)
#ei_net.setVelocityCurrentInput_e()


duration=time.time()-start_time
print "Network setup time:",duration,"seconds"
#                            End Network setup
################################################################################

NSample = int(options.gammaNSample * ei_net.net_Ne)
stateRecF_e = choice(ei_net.E_pop, NSample, replace=False)

stateMonF_params = {
        'withtime' : True,
        'interval' : 0.1,
        'record_from' : ['I_clamp_GABA_A']
}
stateMonF_e = ei_net.getGenericStateMonitor(stateRecF_e, stateMonF_params)


spikeMon_e, spikeMon_i, pc_spikeMon, stateMon_e, stateMon_i = ei_net.getMonitors()


################################################################################
#                              Main cycle
################################################################################
print "Simulation running..."
start_time=time.time()
    
ei_net.simulate(options.time, printTime=True)
duration=time.time()-start_time
print "Simulation time:",duration,"seconds"


output_fname = "{0}/{1}job{2:05}_output".format(options.output_dir, options.fileNamePrefix, options.job_num)
out = DataStorage.open(output_fname + ".h5", 'w')


################################################################################
#                               Save data
out['options'] = options._einet_optdict
out['ei_net']  = ei_net.getAttrDictionary()

# Save spikes
out['spikeMon_e'] = {
        'status'    : nest.GetStatus(spikeMon_e)[0],
        'gid_start' : ei_net.E_pop[0]
    }
out['spikeMon_i'] = {
        'status'    : nest.GetStatus(spikeMon_i)[0],
        'gid_start' : ei_net.I_pop[0]
    }

if (len(ei_net.PC) != 0):
    out['pc_spikeMon'] = {
            'status'    : nest.GetStatus(pc_spikeMon)[0],
            'gid_start' : ei_net.PC[0]
        }

#Save state variables
out['stateMon_e'] = nest.GetStatus(stateMon_e)
out['stateMon_i'] = nest.GetStatus(stateMon_i)
out['stateMonF_e']  = nest.GetStatus(stateMonF_e)

out.close()
#                             End Save data
################################################################################




#                            End main cycle
################################################################################

total_time = time.time()-total_start_t
print "Overall time: ", total_time, " seconds"

