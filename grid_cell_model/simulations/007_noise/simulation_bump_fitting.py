#
#   simulation_bump_fitting.py
#
#   Main simulation run: Fitting a Gaussian to the bump and frequency analysis.
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

from matplotlib.pyplot  import *

from numpy.random       import choice
from optparse           import OptionParser
from matplotlib.mlab    import detrend_mean, psd, window_hanning

from models.parameters  import *
from models.gc_net_nest import *
from data_storage       import DataStorage

import time
import nest
import nest.raster_plot

mV = 1e-3


lg.basicConfig(level=lg.DEBUG)

parser          = getOptParser()
parser.add_option("--ndumps",         type="int",    help="Number of data output dumps during the simulation")
parser.add_option("--gammaNSample",   type="float",   help="Fraction of neurons in the network to sample from, for the frequency analysis.")
parser.add_option("--gammaRangeLow",  type="float",  help="Relative gamma power analysis range low (Hz)")
parser.add_option("--gammaRangeHigh", type="float",  help="Relative gamma power analysis range high (Hz)")

(options, args) = parser.parse_args()
options         = setOptionDictionary(parser, options)

# Other
figSize = (12,8)
rcParams['font.size'] = 16


def getSenders(nest_spikeMon, id, gid_group_start):
    return nest.GetStatus(nest_spikeMon)[id]['events']['senders'] - \
            gid_group_start

def getSpikeTimes(nest_spikeMon, id):
    return nest.GetStatus(nest_spikesMon)[id]['events']['times']




################################################################################
#                              Network setup
################################################################################
print "Starting network and connections initialization..."
start_time=time.time()
total_start_t = time.time()

ei_net = NestGridCellNetwork(options, simulationOpts=None)

const_v = [00.0, 0.0]
ei_net.setConstantVelocityCurrent_e(const_v)
#ei_net.setVelocityCurrentInput_e()

#ei_net.setStartCurrent()
#ei_net.setPlaceCells()


duration=time.time()-start_time
print "Network setup time:",duration,"seconds"
#                            End Network setup
################################################################################

rec_all_spikes = True
if rec_all_spikes:
    nrecSpike_e = ei_net.Ne_x*ei_net.Ne_y
    nrecSpike_i = ei_net.Ni_x*ei_net.Ni_y
else:
    nrecSpike_e = 200
    nrecSpike_i = 50

state_record_e = [ei_net.Ne_x/2 -1 , ei_net.Ne_y/2*ei_net.Ne_x + ei_net.Ne_x/2 - 1]
state_record_i = [ei_net.Ni_x/2 - 1, ei_net.Ni_y/2*ei_net.Ni_x + ei_net.Ni_x/2 - 1]

NSample = int(options.gammaNSample * ei_net.net_Ne)
stateRecF_e = choice(ei_net.E_pop, NSample, replace=False)

spikeMon_e = ei_net.getSpikeDetector("E")
spikeMon_i = ei_net.getSpikeDetector("I")
pc_spikemon = ei_net.getGenericSpikeDetector(ei_net.PC, "Place cells")

stateMon_params = {
        'withtime' : True,
        'interval' : 0.1,
        'record_from' : ['V_m', 'I_clamp_AMPA', 'I_clamp_NMDA',
            'I_clamp_GABA_A', 'I_stim']
}
stateMonF_params = {
        'withtime' : True,
        'interval' : 0.1,
        'record_from' : ['I_clamp_GABA_A']
}
stateMon_e  = ei_net.getStateMonitor("E", state_record_e, stateMon_params)
stateMon_i  = ei_net.getStateMonitor("I", state_record_i, stateMon_params)
stateMonF_e = ei_net.getGenericStateMonitor(stateRecF_e, stateMonF_params)



x_lim = [options.time-2e3, options.time]
#x_lim = [0, options.time]



################################################################################
#                              Main cycle
################################################################################
print "Simulation running..."
start_time=time.time()
    
ei_net.simulate(options.time, printTime=True)
duration=time.time()-start_time
print "Simulation time:",duration,"seconds"


output_fname = "{0}/{1}job{2:05}".format(options.output_dir, options.fileNamePrefix, options.job_num)
out = DataStorage.open(output_fname + ".h5")


################################################################################
#                               Save data
# Save spikes
out['spikeMon_e'] = {
        'status'    = nest.GetStatus(spikeMon_e)[0]
        'gid_start' = ei_net.E_pop[0]
    }
out['spikeMon_i'] = {
        'status'    = nest.GetStatus(spikeMon_i)[0]
        'gid_start' = ei_net.I_pop[0]
    }

if (len(ei_net.PC) != 0):
    out['pc_spikeMon'] = {
            'status'    = nest.GetStatus(pc_spikeMon)[0]
            'gid_start' = ei_net.PC[0]
        }

#Save state variables
out['stateMon_i_0'] = nest.GetStatus(stateMon_i)[0]
out['stateMon_i_1'] = nest.GetStatus(stateMon_i)[1]

out.close()
#                             End Save data
################################################################################




#                            End main cycle
################################################################################

total_time = time.time()-total_start_t
print "Overall time: ", total_time, " seconds"

