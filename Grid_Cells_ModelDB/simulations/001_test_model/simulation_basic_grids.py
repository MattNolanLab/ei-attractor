#
#   simulation_basic_grids.py
#
#   Main simulation run: grid fields with theta input and all the inhibition
#   (for gamma) and place input.
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

from scipy.io   import loadmat
from scipy.io   import savemat
from optparse   import OptionParser

from parameters              import *
from grid_cell_network_nest  import *
from spike_analysis          import firingRateFromPairs

import time
import numpy as np
import nest
import nest.raster_plot

mV = 1e-3


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

ei_net = NestGridCellNetwork(options, simulationOpts=None)

#ei_net.uniformInhibition()
#ei_net.setConstantCurrent()
#ei_net.setStartCurrent()
#ei_net.setThetaCurrentStimulation()
#ei_net.setPlaceCurrentInput()

#const_v = [0.0, 1.0]
#ei_net.setConstantVelocityCurrent_e(const_v)
#ei_net.setVelocityCurrentInput_e()

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

spikeMon_e = ei_net.getSpikeDetector("E")
spikeMon_i = ei_net.getSpikeDetector("I")

stateMon_params = {
        'withtime' : True,
        'interval' : 0.1,
        'record_from' : ['V_m', 'I_clamp_AMPA', 'I_clamp_NMDA',
            'I_clamp_GABA_A', 'I_stim']
}
stateMon_e = ei_net.getStateMonitor("E", state_record_e, stateMon_params)
stateMon_i = ei_net.getStateMonitor("I", state_record_i, stateMon_params)



#x_lim = [options.time-0.5, options.time]
x_lim = [0, options.time]



################################################################################
#                              Main cycle
################################################################################
print "Simulation running..."
start_time=time.time()
    
ei_net.simulate(options.time)
duration=time.time()-start_time
print "Simulation time:",duration,"seconds"

output_fname = "{0}/{1}job{2:04}_".format(options.output_dir, options.fileNamePrefix, options.job_num)



# Raster plot
nest.raster_plot.from_device(spikeMon_e, hist=False)
nest.raster_plot.from_device(spikeMon_i, hist=False)

events_e = nest.GetStatus(stateMon_e)[0]['events']
events_i = nest.GetStatus(stateMon_i)[0]['events']

# External currents
figure()
ax = subplot(211)
plot(events_e['times'], events_e['I_stim'])
ylabel('E cell $I_{stim}$')
subplot(212)
plot(events_i['times'], events_i['I_stim'])
ylabel('I cell $I_{stim}$')
xlabel('Time (ms)')


# E/I Vm
figure()
ax = subplot(211)
plot(events_e['times'], events_e['V_m'])
ylabel('E cell $V_m$')
subplot(212)
plot(events_i['times'], events_i['V_m'])
ylabel('I cell $V_m$')
xlabel('Time (ms)')



figure()
ax = subplot(211)
plot(events_e['times'], events_e['I_clamp_GABA_A'])
ylabel('E synaptic current (pA)')
subplot(212, sharex=ax)
plot(events_i['times'], events_i['I_clamp_AMPA'] + events_i['I_clamp_NMDA'])
xlabel('Time (s)')
ylabel('I synaptic current (pA)')
xlim(x_lim)
savefig(output_fname + '_Isyn.pdf')

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
#


# Firing rate of E cells on the twisted torus
F_tstart = ei_net.no.time - 1e3
F_tend = ei_net.no.time

figure()
senders = nest.GetStatus(spikeMon_e)[0]['events']['senders'] - ei_net.E_pop[0]
spikeTimes   = nest.GetStatus(spikeMon_e)[0]['events']['times']
Fe = firingRateFromPairs(ei_net.net_Ne, senders, spikeTimes, F_tstart, F_tend)
pcolormesh(np.reshape(Fe*1e3, (ei_net.Ne_y, ei_net.Ne_x)))
xlabel('E neuron no.')
ylabel('E neuron no.')
colorbar()
axis('equal')
title('Firing rates (torus) of E cells')
savefig(output_fname + '_firing_snapshot_e.png')


## Firing rate of I cells on the twisted torus
figure()
senders = nest.GetStatus(spikeMon_i)[0]['events']['senders'] - ei_net.I_pop[0]
spikeTimes   = nest.GetStatus(spikeMon_i)[0]['events']['times']
Fi = firingRateFromPairs(ei_net.net_Ni, senders, spikeTimes, F_tstart, F_tend)
pcolormesh(np.reshape(Fi*1e3, (ei_net.Ni_y, ei_net.Ni_x)))
xlabel('I neuron no.')
ylabel('I neuron no.')
colorbar()
axis('equal')
title('Firing rates (torus) of I cells')
savefig(output_fname + '_firing_snapshot_i.png')


show()

        
#                            End main cycle
################################################################################

total_time = time.time()-total_start_t
print "Overall time: ", total_time, " seconds"

