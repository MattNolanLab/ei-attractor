#import brian_no_units
from brian import *

from brian import *
from brian.library.IF import *
from brian.library.synapses import *

from scipy import linspace
from scipy.io import loadmat
from scipy.io import savemat
from optparse import OptionParser
from datetime import datetime

import time
import math
import sys
import numpy as np
import logging as lg

from EI_network import *
from EI_network_sim_mod import *
from custombrian import *

lg.basicConfig(level=lg.DEBUG)


(options, args) = getOptParser().parse_args()
print "Options:"
print options

# Clock definitions
sim_dt = options.sim_dt*second
simulationClock = Clock(dt=sim_dt)
stimClock = Clock(50*msecond)
rateClock = Clock(500*msecond)


################################################################################
#                              Network setup
################################################################################
print "Starting network and connections initialization..."
start_time=time.time()

options.ndim = 2
ei_net = EI_Network(options, simulationClock)

# Mexican hat properties and AMPA/GABA connections
pAMPA_mu = 0.5
pAMPA_sigma = 0.25/6
pGABA_sigma = 0.5/6
ei_net.connMexicanHat(pAMPA_mu, pAMPA_sigma, pGABA_sigma)


duration=time.time()-start_time
print "Network setup time:",duration,"seconds"
#                            End Network setup
################################################################################

stim_start = int(0.4*ei_net.o.Ne)
stim_range = int(0.2*ei_net.o.Ne)
stim_current = 1200*pA

@network_operation(stimClock)
def stimulateSubPopulation():
    if simulationClock.t > 500*msecond and simulationClock.t < 650*msecond:
        tmp = ei_net.E_pop.Iext.reshape((options.Ne, options.Ne))
        tmp[stim_start:stim_start+stim_range, stim_start:stim_start+stim_range] =\
            linspace(stim_current, stim_current, stim_range**2).reshape((stim_range, stim_range))
        ei_net.E_pop.Iext = tmp.ravel()
        print "Stimulation..."
    else:
        ei_net.E_pop.Iext = [ei_net.E_pop.Iext[0]] * len(ei_net.E_pop)


state_record_e = [465]
state_record_i = [465]

spikeMon_e = ExtendedSpikeMonitor(ei_net.E_pop)
spikeMon_i = ExtendedSpikeMonitor(ei_net.I_pop)
stateMon_e = StateMonitor(ei_net.E_pop, 'vm', record = state_record_e, clock=simulationClock)
stateMon_i = StateMonitor(ei_net.I_pop, 'vm', record = state_record_i, clock=simulationClock)
stateMon_Isyn_e = StateMonitor(ei_net.E_pop, 'Isyn', record = state_record_e, clock=simulationClock)
stateMon_Isyn_i = StateMonitor(ei_net.I_pop, 'Isyn', record = state_record_i, clock=simulationClock)
stateMon_g_ad_e = StateMonitor(ei_net.E_pop, 'g_ad', record = state_record_e, clock=simulationClock)
stateMon_g_ad_i = StateMonitor(ei_net.I_pop, 'g_ad', record = state_record_i, clock=simulationClock)

ei_net.net.add(spikeMon_e, spikeMon_i, stateMon_e, stateMon_i, stateMon_Isyn_e,
        stateMon_Isyn_i, stateMon_g_ad_e, stateMon_g_ad_i)
ei_net.net.add(stimulateSubPopulation)

print "Simulation running..."
start_time=time.time()

ei_net.net.run(options.time*second)
duration=time.time()-start_time
print "Simulation time:",duration,"seconds"

figSize = (12,8)

figure(figsize=figSize)
subplot(211)
plot(stateMon_e.times, stateMon_e.values[0]*1e3)
xlabel('Time (s)')
ylabel('$V_m$ (mV)')
title('E membrane voltage')

subplot(212)
plot(stateMon_Isyn_e.times, stateMon_Isyn_e.values[0]*1e12)
xlabel('Time (s)')
ylabel('Current (pA)')
#figure()
#plot(stateMon_g_ad_e.times, stateMon_g_ad_e.values[0]*1e9)
#xlabel('Time (ms)')
#ylabel('Conductance (nS)')

figure()
subplot(211)
plot(stateMon_i.times, stateMon_i.values[0]*1e3)
xlabel('Time (s)')
ylabel('$V_m$ (mV)')
title('I membrane voltage')
subplot(212)
plot(stateMon_Isyn_i.times, stateMon_Isyn_i.values[0]*1e12)
xlabel('Time (s)')
ylabel('Current (pA)')


#f = figure(figsize=figSize)
#subplot2grid((5,1), (0, 0), rowspan=3)
#raster_plot(spikeMon_e)
#title('Network size: ' + str(ei_net.o.Ne) + '(E), ' + str(ei_net.o.Ni) + '(I)')
#ylabel('Neuron no. (E)')
#xlim((0, ei_net.o.time*1000))
#xlabel('')
#subplot2grid((5,1), (3, 0), rowspan=2)
#raster_plot(spikeMon_i)
#ylabel('Neuron no. (I)')
#xlim((0, ei_net.o.time*1000))


F_tstart = 0
F_tend = options.time
F_dt = 0.2
F_winLen = 1.
f = figure(figsize=figSize)
title('Network size: ' + str(ei_net.o.Ne) + '(E), ' + str(ei_net.o.Ni) + '(I)')
subplot2grid((5,1), (0, 0), rowspan=3)
Fe, Fe_t = spikeMon_e.getFiringRate(F_tstart, F_tend, F_dt, F_winLen) 
X, Y = meshgrid(Fe_t, np.arange(len(Fe)))
pcolor(X, Y, Fe)
colorbar()
ylabel('Neuron no. (E)')
xlabel('')
subplot2grid((5,1), (3, 0), rowspan=2)
Fi, Fi_t = spikeMon_i.getFiringRate(F_tstart, F_tend, F_dt, F_winLen)
X, Y = meshgrid(Fi_t, np.arange(len(Fi)))
pcolor(X, Y, Fi)
colorbar()
ylabel('Neuron no. (I)')
xlabel('Time (s)')


# Print a plot of bump position
(pos, times) = spikeMon_e.torusPopulationVector(ei_net.o.Ne, F_tstart, F_tend, F_dt,
        F_winLen)
figure(figsize=figSize)
plot(times, pos)
xlabel('Time (s)')
ylabel('Bump position (neurons)')


timeSnapshot = datetime.now().strftime("%Y-%m-%dT%H-%M-%S")
dirName = options.output_dir

output_fname = dirName
if options.job_num != -1:
    output_fname = output_fname + '/job' + str(options.job_num)
output_fname += '_Ne_' + str(ei_net.o.Ne)
output_fname +=  '_' + timeSnapshot + '_raster_plot.pdf'

#f.savefig(output_fname)
show()

