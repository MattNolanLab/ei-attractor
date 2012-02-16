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

from EI_network import *
from EI_network_sim_mod import *


(options, args) = getOptParser().parse_args()
print "Options:"
print options

# Clock definitions
sim_dt = options.sim_dt*second
simulationClock = Clock(dt=sim_dt)
stimClock = Clock(50*msecond)


################################################################################
#                              Network setup
################################################################################
print "Starting network and connections initialization..."
start_time=time.time()

ei_net = EI_Network(options, simulationClock)


duration=time.time()-start_time
print "Network setup time:",duration,"seconds"
#                            End Network setup
################################################################################

stim_start = int(0.4*ei_net.o.Ne)
stim_range = int(0.2*ei_net.o.Ne)
stim_current = 0*pA

@network_operation(stimClock)
def stimulateSubPopulation():
    if simulationClock.t > 500*msecond and simulationClock.t < 550*msecond:
        ei_net.E_pop.Iext[stim_start:stim_start+stim_range] = linspace(stim_current, stim_current, stim_range)
        print "Stimulation..."
    else:
        ei_net.E_pop.Iext = [ei_net.E_pop.Iext[0]] * len(ei_net.E_pop)


state_record = [0]

spikeMon_e = SpikeMonitor(ei_net.E_pop)
spikeMon_i = SpikeMonitor(ei_net.I_pop)
stateMon_e = StateMonitor(ei_net.E_pop, 'vm', record = state_record, clock=simulationClock)
stateMon_i = StateMonitor(ei_net.I_pop, 'vm', record = state_record, clock=simulationClock)
stateMon_Isyn_e = StateMonitor(ei_net.E_pop, 'Isyn', record = state_record, clock=simulationClock)
stateMon_Isyn_i = StateMonitor(ei_net.I_pop, 'Isyn', record = state_record, clock=simulationClock)
stateMon_g_ad_e = StateMonitor(ei_net.E_pop, 'g_ad', record = state_record, clock=simulationClock)
stateMon_g_ad_i = StateMonitor(ei_net.I_pop, 'g_ad', record = state_record, clock=simulationClock)

ei_net.net.add(spikeMon_e, spikeMon_i, stateMon_e, stateMon_i, stateMon_Isyn_e,
        stateMon_Isyn_i, stateMon_g_ad_e, stateMon_g_ad_i)
ei_net.net.add(stimulateSubPopulation)

print "Simulation running..."
start_time=time.time()

ei_net.net.run(options.time*second)
duration=time.time()-start_time
print "Simulation time:",duration,"seconds"

figSize = (12,8)

#figure(figsize=figSize)
#stateMon_e.plot()
#figure()
#stateMon_Isyn_e.plot()
#figure()
#stateMon_g_ad_e.plot()
#figure()
#stateMon_i.plot()
#figure()
#stateMon_Isyn_i.plot()
#figure()
#stateMon_g_ad_i.plot()

f = figure(figsize=figSize)
subplot2grid((5,1), (0, 0), rowspan=4)
raster_plot(spikeMon_e)
title('Network size: ' + str(ei_net.o.Ne) + '(E), ' + str(ei_net.o.Ni) + '(I)')
ylabel('Neuron no. (E)')
xlim((0, ei_net.o.time*1000))
xlabel('')
subplot2grid((5,1), (4, 0), rowspan=1)
raster_plot(spikeMon_i)
ylabel('Neuron no. (I)')
xlim((0, ei_net.o.time*1000))

timeSnapshot = datetime.now().strftime("%Y-%m-%dT%H-%M-%S")
dirName = options.output_dir

output_fname = dirName
if options.job_num != -1:
    output_fname = output_fname + '/job' + str(options.job_num)
output_fname += '_Ne_' + str(ei_net.o.Ne)
output_fname +=  '_' + timeSnapshot + '_raster_plot.pdf'

f.savefig(output_fname)

