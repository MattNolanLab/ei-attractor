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
sim_dt = options.sim_dt*ms
simulationClock = Clock(dt=sim_dt)
#SNClock = Clock(dt=10*ms)




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

spikeMon_e = SpikeMonitor(ei_net.E_pop)
spikeMon_i = SpikeMonitor(ei_net.I_pop)
stateMon_e = StateMonitor(ei_net.E_pop, 'vm', record = [1, 2, 3])
stateMon_i = StateMonitor(ei_net.I_pop, 'vm', record = [1, 2, 3])
stateMon_Isyn_e = StateMonitor(ei_net.E_pop, 'Isyn', record = [1, 2, 3])
stateMon_Isyn_i = StateMonitor(ei_net.I_pop, 'Isyn', record = [1, 2, 3])

ei_net.net.add(spikeMon_e, spikeMon_i, stateMon_e, stateMon_i, stateMon_Isyn_e,
        stateMon_Isyn_i)

print "Simulation running..."
start_time=time.time()

ei_net.net.run(options.time*second)
duration=time.time()-start_time
print "Simulation time:",duration,"seconds"

figure()
raster_plot(spikeMon_e, spikeMon_i)

figure()
stateMon_e.plot()
figure()
stateMon_Isyn_e.plot()
figure()
stateMon_i.plot()
figure()
stateMon_Isyn_i.plot()

show()

#saveResultsToMat(options, None, spikeMonitor, None, None)
