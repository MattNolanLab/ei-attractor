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
SNClock = Clock(dt=10*ms)
printStatusClock = Clock(dt=options.update_interval*second)




################################################################################
#                              Network setup
################################################################################
print "Starting network and connections initialization..."
start_time=time.time()

net = EI_Network(options, simulationClock)


duration=time.time()-start_time
print "Network setup time:",duration,"seconds"
#                            End Network setup
################################################################################


## Velocity inputs - for now zero velocity
#input = options.input
#sheetGroup.B = linspace(input*namp, input*namp, sheet_size**2)
#
#spikeMonitor = SpikeMonitor(sheetGroup)
#
##printConn(sheet_size, inhibConn, options.write_data, options.print_only_conn)
#
#print "Simulation running..."
#start_time=time.time()
#net = Network(sheetGroup, inhibConn, printStatus, spikeMonitor)
#
#net.run(options.time*second)
#duration=time.time()-start_time
#print "Simulation time:",duration,"seconds"
#
#saveResultsToMat(options, None, spikeMonitor, None, None)
