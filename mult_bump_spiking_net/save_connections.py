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

from mult_bump_network import *
from mult_bump_simulation_mod import *

# Perform the network setup, but only export connection weights data to
# a specified file

optParser = getOptParser()
(options, args) = optParser.parse_args()
print "Options:"
print options


simulationClock = Clock()


################################################################################
#                              Network setup
################################################################################
print "Starting network and connections initialization..."
start_time=time.time()

sheet_size  = options.sheet_size  # Total no. of neurons will be sheet_size^2
W           = None


[sheetGroup, inhibConn] = createNetwork(options, simulationClock, W)


duration=time.time()-start_time
print "Network setup time:",duration,"seconds"
#                            End Network setup
################################################################################


# Save connections if necessary and exit
saveConnectionsToMat("results/connection_matrix.mat", inhibConn, options)
