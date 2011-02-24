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

from mult_bump_network import *
from mult_bump_simulation_mod import *



(options, args) = getOptParser().parse_args()
print "Options:"
print options

# Clock definitions
sim_dt = options.sim_dt*ms
vel_dt = 0.02*second
simulationClock = Clock(dt=sim_dt)
SNClock = Clock(dt=10*ms)
velocityClock = Clock(dt=vel_dt)
printStatusClock = Clock(dt=options.update_interval*second)

# Provisional loading of connection weights and initial conditions from a mat file (settings)
connFileName = '../../central_data_store/data/connections_96.mat'
initCondFileName = '../../data/initial_conditions/initial_conditions_96.mat'
vm_init_rand = 10.0*mV




################################################################################
#                              Network operations
################################################################################
@network_operation(velocityClock)
def updateVelocity():
    updateVelocityZero()


@network_operation(printStatusClock)
def printStatus():
    print "Simulated " + str(printStatusClock.t) + " seconds."
    sys.stdout.flush()
################################################################################





################################################################################
#                              Network setup
################################################################################
print "Starting network and connections initialization..."
start_time=time.time()

sheet_size  = options.sheet_size  # Total no. of neurons will be sheet_size^2
W           = loadConnectionsFromMat(connFileName)
initCond    = loadInitialConditions(initCondFileName)


[sheetGroup, inhibConn] = createNetwork(options, simulationClock, W)
sheetGroup.vm = initCond['membrane_potentials']


duration=time.time()-start_time
print "Network setup time:",duration,"seconds"
#                            End Network setup
################################################################################




# Print connections if necessary and exit
if (options.save_conn):
    saveConnectionsToMat("results/connection_matrix.mat", inhibConn, options)
    exit();

ratData = loadmat("../../data/hafting_et_al_2005/Hafting_Fig2c_Trial1_preprocessed.mat")
#print ratData['pos_timeStamps']

# Velocity inputs - for now zero velocity
input = options.input
sheetGroup.B = linspace(input*namp, input*namp, sheet_size**2)
vIndex = 0  # Bad habit, but there are no static variables in python,
rat_pos_x = ratData['pos_x']
rat_pos_y = ratData['pos_y']



# Record the number of spikes
if options.record_sn_row == True:
    #SNList = range(sheet_size**2 / 2, sheet_size**2 / 2 + sheet_size)
    SNList = range(0, sheet_size**2, sheet_size+1)
else:
    #SNList = [sheet_size**2/4, sheet_size**2/2, (sheet_size**2)*3/4]
    SNList = True
SNMonitor = StateMonitor(sheetGroup, 'vm', record = SNList,
        clock=SNClock)
SNgMonitor = StateMonitor(sheetGroup, 'gi', record = SNList,
        clock=SNClock)

spikeGroupStart = int((sheet_size**2)*0.45)
spikeGroupEnd = int(spikeGroupStart + (sheet_size**2)/10)
spikeMonitorG = sheetGroup[spikeGroupStart:spikeGroupEnd]
#spikeMonitor = SpikeMonitor(spikeMonitorG)
spikeMonitor = SpikeMonitor(sheetGroup)

#printConn(sheet_size, inhibConn, options.write_data, options.print_only_conn)

print "Simulation running..."
start_time=time.time()
net = Network(sheetGroup, inhibConn, printStatus, updateVelocity, spikeMonitor)
if options.record_sn == True:
    net.add(SNMonitor)
    net.add(SNgMonitor)

net.run(options.time*second)
duration=time.time()-start_time
print "Simulation time:",duration,"seconds"


# Directory and filenames constants
timeSnapshot = datetime.now().strftime("%Y-%m-%dT%H-%M-%S")
dirName = options.output_dir
population_fname = dirName + timeSnapshot + '_spacePlot.eps'
options_fname = dirName + timeSnapshot + '.params'
count_fname = dirName + timeSnapshot + '_linePlot.eps'
conn_fname = dirName + timeSnapshot + "_conn.eps"

output_fname = dirName
if options.job_num != -1:
    output_fname = output_fname + 'job' + str(options.job_num)
output_fname +=  '_' + timeSnapshot + '_output.mat'



# write recorded data into .mat file
if (options.write_data):
    #f = open(options_fname, 'w')
    #f.write(str(options))
    #f.close()
    saveResultsToMat(output_fname, options, ratData, spikeMonitor, SNMonitor,
            SNgMonitor)

