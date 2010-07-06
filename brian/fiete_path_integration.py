import brian_no_units
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

from fiete_network import *


# Network parameters definition
optParser = OptionParser()
optParser.add_option("--lambda-net", type="float", default=13,
        dest="lambda_net")
optParser.add_option("-s", "--sheet-size", type="int", default=40, dest="sheet_size")
optParser.add_option("-l", type="float", default=2.0, dest="l")
optParser.add_option("-a", type="float", default=1.0, dest="a");
optParser.add_option("-c", "--conn-mult", type="float", default=10.0,
        dest="connMult")
optParser.add_option("-w", "--write-data", action="store_true",
        dest="write_data", default=False)
optParser.add_option("-i", "--input", type="float", default=0.3, dest="input",
        help="Feedforward constant input to all neurons")
optParser.add_option("-p", "--print-only-conn", action="store_true",
        dest="print_only_conn", default=False)
optParser.add_option("--alpha", type="float", default=0.10315, dest="alpha",
        help="Network alpha parameter (see Burak&Fiete, 2009)")
optParser.add_option("-t", "--time", type="float", default=5.0, dest="time",
        help="Total simulation time [seconds]")
optParser.add_option("-u", '--update-interval', type=float, default=5.0,
        dest="update_interval", help="Duration between simulation status printouts");

(options, args) = optParser.parse_args()
print "Options:"
print options

# Clock definitions
sim_dt = 0.1*ms
vel_dt = 20*ms
simulationClock = Clock(dt=sim_dt)
#SNClock = Clock(dt=10*sim_dt)
velocityClock = Clock(dt=vel_dt)
printStatusClock = Clock(dt=options.update_interval*second)

# Time for the network to settle down into an attractor state at the beginning
network_setup_time = 5*second


# Save the results of the simulation in a .mat file
def saveResultsToMat(fileName, spikeMonitor, ratData, options):
    # SNMonitor - monitor of single neurons, defined by list
    # SNList - list of neuron numbers monitored (SN response)
    # spikeMonitor - monitor of spikes of all neurons
    # rateMonitor - rate monitor of several neurons (TODO)
    outData = dict()

    for k,v in spikeMonitor.spiketimes.iteritems():
        outData['spikeMonitor_times_n' + str(k)] = v

    #outData['SNMonitor_values'] = SNMonitor.values_
    #outData['SNMonitor_times'] = SNMonitor.times_
    #outData['SNList'] = SNList

    # Merge the rat position and timestamp data into the output so it is self
    # contained
    for k,v in ratData.iteritems():
        outData['ratData_' + k] = v

    # Save options as string - couldn't find any other way how to convert it to
    # variables
    outData['options'] = str(options)
    outData['sheet_size'] = options.sheet_size

    savemat(fileName, outData, do_compression=True)


# Definition of 2d topography on a sheet and connections
sheet_size = options.sheet_size  # Total no. of neurons will be sheet_size^2

# Definition of connection parameters and setting the connection iteratively:
# Mexican hat connectivity

start_time=time.time()

[sheetGroup, inhibConn] = createNetwork(sheet_size, options.lambda_net,
        options.l, options.a, options.connMult, simulationClock)

duration=time.time()-start_time
print "Connection setup time:",duration,"seconds"

ratData = loadmat("../../data/hafting_et_al_2005/Hafting_Fig2c_Trial2.mat")
#print ratData['pos_timeStamps']

# Velocity inputs - for now zero velocity
input = options.input
sheetGroup.B = linspace(input*namp, input*namp, sheet_size**2)
vIndex = 0  # Bad habit, but there are no static variables in python,
rat_pos_x = ratData['pos_x'][0]
rat_pos_y = ratData['pos_y'][0]#
@network_operation(velocityClock)
def updateVelocity():
    #updateVelocityLinear()
    updateVelocityRat()


def updateVelocityRat():
    global vIndex
    global input
    # the original data are in cm/s, however, we rather want m/s 
    vel_x = (rat_pos_x[vIndex + 1] - rat_pos_x[vIndex])/vel_dt/100
    vel_y = (rat_pos_y[vIndex + 1] - rat_pos_y[vIndex])/vel_dt/100
    #vel_x = 0
    #vel_y = 0

    i = 0
    for i_y in range(sheet_size):
        for i_x in range(sheet_size):
            prefDir = getPreferredDirection(i_x, i_y)
            sheetGroup.B[i] = (input + options.alpha*(prefDir[0]*vel_x +
                prefDir[1]*vel_y))*namp
            i+=1

    vIndex+=1

linSize = 0.9   # m
ratSpeed = 0.05 # m/s
dts_side = linSize/ratSpeed/vel_dt # Number of vel_dts to get from one side to
                                   # the other
curr_pos_x = 0                                     

def updateVelocityLinear():
    global vIndex
    global input
    global curr_pos_x

    dir = (vIndex // dts_side) % 2
    if dir == 0:
        vel_x = ratSpeed
    else:
        vel_x = -ratSpeed

    vel_y = 0

    i = 0
    for i_y in range(sheet_size):
        for i_x in range(sheet_size):
            prefDir = getPreferredDirection(i_x, i_y)
            sheetGroup.B[i] = (input + options.alpha*(prefDir[0]*vel_x +
                prefDir[1]*vel_y))*namp
            i+=1

    ratData['pos_x'][0][vIndex] = curr_pos_x*100 # m --> cm
    ratData['pos_y'][0][vIndex] = 0
    curr_pos_x += vel_x*vel_dt
    vIndex+=1


@network_operation(printStatusClock)
def printStatus():
    print "Simulated " + str(printStatusClock.t) + " seconds."

# Record the number of spikes
#SNList = [sheet_size**2/2]
#Monitor = SpikeCounter(sheetGroup)
#SNMonitor = StateMonitor(sheetGroup, 'vm', record = SNList,
#        clock=SNClock)

#spikeGroupStart = (sheet_size**2)/4
#spikeGroupEnd = spikeGroupStart + (sheet_size**2)/2
#spikeMonitorG = sheetGroup[spikeGroupStart:spikeGroupEnd]
spikeMonitor = SpikeMonitor(sheetGroup)

#printConn(sheet_size, inhibConn, options.write_data, options.print_only_conn)

print "Simulation running..."
start_time=time.time()
#net = Network(sheetGroup, inhibConn, printStatus)
net = Network(sheetGroup, inhibConn, printStatus, updateVelocity, spikeMonitor)
net.run(options.time*second)
duration=time.time()-start_time
print "Simulation time:",duration,"seconds"


# Directory and filenames constants
timeSnapshot = datetime.now().strftime("%Y-%m-%dT%H-%M-%S")
dirName = 'results/'
population_fname = dirName + timeSnapshot + '_spacePlot.eps'
options_fname = dirName + timeSnapshot + '.params'
count_fname = dirName + timeSnapshot + '_linePlot.eps'
conn_fname = dirName + timeSnapshot + "_conn.eps"

output_fname = dirName + timeSnapshot + '_output.mat'



# write recorded data into .mat file
if (options.write_data):
    #f = open(options_fname, 'w')
    #f.write(str(options))
    #f.close()
    saveResultsToMat(output_fname, spikeMonitor, ratData, options)

