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

#from mult_bump_network import *


def getOptParser():
    # Network parameters definition
    optParser = OptionParser()
    optParser.add_option("--lambda-net", type="float", default=13,
            dest="lambda_net")
    optParser.add_option("-s", "--sheet-size", type="int", default=40, dest="sheet_size")
    optParser.add_option("-l", type="float", default=2.0, dest="l")
    optParser.add_option("-a", type="float", default=1.0, dest="a");
    optParser.add_option("-c", "--conn-mult", type="float", default=10.0,
            dest="connMult")
    optParser.add_option("-i", "--input", type="float", default=0.3, dest="input",
            help="Feedforward constant input to all neurons")
    optParser.add_option("--noise-sigma", type="float", default=0, dest="noise_sigma", help="Input current Gaussian noise sigma (mV)")
    optParser.add_option("--alpha", type="float", default=0.10315, dest="alpha",
            help="Network alpha parameter (see Burak&Fiete, 2009)")
    optParser.add_option("-t", "--time", type="float", default=5.0, dest="time",
            help="Total simulation time [seconds]")
    optParser.add_option("-u", '--update-interval', type="float", default=5.0,
            dest="update_interval", help="Duration between simulation status printouts")
    optParser.add_option("--sim_dt", type=float, default=0.1, dest="sim_dt",
            help="Simulation time step (ms)")
    optParser.add_option("-n", "--job-num", type="int", default=-1, dest="job_num",
            help="Use argument of this option to specify the output file name number, instead of using time")
    optParser.add_option("--taum", type="float", default=9.3, dest="taum",
            help="Neuron membrane time constant (ms)")
    optParser.add_option("--taui", type="float", default=10, dest="taui",
            help="Inhibitory synaptic time constant (ms)")
    optParser.add_option("--EL", type="float", default=-68.5, dest="EL",
            help="Resting membrane potential of neurons (mV)")
    optParser.add_option("--threshold", type="float", default=-55, dest="threshold",
            help="Integrate and fire spiking threshold (mV)")
    optParser.add_option("--Rm", type=float, default=44, dest="Rm",
            help="Membrane resistance (MOhm)")
    optParser.add_option("--tau_ad", type="float", default=40, dest="tau_ad",
            help="Integrate and fire adaptation time constant (msec)")
    optParser.add_option("--output-dir", type="string", default="results/", dest="output_dir", help="Output directory path.")

    optParser.add_option("--Ilinresp", type="float", default=0, dest="Ilinresp",
            help="Amount of current to inject in addition to 'input' in network response simulations (nA)")
    optParser.add_option("--initT", type="float", default=0, dest="initT",
            help="Initialisation time, during which no movement of rat is forced (seconds)")

    return optParser


# Save the results of the simulation in a .mat file
def saveResultsToMat(options, ratData, spikeMonitor, SNList, SNMonitor,
        SNgMonitor, SNg_adMonitor):
    # SNMonitor - monitor of single neurons, defined by list
    # SNList - list of neuron numbers monitored (SN response)
    # spikeMonitor - monitor of spikes of all neurons
    # rateMonitor - rate monitor of several neurons (TODO)

    # Directory and filenames
    timeSnapshot = datetime.now().strftime("%Y-%m-%dT%H-%M-%S")
    dirName = options.output_dir
    
    output_fname = dirName
    if options.job_num != -1:
        output_fname = output_fname + 'job' + str(options.job_num)
    output_fname +=  '_' + timeSnapshot + '_output.mat'

    # Start saving everything
    if ratData == None:
        outData = {}
    else:
        outData = ratData

    spikeCell = empty((options.sheet_size**2), dtype=object);

    if spikeMonitor != None:
        for k,v in spikeMonitor.spiketimes.iteritems():
            spikeCell[k] = v
        outData['spikeCell'] = spikeCell

    if SNMonitor != None:
        outData['SNMonitor_values'] = SNMonitor.values_
        outData['SNMonitor_times'] =  SNMonitor.times_
        outData['SNList'] = SNList

    if SNgMonitor != None:
        outData['SNgMonitor_times'] = SNgMonitor.times_
        outData['SNgMonitor_values'] =SNgMonitor.values_ 
        outData['SNList'] = SNList

    if SNg_adMonitor != None:
        outData['SNg_adMonitor_times']  = SNg_adMonitor.times_
        outData['SNg_adMonitor_values'] = SNg_adMonitor.values_
        outData['SNList'] = SNList

    # Merge the rat position and timestamp data into the output so it is self
    # contained
    #for k,v in ratData.iteritems():
    #    outData['ratData_' + k] = v

    # Save options as string - couldn't find any other way how to convert it to
    # variables
    outData['options'] = str(options)
    outData['sheet_size'] = options.sheet_size

    savemat(output_fname, outData, do_compression=True)

def saveConnectionsToMat(fileName, conn, options):
    outData = {};
    outData['connections'] = asarray(conn.W)
    outData['options'] = str(options)
    savemat(fileName, outData, do_compression=True)


def loadConnectionsFromMat(fileName):
    return loadmat(fileName)['connections']

def loadInitialConditions(fileName):
    return loadmat(fileName)



#def updateVelocityRat(vIndex, input, rat_pos_x, rat_pos_y, sheetGroup, vel_dt):
#    # the original data are in cm/s, however, we rather want m/s 
#    vel_x = (rat_pos_x[vIndex + 1] - rat_pos_x[vIndex])/vel_dt/100*second
#    vel_y = (rat_pos_y[vIndex + 1] - rat_pos_y[vIndex])/vel_dt/100*second
#
#    i = 0
#    for i_y in xrange(sheet_size):
#        for i_x in xrange(sheet_size):
#            prefDir = getPreferredDirection(i_x, i_y)
#            sheetGroup.B[i] = (input + options.alpha*(prefDir[0]*vel_x +
#                prefDir[1]*vel_y))*namp
#            i+=1


def updateVelocityZero():
    # Velocity is zero, do nothing
    nothing = 5

##linSize = 0.2   # m
#ratSpeed = 0.1 # m/s
##dts_side = linSize/ratSpeed/vel_dt # Number of vel_dts to get from one side to
#                                   # the other
#curr_pos_x = 0                                     
#
#
#
#def updateVelocityUp():
#    # Move the rat always to the right
#    vel_x = 0
#    vel_y = ratSpeed
#
#    i = 0
#    for i_y in range(sheet_size):
#        for i_x in range(sheet_size):
#            prefDir = getPreferredDirection(i_x, i_y)
#            sheetGroup.B[i] = (input + options.alpha*(prefDir[0]*vel_x +
#                prefDir[1]*vel_y))*namp
#            i+=1
#
#
#def updateVelocityLinear():
#    global vIndex
#    global input
#    global curr_pos_x
#
#    if vIndex == 0:
#        vel_x = 0
#        vel_y = ratSpeed
#    elif vIndex == 1:
#        vel_x = 0
#        vel_y = 0
#    elif vIndex == 2:
#        vel_x = ratSpeed/math.sqrt(2)
#        vel_y = ratSpeed/math.sqrt(2)
#    elif vIndex == 3:
#        vel_x = 0
#        vel_y = 0
#    elif vIndex == 4:
#        vel_x = ratSpeed
#        vel_y = 0
#    elif vIndex == 5:
#        vel_x = 0
#        vel_y = 0
#    elif vIndex == 6:
#        vel_x = ratSpeed/math.sqrt(2)
#        vel_y = -ratSpeed/math.sqrt(2)
#    elif vIndex == 7:
#        vel_x = 0
#        vel_y = 0
#    elif vIndex == 8:
#        vel_x = 0
#        vel_y = -ratSpeed
#    elif vIndex == 9:
#        vel_x = 0
#        vel_y = 0
#    elif vIndex == 10:
#        vel_x = -ratSpeed/math.sqrt(2)
#        vel_y = -ratSpeed/math.sqrt(2)
#    elif vIndex == 11:
#        vel_x = 0
#        vel_y = 0
#    elif vIndex == 12:
#        vel_x = -ratSpeed
#        vel_y = 0
#    elif vIndex == 13:
#        vel_x = 0
#        vel_y = 0
#    elif vIndex == 14:
#        vel_x = -ratSpeed/math.sqrt(2)
#        vel_y = ratSpeed/math.sqrt(2)
#    else: 
#        vel_x = 0
#        vel_y = 0
#
#    i = 0
#    for i_y in range(sheet_size):
#        for i_x in range(sheet_size):
#            prefDir = getPreferredDirection(i_x, i_y)
#            sheetGroup.B[i] = (input + options.alpha*(prefDir[0]*vel_x +
#                prefDir[1]*vel_y))*namp
#            i+=1
#
#    ratData['pos_x'][0][vIndex] = curr_pos_x*100 # m --> cm
#    ratData['pos_y'][0][vIndex] = 0
#    curr_pos_x += vel_x*vel_dt
#    vIndex+=1
#    vIndex = vIndex % 16


