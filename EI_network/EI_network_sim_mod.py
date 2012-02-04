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



def getOptParser():
    # Network parameters definition
    optParser = OptionParser()
    optParser.add_option("--Ne", type="int", default=1000, dest="Ne",
            help="Number of excitatory neurons")
    optParser.add_option("--Ni", type="int", default=100, dest="Ni",
            help="Number of inhibitory neurons")

    optParser.add_option("--AMPA_density", type="float", default=0.4,
            dest="AMPA_density", help="Density of E-->I connections")
    optParser.add_option("--GABA_density", type="float", default=0.8,
            dest="i_density", help="Density of I-->E connections")

    optParser.add_option("--taum_e", type="float", default=9.3, dest="taum_e",
            help="Mean of excitatory neuron membrane time constant (ms)")
    optParser.add_option("--EL_e", type="float", default=-68.5, dest="EL_e",
            help="Mean resting membrane potential of excitatory neurons (mV)")
    optParser.add_option("--Vt_e", type="float", default=-50, dest="Vt_e",
            help="Mean of the excitatory integrate and fire spiking threshold (mV)")
    optParser.add_option("--Vr_e", type="float", default=-50, dest="Vr_e",
            help="Mean of the excitatory integrate and fire reset potential (mV)")
    optParser.add_option("--Rm_e", type=float, default=44, dest="Rm_e",
            help="Mean of excitatory membrane resistance (MOhm)")
    optParser.add_option("--ad_tau_e_mean", type=float, default=40,
            dest="ad_tau_e_mean", help="Mean of excitatory adaptation time constant (msec)")
    optParser.add_option("--ad_tau_e_std", type=float, default=5,
            dest="ad_tau_e_std", help="Std. deviation of excitatory adaptation time constant (msec)")
    optParser.add_option("--ad_e_g_inc", type=float, default=1.136e-8,
            dest="ad_e_g_inc", help="After-spike excitatory increase of leak conductance (Siemens)")
    optParser.add_option("--deltaT_e", type=float, default=3,
            dest="deltaT_e", help="Sharpness of exponential I&F neuron spike
            initiation (excitatory, mV)")

    optParser.add_option("--taum_i", type="float", default=10, dest="taum_i",
            help="Mean of inhibitory neuron membrane time constant (ms)")
    optParser.add_option("--EL_i", type="float", default=-60, dest="EL_i",
            help="Mean resting membrane potential of inhibitory neurons (mV)")
    optParser.add_option("--Vt_i", type="float", default=-50, dest="Vt_i",
            help="Mean of the excitatory integrate and fire spiking threshold (mV)")
    optParser.add_option("--Vr_i", type="float", default=-50, dest="Vr_i",
            help="Mean of the inhibitory integrate and fire reset potential (mV)")
    optParser.add_option("--Rm_i", type=float, default=44, dest="Rm_i",
            help="Mean of inhibitory membrane resistance (MOhm)")
    optParser.add_option("--ad_tau_i_mean", type=float, default=7.5,
            dest="ad_tau_i_mean", help="Mean of inhibitory adaptation time constant (msec)")
    optParser.add_option("--ad_tau_i_std", type=float, default=0.5,
            dest="ad_tau_i_std", help="Std. deviation of inhibitory adaptation time constant (msec)")
    optParser.add_option("--ad_i_g_inc", type=float, default=1.136e-8,
            dest="ad_i_g_inc", help="After-spike inhibitory increase of leak
            conductance (S)")
    optParser.add_option("--deltaT_i", type=float, default=3,
            dest="deltaT_i", help="Sharpness of exponential I&F neuron spike
            initiation (inhibitory, mV)")

    optParser.add_option("--tau_AMPA", type="float", default=1e-3,
            dest="tau_AMPA", help="Mean of AMPA synaptic conductance time constant (ms)")
    optParser.add_option("--g_AMPA_mean", type="float", dest="g_AMPA_mean",
            help="Mean of AMPA connection synaptic conductance (nS)")
    optParser.add_option("--g_AMPA_std", type="float", dest="g_AMPA_std",
            help="Std. deviation of AMPA connections synaptic conductance (nS)")
    optParser.add_option("--tau_GABA_rise", type="float", default=1e-3,
            dest="tau_GABA_rise", help="Mean of GABA rising time constant (ms)")
    optParser.add_option("--tau_GABA_fall", type="float", default=1e-3,
            dest="tau_GABA_fall", help="Mean of GABA fall time constant (ms)")
    optParser.add_option("--g_GABA_mean", type="float", dest="g_GABA_mean",
            help="Mean of GABA connections synaptic conductance (nS)")

    optParser.add_option("--Vrev_AMPA", type="float", default=0, dest="Vrev_AMPA",
            help="AMPA reversal potential (mV)")
    optParser.add_option("--Vrev_GABA", type="float", default=-75,
            dest="Vrev_GABA", help="GABA reversal potential (mV)")

    optParser.add_option("--noise_sigma", type="float", default=0.02,
            dest="noise_sigma", help="Std. dev of neural noise (mV)")
    optParser.add_option("--sigma_init_cond", type="float", default=10,
            dest="sigma_init_cond", help="Std. dev distribution of initial
            membrane voltages (mV)")

    optParser.add_option("--refrac_abs", type="float", default=1,
            dest="refrac_abs", help="Absolute refractory period (msec)")

    optParser.add_option("-t", "--time", type="float", default=5.0, dest="time",
            help="Total simulation time [seconds]")
    optParser.add_option("--sim_dt", type=float, default=0.1, dest="sim_dt",
            help="Simulation time step (ms)")
    optParse.add_options("--spike_detect_th", type=float, default=40,
            dest="spike_detect_th", help="Spike detection threshold during
            numerical simulation (mV)")


    optParser.add_option("--output-dir", type="string", default="results/", dest="output_dir", help="Output directory path.")
    optParser.add_option("-u", '--update-interval', type="float", default=5.0,
            dest="update_interval", help="Duration between simulation status printouts")
    optParser.add_option("-n", "--job-num", type="int", default=-1, dest="job_num",
            help="Use argument of this option to specify the output file name number, instead of using time")

    return optParser


## Save the results of the simulation in a .mat file
#def saveResultsToMat(options, ratData, spikeMonitor, SNMonitor,
#        SNgMonitor, SNList, spikeMonList):
#    # SNMonitor - monitor of single neurons, defined by list
#    # SNList - list of neuron numbers monitored (SN response)
#    # spikeMonitor - monitor of spikes of all neurons
#    # rateMonitor - rate monitor of several neurons (TODO)
#
#    # Directory and filenames
#    timeSnapshot = datetime.now().strftime("%Y-%m-%dT%H-%M-%S")
#    dirName = options.output_dir
#    
#    output_fname = dirName
#    if options.job_num != -1:
#        output_fname = output_fname + 'job' + str(options.job_num)
#    output_fname +=  '_' + timeSnapshot + '_output.mat'
#
#    # Start saving everything
#    if ratData == None:
#        outData = {}
#    else:
#        outData = ratData
#
#    spikeCell = empty(len(spikeMonList), dtype=object);
#
#    if spikeMonitor != None:
#        for k,v in spikeMonitor.spiketimes.iteritems():
#            spikeCell[k] = v
#        outData['spikeCell'] = spikeCell
#
#    if SNMonitor != None:
#        outData['SNMonitor_values'] = SNMonitor.values_
#        outData['SNMonitor_times'] =  SNMonitor.times_
#        outData['SNList'] = SNList
#
#    if SNgMonitor != None:
#        outData['SNgMonitor_times'] = SNgMonitor.times_
#        outData['SNgMonitor_values'] =SNgMonitor.values_ 
#        outData['SNList'] = SNList
#
#    if SNg_adMonitor != None:
#        outData['SNg_adMonitor_times']  = SNg_adMonitor.times_
#        outData['SNg_adMonitor_values'] = SNg_adMonitor.values_
#        outData['SNList'] = SNList
#
#    # Merge the rat position and timestamp data into the output so it is self
#    # contained
#    #for k,v in ratData.iteritems():
#    #    outData['ratData_' + k] = v
#
#    # Save options as string - couldn't find any other way how to convert it to
#    # variables
#    outData['options'] = str(options)
#    outData['sheet_size'] = options.sheet_size
#
#    savemat(output_fname, outData, do_compression=True)
#
