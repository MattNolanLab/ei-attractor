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
    optParser.add_option("--Ne", type="int", dest="Ne",
            help="Number of excitatory neurons")
    optParser.add_option("--Ni", type="int", dest="Ni",
            help="Number of inhibitory neurons")
    optParser.add_option("--ntrials", type=int, default=1,
            help="Number of trials for the parameter set")

    optParser.add_option("--AMPA_density", type="float", 
            dest="AMPA_density", help="Density of E-->I connections")
    optParser.add_option("--GABA_density", type="float", 
            dest="GABA_density", help="Density of I-->E connections")

    optParser.add_option("--Iext_e", type="float",
            help="External background current into the excitatory population (A)")
    optParser.add_option("--Iext_i", type="float",
            help="External background current into the inhibitory population (A)")

    optParser.add_option("--taum_e", type="float", dest="taum_e",
            help="Mean of excitatory neuron membrane time constant (s)")
    optParser.add_option("--taum_e_spread", type="float",
            help="Spread of --taum_e (see --taum_i_spread option) (s)")
    optParser.add_option("--EL_e", type="float", dest="EL_e",
            help="Mean resting membrane potential of excitatory neurons (V)")
    optParser.add_option("--EL_e_spread", type="float",
            help="EL_e spread (see --taum_e_spread options for details) (V)")
    optParser.add_option("--Vt_e", type="float", dest="Vt_e",
            help="Mean of the excitatory integrate and fire spiking threshold (V)")
    optParser.add_option("--Vr_e", type="float", dest="Vr_e",
            help="Mean of the excitatory integrate and fire reset potential (V)")
    optParser.add_option("--Rm_e", type="float", dest="Rm_e",
            help="Mean of excitatory membrane resistance (Ohm)")
    optParser.add_option("--ad_tau_e_mean", type="float",
            dest="ad_tau_e_mean", help="Mean of excitatory adaptation time constant (sec)")
    optParser.add_option("--ad_tau_e_std", type="float",
            dest="ad_tau_e_std", help="Std. deviation of excitatory adaptation time constant (sec)")
    optParser.add_option("--ad_e_g_inc", type="float",
            dest="ad_e_g_inc", help="After-spike excitatory increase of leak conductance (Siemens)")
    optParser.add_option("--deltaT_e", type="float",
            dest="deltaT_e", help="Sharpness of exponential I&F neuron spike"
            "initiation (excitatory, V)")
    optParser.add_option("--Eahp_e", type="float",
            help="Excitatory AHP reversal potential (Volt)")
    optParser.add_option("--g_ahp_e", type="float",
            help="Maximal AHP excitatory conductance (Siemens)")
    optParser.add_option("--tau_ahp_e", type="float",
            help="Excitatory AHP time decay constant (seconds)")

    optParser.add_option("--taum_i", type="float", dest="taum_i",
            help="Mean of inhibitory neuron membrane time constant (s)")
    optParser.add_option("--taum_i_spread", type="float",
            help="Spread of --taum_i: drawn from uniform distrib. "
            "<taum_i-spread/2, taum_i+spread/2> (s)")
    optParser.add_option("--EL_i", type="float", dest="EL_i",
            help="Mean resting membrane potential of inhibitory neurons (V)")
    optParser.add_option("--EL_i_spread", type="float",
            help="ELs will be uniformly generated from <EL-spread, EL+spread> (V)")
    optParser.add_option("--Vt_i", type="float", dest="Vt_i",
            help="Mean of the excitatory integrate and fire spiking threshold (V)")
    optParser.add_option("--Vr_i", type="float", dest="Vr_i",
            help="Mean of the inhibitory integrate and fire reset potential (V)")
    optParser.add_option("--Rm_i", type="float", dest="Rm_i",
            help="Mean of inhibitory membrane resistance (Ohm)")
    optParser.add_option("--ad_tau_i_mean", type="float",
            dest="ad_tau_i_mean", help="Mean of inhibitory adaptation time constant (sec)")
    optParser.add_option("--ad_tau_i_std", type="float",
            dest="ad_tau_i_std", help="Std. deviation of inhibitory adaptation time constant (sec)")
    optParser.add_option("--ad_i_g_inc", type="float",
            dest="ad_i_g_inc", help="After-spike inhibitory increase of leak"
            "conductance (S)")
    optParser.add_option("--deltaT_i", type="float",
            dest="deltaT_i", help="Sharpness of exponential I&F neuron spike"
            "initiation (inhibitory, V)")

    optParser.add_option("--tau_AMPA", type="float",
            dest="tau_AMPA", help="Mean of AMPA synaptic conductance time constant (s)")
    optParser.add_option("--g_AMPA_total", type="float",
            help="Mean of total AMPA connection synaptic conductance (S)")
    optParser.add_option("--g_AMPA_std", type="float", dest="g_AMPA_std",
            help="Std. deviation of AMPA connections synaptic conductance (S)")
    optParser.add_option("--tau_GABA_rise", type="float",
            dest="tau_GABA_rise", help="Mean of GABA rising time constant (s)")
    optParser.add_option("--tau_GABA_fall", type="float",
            dest="tau_GABA_fall", help="Mean of GABA fall time constant (s)")
    optParser.add_option("--g_GABA_total", type="float",
            help="Mean of GABA connections synaptic conductance (S)")

    optParser.add_option("--Vrev_AMPA", type="float", dest="Vrev_AMPA",
            help="AMPA reversal potential (V)")
    optParser.add_option("--Vrev_GABA", type="float",
            dest="Vrev_GABA", help="GABA reversal potential (V)")

    optParser.add_option("--noise_sigma", type="float",
            dest="noise_sigma", help="Std. dev of neural noise (V)")
    optParser.add_option("--sigma_init_cond", type="float",
            dest="sigma_init_cond", help="Std. dev distribution of initial"
            "membrane voltages (V)")

    optParser.add_option("--refrac_abs", type="float",
            dest="refrac_abs", help="Absolute refractory period (sec)")

    optParser.add_option("-t", "--time", type="float", dest="time",
            help="Total simulation time [seconds]")
    optParser.add_option("--sim_dt", type="float", dest="sim_dt",
            help="Simulation time step (s)")
    optParser.add_option("--spike_detect_th", type="float",
            dest="spike_detect_th", help="Spike detection threshold during"
            "numerical simulation (V)")

    optParser.add_option('--Vclamp', type=float,
            help="Clamp potential (for simulated voltage clamp, V)")

    optParser.add_option("--output_dir", type="string",
            dest="output_dir", help="Output directory path.")
    optParser.add_option("--fileNamePrefix", type="string", default='',
            help="Prefix to include for each output file")
    optParser.add_option("-u", '--update_interval', type="float",
            dest="update_interval", help="Duration between simulation status printouts (s)")
    optParser.add_option("-n", "--job_num", type="int", dest="job_num",
            help="Use argument of this option to specify the output file name number, instead of using time")

    return optParser

def setOptionDictionary(parser, options):
    d = {}
    attribs =  [x.dest for x in parser._get_all_options()[1:]]

    for attr_name in attribs:
        a = getattr(options, attr_name)
        if a is not None:
            d[attr_name] = a

    options._einet_optdict = d
    return options


# Save the results of the simulation in a .mat file
def saveResultsToMat(options, output_fname=None, spikeMon_e=None,
        spikeMon_i=None,  ratData=None, do_compression=True):

    # Directory and filenames
    
    if output_fname is None:
        timeSnapshot = datetime.now().strftime("%Y-%m-%dT%H-%M-%S")
        output_fname = options.output_dir
        if options.job_num != -1:
            output_fname = output_fname + '/job' + str(options.job_num)
        output_fname +=  '_' + timeSnapshot + '_output.mat'

    # Start saving everything
    if ratData == None:
        outData = {}
    else:
        outData = ratData

    if spikeMon_e is not None:
        outData['spikeCell_e'] = spikeMon_e.aspikes
    if spikeMon_i is not None:
        outData['spikeCell_i'] = spikeMon_i.aspikes

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

    outData['options'] = options._einet_optdict

    savemat(output_fname, outData, do_compression=do_compression)

