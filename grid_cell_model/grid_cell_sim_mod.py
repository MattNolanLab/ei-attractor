#
#   grid_cell_sim_mod.py
#   
#   This file defines some basic simulation functions, independent of modelling.
#   
#
#       Copyright (C) 2012  Lukas Solanka <l.solanka@sms.ed.ac.uk>
#       
#       This program is free software: you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation, either version 3 of the License, or
#       (at your option) any later version.
#       
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#       
#       You should have received a copy of the GNU General Public License
#       along with this program.  If not, see <http://www.gnu.org/licenses/>.
#


from scipy.io import savemat
from optparse import OptionParser
from datetime import datetime



def getOptParser():
    # Network parameters definition
    optParser = OptionParser()
    optParser.add_option("--Ne",              type="int",    help="Number of excitatory neurons")
    optParser.add_option("--Ni",              type="int",    help="Number of inhibitory neurons")
    optParser.add_option("--ntrials",         type="int",    help="Number of trials for the parameter set")

    optParser.add_option("--AMPA_density",    type="float",  help="Density of E-->I connections")
    optParser.add_option("--GABA_density",    type="float",  help="Density of I-->E connections")

    optParser.add_option("--Iext_e",          type="float",  help="External background current into the excitatory population (pA)")
    optParser.add_option("--Iext_i",          type="float",  help="External background current into the inhibitory population (pA)")

    optParser.add_option("--taum_e",          type="float",  help="Mean of excitatory neuron membrane time constant (ms)")
    optParser.add_option("--taum_e_spread",   type="float",  help="Spread of --taum_e (see --taum_i_spread option) (ms)")
    optParser.add_option("--EL_e",            type="float",  help="Mean resting membrane potential of excitatory neurons (mV)")
    optParser.add_option("--EL_e_spread",     type="float",  help="EL_e spread (see --taum_e_spread options for details) (mV)")
    optParser.add_option("--Vt_e",            type="float",  help="Mean of the excitatory integrate and fire spiking threshold (mV)")
    optParser.add_option("--Vr_e",            type="float",  help="Mean of the excitatory integrate and fire reset potential (mV)")
    optParser.add_option("--gL_e",            type="float",  help="Mean of excitatory membrane resistance (nS)")
    optParser.add_option("--deltaT_e",        type="float",  help="Sharpness of exponential I&F neuron spike initiation (excitatory, mV)")
    optParser.add_option("--E_AHP_e",         type="float",  help="Excitatory AHP reversal potential (mV)")
    optParser.add_option("--g_AHP_e",         type="float",  help="Maximal AHP excitatory conductance (nS)")
    optParser.add_option("--tau_ahp_e",       type="float",  help="Excitatory AHP time decay constant (ms)")

    optParser.add_option("--taum_i",          type="float",  help="Mean of inhibitory neuron membrane time constant (ms)")
    optParser.add_option("--taum_i_spread",   type="float",  help="Spread of --taum_i: drawn from uniform distrib. <taum_i-spread/2, taum_i+spread/2> (ms)")
    optParser.add_option("--EL_i",            type="float",  help="Mean resting membrane potential of inhibitory neurons (mV)")
    optParser.add_option("--EL_i_spread",     type="float",  help="ELs will be uniformly generated from <EL-spread, EL+spread> (mV)")
    optParser.add_option("--Vt_i",            type="float",  help="Mean of the excitatory integrate and fire spiking threshold (mV)")
    optParser.add_option("--Vr_i",            type="float",  help="Mean of the inhibitory integrate and fire reset potential (mV)")
    optParser.add_option("--gL_i",            type="float",  help="Mean of inhibitory membrane resistance (nS)")
    optParser.add_option("--ad_tau_i_mean",  type="float",  help="Mean of inhibitory adaptation time constant (ms)")
    optParser.add_option("--ad_tau_i_std",   type="float",  help="Std. deviation of inhibitory adaptation time constant (ms)")
    optParser.add_option("--ad_i_g_inc",     type="float",  help="After-spike inhibitory increase of leak conductance (nS)")
    optParser.add_option("--deltaT_i",        type="float",  help="Sharpness of exponential I&F neuron spike initiation (mV)")

    optParser.add_option("--tau_AMPA",        type="float",  help="Mean of AMPA synaptic conductance time constant (ms)")
    optParser.add_option("--g_AMPA_total",    type="float",  help="Total AMPA connection synaptic conductance (nS)")
#    optParser.add_option("--g_AMPA_std",     type="float",  help="Std. deviation of AMPA connections synaptic conductance (S)")
    optParser.add_option("--g_GABA_total",    type="float",  help="Total GABA connections synaptic conductance (S)")
    optParser.add_option("--NMDA_percent",    type="float",  help="Percentage of NMDA conductance the excitatory synapse contains (%)")
    optParser.add_option("--tau_NMDA_rise",   type="float",  help="NMDA rise time constant (ms)")
    optParser.add_option("--tau_NMDA_fall",   type="float",  help="NMDA fall time constant (ms)")
    optParser.add_option("--tau_GABA_rise",   type="float",  help="Mean of GABA rising time constant (ms)")
    optParser.add_option("--tau_GABA_fall",   type="float",  help="Mean of GABA fall time constant (s)")

    optParser.add_option("--Vrev_AMPA",       type="float",  help="AMPA reversal potential (V)")
    optParser.add_option("--Vrev_GABA",       type="float",  help="GABA reversal potential (V)")

    optParser.add_option("--noise_sigma",     type="float",  help="Std. dev of neural noise (V)")
    optParser.add_option("--sigma_init_cond", type="float",  help="Std. dev distribution of initial membrane voltages (V)")

    optParser.add_option("--refrac_abs",      type="float",  help="Absolute refractory period (sec)")

    optParser.add_option("-t", "--time",      type="float",  help="Total simulation time [seconds]")
    optParser.add_option("--sim_dt",          type="float",  help="Simulation time step (s)")
    optParser.add_option("--spike_detect_th", type="float",  help="Spike detection threshold during numerical simulation (V)")

    optParser.add_option('--Vclamp',          type=float,    help="Clamp potential (for simulated voltage clamp, V)")

    optParser.add_option("--output_dir",      type="string", help="Output directory path.")
    optParser.add_option("--fileNamePrefix",  type="string", default='', help="Prefix to include for each output file")
    optParser.add_option('--update_interval', type="float",  help="Duration between simulation status printouts (s)")
    optParser.add_option("--job_num",         type="int",    help="Use argument of this option to specify the output file name number, instead of using time")

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

