'''Command line options.'''
from __future__ import absolute_import, print_function, division

import copy
from ..submitting import flagparse


class DOptionParser(flagparse.FlagParser):
    '''
    Default EI network option parser. This parser adds a special flag
    ``--verbosity`` into the list of options, and sets the logging level
    according to this option (as specified in the python logging module).
    '''

    def __init__(self):
        super(DOptionParser, self).__init__()


    def parse_args(self, args=None, namespace=None):
        options = super(DOptionParser, self).parse_args(args, namespace)
        options = self.setOptionDictionary(options)
        return options, None


    def setOptionDictionary(self, options):
        optDict = copy.deepcopy(vars(options))
        for key in optDict.keys():
            if optDict[key] is None:
                del optDict[key]
        options._einet_optdict = optDict
        return options



def getOptParser():
    # Network parameters definition
    optParser = DOptionParser()
    optParser.add_argument("--master_seed",           type=int,    help="Master random number generator seed")
    optParser.add_argument("--Ne",                    type=int,    help="Number of excitatory neurons")
    optParser.add_argument("--Ni",                    type=int,    help="Number of inhibitory neurons")
    optParser.add_argument("--ntrials",               type=int,    help="Number of trials for the parameter set")
    optParser.add_argument("--delay",                 type=float,  help="Synaptic delay (ms)")
    optParser.add_argument("--nthreads",              type=int,    help="Number of threads (NEST)")
    optParser.add_argument("--printTime",             type=int,    help="Whether to print time (0=False, 1=True)")

    optParser.add_argument("--Ivel",                  type=float,  help="Velocity input (In constant velocity currents, pA)")
    optParser.add_argument("--bumpCurrentSlope",      type=float,  help="Slope of the bump-current linear relationsip (starts at zero, pA/(cm/s)")
    optParser.add_argument("--ratVelFName",           type=str,    help="Positional input file name (matlab)")

    optParser.add_argument("--Iext_e_const",          type=float,  help="Constant part of externals E current (pA)")
    optParser.add_argument("--Iext_i_const",          type=float,  help="Constant part of externals I current (pA)")
    optParser.add_argument("--Iext_e_theta",          type=float,  help="Theta E external current aplitude (pA)")
    optParser.add_argument("--Iext_i_theta",          type=float,  help="Theta I external current aplitude (pA)")

    optParser.add_argument("--theta_start_t",         type=float,  help="Start time of theta stimulation (ms)")
    optParser.add_argument("--theta_freq",            type=float,  help="Theta oscillation frequency (Hz)")

    optParser.add_argument("--AMPA_gaussian",         type=float, help="If 1, AMPA profile will be gauusian, if 0, ring-like. Can be used to swap connectivity types")
    optParser.add_argument("--pAMPA_mu",              type=float, help="AMPA profile center (normalised)")
    optParser.add_argument("--pAMPA_sigma",           type=float, help="AMPA profile spread (normalised)")
    optParser.add_argument("--pGABA_mu",              type=float, help="GABA A profile center (normalised)")
    optParser.add_argument("--pGABA_sigma",           type=float, help="GABA A profile spread (normalised)")
    optParser.add_argument("--NMDA_amount",           type=float, help="NMDA portion relative to AMPA (%%)")
    optParser.add_argument("--C_Mg",                  type=float, help="Mg2+ concentration; used for NMDA voltage dependence")

    optParser.add_argument("--prefDirC_e",            type=float, help="Excitatory (E-->I) preferred direction multiplier")
    optParser.add_argument("--prefDirC_i",            type=float, help="Inhibitory (I-->E) preferred direction multiplier")
    optParser.add_argument("--arenaSize",             type=float, help="Size of the arena where the rat runs (cm)")
    optParser.add_argument("--gridSep",               type=float, help="Distance between grid field peaks (cm)")

    optParser.add_argument("--taum_e",                type=float,  help="Mean of excitatory neuron membrane time constant (ms)")
    optParser.add_argument("--taum_e_spread",         type=float,  help="Spread of --taum_e (see --taum_i_spread option) (ms)")
    optParser.add_argument("--EL_e",                  type=float,  help="Mean resting membrane potential of excitatory neurons (mV)")
    optParser.add_argument("--EL_e_spread",           type=float,  help="EL_e spread (see --taum_e_spread options for details) (mV)")
    optParser.add_argument("--Vt_e",                  type=float,  help="Mean of the excitatory integrate and fire spiking threshold (mV)")
    optParser.add_argument("--Vr_e",                  type=float,  help="Mean of the excitatory integrate and fire reset potential (mV)")
    optParser.add_argument("--gL_e",                  type=float,  help="Mean of excitatory membrane resistance (nS)")
    optParser.add_argument("--deltaT_e",              type=float,  help="Sharpness of exponential I&F neuron spike initiation (excitatory, mV)")
    optParser.add_argument("--E_AHP_e",               type=float,  help="Excitatory AHP reversal potential (mV)")
    optParser.add_argument("--g_AHP_e_max",           type=float,  help="Maximal AHP excitatory conductance (nS)")
    optParser.add_argument("--tau_AHP_e",             type=float,  help="Excitatory AHP time decay constant (ms)")
    optParser.add_argument("--t_ref_e",               type=float,  help="Excitatory absolute refractory period (msec)")
    optParser.add_argument("--V_peak_e",              type=float,  help="Excitatory spike detection threshold (mV)")

    optParser.add_argument("--taum_i",                type=float,  help="Mean of inhibitory neuron membrane time constant (ms)")
    optParser.add_argument("--taum_i_spread",         type=float,  help="Spread of --taum_i: drawn from uniform distrib. <taum_i-spread/2, taum_i+spread/2> (ms)")
    optParser.add_argument("--EL_i",                  type=float,  help="Mean resting membrane potential of inhibitory neurons (mV)")
    optParser.add_argument("--EL_i_spread",           type=float,  help="ELs will be uniformly generated from <EL-spread, EL+spread> (mV)")
    optParser.add_argument("--Vt_i",                  type=float,  help="Mean of the excitatory integrate and fire spiking threshold (mV)")
    optParser.add_argument("--Vr_i",                  type=float,  help="Mean of the inhibitory integrate and fire reset potential (mV)")
    optParser.add_argument("--gL_i",                  type=float,  help="Mean of inhibitory membrane resistance (nS)")
    optParser.add_argument("--ad_tau_i_mean",         type=float,  help="Mean of inhibitory adaptation time constant (ms)")
    optParser.add_argument("--ad_tau_i_std",          type=float,  help="Std. deviation of inhibitory adaptation time constant (ms)")
    optParser.add_argument("--ad_i_g_inc",            type=float,  help="After-spike inhibitory increase of leak conductance (nS)")
    optParser.add_argument("--deltaT_i",              type=float,  help="Sharpness of exponential I&F neuron spike initiation (mV)")
    optParser.add_argument("--t_ref_i",               type=float,  help="Inhibitory absolute refractory period (msec)")
    optParser.add_argument("--E_AHP_i",               type=float,  help="Inhibitory AHP reversal potential (mV)")
    optParser.add_argument("--g_AHP_i_max",           type=float,  help="Maximal AHP inhibitory conductance (nS)")
    optParser.add_argument("--tau_AHP_i",             type=float,  help="Inhibitory AHP time decay constant (ms)")
    optParser.add_argument("--V_peak_i",              type=float,  help="Inhibitory spike detection threshold (mV)")

    optParser.add_argument("--tau_AMPA",              type=float,  help="Mean of AMPA synaptic conductance time constant (ms)")
    optParser.add_argument("--g_AMPA_total",          type=float,  help="Total AMPA connection synaptic conductance (nS)")
    optParser.add_argument("--g_AMPA_std",            type=float,  help="Std. deviation of AMPA connections synaptic conductance (nS)")
    optParser.add_argument("--g_uni_AMPA_total",      type=float,  help="Total AMPA connections synaptic conductance (nS)")
    optParser.add_argument("--uni_AMPA_density",      type=float,  help="Density of uniform AMPA connections (fraction)")
    optParser.add_argument("--g_GABA_total",          type=float,  help="Total GABA connections synaptic conductance (nS)")
    optParser.add_argument("--g_uni_GABA_frac",       type=float,  help="Total uniform GABA A connections synaptic conductance (fraction of g_GABA_total)")
    optParser.add_argument("--uni_GABA_density",      type=float,  help="Density of uniform GABA A connections")
    optParser.add_argument("--tau_NMDA_rise",         type=float,  help="NMDA rise time constant (ms)")
    optParser.add_argument("--tau_NMDA_fall",         type=float,  help="NMDA fall time constant (ms)")
    optParser.add_argument("--tau_GABA_A_rise",       type=float,  help="Mean of GABA A rising time constant (ms)")
    optParser.add_argument("--tau_GABA_A_fall",       type=float,  help="Mean of GABA A fall time constant (s)")

    optParser.add_argument("--E_AMPA",                type=float,  help="AMPA reversal potential (V)")
    optParser.add_argument("--E_GABA_A",              type=float,  help="GABA A reversal potential (V)")

    optParser.add_argument("--N_place_cells",         type=int,    help="Number of place cells connected to grid cells in one dimension. The total number of place cells will be the square of this parameter.")
    optParser.add_argument("--pc_max_rate",           type=float,  help="Place cell maximal firing rate (Hz)")
    optParser.add_argument("--pc_field_std",          type=float,  help="Std. deviation of the Gaussian specifying the place field size (cm)")
    optParser.add_argument("--pc_conn_weight",        type=float,  help="Connection weight from each place to grid cells (nS)")
    optParser.add_argument("--pc_start_max_rate",     type=float,  help="Init place cell maximal firing rate (Hz)")
    optParser.add_argument("--pc_start_conn_weight",  type=float,  help="Connection weight from each initialisation place cell to grid cells (nS)")

    optParser.add_argument("--noise_sigma",           type=float,  help="Std. dev of neural noise (V)")


    optParser.add_argument("--time",                  type=float,  help="Total simulation time (ms)")
    optParser.add_argument("--sim_dt",                type=float,  help="Simulation time step (ms)")

    optParser.add_argument('--Vclamp',                type=float,  help="Clamp potential (for simulated voltage clamp, V)")

    optParser.add_argument("--output_dir",            type=str,    help="Output directory path.")
    optParser.add_argument("--fileNamePrefix",        type=str,    default='', help="Prefix to include for each output file")
    optParser.add_argument("--stateMonDur",           type=float,  help="State monitors window duration (ms)")
    optParser.add_argument("--job_num",               type=int,    help="Use argument of this option to specify the output file name number, instead of using time")

    optParser.add_argument("--gammaNSample",          type=int,    help="Number of neurons in the network to sample from, for the frequency analysis.")
    optParser.add_argument("--connNE",                type=int,    help="Number of E neurons to sample from, to save input connection weights.")
    optParser.add_argument("--connNI",                type=int,    help="Number of I neurons to sample from, to save input connection weights.")
    return optParser


