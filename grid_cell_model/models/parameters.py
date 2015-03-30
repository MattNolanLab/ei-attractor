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


class ParameterSelector(object):
    '''Allow the client to select portions of the EI network parameters.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        The parser which will receive the argument definitions.
    '''
    def __init__(self, parser):
        self._parser = parser

    @property
    def parser(self):
        '''Return the parser associated with this selector.'''
        return self._parser

    def simulation_params(self):
        '''Simulation parameters.'''
        self.parser.add_argument("--master_seed", type=int,   help="Master random number generator seed")
        self.parser.add_argument("--ntrials",     type=int,   help="Number of trials for the parameter set")
        self.parser.add_argument("--delay",       type=float, help="Synaptic delay (ms)")
        self.parser.add_argument("--nthreads",    type=int,   help="Number of threads (NEST)")
        self.parser.add_argument("--printTime",   type=int,   help="Whether to print time (0=False, 1=True)")
        self.parser.add_argument("--time",        type=float, help="Total simulation time (ms)")
        self.parser.add_argument("--sim_dt",      type=float, help="Simulation time step (ms)")

        self.parser.add_argument("--output_dir",     type=str,   help="Output directory path.")
        self.parser.add_argument("--fileNamePrefix", type=str,   default='', help="Prefix to include for each output file")
        self.parser.add_argument("--stateMonDur",    type=float, help="State monitors window duration (ms)")
        self.parser.add_argument("--job_num",        type=int,   help="Use argument of this option to specify the output file name number, instead of using time")

    def external_currents(self):
        '''External currents (theta) parameters.'''
        self.parser.add_argument("--Iext_e_const",  type=float,  help="Constant part of externals E current (pA)")
        self.parser.add_argument("--Iext_i_const",  type=float,  help="Constant part of externals I current (pA)")
        self.parser.add_argument("--Iext_e_theta",  type=float,  help="Theta E external current aplitude (pA)")
        self.parser.add_argument("--Iext_i_theta",  type=float,  help="Theta I external current aplitude (pA)")

        self.parser.add_argument("--theta_start_t", type=float,  help="Start time of theta stimulation (ms)")
        self.parser.add_argument("--theta_freq",    type=float,  help="Theta oscillation frequency (Hz)")


    def e_cell_params(self):
        '''E cell parameters.'''
        self.parser.add_argument("--taum_e",        type=float, help="Mean of excitatory neuron membrane time constant (ms)")
        self.parser.add_argument("--taum_e_spread", type=float, help="Spread of --taum_e (see --taum_i_spread option) (ms)")
        self.parser.add_argument("--EL_e",          type=float, help="Mean resting membrane potential of excitatory neurons (mV)")
        self.parser.add_argument("--EL_e_spread",   type=float, help="EL_e spread (see --taum_e_spread options for details) (mV)")
        self.parser.add_argument("--Vt_e",          type=float, help="Mean of the excitatory integrate and fire spiking threshold (mV)")
        self.parser.add_argument("--Vr_e",          type=float, help="Mean of the excitatory integrate and fire reset potential (mV)")
        self.parser.add_argument("--gL_e",          type=float, help="Mean of excitatory membrane resistance (nS)")
        self.parser.add_argument("--deltaT_e",      type=float, help="Sharpness of exponential I&F neuron spike initiation (excitatory, mV)")
        self.parser.add_argument("--E_AHP_e",       type=float, help="Excitatory AHP reversal potential (mV)")
        self.parser.add_argument("--g_AHP_e_max",   type=float, help="Maximal AHP excitatory conductance (nS)")
        self.parser.add_argument("--tau_AHP_e",     type=float, help="Excitatory AHP time decay constant (ms)")
        self.parser.add_argument("--t_ref_e",       type=float, help="Excitatory absolute refractory period (msec)")
        self.parser.add_argument("--V_peak_e",      type=float, help="Excitatory spike detection threshold (mV)")

    def i_cell_params(self):
        '''I cell parameters.'''
        self.parser.add_argument("--taum_i",        type=float, help="Mean of inhibitory neuron membrane time constant (ms)")
        self.parser.add_argument("--taum_i_spread", type=float, help="Spread of --taum_i: drawn from uniform distrib. <taum_i-spread/2, taum_i+spread/2> (ms)")
        self.parser.add_argument("--EL_i",          type=float, help="Mean resting membrane potential of inhibitory neurons (mV)")
        self.parser.add_argument("--EL_i_spread",   type=float, help="ELs will be uniformly generated from <EL-spread, EL+spread> (mV)")
        self.parser.add_argument("--Vt_i",          type=float, help="Mean of the excitatory integrate and fire spiking threshold (mV)")
        self.parser.add_argument("--Vr_i",          type=float, help="Mean of the inhibitory integrate and fire reset potential (mV)")
        self.parser.add_argument("--gL_i",          type=float, help="Mean of inhibitory membrane resistance (nS)")
        self.parser.add_argument("--ad_tau_i_mean", type=float, help="Mean of inhibitory adaptation time constant (ms)")
        self.parser.add_argument("--ad_tau_i_std",  type=float, help="Std. deviation of inhibitory adaptation time constant (ms)")
        self.parser.add_argument("--ad_i_g_inc",    type=float, help="After-spike inhibitory increase of leak conductance (nS)")
        self.parser.add_argument("--deltaT_i",      type=float, help="Sharpness of exponential I&F neuron spike initiation (mV)")
        self.parser.add_argument("--t_ref_i",       type=float, help="Inhibitory absolute refractory period (msec)")
        self.parser.add_argument("--E_AHP_i",       type=float, help="Inhibitory AHP reversal potential (mV)")
        self.parser.add_argument("--g_AHP_i_max",   type=float, help="Maximal AHP inhibitory conductance (nS)")
        self.parser.add_argument("--tau_AHP_i",     type=float, help="Inhibitory AHP time decay constant (ms)")
        self.parser.add_argument("--V_peak_i",      type=float, help="Inhibitory spike detection threshold (mV)")

    def network_size(self):
        '''Network size parameters.'''
        self.parser.add_argument("--Ne", type=int,    help="Number of excitatory neurons")
        self.parser.add_argument("--Ni", type=int,    help="Number of inhibitory neurons")

    def noise(self):
        '''Noise levels.'''
        self.parser.add_argument("--noise_sigma", type=float,  help="Std. dev of neural noise (V)")

    def velocity_inputs(self):
        '''Velocity input parameters.'''
        self.parser.add_argument("--Ivel",             type=float,  help="Velocity input (In constant velocity currents, pA)")
        self.parser.add_argument("--bumpCurrentSlope", type=float,  help="Slope of the bump-current linear relationsip (starts at zero, pA/(cm/s)")
        self.parser.add_argument("--ratVelFName",      type=str,    help="Positional input file name (matlab)")

    def place_cells(self):
        '''Place cell parameters.'''
        self.parser.add_argument("--N_place_cells",        type=int,    help="Number of place cells connected to grid cells in one dimension. The total number of place cells will be the square of this parameter.")
        self.parser.add_argument("--pc_max_rate",          type=float,  help="Place cell maximal firing rate (Hz)")
        self.parser.add_argument("--pc_field_std",         type=float,  help="Std. deviation of the Gaussian specifying the place field size (cm)")
        self.parser.add_argument("--pc_conn_weight",       type=float,  help="Connection weight from each place to grid cells (nS)")
        self.parser.add_argument("--pc_start_max_rate",    type=float,  help="Init place cell maximal firing rate (Hz)")
        self.parser.add_argument("--pc_start_conn_weight", type=float,  help="Connection weight from each initialisation place cell to grid cells (nS)")

    def synapse_properties(self):
        '''Properties of synapses.'''
        self.parser.add_argument("--tau_AMPA",        type=float,  help="Mean of AMPA synaptic conductance time constant (ms)")
        self.parser.add_argument("--tau_NMDA_rise",   type=float,  help="NMDA rise time constant (ms)")
        self.parser.add_argument("--tau_NMDA_fall",   type=float,  help="NMDA fall time constant (ms)")
        self.parser.add_argument("--tau_GABA_A_rise", type=float,  help="Mean of GABA A rising time constant (ms)")
        self.parser.add_argument("--tau_GABA_A_fall", type=float,  help="Mean of GABA A fall time constant (s)")

        self.parser.add_argument("--E_AMPA",   type=float,  help="AMPA reversal potential (V)")
        self.parser.add_argument("--E_GABA_A", type=float,  help="GABA A reversal potential (V)")

        self.parser.add_argument("--C_Mg",        type=float, help="Mg2+ concentration; used for NMDA voltage dependence")
        self.parser.add_argument("--NMDA_amount", type=float, help="NMDA portion relative to AMPA (%%)")

    def preferred_directions(self):
        '''Preferred directions.'''
        self.parser.add_argument("--prefDirC_e",  type=float, help="Excitatory (E-->I) preferred direction multiplier")
        self.parser.add_argument("--prefDirC_ee", type=float, help="Excitatory (E-->E) preferred direction multiplier")
        self.parser.add_argument("--prefDirC_i",  type=float, help="Inhibitory (I-->E) preferred direction multiplier")

    def spatial_properties(self):
        '''Properties of grid fields and arenas.'''
        self.parser.add_argument("--arenaSize", type=float, help="Size of the arena where the rat runs (cm)")
        self.parser.add_argument("--gridSep",   type=float, help="Distance between grid field peaks (cm)")

    def recording_parameters(self):
        '''Simulation result recordings.'''
        self.parser.add_argument('--Vclamp',       type=float,  help="Clamp potential (for simulated voltage clamp, V)")

        self.parser.add_argument("--gammaNSample", type=int,    help="Number of neurons in the network to sample from, for the frequency analysis.")
        self.parser.add_argument("--connNE",       type=int,    help="Number of E neurons to sample from, to save input connection weights.")
        self.parser.add_argument("--connNI",       type=int,    help="Number of I neurons to sample from, to save input connection weights.")

    def ei_profile(self):
        '''Properties of E --> I synaptic profiles.'''
        self.parser.add_argument("--EI_flat",
                                 type=int, choices=[0, 1], required=True,
                                 help="Whether the E-->I profile is flat or distance dependent.")
        self.parser.add_argument("--AMPA_gaussian",    type=float, help="If 1, AMPA profile will be Gaussian, if 0, ring-like. Can be used to swap connectivity types")
        self.parser.add_argument("--pAMPA_mu",         type=float, help="AMPA profile center (normalised)")
        self.parser.add_argument("--pAMPA_sigma",      type=float, help="AMPA profile spread (normalised)")
        self.parser.add_argument("--g_AMPA_total",     type=float, help="Total AMPA connection synaptic conductance (nS)")
        self.parser.add_argument("--g_uni_AMPA_total", type=float, help="Total AMPA connections synaptic conductance (nS)")
        self.parser.add_argument("--uni_AMPA_density", type=float, help="Density of uniform AMPA connections (fraction)")

    def ie_profile(self):
        '''Properties of I --> E synaptic profiles.'''
        self.parser.add_argument("--IE_flat",
                                 type=int, choices=[0, 1], required=True,
                                 help="Whether the I-->E profile is flat or distance dependent.")
        self.parser.add_argument("--pGABA_mu",         type=float, help="GABA A profile center (normalised)")
        self.parser.add_argument("--pGABA_sigma",      type=float, help="GABA A profile spread (normalised)")

        self.parser.add_argument("--g_GABA_total",     type=float, help="Total GABA connections synaptic conductance (nS)")
        self.parser.add_argument("--g_uni_GABA_frac",  type=float, help="Total uniform GABA A connections synaptic conductance (fraction of g_GABA_total)")
        self.parser.add_argument("--uni_GABA_density", type=float, help="Density of uniform GABA A connections")

    def ee_profile(self):
        '''Properties of E --> E synaptic profiles.

        In the current network model, E --> E connections are always
        Gaussian-like.
        '''
        self.parser.add_argument("--use_EE",     type=int, choices=[0, 1], help="Whether to use the E-->E connectivity profiles.")
        self.parser.add_argument("--pEE_sigma",  type=float, help="E-->E profile spread (normalised).")
        self.parser.add_argument("--g_EE_total", type=float, help="Total AMPA amount for the E-->E connections.")
        self.parser.add_argument("--g_EI_uni_density", type=float, help="Probability of an E-->I connection.")
        self.parser.add_argument("--g_IE_uni_density", type=float, help="Probability of an I-->E connection.")


def getOptParser():
    '''Create a generic option parser.

    .. note::

        This is only for compatibility with old submission code. Do not use
        this otherwise.
    '''
    s = ParameterSelector(DOptionParser())
    s.simulation_params()
    s.external_currents()
    s.e_cell_params()
    s.i_cell_params()
    s.network_size()
    s.noise()
    s.velocity_inputs()
    s.place_cells()
    s.synapse_properties()
    s.preferred_directions()
    s.spatial_properties()
    s.recording_parameters()
    s.ei_profile()
    s.ie_profile()
    s.ee_profile()

    return s.parser
