#
#   parameters.py
#   
#   Command line options
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


from optparse import OptionParser


class DOptionParser(OptionParser):

    def __init__(self):
        OptionParser.__init__(self)


    def parse_args(self):
        options, args = OptionParser.parse_args(self)
        options = self.setOptionDictionary(options)
        return options, args


    def setOptionDictionary(self, options):
        d = {}
        attribs =  [x.dest for x in self._get_all_options()[1:]]
    
        for attr_name in attribs:
            a = getattr(options, attr_name)
            if a is not None:
                d[attr_name] = a
    
        options._einet_optdict = d
        return options



def getOptParser():
    # Network parameters definition
    optParser = DOptionParser()
    optParser.add_option("--Ne",                   type="int",    help="Number of excitatory neurons")
    optParser.add_option("--Ni",                   type="int",    help="Number of inhibitory neurons")
    optParser.add_option("--ntrials",              type="int",    help="Number of trials for the parameter set")
    optParser.add_option("--delay",                type="float",  help="Synaptic delay (ms)")
    optParser.add_option("--nthreads",             type="int",    help="Number of threads (NEST)")
    optParser.add_option("--printTime",            type="int",    help="Whether to print time (0=False, 1=True)")

    optParser.add_option("--Ivel",                 type="float",  help="Velocity input (In constant velocity currents, pA)")
    optParser.add_option("--bumpCurrentSlope",     type="float",  help="Slope of the bump-current linear relationsip (starts at zero, pA/(cm/s)")
    optParser.add_option("--ratVelFName",          type="string", help="Positional input file name (matlab)")

    optParser.add_option("--Iext_e_const",         type="float",  help="Constant part of externals E current (pA)")
    optParser.add_option("--Iext_i_const",         type="float",  help="Constant part of externals I current (pA)")
    optParser.add_option("--Iext_e_theta",         type="float",  help="Theta E external current aplitude (pA)")
    optParser.add_option("--Iext_i_theta",         type="float",  help="Theta I external current aplitude (pA)")

    optParser.add_option("--theta_start_t",        type="float",  help="Start time of theta stimulation (ms)")
    optParser.add_option("--theta_freq",           type="float",  help="Theta oscillation frequency (Hz)")

    optParser.add_option("--AMPA_gaussian",     type="float", help="If 1, AMPA profile will be gauusian, if 0, ring-like. Can be used to swap connectivity types")
    optParser.add_option("--pAMPA_mu",          type="float", help="AMPA profile center (normalised)")
    optParser.add_option("--pAMPA_sigma",       type="float", help="AMPA profile spread (normalised)")
    optParser.add_option("--pGABA_mu",          type="float", help="GABA A profile center (normalised)")
    optParser.add_option("--pGABA_sigma",       type="float", help="GABA A profile spread (normalised)")
    optParser.add_option("--NMDA_amount",       type="float", help="NMDA portion relative to AMPA (%)")

    optParser.add_option("--prefDirC_e",        type="float", help="Excitatory (E-->I) preferred direction multiplier")
    optParser.add_option("--prefDirC_i",        type="float", help="Inhibitory (I-->E) preferred direction multiplier")
    optParser.add_option("--arenaSize",         type="float", help="Size of the arena where the rat runs (cm)")
    optParser.add_option("--gridSep",           type="float", help="Distance between grid field peaks (cm)")

    optParser.add_option("--taum_e",            type="float",  help="Mean of excitatory neuron membrane time constant (ms)")
    optParser.add_option("--taum_e_spread",     type="float",  help="Spread of --taum_e (see --taum_i_spread option) (ms)")
    optParser.add_option("--EL_e",              type="float",  help="Mean resting membrane potential of excitatory neurons (mV)")
    optParser.add_option("--EL_e_spread",       type="float",  help="EL_e spread (see --taum_e_spread options for details) (mV)")
    optParser.add_option("--Vt_e",              type="float",  help="Mean of the excitatory integrate and fire spiking threshold (mV)")
    optParser.add_option("--Vr_e",              type="float",  help="Mean of the excitatory integrate and fire reset potential (mV)")
    optParser.add_option("--gL_e",              type="float",  help="Mean of excitatory membrane resistance (nS)")
    optParser.add_option("--deltaT_e",          type="float",  help="Sharpness of exponential I&F neuron spike initiation (excitatory, mV)")
    optParser.add_option("--E_AHP_e",           type="float",  help="Excitatory AHP reversal potential (mV)")
    optParser.add_option("--g_AHP_e_max",       type="float",  help="Maximal AHP excitatory conductance (nS)")
    optParser.add_option("--tau_AHP_e",         type="float",  help="Excitatory AHP time decay constant (ms)")
    optParser.add_option("--t_ref_e",           type="float",  help="Excitatory absolute refractory period (msec)")
    optParser.add_option("--V_peak_e",          type="float",  help="Excitatory spike detection threshold (mV)")

    optParser.add_option("--taum_i",            type="float",  help="Mean of inhibitory neuron membrane time constant (ms)")
    optParser.add_option("--taum_i_spread",     type="float",  help="Spread of --taum_i: drawn from uniform distrib. <taum_i-spread/2, taum_i+spread/2> (ms)")
    optParser.add_option("--EL_i",              type="float",  help="Mean resting membrane potential of inhibitory neurons (mV)")
    optParser.add_option("--EL_i_spread",       type="float",  help="ELs will be uniformly generated from <EL-spread, EL+spread> (mV)")
    optParser.add_option("--Vt_i",              type="float",  help="Mean of the excitatory integrate and fire spiking threshold (mV)")
    optParser.add_option("--Vr_i",              type="float",  help="Mean of the inhibitory integrate and fire reset potential (mV)")
    optParser.add_option("--gL_i",              type="float",  help="Mean of inhibitory membrane resistance (nS)")
    optParser.add_option("--ad_tau_i_mean",     type="float",  help="Mean of inhibitory adaptation time constant (ms)")
    optParser.add_option("--ad_tau_i_std",      type="float",  help="Std. deviation of inhibitory adaptation time constant (ms)")
    optParser.add_option("--ad_i_g_inc",        type="float",  help="After-spike inhibitory increase of leak conductance (nS)")
    optParser.add_option("--deltaT_i",          type="float",  help="Sharpness of exponential I&F neuron spike initiation (mV)")
    optParser.add_option("--t_ref_i",           type="float",  help="Inhibitory absolute refractory period (msec)")
    optParser.add_option("--E_AHP_i",           type="float",  help="Inhibitory AHP reversal potential (mV)")
    optParser.add_option("--g_AHP_i_max",       type="float",  help="Maximal AHP inhibitory conductance (nS)")
    optParser.add_option("--tau_AHP_i",         type="float",  help="Inhibitory AHP time decay constant (ms)")
    optParser.add_option("--V_peak_i",          type="float",  help="Inhibitory spike detection threshold (mV)")

    optParser.add_option("--tau_AMPA",          type="float",  help="Mean of AMPA synaptic conductance time constant (ms)")
    optParser.add_option("--g_AMPA_total",      type="float",  help="Total AMPA connection synaptic conductance (nS)")
    optParser.add_option("--g_AMPA_std",        type="float",  help="Std. deviation of AMPA connections synaptic conductance (nS)")
    optParser.add_option("--g_uni_AMPA_total",  type="float",  help="Total AMPA connections synaptic conductance (nS)")
    optParser.add_option("--uni_AMPA_density",  type="float",  help="Density of uniform AMPA connections (fraction)")
    optParser.add_option("--g_GABA_total",      type="float",  help="Total GABA connections synaptic conductance (nS)")
    optParser.add_option("--g_uni_GABA_frac",   type="float",  help="Total uniform GABA A connections synaptic conductance (fraction of g_GABA_total)")
    optParser.add_option("--uni_GABA_density",  type="float",  help="Density of uniform GABA A connections")
    optParser.add_option("--NMDA_percent",      type="float",  help="Percentage of NMDA conductance the excitatory synapse contains (%)")
    optParser.add_option("--tau_NMDA_rise",     type="float",  help="NMDA rise time constant (ms)")
    optParser.add_option("--tau_NMDA_fall",     type="float",  help="NMDA fall time constant (ms)")
    optParser.add_option("--tau_GABA_A_rise",   type="float",  help="Mean of GABA A rising time constant (ms)")
    optParser.add_option("--tau_GABA_A_fall",   type="float",  help="Mean of GABA A fall time constant (s)")

    optParser.add_option("--E_AMPA",            type="float",  help="AMPA reversal potential (V)")
    optParser.add_option("--E_GABA_A",          type="float",  help="GABA A reversal potential (V)")

    optParser.add_option("--N_place_cells",     type="int",  help="Number of place cells connected to grid cells in one dimension. The total number of place cells will be the square of this parameter.")
    optParser.add_option("--pc_max_rate",       type="float",  help="Place cell maximal firing rate (Hz)")
    optParser.add_option("--pc_field_std",      type="float",  help="Std. deviation of the Gaussian specifying the place field size (cm)")
    optParser.add_option("--pc_conn_weight",    type="float",  help="Connection weight from each place to grid cells (nS)")
    optParser.add_option("--pc_start_max_rate", type="float",  help="Init place cell maximal firing rate (Hz)")
    optParser.add_option("--pc_start_conn_weight", type="float",  help="Connection weight from each initialisation place cell to grid cells (nS)")

    optParser.add_option("--noise_sigma",       type="float",  help="Std. dev of neural noise (V)")


    optParser.add_option("--time",              type="float",  help="Total simulation time (ms)")
    optParser.add_option("--sim_dt",            type="float",  help="Simulation time step (ms)")

    optParser.add_option('--Vclamp',            type=float,    help="Clamp potential (for simulated voltage clamp, V)")

    optParser.add_option("--output_dir",        type="string", help="Output directory path.")
    optParser.add_option("--fileNamePrefix",    type="string", default='', help="Prefix to include for each output file")
    optParser.add_option("--stateMonDur",       type="float",  help="State monitors window duration (ms)")
    optParser.add_option("--job_num",           type="int",    help="Use argument of this option to specify the output file name number, instead of using time")

    optParser.add_option("--gammaNSample",      type="int",    help="Number of neurons in the network to sample from, for the frequency analysis.")
    return optParser


