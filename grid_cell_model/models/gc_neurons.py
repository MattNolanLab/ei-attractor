'''Grid cell network neuron definitions.

.. currentmodule:: grid_cell_model.models.gc_neurons

A helper module to set up parameters of E and I populations.

Functions
---------

.. autosummary::

    getENeuronParams
    getINeuronParams

'''
from __future__ import absolute_import, print_function, division

import numpy as np

__all__ = ['getENeuronParams', 'getINeuronParams']


def getENeuronParams(no):
    '''
    Return a dictionary of E neuron parameters, using `no`: neuron options
    object.
    '''
    return {
            "V_m"              : no.EL_e,
            "C_m"              : no.taum_e * no.gL_e,
            "t_ref"            : no.t_ref_e,
            "V_peak"           : no.V_peak_e,
            "V_reset"          : no.Vr_e,
            "E_L"              : no.EL_e,
            "g_L"              : no.gL_e,
            "Delta_T"          : no.deltaT_e,
            "V_th"             : no.Vt_e,
            "E_AMPA"           : no.E_AMPA,
            "E_GABA_A"         : no.E_GABA_A,
            "tau_AMPA_fall"    : no.tau_AMPA,
            "tau_NMDA_fall"    : no.tau_NMDA_fall,
            "tau_GABA_A_fall"  : no.tau_GABA_A_fall,
            "tau_AHP"          : no.tau_AHP_e,
            "E_AHP"            : no.E_AHP_e,
            "g_AHP_max"        : no.g_AHP_e_max,
            "g_AHP_ad"         : False,
            "I_const"          : no.Iext_e_const,
            "I_ac_amp"         : no.Iext_e_theta,
            "I_ac_freq"        : no.theta_freq,
            "I_ac_phase"       : -np.pi/2,
            "I_ac_start_t"     : no.theta_start_t,
            "I_noise_std"      : no.noise_sigma,
            "V_clamp"          : no.Vclamp,
            "C_Mg"             : no.C_Mg,
            "rat_pos_x"        : [],
            "rat_pos_y"        : []}


def getINeuronParams(no):
    '''
    Return a dictionary of I neuron parameters, using `no`: neuron options
    object.
    '''
    return {
            "V_m"              : no.EL_i,
            "C_m"              : no.taum_i * no.gL_i,
            "t_ref"            : no.t_ref_i,
            "V_peak"           : no.V_peak_i,
            "V_reset"          : no.Vr_i,
            "E_L"              : no.EL_i,
            "g_L"              : no.gL_i,
            "Delta_T"          : no.deltaT_i,
            "V_th"             : no.Vt_i,
            "E_AMPA"           : no.E_AMPA,
            "E_GABA_A"         : no.E_GABA_A,
            "tau_AMPA_fall"    : no.tau_AMPA,
            "tau_NMDA_fall"    : no.tau_NMDA_fall,
            "tau_GABA_A_fall"  : no.tau_GABA_A_fall,
            "tau_AHP"          : no.ad_tau_i_mean,
            "E_AHP"            : no.EL_i,  # AHP has a role of adaptation here
            "g_AHP_max"        : no.ad_i_g_inc,
            "g_AHP_ad"         : True,
            "I_const"          : no.Iext_i_const,
            "I_ac_amp"         : no.Iext_i_theta,
            "I_ac_freq"        : no.theta_freq,
            "I_ac_phase"       : -np.pi/2,
            "I_ac_start_t"     : no.theta_start_t,
            "I_noise_std"      : no.noise_sigma,
            "V_clamp"          : no.Vclamp,
            "g_NMDA_fraction"  : no.NMDA_amount,
            "C_Mg"             : no.C_Mg,
            "rat_pos_x"        : [],
            "rat_pos_y"        : []}

