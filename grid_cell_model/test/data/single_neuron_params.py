'''Parameters for tests of single neurons'''
import numpy as np

__all__ = ['defaultParameters']

defaultParameters = {
        "Ne"                    :     1,
        "Ni"                    :     1,
        "delay"                 :   0.1,          # ms
        "nthreads"              :   1,
        "printTime"             :   0,            # This is boolean

        "NMDA_amount"           :   100.0,        # %

        "Iext_e_const"          :   300.0,        # pA
        "Iext_i_const"          :   200.0,        # pA
        "Iext_e_theta"          :   375.0,        # pA
        "Iext_i_theta"          :    25.0,        # pA

        "theta_start_t"         :   0.5e3,        # ms
        "theta_freq"            :   8.,            # Hz

        "taum_e"                :   9.3,          # ms
        "taum_e_spread"         :   0.31,         # ms
        "EL_e"                  :   -68.5,        # mV
        "EL_e_spread"           :   0.20,         # mV
        "Vt_e"                  :   -50.0,        # mV
        "Vr_e"                  :   -68.5,        # mV
        "gL_e"                  :   22.73,        # nS
        "deltaT_e"              :   0.4,          # mV
        "E_AHP_e"               :   -80.0,        # mV
        "tau_AHP_e"             :   20.0,         # ms
        "g_AHP_e_max"           :   5.0,          # nS
        "t_ref_e"               :   0.1,          # ms
        "V_peak_e"              :   -40.0,        # mV
        
        "taum_i"                :   10.,          # ms
        "taum_i_spread"         :   0.,           # ms
        "EL_i"                  :   -60.,         # mV
        "EL_i_spread"           :   0.,           # mV
        "Vt_i"                  :   -45.,         # mV
        "Vr_i"                  :   -60.,         # mV
        "gL_i"                  :   22.73,        # nS
        "t_ref_i"               :   0.1,          # ms
        "deltaT_i"              :   0.4,          # mV
        "ad_tau_i_mean"         :   7.5,          # ms
        "ad_tau_i_std"          :   0.5,          # ms, Unused in the simulation for now
        "ad_i_g_inc"            :   22.73,        # nS
        "V_peak_i"              :  -35.,          # mV
        
        "tau_AMPA"              :   1.,           # ms
        "tau_NMDA_fall"         :    100.,        # ms, only a single exponential used here
        "tau_GABA_A_rise"       :   0.1,          # ms
        "tau_GABA_A_fall"       :   5.,           # ms

        "C_Mg"                  :    0.1,         # mM

        "E_AMPA"                :   0.,           # mV
        "E_GABA_A"              :   -75.,         # mV
        
        "noise_sigma"           :   0.,           # pA            
        "sim_dt"                :   0.1,          # ms
        "Vclamp"                :   -50.,         # mV
        
        "stateMonDur"           :  20e3,          # ms
}


