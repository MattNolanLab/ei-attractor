'''Default model parameters; E-->E structured connections.'''
import numpy as np

__all__ = ['defaultParameters']

_defaultOutputDir = "output/"

y_dim       = np.sqrt(3)/2.

defaultParameters = {
        "Ne"                    :   34,
        "Ni"                    :   34,
        "delay"                 :   0.1,          # ms
        "nthreads"              :   1,
        "printTime"             :   0,            # This is boolean

        "ratVelFName"           : '../../../data/hafting_et_al_2005/rat_trajectory_lowpass.mat',
        #"ratVelFName"           : '../../../data/hafting_et_al_2005/rat_trajectory_original.mat',

        "gridSep"               :   60,             # cm

        "EI_flat"               :   1,            # bool
        "IE_flat"               :   1,            # bool
        "use_EE"                :   1,            # bool
        "AMPA_gaussian"         :   0,            # bool
        "pEE_sigma"             :   .05,
        "pAMPA_mu"              :   y_dim/2.0,
        "pAMPA_sigma"           :   0.5/6,
        "pGABA_mu"              :   y_dim/2.0,
        "pGABA_sigma"           :   0.5/6,
        "prefDirC_e"            :   4.0,
        "prefDirC_ee"           :   0.0,
        "prefDirC_i"            :   0.0,
        "arenaSize"             :   180.0,        # cm

        "NMDA_amount"           :   2.0,          # %
        "C_Mg"                  :      .0,        # mM; def is no Vm dependence

        "Iext_e_const"          :   300.0,        # pA
        "Iext_i_const"          :   200.0,        # pA
        "Iext_e_theta"          :   375.0,        # pA
        "Iext_i_theta"          :    25.0,        # pA

        "theta_start_t"         :   0.5e3,        # ms
        "theta_freq"            :   8,            # Hz

        "taum_e"                :   9.3,          # ms
        "taum_e_spread"         :   0.31,         # ms
        "EL_e"                  :   -68.5,        # mV
        "EL_e_spread"           :   0.20,         # mV
        "Vt_e"                  :   -50,          # mV
        "Vr_e"                  :   -68.5,        # mV
        "gL_e"                  :   22.73,        # nS
        "deltaT_e"              :   0.4,          # mV
        "E_AHP_e"               :   -80,          # mV
        "tau_AHP_e"             :   20,           # ms
        "g_AHP_e_max"           :   5.0,          # nS
        "t_ref_e"               :   0.1,          # ms
        "V_peak_e"              :   -40,          # mV

        "taum_i"                :   10,           # ms
        "taum_i_spread"         :   0,            # ms
        "EL_i"                  :   -60,          # mV
        "EL_i_spread"           :   0,            # mV
        "Vt_i"                  :   -45,          # mV
        "Vr_i"                  :   -60,          # mV
        "gL_i"                  :   22.73,        # nS
        "t_ref_i"               :   0.1,          # ms
        "deltaT_i"              :   0.4,          # mV
        "ad_tau_i_mean"         :   7.5,          # ms
        "ad_tau_i_std"          :   0.5,          # ms, Unused in the simulation for now
        "ad_i_g_inc"            :   22.73,        # nS
        "V_peak_i"              :  -35,           # mV

        "tau_AMPA"              :   1,            # ms
        "tau_NMDA_fall"         :    100,         # ms, only a single exponential used here
        "g_EE_total"            :   4000.,        # nS
        "g_AMPA_total"          :   1400,         # nS
        "g_uni_AMPA_total"      :    0,           # nS
        "uni_AMPA_density"      :   0.001,        # fraction
        "tau_GABA_A_rise"       :   0.1,          # ms
        "tau_GABA_A_fall"       :   5,            # ms
        "g_GABA_total"          :   2160,         # nS
        "g_uni_GABA_frac"       :   0.013,        # fraction of g_GABA_total
        "uni_GABA_density"      :   0.4,

        "g_EI_uni_density"      :     .1,         # Probability
        "g_IE_uni_density"      :     .1,         # Probability

        "E_AMPA"                :   0,            # mV
        "E_GABA_A"              :   -75,          # mV

        "N_place_cells"         :   30,           # sqrt(total PC number)
        "pc_max_rate"           :   50.0,         # Hz
        "pc_conn_weight"        :   0.5,          # nS
        "pc_field_std"          :   20.0,         # cm
        "bumpCurrentSlope"      :   0.53,         # neurons/s/pA, !! this will depend on prefDirC !!
        "pc_start_max_rate"     :   100.0,        # Hz
        "pc_start_conn_weight"  :   5.0,          # nS

        "noise_sigma"           :   150.0,        # pA
        "gammaNSample"          :   25,           # No. of neurons

        "sim_dt"                :   0.1,          # ms
        "Vclamp"                :   -50,          # mV

        "ntrials"               :   1,

        "output_dir"            :   _defaultOutputDir,
        "stateMonDur"           :  20e3,          # ms
}


