#
#   default_params.py
#
#   Default neuron and network parameters
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

import numpy as np

__all__ = ['defaultParameters']

_defaultOutputDir = "output/"

y_dim       = np.sqrt(3)/2.

defaultParameters = {
        "Ne"                    :   68,
        "Ni"                    :   34,
        "delay"                 :   1.0,          # ms

        "rat_dt"                :   20.0,         # ms
        "ratVelFName"           : './rat_trajectory_lowpass.mat',

        "gridSep"               : 70,             # cm

        "pAMPA_mu"              :   y_dim/2.0,
        "pAMPA_sigma"           :   0.5/6,
        "pGABA_mu"              :   y_dim/2.0,
        "pGABA_sigma"           :   0.5/6,
        "AMPA_gaussian"         :   0,            # bool
        "prefDirC_e"            :   4.0,
        "prefDirC_i"            :   0.0,
        "arenaSize"             :   180.0,        # cm
        "Iplace"                :   250.0,        # pA
        "placeDur"              :    100,         # ms
        "placeSigma"            :      7,         # cm (CHECK)

        "NMDA_amount"           :   2.0,          # %

        "Iext_e_const"          :   300.0,        # pA
        "Iext_i_const"          :   200.0,        # pA
        "Iext_start"            :   300.0,        # pA
        "Iext_start_dur"        :   100.0,        # ms
        "Iext_e_theta"          :   375.0,        # pA
        "Iext_i_theta"          :    25.0,        # pA
        "sigmaIextGaussian"     :     0.5,        # Normalised to <0, 1>
        "shuffleIextGaussian"   :       0,        # 0 or 1
        "theta_start_t"         :   0.5e3,        # ms
        "theta_freq"            :   8,            # Hz
        "theta_noise_sigma"     :   0.0,            # pA
        "theta_ph_jit_mean_e"   :   0.0,          # rad
        "theta_ph_jit_spread_e" :   0.0,            # rad
        "theta_ph_jit_mean_i"   :   0.0,          # rad
        "theta_ph_jit_spread_i" :   0.0,            # rad

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
        "g_AMPA_total"          :   1400,         # nS
        "g_uni_AMPA_total"      :    0,           # nS
        "uni_AMPA_density"      :   0.001,        # fraction
        "tau_GABA_A_rise"       :   0.1,          # ms
        "tau_GABA_A_fall"       :   5,            # ms
        "g_GABA_total"          :   2160,         # nS
        "g_uni_GABA_total"      :   28,           # nS
        "uni_GABA_density"      :   0.4,

        "condAddPercSynapses_e" :   0.0,          # %
        "condAdd_e"             :   0.0,          # nS
    
        "E_AMPA"                :   0,            # mV
        "E_GABA_A"              :   -75,          # mV
        
        "noise_sigma"           :   2,            # mV
        "sigma_init_cond"       :   10,           # mV
        
        "sim_dt"                :   0.1,          # ms
        "Vclamp"                :   -50,          # mV
        
        "ntrials"               :   1,
        
        "output_dir"            :   _defaultOutputDir,
        "stateMonDur"           :  20e3,          # ms

        "stim_spread"           :  1.0}


