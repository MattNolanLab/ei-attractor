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

_defaultOutputDir = "output_local/"

y_dim = np.sqrt(3)/2.
pAMPA_mu = y_dim/2.
pAMPA_sigma = 0.5/6
pGABA_sigma = 0.5/6

defaultParameters = {
        "Ne"                  :   68,
        "Ni"                  :   34,
        "delay"               :   1.0,          # ms

        "rat_dt"              :   20.0,         # ms
        "Ivel"                :   50.0,         # pA
        "Ivel_max"            :   50.0,         # pA
        "ratVelFName"         : '../../../../data/hafting_et_al_2005/rat_trajectory_lowpass.mat',
        #"ratVelFName"         : '../../../../data/hafting_et_al_2005/rat_data_original.mat',

        "pAMPA_mu"            :   y_dim/2.0,
        "pAMPA_sigma"         :   0.5/6,
        "pGABA_sigma"         :   0.5/6,
        "prefDirC"            :   4.0,
        "arenaSize"           :   180.0,        # cm
        "gridsPerArena"       :   2.5,

        "NMDA_amount"         :   2.0,          # %

        "Iext_e_const"        :   425.0,        # pA
        "Iext_i_const"        :   200.0,        # pA
        "Iext_start"          :   300.0,        # pA
        "Iext_start_dur"      :   100.0,        # ms
        "theta_start_t"       :   500.0,        # ms

        "taum_e"              :   9.3,          # ms
        "taum_e_spread"       :   0.31,         # ms
        "EL_e"                :   -68.5,        # mV
        "EL_e_spread"         :   0.20,         # mV
        "Vt_e"                :   -50,          # mV
        "Vr_e"                :   -68.5,        # mV
        "gL_e"                :   22.73,        # nS
        "deltaT_e"            :   0.4,          # mV
        "E_AHP_e"             :   -80,          # mV
        "tau_AHP_e"           :   20,           # ms
        "g_AHP_e_max"         :   5.0,          # nS
        "t_ref_e"             :   0.1,          # ms
        "V_peak_e"            :   -40,          # mV
        
        "taum_i"              :   10,           # ms
        "taum_i_spread"       :   0,            # ms
        "EL_i"                :   -60,          # mV
        "EL_i_spread"         :   0,            # mV
        "Vt_i"                :   -45,          # mV
        "Vr_i"                :   -60,          # mV
        "gL_i"                :   22.73,        # nS
        "t_ref_i"             :   0.1,          # ms
        "deltaT_i"            :   0.4,          # mV
        "ad_tau_i_mean"       :   7.5,          # ms
        "ad_tau_i_std"        :   0.5,          # ms, Unused in the simulation for now
        "ad_i_g_inc"          :   22.73,        # nS
        "V_peak_i"            :  -35,           # mV
        
        "tau_AMPA"            :   1,            # ms
        "g_AMPA_total"        :   1400,         # nS
        "tau_GABA_A_rise"     :   0.1,          # ms
        "tau_GABA_A_fall"     :   5,            # ms
        "g_GABA_total"        :   2160,         # nS
        "g_uni_GABA_total"    :   28,           # nS
        "uni_GABA_density"    :   0.4,
    
        "E_AMPA"              :   0,            # mV
        "E_GABA_A"            :   -75,          # mV
        
        "noise_sigma"         :   2,            # mV
        "sigma_init_cond"     :   10,           # mV
        
        "sim_dt"              :   0.1,          # ms
        "Vclamp"              :   -50,          # mV
        
        "ntrials"             :   1,
        
        "output_dir"          :   _defaultOutputDir,

        "stim_spread"         :  1.0}


