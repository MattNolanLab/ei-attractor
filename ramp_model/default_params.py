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

__all__ = ['defaultParameters']

_defaultOutputDir = "output/"

defaultParameters = {
        "Ne"                  :   500,
        "Ni"                  :   125,
        "Iext_e"              :   500.0,        # pA
        "Iext_i"              :   300.0,        # pA
        "AMPA_density"        :   0.4,
        "GABA_density"        :   0.4,
        "delay"               :   1.0,          # ms

        "taum_e"              :   9.3,          # ms
        "taum_e_spread"       :   3.1,          # ms
        "EL_e"                :   -68.5,        # mV
        "EL_e_spread"         :   2.0,          # mV
        "Vt_e"                :   -50,          # mV
        "Vr_e"                :   -68.5,        # mV
        "gL_e"                :   22.73,        # nS
        "deltaT_e"            :   0.5,          # mV
        "E_AHP_e"             :   -80,          # mV
        "tau_AHP_e"           :   20,           # ms
        "g_AHP_e_max"         :   30.0,         # nS
        "t_ref_e"             :   0.1,          # ms
        "V_peak_e"            :   40,           # mV
        
        "taum_i"              :   10,           # ms
        "taum_i_spread"       :   4,            # ms
        "EL_i"                :   -60,          # mV
        "EL_i_spread"         :   20,           # mV
        "Vt_i"                :   -45,          # mV
        "Vr_i"                :   -60,          # mV
        "gL_i"                :   22.73,        # nS
        "t_ref_i"             :   0.1,          # ms
        "deltaT_i"            :   0.4,          # mV
        "E_AHP_i"             :   -80,          # mV
        "tau_AHP_i"           :   5,           # ms
        "g_AHP_i_max"         :   30.0,         # nS
        "V_peak_i"            :   40,           # mV
        
        "tau_AMPA"            :   1,            # ms
        "g_AMPA_total"        :   70,           # nS
        "g_AMPA_std"          :   0.6,          # nS
        "tau_GABA_A_rise"     :   1,            # ms
        "tau_GABA_A_fall"     :   5,            # ms
        "g_GABA_total"        :   40,           # nS
    
        "E_AMPA"              :   0,            # mV
        "E_GABA_A"            :   -75,          # mV
        
        "noise_sigma"         :   2,            # mV
        "sigma_init_cond"     :   10,           # mV
        
        "sim_dt"              :   0.1,          # ms
        "Vclamp"              :   -50,          # mV
        
        "ntrials"             :   1,
        
        "output_dir"          :   _defaultOutputDir,
        "job_num"             :   0,

        "stim_spread"         :  1.0}


