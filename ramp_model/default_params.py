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
        "taum_e"              :   9.3,          # ms
        "taum_e_spread"       :   3.1,          # ms
        "EL_e"                :   -68.5,        # mV
        "EL_e_spread"         :   2.0,          # mV
        "Vt_e"                :   -50,          # mV
        "Vr_e"                :   -68.5,        # mV
        "gL_e"                :   None,
        "deltaT_e"            :   0.4,          # mV
        "Eahp_e"              :   -80,          # mV
        "g_ahp_e"             :   5,            # nS
        "tau_ahp_e"           :   20,           # ms
        
        "taum_i"              :   10,           # ms
        "taum_i_spread"       :   4,            # ms
        "EL_i"                :   -60,          # mV
        "EL_i_spread"         :   20,           # mV
        "Vt_i"                :   -45,          # mV
        "Vr_i"                :   -60,          # mV
        "Rm_i"                :   44e6, 
        "ad_tau_i_mean"       :   7.5,          # ms
        "ad_tau_i_std"        :   0.5,          # ms # Unused in the simulation for now
        "ad_i_g_inc"          :   22.73,        # nS
        "deltaT_i"            :   0.4,          # mV
        
        "tau_AMPA"            :   1,            # ms
        "g_AMPA_total"        :   35,           # nS
        "g_AMPA_std"          :   0.6,          # nS
        "tau_GABA_rise"       :   1,            # ms
        "tau_GABA_fall"       :   5,            # ms
        "g_GABA_total"        :   40,           # nS
    
        "g_extraGABA_total"   :   40,           # nS
        "extraGABA_density"   :   0.4,
        
        "Vrev_AMPA"           :   0,            # mV
        "Vrev_GABA"           :   -75,          # mV
        
        "noise_sigma"         :   2,            # mV
        "sigma_init_cond"     :   10,           # mV
        
        "refrac_abs"          :   0.1,          # ms
        
        "time"                :   0.6,          # ms
        "sim_dt"              :   0.1,          # ms
        "spike_detect_th"     :   40,           # mV
        "Vclamp"              :   -50,          # mV
        
        "ntrials"             :   1,
        
        "output_dir"          :   _defaultOutputDir,
        "job_num"             :   0}


