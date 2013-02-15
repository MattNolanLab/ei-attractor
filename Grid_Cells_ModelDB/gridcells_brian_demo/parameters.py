#
#   parameters.py
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

__all__ = ['defaultParameters', 'defaultParameters']


y_dim       = np.sqrt(3)/2.


class DefaultParameters(object):
    def __init__(s):
        s.time                  =   1e3          # ms
        s.output_dir            =   "./"
        
        s.Ne                    =   68
        s.Ni                    =   34
        
        s.gridSep               = 70             # cm
        
        s.pAMPA_mu              =   y_dim/2.0
        s.pAMPA_sigma           =   0.5/6
        s.pGABA_mu              =   y_dim/2.0
        s.pGABA_sigma           =   0.5/6
        s.AMPA_gaussian         =   0            # bool
        s.prefDirC_e            =   0            # neurons: not normalized, i.e. 4
        s.prefDirC_i            =   0            # neurons: not normalized
        s.arenaSize             =   180.0        # cm
        
        s.NMDA_amount           =   2.0          # %
        
        s.Iext_e_const          =   500.0        # pA
        s.Iext_i_const          =   300.0        # pA
        s.Iext_start            =   100.0        # pA
        s.Iext_start_dur        =   200.0        # ms
        s.Iext_start_size       =   1./3         # fraction
        
        
        s.taum_e                =   9.3          # ms
        s.taum_e_spread         =   2.0          # ms
        s.EL_e                  =   -68.5        # mV
        s.EL_e_spread           =   5.0          # mV
        s.Vt_e                  =   -50          # mV
        s.Vr_e                  =   -68.5        # mV
        s.gL_e                  =   22.73        # nS
        s.deltaT_e              =   0.4          # mV
        s.E_AHP_e               =   -80          # mV
        s.tau_AHP_e             =   20           # ms
        s.g_AHP_e_max           =   5.0          # nS
        s.t_ref_e               =   0.1          # ms
        s.V_peak_e              =   -40          # mV
        
        s.taum_i                =   10           # ms
        s.taum_i_spread         =   2.0          # ms
        s.EL_i                  =   -60          # mV
        s.EL_i_spread           =   2.0          # mV
        s.Vt_i                  =   -45          # mV
        s.Vr_i                  =   -60          # mV
        s.gL_i                  =   22.73        # nS
        s.t_ref_i               =   0.1          # ms
        s.deltaT_i              =   0.4          # mV
        s.ad_tau_i_mean         =   7.5          # ms
        s.ad_tau_i_std          =   0.5          # ms Unused in the simulation for now
        s.ad_i_g_inc            =   22.73        # nS
        s.V_peak_i              =  -35           # mV
        
        s.tau_AMPA              =   1            # ms
        s.tau_NMDA_fall         =    100         # ms only a single exponential used here
        s.g_AMPA_total          =   1200.0        # nS
        s.tau_GABA_A_rise       =   0.1          # ms
        s.tau_GABA_A_fall       =   5            # ms
        s.g_GABA_total          =   400.0        # nS
        
        s.E_AMPA                =   0            # mV
        s.E_GABA_A              =   -75          # mV
        
        s.noise_sigma           =   2            # mV
        s.sigma_init_cond       =   10           # mV
        
        s.sim_dt                =   0.1          # ms
        s.Vclamp                =   -50          # mV
        
        
        
        s.stateMonDur           =  20e3          # ms


defaultParameters = DefaultParameters()

