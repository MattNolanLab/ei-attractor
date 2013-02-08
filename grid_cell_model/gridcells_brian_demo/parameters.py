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

__all__ = ['defaultParameters']


y_dim       = np.sqrt(3)/2.


class Parameters(object):
    pass


defaultParameters = Parameters()

defaultParameters.time                  =   1e3          # ms

defaultParameters.Ne                    =   20
defaultParameters.Ni                    =   20

defaultParameters.gridSep               = 70             # cm

defaultParameters.pAMPA_mu              =   y_dim/2.0
defaultParameters.pAMPA_sigma           =   0.5/6
defaultParameters.pGABA_mu              =   y_dim/2.0
defaultParameters.pGABA_sigma           =   0.5/6
defaultParameters.AMPA_gaussian         =   0            # bool
defaultParameters.prefDirC_e            =   0.0
defaultParameters.prefDirC_i            =   0.0
defaultParameters.arenaSize             =   180.0        # cm

defaultParameters.NMDA_amount           =   2.0          # %

defaultParameters.Iext_e_const          =   500.0        # pA
defaultParameters.Iext_i_const          =   300.0        # pA
defaultParameters.Iext_start            =   100.0        # pA
defaultParameters.Iext_start_dur        =   200.0        # ms
defaultParameters.Iext_start_size       =   1./3         # fraction


defaultParameters.taum_e                =   9.3          # ms
defaultParameters.taum_e_spread         =   2.0          # ms
defaultParameters.EL_e                  =   -68.5        # mV
defaultParameters.EL_e_spread           =   5.0          # mV
defaultParameters.Vt_e                  =   -50          # mV
defaultParameters.Vr_e                  =   -68.5        # mV
defaultParameters.gL_e                  =   22.73        # nS
defaultParameters.deltaT_e              =   0.4          # mV
defaultParameters.E_AHP_e               =   -80          # mV
defaultParameters.tau_AHP_e             =   20           # ms
defaultParameters.g_AHP_e_max           =   5.0          # nS
defaultParameters.t_ref_e               =   0.1          # ms
defaultParameters.V_peak_e              =   -40          # mV

defaultParameters.taum_i                =   10           # ms
defaultParameters.taum_i_spread         =   2.0          # ms
defaultParameters.EL_i                  =   -60          # mV
defaultParameters.EL_i_spread           =   2.0          # mV
defaultParameters.Vt_i                  =   -45          # mV
defaultParameters.Vr_i                  =   -60          # mV
defaultParameters.gL_i                  =   22.73        # nS
defaultParameters.t_ref_i               =   0.1          # ms
defaultParameters.deltaT_i              =   0.4          # mV
defaultParameters.ad_tau_i_mean         =   7.5          # ms
defaultParameters.ad_tau_i_std          =   0.5          # ms Unused in the simulation for now
defaultParameters.ad_i_g_inc            =   22.73        # nS
defaultParameters.V_peak_i              =  -35           # mV

defaultParameters.tau_AMPA              =   1            # ms
defaultParameters.tau_NMDA_fall         =    100         # ms only a single exponential used here
defaultParameters.g_AMPA_total          =   1200.0        # nS
defaultParameters.tau_GABA_A_rise       =   0.1          # ms
defaultParameters.tau_GABA_A_fall       =   5            # ms
defaultParameters.g_GABA_total          =   400.0        # nS

defaultParameters.E_AMPA                =   0            # mV
defaultParameters.E_GABA_A              =   -75          # mV

defaultParameters.noise_sigma           =   2            # mV
defaultParameters.sigma_init_cond       =   10           # mV

defaultParameters.sim_dt                =   0.1          # ms
defaultParameters.Vclamp                =   -50          # mV


defaultParameters.output_dir            =   "output_local/"

defaultParameters.stateMonDur           =  20e3          # ms



