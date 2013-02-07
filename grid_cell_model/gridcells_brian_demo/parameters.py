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

_defaultOutputDir = "output/"

y_dim       = np.sqrt(3)/2.


class Parameters(object):
    pass


defaultParameters = Parameters()

defaultParameters.Ne                    =   68
defaultParameters.Ni                    =   34
defaultParameters.delay                 =   1.0          # ms

defaultParameters.rat_dt                =   20.0         # ms
defaultParameters.ratVelFName           = '../../../../data/hafting_et_al_2005/rat_trajectory_lowpass.mat'

defaultParameters.gridSep               = 70             # cm

defaultParameters.pAMPA_mu              =   y_dim/2.0
defaultParameters.pAMPA_sigma           =   0.5/6
defaultParameters.pGABA_mu              =   y_dim/2.0
defaultParameters.pGABA_sigma           =   0.5/6
defaultParameters.AMPA_gaussian         =   0            # bool
defaultParameters.prefDirC_e            =   4.0
defaultParameters.prefDirC_i            =   0.0
defaultParameters.arenaSize             =   180.0        # cm
defaultParameters.Iplace                =   250.0        # pA
defaultParameters.placeDur              =    100         # ms
defaultParameters.placeSigma            =      7         # cm (CHECK)

defaultParameters.NMDA_amount           =   2.0          # %

defaultParameters.Iext_e_const          =   300.0        # pA
defaultParameters.Iext_i_const          =   200.0        # pA
defaultParameters.Iext_start            =   300.0        # pA
defaultParameters.Iext_start_dur        =   100.0        # ms
defaultParameters.Iext_e_theta          =   375.0        # pA
defaultParameters.Iext_i_theta          =    25.0        # pA
defaultParameters.theta_start_t         =   0.5e3        # ms
defaultParameters.theta_freq            =   8            # Hz
defaultParameters.theta_noise_sigma     =   0.0          # pA

defaultParameters.taum_e                =   9.3          # ms
defaultParameters.taum_e_spread         =   0.31         # ms
defaultParameters.EL_e                  =   -68.5        # mV
defaultParameters.EL_e_spread           =   0.20         # mV
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
defaultParameters.taum_i_spread         =   0            # ms
defaultParameters.EL_i                  =   -60          # mV
defaultParameters.EL_i_spread           =   0            # mV
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
defaultParameters.g_AMPA_total          =   1400         # nS
defaultParameters.tau_GABA_A_rise       =   0.1          # ms
defaultParameters.tau_GABA_A_fall       =   5            # ms
defaultParameters.g_GABA_total          =   2160         # nS
defaultParameters.g_uni_GABA_total      =   28           # nS
defaultParameters.uni_GABA_density      =   0.4


defaultParameters.E_AMPA                =   0            # mV
defaultParameters.E_GABA_A              =   -75          # mV

defaultParameters.noise_sigma           =   2            # mV
defaultParameters.sigma_init_cond       =   10           # mV

defaultParameters.sim_dt                =   0.1          # ms
defaultParameters.Vclamp                =   -50          # mV


defaultParameters.output_dir            =   _defaultOutputDir
defaultParameters.stateMonDur           =  20e3          # ms



