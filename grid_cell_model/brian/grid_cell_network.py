#
#   grid_cell_network.py
#
#   Brian-specific implementation of the grid cell model
#
#     Copyright (C) 2012  Lukas Solanka <l.solanka@sms.ed.ac.uk>
#     
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#     
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#     
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
#


from common import *

import numpy as np
import logging as lg
from numpy.random import rand, randn
from brian import *
from scipy import linspace


class BrianGridCellNetwork(GridCellNetwork):
    def uniformDistrib(self, mean, spread, N):
        return mean - spread/2.0 * rand(N)

    def _initStates(self):
        # Initialize membrane potential randomly
        self.E_pop.vm       = self.no.EL_e + (self.no.Vt_e-self.no.EL_e) * rand(len(self.E_pop)) * mV
        self.I_pop.vm       = self.no.EL_i + (self.no.Vt_i-self.no.EL_i) * rand(len(self.I_pop)) * mV

        self.E_pop.EL       = uniformDistrib(self.no.EL_e,   self.no.EL_e_spread,   len(self.E_pop)) * mV
        self.E_pop.taum     = uniformDistrib(self.no.taum_e, self.no.taum_e_spread, len(self.E_pop)) * ms
        self.I_pop.EL       = uniformDistrib(self.no.EL_i,   self.no.EL_i_spread,   len(self.I_pop)) * mV
        self.I_pop.taum     = uniformDistrib(self.no.taum_i, self.no.taum_i_spread, len(self.I_pop)) * ms

        self.I_pop.tau_ad   = (self.no.ad_tau_i_mean + self.no.ad_tau_i_std * randn(len(self.I_pop.tau_ad))) * ms

        self._setBackgroundInput(self.no.Iext_e * pA, self.no.Iext_i * pA)

    def reinit(self):
        self.net.reinit(states=True)
        self._initStates()

    def __init__(self, neuronOpts, simulationOpts):
        GridCellNetwork.__init__(self, neuronOpts, simulationOpts)

        no = neuronOpts
        so = simulationOpts


        tau1_GABA = no.tau_GABA_fall
        tau2_GABA = no.tau_GABA_rise * no.tau_GABA_fall / (no.tau_GABA_rise + no.tau_GABA_fall);
        self.B_GABA = 1/((tau2_GABA/tau1_GABA)**(tau_GABA_rise/tau1_GABA) - 
                (tau2_GABA/tau1_GABA)**(tau_GABA_rise/tau2_GABA))

        self.eqs_e = Equations('''
            dvm/dt      = 1/C*Im + (noise_sigma*xi/taum_mean**.5)                : volt
            Ispike      = gL*deltaT*exp((vm-Vt)/deltaT)                          : amp
            Im          = gL*(EL-vm) + g_ahp*(Eahp - vm) + Ispike + Isyn + Iext  : amp
            Isyn        = B_GABA*(gi1 - gi2)*(Esyn - vm)                         : amp
            Iclamp      = -(gi1 - gi2)*(Esyn - Vclamp)                           : amp
            dgi1/dt     = -gi1/syn_tau1                                          : siemens
            dgi2/dt     = -gi2/syn_tau2                                          : siemens
            dg_ahp/dt   = -g_ahp/tau_ahp                                         : siemens
            Iext                                                                 : amp
            EL                                                                   : volt
            taum                                                                 : second
            ''',
            C           = no.taum_e * no.gL_e * pF,
            gL          = no.gL_e * nS,
            noise_sigma = no.noise_sigma * mV,
            deltaT      = no.deltaT_e * mV,
            Vt          = no.Vt_e * mV,
            Esyn        = no.Vrev_GABA * mV,
            Vclamp      = no.Vclamp * mV,
            syn_tau1    = tau1_GABA * ms,
            syn_tau2    = tau2_GABA * ms,
            B_GABA      = self.B_GABA,
            taum_mean   = no.taum_e * ms,
            tau_ahp     = no.tau_ahp_e * ms,
            Eahp        = no.E_AHP_e * mV)


        self.eqs_i = Equations('''
            dvm/dt      = 1/C*Im + (noise_sigma*xi/taum_mean**.5)        : volt
            Ispike      = gL*deltaT*exp((vm-Vt)/deltaT)                  : amp
            Im          = gL*(EL-vm)*(1+g_ad/gL) + Ispike + Isyn + Iext  : amp
            Isyn        = ge*(Esyn - vm) + gNMDA*(Esyn - vm)             : amp
            Iclamp      = -(ge*(Esyn - Vclamp) + gNMDA*(Esyn - Vclamp))  : amp
            dge/dt      = -ge/syn_tau                                    : siemens
            dg_ad/dt    = -g_ad/tau_ad                                   : siemens
            dgNMDA/dt   = -gNMDA/(100*msecond)                           : siemens
            tau_ad                                                       : second
            Iext                                                         : amp
            EL                                                           : volt
            taum                                                         : second
            ''',
            C           = no.taum_i * no.gL_i * pF,
            gL          = no.gL_i * nS,
            noise_sigma = no.noise_sigma * mV,
            deltaT      = no.deltaT_i * mV,
            Vt          = no.Vt_i * mV,
            Esyn        = no.Vrev_AMPA * mV,
            Vclamp      = no.Vclamp * mV,
            syn_tau     = no.tau_AMPA * ms,
            taum_mean   = no.taum_i * ms)


        # Other constants
        refrac_abs      = no.refrac_abs * ms
        spike_detect_th = no.spike_detect_th * mV
        g_AHP_e         = no.g_AHP_e * nS


        # Setup neuron groups and connections
        self.E_pop = NeuronGroup(
                N = self.net_Ne,
                model=self.eqs_e,
                threshold=spike_detect_th,
                reset="vm=Vr_e; g_ahp=g_AHP_e",
                refractory=refrac_abs,
                clock=clk)

        self.I_pop = NeuronGroup(
                N = self.net_Ni,
                model=self.eqs_i,
                threshold=spike_detect_th,
                reset=Vr_i,
                refractory=refrac_abs,
                clock=clk)

        self.net = Network(self.E_pop, self.I_pop)

        # Setup adaptation connections: neuron on itself
        if no.ad_i_g_inc != 0.0:
            self.adaptConn_i = IdentityConnection(self.I_pop, self.I_pop, 'g_ad',
                    weight=o.ad_i_g_inc*nS)
            self.net.add(self.adaptConn_i)

        # Connect E-->I and I-->E
        self.AMPA_conn = Connection(self.E_pop, self.I_pop, 'ge',
            structure='dense')
        self.NMDA_conn = Connection(self.E_pop, self.I_pop, 'gNMDA',
            structure='dense')
        self.GABA_conn1 = Connection(self.I_pop, self.E_pop, 'gi1')
        self.GABA_conn2 = Connection(self.I_pop, self.E_pop, 'gi2')

        self._centerSurroundConnection(no.pAMPA_mu, no.pAMPA_sigma, no.pGABA_sigma)

        # Now simply copy AMPA --> NMDA and GABA_conn1 --> GABA_conn2
        self.NMDA_conn.connect(self.E_pop, self.I_pop, self.AMPA_conn.W * .01 * no.NMDA_amount)
        self.GABA_conn2.connect(self.I_pop, self.E_pop, self.GABA_conn1.W)

        self.net.add(self.AMPA_conn, self.NMDA_conn, self.GABA_conn1, self.GABA_conn2)

        self._initStates()


    def _setBackgroundInput(self, Iext_e, Iext_i):
        self.E_pop.Iext = linspace(Iext_e, Iext_e, len(self.E_pop))
        self.I_pop.Iext = linspace(Iext_i, Iext_i, len(self.I_pop))

