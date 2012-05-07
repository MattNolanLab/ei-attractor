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


import numpy    as np
import logging  as lg

from numpy.random       import rand, randn
from brian              import *
from scipy              import linspace

from common             import *
from grid_cell_network  import *
from place_input        import *

class IextManager(object):
    '''
    Simple interface to manage external currents.
    '''
    def __init__(self):
        self._Iext_e = None
        self._Iext_i = None

    def set_Iext(self, Iext_e, Iext_i):
        self._Iext_e = np.array(Iext_e)
        self._Iext_i = np.array(Iext_i)

    def add_Iext(self, Iext_e, Iext_i):
        self._Iext_e += Iext_e
        self._Iext_i += Iext_i



class BrianGridCellNetwork(GridCellNetwork):
    def uniformDistrib(self, mean, spread, N):
        return mean - spread/2.0 * rand(N)

    def _initStates(self):
        # Initialize membrane potential randomly
        self.E_pop.vm       = (self.no.EL_e + (self.no.Vt_e-self.no.EL_e) * rand(len(self.E_pop))) * mV
        self.I_pop.vm       = (self.no.EL_i + (self.no.Vt_i-self.no.EL_i) * rand(len(self.I_pop))) * mV

        self.E_pop.EL       = self.uniformDistrib(self.no.EL_e,   self.no.EL_e_spread,   len(self.E_pop)) * mV
        self.E_pop.taum     = self.uniformDistrib(self.no.taum_e, self.no.taum_e_spread, len(self.E_pop)) * ms
        self.I_pop.EL       = self.uniformDistrib(self.no.EL_i,   self.no.EL_i_spread,   len(self.I_pop)) * mV
        self.I_pop.taum     = self.uniformDistrib(self.no.taum_i, self.no.taum_i_spread, len(self.I_pop)) * ms

        self.I_pop.tau_ad   = (self.no.ad_tau_i_mean + self.no.ad_tau_i_std * randn(len(self.I_pop.tau_ad))) * ms

    def reinit(self):
        self.net.reinit(states=True)
        self._initStates()

    def __init__(self, neuronOpts, simulationOpts):
        GridCellNetwork.__init__(self, neuronOpts, simulationOpts)

        self._simulationClock = Clock(dt = self.no.sim_dt * msecond)

        # Place cell input settings
        self._gridsep = self.no.arenaSize/self.no.gridsPerArena
        self._gridCenter = [0, 0]
        self._place_current = 250 * pA
        self._pc = PlaceCellInput(self.Ne_x, self.Ne_y, self.no.arenaSize, self._gridsep, self._gridCenter)

        tau1_GABA = self.no.tau_GABA_A_fall
        tau2_GABA = self.no.tau_GABA_A_rise * self.no.tau_GABA_A_fall / \
                (self.no.tau_GABA_A_rise + self.no.tau_GABA_A_fall);
        self.B_GABA = 1/((tau2_GABA/tau1_GABA)**(self.no.tau_GABA_A_rise/tau1_GABA) - 
                (tau2_GABA/tau1_GABA)**(self.no.tau_GABA_A_rise/tau2_GABA))

        self.eqs_e = Equations('''
            dvm/dt      = 1/C*Im + (noise_sigma*xi/taum_mean**.5)                : volt
            Ispike      = gL*deltaT*exp((vm-Vt)/deltaT)                          : amp
            Im          = gL*(EL-vm) + g_ahp*(Eahp - vm) + Ispike + Isyn + Iext  : amp
            Isyn        = B_GABA*(gi1 - gi2)*(Esyn - vm)                         : amp
            Iclamp      = -(gi1 - gi2)*(Esyn - Vclamp)                           : amp
            dgi1/dt     = -gi1/syn_tau1                                          : siemens
            dgi2/dt     = -gi2/syn_tau2                                          : siemens
            dg_ahp/dt   = -g_ahp/tau_ahp                                         : siemens
            Iext        = Iext_const + Iext_theta + Iext_vel + Iext_start        : amp
            Iext_const                                                           : amp
            Iext_theta                                                           : amp
            Iext_vel                                                             : amp
            Iext_start                                                           : amp
            EL                                                                   : volt
            taum                                                                 : second
            ''',
            C           = self.no.taum_e * self.no.gL_e * pF,
            gL          = self.no.gL_e * nS,
            noise_sigma = self.no.noise_sigma * mV,
            deltaT      = self.no.deltaT_e * mV,
            Vt          = self.no.Vt_e * mV,
            Esyn        = self.no.E_GABA_A * mV,
            Vclamp      = self.no.Vclamp * mV,
            syn_tau1    = tau1_GABA * ms,
            syn_tau2    = tau2_GABA * ms,
            B_GABA      = self.B_GABA,
            taum_mean   = self.no.taum_e * ms,
            tau_ahp     = self.no.tau_AHP_e * ms,
            Eahp        = self.no.E_AHP_e * mV)


        self.eqs_i = Equations('''
            dvm/dt      = 1/C*Im + (noise_sigma*xi/taum_mean**.5)           : volt
            Ispike      = gL*deltaT*exp((vm-Vt)/deltaT)                     : amp
            Im          = gL*(EL-vm)*(1+g_ad/gL) + Ispike + Isyn + Iext     : amp
            Isyn        = ge*(Esyn - vm) + gNMDA*(Esyn - vm)                : amp
            Iclamp      = -(ge*(Esyn - Vclamp) + gNMDA*(Esyn - Vclamp))     : amp
            dge/dt      = -ge/syn_tau                                       : siemens
            dg_ad/dt    = -g_ad/tau_ad                                      : siemens
            dgNMDA/dt   = -gNMDA/(100*msecond)                              : siemens
            tau_ad                                                          : second
            Iext        = Iext_const + Iext_theta + Iext_vel                : amp
            Iext_const                                                      : amp
            Iext_theta                                                      : amp
            Iext_vel                                                        : amp
            Iext_start                                                      : amp
            EL                                                              : volt
            taum                                                            : second
            ''',
            C           = self.no.taum_i * self.no.gL_i * pF,
            gL          = self.no.gL_i * nS,
            noise_sigma = self.no.noise_sigma * mV,
            deltaT      = self.no.deltaT_i * mV,
            Vt          = self.no.Vt_i * mV,
            Esyn        = self.no.E_AMPA * mV,
            Vclamp      = self.no.Vclamp * mV,
            syn_tau     = self.no.tau_AMPA * ms,
            taum_mean   = self.no.taum_i * ms)


        # Other constants
        g_AHP_e = self.no.g_AHP_e_max * nS
        Vr_e    = self.no.Vr_e * mV    


        # Setup neuron groups and connections
        self.E_pop = NeuronGroup(
                N = self.net_Ne,
                model=self.eqs_e,
                threshold=self.no.V_peak_e * mV,
                reset="vm=Vr_e; g_ahp=g_AHP_e",
                refractory=self.no.t_ref_e * msecond,
                clock=self._simulationClock)

        self.I_pop = NeuronGroup(
                N = self.net_Ni,
                model=self.eqs_i,
                threshold=self.no.V_peak_i * mV,
                reset=self.no.Vr_i * mV,
                refractory=self.no.t_ref_i * msecond,
                clock=self._simulationClock)

        self.net = Network(self.E_pop, self.I_pop)

        # Setup adaptation connections: neuron on itself
        if self.no.ad_i_g_inc != 0.0:
            self.adaptConn_i = IdentityConnection(self.I_pop, self.I_pop, 'g_ad',
                    weight=self.no.ad_i_g_inc*nS)
            self.net.add(self.adaptConn_i)

        # Connect E-->I and I-->E
        self.AMPA_conn = Connection(self.E_pop, self.I_pop, 'ge',
            structure='dense')
        self.NMDA_conn = Connection(self.E_pop, self.I_pop, 'gNMDA',
            structure='dense')
        self.GABA_conn1 = Connection(self.I_pop, self.E_pop, 'gi1')
        self.GABA_conn2 = Connection(self.I_pop, self.E_pop, 'gi2')

        # Weight matrices which are used in _divergentConnectXY() functions
        self._E_W = np.asarray(self.AMPA_conn.W)

        self._centerSurroundConnection(self.no.pAMPA_mu, self.no.pAMPA_sigma, self.no.pGABA_sigma)

        # Now simply copy AMPA --> NMDA and GABA_conn1 --> GABA_conn2
        self.NMDA_conn.connect(self.E_pop, self.I_pop, self.AMPA_conn.W * .01 * self.no.NMDA_amount)
        self.GABA_conn2.connect(self.I_pop, self.E_pop, self.GABA_conn1.W)

        self.net.add(self.AMPA_conn, self.NMDA_conn, self.GABA_conn1, self.GABA_conn2)

        self._initStates()


    def uniformInhibition(self):
        g_GABA_mean = self.no.g_uni_GABA_total / self.net_Ni / self.no.uni_GABA_density * nS

        self.extraGABA_conn1 = Connection(self.I_pop, self.E_pop, 'gi1')
        self.extraGABA_conn2 = Connection(self.I_pop, self.E_pop, 'gi2')

        self.extraGABA_conn1.connect_random(self.I_pop, self.E_pop, self.no.uni_GABA_density,
                weight=self.B_GABA*g_GABA_mean)
        self.extraGABA_conn2.connect(self.I_pop, self.E_pop, self.extraGABA_conn1.W) 

        self.net.add(self.extraGABA_conn1, self.extraGABA_conn2)


    def _divergentConnectEI(self, pre, post, weights):
        self._E_W[pre, post] = weights * nS


    def _divergentConnectIE(self, pre, post, weights):
        '''The inhibitory matrix is sparse!'''
        self.GABA_conn1.W.rows[pre] = post
        self.GABA_conn1.W.data[pre] = weights * nS



    ############################################################################ 
    #                     External sources definitions
    ############################################################################ 
    def setConstantCurrent(self):
        self.E_pop.Iext_const = self.no.Iext_e_const * pA
        self.I_pop.Iext_const = self.no.Iext_i_const * pA


    def setStartCurrent(self):
        self._startCurrentClock = Clock(dt=50*ms)

        @network_operation(self._startCurrentClock)
        def startCurrentFun():
            if self._simulationClock.t >= 0*msecond and self._simulationClock.t < self.no.Iext_start_dur*msecond:
                self.E_pop.Iext_start = self._pc.getSheetInput(0.0, 0.0).ravel() * self.no.Iext_start * pA
                print "Bump initialisation..."
            else:
                self.E_pop.Iext_start = 0.0

        self.net.add(startCurrentFun)


    def setThetaCurrentStimulation(self):
        '''
        Create procedures for theta stimulation:
            theat_start_t   Start time of theta stimulation (ms)
        '''
        self.stim_omega = 2*np.pi*self.no.theta_freq*Hz
        self.stim_e_A  = self.no.Iext_e_theta/2 * pA
        self.stim_i_A  = self.no.Iext_i_theta/2 * pA

        @network_operation(self._simulationClock)
        def thetaStimulationFun():
            global place_flag
            global place_I
            if self._simulationClock.t >= self.no.theta_start_t*msecond and \
                    self._simulationClock.t < self.no.time*msecond:
                ph = self.stim_omega*self._simulationClock.t
                self.E_pop.Iext_theta = self.stim_e_A + self.stim_e_A*np.sin(ph - np.pi/2)
                self.I_pop.Iext_theta = self.stim_i_A + self.stim_i_A*np.sin(ph - np.pi/2)
            else:
                self.E_pop.Iext_theta = 0.0
                self.I_pop.Iext_theta = 0.0

        self.net.add(thetaStimulationFun)


    ############################################################################ 
    #                                   Other
    ############################################################################ 
    def _getSimulationClock(self):
        return self._simulationClock


