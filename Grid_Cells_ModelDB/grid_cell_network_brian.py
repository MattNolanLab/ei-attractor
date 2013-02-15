#
#   grid_cell_network_brian.py
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

#from numpy.random       import rand, randn
from brian              import *
from scipy              import linspace

from common             import *
from grid_cell_network  import *
from place_input        import *

import random


class BrianGridCellNetwork(GridCellNetwork):
    def uniformDistrib(self, mean, spread, N):
        return mean - spread/2.0 * rand(N)

    def _initStates(self):
        # Initialize membrane potential randomly
        self.E_pop.vm       = (self.no.EL_e + (self.no.Vt_e-self.no.EL_e) * np.random.rand(len(self.E_pop))) * mV
        self.E_pop.gi1      = 0.0
        self.E_pop.gi2      = 0.0
        self.E_pop.g_ahp    = 0.0

        self.I_pop.vm       = (self.no.EL_i + (self.no.Vt_i-self.no.EL_i) * np.random.rand(len(self.I_pop))) * mV
        self.I_pop.ge       = 0.0
        self.I_pop.g_ad     = 0.0
        self.I_pop.gNMDA    = 0.0


    def _initCellularProperties(self):
        self.E_pop.EL       = self.uniformDistrib(self.no.EL_e,   self.no.EL_e_spread,   len(self.E_pop)) * mV
        self.E_pop.taum     = self.uniformDistrib(self.no.taum_e, self.no.taum_e_spread, len(self.E_pop)) * ms
        self.I_pop.EL       = self.uniformDistrib(self.no.EL_i,   self.no.EL_i_spread,   len(self.I_pop)) * mV
        self.I_pop.taum     = self.uniformDistrib(self.no.taum_i, self.no.taum_i_spread, len(self.I_pop)) * ms

        self.I_pop.tau_ad   = (self.no.ad_tau_i_mean + self.no.ad_tau_i_std * np.random.randn(len(self.I_pop.tau_ad))) * ms

    def _initClocks(self):
        for clk in self._clocks:
            clk.reinit()

    def reinit(self):
        self.net.reinit(states=False)
        self._initStates()
        self._initClocks()


    def setGaussianIextEnv(self):
        if self.no.sigmaIextGaussian is None:
            return 1.0, 1.0

        # E cells, distance normalised
        X_e, Y_e = np.meshgrid(range(self.Ne_x), range(self.Ne_y))
        X_e = (X_e - self.Ne_x/2.0)
        Y_e = (Y_e - self.Ne_y/2.0)
        d2_e = (X_e**2 + Y_e**2) / self.Ne_x**2 
        env_e = (np.exp(-d2_e / 2.0 / self.no.sigmaIextGaussian**2)).ravel()

        # I cells
        X_i, Y_i = np.meshgrid(range(self.Ni_x), range(self.Ni_y))
        X_i = (X_i - self.Ni_x/2.0)
        Y_i = (Y_i - self.Ni_y/2.0)
        d2_i = (X_i**2 + Y_i**2) / self.Ni_x**2
        env_i = (np.exp(-d2_i / 2.0 / self.no.sigmaIextGaussian**2)).ravel()

        if self.no.shuffleIextGaussian == 1:
            random.shuffle(env_e)
            random.shuffle(env_i)

        return env_e, env_i


    def __init__(self, neuronOpts, simulationOpts):
        GridCellNetwork.__init__(self, neuronOpts, simulationOpts)

        self._ratVelocitiesLoaded = False

        self._clocks = []
        self._simulationClock = Clock(dt = self.no.sim_dt * msecond)
        self._ratClock = Clock(self.no.rat_dt*msecond)
        self._clocks.append(self._simulationClock)
        self._clocks.append(self._ratClock)

        # Place cell input settings
        self._gridsep = self.no.gridSep
        self._gridCenter = [0, 0]
        self._pc = PlaceCellInput(self.Ne_x, self.Ne_y, self.no.arenaSize,
                self._gridsep, self._gridCenter, self.no.placeSigma)
        # Has place cell input been set up?
        self._placeCellInputOn = False

        # Set up the Gaussian envelope
        self.gaussianEnv_e, self.gaussianEnv_i = self.setGaussianIextEnv()


        tau1_GABA = self.no.tau_GABA_A_fall
        tau2_GABA = self.no.tau_GABA_A_rise * self.no.tau_GABA_A_fall / \
                (self.no.tau_GABA_A_rise + self.no.tau_GABA_A_fall);
        self.B_GABA = 1/((tau2_GABA/tau1_GABA)**(self.no.tau_GABA_A_rise/tau1_GABA) - 
                (tau2_GABA/tau1_GABA)**(self.no.tau_GABA_A_rise/tau2_GABA))

        self.eqs_e = Equations('''
            dvm/dt      = 1/C*Im + (noise_sigma*xi/taum_mean**.5)                      : volt
            Ispike      = gL*deltaT*exp((vm-Vt)/deltaT)                                : amp
            Im          = gL*(EL-vm) + g_ahp*(Eahp - vm) + Ispike + Isyn + Iext        : amp
            Isyn        = B_GABA*(gi1 - gi2)*(Esyn_i - vm) + ge*(Esyn_e - vm) + gNMDA*(Esyn_e - vm) : amp
            Iclamp      = -(B_GABA*(gi1 - gi2)*(Esyn_i - Vclamp) + ge*(Esyn_e - Vclamp) + gNMDA*(Esyn_e - Vclamp))                        : amp
            Iclamp_all  = -( gL*(EL-Vclamp) + gL*deltaT*exp((Vclamp-Vt)/deltaT) + Iext ) + Iclamp : amp
            dge/dt      = -ge/syn_tau_e                                                : siemens
            dgNMDA/dt   = -gNMDA/tau_NMDA_fall                                         : siemens
            dgi1/dt     = -gi1/syn_tau1                                                : siemens
            dgi2/dt     = -gi2/syn_tau2                                                : siemens
            dg_ahp/dt   = -g_ahp/tau_ahp                                               : siemens
            Iext        = Iext_const + Iext_theta + Iext_vel + Iext_start + Iext_place : amp
            Iext_const                                                                 : amp
            Iext_theta                                                                 : amp
            Iext_vel                                                                   : amp
            Iext_start                                                                 : amp
            Iext_place                                                                 : amp
            EL                                                                         : volt
            taum                                                                       : second
            ''',
            C             = self.no.taum_e * self.no.gL_e * pF,
            gL            = self.no.gL_e * nS,
            noise_sigma   = self.no.noise_sigma * mV,
            deltaT        = self.no.deltaT_e * mV,
            Vt            = self.no.Vt_e * mV,
            Esyn_i        = self.no.E_GABA_A * mV,
            Esyn_e        = self.no.E_AMPA * mV,
            Vclamp        = self.no.Vclamp * mV,
            syn_tau_e     = self.no.tau_AMPA * ms,
            tau_NMDA_fall = self.no.tau_NMDA_fall * ms,
            syn_tau1      = tau1_GABA * ms,
            syn_tau2      = tau2_GABA * ms,
            B_GABA        = self.B_GABA,
            taum_mean     = self.no.taum_e * ms,
            tau_ahp       = self.no.tau_AHP_e * ms,
            Eahp          = self.no.E_AHP_e * mV)


        self.eqs_i = Equations('''
            dvm/dt      = 1/C*Im + (noise_sigma*xi/taum_mean**.5)           : volt
            Ispike      = gL*deltaT*exp((vm-Vt)/deltaT)                     : amp
            Im          = gL*(EL-vm)*(1+g_ad/gL) + Ispike + Isyn + Iext     : amp
            Isyn        = ge*(Esyn - vm) + gNMDA*(Esyn - vm)                : amp
            Iclamp      = -(ge*(Esyn - Vclamp) + gNMDA*(Esyn - Vclamp))     : amp
            Iclamp_all  = -( gL*(EL-Vclamp) + gL*deltaT*exp((Vclamp-Vt)/deltaT) + ge*(Esyn - Vclamp) + gNMDA*(Esyn - Vclamp) + Iext )     : amp
            dge/dt      = -ge/syn_tau                                       : siemens
            dg_ad/dt    = -g_ad/tau_ad                                      : siemens
            dgNMDA/dt   = -gNMDA/tau_NMDA_fall                              : siemens
            tau_ad                                                          : second
            Iext        = Iext_const + Iext_theta + Iext_vel                : amp
            Iext_const                                                      : amp
            Iext_theta                                                      : amp
            Iext_vel                                                        : amp
            Iext_start                                                      : amp
            EL                                                              : volt
            taum                                                            : second
            ''',
            C             = self.no.taum_i * self.no.gL_i * pF,
            gL            = self.no.gL_i * nS,
            noise_sigma   = self.no.noise_sigma * mV,
            deltaT        = self.no.deltaT_i * mV,
            Vt            = self.no.Vt_i * mV,
            Esyn          = self.no.E_AMPA * mV,
            Vclamp        = self.no.Vclamp * mV,
            syn_tau       = self.no.tau_AMPA * ms,
            tau_NMDA_fall = self.no.tau_NMDA_fall * ms,
            taum_mean     = self.no.taum_i * ms)


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

        self._centerSurroundConnection(self.no.AMPA_gaussian, self.no.pAMPA_mu,
                self.no.pAMPA_sigma, self.no.pGABA_mu, self.no.pGABA_sigma)

        # Now simply copy AMPA --> NMDA and GABA_conn1 --> GABA_conn2
        self.NMDA_conn.connect(self.E_pop, self.I_pop, self.AMPA_conn.W * .01 * self.no.NMDA_amount)
        self.GABA_conn2.connect(self.I_pop, self.E_pop, self.GABA_conn1.W)

        self.net.add(self.AMPA_conn, self.NMDA_conn, self.GABA_conn1, self.GABA_conn2)

        self._initStates()
        self._initCellularProperties()


    def uniformInhibition(self):
        g_GABA_mean = self.no.g_uni_GABA_total / self.net_Ni / self.no.uni_GABA_density * nS

        self.extraGABA_conn1 = Connection(self.I_pop, self.E_pop, 'gi1')
        self.extraGABA_conn2 = Connection(self.I_pop, self.E_pop, 'gi2')

        self.extraGABA_conn1.connect_random(self.I_pop, self.E_pop, self.no.uni_GABA_density,
                weight=self.B_GABA*g_GABA_mean)
        self.extraGABA_conn2.connect(self.I_pop, self.E_pop, self.extraGABA_conn1.W) 

        self.net.add(self.extraGABA_conn1, self.extraGABA_conn2)


    def uniformExcitation(self):
        # Uniform E-E connections (only E-E, this does not consider E-I connections)
        g_AMPA_mean = self.no.g_uni_AMPA_total / self.net_Ne / self.no.uni_AMPA_density * nS

        self.uniAMPA_conn = Connection(self.E_pop, self.E_pop, 'ge')
        self.uniAMPA_conn.connect_random(self.E_pop, self.E_pop, self.no.uni_AMPA_density,
                weight=g_AMPA_mean)

        self.net.add(self.uniAMPA_conn)


    #def _divergentConnectEE(self, pre, post, weights):

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
        '''
        Constant uniform excitatory drive. All the populations will receive this
        current throughout the simulation, given by:
          Iext_e_const  E population
          Iext_i_const  I population
        '''
        self.E_pop.Iext_const = self.no.Iext_e_const * pA
        self.I_pop.Iext_const = self.no.Iext_i_const * pA


    def setConstantGaussianCurrent(self):
        '''
        The constant current will be modulated by a Gaussian with spread sigma,
        in normalised units (normalised against the X size of the twisted torus).
        The Gaussian current is implemented by providing an extra multiplier into
        the external current sources of the model (where appropriate)
        '''
        self.E_pop.Iext_const = self.no.Iext_e_const * self.gaussianEnv_e * pA
        self.I_pop.Iext_const = self.no.Iext_i_const * self.gaussianEnv_i * pA


    def setStartCurrent(self, force_pos=None):
        self._startCurrentClock = Clock(dt=50*ms)
        self._clocks.append(self._startCurrentClock)

        if force_pos is not None:
            self.rat_pos_x = [force_pos[0]]
            self.rat_pos_y = [force_pos[1]]
            self.Ivel_it = 0
        # If no one has overwritten the velocity inputs, init the positions from [0.0, 0.0]
        # Otherwise they will be set by the place cell input methods
        elif not self._placeCellInputOn:
            self.rat_pos_x = [0.0]
            self.rat_pos_y = [0.0]
            self.Ivel_it = 0
        else:
            pass

        @network_operation(self._startCurrentClock)
        def startCurrentFun():
            if self._simulationClock.t < self.no.Iext_start_dur*msecond:
                #self.E_pop.Iext_start = self._pc.getSheetInput(0.0, 0.0).ravel() * self.no.Iext_start * pA
                self.E_pop.Iext_start = self._pc.getSheetInput(self.rat_pos_x[self.Ivel_it],
                        self.rat_pos_y[self.Ivel_it]).ravel() * self.gaussianEnv_e * self.no.Iext_start * pA
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
        self.theta_ph_jitter_e = self.uniformDistrib(self.no.theta_ph_jit_mean_e, self.no.theta_ph_jit_spread_e, len(self.E_pop))
        self.theta_ph_jitter_i = self.uniformDistrib(self.no.theta_ph_jit_mean_i, self.no.theta_ph_jit_spread_i, len(self.I_pop))

        @network_operation(self._simulationClock)
        def thetaStimulationFun():
            global place_flag
            global place_I
            if self._simulationClock.t < self.no.theta_start_t*msecond:
                self.E_pop.Iext_theta = 2 * self.stim_e_A
                self.I_pop.Iext_theta = 2 * self.stim_i_A
            elif self._simulationClock.t >= self.no.theta_start_t*msecond:
                ph = self.stim_omega*self._simulationClock.t
                self.E_pop.Iext_theta = self.stim_e_A + self.stim_e_A*np.sin(ph - np.pi/2 - self.theta_ph_jitter_e) + self.no.theta_noise_sigma*np.random.randn(self.net_Ne)*pA
                self.I_pop.Iext_theta = self.stim_i_A + self.stim_i_A*np.sin(ph - np.pi/2 - self.theta_ph_jitter_i) + self.no.theta_noise_sigma*np.random.randn(self.net_Ni)*pA
            else:
                self.E_pop.Iext_theta = 0.0
                self.I_pop.Iext_theta = 0.0

        self.net.add(thetaStimulationFun)

    def setGaussianThetaCurrent(self):
        '''
        Theta current external stimulation, but multiplicatively modulated
        by a Gaussian function
        '''
        self.stim_omega = 2*np.pi*self.no.theta_freq*Hz
        self.stim_e_A  = self.no.Iext_e_theta/2 * pA
        self.stim_i_A  = self.no.Iext_i_theta/2 * pA

        @network_operation(self._simulationClock)
        def thetaStimulationFun():
            global place_flag
            global place_I
            if self._simulationClock.t < self.no.theta_start_t*msecond:
                self.E_pop.Iext_theta = 2 * self.stim_e_A * self.gaussianEnv_e
                self.I_pop.Iext_theta = 2 * self.stim_i_A * self.gaussianEnv_i
            elif self._simulationClock.t >= self.no.theta_start_t*msecond:
                ph = self.stim_omega*self._simulationClock.t
                self.E_pop.Iext_theta = (self.stim_e_A + self.stim_e_A*np.sin(ph - np.pi/2)) * self.gaussianEnv_e + self.no.theta_noise_sigma*np.random.randn(self.net_Ne)*pA
                self.I_pop.Iext_theta = (self.stim_i_A + self.stim_i_A*np.sin(ph - np.pi/2)) * self.gaussianEnv_i + self.no.theta_noise_sigma*np.random.randn(self.net_Ni)*pA
            else:
                self.E_pop.Iext_theta = 0.0
                self.I_pop.Iext_theta = 0.0

        self.net.add(thetaStimulationFun)



    def _loadRatVelocities(self):
        if self._ratVelocitiesLoaded:
            return

        self.ratData    = loadmat(self.no.ratVelFName)
        self.rat_dt     = self.ratData['dt'][0][0]

        self.rat_pos_x  = self.ratData['pos_x'].ravel()
        self.rat_pos_y  = self.ratData['pos_y'].ravel()
        self.rat_vel_x  = np.diff(self.rat_pos_x)/(self.no.rat_dt * ms)
        self.rat_vel_y  = np.diff(self.rat_pos_y)/(self.no.rat_dt * ms)
        
        # Map velocities to currents: we use the slope of bump speed vs. rat speed and
        # inter-peak grid field distance to remap
        # Bump current slope must be estimated
        velC = self.Ne_x / self.no.gridSep / self.no.bumpCurrentSlope
        self.rat_Ivel_x = self.rat_vel_x * velC * pA
        self.rat_Ivel_y = self.rat_vel_y * velC * pA
        
        print "mean rat_vel_x:  " + str(np.mean(np.abs(self.rat_vel_x))) + " cm/s"
        print "mean rat_vel_y:  " + str(np.mean(np.abs(self.rat_vel_y))) + " cm/s"
        print "max  rat_vel_x:  " + str(np.max(np.abs(self.rat_vel_x))) + " cm/s"
        print "max  rat_vel_y:  " + str(np.max(np.abs(self.rat_vel_y))) + " cm/s"
        print "mean rat_Ivel_x: " + str(np.mean(np.abs(self.rat_Ivel_x))/pA) + " pA"
        print "mean rat_Ivel_y: " + str(np.mean(np.abs(self.rat_Ivel_y))/pA) + " pA"
        print "max rat_Ivel_x: "  + str(np.max(np.abs(self.rat_Ivel_x))/pA) + " pA"
        print "max rat_Ivel_y: "  + str(np.max(np.abs(self.rat_Ivel_y))/pA) + " pA"

        self.Ivel_it   = 0

        self._ratVelocitiesLoaded = True


    def setVelocityCurrentInput_e(self, prefDirs_mask=None):
        '''
        Set up movement simulation, based on preferred directions of neurons.
        prefDirs_mask can be used to manipulate velocity input strength
        for each neuron.
        '''
        self._loadRatVelocities()

        if prefDirs_mask is None:
            self._prefDirs_mask_e = np.ndarray((len(self.E_pop), 2))
            self._prefDirs_mask_e[:, :] = 1.0
        else:
            self._prefDirs_mask_e = prefDirs_mask

        @network_operation(self._ratClock)
        def velocityChange():
            if self._simulationClock.t >= self.no.theta_start_t*msecond:
                v = np.array([[self.rat_Ivel_x[self.Ivel_it], self.rat_Ivel_y[self.Ivel_it]]]).T
                self.E_pop.Iext_vel = np.dot(self.prefDirs_e*self._prefDirs_mask_e, v).ravel()
                self.Ivel_it += 1
            else:
                self.E_pop.Iext_vel = 0.0

        self.net.add(velocityChange)



    def setVelocityCurrentInput_i(self):
        self._loadRatVelocities()

        @network_operation(self._ratClock)
        def velocityChange():
            if self._simulationClock.t >= self.no.theta_start_t*msecond:
                v = np.array([[self.rat_Ivel_x[self.Ivel_it], self.rat_Ivel_y[self.Ivel_it]]]).T
                self.I_pop.Iext_vel = np.dot(self.prefDirs_i, v).ravel()
                self.Ivel_it += 1
            else:
                self.I_pop.Iext_vel = 0.0

        self.net.add(velocityChange)



    def setConstantVelocityCurrent_e(self, vel, start_t=None, end_t=None):
        if start_t == None:
            self._Iext_const_start_t = self.no.theta_start_t
        else:
            self._Iext_const_start_t = start_t
        if end_t == None:
            self._Iext_const_end_t = self.no.time
        else:
            self._Iext_const_end_t = end_t

        @network_operation(self._ratClock)
        def velocityChange():
            if self._simulationClock.t >= self._Iext_const_start_t*msecond and self._simulationClock.t <= self._Iext_const_end_t*msecond:
                v = np.array([vel]).T
                self.E_pop.Iext_vel = (np.dot(self.prefDirs_e, v)).ravel() * self.no.Ivel * pA
            else:
                self.E_pop.Iext_vel = 0.0

        self.net.add(velocityChange)


    def setConstantVelocityCurrent_i(self, vel, start_t=None, end_t=None):
        if start_t == None:
            self._Iext_const_start_t = self.no.theta_start_t
        else:
            self._Iext_const_start_t = start_t
        if end_t == None:
            self._Iext_const_end_t = self.no.time
        else:
            self._Iext_const_end_t = end_t

        @network_operation(self._ratClock)
        def velocityChange():
            if self._simulationClock.t >= self._Iext_const_start_t*msecond and self._simulationClock.t <= self._Iext_const_end_t*msecond:
                v = np.array([vel]).T
                self.I_pop.Iext_vel = (np.dot(self.prefDirs_i, v)).ravel() * self.no.Ivel * pA
            else:
                self.I_pop.Iext_vel = 0.0

        self.net.add(velocityChange)


    def setPlaceCurrentInput(self):
        self._placeClock  = Clock(self.no.placeDur*msecond)
        self._clocks.append(self._placeClock)
        self._placeCellInputOn = True

        self._loadRatVelocities()

        @network_operation(self._placeClock)
        def switchPlaceInputOn():
            t_within = np.mod(self._simulationClock.t,
                    self.no.placeT*msecond)*second
            if self._simulationClock.t >= self.no.theta_start_t*msecond and \
                t_within > 0*ms and t_within <= self.no.placeDur*msecond:
                self.E_pop.Iext_place = self._pc.getSheetInput(self.rat_pos_x[self.Ivel_it],
                        self.rat_pos_y[self.Ivel_it]).ravel() * self.no.Iplace * pA
                print "  Place cell input on; pos_x = " + \
                        str(self.rat_pos_x[self.Ivel_it]) + "; pos_y = " + \
                        str(self.rat_pos_y[self.Ivel_it]) + "; t = " + \
                        str(self._simulationClock.t)
            else:
                self.E_pop.Iext_place = 0.0
                #print "  Place cell input off; t = " + str(self._simulationClock.t)
        
        self.net.add(switchPlaceInputOn)


    def setPlaceThetaCurrentInput(self, const=False):
        '''
        Place cell input which is more over modulated by theta. The gaussian
        envelope is also multiplied by the theta modulator
        We define the frequency and phase of the theta place cell input
        thetaPlaceFreq
        thetaPlacePhase
        
        If const == True, the input points to (0, 0) and doesn't change over time
        '''
        self._placeCellInputOn = True
        self._thetaPlaceInOmega = 2*np.pi*self.no.thetaPlaceFreq*Hz

        if const == False:
            self._loadRatVelocities()
            @network_operation(self._simulationClock)
            def placeThetaCurrentStimulationVar():
                if self._simulationClock.t >= self.no.theta_start_t*msecond:
                    theta_envelope = .5 + .5 * np.sin(self._thetaPlaceInOmega * self._simulationClock.t -
                            np.pi/2. - self.no.thetaPlacePhase)
                    self.E_pop.Iext_place = self._pc.getSheetInput(self.rat_pos_x[self.Ivel_it],
                            self.rat_pos_y[self.Ivel_it]).ravel() * self.no.Iplace * \
                            theta_envelope * pA
                else:
                    self.E_pop.Iext_place = 0.0
            
            self.net.add(placeThetaCurrentStimulationVar)
        else:
            @network_operation(self._simulationClock)
            def placeThetaCurrentStimulationConst():
                if self._simulationClock.t >= self.no.theta_start_t*msecond:
                    theta_envelope = .5 + .5 * np.sin(self._thetaPlaceInOmega * self._simulationClock.t -
                            np.pi/2. - self.no.thetaPlacePhase)
                    self.E_pop.Iext_place = self._pc.getSheetInput(0, 0).ravel() * self.no.Iplace * \
                            theta_envelope * pA
                else:
                    self.E_pop.Iext_place = 0.0
            
            self.net.add(placeThetaCurrentStimulationConst)




    ############################################################################ 
    #                                   Other
    ############################################################################ 
    def _getSimulationClock(self):
        return self._simulationClock

    def getRatData(self):
        return self.ratData
