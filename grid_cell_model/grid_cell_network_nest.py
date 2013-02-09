#
#   grid_cell_network_nest.py
#
#   Nest-specific implementation of the grid cell model
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

from scipy              import linspace
from numpy.random       import rand, randn

from common             import *
from grid_cell_network  import *
from place_input        import *

import random
import nest


class NestGridCellNetwork(GridCellNetwork):
    def uniformDistrib(self, mean, spread, N):
        return mean - spread/2.0 * rand(N)

    def _initStates(self):
        # Initialize membrane potential randomly
        nest.SetStatus(self.E_pop, 'V_m', self.no.EL_e + (self.no.Vt_e-self.no.EL_e) * np.random.rand(len(self.E_pop)))
        nest.SetStatus(self.I_pop, 'V_m', self.no.EL_i + (self.no.Vt_i-self.no.EL_i) * np.random.rand(len(self.I_pop)))


    def _initCellularProperties(self):
        EL_e    = self.uniformDistrib(self.no.EL_e,   self.no.EL_e_spread,   len(self.E_pop))
        taum_e  = self.uniformDistrib(self.no.taum_e, self.no.taum_e_spread, len(self.E_pop))
        EL_i    = self.uniformDistrib(self.no.EL_i,   self.no.EL_i_spread,   len(self.I_pop))
        taum_i  = self.uniformDistrib(self.no.taum_i, self.no.taum_i_spread, len(self.I_pop))
        nest.SetStatus(self.E_pop, 'E_L',  EL_e)
        nest.SetStatus(self.E_pop, 'C_m',  taum_e * self.no.gL_e)
        nest.SetStatus(self.I_pop, 'E_L',  EL_i)
        nest.SetStatus(self.I_pop, 'C_m', taum_i * self.no.gL_i)

        #self.I_pop.tau_ad   = (self.no.ad_tau_i_mean + self.no.ad_tau_i_std * np.random.randn(len(self.I_pop.tau_ad))) * ms

    def _initClocks(self):
        for clk in self._clocks:
            clk.reinit()

    def reinit(self):
        self.net.reinit(states=False)
        self._initStates()
        self._initClocks()



    def __init__(self, neuronOpts, simulationOpts):
        GridCellNetwork.__init__(self, neuronOpts, simulationOpts)

        self.spikeMon_e = None
        self.spikeMon_i = None
        self.stateMon_e = None
        self.stateMon_i = None

        nest.Install('gridcellsmodule')
        nest.ResetKernel()
        nest.SetKernelStatus({"resolution" : self.no.sim_dt, "print_time": True})
        nest.SetKernelStatus({"local_num_threads" : 4})

        self.e_neuron_params = {
                "V_m"              : self.no.EL_e,
                "C_m"              : self.no.taum_e * self.no.gL_e,
                "t_ref"            : self.no.t_ref_e,
                "V_peak"           : self.no.V_peak_e,
                "V_reset"          : self.no.Vr_e,
                "E_L"              : self.no.EL_e,
                "g_L"              : self.no.gL_e,
                "Delta_T"          : self.no.deltaT_e,
                "V_th"             : self.no.Vt_e,
                "E_AMPA"           : self.no.E_AMPA,
                "E_GABA_A"         : self.no.E_GABA_A,
                "tau_AMPA_fall"    : self.no.tau_AMPA,
                "tau_NMDA_fall"    : self.no.tau_NMDA_fall,
                "tau_GABA_A_fall"  : self.no.tau_GABA_A_fall,
                "tau_AHP"          : self.no.tau_AHP_e,
                "E_AHP"            : self.no.E_AHP_e,
                "g_AHP_max"        : self.no.g_AHP_e_max,
                "I_e"              : self.no.Iext_e_const,
                "V_clamp"          : self.no.Vclamp}

        self.i_neuron_params = {
                "V_m"              : self.no.EL_i,
                "C_m"              : self.no.taum_i * self.no.gL_i,
                "t_ref"            : self.no.t_ref_i,
                "V_peak"           : self.no.V_peak_i,
                "V_reset"          : self.no.Vr_i,
                "E_L"              : self.no.EL_i,
                "g_L"              : self.no.gL_i,
                "Delta_T"          : self.no.deltaT_i,
                "V_th"             : self.no.Vt_i,
                "E_AMPA"           : self.no.E_AMPA,
                "E_GABA_A"         : self.no.E_GABA_A,
                "tau_AMPA_fall"    : self.no.tau_AMPA,
                "tau_NMDA_fall"    : self.no.tau_NMDA_fall,
                "tau_GABA_A_fall"  : self.no.tau_GABA_A_fall,
                "tau_AHP"          : self.no.ad_tau_i_mean,
                "E_AHP"            : self.no.EL_i,  # !!! Here there is no adaptaion
                "g_AHP_max"        : self.no.ad_i_g_inc,
                "I_e"              : self.no.Iext_i_const,
                "V_clamp"          : self.no.Vclamp,
                "g_NMDA_fraction"  : self.no.NMDA_amount}


        #tau1_GABA = self.no.tau_GABA_A_fall
        #tau2_GABA = self.no.tau_GABA_A_rise * self.no.tau_GABA_A_fall / \
        #        (self.no.tau_GABA_A_rise + self.no.tau_GABA_A_fall);
        #self.B_GABA = 1/((tau2_GABA/tau1_GABA)**(self.no.tau_GABA_A_rise/tau1_GABA) - 
        #        (tau2_GABA/tau1_GABA)**(self.no.tau_GABA_A_rise/tau2_GABA))
        self.B_GABA = 1.0

        self.e_model_name = "iaf_gridcells"
        self.i_model_name = "iaf_gridcells"
        self.e_receptors = nest.GetDefaults(self.e_model_name)['receptor_types']
        self.i_receptors = nest.GetDefaults(self.i_model_name)['receptor_types']
        self.E_pop = nest.Create(self.e_model_name, self.net_Ne,
                params=self.e_neuron_params)
        self.I_pop = nest.Create(self.i_model_name, self.net_Ni, params =
                self.i_neuron_params)

        nest.CopyModel('static_synapse', 'I_AMPA_NMDA',
                params={'receptor_type' : self.i_receptors['AMPA_NMDA']})
        nest.CopyModel('static_synapse', 'E_GABA_A',
                params={'receptor_type' : self.e_receptors['GABA_A']})


        # Connect E-->I and I-->E
        self._centerSurroundConnection(self.no.AMPA_gaussian, self.no.pAMPA_mu,
                self.no.pAMPA_sigma, self.no.pGABA_mu, self.no.pGABA_sigma)

        # Noise generators
        self.noise_gen_e = nest.Create('noise_generator', self.net_Ne,
                params={'mean' : 0.0, 'std' : self.no.noise_sigma})
        self.noise_gen_i = nest.Create('noise_generator', self.net_Ni,
                params={'mean' : 0.0, 'std' : self.no.noise_sigma})

        nest.Connect(self.noise_gen_e, self.E_pop)
        nest.Connect(self.noise_gen_i, self.I_pop)

        self._initStates()
        self._initCellularProperties()


    def simulate(self, time):
        '''Run the simulation'''
        nest.Simulate(time)


    def getSpikeDetector(self, type, N_ids=None):
        '''
        Get a spike detector that records from neurons given N_ids and from the
        population type given by type
        '''
        if (type == "E"):
            if self.spikeMon_e is not None:
                return self.spikeMon_e
            else:
                self.spikeMon_e = nest.Create('spike_detector')
                nest.SetStatus(self.spikeMon_e, {
                    "label"     : "E spikes",
                    'withtime'  : True,
                    'withgid'   : True})
                nest.ConvergentConnect(self.E_pop, self.spikeMon_e)
                return self.spikeMon_e
        elif (type == "I"):
            if (self.spikeMon_i is not None):
                return self.spikeMon_i
            else:
                self.spikeMon_i = nest.Create('spike_detector')
                nest.SetStatus(self.spikeMon_i, {
                    "label"     : "I spikes",
                    'withtime'  : True,
                    'withgid'   : True})
                nest.ConvergentConnect(self.I_pop, self.spikeMon_i)
                return self.spikeMon_i
        else:
            raise ValueError("Unsupported type of spike detector: " + type)


    def getStateMonitor(self, type, N_ids, params):
        '''
        Return a state monitor for a given population (type) and relative
        indexes of neurons (N_ids), with parameters given by params
        '''
        if (len(N_ids) == 0):
            raise ValueError("State monitor needs to record from at least one neuron")

        N = len(N_ids)

        if (type == "E"):
            if self.stateMon_e is None:
                self.stateMon_e = nest.Create('multimeter', N, params=params)
                nest.Connect(self.stateMon_e, self.E_pop[0] + np.array(N_ids))
            return self.stateMon_e
        elif (type == "I"):
            if (self.stateMon_i is None):
                self.stateMon_i = nest.Create('multimeter', N, params=params)
                nest.Connect(self.stateMon_i, self.I_pop[0] + np.array(N_ids))
            return self.stateMon_i


    def uniformInhibition(self):
        g_GABA_mean = self.no.g_uni_GABA_total / self.net_Ni / self.no.uni_GABA_density * nS

        self.extraGABA_conn1 = Connection(self.I_pop, self.E_pop, 'gi1')
        self.extraGABA_conn2 = Connection(self.I_pop, self.E_pop, 'gi2')



    #def _divergentConnectEE(self, pre, post, weights):

    def _divergentConnectEI(self, pre, post, weights):
        post_global = list(self.I_pop[0] + np.array(post))
        nest.DivergentConnect([self.E_pop[0] + pre], post_global, model='I_AMPA_NMDA',
                weight=list(weights), delay=[self.no.delay]*len(weights))

    def _divergentConnectIE(self, pre, post, weights):
        post_global = list(self.E_pop[0] + np.array(post))
        nest.DivergentConnect([self.I_pop[0] + pre], post_global, model='E_GABA_A',
                weight=list(weights), delay=[self.no.delay]*len(weights))



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


    def setStartCurrent(self, force_pos=None):
        # This will either need to instantiate the start current model or set
        # parameters of the external current model
        pass



    def setThetaCurrentStimulation(self):
        pass


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




    def setPlaceCurrentInput(self):
        pass




    ############################################################################ 
    #                                   Other
    ############################################################################ 
    def getRatData(self):
        return self.ratData
