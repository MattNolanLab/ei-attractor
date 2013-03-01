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
from place_cells        import *

import random
import nest


class PosInputs(object):
    def __init__(self, pos_x, pos_y, pos_dt):
        self.pos_x = pos_x
        self.pos_y = pos_y
        self.pos_dt = pos_dt



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

        self._ratVelocitiesLoaded = False
        self._placeCellsLoaded = False

        self.PC = []

        nest.Install('gridcellsmodule')
        nest.ResetKernel()
        nest.SetKernelStatus({"resolution" : self.no.sim_dt, "print_time": False})
        nest.SetKernelStatus({"local_num_threads" : self.no.nthreads})

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
                "g_AHP_ad"         : False,
                "I_const"          : self.no.Iext_e_const,
                "I_ac_amp"         : self.no.Iext_e_theta,
                "I_ac_freq"        : self.no.theta_freq,
                "I_ac_phase"       : -np.pi/2,
                "I_ac_start_t"     : self.no.theta_start_t,
                "I_noise_std"      : self.no.noise_sigma,
                "V_clamp"          : self.no.Vclamp,
                "rat_pos_x"        : [],
                "rat_pos_y"        : []}

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
                "E_AHP"            : self.no.EL_i,  # AHP has a role of adaptation here 
                "g_AHP_max"        : self.no.ad_i_g_inc,
                "g_AHP_ad"         : True,
                "I_const"          : self.no.Iext_i_const,
                "I_ac_amp"         : self.no.Iext_i_theta,
                "I_ac_freq"        : self.no.theta_freq,
                "I_ac_phase"       : -np.pi/2,
                "I_ac_start_t"     : self.no.theta_start_t,
                "I_noise_std"      : self.no.noise_sigma,
                "V_clamp"          : self.no.Vclamp,
                "g_NMDA_fraction"  : self.no.NMDA_amount,
                "rat_pos_x"        : [],
                "rat_pos_y"        : []}


        self.B_GABA = 1.0   # Must be here for compatibility with brian code

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
        nest.CopyModel('static_synapse', 'PC_AMPA',
                params={'receptor_type' : self.e_receptors['AMPA']})


        # Connect E-->I and I-->E
        self._centerSurroundConnection(self.no.AMPA_gaussian, self.no.pAMPA_mu,
                self.no.pAMPA_sigma, self.no.pGABA_mu, self.no.pGABA_sigma)

        self._initStates()
        self._initCellularProperties()


    def simulate(self, time, printTime=True):
        '''Run the simulation'''
        nest.SetKernelStatus({"print_time": printTime})
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



    def getGenericSpikeDetector(self, gids, label=""):
        '''
        NEST specific function to get a spike detector that monitors a
        population of neurons with global id set to gids

        gids    A list of global ids of the neurons to monitor
        '''
        mon = nest.Create('spike_detector')
        nest.SetStatus(mon, {
            "label"     : label,
            'withtime'  : True,
            'withgid'   : True})
        nest.ConvergentConnect(gids, mon)
        return mon



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

    def setStartCurrent(self, force_pos=None):
        # This will either need to instantiate the start current model or set
        # parameters of the external current model
        pass



    def _loadRatVelocities(self):
        '''
        Load rat velocities (in this case positions only)
        '''
        if self._ratVelocitiesLoaded:
            return

        self.ratData    = loadmat(self.no.ratVelFName)
        self.rat_dt     = self.ratData['dt'][0][0]*1e3      # units: ms

        self.rat_pos_x  = self.ratData['pos_x'].ravel()
        self.rat_pos_y  = self.ratData['pos_y'].ravel()

        # Map velocities to currents: we use the slope of bump speed vs. rat speed and
        # inter-peak grid field distance to remap
        # Bump current slope must be estimated
        self.velC = self.Ne_x / self.no.gridSep / self.no.bumpCurrentSlope

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
            raise NotImplementedException("NestGridCellNetwork.setVelocityCurrentInput_e.prefDirs_mask")


        # Load velocities into nest: they are all shared among all iaf_gridcells
        # nodes so only one neuron needs setting the actual values
        npos = int(self.no.time / self.rat_dt)
        nest.SetStatus([self.E_pop[0]], {
            "rat_pos_x" :  self.rat_pos_x[0:npos],
            "rat_pos_y" :  self.rat_pos_y[0:npos],
            "rat_pos_dt":  self.rat_dt}) # s --> ms

        nest.SetStatus(self.E_pop, "pref_dir_x", self.prefDirs_e[:, 0]);
        nest.SetStatus(self.E_pop, "pref_dir_y", self.prefDirs_e[:, 1]);
        nest.SetStatus(self.E_pop, "velC"      , self.velC);


    def setConstantVelocityCurrent_e(self, vel, start_t=None, end_t=None):
        '''
        Set the model so that there is only a constant velocity current input.
        '''
        if start_t is not None:
            raise Exception("Const velocity start time cannot be overridden in this model!")

        start_t = self.no.theta_start_t

        if end_t is None:
            end_t = self.no.time


        self.rat_dt = 20.0 # ms
        nVel = int((end_t - start_t) / self.rat_dt)
        self.rat_pos_x = np.cumsum(np.array([vel[0]] * nVel)) * (self.rat_dt*1e-3)
        self.rat_pos_y = np.cumsum(np.array([vel[1]] * nVel)) * (self.rat_dt*1e-3)

        self._ratVelocitiesLoaded = True # Force this velocities, not the animal

        # Load velocities into nest: they are all shared among all iaf_gridcells
        # nodes so only one neuron needs setting the actual values
        nest.SetStatus([self.E_pop[0]], {
            "rat_pos_x" :  self.rat_pos_x,
            "rat_pos_y" :  self.rat_pos_y,
            "rat_pos_dt":  self.rat_dt}) # s --> ms

        print self.rat_pos_x
        print self.rat_pos_y

        # Map velocities to currents: we use the slope of bump speed vs. rat speed and
        # inter-peak grid field distance to remap
        # Bump current slope must be estimated
        self.velC = self.Ne_x / self.no.gridSep / self.no.bumpCurrentSlope

        nest.SetStatus(self.E_pop, "pref_dir_x", self.prefDirs_e[:, 0]);
        nest.SetStatus(self.E_pop, "pref_dir_y", self.prefDirs_e[:, 1]);
        nest.SetStatus(self.E_pop, "velC"      , self.velC);

        self.setPlaceCells(start=0.0, end=self.no.theta_start_t,
                posIn=PosInputs([.0], [.0], self.rat_dt))


    def setPlaceCells(self, start=None, end=None, posIn=None):
        '''
        Generate place cells and connect them to grid cells. The wiring is
        fixed, and there is no plasticity.
        '''
        if start is None:
            start = self.no.theta_start_t
        if end is None:
            end = self.no.time
        if (posIn is None):
            self._loadRatVelocities()
            posIn = PosInputs(self.rat_pos_x, self.rat_pos_y, self.rat_dt)

        if (self.no.N_place_cells != 0):
            print "Setting up place cells"

            boxSize = [self.no.arenaSize, self.no.arenaSize]
            N_pc_size = int(np.sqrt(self.no.N_place_cells))
            N = [N_pc_size, N_pc_size]
            self.N_pc_created = N_pc_size**2
            self.PCHelper = UniformBoxPlaceCells(boxSize, N,
                    self.no.pc_max_rate, self.no.pc_field_std, random=False)

            self.PC = nest.Create('place_cell_generator', self.N_pc_created,
                    params={'rate'       : self.no.pc_max_rate,
                            'field_size' : self.no.pc_field_std,
                            'start'      : start,
                            'stop'       : end})
            nest.SetStatus(self.PC, 'ctr_x', self.PCHelper.centers[:, 0])
            nest.SetStatus(self.PC, 'ctr_y', self.PCHelper.centers[:, 1])

            #print "ctr_x", self.PCHelper.centers[:, 0]
            #print "ctr_y", self.PCHelper.centers[:, 1]

            npos = int(self.no.time / self.rat_dt)
            nest.SetStatus([self.PC[0]], params={
                'rat_pos_x' : posIn.pos_x[0:npos],
                'rat_pos_y' : posIn.pos_y[0:npos],
                'rat_pos_dt': posIn.pos_dt})

            test_x = nest.GetStatus([self.PC[0]], 'rat_pos_x')
            test_y = nest.GetStatus([self.PC[0]], 'rat_pos_y')
            #print test_x, test_y


            # Connections
            # Here we extract connections from the PlaceCellInput class that was
            # originaly used as a current input generator for place cell
            # resetting mechanism. The output of this class perfectly matches
            # how divergent connections from a single place cell should be
            # mapped onto the twisted torus grid cell sheet

            # how divergent the connections are, 3sigma rule --> division by 6.
            connStdDev          = self.no.gridSep / 2. / 6.
            pc_weight_threshold = 0.1

            pc_input = PlaceCellInput(self.Ne_x, self.Ne_y, self.no.arenaSize,
                    self.no.gridSep, [.0, .0], fieldSigma=connStdDev)
            ctr_x = nest.GetStatus(self.PC, 'ctr_x')
            ctr_y = nest.GetStatus(self.PC, 'ctr_y')
            #print ctr_x, ctr_y
            for pc_id in xrange(self.N_pc_created):
                w = pc_input.getSheetInput(ctr_x[pc_id], ctr_y[pc_id]).flatten()
                gt_th = w > pc_weight_threshold
                post = np.array(self.E_pop)[gt_th]
                w    = w[gt_th]
                #print post, w
                nest.DivergentConnect(
                        [self.PC[pc_id]],
                        list(post),
                        weight=list(w * self.no.pc_conn_weight),
                        delay=[self.no.delay] * len(w),
                        model='PC_AMPA')



        else:
            print "Warning: trying to set up place cells with N_place_cells == 0"

        self._placeCellsLoaded = True


    ############################################################################ 
    #                                   Other
    ############################################################################ 
    def getRatData(self):
        return self.ratData
