#
#   gc_net_nest.py
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

from numpy.random import rand, randn

from gc_net       import GridCellNetwork
from place_input  import PlaceCellInput
from place_cells  import UniformBoxPlaceCells
from data_storage import DataStorage

import nest


nest.Install('gridcellsmodule')


class PosInputs(object):
    def __init__(self, pos_x, pos_y, pos_dt):
        self.pos_x = pos_x
        self.pos_y = pos_y
        self.pos_dt = pos_dt



class NestGridCellNetwork(GridCellNetwork):
    def __init__(self, neuronOpts, simulationOpts):
        GridCellNetwork.__init__(self, neuronOpts, simulationOpts)

        self.spikeMon_e = None
        self.spikeMon_i = None
        self.stateMon_e = None
        self.stateMon_i = None

        # Extra monitors
        self._extraSpikeMons = {}   # Extra spike monitors
        self._extraStateMons = {}   # Extra state monitors

        self._ratVelocitiesLoaded = False
        self._placeCellsLoaded = False

        self.PC = []

        self._initNESTKernel()
        self._constructNetwork()
        self._initStates()
        self._initCellularProperties()


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

    def _initNESTKernel(self):
        nest.ResetKernel()
        nest.SetKernelStatus({"resolution" : self.no.sim_dt, "print_time": False})
        nest.SetKernelStatus({"local_num_threads" : self.no.nthreads})

    #def reinit(self):
    #    self._initNESTKernel()
    #    self._initStates()
    #    self._initClocks()



    def _constructNetwork(self):
        '''Construct the E/I network'''
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


    def simulate(self, time, printTime=True):
        '''Run the simulation'''
        self.endConstruction()
        self.beginSimulation()
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



    def getGenericSpikeDetector(self, gids, label):
        '''
        NEST specific function to get a spike detector that monitors a
        population of neurons with global id set to gids.

        Parameters
        ----------
        gids : list
            A list of global ids of the neurons to monitor.  There is one
            profound limitation of this method: the gids must be a list of
            increasing integers without gaps. Otherwise the local translation
            of neuron numbers when saving the data will not work!.
        '''
        mon = nest.Create('spike_detector')
        nest.SetStatus(mon, {
            "label"     : label,
            'withtime'  : True,
            'withgid'   : True})
        nest.ConvergentConnect(gids, mon)

        self._extraSpikeMons[label] = (mon, gids[0])
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


    ## Create a user defined state monitor (multimeter)
    #
    # @param gids   Global IDs of the neurons.
    # @param params Parameter of the state monitors
    # @return Global IDs for manipulation in the nest space
    #
    def getGenericStateMonitor(self, gids, params, label):
        mon = nest.Create('multimeter', len(gids), params=params)
        nest.Connect(mon, gids)
        self._extraStateMons[label] = (mon)
        return mon



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
            raise NotImplementedError()


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

        # Load velocities into nest: they are all shared among all
        # iaf_gridcells nodes so only one neuron needs setting the actual
        # values
        nest.SetStatus([self.E_pop[0]], {
            "rat_pos_x" :  self.rat_pos_x,
            "rat_pos_y" :  self.rat_pos_y,
            "rat_pos_dt":  self.rat_dt}) # s --> ms

        print self.rat_pos_x
        print self.rat_pos_y

        # Map velocities to currents: Here the mapping is 1:1, i.e. the
        # velocity dictates the current
        self.velC = 1.

        nest.SetStatus(self.E_pop, "pref_dir_x", self.prefDirs_e[:, 0]);
        nest.SetStatus(self.E_pop, "pref_dir_y", self.prefDirs_e[:, 1]);
        nest.SetStatus(self.E_pop, "velC"      , self.velC);

        self.setPlaceCells(start=0.0, end=self.no.theta_start_t,
                posIn=PosInputs([0.], [.0], self.rat_dt))


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


    def getAttrDictionary(self):
        d = {}

        d['e_neuron_params'] = self.e_neuron_params
        d['i_neuron_params'] = self.i_neuron_params
        d['B_GABA'         ] = self.B_GABA
        d['E_pop'          ] = np.array(self.E_pop)
        d['I_pop'          ] = np.array(self.I_pop)
        d['PC'             ] = np.array(self.PC)
        d['net_Ne'         ] = self.net_Ne
        d['net_Ni'         ] = self.net_Ni
        d['rat_pos_x'      ] = self.rat_pos_x
        d['rat_pos_y'      ] = self.rat_pos_y
        d['rat_dt'         ] = self.rat_dt
        d['velC'           ] = self.velC
        d['Ne_x'           ] = self.Ne_x
        d['Ne_y'           ] = self.Ne_y
        d['Ni_x'           ] = self.Ni_x
        d['Ni_y'           ] = self.Ni_y
        d['prefDirs_e'     ] = self.prefDirs_e
        d['N_pc_created'   ] = self.N_pc_created

        return d
            



class BasicGridCellNetwork(NestGridCellNetwork):
    '''
    A grid cell network that generates the common network and creates a basic
    set of spike monitors and state monitors which are generically usable in
    most of the simulation setups.
    '''

    def getDefaultStateMonParams(self):
        return {
            'withtime' : True,
            'interval' : 10.0 * self.no.sim_dt,
            'record_from' : ['V_m', 'I_clamp_AMPA', 'I_clamp_NMDA',
                'I_clamp_GABA_A', 'I_stim']
        }


    def fillParams(self, dest, src):
        for key, value in src:
            dest[key] = value
        return dest


    def __init__(self, options, simulationOpts=None,
            nrec_spikes      = (None, None),
            stateRecord_type = 'middle-center',
            stateRecParams   = (None, None)):
        '''
        TODO
        '''
        NestGridCellNetwork.__init__(self, options, simulationOpts)


        # Spikes
        self.nrecSpikes_e = nrec_spikes[0]
        self.nrecSpikes_i = nrec_spikes[1]

        if (self.nrecSpikes_e is None):
            self.nrecSpikes_e = self.Ne_x*self.Ne_y
        if (self.nrecSpikes_i is None):
            self.nrecSpike_i = self.Ni_x*self.Ni_y

        self.spikeMon_e  = self.getSpikeDetector("E")
        self.spikeMon_i  = self.getSpikeDetector("I")


        # States
        if (stateRecord_type == 'middle-center'):
            self.state_record_e = [self.Ne_x/2 -1 , self.Ne_y/2*self.Ne_x +
                    self.Ne_x/2 - 1]
            self.state_record_i = [self.Ni_x/2 - 1, self.Ni_y/2*self.Ni_x +
                    self.Ni_x/2 - 1]
        else:
            raise ValueError("Currently stateRecordType must be 'middle-center'")
        
        self.stateMonParams_e = self.getDefaultStateMonParams()
        self.stateMonParams_i = self.getDefaultStateMonParams()

        stRecp_e = stateRecParams[0]
        stRecp_i = stateRecParams[1]
        if (stRecp_e is not None):
            self.fillParams(self.stateMonParams_e, stRecp_e);
        if (stRecp_i is not None):
            self.fillParams(self.stateMonParams_i, stRecp_i);


        self.stateMon_e  = self.getStateMonitor("E",
                self.state_record_e,
                self.stateMonParams_e)
        self.stateMon_i  = self.getStateMonitor("I",
                self.state_record_i,
                self.stateMonParams_i)

    def getMonitors(self):
        return (
            self.spikeMon_e,
            self.spikeMon_i,
            self.stateMon_e,
            self.stateMon_i
        )

  
    def _getSpikeMonData(self, mon, gidStart):
        '''
        Generate a dictionary of a spike data from the monitor ``mon``
        '''
        st = nest.GetStatus(mon)[0]
        st['events']['senders'] -= gidStart
        return st


    def getAllData(self):
        '''
        Save all the simulated data into a dictionary and return it.
        '''
        out = {}
        out['options']  = self.no._einet_optdict
        out['net_attr'] = self.getAttrDictionary()
        
        # Save spikes
        if (self.spikeMon_e is not None):
            out['spikeMon_e'] = self._getSpikeMonData(self.spikeMon_e,
                    self.E_pop[0])
        if (self.spikeMon_i is not None):
            out['spikeMon_i'] = self._getSpikeMonData(self.spikeMon_i,
                    self.I_pop[0])

        for label, vals in self._extraSpikeMons.iteritems():
            assert(label not in out.keys())
            out[label] = self._getSpikeMonData(vals[0], vals[1])


        
        #Save state variables
        out['stateMon_e']   = nest.GetStatus(self.stateMon_e)
        out['stateMon_i']   = nest.GetStatus(self.stateMon_i)
        for label, val in self._extraStateMons.iteritems():
            assert(label not in out.keys())
            out[label] = nest.GetStatus(val)

        return out


    def saveAll(self, fileName):
        '''
        Save all the simulated data that has been recorded into a file.

        Parameters
        ----------
        fileName : string
            Path and name of the file
        '''
        out = DataStorage.open(fileName, 'w')
        d = self.getAllData()
        for key, val in d.iteritems():
            out[key] = val
        out.close()



