'''Nest-specific implementation of the grid cell model.'''
from __future__ import absolute_import, print_function

import logging
import collections

import numpy as np
from scipy.io import loadmat
import nest

from . import gc_neurons
from .gc_net import GridCellNetwork
from .place_input import PlaceCellInput
from .place_cells import UniformBoxPlaceCells
from ..data_storage import DataStorage

logger = logging.getLogger(__name__)
gcnLogger = logging.getLogger('{0}.{1}'.format(__name__,
                                               "NestGridCellNetwork"))

nest.Install('gridcellsmodule')


class PosInputs(object):
    '''Data representing animal position input.'''
    def __init__(self, pos_x, pos_y, pos_dt):
        self.pos_x = pos_x
        self.pos_y = pos_y
        self.pos_dt = pos_dt

    def __str__(self):
        res = ("PosInputs:\n  pos_x: {0}\n  pos_y: {1}\n  "
               "pos_dt: {2}".format(self.pos_x, self.pos_y, self.pos_dt))
        return res


class ConstPosInputs(PosInputs):
    '''Data representing constant position of the animal.'''
    def __init__(self, pos_x, pos_y):
        # dt here is irrelevant (say 1e3). This data will never get advanced
        super(ConstPosInputs, self).__init__([float(pos_x)], [float(pos_y)],
                                             1e3)

    def __str__(self):
        res = ('ConstPosInputs:\n  pos_x: {0}\n  pos_y: {1}\n  '
               'pos_dt: {2}'.format(self.pos_x, self.pos_y, self.pos_dt))
        return res


class NestGridCellNetwork(GridCellNetwork):
    '''Grid cell network implemented in NEST simulator.'''
    def __init__(self, neuronOpts, simulationOpts):
        GridCellNetwork.__init__(self, neuronOpts, simulationOpts)
        self.velocityInputInitialized = False

        self.spikeMon_e = None
        self.spikeMon_i = None
        self.stateMon_e = None
        self.stateMon_i = None

        # Extra monitors
        self._extraSpikeMons = {}   # Extra spike monitors
        self._extraStateMons = {}   # Extra state monitors

        self._ratVelocitiesLoaded = False
        self._placeCellsLoaded = False
        self._i_placeCellsLoaded = False

        self.PC = []
        self.PC_start = []

        self.IPC = []
        self.IPCHelper = None
        self.NIPC = None

        self._initNESTKernel()
        self._constructNetwork()
        self._initStates()
        self._initCellularProperties()

    def uniformDistrib(self, mean, spread, N):
        '''Generate a uniform distribution of neurons parameters.

        Parameters
        ----------
        mean : float
            Mean of the distribution
        spread : float
            Width of the distribution around mean.
        N : float
            Number of numbers to generate.

        Returns
        -------
        An array of numbers drawn from this distribution.
        '''
        return mean - spread / 2.0 * np.random.rand(N)

    def _initStates(self):
        '''Initialise states of E and I neurons randomly.'''
        nest.SetStatus(self.E_pop, 'V_m', (self.no.EL_e +
                                           (self.no.Vt_e - self.no.EL_e) *
                                           np.random.rand(len(self.E_pop))))
        nest.SetStatus(self.I_pop, 'V_m', (self.no.EL_i +
                                           (self.no.Vt_i - self.no.EL_i) *
                                           np.random.rand(len(self.I_pop))))

    def _initCellularProperties(self):
        '''Initialise the cellular properties of neurons in the network.'''
        EL_e    = self.uniformDistrib(self.no.EL_e, self.no.EL_e_spread,
                                      len(self.E_pop))
        taum_e  = self.uniformDistrib(self.no.taum_e, self.no.taum_e_spread,
                                      len(self.E_pop))
        EL_i    = self.uniformDistrib(self.no.EL_i, self.no.EL_i_spread,
                                      len(self.I_pop))
        taum_i  = self.uniformDistrib(self.no.taum_i, self.no.taum_i_spread,
                                      len(self.I_pop))
        nest.SetStatus(self.E_pop, 'E_L', EL_e)
        nest.SetStatus(self.E_pop, 'C_m', taum_e * self.no.gL_e)
        nest.SetStatus(self.I_pop, 'E_L', EL_i)
        nest.SetStatus(self.I_pop, 'C_m', taum_i * self.no.gL_i)

    # def _initClocks(self):
    #     for clk in self._clocks:
    #         clk.reinit()

    def _initNESTKernel(self):
        '''Initialise the NEST kernel.'''
        gcnLogger.debug('Initializing NEST kernel: no. of threads: %d',
                        self.no.nthreads)
        nest.ResetKernel()
        nest.SetKernelStatus({"resolution": self.no.sim_dt,
                              "print_time": False})
        nest.SetKernelStatus({"local_num_threads": self.no.nthreads})

    # def reinit(self):
    #     self._initNESTKernel()
    #     self._initStates()
    #     self._initClocks()

    def _constructNetwork(self):
        '''Construct the E/I network'''
        self.e_neuron_params = gc_neurons.getENeuronParams(self.no)
        self.i_neuron_params = gc_neurons.getINeuronParams(self.no)

        self.B_GABA = 1.0   # Must be here for compatibility with brian code

        self.e_model_name = "iaf_gridcells"
        self.i_model_name = "iaf_gridcells"
        self.e_receptors = \
            nest.GetDefaults(self.e_model_name)['receptor_types']
        self.i_receptors = \
            nest.GetDefaults(self.i_model_name)['receptor_types']
        self.E_pop = nest.Create(self.e_model_name, self.net_Ne,
                                 params=self.e_neuron_params)
        self.I_pop = nest.Create(self.i_model_name, self.net_Ni,
                                 params=self.i_neuron_params)

        nest.CopyModel(
            'static_synapse', 'I_AMPA_NMDA',
            params={'receptor_type': self.i_receptors['AMPA_NMDA']})
        nest.CopyModel(
            'static_synapse', 'E_GABA_A',
            params={'receptor_type': self.e_receptors['GABA_A']})
        nest.CopyModel(
            'static_synapse', 'PC_AMPA',
            params={'receptor_type': self.e_receptors['AMPA']})

        # Connect E-->I and I-->E
        self._connect_network()

    def simulate(self, time, printTime=True):
        '''Run the simulation'''
        self.endConstruction()
        self.beginSimulation()
        nest.SetKernelStatus({"print_time": bool(printTime)})

        if not self.velocityInputInitialized:
            velMsg = ("Velocity input has not been initialized. Make sure "
                      "this is the desired behavior. If you have set the "
                      "'velON' parameter to 1, then this message probably "
                      "indicates a bug in the simulation code.")
            gcnLogger.warn(velMsg)
        nest.Simulate(time)

    def getSpikeDetector(self, type, N_ids=None):
        '''
        Get a spike detector that records from neurons given N_ids and from the
        population type given by type
        '''
        if type == "E":
            if self.spikeMon_e is not None:
                return self.spikeMon_e
            else:
                if N_ids is None:
                    N_ids = np.arange(len(self.E_pop))
                src = list(np.array(self.E_pop)[N_ids])
                self.spikeMon_e = nest.Create('spike_detector')
                nest.SetStatus(self.spikeMon_e, {
                    "label"   : "E spikes",
                    'withtime': True,
                    'withgid' : True})
                nest.ConvergentConnect(src, self.spikeMon_e)
                return self.spikeMon_e
        elif type == "I":
            if self.spikeMon_i is not None:
                return self.spikeMon_i
            else:
                if N_ids is None:
                    N_ids = np.arange(len(self.I_pop))
                src = list(np.array(self.I_pop)[N_ids])
                self.spikeMon_i = nest.Create('spike_detector')
                nest.SetStatus(self.spikeMon_i, {
                    "label"   : "I spikes",
                    'withtime': True,
                    'withgid' : True})
                # print src
                nest.ConvergentConnect(src, self.spikeMon_i)
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
            "label"   : label,
            'withtime': True,
            'withgid' : True})
        nest.ConvergentConnect(gids, mon)

        self._extraSpikeMons[label] = (mon, gids[0])
        return mon

    def getStateMonitor(self, type, N_ids, params):
        '''
        Return a state monitor for a given population (type) and relative
        indexes of neurons (N_ids), with parameters given by params
        '''
        if len(N_ids) == 0:
            raise ValueError("State monitor needs to record from at least one "
                             "neuron")

        N = len(N_ids)

        if type == "E":
            if self.stateMon_e is None:
                self.stateMon_e = nest.Create('multimeter', N, params=params)
                nest.Connect(self.stateMon_e, self.E_pop[0] + np.array(N_ids))
            return self.stateMon_e
        elif type == "I":
            if self.stateMon_i is None:
                self.stateMon_i = nest.Create('multimeter', N, params=params)
                nest.Connect(self.stateMon_i, self.I_pop[0] + np.array(N_ids))
            return self.stateMon_i

    def getGenericStateMonitor(self, gids, params, label):
        '''Create a user defined state monitor (multimeter)

        @param gids   Global IDs of the neurons.
        @param params Parameter of the state monitors
        @return Global IDs for manipulation in the nest space
        '''
        if len(gids) == 0:
            logger.warn('Requested to create 0 state monitors. Ignoring...')
            return []
        mon = nest.Create('multimeter', len(gids), params=params)
        nest.Connect(mon, gids)
        self._extraStateMons[label] = (mon)
        return mon

    def _divergentConnectEE(self, pre, post, weights):
        post_global = list(self.E_pop[0] + np.asanyarray(post))
        nest.DivergentConnect([self.E_pop[0] + pre], post_global,
                              model='I_AMPA_NMDA', weight=list(weights),
                              delay=[self.no.delay] * len(weights))

    def _divergentConnectEI(self, pre, post, weights):
        post_global = list(self.I_pop[0] + np.array(post))
        nest.DivergentConnect([self.E_pop[0] + pre], post_global,
                              model='I_AMPA_NMDA', weight=list(weights),
                              delay=[self.no.delay] * len(weights))

    def _randomDivergentConnectEI(self, pre, post, n, weights):
        '''Connect each neuron in ``pre`` (E population) to n randomly selected
        neurons in ``post`` (I population), with weights specified in
        ``weights``. If weights is a float then all the weights are constant.
        '''
        if isinstance(weights, collections.Iterable):
            delay = [self.no.delay] * len(weights)
        else:
            delay = self.no.delay
        nest.RandomDivergentConnect(
            (self.E_pop[0] + np.asanyarray(pre)).tolist(),
            (self.I_pop[0] + np.asanyarray(post)).tolist(),
            n,
            weight=weights,
            model='I_AMPA_NMDA',
            delay=delay)

    def _randomDivergentConnectIE(self, pre, post, n, weights):
        '''Connect each neuron in ``pre`` (I population) to n randomly selected
        neurons in ``post`` (E population), with weights specified in
        ``weights``. If weights is a float then all the weights are constant.
        '''
        if isinstance(weights, collections.Iterable):
            delay = [self.no.delay] * len(weights)
        else:
            delay = self.no.delay
        nest.RandomDivergentConnect(
            (self.I_pop[0] + np.asanyarray(pre)).tolist(),
            (self.E_pop[0] + np.asanyarray(post)).tolist(),
            n,
            weight=weights,
            model='E_GABA_A',
            delay=delay)

    def _randomDivergentConnectII(self, pre, post, n, weights,
                                  allow_autapses=False, allow_multapses=False):
        '''Connect each neuron in ``pre`` (I population) to n randomly selected
        neurons in ``post`` (I population), with weights specified in
        ``weights``. If weights is a float then all the weights are constant.
        '''
        if isinstance(weights, collections.Iterable):
            delay = [self.no.delay] * len(weights)
        else:
            delay = self.no.delay
        nest.RandomDivergentConnect(
            (self.I_pop[0] + np.asanyarray(pre)).tolist(),
            (self.I_pop[0] + np.asanyarray(post)).tolist(),
            n,
            weight=weights,
            model='E_GABA_A',
            delay=delay,
            options={'allow_autapses': allow_autapses,
                     'allow_multapses': allow_multapses})

    def _divergentConnectIE(self, pre, post, weights):
        post_global = list(self.E_pop[0] + np.array(post))
        nest.DivergentConnect([self.I_pop[0] + pre], post_global,
                              model='E_GABA_A', weight=list(weights),
                              delay=[self.no.delay] * len(weights))

    def getConnMatrix(self, popType):
        '''
        Return all *input* connections to neuron with index post from the
        specified popType.

        Parameters
        ----------
        popType : string, 'E' or 'I'
            Type of the population. If popType == 'E', return connection
            weights for AMPA connections only. The NMDA connections will be a
            fraction of the AMPA connection strength specified by the
            NMDA_amount parameter.

            If popType == 'I' the connection weights returned will be for
            GABA_A connections.
        output : a 2D numpy array
            An array containing the connections. The shape is (post,
            pre)/(target, source).
        '''
        EStart = np.min(self.E_pop)
        IStart = np.min(self.I_pop)
        print("EStart: {0}".format(EStart))
        print("IStart: {0}".format(IStart))
        print("len(self.E_pop): {0}".format(len(self.E_pop)))
        print("len(self.I_pop): {0}".format(len(self.I_pop)))
        if popType == 'E':
            W_IE = np.zeros((len(self.I_pop), len(self.E_pop)))
            for e in xrange(len(self.E_pop)):
                print("E neuron {0} --> I neurons".format(e))
                conns = nest.FindConnections([self.E_pop[e]])
                for i in xrange(len(conns)):
                    target = nest.GetStatus([conns[i]], 'target')
                    if target[0] in self.I_pop:
                        W_IE[target[0] - IStart, e] = \
                            nest.GetStatus([conns[i]], 'weight')[0]
            return W_IE
        elif popType == 'I':
            W_EI = np.zeros((len(self.E_pop), len(self.I_pop)))
            for i in xrange(len(self.I_pop)):
                print("I neuron {0} --> E neurons".format(i))
                conns = nest.FindConnections([self.I_pop[i]])
                for e in xrange(len(conns)):
                    target = nest.GetStatus([conns[e]], 'target')
                    if target[0] in self.E_pop:
                        W_EI[target[0] - EStart, i] = \
                            nest.GetStatus([conns[e]], 'weight')[0]
            return W_EI

        else:
            msg = 'popType must be either \'E\' or \'I\'. Got {0}'
            raise ValueError(msg.format(popType))

    ###########################################################################
    #                     External sources definitions
    ###########################################################################
    def _loadRatVelocities(self):
        '''
        Load rat velocities (in this case positions only)
        '''
        if self._ratVelocitiesLoaded:
            return

        logger.info('Loading rat velocities')

        self.ratData    = loadmat(self.no.ratVelFName)
        self.rat_dt     = self.ratData['dt'][0][0] * 1e3      # units: ms

        self.rat_pos_x  = self.ratData['pos_x'].ravel()
        self.rat_pos_y  = self.ratData['pos_y'].ravel()

        # Map velocities to currents: we use the slope of bump speed vs. rat
        # speed and inter-peak grid field distance to remap
        # Bump speed-current slope must be estimated
        self.velC = self.Ne_x / self.no.gridSep / self.no.bumpCurrentSlope

        self._ratVelocitiesLoaded = True

        gcnLogger.debug('velC: %f, bumpCurrentSlope: %f, gridSep: %f',
                        self.velC, self.no.bumpCurrentSlope, self.no.gridSep)

    def setVelocityCurrentInput_e(self, prefDirs_mask=None):
        '''
        Set up movement simulation, based on preferred directions of neurons.
        prefDirs_mask can be used to manipulate velocity input strength
        for each neuron.
        '''
        logger.info("Setting up velocity input current.")
        self._loadRatVelocities()

        if prefDirs_mask is None:
            self._prefDirs_mask_e = np.ndarray((len(self.E_pop), 2))
            self._prefDirs_mask_e[:, :] = 1.0
        else:
            raise NotImplementedError()

        # Load velocities into nest: they are all shared among all
        # iaf_gridcells nodes so only one neuron needs setting the actual
        # values
        npos = int(self.no.time / self.rat_dt)
        nest.SetStatus([self.E_pop[0]], {
            "rat_pos_x" : self.rat_pos_x[0:npos].tolist(),
            "rat_pos_y" : self.rat_pos_y[0:npos].tolist(),
            "rat_pos_dt": self.rat_dt})  # s --> ms

        nest.SetStatus(self.E_pop, "pref_dir_x", self.prefDirs_e[:, 0])
        nest.SetStatus(self.E_pop, "pref_dir_y", self.prefDirs_e[:, 1])
        nest.SetStatus(self.E_pop, "velC", self.velC)

        self.velocityInputInitialized = True

    def setConstantVelocityCurrent_e(self, vel, start_t=None, end_t=None):
        '''
        Set the model so that there is only a constant velocity current input.
        '''
        gcnLogger.info('Setting up constant velocity current '
                       'input: {0}'.format(vel))
        if start_t is not None:
            raise Exception("Const velocity start time cannot be overridden "
                            "in this model!")

        start_t = self.no.theta_start_t

        if end_t is None:
            end_t = self.no.time

        self.rat_dt = 20.0  # ms
        nVel = int((end_t - start_t) / self.rat_dt)
        self.rat_pos_x = np.cumsum(np.array([vel[0]] * nVel)) * (self.rat_dt *
                                                                 1e-3)
        self.rat_pos_y = np.cumsum(np.array([vel[1]] * nVel)) * (self.rat_dt *
                                                                 1e-3)

        # Force these velocities, not the animal velocitites
        self._ratVelocitiesLoaded = True

        # Load velocities into nest: they are all shared among all
        # iaf_gridcells nodes so only one neuron needs setting the actual
        # values
        nest.SetStatus([self.E_pop[0]], {
            "rat_pos_x" : self.rat_pos_x.tolist(),
            "rat_pos_y" : self.rat_pos_y.tolist(),
            "rat_pos_dt": self.rat_dt})  # s --> ms

        print(self.rat_pos_x)
        print(self.rat_pos_y)

        # Map velocities to currents: Here the mapping is 1:1, i.e. the
        # velocity dictates the current
        self.velC = 1.

        nest.SetStatus(self.E_pop, "pref_dir_x", self.prefDirs_e[:, 0])
        nest.SetStatus(self.E_pop, "pref_dir_y", self.prefDirs_e[:, 1])
        nest.SetStatus(self.E_pop, "velC", self.velC)

        self.setStartPlaceCells(PosInputs([0.], [.0], self.rat_dt))

        self.velocityInputInitialized = True

    def setStartPlaceCells(self, posIn):
        '''Create and connect the initialisation place cells.'''
        if len(self.PC_start) == 0:
            gcnLogger.info("Setting up initialization place cells")
            gcnLogger.debug("Init place cell positional input: {0}".format(
                str(posIn)))
            gcnLogger.debug("Init place cells: start: {0}, end: {1}".format(
                0, self.no.theta_start_t))
            self.PC_start, _, _ = self.createGenericPlaceCells(
                self.no.N_place_cells,
                self.no.pc_start_max_rate,
                self.no.pc_start_conn_weight,
                start=0.0,
                end=self.no.theta_start_t,
                posIn=posIn)
        else:
            gcnLogger.info('Initialization place cells already set. Skipping '
                           'the set up')

    def setPlaceCells(self, start=None, end=None, posIn=None):
        '''Place cells to initialize the bump.

        It should be initialized onto the correct position, i.e. the bump must
        be at the correct starting position, which matches the actual velocity
        simulation place cell input.
        '''
        if posIn is None:
            self._loadRatVelocities()
            startPos = ConstPosInputs(self.rat_pos_x[0], self.rat_pos_y[0])
        else:
            startPos = ConstPosInputs(posIn.pos_x[0], posIn.pos_y[0])
        self.setStartPlaceCells(startPos)

        # Here the actual velocity place cells
        gcnLogger.info("Setting up place cells. User defined positional "
                       "data: {0}".format('no' if posIn is None else 'yes'))
        gcnLogger.debug("Place cell positional input: {0}".format(str(posIn)))

        self.PC, _, _ = self.createGenericPlaceCells(self.no.N_place_cells,
                                                     self.no.pc_max_rate,
                                                     self.no.pc_conn_weight,
                                                     start, end, posIn)

    def createGenericPlaceCells(self, N, maxRate, weight, start=None, end=None,
                                posIn=None):
        '''
        Generate place cells and connect them to grid cells. The wiring is
        fixed, and there is no plasticity. This method can be used more than
        once, to set up different populations of place cells.
        '''
        if start is None:
            start = self.no.theta_start_t
        if end is None:
            end = self.no.time
        if posIn is None:
            self._loadRatVelocities()
            posIn = PosInputs(self.rat_pos_x, self.rat_pos_y, self.rat_dt)

        if N != 0:
            gcnLogger.info('Setting up generic place cells')
            NTotal = N * N

            boxSize = [self.no.arenaSize, self.no.arenaSize]
            PCHelper = UniformBoxPlaceCells(boxSize, (N, N), maxRate,
                                            self.no.pc_field_std, random=False)

            PC = nest.Create('place_cell_generator', NTotal,
                             params={'rate'      : maxRate,
                                     'field_size': self.no.pc_field_std,
                                     'start'     : start,
                                     'stop'      : end})
            nest.SetStatus(PC, 'ctr_x', PCHelper.centers[:, 0])
            nest.SetStatus(PC, 'ctr_y', PCHelper.centers[:, 1])

            npos = int(self.no.time / posIn.pos_dt)
            nest.SetStatus([PC[0]], params={
                'rat_pos_x' : list(posIn.pos_x[0:npos]),
                'rat_pos_y' : list(posIn.pos_y[0:npos]),
                'rat_pos_dt': posIn.pos_dt})

            # test_x = nest.GetStatus([PC[0]], 'rat_pos_x')
            # test_y = nest.GetStatus([PC[0]], 'rat_pos_y')
            # print test_x, test_y

            # Connections
            # Here we extract connections from the PlaceCellInput class that
            # was originaly used as a current input generator for place cell
            # resetting mechanism. The output of this class perfectly matches
            # how divergent connections from a single place cell should be
            # mapped onto the twisted torus grid cell sheet

            # how divergent the connections are, 3sigma rule --> division by 6.
            connStdDev          = self.no.gridSep / 2. / 6.
            pc_weight_threshold = 0.1

            pc_input = PlaceCellInput(self.Ne_x, self.Ne_y, self.no.arenaSize,
                                      self.no.gridSep, [.0, .0],
                                      fieldSigma=connStdDev)
            ctr_x = nest.GetStatus(PC, 'ctr_x')
            ctr_y = nest.GetStatus(PC, 'ctr_y')
            for pc_id in xrange(NTotal):
                w = pc_input.getSheetInput(ctr_x[pc_id],
                                           ctr_y[pc_id]).flatten()
                gt_th = w > pc_weight_threshold
                post = np.array(self.E_pop)[gt_th]
                w    = w[gt_th]
                # print post, w
                nest.DivergentConnect(
                    [PC[pc_id]],
                    list(post),
                    weight=list(w * weight),
                    delay=[self.no.delay] * len(w),
                    model='PC_AMPA')

            return PC, PCHelper, NTotal

        else:
            gcnLogger.warn("trying to set up place cells with N_place_cells "
                           "== 0")

        self._placeCellsLoaded = True

    def setIPlaceCells(self):
        self._createIPlaceCells(self.no.ipc_N,
                                int(self.no.ipc_nconn),
                                self.no.ipc_max_rate,
                                self.no.ipc_weight,
                                self.no.ipc_field_std)

    def _createIPlaceCells(self, N, Nconn_pcs, maxRate, weight, field_std,
                           start=None, end=None, posIn=None):
        '''
        Generate place cells and connect them to I cells. The wiring is
        fixed, and there is no plasticity. This method can be used more than
        once, to set up different populations of place cells.

        Here the widths of the place fields are the same as in the case of the
        generic place cells.

        Parameters
        ----------
        Nconn_pcs : int
            Number of place cells connected to each I neurons.
        '''
        if start is None:
            start = self.no.theta_start_t
        if end is None:
            end = self.no.time
        if posIn is None:
            self._loadRatVelocities()
            posIn = PosInputs(self.rat_pos_x, self.rat_pos_y, self.rat_dt)

        NTotal = N * N
        PC = None
        PCHelper = None

        if N != 0:
            gcnLogger.info('Setting up place cells connected to I cells')
            gcnLogger.info("N: %d, Nconn_pcs: %d, maxRate: %f, weight: %f, field_std: %f", N,
                           int(Nconn_pcs), maxRate, weight, field_std)

            boxSize = [self.no.arenaSize, self.no.arenaSize]
            PCHelper = UniformBoxPlaceCells(boxSize, (N, N), maxRate,
                                            field_std, random=False)

            PC = nest.Create('place_cell_generator', NTotal,
                             params={'rate'      : maxRate,
                                     'field_size': field_std,
                                     'start'     : start,
                                     'stop'      : end})
            nest.SetStatus(PC, 'ctr_x', PCHelper.centers[:, 0])
            nest.SetStatus(PC, 'ctr_y', PCHelper.centers[:, 1])

            npos = int(self.no.time / posIn.pos_dt)
            nest.SetStatus([PC[0]], params={
                'rat_pos_x' : list(posIn.pos_x[0:npos]),
                'rat_pos_y' : list(posIn.pos_y[0:npos]),
                'rat_pos_dt': posIn.pos_dt})

            # Connections
            # I-PCs are connected with a constant connection weight to I cells
            for i_idx in self.I_pop:
                pc_idx = np.random.choice(PC, Nconn_pcs, replace=False)
                nest.ConvergentConnect(pc_idx.tolist(), [i_idx],
                                       weight=weight,
                                       delay=self.no.delay,
                                       model='PC_AMPA')
        else:
            gcnLogger.warn("Trying to set up I place cells with 0 place cells.")

        self._i_placeCellsLoaded = True
        self.IPC = PC
        self.IPCHelper = PCHelper
        self.NIPC = NTotal
        #self._getIPCConnections()

    def _getIPCConnections(self):
        IStart = self.I_pop[0]
        W = np.zeros((len(self.I_pop), len(self.IPC)))
        for pcn in xrange(len(self.IPC)):
            print("IPC {0} --> I neurons".format(pcn))
            conns = nest.FindConnections([self.IPC[pcn]])
            for i in xrange(len(conns)):
                target = nest.GetStatus([conns[i]], 'target')
                if target[0] in self.I_pop:
                    W[target[0] - IStart, pcn] = nest.GetStatus([conns[i]],
                                                                'weight')[0]
                else:
                    print("Target not in I_pop!")
        return W

    ###########################################################################
    #                                   Other
    ###########################################################################
    def getRatData(self):
        '''Return the data representing the animal (rat).'''
        return self.ratData

    def getAttrDictionary(self):
        d = {}

        d['e_neuron_params'] = self.e_neuron_params
        d['i_neuron_params'] = self.i_neuron_params
        d['B_GABA'         ] = self.B_GABA
        d['E_pop'          ] = np.array(self.E_pop)
        d['I_pop'          ] = np.array(self.I_pop)
        d['PC'             ] = np.array(self.PC)
        d['PC_start'       ] = np.array(self.PC_start)
        d['net_Ne'         ] = self.net_Ne
        d['net_Ni'         ] = self.net_Ni
        d['rat_pos_x'      ] = getattr(self, 'rat_pos_x', np.nan)
        d['rat_pos_y'      ] = getattr(self, 'rat_pos_y', np.nan)
        d['rat_dt'         ] = getattr(self, 'rat_dt', np.nan)
        d['velC'           ] = getattr(self, 'velC', np.nan)
        d['Ne_x'           ] = self.Ne_x
        d['Ne_y'           ] = self.Ne_y
        d['Ni_x'           ] = self.Ni_x
        d['Ni_y'           ] = self.Ni_y
        d['prefDirs_e'     ] = self.prefDirs_e

        return d


class BasicGridCellNetwork(NestGridCellNetwork):
    '''
    A grid cell network that generates the common network and creates a basic
    set of spike monitors and state monitors which are generically usable in
    most of the simulation setups.
    '''
    def getDefaultStateMonParams(self):
        '''Generate default, pre-set state monitor parameters.'''
        return {
            'withtime': True,
            'interval': 10.0 * self.no.sim_dt,
            'record_from': ['V_m', 'I_clamp_AMPA', 'I_clamp_NMDA',
                            'I_clamp_GABA_A', 'I_stim']
        }

    def fillParams(self, dest, src):
        '''Fill properties into a dictionary.'''
        for key, value in src.iteritems():
            dest[key] = value
        return dest

    def __init__(self, options, simulationOpts=None,
                 nrec_spikes=(None, None),
                 stateRecord_type='middle-center',
                 stateRecParams=(None, None),
                 rec_spikes_probabilistic=False):
        '''
        TODO
        '''
        NestGridCellNetwork.__init__(self, options, simulationOpts)

        # Spikes
        self.nrecSpikes_e = nrec_spikes[0]
        self.nrecSpikes_i = nrec_spikes[1]

        if self.nrecSpikes_e is None:
            self.nrecSpikes_e = self.Ne_x * self.Ne_y
        if self.nrecSpikes_i is None:
            self.nrecSpikes_i = self.Ni_x * self.Ni_y

        if rec_spikes_probabilistic == False:
            self.spikeMon_e  = self.getSpikeDetector("E",
                                                     np.arange(self.nrecSpikes_e))
            self.spikeMon_i  = self.getSpikeDetector("I",
                                                     np.arange(self.nrecSpikes_i))
        else:
            self.spikeMon_e  = self.getSpikeDetector(
                "E", np.sort(np.random.choice(len(self.E_pop), self.nrecSpikes_e,
                                      replace=False)))
            self.spikeMon_i  = self.getSpikeDetector(
                "I", np.sort(np.random.choice(len(self.I_pop),
                                              self.nrecSpikes_i,
                                              replace=False)))

        # States
        if stateRecord_type == 'middle-center':
            self.state_record_e = [self.Ne_x / 2 - 1,
                                   (self.Ne_y / 2 * self.Ne_x +
                                    self.Ne_x / 2 - 1)]
            self.state_record_i = [self.Ni_x / 2 - 1,
                                   (self.Ni_y / 2 * self.Ni_x +
                                    self.Ni_x / 2 - 1)]
        else:
            raise ValueError("Currently stateRecordType must be "
                             "'middle-center'")

        self.stateMonParams_e = self.getDefaultStateMonParams()
        self.stateMonParams_i = self.getDefaultStateMonParams()

        stRecp_e = stateRecParams[0]
        stRecp_i = stateRecParams[1]
        if stRecp_e is not None:
            self.fillParams(self.stateMonParams_e, stRecp_e)
        if stRecp_i is not None:
            self.fillParams(self.stateMonParams_i, stRecp_i)

        self.stateMon_e  = self.getStateMonitor("E",
                                                self.state_record_e,
                                                self.stateMonParams_e)
        self.stateMon_i  = self.getStateMonitor("I",
                                                self.state_record_i,
                                                self.stateMonParams_i)

    def getMonitors(self):
        '''Return the main spike and state monitors.'''
        return (
            self.spikeMon_e,
            self.spikeMon_i,
            self.stateMon_e,
            self.stateMon_i
        )

    def getSpikeMonData(self, mon, gidStart):
        '''
        Generate a dictionary of a spike data from the monitor ``mon``

        Notes
        -----
        NEST has some troubles with consistency in returning data in a correct
        format on OSX and Linux. On Linux, sequential data are apparently
        returned as np.ndarray, while on OSX the data are returned as lists.
        This obviously causes huge data files on OSX since the data are stored
        as lists into HDF5.
        '''
        st = nest.GetStatus(mon)[0]
        events = st['events']
        for key in events.keys():
            events[key] = np.asanyarray(events[key])
        events['senders'] -= gidStart
        return st

    def getStateMonData(self, mon):
        '''
        Generate a dictionary of state monitor data from the monitor ``mon``

        Notes
        -----
        NEST has some troubles with consistency in returning data in a correct
        format on OSX and Linux. On Linux, sequential data are apparently
        returned as np.ndarray, while on OSX the data are returned as lists.
        This obviously causes huge data files on OSX since the data are
        serialized as lists into HDF5.
        '''
        out = nest.GetStatus(mon)
        for mon_idx in range(len(out)):
            events = out[mon_idx]['events']
            for key in events.keys():
                events[key] = np.asanyarray(events[key])
        return out

    def getSpikes(self, **kw):
        '''
        Return a dictionary of spike monitor data.
        '''
        out = {}

        if self.spikeMon_e is not None:
            out['spikeMon_e'] = self.getSpikeMonData(self.spikeMon_e,
                                                     self.E_pop[0])
        if self.spikeMon_i is not None:
            out['spikeMon_i'] = self.getSpikeMonData(self.spikeMon_i,
                                                     self.I_pop[0])

        for label, vals in self._extraSpikeMons.iteritems():
            assert label not in out.keys()
            out[label] = self.getSpikeMonData(vals[0], vals[1])

        return out

    def getNetParams(self):
        '''Get network and derived network parameters.'''
        out = {}
        out['options'] = self.no._einet_optdict
        out['net_attr'] = self.getAttrDictionary()
        return out

    def getAllData(self):
        '''
        Save all the simulated data into a dictionary and return it.
        '''
        out = self.getNetParams()

        # Spike monitors
        # Note that getSpikes() is overridden in child classes and requires the
        # espikes and ispikes arguments.
        out.update(self.getSpikes(espikes=True, ispikes=True))

        # Save state variables
        out['stateMon_e'] = self.getStateMonData(self.stateMon_e)
        out['stateMon_i'] = self.getStateMonData(self.stateMon_i)
        for label, val in self._extraStateMons.iteritems():
            assert label not in out.keys()
            out[label] = self.getStateMonData(val)

        return out

    def saveSpikes(self, fileName):
        '''
        Save all the simulated spikes that have been recorded into a file.

        Parameters
        ----------
        fileName : string
            Path and name of the file
        '''
        out = DataStorage.open(fileName, 'w')
        d = self.getSpikes()  # FIXME: what is d?
        out.close()

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


class ConstantVelocityNetwork(BasicGridCellNetwork):
    '''
    A grid cell network that simulates a constant velocity in a specified
    direction.
    '''
    def __init__(self, options, simulationOpts=None,
                 vel=[0.0, 0.0],
                 nrec_spikes=(None, None),
                 stateRecord_type='middle-center',
                 stateRecParams=(None, None)):
        '''
        Generate the network.

        Parameters
        ----------
        vel : a pair [x, y]
            Velocity input vector, i.e. it specifies the direction and
            magnitude of the velocity current.
        '''
        BasicGridCellNetwork.__init__(self,
                                      options, simulationOpts,
                                      nrec_spikes,
                                      stateRecord_type,
                                      stateRecParams)

        self.setConstantVelocityCurrent_e(vel)

    def getSpikes(self, **kw):
        '''
        Return a dictionary of spike monitor data.
        For keyword arguments description, see
        :meth:`~ConstantVelocityNetwork.getMinimalSaveData`
        '''
        espikes = kw.get('espikes', True)
        ispikes = kw.get('ispikes', False)

        out = {}
        if espikes:
            out['spikeMon_e'] = self.getSpikeMonData(self.spikeMon_e,
                                                     self.E_pop[0])
        if ispikes:
            out['spikeMon_i'] = self.getSpikeMonData(self.spikeMon_i,
                                                     self.I_pop[0])
        return out

    def getMinimalSaveData(self, **kw):
        '''
        Keyword arguments:
        ``espikes`` : bool
            Whether to return spikes from the E population. Defaults to True.
        ``ispikes`` : bool
            Whether to return spikes from the I population. Defaults to False.
        '''
        out = self.getNetParams()
        out.update(self.getSpikes(**kw))
        return out
