'''Simulator independent grid cell network code.

  Grid cell network set up
  ------------------------

  This file is a module for the grid cell network. It allows you to create a
  network of Exponential integrate and fire neurons. There are two populations:
    Stellate cells -- excitatory neurons (E)
    Interneurons   -- inhibitory (I)

  In the current model, the connections are only E --> I and I --> E. That is
  the main idea behind the model, which is supported by experimental evidence.

  Both neuron types can receive AMPA, NMDA and GABA_A events. NMDA is only in
  E --> I direction.

  The topology of connections has several characteristics:

   - Both neuron models are uniformly placed on a surface of a twisted torus,
     the sides of which are scaled as X:Y = 1:sqrt(3)/2, to accomodate the
     hexagonal nature of grid-cell receptive fields.

   - The user can set the number of neurons in the larger side of the toroidal
     sheet (always X dimension)

   - The connections follow a center-surround profile, i.e. either E-->I
     connections have a surround profile and I-->E connections have a center
     profile, or vice-versa. This can be used to test the effect of the type of
     excitatory or inhibitory profile on the stability of the attractor

   - GABA_A connections (I-->E) can also contain extra, randomly generated
     inhibitory synapses onto stellate cells in order to allow generation of
     gamma oscillations.

   - Technically, the basic functionality of the model (attractor emergence,
     oscillations), shouldn't be very dependent on the spiking neuron type.
     After some parameter setups/changes, one should be able to set the
     simulation with any kind of spiking model (leaky IaF, Hodgkin-Huxley,
     etc.)
'''
from __future__ import absolute_import, print_function

import logging
import numpy as np
import time
import copy

from ..analysis.image import Position2D, remapTwistedTorus
from .construction.weights import (IsomorphicConstructor,
                                   ProbabilisticConstructor)


__all__ = ['GridCellNetwork']


gcnLogger = logging.getLogger('{0}.{1}'.format(__name__,
                                               "NestGridCellNetwork"))


class GridCellNetwork(object):
    '''
    This is an interface to the grid cell network. One should be able to set
    it up, simulate and record from neurons, quite independently of a specific
    simulator. We don't use PyNN here because its generic interface is not
    suitable for this purpose

    The GridCellNetwork creates two separate populations and connects them
    according to the specified connectivity rules.
    '''
    def __init__(self, neuronOpts, simulationOpts):
        # timers
        self._startT = time.time()
        self._constrEndT   = None
        self._simStartT    = None
        self._simEndT      = None
        self.beginConstruction()

        self.no = copy.deepcopy(neuronOpts)
        self.so = copy.deepcopy(simulationOpts)

        # Setup neuron numbers for each dimension (X, Y)
        # We have a single bump and to get hexagonal receptive fields the X:Y
        # size ratio must be 1:sqrt(3)/2
        self.y_dim = np.sqrt(3) / 2.0
        self.Ne_x = self.no.Ne
        self.Ne_y = int(np.ceil(self.no.Ne * self.y_dim)) // 2 * 2

        self.Ni_x = self.no.Ni
        self.Ni_y = int(np.ceil(self.no.Ni * self.y_dim)) // 2 * 2

        self.net_Ne = self.Ne_x * self.Ne_y
        self.net_Ni = self.Ni_x * self.Ni_y

        self.prefDirs_e = None
        self.prefDirs_i = None

        self._weight_constructor = self._select_weight_constructor(neuronOpts)

    def simulate(self, t, printTime):
        '''Simulate the network, after being set up.'''
        raise NotImplementedError()

    @staticmethod
    def _select_weight_constructor(options):
        '''Create an instance of the constructor, based on ``options``.

        Parameters
        ----------
        options : dict-like
            A mapping that contains options necessary to select the appropriate
            constructor.

        Returns
        -------
        constructor : WeightConstructor
            Constructor instance.
        '''
        if options.probabilistic_synapses:
            gcnLogger.debug('Selecting probabilistic contructor for weights.')
            return ProbabilisticConstructor()
        else:
            gcnLogger.debug('Selecting isomorphic contructor for weights.')
            return IsomorphicConstructor()

    def _divergentConnectEE(self, pre, post, weights):
        '''Connect a ``pre`` neuron in the E population to all neurons in the E
        population in the ``post``, with ``weights``.
        '''
        raise NotImplementedError()

    def _divergentConnectEI(self, pre, post, weights):
        '''
        Simply connect a 'pre' neuron in the E population to all neurons in
        the I population from post, with given weights
        '''
        raise NotImplementedError()

    def _divergentConnectIE(self, pre, post, weights):
        '''
        Simply connect a 'pre' neuron in the I population to all neurons in
        the E population from post, with given weights
        '''
        raise NotImplementedError()

    def _shiftOnTwistedTorus(self, val, shift, dim):
        '''Shift a pair of X and Y coordinates on a twisted torus in a specified
        direction.

        Parameters
        ----------
        val : Position2D
            The original coordinates.
        shift : Position2D
            The vector that determines the shift.
        dim : Position2D
            Dimensions of the twisted torus.

        Returns
        -------
        new_coord : Position2D
            New coordinates on the twisted torus.
        '''
        ret = Position2D(val.x, val.y)
        ret.x += shift.x
        ret.y += shift.y

        if ret.y < 0 or ret.y >= dim.y:
            ret.x += dim.x / 2.0

        ret.x %= dim.x
        ret.y %= dim.y

        return ret

    def _generateRinglikeWeights(self, a, others, mu, sigma, prefDir,
                                 prefDirC):
        '''Generate ring-like weights.

        Here we assume that X coordinates are normalised to <0, 1), and Y
        coordinates are normalised to <0, sqrt(3)/2) Y coordinates are twisted,
        i.e. X will have additional position shifts when determining minimum.

        @param a        Neuron center, normalised. A Position2D object.
        @param others   Positions of postsynaptic neurons. A Position2D object.
        @param mu       Radius of the circular function
        @param sigma    Width of the circular function
        @param prefDir  Preferred direction of the cell. A Position2D object.
        @param prefDirC A preferred direction coefficient. A multiplier.
        @return An array (1D) of normalized weights.
        '''
        dim = Position2D()
        dim.x = 1.0
        dim.y = self.y_dim

        # a.x -= prefDirC*prefDir.x
        # a.y -= prefDirC*prefDir.y
        shift = Position2D(-prefDirC * prefDir.x, -prefDirC * prefDir.y)
        a = self._shiftOnTwistedTorus(a, shift, dim)

        d = remapTwistedTorus(a, others, dim)
        return np.exp(-(d - mu)**2 / 2 / sigma**2)

    def _generateGaussianWeights(self, a, others, sigma, prefDir, prefDirC):
        '''Generate Gaussian-like weights, i.e. local connections

        Here we assume that X coordinates are normalised to <0, 1), and Y
        coordinates are normalised to <0, sqrt(3)/2) Y coordinates are twisted,
        i.e. X will have additional position shifts when determining minimum.

        @param a        Neuron center, normalised. A Position2D object.
        @param others   Positions of postsynaptic neurons. A Position2D object.
        @param sigma    Std. dev. of the Gaussian (normalised)
        @param prefDir  Preferred direction of the cell. A Position2D object.
        @param prefDirC A preferred direction coefficient. A multiplier.
        @return An array (1D) of normalized weights.
        '''
        dim = Position2D()
        dim.x = 1.0
        dim.y = self.y_dim

        a.x -= prefDirC * prefDir.x
        a.y -= prefDirC * prefDir.y

        d = remapTwistedTorus(a, others, dim)
        return np.exp(-d**2 / 2. / sigma**2)

    def _addToConnections(self, conductances, perc_synapses, h):
        '''
        Picks perc_synapses% of connections from the array and adds h to them
        '''
        indexes = np.random.choice(
            np.arange(len(conductances)),
            size=int(perc_synapses / 100.0 * len(conductances)),
            replace=False)
        conductances[indexes] += h
        return conductances

    def _connect_network(self):
        '''Make network connections according to parameter settings.'''
        if self.no.EI_flat:
            self._connect_ei_flat()
        else:
            self._connect_ei_distance(self.no.AMPA_gaussian, self.no.pAMPA_mu,
                                      self.no.pAMPA_sigma)

        if self.no.IE_flat:
            self._connect_ie_flat()
        else:
            self._connect_ie_distance(self.no.AMPA_gaussian, self.no.pGABA_mu,
                                      self.no.pGABA_sigma)

        if self.no.use_EE:
            self._connect_ee(self.no.pEE_sigma)

        if self.no.use_II:
            self._connect_ii_flat()

    def _connect_ee(self, pEE_sigma):
        '''Make E-->E connections, according to network options.'''
        gcnLogger.info('Connecting E-->E (distance-dependent).')
        g_EE_mean = self.no.g_EE_total / self.net_Ne
        print("g_EE_mean: %f nS" % g_EE_mean)

        others_e  = Position2D()
        pd_norm_e = Position2D()
        a         = Position2D()

        X, Y = np.meshgrid(np.arange(self.Ne_x), np.arange(self.Ne_y))
        X = 1. * X / self.Ne_x
        Y = 1. * Y / self.Ne_y * self.y_dim
        others_e.x = X.ravel()
        others_e.y = Y.ravel()

        self.prefDirs_e = np.ndarray((self.net_Ne, 2))
        for y in xrange(self.Ne_y):
            y_e_norm = float(y) / self.Ne_y * self.y_dim

            for x in xrange(self.Ne_x):
                it = y * self.Ne_x + x

                x_e_norm = float(x) / self.Ne_x
                a.x = x_e_norm
                a.y = y_e_norm

                pd_e = self.getPreferredDirection(x, y)
                self.prefDirs_e[it, :] = pd_e

                pd_norm_e.x = 1. * pd_e[0] / self.Ne_x
                pd_norm_e.y = 1. * pd_e[1] / self.Ne_y * self.y_dim

                tmp_templ = self._generateGaussianWeights(
                    a, others_e, pEE_sigma, pd_norm_e, self.no.prefDirC_ee)

                # tmp_templ down here must be in the proper units (e.g. nS)
                tmp_templ *= g_EE_mean
                tmp_templ[it] = 0.  # do not allow autapses
                self._divergentConnectEE(it, range(self.net_Ne), tmp_templ)

    def _connect_ei_distance(self, AMPA_gaussian, pAMPA_mu, pAMPA_sigma):
        '''Make E-->I connections, according to network options.

        This doc applies to both connect_ei and connect_ie.

        The connections are remapped to [1.0, sqrt(3)/2], whether the topology
        is a twisted torus or just a regular torus.

        Parameters
        ----------
        AMPA_gaussian : bool
            AMPA_gaussian switches between two cases:
                true    Each exciatory neuron has a 2D excitatory gaussian
                        profile, while each inhibitory neuron has a ring-like
                        profile pAMPA_mu, pAMPA_sigma, pGABA_sigma are used,
                        pGABA_mu is discarded
                false   Each excitatory neuron has a ring-like profile, while
                        each inhibitory neuron has a gaussian profile.
                        pAMPA_sigma, pGABA_mu, pGABA_sigma are used, pAMPA_mu
                        is discarded
        '''
        gcnLogger.info('Connecting E-->I (distance-dependent).')
        g_AMPA_mean = self.no.g_AMPA_total / self.net_Ne

        others_e  = Position2D()
        pd_norm_e = Position2D()
        a         = Position2D()

        X, Y = np.meshgrid(np.arange(self.Ni_x), np.arange(self.Ni_y))
        X = 1. * X / self.Ni_x
        Y = 1. * Y / self.Ni_y * self.y_dim
        others_e.x = X.ravel()
        others_e.y = Y.ravel()

        self.prefDirs_e = np.ndarray((self.net_Ne, 2))
        for y in xrange(self.Ne_y):
            y_e_norm = float(y) / self.Ne_y * self.y_dim

            for x in xrange(self.Ne_x):
                it = y * self.Ne_x + x

                x_e_norm = float(x) / self.Ne_x
                a.x = x_e_norm
                a.y = y_e_norm

                pd_e = self.getPreferredDirection(x, y)
                self.prefDirs_e[it, :] = pd_e

                pd_norm_e.x = 1. * pd_e[0] / self.Ni_x
                pd_norm_e.y = 1. * pd_e[1] / self.Ni_y * self.y_dim

                if AMPA_gaussian == 1:
                    tmp_templ = self._generateGaussianWeights(
                        a, others_e, pAMPA_sigma, pd_norm_e,
                        self.no.prefDirC_e)
                elif AMPA_gaussian == 0:
                    tmp_templ = self._generateRinglikeWeights(
                        a, others_e, pAMPA_mu, pAMPA_sigma, pd_norm_e,
                        self.no.prefDirC_e)
                else:
                    raise Exception('AMPA_gaussian parameters must be 0 or 1')


                tmp_templ = self._weight_constructor.generate_weights(
                    tmp_templ, g_AMPA_mean)
                # tmp_templ down here must be in the proper units (e.g. nS)
                self._divergentConnectEI(it, range(self.net_Ni), tmp_templ)

    def _connect_ei_flat(self):
        '''Make E-->I connections that are distance-independent.'''
        gcnLogger.info('Connecting E-->I (flat).')
        g_EI_mean = (self.no.g_AMPA_total / self.net_Ne /
                     self.no.g_EI_uni_density)
        n = int(float(self.net_Ni) * self.no.g_EI_uni_density)
        self._randomDivergentConnectEI(range(self.net_Ne),
                                       range(self.net_Ni),
                                       n,
                                       g_EI_mean)

    def _connect_ie_distance(self, AMPA_gaussian, pGABA_mu, pGABA_sigma):
        '''Make I-->E connections, according to network options.

        This doc applies to both connect_ei and connect_ie.

        The connections are remapped to [1.0, sqrt(3)/2], whether the topology
        is a twisted torus or just a regular torus.

        Parameters
        ----------
        AMPA_gaussian : bool
            AMPA_gaussian switches between two cases:
                true    Each exciatory neuron has a 2D excitatory gaussian
                        profile, while each inhibitory neuron has a ring-like
                        profile pAMPA_mu, pAMPA_sigma, pGABA_sigma are used,
                        pGABA_mu is discarded
                false   Each excitatory neuron has a ring-like profile, while
                        each inhibitory neuron has a gaussian profile.
                        pAMPA_sigma, pGABA_mu, pGABA_sigma are used, pAMPA_mu
                        is discarded
        '''
        gcnLogger.info('Connecting I-->E (distance-dependent).')
        g_GABA_mean = self.no.g_GABA_total / self.net_Ni
        g_uni_GABA_total = self.no.g_GABA_total * self.no.g_uni_GABA_frac
        g_uni_GABA_mean = (g_uni_GABA_total / self.net_Ni /
                           self.no.uni_GABA_density)
        print("g_uni_GABA_total: ", g_uni_GABA_total)
        print("g_uni_GABA_mean: ", g_uni_GABA_mean)

        others_i  = Position2D()
        pd_norm_i = Position2D()
        a         = Position2D()

        X, Y = np.meshgrid(np.arange(self.Ne_x), np.arange(self.Ne_y))
        X = 1. * X / self.Ne_x
        Y = 1. * Y / self.Ne_y * self.y_dim
        others_i.x = X.ravel()
        others_i.y = Y.ravel()

        conn_th = 1e-5
        self.prefDirs_i = np.ndarray((self.net_Ni, 2))
        for y in xrange(self.Ni_y):
            y_i_norm = float(y) / self.Ni_y * self.y_dim
            for x in xrange(self.Ni_x):
                it = y * self.Ni_x + x
                x_i_norm = float(x) / self.Ni_x

                a.x = x_i_norm
                a.y = y_i_norm

                pd_i = self.getPreferredDirection(x, y)
                self.prefDirs_i[it, :] = pd_i

                pd_norm_i.x = 1. * pd_i[0] / self.Ne_x
                pd_norm_i.y = 1. * pd_i[1] / self.Ne_y * self.y_dim

                if AMPA_gaussian == 1:
                    tmp_templ = self._generateRinglikeWeights(
                        a, others_i, pGABA_mu, pGABA_sigma, pd_norm_i,
                        self.no.prefDirC_i)
                elif AMPA_gaussian == 0:
                    tmp_templ = self._generateGaussianWeights(
                        a, others_i, pGABA_sigma, pd_norm_i,
                        self.no.prefDirC_i)
                else:
                    raise Exception('AMPA_gaussian parameters must be 0 or 1')

                # FIXME: ugly: B_GABA is defined only in child classes
                tmp_templ = self._weight_constructor.generate_weights(
                    tmp_templ, self.B_GABA * g_GABA_mean)
                self._addToConnections(
                    tmp_templ, self.no.uni_GABA_density * 100.0,
                    g_uni_GABA_mean)
                E_nid = (tmp_templ > conn_th).nonzero()[0]
                self._divergentConnectIE(it, E_nid, tmp_templ[E_nid])

    def _connect_ie_flat(self):
        '''Make I-->E connections that are distance independent.'''
        gcnLogger.info('Connecting I-->E (flat).')
        g_IE_mean = (self.no.g_GABA_total / self.net_Ni /
                     self.no.g_IE_uni_density)
        n = int(float(self.net_Ne) * self.no.g_IE_uni_density)
        self._randomDivergentConnectIE(range(self.net_Ni),
                                       range(self.net_Ne),
                                       n,
                                       g_IE_mean)

    def _connect_ii_flat(self):
        '''Make I-->I connections that are distance independent.'''
        gcnLogger.info('Connecting I-->I (flat).')
        g_II_mean = (self.no.g_II_total / self.net_Ni /
                     self.no.g_II_uni_density)
        gcnLogger.debug('g_II_total: %f, g_II_mean: %f', self.no.g_II_total,
                        g_II_mean)
        n = int(float(self.net_Ne) * self.no.g_IE_uni_density)
        self._randomDivergentConnectII(range(self.net_Ni),
                                       range(self.net_Ni),
                                       n,
                                       g_II_mean)

    ###########################################################################
    #                     External sources definitions
    ###########################################################################
    def setVelocityCurrentInput_e(self, prefDirs_mask=None):
        '''
        Setup a velocity input to the excitatory population. Current input.
        '''
        raise NotImplementedError()

    def setVelocityCurrentInput_i(self):
        '''
        Setup a velocity input to the inhibitory population. Current input.
        '''
        raise NotImplementedError()

    def setConstantVelocityCurrent_e(self, vel):
        '''
        Setup a constant velocity current onto E poputlaion, where vel must be
        a list of numbers:
            vel = [vel_x, vel_y]
        '''
        raise NotImplementedError()

    def setConstantVelocityCurrent_i(self, vel):
        '''
        Setup a constant velocity current onto I population, where vel must be
        a list of numbers:
            vel = [vel_x, vel_y]
        '''
        raise NotImplementedError()

    def getAttrDictionary(self):
        '''
        Get a dictionary containing all the necessary attributes the user might
        need in order to work with data produced by the simulation.
        '''
        raise NotImplementedError()

    ###########################################################################
    #                           Other
    ###########################################################################
    def beginConstruction(self):
        '''
        Mark the beginning of network construction.
        '''
        self._constrStartT = time.time()
        print("Starting network construction")

    def endConstruction(self):
        '''
        Mark the end of network construction.
        '''
        self._constrEndT = time.time()
        print("Network construction finished.")

    def constructionTime(self):
        '''Compute network construction time'''
        assert self._constrStartT is not None
        if self._constrEndT is None:
            raise RuntimeError("Cannot compute contruction time. End time has "
                               "not been marked yet.")
        else:
            return self._constrEndT - self._constrStartT

    def beginSimulation(self):
        '''Mark beginning of the simulation'''
        self._simStartT = time.time()
        print("Simulation has started...")

    def endSimulation(self):
        '''Mark end of the simulation'''
        self._simEndT = time.time()
        print("Simulation finished")

    def simulationTime(self):
        '''Compute simulation time'''
        assert self._simStartT is not None
        if self._simEndT is None:
            raise RuntimeError("Cannot compute simulation time. End time has "
                               "not been marked yet (no simulation has been "
                               "run?).")
        else:
            return self._simEndT - self._simStartT

    def totalTime(self):
        '''
        Return elapsed time in seconds, since the network construction start.
        '''
        return time.time() - self._startT

    def _getTimes(self):
        '''Get simulation times'''
        return (self.constructionTime(), self.simulationTime(),
                self.totalTime())

    def printTimes(self, constrT=True, simT=True, totalT=True):
        '''Print the different elapsed simulation times.'''
        constrT, simT, totalT = self._getTimes()
        print("Timer statistics:")
        print("   Construction: {0} s".format(constrT))
        print("   Simulation  : {0} s".format(simT))
        print("   Total       : {0} s".format(totalT))
        return constrT, simT, totalT

    def getPreferredDirection(self, pos_x, pos_y):
        '''
        Get a preferred direction for a neuron.

        Parameters
        ----------
        pos_x/y : int
            Position of neuron in 2d sheet
        '''
        pos4_x = pos_x % 2
        pos2_y = pos_y % 2
        if pos4_x == 0:
            if pos2_y == 0:
                return [-1, 0]  # Left
            else:
                return [0, -1]  # Down
        else:
            if pos2_y == 0:
                return [0, 1]  # up
            else:
                return [1, 0]  # Right
