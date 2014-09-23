'''Main simulation run: Single neuron simulations (figure1).'''
from __future__ import absolute_import, print_function, division

from grid_cell_model.models.parameters       import getOptParser
from grid_cell_model.models.gc_single_neuron import OneNeuronNetwork


parser          = getOptParser()
(options, args) = parser.parse_args()


################################################################################
ei_net = OneNeuronNetwork(options, simulationOpts=None)

ei_net.simulate(options.time, printTime=options.printTime)
ei_net.saveData()

################################################################################

