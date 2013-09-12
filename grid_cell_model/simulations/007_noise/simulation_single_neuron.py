#
#   simulation_single_neuron.py
#
#   Main simulation run: Single neuron simulations (figure1)
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
from models.parameters       import getOptParser
from models.gc_single_neuron import OneNeuronNetwork


parser          = getOptParser()
(options, args) = parser.parse_args()


################################################################################
ei_net = OneNeuronNetwork(options, simulationOpts=None)

ei_net.simulate(options.time, printTime=options.printTime)
ei_net.saveData()

################################################################################

