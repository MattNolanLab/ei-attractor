#
#   simulation_connections.py
#
#   Main simulation run: Only export E and I connections
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
from numpy.random       import choice

from models.parameters  import getOptParser
from models.gc_net_nest import BasicGridCellNetwork
from data_storage       import DataStorage


parser          = getOptParser()
(o, args) = parser.parse_args()


out = []
overalT = 0.
################################################################################
for trial_idx in range(o.ntrials):
    print("\n\t\tStarting trial no. {0}\n".format(trial_idx))
    ei_net = BasicGridCellNetwork(o, simulationOpts=None)
    
    ei_net.endConstruction()
    ei_net.beginSimulation() 

    data = ei_net.getNetParams()
    # E --> I neurons
    data['g_IE'] = ei_net.getConnMatrix("E")
    # I --> E neurons
    data['g_EI'] = ei_net.getConnMatrix("I")

    ei_net.endSimulation()

    
    out.append(data)
    constrT, simT, totalT = ei_net.printTimes()
    overalT += totalT

output_fname = "{0}/{1}job{2:05}_output.h5".format(o.output_dir,
        o.fileNamePrefix, o.job_num)
d = DataStorage.open(output_fname, 'w')
d["trials"] = out
d.close()
print "Script total run time: {0} s".format(overalT)
################################################################################

