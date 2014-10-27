'''Main simulation run: Only export E and I connections.'''
from __future__ import absolute_import, print_function, division

from numpy.random import choice

from grid_cell_model.models.parameters  import getOptParser
from grid_cell_model.models.gc_net_nest import BasicGridCellNetwork
from grid_cell_model.data_storage       import DataStorage


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

