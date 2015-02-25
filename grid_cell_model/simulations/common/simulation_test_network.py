'''Main simulation run: Run a test simulation run of a bump with constant speed.

This simulation will overwrite any old files that are present in the output
directory.
'''
from __future__ import absolute_import, print_function, division

from nest.hl_api import NESTError

from grid_cell_model.models.parameters import getOptParser
from grid_cell_model.models.gc_net_nest import ConstantVelocityNetwork
from grid_cell_model.models.seeds import TrialSeedGenerator
from grid_cell_model.data_storage import DataStorage

parser = getOptParser()
(o, args) = parser.parse_args()

output_fname = "{0}/{1}job{2:05}_output.h5".format(o.output_dir,
                                                   o.fileNamePrefix, o.job_num)
d = DataStorage.open(output_fname, 'w')
d['trials'] = []

overalT = 0.
stop = False
###############################################################################
seed_gen = TrialSeedGenerator(o.master_seed)
for trial_idx in range(o.ntrials):
    seed_gen.set_generators(trial_idx)  # Each trial is reproducible
    const_v = [0.0, -o.Ivel]
    ei_net = ConstantVelocityNetwork(o, simulationOpts=None, vel=const_v)
    d['net_params'] = ei_net.getNetParams()  # Common settings will stay

    try:
        ei_net.simulate(o.time, printTime=o.printTime)
    except NESTError as e:
        print("Simulation interrupted. Message: {0}".format(str(e)))
        print("Not saving the data. Trying to clean up if possible...")
        stop = True
    ei_net.endSimulation()
    d['trials'].append(ei_net.getAllData())
    d.flush()
    constrT, simT, totalT = ei_net.printTimes()
    overalT += totalT
    if stop:
        break

d.close()
print("Script total run time: {0} s".format(overalT))
###############################################################################
