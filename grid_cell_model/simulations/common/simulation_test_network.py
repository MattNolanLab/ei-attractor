'''Main simulation run: Run network with constant speed.'''
from __future__ import absolute_import, print_function, division

from nest.hl_api import NESTError

from grid_cell_model.models.parameters import getOptParser
from grid_cell_model.models.gc_net_nest import ConstantVelocityNetwork
from grid_cell_model.models.seeds import TrialSeedGenerator
from grid_cell_model.data_storage import DataStorage

parser = getOptParser()
parser.add_argument("--trial", type=int, required=True,
                    help='Trial number.')
parser.add_argument("--ispikes", type=int, choices=[0, 1], default=0,
                    help="Whether to save spikes from the I population")
(o, args) = parser.parse_args()

output_fname = "{0}/{1}job{2:05}_output.h5".format(o.output_dir,
                                                   o.fileNamePrefix, o.job_num)
d = DataStorage.open(output_fname, 'w')

overalT = 0.
###############################################################################
try:
    seed_gen = TrialSeedGenerator(o.master_seed)
    seed_gen.set_generators(o.trial)  # Each trial is reproducible
    const_v = [0.0, -o.Ivel]
    ei_net = ConstantVelocityNetwork(o, simulationOpts=None, vel=const_v)

    ei_net.simulate(o.time, printTime=o.printTime)
    ei_net.endSimulation()
    d['Data'] = ei_net.getMinimalSaveData(ispikes=o.ispikes)
    d.flush()
    constrT, simT, totalT = ei_net.printTimes()
    overalT += totalT
except NESTError as e:
    print("Simulation interrupted. Message: {0}".format(str(e)))
    print("Not saving the data. Trying to clean up if possible...")

d.close()
print("Script total run time: {0} s".format(overalT))
###############################################################################
