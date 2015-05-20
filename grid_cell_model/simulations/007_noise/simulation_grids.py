'''
Simulation run with place cells
Depending on the parameters, this simulation performs either the full grid field
simulation with initialisation and velocity place cells active (velOn==True and
constantPosition=False). Or a simulation in which the animal holds still at a
specified position (constantPosition==True).
'''
from __future__ import absolute_import, print_function, division

import numpy as np
from nest.hl_api        import NESTError

from grid_cell_model.models.parameters  import getOptParser
from grid_cell_model.models.gc_net_nest import (BasicGridCellNetwork,
                                                ConstPosInputs)
from grid_cell_model.models.seeds import TrialSeedGenerator
from grid_cell_model.data_storage import DataStorage

import logging
logger = logging.getLogger(__name__)

parser = getOptParser()
parser.add_argument("--velON",            type=int,   choices=[0, 1], required=True, help="Velocity input ON?")
parser.add_argument("--pcON",             type=int,   choices=[0, 1], required=True, help="Place cell input ON?")
parser.add_argument("--constantPosition", type=int,   choices=[0, 1], required=True, help="Should the animat move?")
parser.add_argument("--staticPos_x",      type=float, default=0.0,    help="Static position X coordinate")
parser.add_argument("--staticPos_y",      type=float, default=0.0,    help="Static position Y coordinate")
parser.add_argument("--nrec_spikes_i",    type=int,   default=10,     help="Number of I cell to record spikes from. Chosen randomly.")
parser.add_argument("--rec_spikes_probabilistic", type=int, choices=[0, 1], default=0)
o, _ = parser.parse_args()

# Do nothing when bumpCurrentSlope is NaN
if (np.isnan(o.bumpCurrentSlope)):
    logger.info('bumpCurrentSlope is NaN. Not performing the simulation')
    exit(0)
else:
    logger.info("bumpCurrentSlope: {0}".format(o.bumpCurrentSlope))


stateMonParams = {
        'start' : o.time - o.stateMonDur
}
nrec_spikes_e = None # all neurons
nrec_spikes_i = o.nrec_spikes_i

output_fname = "{0}/{1}job{2:05}_output.h5".format(o.output_dir,
        o.fileNamePrefix, o.job_num)
d = DataStorage.open(output_fname, 'a')
if ("trials" not in d.keys()):
    d['trials'] = []

seed_gen = TrialSeedGenerator(int(o.master_seed))
if len(d['trials']) == 0:
    d['master_seed'] = int(o.master_seed)
else:
    try:
        seed_gen.check_master_seed(d['master_seed'], int(o.master_seed))
    except ValueError as e:
        d.close()
        raise e

overalT = 0.
stop = False
###############################################################################
for trial_idx in range(len(d['trials']), o.ntrials):
    print("\n\t\tStarting trial no. {0}\n".format(trial_idx))
    seed_gen.set_generators(trial_idx)
    d['invalidated'] = 1
    ei_net = BasicGridCellNetwork(o, simulationOpts=None,
                                  nrec_spikes=(nrec_spikes_e, nrec_spikes_i),
                                  stateRecParams=(stateMonParams,
                                                  stateMonParams),
                                  rec_spikes_probabilistic=o.rec_spikes_probabilistic)

    if o.velON and not o.constantPosition:
        ei_net.setVelocityCurrentInput_e()
    if o.pcON:
        # This also sets the start PCs
        posIn = ConstPosInputs(0, 0) if o.constantPosition else None
        ei_net.setPlaceCells(posIn=posIn)
    else:
        # Here the start PCs must be set explicitly
        ei_net.setStartPlaceCells(ConstPosInputs(0, 0))
    if o.ipc_ON:
        if o.constantPosition:
            raise RuntimeError("Place cells connected to I cells cannot be "
                               "used when the constantPosition parameter is "
                               "ON.")
        ei_net.setIPlaceCells()

    d['net_params'] = ei_net.getNetParams()  # Common settings will stay
    d.flush()

    try:
        ei_net.simulate(o.time, printTime=o.printTime)
    except NESTError as e:
        print("Simulation interrupted. Message: {0}".format(str(e)))
        print("Trying to save the simulated data if possible...")
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
################################################################################


