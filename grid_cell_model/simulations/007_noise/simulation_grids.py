'''
Simulation run with place cells
Depending on the parameters, this simulation performs either the full grid field
simulation with initialisation and velocity place cells active (velOn==True and
constantPosition=False). Or a simulation in which the animal holds still at a
specified position (constantPosition==True).
'''
import numpy as np
from os.path            import exists
from numpy.random       import choice

from models.parameters  import getOptParser
from models.gc_net_nest import BasicGridCellNetwork, ConstPosInputs
from data_storage       import DataStorage
from nest.hl_api        import NESTError

import logging
logger = logging.getLogger(__name__)

parser = getOptParser()
parser.add_argument("--velON",            type=int,   choices=[0, 1], help="Velocity input ON?")
parser.add_argument("--constantPosition", type=int,   choices=[0, 1], help="Should the animat move?")
parser.add_argument("--staticPos_x",      type=float, default=0.0,    help="Static position X coordinate")
parser.add_argument("--staticPos_y",      type=float, default=0.0,    help="Static position Y coordinate")
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
nrec_spikes_i = 10


output_fname = "{0}/{1}job{2:05}_output.h5".format(o.output_dir,
        o.fileNamePrefix, o.job_num)
d = DataStorage.open(output_fname, 'a')
if ("trials" not in d.keys()):
    d['trials'] = []

overalT = 0.
################################################################################
for trial_idx in range(len(d['trials']), o.ntrials):
    print("\n\t\tStarting trial no. {0}\n".format(trial_idx))
    d['invalidated'] = 1
    ei_net = BasicGridCellNetwork(o, simulationOpts=None,
            nrec_spikes=(nrec_spikes_e, nrec_spikes_i),
            stateRecParams=(stateMonParams, stateMonParams))
    if o.velON and not o.constantPosition:
        ei_net.setVelocityCurrentInput_e()
    posIn = ConstPosInputs(0, 0) if o.constantPosition else None
    ei_net.setPlaceCells(posIn=posIn)
    
    try:
        ei_net.simulate(o.time, printTime=o.printTime)
    except NESTError as e:
        print("Simulation interrupted. Message: {0}".format(str(e)))
        print("Trying to save the simulated data if possible...")
    ei_net.endSimulation()
    d['trials'].append(ei_net.getAllData())
    d.flush()
    constrT, simT, totalT = ei_net.printTimes()
    overalT += totalT

d.close()
print "Script total run time: {0} s".format(overalT)
################################################################################


