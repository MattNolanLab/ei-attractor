'''Main simulation run: Bump velocity estimation.'''
from __future__ import absolute_import, print_function, division

import numpy as np
from os.path import exists
from nest.hl_api        import NESTError

from grid_cell_model.models.parameters  import getOptParser
from grid_cell_model.models.gc_net_nest import ConstantVelocityNetwork
from grid_cell_model.models.seeds import TrialSeedGenerator
from grid_cell_model.data_storage       import DataStorage


parser = getOptParser()
parser.add_argument("--IvelMax", type=float, required=True, help="Max constant velocity current input (pA)")
parser.add_argument("--dIvel",   type=float, required=True, help="Constant velocity current input step (pA)")
parser.add_argument("--ispikes", type=int,   choices=[0, 1], default=0, help="Whether to save spikes from the I population")

(o, args) = parser.parse_args()


output_fname = "{0}/{1}job{2:05}_output.h5".format(o.output_dir,
        o.fileNamePrefix, o.job_num)
d = DataStorage.open(output_fname, 'a')
if ("trials" not in d.keys()):
    d['trials'] = []

# Initialise seeds and check their consistency with previously saved data
seed_gen = TrialSeedGenerator(o.master_seed)
if len(d['trials']) == 0:
    d['master_seed'] = o.master_seed
else:
    try:
        seed_gen.check_master_seed(d['master_seed'], o.master_seed)
    except ValueError as e:
        d.close()
        raise e

overalT = 0.
oldNTrials = len(d['trials'])
################################################################################
for trial_idx in range(o.ntrials):
    print("\n\t\tStarting/appending to trial no. {0}\n".format(trial_idx))
    if trial_idx >= oldNTrials:  # Create new trial
        d['trials'].append({})
    trialOut = d['trials'][trial_idx]

    # Now check if there is data in the trial and append
    if 'IvelVec' not in trialOut:
        oldNIvel = 0
        trialOut['IvelData'] = []
    else:
        oldNIvel = len(trialOut['IvelVec'])

    try:
        IvelVecAppend = np.arange(oldNIvel*o.dIvel, o.IvelMax + o.dIvel, o.dIvel)
        for Ivel in IvelVecAppend:
            seed_gen.set_generators(trial_idx)  # Each trial is reproducible
            const_v = [0.0, -Ivel]
            ei_net = ConstantVelocityNetwork(o, simulationOpts=None, vel=const_v)

            ei_net.simulate(o.time, printTime=o.printTime)
            ei_net.endSimulation()
            trialOut['IvelData'].append(ei_net.getMinimalSaveData(ispikes=o.ispikes))
            d.flush()
            constrT, simT, totalT = ei_net.printTimes()
            overalT += totalT
        trialOut['IvelVec'] = np.arange(.0, o.IvelMax + o.dIvel, o.dIvel)
        d.flush()
    except NESTError as e:
        print("Simulation interrupted. Message: {0}".format(str(e)))
        print("Not saving the last trial. Trying to clean up if possible...")
        break

d.close()
print("Script total run time: {0} s".format(overalT))
################################################################################

