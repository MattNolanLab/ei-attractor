#
#   simulation_velocity.py
#
#   Main simulation run: Bump velocity estimation.
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
import numpy as np
from os.path import exists
from numpy.random       import choice
from nest.hl_api        import NESTError

from models.parameters  import getOptParser
from models.gc_net_nest import ConstantVelocityNetwork
from data_storage       import DataStorage


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

overalT = 0.
oldNTrials = len(d['trials'])
################################################################################
for trial_idx in range(o.ntrials):
    print("\n\t\tStarting/appending to trial no. {0}\n".format(trial_idx))
    if trial_idx < oldNTrials: # Trial exists
        trialOut = d['trials'][trial_idx]
    else:
        trialOut = {}

    # Now check if there is data in the trial and append
    if 'IvelVec' not in trialOut:
        oldNIvel = 0
        trialOut['IvelData'] = []
    else:
        oldNIvel = len(trialOut['IvelVec'])

    try:
        IvelVecAppend = np.arange(oldNIvel*o.dIvel, o.IvelMax + o.dIvel, o.dIvel)
        for Ivel in IvelVecAppend:
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
print "Script total run time: {0} s".format(overalT)
################################################################################

