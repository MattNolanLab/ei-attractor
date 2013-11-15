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
parser.add_option("--IvelMax", type="float", help="Max constant velocity current input (pA)")
parser.add_option("--dIvel",   type="float", help="Constant velocity current input step (pA)")
parser.add_option("--ispikes", type="int",   help="Whether to save spikes from the I population")

(options, args) = parser.parse_args()


output_fname = "{0}/{1}job{2:05}_output.h5".format(options.output_dir,
        options.fileNamePrefix, options.job_num)
d = DataStorage.open(output_fname, 'a')
if ("trials" not in d.keys()):
    d['trials'] = []

overalT = 0.
################################################################################
for trial_idx in range(len(d['trials']), options.ntrials):
    print("\n\t\tStarting trial no. {0}\n".format(trial_idx))
    try:
        trialOut = {}
        IvelVec = np.arange(0.0, options.IvelMax + options.dIvel, options.dIvel)
        trialOut['IvelVec'] = IvelVec
        trialOut['IvelData'] = []
        for Ivel in IvelVec:
            const_v = [Ivel, 0.0]
            ei_net = ConstantVelocityNetwork(options, simulationOpts=None, vel=const_v)

            ei_net.simulate(options.time, printTime=options.printTime)
            ei_net.endSimulation()
            trialOut['IvelData'].append(ei_net.getMinimalSaveData(ispikes=options.ispikes))
            constrT, simT, totalT = ei_net.printTimes()
            overalT += totalT
        d['trials'].append(trialOut)
        d.flush()
    except NESTError as e:
        print("Simulation interrupted. Message: {0}".format(str(e)))
        print("Not saving the last trial. Trying to clean up if possible...")
        break

d.close()
print "Script total run time: {0} s".format(overalT)
################################################################################

