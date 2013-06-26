#
#   simulation_bump_fitting.py
#
#   Main simulation run: Fitting a Gaussian to the bump and frequency analysis.
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
from numpy.random       import choice

from models.parameters  import getOptParser
from models.gc_net_nest import ConstantVelocityNetwork
from data_storage       import DataStorage


parser          = getOptParser()
parser.add_option("--IvelMax",      type="float", help="Max constant velocity current input (pA)")
parser.add_option("--dIvel",        type="float", help="Constant velocity current input step (pA)")

(options, args) = parser.parse_args()


out = []
overalT = 0.
################################################################################
for trial_idx in range(options.ntrials):
    print("\n\t\tStarting trial no. {0}\n".format(trial_idx))
    trialOut = {}
    IvelVec = np.arange(0.0, options.IvelMax + options.dIvel, options.dIvel)
    trialOut['IvelVec'] = IvelVec
    trialOut['IvelData'] = []
    for Ivel in IvelVec:
        const_v = [Ivel, 0.0]
        ei_net = ConstantVelocityNetwork(options, simulationOpts=None, vel=const_v)

        ei_net.simulate(options.time, printTime=options.printTime)
        ei_net.endSimulation()
        trialOut['IvelData'].append(ei_net.getSpikes())
        constrT, simT, totalT = ei_net.printTimes()
        overalT += totalT
    out.append(trialOut)

output_fname = "{0}/{1}job{2:05}_output.h5".format(options.output_dir,
        options.fileNamePrefix, options.job_num)
d = DataStorage.open(output_fname, 'w')
d["trials"] = out
d.close()
print "Script total run time: {0} s".format(overalT)
################################################################################

