#
#   simulation_stationary.py
#
#   Main simulation run: Simulation of a stationary bump.
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
from nest.hl_api        import NESTError

from models.parameters  import getOptParser
from models.gc_net_nest import BasicGridCellNetwork
from data_storage       import DataStorage


parser          = getOptParser()
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
    d['invalidated'] = 1
    try:
        ei_net = BasicGridCellNetwork(options, simulationOpts=None)
        
        const_v = [0.0, 0.0]
        ei_net.setConstantVelocityCurrent_e(const_v)
        
        ei_net.simulate(options.time, printTime=options.printTime)
        ei_net.endSimulation()
        d['trials'].append(ei_net.getAllData())
        d.flush()
        constrT, simT, totalT = ei_net.printTimes()
        overalT += totalT
    except NESTError as e:
        print("Simulation interrupted. Message: {0}".format(str(e)))
        print("Trying to save the simulated data if possible...")
        break

d.close()
print "Script total run time: {0} s".format(overalT)
################################################################################
