#
#   simulation_grids.py
#
#   Main simulation run: Simulation of an animal movement (grid fields)
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
from os.path            import exists
from numpy.random       import choice

from models.parameters  import getOptParser
from models.gc_net_nest import BasicGridCellNetwork
from data_storage       import DataStorage
from nest.hl_api        import NESTError


parser          = getOptParser()
parser.add_option("--velON", type="int", help="Velocity input ON?")
(options, args) = parser.parse_args()


# Do nothing when bumpCurrentSlope is NaN
if (np.isnan(options.bumpCurrentSlope)):
    print('bumpCurrentSlope is NaN. Not performing the simulation')
    exit(0)
else:
    print("bumpCurrentSlope: {0}".format(options.bumpCurrentSlope))


stateMonParams = {
        'start' : options.time - options.stateMonDur
}
nrec_spikes_e = None # all neurons
nrec_spikes_i = 10


output_fname = "{0}/{1}job{2:05}_output.h5".format(options.output_dir,
        options.fileNamePrefix, options.job_num)
if (exists(output_fname)):
    print('Output file exists, doing nothing')
    exit(0)
d = DataStorage.open(output_fname, 'w')

out = []
overalT = 0.
################################################################################
for trial_idx in range(options.ntrials):
    print("\n\t\tStarting trial no. {0}\n".format(trial_idx))
    ei_net = BasicGridCellNetwork(options, simulationOpts=None,
            nrec_spikes=(nrec_spikes_e, nrec_spikes_i),
            stateRecParams=(stateMonParams, stateMonParams))

    if (options.velON == 1):
        ei_net.setVelocityCurrentInput_e()
    ei_net.setPlaceCells()
    
    try:
        ei_net.simulate(options.time, printTime=options.printTime)
    except NESTError as e:
        print("Simulation interrupted. Message: {0}".format(str(e)))
        print("Trying to save the simulated data if possible...")
    ei_net.endSimulation()
    out.append(ei_net.getAllData())
    constrT, simT, totalT = ei_net.printTimes()
    overalT += totalT

d["trials"] = out
d.close()
print "Script total run time: {0} s".format(overalT)
################################################################################


