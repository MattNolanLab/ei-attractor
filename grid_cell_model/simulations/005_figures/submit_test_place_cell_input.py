#
#   submit_test_place_cell_input.py
#
#   Submit job(s) to the cluster/workstation
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

from default_params import defaultParameters
from common         import *

import logging as lg


lg.basicConfig(level=lg.DEBUG)


EDDIE = False  # if eddie, submit on a cluster using qsub


parameters = defaultParameters

parameters['time']              = 6e3  # ms
parameters['ndumps']            = 1

#parameters['placeT']            = 1e3      # ms

parameters['stateMonDur']       = parameters['time']

parameters['bumpCurrentSlope']  = 1.05     # pA/(cm/s), !! this will depend on prefDirC !!
parameters['gridSep']           = 70        # cm, grid field inter-peak distance

parameters['theta_noise_sigma'] = 0         # pA

#parameters['Ivel']              = 40        # pA

parameters['output_dir']        = 'output'

startJobNum = 2
numRepeat = 1

# Workstation parameters
programName         = 'python2.6 simulation_test_place_cell_input.py'
blocking            = False

# Cluster parameters
eddie_scriptName    = 'eddie_submit.sh simulation_test_place_cell_input.py'
qsub_params         = "-P inf_ndtc -cwd -j y -l h_rt=13:00:00 -pe memory-2G 2"
qsub_output_dir     = parameters['output_dir']

ac = ArgumentCreator(parameters)

iterparams = {
#        'bumpCurrentSlope'  : [1.15, 1.175, 1.2]
        'placeT'    : [1e3]
}
ac.insertDict(iterparams, mult=False)

if EDDIE:
    submitter = QsubSubmitter(ac, eddie_scriptName, qsub_params, qsub_output_dir)
else:
    submitter = GenericSubmitter(ac, programName, blocking=blocking)
submitter.submitAll(startJobNum, numRepeat, dry_run=False)
