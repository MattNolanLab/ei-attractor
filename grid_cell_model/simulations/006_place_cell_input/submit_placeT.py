#
#   submit_basic_grids.py
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
from submitting.submitters         import *

import logging as lg


lg.basicConfig(level=lg.DEBUG)


CLUSTER = True  # if True, submit on a cluster using qsub


parameters = defaultParameters

parameters['time']              = 1199.9e3  # ms
#parameters['time']              = 3e3  # ms
parameters['ndumps']            = 10

parameters['prefDirC_e']        = 4
parameters['prefDirC_i']        = 0

#parameters['placeT']            = 10e3      # ms
parameters['placeDur']          = 125       # ms

parameters['bumpCurrentSlope']  = 1.05      # pA/(cm/s), !! this will depend on prefDirC !!
parameters['gridSep']           = 60        # cm, grid field inter-peak distance
parameters['theta_noise_sigma'] = 0         # pA

parameters['output_dir']        = 'output'


startJobNum =100
numRepeat = 5

# Workstation parameters
programName         = 'python2.6 simulation_placeT.py'
blocking            = False

# Cluster parameters
cluster_scriptName  = 'eddie_submit.sh simulation_placeT.py'
qsub_params         = "-P inf_ndtc -cwd -j y -l h_rt=14:00:00 -pe memory-2G 2"
qsub_output_dir     = parameters['output_dir']

ac = ArgumentCreator(parameters)

iterparams = {
#    'Iplace'    :   [50, 100, 150, 200, 250]
#    'gridSep' : [50, 60]
    'placeT'    :  [125, 250, 375, 500, 625,750, 875, 1e3, 2e3, 3e3, 4e3, 5e3, 6e3, 7e3, 8e3, 9e3, 10e3, 12e3, 14e3, 16e3, 18e3, 20e3, 30e3, 40e3, 50e3, 60e3, 120e3, 180e3, 240e3, 300e3]
}
ac.insertDict(iterparams, mult=False)

if CLUSTER:
    submitter = QsubSubmitter(ac, cluster_scriptName, qsub_params, qsub_output_dir)
else:
    submitter = GenericSubmitter(ac, programName, blocking=blocking)
submitter.submitAll(startJobNum, numRepeat, dry_run=False)
