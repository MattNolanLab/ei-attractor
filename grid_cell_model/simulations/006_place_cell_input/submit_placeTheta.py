#
#   submit_placeTheta.py
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
import numpy as np

from default_params import defaultParameters
from submitting.submitters         import *

import logging as lg


lg.basicConfig(level=lg.DEBUG)


CLUSTER = True  # if True, submit on a cluster using qsub


parameters = defaultParameters

parameters['time']              = 1199.9e3  # ms
#parameters['time']              = 2e3  # ms
parameters['ndumps']            = 10

parameters['thetaPlaceFreq']    = 8 # Hz
#parameters['thetaPlacePhase']   = 0 # rad

parameters['Iplace']            = 10 # pA

parameters['bumpCurrentSlope']  = 1.05      # pA/(cm/s), !! this will depend on prefDirC !!
parameters['gridSep']           = 60        # cm, grid field inter-peak distance
parameters['theta_noise_sigma'] = 0         # pA

parameters['output_dir']        = 'output'

startJobNum = 2100
numRepeat = 10

# Workstation parameters
programName         = 'python2.6 simulation_placeTheta.py'
blocking            = False

# Cluster parameters
cluster_scriptName  = 'eddie_submit.sh simulation_placeTheta.py'
qsub_params         = "-P inf_ndtc -cwd -j y -l h_rt=14:00:00 -pe memory-2G 2"
qsub_output_dir     = parameters['output_dir']

ac = ArgumentCreator(parameters)

iterparams = {
    "thetaPlacePhase" : np.arange(0, np.pi+np.pi/10, np.pi/10)
#    'Iplace'    :   range(0, 50, 10)
#    'gridSep' : [50, 60]
}
ac.insertDict(iterparams, mult=False)

if CLUSTER:
    submitter = QsubSubmitter(ac, cluster_scriptName, qsub_params, qsub_output_dir)
else:
    submitter = GenericSubmitter(ac, programName, blocking=blocking)
submitter.submitAll(startJobNum, numRepeat, dry_run=False)
