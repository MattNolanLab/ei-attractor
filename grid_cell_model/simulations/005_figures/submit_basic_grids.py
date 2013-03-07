#
#   submit_job.py
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
from common         import *

import logging as lg


lg.basicConfig(level=lg.DEBUG)


CLUSTER = False  # if true, submit on a cluster using qsub


parameters = defaultParameters

parameters['time']              = 5e3      # ms
parameters['delay']             = 0.1       # ms

parameters['prefDirC_e']        = 4
parameters['prefDirC_i']        = 0

parameters['theta_noise_sigma'] = 0          # pA
parameters['noise_sigma']       = 150.       # pA

parameters['Ne']                = 34
parameters['Ni']                = 34
parameters['N_place_cells']     = 30*30
parameters['gridSep']           = 60.0      # cm, grid field inter-peak distance
parameters['bumpCurrentSlope']  = 0.53      # pA/(cm/s), !! this will depend on prefDirC !!

#parameters['pc_max_rate']       = 100.0     # Hz

parameters['output_dir']        = 'output_local'
parameters['nthreads']          = 8
parameters['ndumps']            = 1

startJobNum = 0
numRepeat = 1

# Workstation parameters
programName         = 'python2.6 simulation_basic_grids.py'
blocking            = False

# Cluster parameters
eddie_scriptName    = 'cluster_submit.sh simulation_basic_grids.py'
qsub_params         = "-P inf_ndtc -cwd -j y -l h_rt=13:00:00 -pe memory-2G 2"
qsub_output_dir     = parameters['output_dir']

ac = ArgumentCreator(parameters)

iterparams = {
        'pc_max_rate'   : np.arange(0, 90, 10)
}
ac.insertDict(iterparams, mult=False)

if CLUSTER:
    submitter = QsubSubmitter(ac, eddie_scriptName, qsub_params, qsub_output_dir)
else:
    submitter = GenericSubmitter(ac, programName, blocking=blocking)
submitter.submitAll(startJobNum, numRepeat, dry_run=False)
