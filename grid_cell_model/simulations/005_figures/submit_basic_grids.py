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
import numpy as np
from scipy.io       import savemat

from default_params import defaultParameters
from submitting.submitters         import *

import logging as lg


lg.basicConfig(level=lg.DEBUG)


CLUSTER = True  # if true, submit on a cluster using qsub


parameters = defaultParameters

parameters['time']              = 1200e3      # ms
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

parameters['output_dir']        = 'output'
parameters['nthreads']          = 8
parameters['ndumps']            = 1

startJobNum = 2000
numRepeat = 10

# Workstation parameters
programName         = 'python2.6 simulation_basic_grids.py'
blocking            = False

# Cluster parameters
cluster_scriptName  = 'cluster_submit.sh simulation_basic_grids.py'
qsub_params         = "-R y -P inf_ndtc -cwd -j y -l h_rt=01:30:00 -pe OpenMP 8"
qsub_output_dir     = parameters['output_dir']

ac = ArgumentCreator(parameters, printout=True)

iterparams = {
        'pc_max_rate'   : np.arange(0, 100, 10),
        'pc_conn_weight': np.arange(0.5, 5.5, 0.5)
}
ac.insertDict(iterparams, mult=True)

if CLUSTER:
    submitter = QsubSubmitter(ac, cluster_scriptName, qsub_params, qsub_output_dir)
else:
    submitter = GenericSubmitter(ac, programName, blocking=blocking)
submitter.submitAll(startJobNum, numRepeat, dry_run=False)

# Export the iterparams
iterparams_fname = "{0}/job{1:05}_iterparams".format(parameters['output_dir'], startJobNum)
submitter.exportIterParams(iterparams_fname)


