#
#   submit_basic_grids_velocity.py
#
#   Submit job(s) to the cluster/workstation: velocity estimation theta+gamma
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

parameters['time']              = 10e3      # ms
parameters['delay']             = 0.1       # ms

parameters['ngenerations']      = 10
parameters['velModulationType'] = 'excitatory'
parameters['prefDirC_e']        = 4
parameters['prefDirC_i']        = 0

parameters['theta_noise_sigma'] = 0          # pA
parameters['noise_sigma']       = 150.       # pA

parameters['Ne']                = 34
parameters['Ni']                = 34
parameters['N_place_cells']     = 30*30
parameters['gridSep']           = 70.0      # cm, grid field inter-peak distance

#parameters['Ivel']              = 40        # pA

parameters['output_dir']        = 'output'
parameters['nthreads']          = 8

startJobNum = 1000
numRepeat = 1

# Workstation parameters
programName         = 'python2.6 simulation_basic_grids_velocity.py'
blocking            = False

# Cluster parameters
cluster_scriptName  = 'cluster_submit.sh simulation_basic_grids_velocity.py'
qsub_params         = "-P inf_ndtc -cwd -j y -l h_rt=00:20:00 -pe OpenMP 8"
qsub_output_dir     = parameters['output_dir']

ac = ArgumentCreator(parameters, printout=True)

iterparams = {
#        'Iext_e_theta' : np.arange(0, 360, 30),
        'Ivel'         : np.arange(0, 160, 10)}
ac.insertDict(iterparams, mult=True)

if CLUSTER:
    submitter = QsubSubmitter(ac, cluster_scriptName, qsub_params, qsub_output_dir)
else:
    submitter = GenericSubmitter(ac, programName, blocking=blocking)
submitter.submitAll(startJobNum, numRepeat, dry_run=False)

