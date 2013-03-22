#! /usr/bin/env python
#
#   submit_bump_fitting.py
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


CLUSTER = False  # if True, submit on a cluster using qsub


parameters = defaultParameters

#parameters['time']              = 1199.9e3  # ms
parameters['time']              = 10e3  # ms
parameters['nthreads']          = 8

parameters['output_dir']        = 'output_local'

parameters['Ne']                = 34
parameters['Ni']                = 34

parameters['bumpCurrentSlope']  = 1.175     # pA/(cm/s), !! this will depend on prefDirC !!
parameters['gridSep']           = 70.0      # cm, grid field inter-peak distance
parameters['N_place_cells']     = 30*30

# Gamma analysis parameters
parameters['gammaNSample']      = 0.05      # fraction
parameters['gammaRangeLow']     = 50.0      # Hz
parameters['gammaRangeHigh']    = 100.0     # Hz

parameters['noise_sigma']       = 150.0       # pA
parameters['delay']             = 0.1


startJobNum = 0
numRepeat = 1

# Workstation parameters
programName         = 'python simulation_bump_fitting.py'
blocking            = True

# Cluster parameters
cluster_scriptName  = 'cluster_submit.sh simulation_bump_fitting.py'
qsub_params         = "-P inf_ndtc -cwd -j y -l h_rt=00:05:00 -pe openmpi_smp8_mark2 16"
qsub_output_dir     = parameters['output_dir']

ac = ArgumentCreator(parameters)

iterparams = {
        'noise_sigma' : [0.0, 150.0, 200, 250, 300]    # pA
        #'noise_sigma' : [150.0]
}
ac.insertDict(iterparams, mult=False)

if CLUSTER:
    submitter = QsubSubmitter(ac, cluster_scriptName, qsub_params, qsub_output_dir)
else:
    submitter = GenericSubmitter(ac, programName, blocking=blocking)
submitter.submitAll(startJobNum, numRepeat, dry_run=False)
