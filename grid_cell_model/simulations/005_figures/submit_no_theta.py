#
#   submit_no_theta.py
#
#   Submit job(s) to the cluster/workstation: no theta input
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


CLUSTER = True  # if True, submit on a cluster using qsub


parameters = defaultParameters

parameters['time']              = 1199.9e3       # ms
parameters['ndumps']            = 20

parameters['Iext_e_const']      = 500       # pA
parameters['Iext_i_const']      = 200       # pA

parameters['placeT']            = 10e3      # ms

parameters['prefDirC_e']        = 4
parameters['prefDirC_i']        = 0
parameters['bumpCurrentSlope']  = 1.45     # pA/(cm/s), !! this will depend on prefDirC !!
parameters['gridSep']           = 70        # cm, grid field inter-peak distance
parameters['theta_noise_sigma'] = 0         # pA
startJobNum = 3300
numRepeat = 10

# Workstation parameters
programName         = 'python2.6 simulation_no_theta.py'
blocking            = False

# Cluster parameters
cluster_scriptName  = 'eddie_submit.sh simulation_no_theta.py'
qsub_params         = "-P inf_ndtc -cwd -j y -l h_rt=14:00:00 -pe memory-2G 2"
qsub_output_dir     = parameters['output_dir']

ac = ArgumentCreator(parameters)

#iterparams = {
#        'bumpCurrentSlope'    : [0.821, 0.921, 1.021]}
#ac.insertDict(iterparams, mult=False)

if CLUSTER:
    submitter = QsubSubmitter(ac, cluster_scriptName, qsub_params, qsub_output_dir)
else:
    submitter = GenericSubmitter(ac, programName, blocking=blocking)
submitter.submitAll(startJobNum, numRepeat, dry_run=False)
