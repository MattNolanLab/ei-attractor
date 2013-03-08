#
#   submit_AMPA_gaussian_velocity.py
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
from common         import *

import logging as lg


lg.basicConfig(level=lg.DEBUG)


CLUSTER = True  # if True, submit on a cluster using qsub


parameters = defaultParameters

parameters['time']              = 10e3      # ms
parameters['ngenerations']      = 10
parameters['velModulationType'] = 'excitatory'
parameters['prefDirC_e']        = 4
parameters['prefDirC_i']        = 0

parameters['Iext_e_theta']      = 650       # pA
parameters['Iext_i_theta']      = 50        # pA

parameters['pAMPA_mu']          = 1.2/0.6

parameters['AMPA_gaussian']     = 1         # bool
parameters["g_AMPA_total"]      = 4500      # nS
parameters["g_GABA_total"]      = 400       # nS
parameters["g_uni_GABA_total"]  = 125        # nS

parameters['theta_noise_sigma'] = 0         # pA

#parameters['Ivel']              = 120        # pA

startJobNum = 4300
numRepeat = 1

# Workstation parameters
programName         = 'nice python2.6 simulation_AMPA_gaussian_velocity.py'
blocking            = False

# Cluster parameters
cluster_scriptName  = 'eddie_submit.sh simulation_AMPA_gaussian_velocity.py'
qsub_params         = "-P inf_ndtc -cwd -j y -l h_rt=02:00:00 -pe memory-2G 2"
qsub_output_dir     = parameters['output_dir']

ac = ArgumentCreator(parameters)

iterparams = {
        'Ivel'       : np.arange(0, 113, 10)}
ac.insertDict(iterparams, mult=True)

if CLUSTER:
    submitter = QsubSubmitter(ac, cluster_scriptName, qsub_params, qsub_output_dir)
else:
    submitter = GenericSubmitter(ac, programName, blocking=blocking)
submitter.submitAll(startJobNum, numRepeat, dry_run=False)
