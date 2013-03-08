#
#   submit_job_fig_model.py
#
#   Submit job(s) to the cluster/workstation: basic model figures
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

import logging  as lg
import numpy    as np


lg.basicConfig(level=lg.DEBUG)


CLUSTER = True  # if True, submit on a cluster using qsub


parameters = defaultParameters

parameters['time']              = 2e3   # ms
parameters['theta_start_mon_t'] = 1.0e3   # ms

parameters['prefDirC_e']        = 0
parameters['prefDirC_i']        = 0

surround_type = "I_surround"

##############################################################################
## E-surround (AMPA_ring_like params)
#
if surround_type == "E_surround":
    parameters['theta_noise_sigma'] = 0         # pA
##############################################################################



##############################################################################
# I_surround (AMPA_gaussian_params)
elif surround_type == "I_surround":
#    parameters['Iext_e_theta']      = 650       # pA
    parameters['Iext_i_theta']      = 25        # pA
    
    parameters['pAMPA_mu']          = 1.2/0.6
    
    parameters['AMPA_gaussian']     = 1         # bool
    parameters["g_AMPA_total"]      = 4500      # nS
    parameters["g_GABA_total"]      = 400       # nS
    parameters["g_uni_GABA_total"]  = 125        # nS
    
    parameters['placeT']            = 10e3      # ms
    parameters['placeDur']          = 100       # ms
    
    parameters['bumpCurrentSlope']  = 0.933     # pA/(cm/s), !! this will depend on prefDirC !!
    parameters['gridSep']           = 60        # cm, grid field inter-peak distance
    parameters['theta_noise_sigma'] = 0         # pA
##############################################################################
else:
    print "unknown surround profile type!"
    exit(1)


parameters['noise_sigma']       = 0.0       # mV
#parameters['g_uni_AMPA_total']  = 40         # nS



startJobNum = 0
numRepeat = 1

# Workstation parameters
programName         = 'python2.6 -i simulation_fig_model.py'
blocking            = False

# Cluster parameters
cluster_scriptName  = 'eddie_submit.sh simulation_fig_model.py'
qsub_params         = "-P inf_ndtc -cwd -j y -l h_rt=00:10:00"
qsub_output_dir     = parameters['output_dir']

ac = ArgumentCreator(parameters, printout=True)

iterparams = {
    'theta_noise_sigma' : [0, 50, 100, 150, 200, 250, 300, 350]
    #'noise_sigma'   : np.arange(0, 4.1, 0.2)
}
ac.insertDict(iterparams, mult=True)


if CLUSTER:
    submitter = QsubSubmitter(ac, cluster_scriptName, qsub_params, qsub_output_dir)
else:
    submitter = GenericSubmitter(ac, programName, blocking=blocking)
submitter.submitAll(startJobNum, numRepeat, dry_run=False)
