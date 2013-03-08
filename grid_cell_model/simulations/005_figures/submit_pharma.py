#
#   submit_pharma.py
#
#   Submit job(s) to the cluster/workstation.
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

parameters['time']              = 2e3  # ms
#parameters['time']              = 3e3  # ms
parameters['ndumps']            = 1

parameters['prefDirC_e']        = 4
parameters['prefDirC_i']        = 0

parameters['Iext_e_theta']      = 400       # pA
parameters['Iext_i_theta']      = 850       # pA

parameters['pAMPA_sigma']       = 0.5/6.0

parameters['AMPA_gaussian']     = 1         # bool
parameters["g_AMPA_total"]      = 2250      # nS
#parameters["g_AMPA_total"]      = 0        # nS
parameters["g_GABA_total"]      = 320       # nS
parameters["g_uni_GABA_total"]  = 300       # nS

parameters['sigmaIextGaussian'] = 0.6/6.0
parameters['shuffleIextGaussian'] = 0

parameters['condAddPercSynapses_e'] = 50
parameters['condAdd_e']             = 20.0  # nS

parameters['tau_AMPA']          = 2 # ms

parameters['theta_noise_sigma'] = 0         # pA

startJobNum = 8000
numRepeat = 1

# Workstation parameters
programName         = 'python2.6 simulation_pharma.py'
blocking            = False

# Cluster parameters
cluster_scriptName  = 'eddie_submit.sh simulation_pharma.py'
qsub_params         = "-P inf_ndtc -cwd -j y -l h_rt=00:03:00 -pe memory-2G 1"
qsub_output_dir     = parameters['output_dir']

ac = ArgumentCreator(parameters, printout=True)

iterparams = {
#    'Iplace'    :   [50, 100, 150, 200, 250]
#    'gridSep' : [50, 60]
    'sigmaIextGaussian' : np.arange(0.6, 1.6, 0.1)/6.0,
#    'g_uni_GABA_total'      : [75, 100, 125, 150, 175, 200, 225, 250, 275, 300],
    'Iext_e_theta'  : [750, 775, 800, 825, 850, 875, 900, 925, 950, 1000],
#    'Iext_i_theta'  : [650, 700, 750, 800, 850]
#    'condAddPercSynapses_e' : [ 5, 10, 15, 20, 25, 30, 35, 40, 45, 50],
#    'condAdd_e'             : [20, 25, 30, 35, 40]
}
ac.insertDict(iterparams, mult=True)

if CLUSTER:
    submitter = QsubSubmitter(ac, cluster_scriptName, qsub_params, qsub_output_dir)
else:
    submitter = GenericSubmitter(ac, programName, blocking=blocking)
submitter.submitAll(startJobNum, numRepeat, dry_run=False)
