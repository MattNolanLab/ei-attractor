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


EDDIE = True  # if eddie, submit on a cluster using qsub


parameters = defaultParameters

parameters['time']              = 2e3  # ms
#parameters['time']              = 3e3  # ms
parameters['ndumps']            = 1

parameters['prefDirC_e']        = 4
parameters['prefDirC_i']        = 0

parameters['Iext_e_theta']      = 650       # pA
parameters['Iext_i_theta']      = 50        # pA

parameters['pAMPA_mu']          = 0.5/0.6

parameters['AMPA_gaussian']     = 1         # bool
parameters["g_AMPA_total"]      = 4500      # nS
parameters["g_GABA_total"]      = 400       # nS
parameters["g_uni_GABA_total"]  = 125        # nS

parameters['theta_noise_sigma'] = 0         # pA

startJobNum =9000
numRepeat = 1

# Workstation parameters
programName         = 'python2.6 simulation_pharma.py'
blocking            = False

# Cluster parameters
eddie_scriptName    = 'eddie_submit.sh simulation_pharma.py'
qsub_params         = "-P inf_ndtc -cwd -j y -l h_rt=00:05:00 -pe memory-2G 1"
qsub_output_dir     = parameters['output_dir']

ac = ArgumentCreator(parameters)

iterparams = {
#    'Iplace'    :   [50, 100, 150, 200, 250]
#    'gridSep' : [50, 60]
#    'sigmaIextGaussian' : np.array([0.25, 0.5, 0.75, 1.0])/6.0
    'Iext_e_theta'  : [150, 200, 250, 300, 350, 400, 450, 500, 550, 600]
    'Iext_i_theta'  : [50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190]
}
ac.insertDict(iterparams, mult=True)

if EDDIE:
    submitter = QsubSubmitter(ac, eddie_scriptName, qsub_params, qsub_output_dir)
else:
    submitter = GenericSubmitter(ac, programName, blocking=blocking)
submitter.submitAll(startJobNum, numRepeat, dry_run=False)
