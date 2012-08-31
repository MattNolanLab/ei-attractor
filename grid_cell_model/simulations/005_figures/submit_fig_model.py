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


EDDIE = False  # if eddie, submit on a cluster using qsub


parameters = defaultParameters

parameters['time']              = 6.0e3   # ms
parameters['theta_start_mon_t'] = 1.0e3   # ms

parameters['prefDirC_e']        = 4
parameters['prefDirC_i']        = 0

parameters['placeT']            = 10e3      # ms
parameters['placeDur']          = 100       # ms

parameters['Iext_e_theta']      = 2000.0    # pA
parameters['Iext_i_theta']      = 800.0     # pA

parameters['Iext_i_const']      = 200.0     # pA

parameters['taum_i_spread']     = 2.0       # ms
parameters['EL_i_spread']       = 10.0      # mV

parameters['g_AMPA_total']      = 0         # nS
parameters['g_GABA_total']      = 550.0     # nS
parameters['g_uni_GABA_total']  = 36.0      # nS



startJobNum = 0
numRepeat = 4

# Workstation parameters
programName         = 'python2.6 -i simulation_fig_model.py'
blocking            = False

# Cluster parameters
eddie_scriptName    = 'eddie_submit.sh simulation_fig_model.py'
qsub_params         = "-P inf_ndtc -cwd -j y -l h_rt=13:00:00 -pe memory-2G 2"
qsub_output_dir     = parameters['output_dir']

ac = ArgumentCreator(parameters)

#iterparams = {
##        'Iext_e_theta' : [500, 525, 550, 575, 600, 625, 650, 700]
##        'taum_e_spread' : [0.5, 0.75,   1, 1.25, 1.5, 1.75,   2, 2.25],
##        'EL_e_spread'   : [0.5,    1, 1.5,    2, 2.5,    3, 3.5,    4]
##         'theta_noise_sigma' : [0, 20, 40, 60, 80, 100, 120, 140]
##        'g_GABA_total'  : [500, 525, 550, 575]
##        'g_uni_GABA_total'  : [0, 25, 50, 75, 100, 125, 150, 175]
##    'pAMPA_mu'  : np.array([0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]) / 0.6
#}
#ac.insertDict(iterparams, mult=False)


if EDDIE:
    submitter = QsubSubmitter(ac, eddie_scriptName, qsub_params, qsub_output_dir)
else:
    submitter = GenericSubmitter(ac, programName, blocking=blocking)
submitter.submitAll(startJobNum, numRepeat, dry_run=False)
