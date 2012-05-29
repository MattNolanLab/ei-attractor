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

import logging as lg


lg.basicConfig(level=lg.DEBUG)


EDDIE = False  # if eddie, submit on a cluster using qsub


parameters = defaultParameters

parameters['time']              = 6e3   # ms
parameters['theta_start_mon_t'] = 1e3   # ms

#parameters['Iext_e_theta']      =    # pA
parameters['g_GABA_total']      = 0     # nS

#parameters['Iext_e_const']      = 500       # pA
#parameters['Iext_i_const']      = 200       # pA

parameters['noise_sigma']       = 2     # mV

startJobNum = 50
numRepeat = 1

# Workstation parameters
programName         = 'python2.6 -i simulation_fig_model.py'
blocking            = False

# Cluster parameters
eddie_scriptName    = 'eddie_submit.sh simulation_fig_model.py'
qsub_params         = "-P inf_ndtc -cwd -j y -l h_rt=13:00:00 -pe memory-2G 2"
qsub_output_dir     = parameters['output_dir']

ac = ArgumentCreator(parameters)

iterparams = {
        'g_AMPA_total'  : [150, 175, 200, 225, 250, 275, 300, 325]}
ac.insertDict(iterparams, mult=False)


if EDDIE:
    submitter = QsubSubmitter(ac, eddie_scriptName, qsub_params, qsub_output_dir)
else:
    submitter = GenericSubmitter(ac, programName, blocking=blocking)
submitter.submitAll(startJobNum, numRepeat, dry_run=False)
