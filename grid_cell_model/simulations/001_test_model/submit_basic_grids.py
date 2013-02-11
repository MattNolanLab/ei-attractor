#
#   submit_job.py
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
from common         import *


EDDIE = False  # if eddie, submit on a cluster using qsub


parameters = defaultParameters

#parameters['time']              = 1199.9e3  # ms
parameters['time']              = 5.0e3  # ms

parameters['output_dir']        = 'output'



startJobNum = 0
numRepeat = 1

# Workstation parameters
programName         = 'python2.6 -i simulation_basic_grids.py'
#programName         = 'python2.6 -i -m cProfile -s time -o profile.txt simulation_basic_grids.py'
blocking            = True

# Cluster parameters
eddie_scriptName    = 'eddie_submit.sh simulation_basic_grids.py'
qsub_params         = "-P inf_ndtc -cwd -j y -l h_rt=13:00:00 -pe memory-2G 2"
qsub_output_dir     = parameters['output_dir']

ac = ArgumentCreator(parameters)

#iterparams = {
#        'bumpCurrentSlope'  : [1.15, 1.175, 1.2]
#    'g_AMPA_total' : [defaultParameters['g_AMPA_total'], 0.0]
#    'tau_AHP_e'  :   [20,   30,  40,  50,  60,  70,  80,  90],
#    'Iext_e_theta' : [375, 400, 450, 500, 550, 600, 650, 700]
#}
#ac.insertDict(iterparams, mult=False)

if EDDIE:
    submitter = QsubSubmitter(ac, eddie_scriptName, qsub_params, qsub_output_dir)
else:
    submitter = GenericSubmitter(ac, programName, blocking=blocking)
submitter.submitAll(startJobNum, numRepeat, dry_run=False)
