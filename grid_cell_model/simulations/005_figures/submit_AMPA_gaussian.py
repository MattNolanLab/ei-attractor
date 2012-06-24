#
#   submit_AMPA_gaussian.py
#
#   Submit job(s) to the cluster/workstation: A gaussian excitatory weight profile
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

parameters['time']              = 10e3  # ms
#parameters['time']              = 3e3  # ms
parameters['ndumps']            = 1

parameters['prefDirC_e']        = 4
parameters['prefDirC_i']        = 0

parameters["Iext_e_const"]      = 400.0     # pA

parameters['AMPA_gaussian']     = 1         # bool
parameters["g_AMPA_total"]      = 4200      # nS
parameters["g_GABA_total"]      = 1200      # nS
parameters["g_uni_GABA_total"]  = 240       # nS

parameters['placeT']            = 10e3      # ms
parameters['placeDur']          = 100       # ms

parameters['bumpCurrentSlope']  = 1.108     # pA/(cm/s), !! this will depend on prefDirC !!
parameters['gridSep']           = 70        # cm, grid field inter-peak distance
parameters['theta_noise_sigma'] = 0         # pA


startJobNum =0
numRepeat = 4

# Workstation parameters
programName         = 'python2.6 simulation_AMPA_gaussian.py'
blocking            = False

# Cluster parameters
eddie_scriptName    = 'eddie_submit.sh simulation_AMPA_gaussian.py'
qsub_params         = "-P inf_ndtc -cwd -j y -l h_rt=13:00:00 -pe memory-2G 2"
qsub_output_dir     = parameters['output_dir']

ac = ArgumentCreator(parameters)

#iterparams = {
#    'Iplace'    :   [50, 100, 150, 200, 250]
#    'gridSep' : [50, 60]
#    'g_AMPA_total'  : [3600, 3800, 400, 4200]
#    'Iext_e_const'  : [375, 400, 425, 450]
#    'g_uni_GABA_total'  : [160, 200, 240, 280]
#}
#ac.insertDict(iterparams, mult=False)

if EDDIE:
    submitter = QsubSubmitter(ac, eddie_scriptName, qsub_params, qsub_output_dir)
else:
    submitter = GenericSubmitter(ac, programName, blocking=blocking)
submitter.submitAll(startJobNum, numRepeat, dry_run=False)
