#
#   submit_basic_grids.py
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

import logging  as lg
import numpy    as np


lg.basicConfig(level=lg.DEBUG)


EDDIE = True  # if eddie, submit on a cluster using qsub


parameters = defaultParameters

parameters['time']              = 2e3  # ms
#parameters['time']              = 3e3  # ms
parameters['ndumps']            = 1

parameters['prefDirC_e']        = 4
parameters['prefDirC_i']        = 0

parameters['placeT']            = 10e3      # ms
parameters['placeDur']          = 100       # ms

parameters['bumpCurrentSlope']  = 1.05      # pA/(cm/s), !! this will depend on prefDirC !!
#parameters['gridSep']           = 40        # cm, grid field inter-peak distance
parameters['theta_noise_sigma'] = 0         # pA


parameters['Iext_e_theta']      = 2000.0    # pA
parameters['Iext_i_theta']      = 800.0     # pA

parameters['Iext_i_const']      = 200.0     # pA

parameters['taum_i_spread']     = 2.0       # ms
parameters['EL_i_spread']       = 10.0      # mV

parameters['g_AMPA_total']      = 0         # nS
parameters['g_GABA_total']      = 550.0     # nS
parameters['g_uni_GABA_total']  = 36.0      # nS


startJobNum =9000
numRepeat = 1

# Workstation parameters
programName         = 'python2.6 simulation_basic_grids.py'
blocking            = False

# Cluster parameters
eddie_scriptName    = 'eddie_submit.sh simulation_basic_grids.py'
qsub_params         = "-P inf_ndtc -cwd -j y -l h_rt=00:04:00 -pe memory-2G 1"
qsub_output_dir     = parameters['output_dir']

ac = ArgumentCreator(parameters)

iterparams = {
#    'Iplace'    :   [50, 100, 150, 200, 250]
#    'gridSep' : [50, 60]
#    'Iext_e_theta'      : np.arange(1400, 2000, 50)
    'g_uni_GABA_total'  : np.arange(20, 40, 1)

#    'g_GABA_total'      : np.arange(100, 1100, 50),
#    'taum_i_spread' : np.arange(2.0, 4.0, 0.2),
#    'EL_i_spread'   : np.arange(10, 20,  1)
#    'Iext_i_theta' : [750, 800, 850, 900],
#    'Iext_i_const'  : np.arange(200, 800, 50)
}
ac.insertDict(iterparams, mult=True, printout=True)

if EDDIE:
    submitter = QsubSubmitter(ac, eddie_scriptName, qsub_params, qsub_output_dir)
else:
    submitter = GenericSubmitter(ac, programName, blocking=blocking)
submitter.submitAll(startJobNum, numRepeat, dry_run=False)
