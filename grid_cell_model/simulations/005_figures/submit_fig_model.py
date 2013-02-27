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

parameters['time']              = 6e3   # ms
parameters['theta_start_mon_t'] = 1.0e3   # ms

parameters['output_dir'] = 'output_local'

parameters['prefDirC_e']        = 0
parameters['prefDirC_i']        = 0

parameters['Ne'] = 34
parameters['Ni'] = 34

surround_type = "E_surround"

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


#parameters['noise_sigma']       = 0.0       # mV
#parameters['g_uni_AMPA_total']  = 40         # nS



startJobNum = 0
numRepeat = 1

# Workstation parameters
programName         = 'python2.6 -i simulation_fig_model.py'
blocking            = False

# Cluster parameters
eddie_scriptName    = 'eddie_submit.sh simulation_fig_model.py'
qsub_params         = "-P inf_ndtc -cwd -j y -l h_rt=00:10:00"
qsub_output_dir     = parameters['output_dir']

ac = ArgumentCreator(parameters)

#iterparams = {
#        'Ne'    : [8, 12, 16, 20, 24, 28, 32],
#        'Ni'    : [8, 12, 16, 20, 24, 28, 32]
#        
#        'g_uni_AMPA_total'  : np.arange(40, 320, 20) 
#        'Iext_e_theta' : np.arange(375, 625, 25),
#        'Iext_i_theta' : [25, 50]
#        'taum_e_spread' : [0.5, 0.75,   1, 1.25, 1.5, 1.75,   2, 2.25],
#        'EL_e_spread'   : [0.5,    1, 1.5,    2, 2.5,    3, 3.5,    4]
#         'theta_noise_sigma' : [0, 20, 40, 60, 80, 100, 120, 140]
#        'g_GABA_total'  : [500, 525, 550, 575]
#        'g_uni_GABA_total'  : [0, 25, 50, 75, 100, 125, 150, 175]
#    'pAMPA_mu'  : np.array([0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]) / 0.6
#        'noise_sigma'   : np.arange(0, 4.1, 0.2)
#}
#ac.insertDict(iterparams, mult=False, printout=True)


if EDDIE:
    submitter = QsubSubmitter(ac, eddie_scriptName, qsub_params, qsub_output_dir)
else:
    submitter = GenericSubmitter(ac, programName, blocking=blocking)
submitter.submitAll(startJobNum, numRepeat, dry_run=False)
