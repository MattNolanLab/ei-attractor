#!/usr/bin/env python
'''Submit a simulation to generate raster plot data for the thesis.'''

from default_params import defaultParameters
from common         import *

import logging  as lg
import numpy    as np


lg.basicConfig(level=lg.DEBUG)


EDDIE = False  # if eddie, submit on a cluster using qsub


parameters = defaultParameters

parameters['time']              = 6e3   # ms
parameters['theta_start_mon_t'] = 1.0e3   # ms

parameters['prefDirC_e']        = 0
parameters['prefDirC_i']        = 0
parameters['theta_noise_sigma'] = 0         # pA

##############################################################################

startJobNum = 50
numRepeat = 1

# Workstation parameters
programName         = 'python -i simulation_thesis_rasters.py'
blocking            = False

# Cluster parameters
eddie_scriptName    = 'eddie_submit.sh simulation_fig_model.py'
qsub_params         = "-P inf_ndtc -cwd -j y -l h_rt=13:00:00 -pe memory-2G 2"
qsub_output_dir     = parameters['output_dir']

ac = ArgumentCreator(parameters)

iterparams = {'noise_sigma'   : [2]}
ac.insertDict(iterparams, mult=False)


if EDDIE:
    submitter = QsubSubmitter(ac, eddie_scriptName, qsub_params, qsub_output_dir)
else:
    submitter = GenericSubmitter(ac, programName, blocking=blocking)
submitter.submitAll(startJobNum, numRepeat, dry_run=False)
