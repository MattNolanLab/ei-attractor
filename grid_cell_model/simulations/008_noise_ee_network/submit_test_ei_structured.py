#!/usr/bin/env python
'''One parameter sweep in networks with E->E and E<->I connections structured

Run a network for a short time and save data. These simulations are meant to be
run for interactive tuning of the network.
'''
from __future__ import absolute_import, print_function, division

import numpy as np
from grid_cell_model.submitting.factory import SubmitterFactory
from grid_cell_model.submitting.arguments import ArgumentCreator
from grid_cell_model.submitting.noise import SingleParameterSweepParser
from grid_cell_model.submitting.flagparse import positive_int
from default_params import defaultParameters as dp

parser = SingleParameterSweepParser()
parser.add_argument('--nCPU',    type=positive_int, default=1)
o = parser.parse_args()

for noise_sigma in parser.noise_sigmas:
    p = dp.copy()
    p['noise_sigma'] = noise_sigma  # pA

    # Submitting
    ENV         = o.env
    simRootDir  = o.where
    simLabel    = '{0}pA'.format(int(p['noise_sigma']))
    appName     = '../common/simulation_stationary.py'
    rtLimit     = o.rtLimit
    numCPU      = o.nCPU
    blocking    = True
    timePrefix  = False
    numRepeat   = 1
    dry_run     = o.dry_run

    p['master_seed'] = 123456
    p['time']        = 10e3 if o.time is None else o.time  # ms
    p['nthreads']    = 1
    p['ntrials']     = o.ntrials
    p['verbosity']   = o.verbosity

    p['EI_flat'] = 0
    p['IE_flat'] = 0
    p['use_EE']  = 1

    p['g_AMPA_total'] = 2040.    # nS
    p['g_GABA_total'] = 2040.    # nS
    p['g_EE_total']   = 100.     # nS

    # Here, no PC inputs
    #p['pc_start_max_rate'] = .0  # Hz

    p['prefDirC_e'] = 0.

    iterparams = {
        o.explored_param: np.arange(o.param_start, o.param_stop + o.param_step,
                                    o.param_step)
    }
    ac = ArgumentCreator(p, printout=True)
    ac.insertDict(iterparams, mult=False)

    ###########################################################################
    submitter = SubmitterFactory.getSubmitter(
        ac, appName, envType=ENV, rtLimit=rtLimit, output_dir=simRootDir,
        label=simLabel, blocking=blocking, timePrefix=timePrefix,
        numCPU=numCPU)
    ac.setOption('output_dir', submitter.outputDir())
    startJobNum = 0
    submitter.submitAll(startJobNum, numRepeat, dry_run=dry_run,
                        filter=parser.filter)
    submitter.saveIterParams(iterparams, dry_run=dry_run)
