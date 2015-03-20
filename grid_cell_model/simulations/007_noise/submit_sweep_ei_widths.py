#!/usr/bin/env python
'''A parameter sweep that varies g_AMPA_total and pAMPA_sigma.

2D parameter sweep that simulates a stationary bump and records spiking
activity.
'''
from __future__ import absolute_import, print_function, division

import numpy as np
from grid_cell_model.submitting.noise import SubmissionParser
from grid_cell_model.submitting.factory import SubmitterFactory
from grid_cell_model.submitting.arguments import ArgumentCreator
from default_params import defaultParameters as dp

parser = SubmissionParser()
o = parser.parse_args()

for noise_sigma in parser.noise_sigmas:
    p = dp.copy()
    p['noise_sigma'] = noise_sigma # pA

    # Submitting
    ENV         = o.env
    simRootDir  = o.where
    simLabel    = '{0}pA'.format(int(p['noise_sigma']))
    appName     = '../common/simulation_stationary.py'
    rtLimit     = '01:30:00' if o.rtLimit is None else o.rtLimit
    numCPU      = 1
    blocking    = True
    timePrefix  = False
    numRepeat   = 1
    dry_run     = o.dry_run
    rc          = parser.rowcol

    p['master_seed'] = 123456
    p['time']        = 5e3 if o.time is None else o.time  # ms
    p['nthreads']    = 1
    p['ntrials']     = o.ntrials
    p['verbosity']   = o.verbosity

    p['gammaNSample'] = 0

    # Range of E/I synaptic conductances and pAMPA_sigma
    Nvals_G = 61      # Number of values for g_AMPA_total
    startG  = 0.0     # nS
    endG    = 3060.0  # nS

    Nvals_sigma = 11      # Number of values for pAMPA_sigma
    startSigma  = .5 / 6
    endSigma    = 1.

    ###############################################################################
    ac = ArgumentCreator(p, printout=o.printout)

    GArr = np.linspace(startG, endG, Nvals_G)
    sigmaArr = np.linspace(startSigma, endSigma, Nvals_sigma)
    print(GArr)
    print(sigmaArr)

    g_AMPA_total_arr = []
    pAMPA_sigma_arr = []
    for E_coupling in GArr:
        for pAMPA_sigma in sigmaArr:
            g_AMPA_total_arr.append(E_coupling)
            pAMPA_sigma_arr.append(pAMPA_sigma)

    iterparams = {
        'g_AMPA_total': np.array(g_AMPA_total_arr),
        'pAMPA_sigma': np.array(pAMPA_sigma_arr),
    }
    ac.insertDict(iterparams, mult=False)

    ###############################################################################
    submitter = SubmitterFactory.getSubmitter(
        ac, appName, envType=ENV, rtLimit=rtLimit, output_dir=simRootDir,
        label=simLabel, blocking=blocking, timePrefix=timePrefix, numCPU=numCPU)
    ac.setOption('output_dir', submitter.outputDir())
    startJobNum = 0
    rc_filter = rc[0]*len(sigmaArr) + rc[1] if rc is not None else None
    submitter.submitAll(startJobNum, numRepeat, dry_run=dry_run,
                        filter=rc_filter)
    submitter.saveIterParams(iterparams, dry_run=dry_run)
