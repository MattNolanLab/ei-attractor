#!/usr/bin/env python
'''Submit job(s) to the cluster/workstation: grid field parameter sweeps.'''

import numpy as np

from param_sweep    import (submitParamSweep, getBumpCurrentSlope,
                            SubmissionParser)
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
    appName     = 'simulation_grids.py'
    rtLimit     = o.rtLimit
    numCPU      = 1
    blocking    = True
    timePrefix  = False
    numRepeat   = 1
    dry_run     = o.dry_run

    p['master_seed']      = 123456
    p['time']             = 600e3 if o.time is None else o.time  # ms
    p['nthreads']         = 1
    p['ntrials']          = o.ntrials
    p['velON']            = 1
    p['pcON']             = 1
    p['constantPosition'] = 0
    p['verbosity']        = o.verbosity


    # Range of E/I synaptic conductances
    Nvals  = 31      # Number of values for each dimension
    startG = 0.0     # nS
    endG   = 6120.0  # nS

    extraIterparams = {'bumpCurrentSlope' : getBumpCurrentSlope(p['noise_sigma'],
        threshold=-np.infty)}
    #extraIterparams['bumpCurrentSlope'] = [1.0]

    ###############################################################################
    submitParamSweep(p, startG, endG, Nvals, ENV, simRootDir, simLabel,
                     appName, rtLimit, numCPU, blocking, timePrefix, numRepeat,
                     dry_run, extraIterparams, rc=parser.rowcol,
                     printout=o.printout)
