#!/usr/bin/env python
'''
These are the same simulations as submit_param_sweep_grids, except that the
external theta current is replaced by a current with a constant amplitude that
delivers the same charge per theta cycle.
'''
import numpy as np
from submitting.factory   import SubmitterFactory
from submitting.arguments import ArgumentCreator
from default_params       import defaultParameters as dp
from param_sweep          import submitParamSweep, getBumpCurrentSlope
import logging as lg
#lg.basicConfig(level=lg.DEBUG)
lg.basicConfig(level=lg.INFO)

noise_sigma_all = [0.0, 150.0, 300.0] # pA

for noise_sigma in noise_sigma_all:
    p = dp.copy()
    p['noise_sigma'] = noise_sigma # pA

    # Submitting
    ENV         = 'cluster'
    simRootDir  = 'output/no_theta/grids'
    simLabel    = '{0}pA'.format(int(p['noise_sigma']))
    appName     = 'simulation_grids.py'
    rtLimit     = '05:00:00'
    numCPU      = 1
    blocking    = True
    timePrefix  = False
    numRepeat   = 1
    dry_run     = False

    p['time']              = 600e3  # ms
    p['nthreads']          = 1
    p['ntrials']           = 1
    p['velON']             = 1


    # Range of E/I synaptic conductances
    Nvals  = 31      # Number of values for each dimension
    startG = 0.0     # nS
    endG   = 6120.0  # nS

    bumpCurrentSlope = getBumpCurrentSlope(p['noise_sigma'],
            threshold=-np.infty, type='no_theta')
    extraIterparams = {'bumpCurrentSlope' : bumpCurrentSlope}
    #extraIterparams['bumpCurrentSlope'] = [1.0]

    ###############################################################################

    submitParamSweep(p, startG, endG, Nvals, ENV, simRootDir, simLabel,
            appName, rtLimit, numCPU, blocking, timePrefix, numRepeat, dry_run,
            extraIterparams)

