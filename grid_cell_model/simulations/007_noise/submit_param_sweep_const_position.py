#!/usr/bin/env python
'''
E/I parameter sweeps in which velocity input is turned OFF, place cell input is
on and pointing to a constant position (i.e. the animat does not move).

Initialisation place cell inputs and "velocity" place cell input positions
should be in sync. In this particular simulation, the position is at [0, 0].
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
    p['verbosity'] = 'DEBUG'
    p['noise_sigma'] = noise_sigma # pA

    # Submitting
    ENV         = 'cluster'
    simRootDir  = 'output/even_spacing/const_position'
    simLabel    = '{0}pA'.format(int(p['noise_sigma']))
    appName     = 'simulation_grids.py'
    rtLimit     = '01:00:00'
    numCPU      = 1
    blocking    = True
    timePrefix  = False
    numRepeat   = 1
    dry_run     = False

    p['time']              = 10.5e3  # ms
    p['nthreads']          = 1
    p['ntrials']           = 10
    p['bumpCurrentSlope']  = -1
    p['velON']             = 0
    p['constantPosition']  = 1
    p['stateMonDur']       = 1e3 # ms



    # Range of E/I synaptic conductances
    Nvals  = 31      # Number of values for each dimension
    startG = 0.0     # nS
    endG   = 6120.0  # nS


    ###############################################################################

    submitParamSweep(p, startG, endG, Nvals, ENV, simRootDir, simLabel,
            appName, rtLimit, numCPU, blocking, timePrefix, numRepeat, dry_run)

