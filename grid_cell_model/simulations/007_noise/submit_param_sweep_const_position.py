#!/usr/bin/env python
'''
E/I parameter sweeps in which velocity input is turned OFF, place cell input is
on and pointing to a constant position (i.e. the animal does not move).

Initialisation place cell inputs and "velocity" place cell input positions
should be in sync. In this particular simulation, the position is at [0, 0].
'''
import numpy as np

from grid_cell_model.submitting.noise import SubmissionParser
from param_sweep    import submitParamSweep, getBumpCurrentSlope
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

    p['master_seed']       = 123456
    p['time']              = o.time or 10.5e3  # ms
    p['nthreads']          = 1
    p['ntrials']           = o.ntrials
    p['bumpCurrentSlope']  = -1
    p['velON']             = 0
    p['constantPosition']  = 1
    p['pcON']              = 1
    p['stateMonDur']       = 1e3 # ms
    p['verbosity']         = o.verbosity

    # Range of E/I synaptic conductances
    Nvals  = 31      # Number of values for each dimension
    startG = 0.0     # nS
    endG   = 6120.0  # nS

    ###############################################################################
    submitParamSweep(p, startG, endG, Nvals, ENV, simRootDir, simLabel,
                     appName, rtLimit, numCPU, blocking, timePrefix, numRepeat,
                     dry_run, rc=parser.rowcol, printout=o.printout)
