#!/usr/bin/env python
'''
Grid field parameter sweeps: velocity input turned OFF, but place cell input
left ON.
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
    rtLimit     = o.rtLimit or '05:00:00'
    numCPU      = 1
    blocking    = True
    timePrefix  = False
    numRepeat   = 1
    dry_run     = o.dry_run

    p['master_seed']       = 123456
    p['time']              = o.time or 600e3  # ms
    p['nthreads']          = 1
    p['ntrials']           = o.ntrials
    p['bumpCurrentSlope']  = -1
    p['velON']             = 0
    p['constantPosition']  = 0
    p['pcON']              = 1
    p['verbosity']         = o.verbosity

    # Range of E/I synaptic conductances
    Nvals  = 31      # Number of values for each dimension
    startG = 0.0     # nS
    endG   = 6120.0  # nS

    ###############################################################################
    submitParamSweep(p, startG, endG, Nvals, ENV, simRootDir, simLabel,
                     appName, rtLimit, numCPU, blocking, timePrefix, numRepeat,
                     dry_run, rc=parser.rowcol, printout=o.printout)
