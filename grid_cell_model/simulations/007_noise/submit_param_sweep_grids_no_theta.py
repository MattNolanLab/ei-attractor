#!/usr/bin/env python
'''
These are the same simulations as submit_param_sweep_grids, except that the
external theta current is replaced by a current with a constant amplitude that
delivers the same charge per theta cycle.
'''
import numpy as np

from grid_cell_model.submitting.noise import SubmissionParser
from param_sweep    import submitParamSweep, getBumpCurrentSlope
from default_params import defaultParameters as dp

parser = SubmissionParser()
parser.add_argument('--extra_qsub_params', type=str, default='')
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

    p['master_seed']      = 123456
    p['time']             = o.time or 600e3  # ms
    p['nthreads']         = 1
    p['ntrials']          = o.ntrials
    p['velON']            = 1
    p['constantPosition'] = 0
    p['verbosity']        = o.verbosity

    # No theta parameters
    p['Iext_e_const'] = 487.5   # pA
    p['Iext_i_const'] = 212.5   # pA
    p['Iext_e_theta'] = 0       # pA
    p['Iext_i_theta'] = 0       # pA

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
                     appName, rtLimit, numCPU, blocking, timePrefix, numRepeat,
                     dry_run, extraIterparams, rc=parser.rowcol,
                     extra_qsub_params=o.extra_qsub_params)
