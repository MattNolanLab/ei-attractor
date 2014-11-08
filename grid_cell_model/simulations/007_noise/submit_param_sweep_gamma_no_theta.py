#!/usr/bin/env python
'''
These are the same simulations as submit_param_sweep_gamma, except that the
external theta current is replaced by a current with a constant amplitude that
delivers the same charge per theta cycle.
'''
from default_params import defaultParameters as dp
from param_sweep    import submitParamSweep
import logging as lg
#lg.basicConfig(level=lg.DEBUG)
lg.basicConfig(level=lg.INFO)

noise_sigma_all = [0.0, 150.0, 300.0] # pA

for noise_sigma in noise_sigma_all:
    p = dp.copy()
    p['noise_sigma'] = noise_sigma # pA

    # Submitting
    ENV         = 'cluster'
    simRootDir  = 'output/submission/no_theta/gamma_bump'
    simLabel    = '{0}pA'.format(int(p['noise_sigma']))
    appName     = '../common/simulation_stationary.py'
    rtLimit     = '01:30:00'
    numCPU      = 1
    blocking    = True
    timePrefix  = False
    numRepeat   = 1
    dry_run     = False

    p['time']              = 10e3  # ms
    p['nthreads']          = 1
    p['ntrials']           = 5

    # No theta parameters
    p['Iext_e_const'] = 487.5   # pA
    p['Iext_i_const'] = 212.5   # pA
    p['Iext_e_theta'] = 0       # pA
    p['Iext_i_theta'] = 0       # pA


    # Range of E/I synaptic conductances
    Nvals  = 31      # Number of values for each dimension
    startG = 0.0     # nS
    endG   = 6120.0  # nS

    ###############################################################################

    submitParamSweep(p, startG, endG, Nvals, ENV, simRootDir, simLabel,
            appName, rtLimit, numCPU, blocking, timePrefix, numRepeat, dry_run)
