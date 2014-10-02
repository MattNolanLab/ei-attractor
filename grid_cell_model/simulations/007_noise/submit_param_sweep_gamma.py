#!/usr/bin/env python
'''2D parameter sweep that simulates a stationary bump and records spiking
activity and synaptic currents from selected neurons.'''
from __future__ import absolute_import, print_function, division

import logging as lg
#lg.basicConfig(level=lg.DEBUG)
lg.basicConfig(level=lg.INFO)

from default_params import defaultParameters as dp
from param_sweep    import submitParamSweep

noise_sigma_all = [.0] #, 150.0, 300.0] # pA

for noise_sigma in noise_sigma_all:
    p = dp.copy()
    p['noise_sigma'] = noise_sigma # pA

    # Submitting
    ENV         = 'cluster'
    simRootDir  = 'output/nmda_vm_dep/C_Mg_0p1_mM/gamma_bump'
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

    # Range of E/I synaptic conductances
    Nvals  = 31      # Number of values for each dimension
    startG = 0.0     # nS
    endG   = 6120.0  # nS

    ###############################################################################

    submitParamSweep(p, startG, endG, Nvals, ENV, simRootDir, simLabel,
            appName, rtLimit, numCPU, blocking, timePrefix, numRepeat, dry_run)
