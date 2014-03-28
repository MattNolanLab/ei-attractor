#!/usr/bin/env python
'''
Bump velocity gain estimation.

These simulations do the following:
    - Bump is initialized, in the beginning [0, theta_start_t], by a very
      strong place cell input.
    - The bump moves, with a constant speed, following **vertical** direction.
      for the duration of the simulation.
    - The speed is controlled by IvelMax and dIvel parameters (currently
      [0, 100] pA, with a step of 10 pA
    - At the end of each run, spikes from E and I population are exported to
      the output file.
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
    simRootDir  = 'output/even_spacing/velocity_vertical'
    simLabel    = '{0}pA'.format(int(p['noise_sigma']))
    appName     = 'simulation_velocity.py'
    rtLimit     = '08:00:00'
    numCPU      = 1
    blocking    = True
    timePrefix  = False
    numRepeat   = 1
    dry_run     = False


    p['time']     = 10e3  # ms
    p['nthreads'] = 1
    p['ntrials']  = 10

    p['IvelMax']  = 200
    p['dIvel']    = 10

    p['verbosity'] = 'DEBUG'


    # Range of E/I synaptic conductances
    Nvals  = 31      # Number of values for each dimension
    startG = 0.0     # nS
    endG   = 6120.0  # nS

    ###############################################################################

    submitParamSweep(p, startG, endG, Nvals, ENV, simRootDir, simLabel,
            appName, rtLimit, numCPU, blocking, timePrefix, numRepeat, dry_run)
