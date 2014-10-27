#!/usr/bin/env python
#
'''
Detailed noise simulations of gamma/bump parameter sweep. Noise levels are one
of the parameters in the simulation (0-300 pA in steps of 10 pA).
'''
from default_params    import defaultParameters as dp
from param_sweep_noise import submitNoiseSweep, SweepParams
import logging as lg
lg.basicConfig(level=lg.DEBUG)
#lg.basicConfig(level=lg.INFO)

# We are expecting 1-3 distinct simulation runs
#simLabel = 'EI-1_3'
simLabel = 'EI-3_1'

p = dp.copy()
# Submitting
ENV         = 'cluster'
simRootDir  = 'output/detailed_noise_vertical/gamma_bump'
appName     = '../common/simulation_stationary.py'
rtLimit     = '01:00:00'
numCPU      = 1
blocking    = True
timePrefix  = False
numRepeat   = 1
dry_run     = False

p['time']              = 10e3  # ms
p['nthreads']          = 1
p['ntrials']           = 5


# Range of noise and E/I synaptic conductances
noiseParams = SweepParams(0, 300, 31)
if (simLabel == 'EI-1_3'):
    gEParams = SweepParams(816, 1224, 3)
    gIParams = SweepParams(2856, 3264, 3)
elif (simLabel == 'EI-3_1'):
    gEParams = SweepParams(2856, 3264, 3)
    gIParams = SweepParams(816, 1224, 3)
else:
    raise ValueError('Unknown simLabel!')

###############################################################################
submitNoiseSweep(p, gEParams, gIParams, noiseParams,
        ENV, simRootDir, simLabel, appName, rtLimit, numCPU, blocking,
        timePrefix, numRepeat, dry_run)

