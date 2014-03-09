#!/usr/bin/env python
'''
Parameter sweeps of grid field simulations in which detailed noise is one of
the parameters.
'''
import numpy as np
from default_params         import defaultParameters as dp
from param_sweep_noise      import submitNoiseSweep, SweepParams, \
        getBumpCurrentSlope
from submitting             import flagparse
from submitting.flagparse   import positive_int

allowedPositions = ['EI-1_3', 'EI-3_1']

parser = flagparse.FlagParser()
parser.add_argument("--where",   type=str, required=True)
parser.add_argument('--env',     type=str, choices=['workstation', 'cluster'], required=True)
parser.add_argument('--nCPU',    type=positive_int, default=1)
parser.add_argument('--rtLimit', type=str, default='05:00:00')
parser.add_argument('--position',type=str, choices=allowedPositions,
        required=True)
parser.add_flag("--forceUpdate")
parser.add_flag("--dry-run")
o = parser.parse_args()


p = dp.copy()
# Submitting
ENV         = o.env
simRootDir  = o.where
simLabel    = o.position
appName     = 'simulation_grids.py'
rtLimit     = o.rtLimit
numCPU      = o.nCPU
blocking    = True
timePrefix  = False
numRepeat   = 1
dry_run     = o.dry_run

p['time']              = 600e3  # ms
p['nthreads']          = 1
p['ntrials']           = 1
p['velON']             = 1
p['constantPosition']  = 0
p['verbosity']         = o.verbosity

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

extraIterparams = dict(
        bumpCurrentSlope=getBumpCurrentSlope(simLabel, threshold=-np.infty))
#extraIterparams = dict(bumpCurrentSlope=[1.])

###############################################################################
submitNoiseSweep(p, gEParams, gIParams, noiseParams,
        ENV, simRootDir, simLabel, appName, rtLimit, numCPU, blocking,
        timePrefix, numRepeat, dry_run, extraIterparams=extraIterparams)

