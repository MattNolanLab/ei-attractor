#!/usr/bin/env python
'''
Detailed noise simulations of gamma/bump parameter sweep. Noise levels are one
of the parameters in the simulation (0-300 pA in steps of 10 pA).
'''
from __future__ import absolute_import, print_function, division

from grid_cell_model.submitting.base.parsers import GenericSubmissionParser
from default_params import defaultParameters as dp
from param_sweep_noise import submitNoiseSweep, SweepParams


allowedPositions = ['EI-1_3', 'EI-3_1']

parser = GenericSubmissionParser()
parser.add_argument('--position',type=str, choices=allowedPositions,
                    required=True)

o = parser.parse_args()
p = dp.copy()

# Submitting
ENV         = o.env
simRootDir  = o.where
simLabel    = o.position
appName     = '../common/simulation_stationary.py'
rtLimit     = o.rtLimit or '02:30:00'
numCPU      = o.nCPU
blocking    = True
timePrefix  = False
numRepeat   = 1
dry_run     = o.dry_run

p['master_seed'] = 123456
p['time']        = o.time or 10e3  # ms
p['nthreads']    = 1
p['ntrials']     = o.ntrials or 5


# Range of noise and E/I synaptic conductances
noiseParams = SweepParams(0, 300, 31)
if simLabel == 'EI-1_3':
    gEParams = SweepParams(816, 1224, 3)
    gIParams = SweepParams(2856, 3264, 3)
elif simLabel == 'EI-3_1':
    gEParams = SweepParams(2856, 3264, 3)
    gIParams = SweepParams(816, 1224, 3)
else:
    raise ValueError('Unknown simLabel!')

###############################################################################
submitNoiseSweep(p, gEParams, gIParams, noiseParams,
                 ENV, simRootDir, simLabel, appName, rtLimit, numCPU, blocking,
                 timePrefix, numRepeat, dry_run)
