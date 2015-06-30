#!/usr/bin/env python
'''Bump velocity gain estimation - detailed noise levels.

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
from grid_cell_model.submitting.base.parsers import GenericSubmissionParser
from param_sweep_noise import submitNoiseSweep, SweepParams
from default_params import defaultParameters as dp

parser = GenericSubmissionParser()
parser.add_argument('--EI_type', type=str, required=True,
                    choices=['EI-1_3', 'EI-3_1'],
                    help='Value of E and I coupling (g_AMPA_total and '
                         'g_GABA_total)')
o = parser.parse_args()


p = dp.copy()

# Submitting
ENV         = o.env
simRootDir  = o.where
simLabel    = o.EI_type
appName     = 'simulation_velocity.py'
rtLimit     = o.rtLimit or '12:00:00'
numCPU      = 1
blocking    = True
timePrefix  = False
numRepeat   = 1
dry_run     = o.dry_run

p['master_seed'] = 123456
p['time']        = o.time or 10e3  # ms
p['nthreads']    = 1
p['ntrials']     = 10

p['IvelMax']     = 100
p['dIvel']       = 10

p['verbosity']   = o.verbosity


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
submitNoiseSweep(p, gEParams, gIParams, noiseParams, ENV, simRootDir, simLabel,
                 appName, rtLimit, numCPU, blocking, timePrefix, numRepeat,
                 dry_run)

