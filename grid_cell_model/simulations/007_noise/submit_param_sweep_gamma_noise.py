#!/usr/bin/env python
#
#   submit_param_sweep_gamma_noise.py
#
#   Submit job(s) to the cluster/workstation: gamma parameter sweep; detailed
#   noise levels.
#
#       Copyright (C) 2012  Lukas Solanka <l.solanka@sms.ed.ac.uk>
#       
#       This program is free software: you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation, either version 3 of the License, or
#       (at your option) any later version.
#       
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#       
#       You should have received a copy of the GNU General Public License
#       along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
from default_params    import defaultParameters as dp
from param_sweep_noise import submitNoiseSweep, SweepParams
import logging as lg
#lg.basicConfig(level=lg.DEBUG)
lg.basicConfig(level=lg.INFO)

# We are expecting 1-3 distinct simulation runs
#simLabel = 'EI-1_3'
simLabel = 'EI-3_1'

p = dp.copy()
# Submitting
ENV         = 'cluster'
simRootDir  = 'output/detailed_noise/gamma_bump'
appName     = 'simulation_stationary.py'
rtLimit     = '00:45:00'
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
if (simLabel == 'EI-3_1'):
    gEParams = SweepParams(2856, 3264, 3)
    gIParams = SweepParams(816, 1224, 3)
else:
    raise ValueError('Unknown simLabel!')

###############################################################################
submitNoiseSweep(p, gEParams, gIParams, noiseParams,
        ENV, simRootDir, simLabel, appName, rtLimit, numCPU, blocking,
        timePrefix, numRepeat, dry_run)

