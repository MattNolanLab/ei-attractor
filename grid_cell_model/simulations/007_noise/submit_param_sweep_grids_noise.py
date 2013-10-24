#!/usr/bin/env python
#
#   submit_param_sweep_grids_noise.py
#
#   Submit job(s) to the cluster/workstation: grid field parameter sweeps
#   (detailed noise)
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
import numpy as np
from default_params import defaultParameters as dp
from param_sweep_noise import submitNoiseSweep, SweepParams, \
        getBumpCurrentSlope
import logging as lg
#lg.basicConfig(level=lg.DEBUG)
lg.basicConfig(level=lg.INFO)

# We are expecting 1-3 distinct simulation runs
#simLabel = 'EI-1_3'
simLabel = 'EI-3_1'

p = dp.copy()
# Submitting
ENV         = 'cluster'
simRootDir  = 'output/detailed_noise/grids'
appName     = 'simulation_grids.py'
rtLimit     = '05:00:00'
numCPU      = 1
blocking    = True
timePrefix  = False
numRepeat   = 1
dry_run     = False

p['time']              = 600e3  # ms
p['nthreads']          = 1
p['ntrials']           = 1
p['velON']             = 1

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

###############################################################################
submitNoiseSweep(p, gEParams, gIParams, noiseParams,
        ENV, simRootDir, simLabel, appName, rtLimit, numCPU, blocking,
        timePrefix, numRepeat, dry_run, extraIterparams=extraIterparams)

