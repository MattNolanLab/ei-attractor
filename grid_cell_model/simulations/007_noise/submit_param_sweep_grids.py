#!/usr/bin/env python
#
#   submit_basic_grids.py
#
#   Submit job(s) to the cluster/workstation
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
from submitting.factory   import SubmitterFactory
from submitting.arguments import ArgumentCreator
from default_params       import defaultParameters as p
from param_sweep          import submitParamSweep, getBumpCurrentSlope
import logging as lg
#lg.basicConfig(level=lg.DEBUG)
lg.basicConfig(level=lg.INFO)

p['noise_sigma']       = 150.0     # pA

# Submitting
ENV         = 'cluster'
simRootDir  = 'output/grids'
simLabel    = 'EI_param_sweep_{0}pA'.format(int(p['noise_sigma']))
appName     = 'simulation_grids.py'
rtLimit     = '14:00:00'
numCPU      = 1
blocking    = True
timePrefix  = False
numRepeat   = 1
dry_run     = False

p['time']              = 600e3  # ms
p['nthreads']          = 1
p['ntrials']           = 6
p['velON']             = 1


# Range of parameters around default values
Nvals        = 30    # Number of values for each dimension
startFrac    = 0.
endFrac      = 2.8572

extraIterparams = {'bumpCurrentSlope' : getBumpCurrentSlope(p['noise_sigma'],
    threshold=0.05)}
#extraIterparams['bumpCurrentSlope'] = [1.0]

###############################################################################

submitParamSweep(p, startFrac, endFrac, Nvals, ENV, simRootDir, simLabel,
        appName, rtLimit, numCPU, blocking, timePrefix, numRepeat, dry_run,
        extraIterparams)

