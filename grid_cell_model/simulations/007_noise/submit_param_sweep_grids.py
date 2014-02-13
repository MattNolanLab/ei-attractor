#!/usr/bin/env python
#
#   submit_param_sweep_grids.py
#
#   Submit job(s) to the cluster/workstation: grid field parameter sweeps
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
from default_params       import defaultParameters as dp
from param_sweep          import submitParamSweep, getBumpCurrentSlope
import logging as lg
#lg.basicConfig(level=lg.DEBUG)
lg.basicConfig(level=lg.INFO)

noise_sigma_all = [0.0, 150.0, 300.0] # pA

for noise_sigma in noise_sigma_all:
    p = dp.copy()
    p['noise_sigma'] = noise_sigma # pA

    # Submitting
    ENV         = 'cluster'
    simRootDir  = 'output/even_spacing/grids'
    simLabel    = '{0}pA'.format(int(p['noise_sigma']))
    appName     = 'simulation_grids.py'
    rtLimit     = '05:00:00'
    numCPU      = 1
    blocking    = True
    timePrefix  = False
    numRepeat   = 1
    dry_run     = False

    p['time']              = 600e3  # ms
    p['nthreads']          = 1
    p['ntrials']           = 3
    p['velON']             = 1
    p['constantPosition']  = 0


    # Range of E/I synaptic conductances
    Nvals  = 31      # Number of values for each dimension
    startG = 0.0     # nS
    endG   = 6120.0  # nS

    extraIterparams = {'bumpCurrentSlope' : getBumpCurrentSlope(p['noise_sigma'],
        threshold=-np.infty)}
    #extraIterparams['bumpCurrentSlope'] = [1.0]

    ###############################################################################

    submitParamSweep(p, startG, endG, Nvals, ENV, simRootDir, simLabel,
            appName, rtLimit, numCPU, blocking, timePrefix, numRepeat, dry_run,
            extraIterparams)

