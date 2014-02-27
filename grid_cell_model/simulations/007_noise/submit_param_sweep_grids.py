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
import logging as lg
import numpy as np

from submitting.factory   import SubmitterFactory
from submitting.arguments import ArgumentCreator
from default_params       import defaultParameters as dp
from submitting           import flagparse
from submitting.flagparse import positive_int
from param_sweep          import submitParamSweep, getBumpCurrentSlope

parser = flagparse.FlagParser()
parser.add_argument('--row',     type=int)
parser.add_argument('--col',     type=int)
parser.add_argument("--where",      type=str, required=True)
parser.add_argument("--ns",         type=int, choices=[0, 150, 300])
parser.add_argument('--ntrials',    type=positive_int, default=1)
parser.add_argument('--rtLimit',    type=str, default='05:00:00')
parser.add_argument('--env',        type=str, choices=['workstation', 'cluster'], required=True)
parser.add_flag('--dry_run', help='Do no run anything nor save any meta-data')
o = parser.parse_args()

if (o.row is None) ^ (o.col is None):
    raise ValueError("Specify either both --row and --col or None!")

ns_all = [0.0, 150.0, 300.0] # pA
noise_sigmas = ns_all if o.ns is None  else [o.ns]

for noise_sigma in noise_sigmas:
    p = dp.copy()
    p['noise_sigma'] = noise_sigma # pA

    # Submitting
    ENV         = o.env
    simRootDir  = o.where
    simLabel    = '{0}pA'.format(int(p['noise_sigma']))
    appName     = 'simulation_grids.py'
    rtLimit     = o.rtLimit
    numCPU      = 1
    blocking    = True
    timePrefix  = False
    numRepeat   = 1
    dry_run     = o.dry_run

    p['time']              = 600e3  # ms
    p['nthreads']          = 1
    p['ntrials']           = o.ntrials
    p['velON']             = 1
    p['constantPosition']  = 0
    p['verbosity']         = o.verbosity


    # Range of E/I synaptic conductances
    Nvals  = 31      # Number of values for each dimension
    startG = 0.0     # nS
    endG   = 6120.0  # nS

    extraIterparams = {'bumpCurrentSlope' : getBumpCurrentSlope(p['noise_sigma'],
        threshold=-np.infty)}
    #extraIterparams['bumpCurrentSlope'] = [1.0]

    ###############################################################################
    rc = (o.row, o.col) if o.row is not None else None

    submitParamSweep(p, startG, endG, Nvals, ENV, simRootDir, simLabel,
            appName, rtLimit, numCPU, blocking, timePrefix, numRepeat, dry_run,
            extraIterparams, rc=rc)

