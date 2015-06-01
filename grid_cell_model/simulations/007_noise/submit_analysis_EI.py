#!/usr/bin/env python
#
'''
Submit EI analysis jobs.
'''
from __future__ import absolute_import, print_function

import numpy as np

from grid_cell_model.submitting.factory   import SubmitterFactory
from grid_cell_model.submitting.arguments import ArgumentCreator
from grid_cell_model.submitting           import flagparse
from grid_cell_model.submitting.flagparse import positive_int

import common.analysis as common
from param_sweep          import getSpeedPercentile
from default_params       import defaultParameters as dp



parser = flagparse.FlagParser()
parser.add_argument('env',       type=str, choices=['workstation', 'cluster'])
parser.add_argument("where",     type=str)
parser.add_argument('type',      type=str, choices=common.allowedTypes, nargs="+")
parser.add_argument('--row',     type=int)
parser.add_argument('--col',     type=int)
parser.add_argument("--ns",      type=int, choices=[0, 150, 300])
parser.add_argument('--nCPU',    type=positive_int, default=1)
parser.add_argument('--rtLimit', type=str, default='00:05:00')
parser.add_argument('--shape',   type=positive_int, nargs=2,
                    help="Sweep shape (rows, columns)")
parser.add_flag("--ns_all")
parser.add_flag("--forceUpdate")
parser.add_flag("--ignoreErrors")
o = parser.parse_args()

if not o.ns_all and o.ns is None:
    raise RuntimeError("Must specify either --ns or --ns_all!")


# Submitting
ENV         = o.env
appName     = 'analysis_EI.py'
rtLimit     = o.rtLimit
numCPU      = o.nCPU
blocking    = True
timePrefix  = False
numRepeat   = 1
dry_run     = False


ns_all = [0, 150, 300]
shape = (31, 31) if o.shape is None else o.shape
noise_sigmas = ns_all if o.ns_all  else [o.ns]


for noise_sigma in noise_sigmas:
    p = {}
    simRootDir = o.where
    simLabel   = '{0}pA'.format(int(noise_sigma))
    rowN       = shape[0]
    colN       = shape[1]

    p['type']        = o.type
    p['shapeRows']   = rowN
    p['shapeCols']   = colN
    p['verbosity']   = o.verbosity
    p['forceUpdate'] = int(o.forceUpdate)

    if common.velocityType in o.type:
        percentile = 99.0
        p['bumpSpeedMax'] = getSpeedPercentile(percentile, dp['ratVelFName'],
                dp['gridSep'], dp['Ne'])

    ###############################################################################

    ac = ArgumentCreator(p, printout=True)

    iterparams = {
            'row' : np.arange(rowN) if o.row is None else [o.row],
            'col' : np.arange(colN) if o.col is None else [o.col]
    }
    ac.insertDict(iterparams, mult=True)

    ###############################################################################
    submitter = SubmitterFactory.getSubmitter(ac, appName, envType=ENV,
            rtLimit=rtLimit, output_dir=simRootDir, label=simLabel,
            blocking=blocking, timePrefix=timePrefix, numCPU=numCPU,
            ignoreSubmitErrors=o.ignoreErrors)
    ac.setOption('output_dir', submitter.outputDir())
    startJobNum = 0
    submitter.submitAll(startJobNum, numRepeat, dry_run=dry_run)
