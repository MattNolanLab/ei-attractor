#!/usr/bin/env python
from __future__ import absolute_import, print_function
'''Submit EI analysis jobs - detailed noise levels.'''
import numpy as np

from grid_cell_model.submitting import flagparse
from grid_cell_model.submitting.flagparse import positive_int
from grid_cell_model.submitting.factory   import SubmitterFactory
from grid_cell_model.submitting.arguments import ArgumentCreator

from param_sweep          import getSpeedPercentile
from default_params       import defaultParameters as dp
import common.analysis as common

allowedPositions = ['EI-1_3', 'EI-3_1']

parser = flagparse.FlagParser()
parser.add_argument('--row',     type=int)
parser.add_argument('--col',     type=int)
parser.add_argument("--where",   type=str, required=True)
parser.add_argument('--type',    type=str, choices=common.allowedTypes, required=True)
parser.add_argument('--env',     type=str, choices=['workstation', 'cluster'], required=True)
parser.add_argument('--nCPU',    type=positive_int, default=1)
parser.add_argument('--rtLimit', type=str, default='00:05:00')
parser.add_argument('--position',type=str, choices=allowedPositions)
parser.add_flag("--all-positions")
parser.add_flag("--forceUpdate")
parser.add_flag("--ignoreErrors")
o = parser.parse_args()

if not o.all_positions and o.position is None:
    raise RuntimeError("Must specify either --position or --all-positions")


# Submitting
ENV         = o.env
appName     = 'analysis_EI.py'
rtLimit     = o.rtLimit
numCPU      = o.nCPU
blocking    = True
timePrefix  = False
numRepeat   = 1
dry_run     = False


shape = (31, 9)
positions = allowedPositions if o.all_positions else [o.position]


for position in positions:
    p = {}
    simRootDir = o.where
    simLabel   = position
    rowN       = shape[0]
    colN       = shape[1]

    p['type']        = o.type
    p['shapeRows']   = rowN
    p['shapeCols']   = colN
    p['verbosity']   = o.verbosity
    p['forceUpdate'] = int(o.forceUpdate)

    if (p['type'] == common.velocityType):
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
