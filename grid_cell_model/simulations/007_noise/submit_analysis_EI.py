#!/usr/bin/env python
#
'''
Submit EI analysis jobs.
'''

import numpy as np
from submitting.factory   import SubmitterFactory
from submitting.arguments import ArgumentCreator
from param_sweep          import getSpeedPercentile
from default_params       import defaultParameters as dp
from submitting           import flagparse
from submitting.flagparse import positive_int

velocityType  = 'velocity'

parser = flagparse.FlagParser()
parser.add_argument("--where",      type=str, required=True)
parser.add_argument("--ns",         type=int, choices=[0, 150, 300])
parser.add_argument('--type',       type=str, choices=[velocityType], required=True)
parser.add_argument('--env',        type=str, choices=['workstation', 'cluster'], required=True)
parser.add_argument('--nCPU',       type=positive_int, default=1)
parser.add_argument('--rtLimit',    type=str, default='00:05:00')
parser.add_flag("--ns_all")
parser.add_flag("--forceUpdate")
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
shape = (31, 31)
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

    if (p['type'] == velocityType):
        percentile = 99.0
        p['bumpSpeedMax'] = getSpeedPercentile(percentile, dp['ratVelFName'],
                dp['gridSep'], dp['Ne'])

    ###############################################################################

    ac = ArgumentCreator(p, printout=True)

    iterparams = {
            'row' : np.arange(rowN),
            'col' : np.arange(colN)
            #'row' : [5],
            #'col' : [15]
    }
    ac.insertDict(iterparams, mult=True)

    ###############################################################################
    submitter = SubmitterFactory.getSubmitter(ac, appName, envType=ENV,
            rtLimit=rtLimit, output_dir=simRootDir, label=simLabel,
            blocking=blocking, timePrefix=timePrefix, numCPU=numCPU)
    ac.setOption('output_dir', submitter.outputDir())
    startJobNum = 0
    submitter.submitAll(startJobNum, numRepeat, dry_run=dry_run)
