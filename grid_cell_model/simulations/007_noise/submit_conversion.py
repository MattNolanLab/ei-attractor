#!/usr/bin/env python
#
'''
Submit data conversion jobs for the whole parameter sweep.
'''

import numpy as np
from submitting.factory   import SubmitterFactory
from submitting.arguments import ArgumentCreator
from param_sweep          import getSpeedPercentile
from default_params       import defaultParameters as dp
from submitting import flagparse

import logging as lg
lg.basicConfig(level=lg.DEBUG)

parser = flagparse.FlagParser()
parser.add_argument("--where",      type=str, required=True)
parser.add_argument("--ns",         type=int, choices=[0, 150, 300])
parser.add_argument("--what",       type=str, required=True)
parser.add_argument('--type',       type=str, choices=['velocity'], required=True)
parser.add_argument('--targetType', type=str, required=True)
parser.add_argument('--env',        type=str, choices=['workstation', 'cluster'], required=True)
parser.add_flag("--ns_all")
parser.add_flag("--repack", help='Whether to repack data after the operation')
o = parser.parse_args()

if not o.ns_all and o.ns is None:
    raise RuntimeError("Must specify either --ns or --ns_all!")


# Submitting
ENV         = o.env
appName     = 'convert_data.py'
rtLimit     = '00:05:00'
numCPU      = 1
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

    p['where']      = o.where
    p['ns']         = int(noise_sigma)
    p['what']       = o.what
    p['type']       = o.type
    p['targetType'] = o.targetType
    p['verbosity']  = o.verbosity
    if o.repack:
        p['repack'] = ''


    ###############################################################################

    ac = ArgumentCreator(p, printout=True, emitJobNum=False)

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
    startJobNum = 0
    submitter.submitAll(startJobNum, numRepeat, dry_run=dry_run)
