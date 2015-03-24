#!/usr/bin/env python
#
'''
Submit jobs to print some figures for each parameter sweep job.
'''
from __future__ import absolute_import, print_function

import numpy as np

from grid_cell_model.submitting.factory import SubmitterFactory
from grid_cell_model.submitting.arguments import ArgumentCreator
from grid_cell_model.submitting import flagparse
from grid_cell_model.submitting.flagparse import positive_int

import common

parser = flagparse.FlagParser()
parser.add_argument('env',       type=str, choices=['workstation', 'cluster'])
parser.add_argument("where",     type=str)
parser.add_argument('figure_dir',type=str)
parser.add_argument('type',      type=str, choices=common.allowed_types, nargs="+")
parser.add_argument('--row',     type=int)
parser.add_argument('--col',     type=int)
parser.add_argument("--ns",      type=int, choices=[0, 150, 300])
parser.add_argument('--nCPU',    type=positive_int, default=1)
parser.add_argument('--rtLimit', type=str, default='00:05:00')
parser.add_flag("--ns_all")
parser.add_flag("--ignoreErrors")
o = parser.parse_args()

if not o.ns_all and o.ns is None:
    raise RuntimeError("Must specify either --ns or --ns_all!")

# Submitting
ENV         = o.env
appName     = 'figure_sweeps.py'
rtLimit     = o.rtLimit
numCPU      = o.nCPU
blocking    = True
timePrefix  = False
numRepeat   = 1
dry_run     = False


ns_all = [0, 150, 300]
shape = (61, 11)
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
    p['figure_dir']  = o.figure_dir

    ##########################################################################
    ac = ArgumentCreator(p, printout=True)

    iterparams = {
        'row' : np.arange(rowN) if o.row is None else [o.row],
        'col' : np.arange(colN) if o.col is None else [o.col]
    }
    ac.insertDict(iterparams, mult=True)

    ##########################################################################
    submitter = SubmitterFactory.getSubmitter(
        ac, appName, envType=ENV, rtLimit=rtLimit, output_dir=simRootDir,
        label=simLabel, blocking=blocking, timePrefix=timePrefix,
        numCPU=numCPU, ignoreSubmitErrors=o.ignoreErrors)
    ac.setOption('output_dir', submitter.outputDir())
    startJobNum = 0
    submitter.submitAll(startJobNum, numRepeat, dry_run=dry_run)
