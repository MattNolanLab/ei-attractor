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
import logging as lg
lg.basicConfig(level=lg.DEBUG)

# Submitting
ENV         = 'cluster'
appName     = 'analysis_EI.py'
rtLimit     = '00:05:00'
numCPU      = 1
blocking    = True
timePrefix  = False
numRepeat   = 1
dry_run     = False

gammaBumpType = 'gamma-bump'
velocityType  = 'velocity'
gridsType     = 'grids'
posType       = 'positional'

noise_sigma_all = [0, 150.0, 300.0] # pA
dirs = \
    ('output/even_spacing/velocity_vertical', velocityType,  '{0}pA', (31, 31))
    #('output/even_spacing/const_position',    posType,       '{0}pA', (31, 31))
    #('output/even_spacing/gamma_bump',        gammaBumpType, '{0}pA', (31, 31))
    #('output/even_spacing/grids',             gridsType,     '{0}pA', (31, 31))
    #('output/even_spacing/grids_no_velocity', gridsType,     '{0}pA', (31, 31))
    #('output/even_spacing/velocity',          velocityType,  '{0}pA', (31, 31))
    #('output/no_theta/grids',                 gridsType,     '{0}pA', (31, 31))
    #('output/no_theta/velocity',              velocityType,  '{0}pA', (31, 31))
    #('output/no_theta/gamma_bump',            gammaBumpType, '{0}pA', (31, 31))

for noise_sigma in noise_sigma_all:
    p = {}
    simRootDir = dirs[0]
    p['type']  = dirs[1]
    simLabel   = dirs[2].format(int(noise_sigma))
    rowN       = dirs[3][0]
    colN       = dirs[3][1]
    p['verbosity'] = 'INFO'

    p['shapeRows'] = rowN
    p['shapeCols'] = colN
    p['forceUpdate'] = 0

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
