#!/usr/bin/env python
#
#   submit_analysis_EI.py
#
#   Submit job(s) to the cluster/workstation: parameter sweep analysis (noise)
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
from param_sweep          import getSpeedPercentile
from default_params       import defaultParameters as dp
import logging as lg
lg.basicConfig(level=lg.DEBUG)

# Submitting
ENV         = 'cluster'
appName     = 'analysis_EI.py'
rtLimit     = '00:02:00'
numCPU      = 1
blocking    = True
timePrefix  = False
numRepeat   = 1
dry_run     = False

gammaBumpType = 'gamma-bump'
velocityType = 'velocity'
gridsType = 'grids'

noise_sigma_all = [0.0, 150.0, 300.0] # pA
dirs = \
    ('output/no_theta/gamma_bump', gammaBumpType, '{0}pA', (31, 31))
    #('output/grids',      gridsType,     '{0}pA', (30, 30))
    #('output/velocity',   velocityType,  '{0}pA', (30, 30))

for noise_sigma in noise_sigma_all:
    p = {}
    simRootDir = dirs[0]
    p['type']  = dirs[1]
    simLabel   = dirs[2].format(int(noise_sigma))
    rowN       = dirs[3][0]
    colN       = dirs[3][1]

    p['shapeRows'] = rowN
    p['shapeCols'] = colN
    p['forceUpdate'] = 0

    ###############################################################################

    ac = ArgumentCreator(p, printout=True)

    iterparams = {
            'row' : np.arange(rowN),
            'col' : np.arange(colN)
            #'row' : [10],
            #'col' : [10]
    }
    ac.insertDict(iterparams, mult=True)

    ###############################################################################
    submitter = SubmitterFactory.getSubmitter(ac, appName, envType=ENV,
            rtLimit=rtLimit, output_dir=simRootDir, label=simLabel,
            blocking=blocking, timePrefix=timePrefix, numCPU=numCPU)
    ac.setOption('output_dir', submitter.outputDir())
    startJobNum = 0
    submitter.submitAll(startJobNum, numRepeat, dry_run=dry_run)
