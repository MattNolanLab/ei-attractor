#!/usr/bin/env python
#
#   submit_details_EI.py
#
#   Submit job(s) to the cluster/workstation: parameter sweep detailed plots
#   (noise)
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
from submitting.factory     import SubmitterFactory
from submitting.arguments   import ArgumentCreator
import logging as lg
lg.basicConfig(level=lg.DEBUG)

# Submitting
ENV         = 'cluster'
simRootDir  = 'output/no_theta/gamma_bump'
appName     = 'details_EI.py'
rtLimit     = '00:02:30'
numCPU      = 1
blocking    = False
timePrefix  = False
numRepeat   = 1
dry_run     = False

noise_sigma_all = [0.0, 150.0, 300.0] # pA
dirs = ('{0}pA'         , (31, 31))

for noise_sigma in noise_sigma_all:
    p = {}
    simLabel   = dirs[0].format(int(noise_sigma))
    rowN       = dirs[1][0]
    colN       = dirs[1][1]

    p['shapeRows'] = rowN
    p['shapeCols'] = colN

    ###############################################################################

    ac = ArgumentCreator(p, printout=True)

    iterparams = {
            'row' : np.arange(rowN),
            'col' : np.arange(colN)
            #'row' : np.arange(1),
            #'col' : np.arange(1)
    }
    ac.insertDict(iterparams, mult=True)

    ###############################################################################
    submitter = SubmitterFactory.getSubmitter(ac, appName, envType=ENV,
            rtLimit=rtLimit, output_dir=simRootDir, label=simLabel,
            blocking=blocking, timePrefix=timePrefix, numCPU=numCPU)
    ac.setOption('output_dir', submitter.outputDir())
    startJobNum = 0
    submitter.submitAll(startJobNum, numRepeat, dry_run=dry_run)
