#!/usr/bin/env python
#
#   submit_connections.py
#
#   Submit job(s) to the cluster/workstation: export connections to file
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
from submitting.factory   import SubmitterFactory
from submitting.arguments import ArgumentCreator
from default_params       import decideParams
import logging as lg
#lg.basicConfig(level=lg.DEBUG)
lg.basicConfig(level=lg.INFO)

#simLabel    = 'E_surround'
simLabel    = 'I_surround'

p = decideParams(simLabel)
# Submitting
ENV         = 'workstation'
simRootDir  = 'output_local/connections'
appName     = '../common/simulation_connections.py'
numCPU      = 1
blocking    = True
timePrefix  = False
numRepeat   = 1
dry_run     = False

p['time']     = 0.1 # unused
p['nthreads'] = 1
p['ntrials']  = 1

###############################################################################
ac = ArgumentCreator(p, printout=True)
submitter = SubmitterFactory.getSubmitter(ac, appName, envType=ENV,
        output_dir=simRootDir, label=simLabel, blocking=blocking,
        timePrefix=timePrefix, numCPU=numCPU)
ac.setOption('output_dir', submitter.outputDir())
startJobNum = 0
submitter.submitAll(startJobNum, numRepeat, dry_run=dry_run)
