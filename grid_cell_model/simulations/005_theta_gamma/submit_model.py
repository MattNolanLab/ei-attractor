#!/usr/bin/env python
#
#   submit_model.py
#
#   Submit job(s) to the cluster/workstation: Simple stationary bump.
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

simLabel    = 'E_surround'
#simLabel    = 'I_surround'

p = decideParams(simLabel)
# Submitting
ENV         = 'workstation'
simRootDir  = 'output/model'
appName     = '../common/simulation_stationary.py'
rtLimit     = '00:45:00'
numCPU      = 8
blocking    = True
timePrefix  = False
numRepeat   = 1
dry_run     = False

p['time']              = 6e3  # ms
p['nthreads']          = 8
p['ntrials']           = 5


###############################################################################
ac = ArgumentCreator(p, printout=True)
submitter = SubmitterFactory.getSubmitter(ac, appName, envType=ENV,
        rtLimit=rtLimit, output_dir=simRootDir, label=simLabel,
        blocking=blocking, timePrefix=timePrefix, numCPU=numCPU)
ac.setOption('output_dir', submitter.outputDir())
startJobNum = 0
submitter.submitAll(startJobNum, numRepeat, dry_run=dry_run)
