#!/usr/bin/env python
#
#   submit_basic_grids.py
#
#   Submit job(s) to the cluster/workstation
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
from submitting.arguments import ArgumentCreator
from default_params       import defaultParameters as p
from param_sweep          import submitParamSweep
import logging as lg
#lg.basicConfig(level=lg.DEBUG)
lg.basicConfig(level=lg.INFO)

p['noise_sigma']       = 150.0     # pA

# Submitting
ENV         = 'cluster'
simRootDir  = 'output/grids_no_velocity'
simLabel    = 'grids_velocity_ON_{0}pA'.format(int(p['noise_sigma']))
appName     = 'simulation_grids.py'
rtLimit     = '02:00:00'
numCPU      = 4
blocking    = True
timePrefix  = False
numRepeat   = 10
dry_run     = False

p['time']              = 1200e3  # ms
p['nthreads']          = 4
p['ntrials']           = 1

p['bumpCurrentSlope']  = 0.53 # neurons/s/pA, !! this will depend on prefDirC !!
p['pc_conn_weight']    = 0.5

ac = ArgumentCreator(p, printout=True)

iterparams = {
        'pc_max_rate'   : np.arange(0, 100, 10),
}
ac.insertDict(iterparams, mult=True)

###############################################################################
submitter = SubmitterFactory.getSubmitter(ac, appName, envType=ENV,
        rtLimit=rtLimit, output_dir=simRootDir, label=simLabel,
        blocking=blocking, timePrefix=timePrefix, numCPU=numCPU)
ac.setOption('output_dir', submitter.outputDir())
startJobNum = 0
submitter.submitAll(startJobNum, numRepeat, dry_run=dry_run)
submitter.saveIterParams(iterparams, dry_run=dry_run)
