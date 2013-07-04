#!/usr/bin/env python
#
#   submit_single_neuron.py
#
#   Run a simple, single-neuron-from-each-population simulation
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
from default_params import defaultParameters as p
from submitting.factory     import SubmitterFactory
from submitting.arguments   import ArgumentCreator
import logging as lg
#lg.basicConfig(level=lg.DEBUG)
lg.basicConfig(level=lg.INFO)


# Submitting
ENV         = 'workstation'
simRootDir  = 'output_local'
simLabel    = 'single_neuron'
appName     = 'simulation_single_neuron.py'
rtLimit     = '00:02:00'
numCPU      = 1
blocking    = True
timePrefix  = False
numRepeat   = 1
dry_run     = False

p['time']              = 10e3  # ms
p['nthreads']          = 1


###############################################################################
ac = ArgumentCreator(p, printout=True)


iterparams = {
        'noise_sigma' : [0.0, 150.0, 300.0] # pA
}
ac.insertDict(iterparams, mult=False)

###############################################################################
submitter = SubmitterFactory.getSubmitter(ac, appName, envType=ENV,
        rtLimit=rtLimit, output_dir=simRootDir, label=simLabel,
        blocking=blocking, timePrefix=timePrefix, numCPU=numCPU)
ac.setOption('output_dir', submitter.outputDir())
startJobNum = 0
submitter.submitAll(startJobNum, numRepeat, dry_run=dry_run)
submitter.saveIterParams(iterparams, dry_run=dry_run)
