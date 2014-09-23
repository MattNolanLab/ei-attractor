#!/usr/bin/env python
'''Run a simple, single-neuron-from-each-population simulation.'''
from __future__ import absolute_import, print_function, division

import logging as lg
#lg.basicConfig(level=lg.DEBUG)
lg.basicConfig(level=lg.INFO)

from grid_cell_model.submitting.factory   import SubmitterFactory
from grid_cell_model.submitting.arguments import ArgumentCreator

from default_params import defaultParameters as p

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
