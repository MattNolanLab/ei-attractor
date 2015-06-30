#!/usr/bin/env python
'''Run a simple, single-neuron-from-each-population simulation.'''
from __future__ import absolute_import, print_function, division

from grid_cell_model.submitting.base.parsers import GenericSubmissionParser
from grid_cell_model.submitting.arguments import ArgumentCreator
from grid_cell_model.submitting.factory import SubmitterFactory
from default_params import defaultParameters as dp

parser = GenericSubmissionParser()
o = parser.parse_args()

p = dp.copy()

# Submitting
ENV         = o.env
simRootDir  = o.where
simLabel    = 'single_neuron'
appName     = 'simulation_single_neuron.py'
rtLimit     = '00:02:00' if o.rtLimit is None else o.rtLimit
numCPU      = 1
blocking    = True
timePrefix  = False
numRepeat   = 1
dry_run     = o.dry_run

p['master_seed'] = 123456
p['time']        = 10e3 if o.time is None else o.time  # ms
p['nthreads']    = 1
p['ntrials']     = 5 if o.ntrials is None else o.ntrials
p['verbosity']   = o.verbosity

###############################################################################
ac = ArgumentCreator(p, printout=True)

iterparams = {
    'noise_sigma' : [0.0, 150.0, 300.0] # pA
}
ac.insertDict(iterparams, mult=False)

###############################################################################
submitter = SubmitterFactory.getSubmitter(ac, appName, envType=ENV,
                                          rtLimit=rtLimit,
                                          output_dir=simRootDir,
                                          label=simLabel, blocking=blocking,
                                          timePrefix=timePrefix, numCPU=numCPU)
ac.setOption('output_dir', submitter.outputDir())
startJobNum = 0
submitter.submitAll(startJobNum, numRepeat, dry_run=dry_run)
submitter.saveIterParams(iterparams, ['noise_sigma'], [3], dry_run=dry_run)
