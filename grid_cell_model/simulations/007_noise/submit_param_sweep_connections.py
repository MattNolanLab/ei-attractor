#!/usr/bin/env python
'''Submit job(s) to the cluster/workstation: export connections to file'''
from __future__ import absolute_import, print_function

import numpy as np
from grid_cell_model.submitting.noise import SubmissionParser
from grid_cell_model.submitting.factory import SubmitterFactory
from grid_cell_model.submitting.arguments import ArgumentCreator

from default_params import defaultParameters as dp

parser = SubmissionParser()
parser.add_argument('--probabilistic_synapses', type=int, choices=(0, 1),
                    default=0, help=('Whether the network should generate '
                                     'synapse weights probabilistically.'))
o = parser.parse_args()

p = dp.copy()

# Submitting
ENV         = o.env
simRootDir  = o.where
simLabel    = 'connections'
appName     = '../common/simulation_connections.py'
rtLimit     = o.rtLimit
numCPU      = o.nCPU
blocking    = True
timePrefix  = False
numRepeat   = 1
dry_run     = o.dry_run

p['master_seed'] = 123456
p['time']        = 0.1 # unused
p['nthreads']    = 1
p['ntrials']     = 1
p['verbosity']        = o.verbosity

p['probabilistic_synapses'] = o.probabilistic_synapses

# Range of E/I synaptic conductances
Nvals  = 31       # Number of values (these will be equivalent to go through the
                 # diagonal of the parameter space
startG = 0.0     # nS
endG   = 6120.0  # nS

###############################################################################
ac = ArgumentCreator(p, printout=True)

GArr = np.linspace(startG, endG, Nvals)
g_AMPA_total_arr = []
g_GABA_total_arr = []
for coupling in GArr:
    g_AMPA_total_arr.append(coupling)
    g_GABA_total_arr.append(coupling)

iterparams = {
    'g_AMPA_total'      : np.array(g_AMPA_total_arr),
    'g_GABA_total'      : np.array(g_GABA_total_arr),
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
submitter.saveIterParams(iterparams, ['g_total'], [1], dry_run=dry_run)
