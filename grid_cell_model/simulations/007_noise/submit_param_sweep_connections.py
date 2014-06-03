#!/usr/bin/env python
'''Submit job(s) to the cluster/workstation: export connections to file'''
from __future__ import absolute_import, print_function

import logging as lg

import numpy as np
from grid_cell_model.submitting.factory   import SubmitterFactory
from grid_cell_model.submitting.arguments import ArgumentCreator

from default_params       import defaultParameters as dp

#lg.basicConfig(level=lg.DEBUG)
lg.basicConfig(level=lg.INFO)


p = dp.copy()

# Submitting
ENV         = 'workstation'
simRootDir  = 'output/submission'
simLabel    = 'connections'
appName     = 'simulation_connections.py'
rtLimit     = '00:02:00'
numCPU      = 1
blocking    = True
timePrefix  = False
numRepeat   = 1
dry_run     = False

p['time']     = 0.1 # unused
p['nthreads'] = 1
p['ntrials']  = 1


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
    #'g_AMPA_total'      : [1400],
    #'g_GABA_total'      : [2160]
}
ac.insertDict(iterparams, mult=False)

###############################################################################
submitter = SubmitterFactory.getSubmitter(ac, appName, envType=ENV,
        rtLimit=rtLimit, output_dir=simRootDir, label=simLabel,
        blocking=blocking, timePrefix=timePrefix, numCPU=numCPU)
ac.setOption('output_dir', submitter.outputDir())
startJobNum = 0
#submitter.submitAll(startJobNum, numRepeat, dry_run=dry_run)
submitter.saveIterParams(iterparams, dry_run=dry_run)
