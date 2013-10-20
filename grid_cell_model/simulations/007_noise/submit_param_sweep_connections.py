#!/usr/bin/env python
#
#   submit_param_sweep_connections.py
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
import numpy as np
from submitting.factory   import SubmitterFactory
from submitting.arguments import ArgumentCreator
from default_params       import defaultParameters as dp
import logging as lg
#lg.basicConfig(level=lg.DEBUG)
lg.basicConfig(level=lg.INFO)


p = dp.copy()

# Submitting
ENV         = 'workstation'
simRootDir  = 'output_local/even_spacing'
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
submitter.submitAll(startJobNum, numRepeat, dry_run=dry_run)
submitter.saveIterParams(iterparams, dry_run=dry_run)
