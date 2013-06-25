#!/usr/bin/env python
#
#   submit_param_sweep.py
#
#   Submit job(s) to the cluster/workstation: parameter sweep (noise)
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
from default_params         import defaultParameters as p
from submitting.factory     import SubmitterFactory
from submitting.arguments   import ArgumentCreator
import logging as lg
lg.basicConfig(level=lg.DEBUG)

p['noise_sigma']       = 150.0     # pA

# Submitting
ENV         = 'workstation'
simRootDir  = 'output_local/velocity'
simLabel    = 'EI_param_sweep_{0}pA_big'.format(int(p['noise_sigma']))
appName     = 'simulation_velocity.py'
rtLimit     = '00:20:00'
numCPU      = 4
blocking    = True
timePrefix  = False
numRepeat   = 1
dry_run     = False

p['time']       = 10e3  # ms
p['nthreads']   = 4
p['ntrials']    = 10

p['IvelMax']    = 150
p['dIvel']      = 10



###############################################################################

ac = ArgumentCreator(p, printout=True)

# Range of parameters around default values
# Let's choose a 0.5 - 2 range around the default values
Nvals        = 40    # Number of values for each dimension
startFrac    = 0.
endFrac      = 2.8572


#fracArr = np.linspace(startFrac, endFrac, Nvals)
fracArr = [1.0]
print(fracArr)

g_AMPA_total_arr     = []
g_GABA_total_arr     = []
g_uni_GABA_total_arr = []
for E_coupling in fracArr:
    for I_coupling in fracArr:
        g_AMPA_total_arr.append(E_coupling*p['g_AMPA_total'])
        g_GABA_total_arr.append(I_coupling*p['g_GABA_total'])
        g_uni_GABA_total_arr.append(I_coupling*p['g_uni_GABA_total'])


iterparams = {
    'g_AMPA_total'      : np.array(g_AMPA_total_arr),
    'g_GABA_total'      : np.array(g_GABA_total_arr),
    'g_uni_GABA_total'  : np.array(g_uni_GABA_total_arr)
}
ac.insertDict(iterparams, mult=False)

###############################################################################
submitter = SubmitterFactory.getSubmitter(ac, appName, envType=ENV,
        rtLimit=rtLimit, output_dir=simRootDir, label=simLabel,
        blocking=blocking, timePrefix=timePrefix, numCPU=numCPU)
ac.setOption('output_dir', submitter.outputDir())
startJobNum = 0
submitter.submitAll(startJobNum, numRepeat, dry_run=dry_run)
submitter.saveIterParams(iterparams)
