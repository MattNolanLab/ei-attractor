#!/usr/bin/env python
#
#   param_sweep_noise.py
#
#   Detailed noise level parameter sweeps.
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


def submitNoiseSweep(p, gEp, gIp, noisep,
        ENV, simRootDir, simLabel, appName, rtLimit, numCPU, blocking,
        timePrefix, numRepeat, dry_run, extraIterparams={}):
    '''
    Submit a parameter sweep with an extra dimension: noise_sigma.
    '''
    noise_sigma_arr  = []
    g_AMPA_total_arr = []
    g_GABA_total_arr = []
    for n in noisep.values():
        for E_coupling in gEp.values():
            for I_coupling in gIp.values():
                noise_sigma_arr.append(n)
                g_AMPA_total_arr.append(E_coupling)
                g_GABA_total_arr.append(I_coupling)


    iterparams = {
        'noise_sigma'  : np.array(noise_sigma_arr),
        'g_AMPA_total' : np.array(g_AMPA_total_arr),
        'g_GABA_total' : np.array(g_GABA_total_arr),
        #'noise_sigma'  : [150],
        #'g_AMPA_total' : [1020],
        #'g_GABA_total' : [3060]
    }
    iterparams.update(extraIterparams)
    ac = ArgumentCreator(p, printout=True)
    ac.insertDict(iterparams, mult=False)

    ###############################################################################
    submitter = SubmitterFactory.getSubmitter(ac, appName, envType=ENV,
            rtLimit=rtLimit, output_dir=simRootDir, label=simLabel,
            blocking=blocking, timePrefix=timePrefix, numCPU=numCPU)
    ac.setOption('output_dir', submitter.outputDir())
    startJobNum = 0
    submitter.submitAll(startJobNum, numRepeat, dry_run=dry_run)
    submitter.saveIterParams(iterparams, dry_run=dry_run)


class SweepParams(object):
    '''
    Parameters of a linear sweep in one dimension.
    '''
    def __init__(self, start, end, nVals):
        self.start = start
        self.end   = end
        self.nVals = nVals

    def values(self):
        return np.linspace(self.start, self.end, self.nVals)

