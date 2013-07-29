#
#   param_sweep.py
#
#   Parameter sweep (2D): shared procedures.
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
from submitting.arguments   import ArgumentCreator
from data_storage           import DataStorage


def submitParamSweep(p, startFrac, endFrac, Nvals, ENV, simRootDir, simLabel,
        appName, rtLimit, numCPU, blocking, timePrefix, numRepeat, dry_run,
        extraIterparams={}):
    ac = ArgumentCreator(p, printout=True)

    fracArr = np.linspace(startFrac, endFrac, Nvals)
    #fracArr = [1.0]
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
    iterparams.update(extraIterparams)
    ac.insertDict(iterparams, mult=False)

    ###############################################################################
    submitter = SubmitterFactory.getSubmitter(ac, appName, envType=ENV,
            rtLimit=rtLimit, output_dir=simRootDir, label=simLabel,
            blocking=blocking, timePrefix=timePrefix, numCPU=numCPU)
    ac.setOption('output_dir', submitter.outputDir())
    startJobNum = 0
    submitter.submitAll(startJobNum, numRepeat, dry_run=dry_run)
    submitter.saveIterParams(iterparams, dry_run=dry_run)


###############################################################################

def getBumpCurrentSlope(noise_sigma, threshold=0):
    fileName = 'bump_slope_data/bump_slope_{0}pA.h5'.format(int(noise_sigma))
    ds = DataStorage.open(fileName, 'r')
    slopes = np.abs(ds['lineFitSlope'].flatten())
    ds.close()
    slopes[slopes < threshold] = np.nan
    return slopes

