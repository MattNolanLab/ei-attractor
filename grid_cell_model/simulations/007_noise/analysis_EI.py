#!/usr/bin/env python
#
#   analysis_EI.py
#
#   Theta/gamma analysis using a custom "peak" method - E/I coupling parameter
#   sweep.
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
import matplotlib.pyplot as plt

import numpy.ma as ma
from matplotlib.ticker  import MaxNLocator, LinearLocator

from analysis.visitors    import AutoCorrelationVisitor
from parameters           import JobTrialSpace2D
from plotting.global_defs import globalAxesSettings, createColorbar
import logging as lg
lg.basicConfig(level=lg.WARN)
#lg.basicConfig(level=lg.INFO)


# Other
plt.rcParams['font.size'] = 11

###############################################################################

def plot2DTrial(sp, varName, iterList, trialNum=0, xlabel="", ylabel="",
        colorBar=True, clBarLabel="", title="", clbarNTicks=2, xticks=True,
        yticks=True):
    if (len(iterList) != 2):
        raise ValueError("iterList must contain exactly 2 elements.")
    shape = sp.getShape()
    rows = shape[0]
    cols = shape[1]
    retVar = np.ndarray(shape)
    for r in xrange(rows):
        for c in xrange(cols):
            if (len(sp[r][c]) == 0):
                retVar[r][c] = np.nan
            else:
                data = sp[r][c][trialNum].data['analysis']
                retVar[r][c] = np.mean(data[varName])

    retVar =  ma.MaskedArray(retVar, mask=np.isnan(retVar))


    ax = plt.gca()
    globalAxesSettings(ax)
    Y, X = sp.getIteratedParameters(iterList)
    plt.pcolormesh(X, Y, retVar)
    createColorbar(ax, retVar, clBarLabel, nticks=clbarNTicks,
            orientation='horizontal', pad=0.3)
    if (xlabel != ""):
        plt.xlabel(xlabel)
    if (ylabel != ""):
        #plt.ylabel(ylabel, color='white')
        plt.ylabel(ylabel)
    plt.axis('tight')
    ax.xaxis.set_major_locator(LinearLocator(3))
    ax.yaxis.set_major_locator(LinearLocator(3))
    if (not xticks):
        ax.xaxis.set_ticklabels([])
    if (not yticks):
        ax.yaxis.set_ticklabels([])

    return retVar


###############################################################################

dirs = {
    #'./output_local/2013-04-24T15-27-30_EI_param_sweep_0pA_big'   : (40, 40),
    #'./output_local/2013-04-24T21-37-47_EI_param_sweep_150pA_big' : (40, 40),
    'output_local/2013-04-24T21-43-32_EI_param_sweep_300pA_big' : (40, 40)
}
NTrials = 5

for rootDir, shape in dirs.iteritems():
    sp = JobTrialSpace2D(shape, rootDir)
    
    monName   = 'stateMonF_e'
    stateList = ['I_clamp_GABA_A']
    iterList  = ['g_AMPA_total', 'g_GABA_total']
    forceUpdate = False
    visitor = AutoCorrelationVisitor(monName, stateList, forceUpdate=forceUpdate)
    sp.visit(visitor)

    for trialNum in xrange(NTrials):
        print("Plotting results for trial {0}".format(trialNum))
        
        ###############################################################################
        plt.figure(figsize=(5.1, 2.9))
        N = 2
        plt.subplot(1, N, 1)
        
        acVal = plot2DTrial(sp, 'acVal', iterList,
                trialNum=trialNum,
                xlabel="I coupling strength (nS)",
                ylabel='E coupling strength (nS)',
                clBarLabel = "Correlation",
                clbarNTicks=3)
        
        ###############################################################################
        plt.subplot(1, N, 2)
        freq = plot2DTrial(sp, 'freq', iterList,
                trialNum=trialNum,
                xlabel="I coupling strength (nS)",
                clBarLabel = "Frequency (Hz)",
                clbarNTicks=3,
                yticks=False)
        ###############################################################################
        plt.tight_layout()
        noise_sigma = sp.getParam(sp[0][0][0].data, 'noise_sigma')
        plt.savefig(sp.rootDir +
                '/analysis_EI_{0}pA_trial{1}.png'.format(int(noise_sigma),
                    trialNum))

#plt.show()

