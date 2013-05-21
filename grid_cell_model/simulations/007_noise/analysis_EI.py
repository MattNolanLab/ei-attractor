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

from analysis.visitors    import AutoCorrelationVisitor, BumpFittingVisitor, \
        FiringRateVisitor
from parameters           import JobTrialSpace2D
from plotting.global_defs import globalAxesSettings, createColorbar
import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)


# Other
plt.rcParams['font.size'] = 11

###############################################################################

def getDictData(d, keyList):
    ret = d
    for key in keyList:
        ret = ret[key]
    return ret

def aggregate2DTrial(sp, varList, iterList, trialNum):
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
                retVar[r][c] = np.mean(getDictData(data, varList))

    return retVar

def plot2DTrial(X, Y, C, xlabel="", ylabel="",
        colorBar=True, clBarLabel="", vmax=None, title="", clbarNTicks=2,
        xticks=True, yticks=True):

    ax = plt.gca()
    globalAxesSettings(ax)
    plt.pcolormesh(X, Y, C, vmax=vmax)
    if (clbarNTicks == None):
        createColorbar(ax, None, clBarLabel, orientation='horizontal', pad=0.2)
    else:
        createColorbar(ax, C, clBarLabel, nticks=clbarNTicks,
                orientation='horizontal', pad=0.2)
    if (xlabel != ""):
        plt.xlabel(xlabel, va='top')
        ax.xaxis.set_label_coords(0.5, -0.10)
    if (ylabel != ""):
        plt.ylabel(ylabel, ha='right')
        ax.yaxis.set_label_coords(-0.10, 0.5)
    plt.axis('tight')
    ax.xaxis.set_major_locator(LinearLocator(2))
    ax.yaxis.set_major_locator(LinearLocator(2))
    if (not xticks):
        ax.xaxis.set_ticklabels([])
    if (not yticks):
        ax.yaxis.set_ticklabels([])

    return C


def plotACTrial(sp, varList, iterList, trialNum=0, xlabel="", ylabel="",
        colorBar=True, clBarLabel="", vmax=None, title="", clbarNTicks=2,
        xticks=True, yticks=True):
    C = aggregate2DTrial(sp, varList, iterList, trialNum)
    C = ma.MaskedArray(C, mask=np.isnan(C))
    Y, X = sp.getIteratedParameters(iterList)
    plot2DTrial(X, Y, C, xlabel, ylabel, colorBar, clBarLabel, vmax, title,
            clbarNTicks, xticks, yticks)

def plotBumpSigmaTrial(sp, varList, iterList, thr=np.infty, trialNum=0, xlabel="", ylabel="",
        colorBar=True, clBarLabel="", vmax=None, title="", clbarNTicks=2,
        xticks=True, yticks=True):
    C = aggregate2DTrial(sp, varList, iterList, trialNum)
    C = ma.MaskedArray(C, mask=np.logical_or(np.isnan(C), C > thr))
    Y, X = sp.getIteratedParameters(iterList)
    return plot2DTrial(X, Y, C, xlabel, ylabel, colorBar, clBarLabel, vmax,
            title, clbarNTicks, xticks, yticks)

def plotBumpErrTrial(sp, varList, iterList, thr=np.infty, mask=None,
        trialNum=0, xlabel="", ylabel="", colorBar=True, clBarLabel="",
        vmax=None, title="", clbarNTicks=2, xticks=True, yticks=True):
    C = np.sqrt(aggregate2DTrial(sp, varList, iterList, trialNum))
    if mask is None:
        mask = False
    C = ma.MaskedArray(C, mask=np.logical_or(np.logical_or(np.isnan(C), C >
        thr), mask))
    Y, X = sp.getIteratedParameters(iterList)
    return plot2DTrial(X, Y, C, xlabel, ylabel, colorBar, clBarLabel, vmax,
            title, clbarNTicks, xticks, yticks)


###############################################################################

dirs = {
    './output_local/2013-03-30T19-29-21_EI_param_sweep_0pA_small_sample' : (2, 2),
    #'./output_local/2013-04-24T15-27-30_EI_param_sweep_0pA_big'   : (40, 40),
    #'./output_local/2013-04-24T21-37-47_EI_param_sweep_150pA_big' : (40, 40),
    #'output_local/2013-04-24T21-43-32_EI_param_sweep_300pA_big' : (40, 40)
}
NTrials = 5

for rootDir, shape in dirs.iteritems():
    sp = JobTrialSpace2D(shape, rootDir)
    
    monName   = 'stateMonF_e'
    stateList = ['I_clamp_GABA_A']
    iterList  = ['g_AMPA_total', 'g_GABA_total']
    forceUpdate = False
    ACVisitor = AutoCorrelationVisitor(monName, stateList, forceUpdate=forceUpdate)
    bumpVisitor = BumpFittingVisitor(forceUpdate=forceUpdate)
    FRVisitor = FiringRateVisitor(forceUpdate=forceUpdate)

    sp.visit(ACVisitor)
    sp.visit(bumpVisitor)
    sp.visit(FRVisitor)

    for trialNum in xrange(NTrials):
        print("Plotting results for trial {0}".format(trialNum))
        
        ###############################################################################
        plt.figure(figsize=(5.1, 2.9))
        N = 2
        plt.subplot(1, N, 1)
        
        acVal = plotACTrial(sp, ['acVal'], iterList,
                trialNum=trialNum,
                xlabel="I (nS)",
                ylabel='E (nS)',
                clBarLabel = "Correlation",
                clbarNTicks=3)
        ###############################################################################
        plt.subplot(1, N, 2)
        freq = plotACTrial(sp, ['freq'], iterList,
                trialNum=trialNum,
                xlabel="I (nS)",
                clBarLabel = "Frequency (Hz)",
                clbarNTicks=3,
                yticks=False)
        ###############################################################################
        plt.tight_layout()
        noise_sigma = sp.getParam(sp[0][0][0].data, 'noise_sigma')
        plt.savefig(sp.rootDir +
                '/analysis_EI_{0}pA_trial{1}.png'.format(int(noise_sigma),
                    trialNum))

        ###############################################################################
        ###############################################################################
        plt.figure(figsize=(5.1, 2.9))
        N = 2
        plt.subplot(1, N, 1)
        
        bumpSigmaThreshold = 50
        bump_sigma = plotBumpSigmaTrial(sp, ['bump_e', 'sigma'], iterList,
                thr=bumpSigmaThreshold,
                trialNum=trialNum,
                xlabel="I (nS)",
                ylabel='E (nS)',
                clBarLabel = "Gaussian $\sigma$ (nrns)",
                #vmax  = 100,
                clbarNTicks=4)
        ###############################################################################
        plt.subplot(1, N, 2)
        bump_err2 = plotBumpErrTrial(sp, ['bump_e', 'err2'], iterList,
                trialNum=trialNum,
                mask=bump_sigma.mask,
                xlabel="I (nS)",
                clBarLabel = "Error of fit (Hz)",
                clbarNTicks=3,
                yticks=False)
        ###############################################################################
        plt.tight_layout()
        noise_sigma = sp.getParam(sp[0][0][0].data, 'noise_sigma')
        plt.savefig(sp.rootDir +
                '/analysis_EI_{0}pA_bump_trial{1}.png'.format(int(noise_sigma),
                    trialNum))

#plt.show()

