#!/usr/bin/env python
#
#   figure4.py
#
#   Noise paper figure 4.
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
import matplotlib
#matplotlib.use('cairo')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import LinearLocator

import numpy.ma as ma

from parameters  import JobTrialSpace2D, DataSpace
from EI_plotting import plot2DTrial
from plotting.global_defs import globalAxesSettings

from plot_EI import plotBumpSigmaTrial, plotBumpErrTrial
from figure3 import plot2AxisBar


import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)


# Other
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)
rc('pdf', fonttype=42)

plt.rcParams['font.size'] = 25

###############################################################################
cFreq = 'blue'
cAC = 'green'

###############################################################################

def aggregate2DTrial(sp, varList, trialNumList):
    varList = ['analysis'] + varList
    retVar = sp.aggregateData(varList, trialNumList, funReduce=np.mean,
            saveData=True)
    return np.mean(retVar, 2)



def computeYX(sp, iterList):
    E, I = sp.getIteratedParameters(iterList)
    Ne = DataSpace.getNetParam(sp[0][0][0].data, 'net_Ne')
    Ni = DataSpace.getNetParam(sp[0][0][0].data, 'net_Ni')
    return E/Ne, I/Ni


 
def drawColorbar(drawAx, label):
    pos = drawAx.get_position()
    pos.y0 -= 0.12
    pos.y1 -= 0.12
    pos.y1 = pos.y0 + 0.1*(pos.y1 - pos.y0)
    w = pos.x1 - pos.x0
    pos.x0 += 0.1*w
    pos.x1 -= 0.1*w
    clba = plt.gcf().add_axes(pos)
    globalAxesSettings(clba)
    clba.minorticks_on()
    cb = plt.colorbar(cax=clba, orientation='horizontal',
            ticks=LinearLocator(2))
    cb.set_label(label)


def plot2DNoiseBumps(spList, iterList, trialNumList=[0]):
    NSp = len(spList)
    Sz2D = 1.
    clBarSz = 0.1
    hr = [Sz2D]*NSp + [clBarSz]
    gs = GridSpec(len(spList), 2, wspace=0.15, hspace=0.15, height_ratios=hr,
            bottom=0.22)
    for spIdx in range(NSp):
        if (spIdx == NSp - 1):
            xlabel = "I (nS)"
            xticks = True
            clbar = True
        else:
            xlabel = ""
            xticks = False
            clbar = False


        bumpSigmaThreshold = 20
        ax_sigma = plt.subplot(gs[spIdx, 0])
        bump_sigma = plotBumpSigmaTrial(spList[spIdx], ['bump_e', 'sigma'], iterList,
                thr=bumpSigmaThreshold,
                trialNumList=range(NTrials),
                xlabel=xlabel,
                ylabel='E (nS)',
                colorBar=False,
                vmin=0,
                vmax=bumpSigmaThreshold,
                xticks=xticks,
                clbarNTicks=4)
        if (clbar):
            drawColorbar(ax_sigma, 'Bump $\sigma$ (nrns)')

        ax_err = plt.subplot(gs[spIdx, 1])
        bump_err2 = plotBumpErrTrial(spList[spIdx], ['bump_e', 'err2'], iterList,
                trialNumList=range(NTrials),
                mask=bump_sigma.mask,
                xlabel=xlabel,
                colorBar=False,
                clbarNTicks=3,
                xticks=xticks,
                yticks=False,
                vmin=0,
                vmax=200)
        if (clbar):
            drawColorbar(ax_err, 'Error of fit (Hz)')




###############################################################################
def aggregateBarBumps(spList, varLists, trialNumList, thresholds=(np.infty, \
        np.infty), func=(None, None)):
    # Ugly hack!!!
    # asusming sigma is first, err2 is second
    means = ([], [])
    errs  = ([], [])
    noise_sigma = []
    for idx in xrange(len(spList)):
        mask = False
        for varIdx in range(len(varLists)):
            f = func[varIdx]
            if f is None:
                f = lambda x: x
            var = f(aggregate2DTrial(spList[idx], varLists[varIdx],
                trialNumList).flatten())
            mask = np.logical_or(np.logical_or(np.isnan(var), var >
                thresholds[varIdx]), mask)
            var = ma.MaskedArray(var, mask=mask)
            means[varIdx].append(np.mean(var))
            errs[varIdx].append(np.std(var))
        noise_sigma.append(spList[idx][0][0][0].data['options']['noise_sigma'])

    noise_sigma = np.array(noise_sigma, dtype=int)
    return means, errs, noise_sigma



def plotAggregateBarBumps(spList, trialNumList, ylabels, sigmaTh=20,
        errTh=np.infty):
    varLists = [['bump_e', 'sigma'], ['bump_e', 'err2']]
    (sigmaMean, err2Mean), (sigmaStd, err2Std), noise_sigma = \
            aggregateBarBumps(spList, varLists, trialNumList, thresholds=(sigmaTh,
                errTh), func=(None, np.sqrt))
    

    print sigmaMean, err2Mean

    plot2AxisBar(noise_sigma,
            (sigmaMean, err2Mean),
            (sigmaStd,  err2Std),
            xlabel='$\sigma$ (pA)',
            ylabels=ylabels,
            colors=(cAC, cFreq))
###############################################################################
            
if (__name__ == '__main__'):

    NTrials = 5
    iterList  = ['g_AMPA_total', 'g_GABA_total']
    r = 14
    c = 13

    baseDir = 'output_local/one_to_one'
    dirs = [ \
        ('EI_param_sweep_0pA',    (40, 40)),
        ('EI_param_sweep_150pA',  (40, 40)),
        ('EI_param_sweep_300pA',  (40, 40))
    ]

    dataSpaces = []
    for (dir, shape) in dirs:
        rootDir = "{0}/{1}".format(baseDir, dir)
        dataSpaces.append(JobTrialSpace2D(shape, rootDir))
        

    ###############################################################################
    plt.figure(figsize=(6.5, 8))
    plot2DNoiseBumps(dataSpaces, iterList)
    #plt.tight_layout()
    plt.savefig('{0}/noise_bumps_2d.png'.format(baseDir), transparent=True,
            dpi=300)
    ###############################################################################
    plt.figure(figsize=(5.5, 5))
    plotAggregateBarBumps(dataSpaces,
            trialNumList= range(NTrials),
            ylabels=('Bump $\sigma$ (neurons)', 'Error of fit (Hz)'))
    plt.tight_layout()
    plt.savefig('{0}/aggregateBar_bump.pdf'.format(baseDir), transparent=True,
            dpi=300)
    ################################################################################
