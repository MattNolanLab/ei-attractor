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
from matplotlib.ticker import LinearLocator, NullLocator

import numpy.ma as ma

from parameters      import JobTrialSpace2D, DataSpace
from analysis.visitors.interface import EmptyDSVisitor
from analysis.spikes import PopulationSpikes
from EI_plotting     import plot2DTrial
from plotting.global_defs import globalAxesSettings
from plotting.bumps  import torusFiringRate
from plot_EI         import plotBumpSigmaTrial, plotBumpErrTrial
from figures_shared  import _getSpikeTrain
from figure3         import plot2AxisBar


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


def drawBumpColorbar(drawAx):
    pos = drawAx.get_position()
    pos.x0 += 0.05
    pos.x1 += 0.05
    pos.x0 = pos.x1 - 0.1*(pos.x1 - pos.x0)
    h = pos.y1 - pos.y0
    pos.y0 += 0.1*h
    pos.y1 -= 0.1*h
    clba = plt.gcf().add_axes(pos)
    globalAxesSettings(clba)
    cb = plt.colorbar(cax=clba, orientation='vertical',
            ticks=NullLocator())
    #clba.yaxis.set_ticklabels(['min'], ['max'])
    cb.set_label('F (Hz)', fontsize='small')

def plotBumps(spList, r, c, trialNum=0):
    N = len(spList)
    gs = GridSpec(1, len(spList), wspace=0.15, left=0.01, right=0.85,
            bottom=0.01, top=0.99)
    for spIdx in range(len(spList)):
        ax = plt.subplot(gs[0, spIdx])
        data = spList[spIdx][r][c][trialNum].data
        v = EmptyDSVisitor()
        tStart = v.getOption(data, 'theta_start_t')
        tEnd   = v.getOption(data, 'time')
        ESenders, ETimes, EN = v._getSpikeTrain(data, 'spikeMon_e', ('Ne_x',
            'Ne_y'))
        s = PopulationSpikes(EN, ESenders, ETimes)
        Fe = s.avgFiringRate(tStart, tEnd)
        Ne_x = v.getNetParam(data, 'Ne_x')
        Ne_y = v.getNetParam(data, 'Ne_y')
        bump_e = np.reshape(Fe, (Ne_y, Ne_x))

        if (spIdx == 0):
            yTickLabels = True
        else:
            yTickLabels = False
        torusFiringRate(
                rateMap  = bump_e,
                labelx   = '',
                labely   = '',
                titleStr = '',
                clbar    = False,
                xTicks   = False,
                yTicks   = False,
                xTickLabels = False,
                yTickLabels = False)

    drawBumpColorbar(ax)


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
    plt.figure(figsize=(5, 1.5))
    plotBumps(dataSpaces, r, c)
    plt.savefig('{0}/bumps_example.png'.format(baseDir), transparent=False,
            dpi=300)
    #################################################################################
    rc('text', usetex=True)
    fig = plt.figure(figsize=(5.5, 1.5))
    fig.text(0.5, 0.5,
            '$G(x, y) = A \exp\left(\\frac{(x - \mu_x)^2 + (y - ' +
            ' \mu_y)^2}{2\sigma^2}\\right)$',
            ha='center', va='center', fontsize=plt.rcParams['font.size'] + 3)
    plt.savefig('{0}/gaussian_equation.pdf'.format(baseDir), transparent=True)

