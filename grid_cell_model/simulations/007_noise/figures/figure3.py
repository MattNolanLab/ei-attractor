#!/usr/bin/env python
#
#   figure3.py
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
from matplotlib.pyplot   import figure, plot, savefig
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator, AutoMinorLocator, LinearLocator, MaxNLocator, \
        ScalarFormatter
from matplotlib.transforms import Bbox
from matplotlib.colorbar import make_axes

import numpy.ma as ma

from parameters  import JobTrialSpace2D, DataSpace
from EI_plotting import plot2DTrial, plotACTrial
from plotting.global_defs import globalAxesSettings
from figures_shared import getNoiseRoots, plotOneHist
import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)


# Other
from matplotlib import rc
rc('pdf', fonttype=42)
rc('mathtext', default='regular')

plt.rcParams['font.size'] = 11

###############################################################################
cFreq = 'blue'
cAC = 'green'
cCount = 'red'

outputDir = "."
NTrials = 5
iterList  = ['g_AMPA_total', 'g_GABA_total']

noise_sigmas  = [0, 150, 300]
exampleIdx    = [(0, 0), (0, 0), (0, 0)] # (row, col)
gammaDataRoot = 'output_local/even_spacing/gamma_bump'
gammaShape    = (31, 31)

gammaSweep0   = 1
gammaSweep150 = 1
gammaSweep300 = 1
threshold     = 1
freqHist      = 1

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


def drawGammaPowerSweeps(ax, dataSpace, iterList, noise_sigma, trialNumList,
        r=0, c=0, yLabelOn=True, yticks=True, exRows=[], exCols=[],
        exColor='white', cbar=False):
    xLabelText = '$g_I$ (nS)'
    if (yLabelOn):
        yLabelText = '$g_E$ (nS)'
    else:
        yLabelText = ''

    if (ax is None):
        ax = plt.gca()

    C = plotACTrial(dataSpace, ['acVal'], iterList,
            trialNumList=trialNumList,
            xlabel=xLabelText,
            ylabel=yLabelText,
            colorBar=False,
            yticks=yticks,
            vmin=-0.09,
            vmax=0.675)
    cax, kw = make_axes(ax, orientation='vertical', shrink=0.9,
            pad=0.05)
    globalAxesSettings(cax)
    cb = plt.colorbar(ax=ax, cax=cax, ticks=MultipleLocator(0.2), **kw)
    cax.yaxis.set_minor_locator(AutoMinorLocator(2))
    cb.set_label('Correlation')
    if (cbar == False):
        cax.set_visible(False)
    ax.set_title('$\sigma$ = {0} pA'.format(int(noise_sigma)))

    print("    max(C): {0}".format(np.max(C)))
    print("    min(C): {0}".format(np.min(C)))
    return ax, cax




def extractACExample(sp, r, c, trialNum):
    data = sp[r][c][trialNum].data
    ac = data['analysis']['acVec'][0]
    dt = data['stateMonF_e'][0]['interval']
    freq = data['analysis']['freq'][0]
    acVal = data['analysis']['acVal'][0]
    noise_sigma = data['options']['noise_sigma']
    return ac, dt, freq, acVal, noise_sigma


def aggregateBar2(spList, varLists, trialNumList, func=(None, None)):
    vars = ([], [])
    noise_sigma = []
    for idx in xrange(len(spList)):
        for varIdx in range(len(varLists)):
            f = func[varIdx]
            if f is None:
                f = lambda x: x
            vars[varIdx].append(f(aggregate2DTrial(spList[idx], varLists[varIdx],
                trialNumList).flatten()))
        noise_sigma.append(spList[idx][0][0][0].data['options']['noise_sigma'])

    noise_sigma = np.array(noise_sigma, dtype=int)
    return vars, noise_sigma
 


def getACFreqThreshold(spList, trialNumList, ACThr):
    varLists = [['acVal'], ['freq']]
    vars, noise_sigma = aggregateBar2(spList, varLists, trialNumList)
    AC = vars[0]
    freq = vars[1]
    ACMean   = []
    freqMean = []
    ACStd    = []
    freqStd  = []
    thrCount = []
    for spIdx in range(len(spList)):
        thrIdx = np.logical_and(AC[spIdx] >= ACThr,
                np.logical_not(np.isnan(AC[spIdx])))
        ac_filt = AC[spIdx][thrIdx]
        ACMean.append(np.mean(ac_filt))
        ACStd.append(np.std(ac_filt))
        freq_filt = freq[spIdx][thrIdx]
        freqMean.append(np.mean(freq_filt))
        freqStd.append(np.std(freq_filt))
        thrCount.append(float(len(AC[spIdx][thrIdx])) / len (AC[spIdx]))

    return (ACMean, ACStd), (freqMean, freqStd), thrCount, noise_sigma


def plotThresholdComparison(spList, trialNumList, ACThrList):
    counts = []
    noise_sigma = None
    for ACThr in ACThrList:
        _, _, thrCount, noise_sigma = getACFreqThreshold(spList, trialNumList,
                ACThr)
        counts.append(thrCount)
    counts = np.array(counts)

    print ACThrList, counts

    ax = plt.gca()
    globalAxesSettings(ax)
    plt.plot(ACThrList, counts, 'o-')
    plt.plot([0], [1], linestyle='None', marker='None')
    ax.set_xlabel('Correlation threshold', labelpad=5)
    ax.set_ylabel('Count', labelpad=5)
    leg = []
    for s in noise_sigma:
        leg.append("{0}".format(int(s)))
    ax.legend(leg, loc=(0.8, 0.55), title='$\sigma$ (pA)', frameon=False,
            fontsize='small')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_major_locator(MultipleLocator(0.3))
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(AutoMinorLocator(3))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.margins(0.025)
    

def plotFreqHistogram(spList, trialNumList, ylabelPos=-0.2, CThreshold=0.1):
    FVarList = ['freq']
    CVarList = ['acVal']
    noise_sigma = [0, 150, 300]
    colors = ['red', 'green', 'blue']

    ax = plt.gca()
    plt.hold('on')
    globalAxesSettings(ax)

    for idx, sp in enumerate(spList):
        F = aggregate2DTrial(sp, FVarList, trialNumList).flatten()
        C = aggregate2DTrial(sp, CVarList, trialNumList).flatten()
        filtIdx = np.logical_and(np.logical_not(np.isnan(F)), C > CThreshold)
        plotOneHist(F[filtIdx], bins=20, normed=True)
    leg = []
    for s in noise_sigma:
        leg.append("{0}".format(int(s)))
    l = ax.legend(leg, loc=(0.8, 0.5), title='$\sigma$ (pA)', frameon=False,
            fontsize='x-small', ncol=1)
    plt.setp(l.get_title(), fontsize='x-small')

    ax.set_xlabel("Oscillation frequency (Hz)")
    #ax.text(ylabelPos, 0.5, 'p(F)', rotation=90, transform=ax.transAxes,
    #        va='center', ha='right')
    ax.set_ylabel('p(Frequency)')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_major_locator(MultipleLocator(20))
    ax.yaxis.set_major_locator(MaxNLocator(4))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    f = ScalarFormatter(useMathText=True)
    f.set_scientific(True)
    f.set_powerlimits([0, 3])
    ax.yaxis.set_major_formatter(f)
    #ax.margins(0.01, 0.00)

    thStr = 'Frequencies with C > {0}'.format(CThreshold)
    ax.text(0.99, 1.1, thStr, transform=ax.transAxes, va='bottom',
            ha='right')
    



###############################################################################

gammaRoots = getNoiseRoots(gammaDataRoot, noise_sigmas)
gammaDataSpace0   = JobTrialSpace2D(gammaShape, gammaRoots[0])
gammaDataSpace150 = JobTrialSpace2D(gammaShape, gammaRoots[1])
gammaDataSpace300 = JobTrialSpace2D(gammaShape, gammaRoots[2])


sweepFigSize = (3.2, 2.5)
sweepLeft   = 0.15
sweepBottom = 0.2
sweepRight  = 0.85
sweepTop    = 0.85


if (gammaSweep0):
    # noise_sigma = 0 pA
    fig = figure("sweeps0", figsize=sweepFigSize)
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    ax, cax = drawGammaPowerSweeps(ax, gammaDataSpace0, iterList,
            noise_sigma=noise_sigmas[0], trialNumList=xrange(NTrials), cbar=False)
    fname = outputDir + "/figure3_sweeps0.png"
    fig.savefig(fname, dpi=300, transparent=True)
        
if (gammaSweep150):
    # noise_sigma = 0 pA
    fig = figure("sweeps150", figsize=sweepFigSize)
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    ax, cax = drawGammaPowerSweeps(ax, gammaDataSpace150, iterList,
            noise_sigma=noise_sigmas[1], trialNumList=xrange(NTrials),
            yticks=False, yLabelOn=False, cbar=False)
    fname = outputDir + "/figure3_sweeps150.png"
    fig.savefig(fname, dpi=300, transparent=True)
        
if (gammaSweep150):
    # noise_sigma = 0 pA
    fig = figure("sweeps300", figsize=sweepFigSize)
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    ax, cax = drawGammaPowerSweeps(ax, gammaDataSpace300, iterList,
            noise_sigma=noise_sigmas[2], trialNumList=xrange(NTrials),
            yticks=False, yLabelOn=False, cbar=True)
    fname = outputDir + "/figure3_sweeps300.png"
    fig.savefig(fname, dpi=300, transparent=True)
        

gammaSpList = [gammaDataSpace0, gammaDataSpace150, gammaDataSpace300]
if (threshold):
    ###############################################################################
    plt.figure(figsize=(3.5, 2))
    plotThresholdComparison(gammaSpList,
            trialNumList=range(NTrials),
            ACThrList=np.arange(0, 0.65, 0.05))
    plt.tight_layout()
    plt.savefig('figure3_AC_threshold_comparison.pdf'.format(baseDir), transparent=True,
            dpi=300)


if (freqHist):
    ylabelPos = -0.16
    fig = figure(figsize=(3.7, 2.5))
    plotFreqHistogram(gammaSpList, range(NTrials), ylabelPos=ylabelPos)
    plt.tight_layout()
    fname = outputDir + "/figure3_freq_histograms.pdf"
    savefig(fname, dpi=300, transparent=True)


