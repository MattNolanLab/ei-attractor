#!/usr/bin/env python
#
#   figure2.py
#
#   Noise publication Figure 2.
#
#       Copyright (C) 2013  Lukas Solanka <l.solanka@sms.ed.ac.uk>
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
from matplotlib.pyplot   import figure, subplot, plot, savefig
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator, AutoMinorLocator, LinearLocator, MaxNLocator, \
        ScalarFormatter

from parameters  import JobTrialSpace2D
from EI_plotting import plotGridTrial, computeYX, aggregate2DTrial, aggregate2D
from plotting.grids import plotGridRateMap, plotAutoCorrelation, plotSpikes2D
from plotting.global_defs import globalAxesSettings
from figures_shared import plotOneHist

import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)

from matplotlib import rc
rc('pdf', fonttype=42)
rc('mathtext', default='regular')

plt.rcParams['font.size'] = 12

outputDir = "."

NTrials=10
iterList  = ['g_AMPA_total', 'g_GABA_total']

noise_sigmas = [0, 150, 300]
gridsDataRoot= 'output_local/grids'
velDataRoot = 'output_local/velocity'
shape = (30, 30)

grids     = 0
hists     = 1

##############################################################################

def getNoiseRootDir(prefix, noise_sigma):
    return  "{0}/EI_param_sweep_{1}pA".format(prefix, int(noise_sigma))


def getNoiseRoots(prefix, noise_sigmas):
    roots = []
    for s in noise_sigmas:
        roots.append(getNoiseRootDir(prefix, s))
    return roots



def drawGridSweeps(gs, dataSpace, iterList, NTrials=1, r=0, c=0, xLabelOn=True,
        xticks=True, exRows=[], exCols=[], exLetters=[], exDir=[],
        exColor='white'):
    yLabelText = '$w_E$ (nS)'
    if (xLabelOn):
        xLabelText = '$w_I$ (nS)'
    else:
        xLabelText = ''

    ax0 = subplot(gs[0:2, 0])
    G = plotGridTrial(dataSpace, ['gridnessScore'], iterList,
            trialNumList=range(NTrials),
            r=r,
            c=c,
            xlabel=xLabelText,
            ylabel=yLabelText,
            colorBar=False,
            clBarLabel = "Gridness score",
            clbarNTicks=3,
            xticks=xticks,
            vmin=-0.5,
            vmax=1.01)
    print("    max(G): {0}".format(np.max(G)))
    print("    min(G): {0}".format(np.min(G)))





def plotGridnessThresholdComparison(spList, trialNumList, thrList, r=0, c=0,
        ylabelPos=-0.2):
    varList = ['gridnessScore']
    #G = []
    noise_sigma = [0, 150, 300]

    ax = plt.gca()
    plt.hold('on')
    globalAxesSettings(ax)

    for sp in spList:
        G = aggregate2DTrial(sp, varList, trialNumList).flatten()
        counts = []
        for thr in thrList:
            thrIdx = np.logical_and(G > thr, np.logical_not(np.isnan(G)))
            counts.append(float(len(G[thrIdx])) / len(G))
        plt.plot(thrList, counts, 'o-', markersize=4)


    plt.plot([0], [1], linestyle='None', marker='None')
    ax.set_xlabel('Gridness score threshold')
    ax.text(ylabelPos, 0.5, 'Count', rotation=90, transform=ax.transAxes,
            va='center', ha='right')
    leg = []
    for s in noise_sigma:
        leg.append("{0}".format(int(s)))
    l = ax.legend(leg, loc=(0.8, 0.5), title='$\sigma$ (pA)', frameon=False,
            fontsize='x-small', ncol=1)
    plt.setp(l.get_title(), fontsize='x-small')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_major_locator(MultipleLocator(0.3))
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(AutoMinorLocator(3))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.margins(0.02, 0.07)


def plotGridnessHistogram(spList, trialNumList, ylabelPos=-0.2):
    varList = ['gridnessScore']
    #G = []
    noise_sigma = [0, 150, 300]
    colors = ['red', 'green', 'blue']

    ax = plt.gca()
    plt.hold('on')
    globalAxesSettings(ax)

    for idx, sp in enumerate(spList):
        G = aggregate2DTrial(sp, varList, trialNumList).flatten()
        filtIdx = np.logical_not(np.isnan(G))
        plotOneHist(G[filtIdx]) #, color=colors[idx])
    leg = []
    for s in noise_sigma:
        leg.append("{0}".format(int(s)))
    l = ax.legend(leg, loc=(0.8, 0.5), title='$\sigma$ (pA)', frameon=False,
            fontsize='x-small', ncol=1)
    plt.setp(l.get_title(), fontsize='x-small')

    ax.set_xlabel("Gridness score")
    ax.text(ylabelPos, 0.5, 'Count', rotation=90, transform=ax.transAxes,
            va='center', ha='right')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_major_locator(MultipleLocator(50))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(3))
    ax.margins(0.05, 0.025)
    


def plotGridnessVsFitErr(spListGrids, spListVelocity, trialNumList,
        ylabelPos=-0.2, maxErr=None):
    GVars = ['gridnessScore']
    errVars = ['lineFitErr']
    slopeVars = ['lineFitSlope']
    noise_sigma = [0, 150, 300]

    ax = plt.gca()
    plt.hold('on')
    globalAxesSettings(ax)

    for idx, (spGrids, spVel) in enumerate(zip(spListGrids, spListVelocity)):
        G = aggregate2DTrial(spGrids, GVars, trialNumList).flatten()
        errs = aggregate2D(spVel, errVars, funReduce=np.sum).flatten()
        slopes = np.abs(aggregate2D(spVel, slopeVars,
            funReduce=None).flatten())
        #filtIdx = np.logical_not(np.isnan(G))
        ax.plot(G, errs/slopes, 'o', markersize=2)

    if (maxErr is not None):
        ax.set_ylim([0, maxErr])

    leg = []
    for s in noise_sigma:
        leg.append("{0}".format(int(s)))
    l = ax.legend(leg, loc=(0.8, 0.5), title='$\sigma$ (pA)', frameon=False,
            fontsize='x-small', ncol=1)
    plt.setp(l.get_title(), fontsize='x-small')

    ax.set_xlabel("Gridness score")
    ax.text(ylabelPos, 0.5, 'Error of fit (norm., pA)', rotation=90, transform=ax.transAxes,
            va='center', ha='right')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    #ax.yaxis.set_major_locator(MultipleLocator(50))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(3))
    ax.margins(0.05, 0.025)


th = 1 # top plot height
hr = 0.5
space = 0.5
width = 0.75
left = 0.075
right = left + width
top = 0.93
margin = 0.05
div = 0.1
height = 0.22
hspace = 0
wspace = 0.
gsWidthRatios = [th, space, hr, hr, hr]

letter_top=0.95
letter_left_offset=0.06
letter_va='bottom'
letter_ha='left'
letter_ex_yoffset=0.4
letter_ex_xoffset=0.01
letter_ex2_mult = 0.55

gsRows = 2
gsCols = 5

gridRoots = getNoiseRoots(gridsDataRoot, noise_sigmas)
gridDataSpace0   = JobTrialSpace2D(shape, gridRoots[0])
gridDataSpace150 = JobTrialSpace2D(shape, gridRoots[1])
gridDataSpace300 = JobTrialSpace2D(shape, gridRoots[2])

velRoots = getNoiseRoots(velDataRoot, noise_sigmas)
velDataSpace0   = JobTrialSpace2D(shape, velRoots[0])
velDataSpace150 = JobTrialSpace2D(shape, velRoots[1])
velDataSpace300 = JobTrialSpace2D(shape, velRoots[2])

if (grids):
    figSize = (8.5, 6)
    fig = figure(figsize=figSize)

    top = 1. - margin
    bottom = top - height
    # noise_sigm = 0 pA
    gs = GridSpec(gsRows, gsCols, width_ratios=gsWidthRatios)
    gs.update(left=left, right=right, bottom=bottom, top=top, wspace=wspace)
    exRows = [28, 15]
    exCols = [3, 15]
    exDir  = [(1., -1), (1., 1.)]


    # noise_sigma = 150 pA
    gs = GridSpec(gsRows, gsCols, width_ratios=gsWidthRatios)
    top = bottom - div
    bottom = top - height
    gs.update(left=left, right=right, bottom=bottom, top=top, wspace=wspace)
    exRows = [8, 2]
    exCols = [10, 9]
    exDir  = [(-1., 1), (1., 1.)]
    drawGridSweeps(gs, gridDataSpace150, iterList, NTrials=NTrials, r=11, c=10,
            xLabelOn=False, xticks=False, exRows=exRows, exCols=exCols,
            exLetters=['E', 'F'], exDir=exDir)



    # noise_sigma = 300 pA
    gs = GridSpec(gsRows, gsCols, width_ratios=gsWidthRatios)
    top = bottom - div
    bottom = top - height
    gs.update(left=left, right=right, bottom=bottom, top=top, wspace=wspace)
    exRows = [16, 15]
    exCols = [6, 23]
    exDir  = [(1, 1), (1., 1.)]
    drawGridSweeps(gs, gridDataSpace300, iterList, NTrials=NTrials, r=0, c=5,
            xticks=True, exRows=exRows, exCols=exCols, exLetters=['H', 'I'],
            exDir=exDir, exColor='black')

    fname = outputDir + "/figure2_sweeps.png"
    savefig(fname, dpi=300, transparent=True)


# Stats
if (hists):
    ylabelPos = -0.225
    gridSpList = [gridDataSpace0, gridDataSpace150, gridDataSpace300]
    velSpList = [velDataSpace0, velDataSpace150, velDataSpace300]
    fig = figure(figsize=(4, 6))
    gs = GridSpec(3, 1, height_ratios=[0.8, 0.6, 1])

    ax_hist = subplot(gs[0, 0])
    plotGridnessHistogram(gridSpList, range(NTrials), ylabelPos=ylabelPos)

    ax_threshold = subplot(gs[1, 0])
    plotGridnessThresholdComparison(gridSpList, range(NTrials),
            thrList=np.arange(-0.4, 1.2, 0.05), ylabelPos=ylabelPos)

    ax_grids_vel = subplot(gs[2, 0])
    plotGridnessVsFitErr(gridSpList, velSpList, range(NTrials),
            ylabelPos=ylabelPos, maxErr=1000)

    gs.tight_layout(fig, rect=[0.05, 0, 1, 1], h_pad=3.0)
    fname = outputDir + "/figure2_histograms.pdf"
    savefig(fname, dpi=300, transparent=True)


