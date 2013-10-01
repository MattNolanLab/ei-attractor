#!/usr/bin/env python
#
#   figure1.py
#
#   Noise publication Figure 1.
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
from matplotlib.gridspec   import GridSpec
from matplotlib.ticker     import MultipleLocator, AutoMinorLocator
from matplotlib.colorbar   import make_axes
from matplotlib.transforms import Bbox

from EI_plotting          import plotGridTrial, aggregate2DTrial, plotSquareGridExample
from plotting.global_defs import globalAxesSettings
from figures_shared       import plotOneHist, getNoiseDataSpaces

import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)

from matplotlib import rc
rc('pdf', fonttype=42)
rc('mathtext', default='regular')

plt.rcParams['font.size'] = 11

outputDir = "."

NTrials=10
iterList  = ['g_AMPA_total', 'g_GABA_total']

noise_sigmas = [0, 150, 300]
exampleIdx   = [(1, 22), (1, 22), (1, 22)] # (row, col)
gridsDataRoot= 'output_local/even_spacing/grids'
velDataRoot = 'output_local/velocity'
shape = (31, 31)

grid_examples = 1
grids0        = 1
grids150      = 1
grids300      = 1
hists         = 1

##############################################################################

def drawGridSweeps(ax, dataSpace, iterList, NTrials=1, r=0, c=0, xLabelOn=True,
        xticks=True, exRows=[], exCols=[], exColor='white',
        cbar=False):
    yLabelText = '$g_E$ (nS)'
    if (xLabelOn):
        xLabelText = '$g_I$ (nS)'
    else:
        xLabelText = ''

    if (ax is None):
        ax = plt.gca()
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
            vmax=1.15,
            nansAs0=False)
    cax, kw = make_axes(ax, orientation='horizontal', shrink=0.8,
            pad=0.2)
    globalAxesSettings(cax)
    cb = plt.colorbar(ax=ax, cax=cax, ticks=MultipleLocator(0.5), **kw)
    cb.set_label('Gridness score')
    if (cbar == False):
        cax.set_visible(False)

    print("    max(G): {0}".format(np.max(G)))
    print("    min(G): {0}".format(np.min(G)))
    return ax, cax





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
        plotOneHist(G[filtIdx], normed=True)
    leg = []
    for s in noise_sigma:
        leg.append("{0}".format(int(s)))
    l = ax.legend(leg, loc=(0.8, 0.5), title='$\sigma$ (pA)', frameon=False,
            fontsize='x-small', ncol=1)
    plt.setp(l.get_title(), fontsize='x-small')

    ax.set_xlabel("Gridness score")
    ax.text(ylabelPos, 0.5, 'p(G)', rotation=90, transform=ax.transAxes,
            va='center', ha='right')
    #ax.set_ylabel("p(G)")

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.margins(0.05, 0.025)
    


##############################################################################
velSpaces  = getNoiseDataSpaces(velDataRoot,   noise_sigmas, shape)
gridSpaces = getNoiseDataSpaces(gridsDataRoot, noise_sigmas, shape)

exSz = 4
exMargin = 0.1
exGsCoords = 0, 0.0, 1.0-exMargin, 1.0-exMargin
exWspace=0.2
exHspace=0.2

sweepFigSize = (2.43, 3.33)
sweepLeft   = 0.15
sweepBottom = 0.1
sweepRight  = 0.99
sweepTop    = 0.95
if (grids0):

    # noise_sigma = 0 pA
    fig = plt.figure("sweeps0", figsize=sweepFigSize)
    exRows = [28, 15]
    exCols = [3, 15]
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    ax, cax = drawGridSweeps(ax, gridSpaces[0], iterList, NTrials=NTrials,
            r=exampleIdx[0][0], c=exampleIdx[0][1], xLabelOn=False,
            exRows=exRows, exCols=exCols, xticks=False)
    if (grid_examples):
        exLeft = 2
        exBottom = 24
        fname = outputDir + "/figure1_examples_0pA_0.png"
        plotSquareGridExample(exLeft, exBottom, exSz, fname, exampleIdx[0], ax,
                gridSpaces[0], iterList, exGsCoords, xlabel2=False,
                ylabel2=False, xlabel=False, ylabel=False, wspace=exWspace,
                hspace=exHspace, maxRate=True, plotGScore=False)
        
        exLeft = 18
        exBottom = 14
        fname = outputDir + "/figure1_examples_0pA_1.png"
        plotSquareGridExample(exLeft, exBottom, exSz, fname, exampleIdx[0], ax,
                gridSpaces[0], iterList, exGsCoords, xlabel2=False,
                ylabel2=False, xlabel=False, ylabel=False, wspace=exWspace,
                hspace=exHspace, maxRate=True, plotGScore=False)
    fname = outputDir + "/figure1_sweeps0.png"
    fig.savefig(fname, dpi=300, transparent=True)



if (grids150):
    # noise_sigma = 150 pA
    fig = plt.figure("sweeps150", figsize=sweepFigSize)
    exRows = [8, 2]
    exCols = [10, 9]
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    ax, cax = drawGridSweeps(ax, gridSpaces[1], iterList, NTrials=NTrials,
            r=exampleIdx[1][0], c=exampleIdx[1][1], xLabelOn=False,
            xticks=False, exRows=exRows, exCols=exCols) 
    if (grid_examples):
        exLeft = 2
        exBottom = 24
        fname = outputDir + "/figure1_examples_150pA_0.png"
        plotSquareGridExample(exLeft, exBottom, exSz, fname, exampleIdx[1], ax,
                gridSpaces[1], iterList, exGsCoords, xlabel2=False,
                ylabel2=False, xlabel=False, ylabel=False, wspace=exWspace,
                hspace=exHspace, maxRate=True, plotGScore=False)

        exLeft = 18
        exBottom = 14
        fname = outputDir + "/figure1_examples_150pA_1.png"
        plotSquareGridExample(exLeft, exBottom, exSz, fname, exampleIdx[1], ax,
                gridSpaces[1], iterList, exGsCoords, xlabel2=False,
                ylabel2=False, xlabel=False, ylabel=False, wspace=exWspace,
                hspace=exHspace, maxRate=True, plotGScore=False)
    fname = outputDir + "/figure1_sweeps150.png"
    fig.savefig(fname, dpi=300, transparent=True)



if (grids300):
    # noise_sigma = 300 pA
    fig = plt.figure("sweeps300", figsize=sweepFigSize)
    exRows = [16, 15]
    exCols = [6, 23]
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    _, cax = drawGridSweeps(ax, gridSpaces[2], iterList, NTrials=NTrials,
            r=exampleIdx[2][0], c=exampleIdx[2][1], xticks=True, exRows=exRows,
            exCols=exCols, exColor='black', cbar=True)
    if (grid_examples):
        exLeft = 2
        exBottom = 24
        fname = outputDir + "/figure1_examples_300pA_0.png"
        plotSquareGridExample(exLeft, exBottom, exSz, fname, exampleIdx[2], ax,
                gridSpaces[2], iterList, exGsCoords, xlabel2=False,
                ylabel2=False, xlabel=False, ylabel=False, wspace=exWspace,
                hspace=exHspace, maxRate=True, plotGScore=False)

        exLeft = 18
        exBottom = 14
        fname = outputDir + "/figure1_examples_300pA_1.png"
        plotSquareGridExample(exLeft, exBottom, exSz, fname, exampleIdx[2], ax,
                gridSpaces[2], iterList, exGsCoords, xlabel2=False,
                ylabel2=False, xlabel=False, ylabel=False, wspace=exWspace,
                hspace=exHspace, maxRate=True, plotGScore=False)

    fname = outputDir + "/figure1_sweeps300.png"
    fig.savefig(fname, dpi=300, transparent=True)


# Stats
if (hists):
    ylabelPos = -0.16
    fig = plt.figure(figsize=(3.7, 6))
    gs = GridSpec(3, 1, height_ratios=[0.8, 0.6, 1])

    ax_hist = plt.subplot(gs[0, 0])
    plotGridnessHistogram(gridSpaces, range(NTrials), ylabelPos=ylabelPos)

    ax_threshold = plt.subplot(gs[1, 0])
    plotGridnessThresholdComparison(gridSpaces, range(NTrials),
            thrList=np.arange(-0.4, 1.2, 0.05), ylabelPos=ylabelPos)

    gs.tight_layout(fig, rect=[0.1, 0, 1, 1], h_pad=3.0, pad=0.5)
    fname = outputDir + "/figure1_histograms.pdf"
    plt.savefig(fname, dpi=300, transparent=True)


