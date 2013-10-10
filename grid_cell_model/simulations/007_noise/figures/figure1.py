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

import EI_plotting as EI
from EI_plotting          import plotGridTrial, aggregate2DTrial, plotSquareGridExample
from plotting.global_defs import globalAxesSettings
from figures_shared       import plotOneHist, NoiseDataSpaces

import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)

from matplotlib import rc
rc('pdf', fonttype=42)
rc('mathtext', default='regular')

plt.rcParams['font.size'] = 11

outputDir = "."

NTrials=10
gridTrialNumList = np.arange(NTrials)
iterList  = ['g_AMPA_total', 'g_GABA_total']

noise_sigmas  = [0, 150, 300]
exampleIdx    = [(1, 22), (1, 22), (1, 22)] # (row, col)
bumpDataRoot  = None
velDataRoot   = None
gridsDataRoot = 'output_local/even_spacing/grids'
shape = (31, 31)

grid_examples = 0
grids0        = 1
grids150      = 1
grids300      = 1
hists         = 1
slices        = 1

##############################################################################



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
roots = NoiseDataSpaces.Roots(bumpDataRoot, velDataRoot, gridsDataRoot)
ps    = NoiseDataSpaces(roots, shape, noise_sigmas)

exSz = 4
exMargin = 0.1
exGsCoords = 0, 0.0, 1.0-exMargin, 1.0-exMargin
exWspace=0.2
exHspace=0.2

sweepFigSize = (2.43, 3.33)
sweepLeft   = 0.15
sweepBottom = 0.1
sweepRight  = 0.99
sweepTop    = 0.92

cbar_kwargs = {'label' : 'Gridness score',
    'orientation': 'horizontal',
    'shrink': 0.8,
    'pad' : 0.2,
    'ticks' : MultipleLocator(0.5)}

vmin = -0.5
vmax = 1.1

##############################################################################

varList = ['gridnessScore']

if (grids0):

    # noise_sigma = 0 pA
    fig = plt.figure("sweeps0", figsize=sweepFigSize)
    exRows = [28, 15]
    exCols = [3, 15]
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    plotGridTrial(ps.grids[0], varList, iterList, ps.noise_sigmas[0],
            trialNumList=gridTrialNumList,
            ax=ax,
            r=exampleIdx[0][0], c=exampleIdx[0][1],
            xlabel='', xticks=False,
            cbar=False, cbar_kwargs=cbar_kwargs,
            vmin=vmin, vmax=vmax,
            ignoreNaNs=True)
    if (grid_examples):
        exLeft = 2
        exBottom = 24
        fname = outputDir + "/figure1_examples_0pA_0.png"
        plotSquareGridExample(exLeft, exBottom, exSz, fname, exampleIdx[0], ax,
                ps.grids[0], iterList, exGsCoords, xlabel2=False,
                ylabel2=False, xlabel=False, ylabel=False, wspace=exWspace,
                hspace=exHspace, maxRate=True, plotGScore=False)
        
        exLeft = 18
        exBottom = 14
        fname = outputDir + "/figure1_examples_0pA_1.png"
        plotSquareGridExample(exLeft, exBottom, exSz, fname, exampleIdx[0], ax,
                ps.grids[0], iterList, exGsCoords, xlabel2=False,
                ylabel2=False, xlabel=False, ylabel=False, wspace=exWspace,
                hspace=exHspace, maxRate=True, plotGScore=False)
    fname = outputDir + "/figure1_sweeps0.png"
    fig.savefig(fname, dpi=300, transparent=True)
    plt.close()



if (grids150):
    # noise_sigma = 150 pA
    fig = plt.figure("sweeps150", figsize=sweepFigSize)
    exRows = [8, 2]
    exCols = [10, 9]
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    plotGridTrial(ps.grids[1], varList, iterList, ps.noise_sigmas[1],
            trialNumList=gridTrialNumList,
            ax=ax,
            r=exampleIdx[1][0], c=exampleIdx[1][1],
            xlabel='', xticks=False,
            cbar=False, cbar_kwargs=cbar_kwargs,
            vmin=vmin, vmax=vmax,
            ignoreNaNs=True)
    if (grid_examples):
        exLeft = 2
        exBottom = 24
        fname = outputDir + "/figure1_examples_150pA_0.png"
        plotSquareGridExample(exLeft, exBottom, exSz, fname, exampleIdx[1], ax,
                ps.grids[1], iterList, exGsCoords, xlabel2=False,
                ylabel2=False, xlabel=False, ylabel=False, wspace=exWspace,
                hspace=exHspace, maxRate=True, plotGScore=False)

        exLeft = 18
        exBottom = 14
        fname = outputDir + "/figure1_examples_150pA_1.png"
        plotSquareGridExample(exLeft, exBottom, exSz, fname, exampleIdx[1], ax,
                ps.grids[1], iterList, exGsCoords, xlabel2=False,
                ylabel2=False, xlabel=False, ylabel=False, wspace=exWspace,
                hspace=exHspace, maxRate=True, plotGScore=False)
    fname = outputDir + "/figure1_sweeps150.png"
    fig.savefig(fname, dpi=300, transparent=True)
    plt.close()



if (grids300):
    # noise_sigma = 300 pA
    fig = plt.figure("sweeps300", figsize=sweepFigSize)
    exRows = [16, 15]
    exCols = [6, 23]
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    plotGridTrial(ps.grids[2], varList, iterList, ps.noise_sigmas[2],
            trialNumList=gridTrialNumList,
            ax=ax,
            r=exampleIdx[2][0], c=exampleIdx[2][1],
            cbar_kwargs=cbar_kwargs,
            vmin=vmin, vmax=vmax,
            ignoreNaNs=True)
    if (grid_examples):
        exLeft = 2
        exBottom = 24
        fname = outputDir + "/figure1_examples_300pA_0.png"
        plotSquareGridExample(exLeft, exBottom, exSz, fname, exampleIdx[2], ax,
                ps.grids[2], iterList, exGsCoords, xlabel2=False,
                ylabel2=False, xlabel=False, ylabel=False, wspace=exWspace,
                hspace=exHspace, maxRate=True, plotGScore=False)

        exLeft = 18
        exBottom = 14
        fname = outputDir + "/figure1_examples_300pA_1.png"
        plotSquareGridExample(exLeft, exBottom, exSz, fname, exampleIdx[2], ax,
                ps.grids[2], iterList, exGsCoords, xlabel2=False,
                ylabel2=False, xlabel=False, ylabel=False, wspace=exWspace,
                hspace=exHspace, maxRate=True, plotGScore=False)

    fname = outputDir + "/figure1_sweeps300.png"
    fig.savefig(fname, dpi=300, transparent=True)
    plt.close()


# Stats
sliceFigSize = (3.7, 2)
sliceLeft   = 0.2
sliceBottom = 0.3
sliceRight  = 0.99
sliceTop    = 0.85
if (hists):
    ylabelPos = -0.16
    fig = plt.figure(figsize=sliceFigSize)
    ax = fig.add_axes(Bbox.from_extents(sliceLeft, sliceBottom, sliceRight,
        sliceTop))
    plotGridnessThresholdComparison(ps.grids, range(NTrials),
            thrList=np.arange(-0.4, 1.2, 0.05), ylabelPos=ylabelPos)
    fname = outputDir + "/figure1_threshold_comparison.pdf"
    plt.savefig(fname, dpi=300, transparent=True)


if (slices):
    ylabelPos = -0.16
    idx_horizontal = 15
    fig = plt.figure(figsize=sliceFigSize)
    ax = fig.add_axes(Bbox.from_extents(sliceLeft, sliceBottom, sliceRight,
        sliceTop))
    EI.plotGridnessSlice(ps, idx_horizontal, slice(None), ax=ax)
    ax.yaxis.set_major_locator(MultipleLocator(0.4))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.set_ylim([-0.5, 1.21])
    fname = "figure1_slice_horizontal.pdf"
    plt.savefig(fname, dpi=300, transparent=True)
    plt.close()

    idx_vertical = 15
    fig = plt.figure(figsize=sliceFigSize)
    ax = fig.add_axes(Bbox.from_extents(sliceLeft, sliceBottom, sliceRight,
        sliceTop))
    EI.plotGridnessSlice(ps, slice(None), idx_vertical, ax=ax)
    ax.yaxis.set_major_locator(MultipleLocator(0.4))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.set_ylim([-0.5, 1.21])
    fname = "figure1_slice_vertical.pdf"
    plt.savefig(fname, dpi=300, transparent=True)
    plt.close()



