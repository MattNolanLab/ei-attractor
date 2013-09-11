#!/usr/bin/env python
#
#   figure1.py
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
from matplotlib.colorbar import make_axes
from matplotlib.transforms import Bbox

from parameters  import JobTrialSpace2D
from EI_plotting import plotBumpSigmaTrial, computeYX, aggregate2DTrial, \
        aggregate2D, drawEIRectSelection, drawBumpExamples
from plotting.grids import plotGridRateMap, plotAutoCorrelation, plotSpikes2D
from plotting.global_defs import globalAxesSettings, createColorbar
from figures_shared import plotOneHist

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
exampleIdx   = [(0, 0), (0, 0), (0, 0)] # (row, col)
bumpDataRoot= 'output_local/one_to_one'
velDataRoot = 'output_local/velocity'
bumpShape = (40, 40)

examples      = 1
sweep0        = 1
sweep150      = 1
sweep300      = 1

##############################################################################

def getNoiseRootDir(prefix, noise_sigma):
    return  "{0}/EI_param_sweep_{1}pA".format(prefix, int(noise_sigma))


def getNoiseRoots(prefix, noise_sigmas):
    roots = []
    for s in noise_sigmas:
        roots.append(getNoiseRootDir(prefix, s))
    return roots



def drawSweeps(ax, dataSpace, iterList, noise_sigma, NTrials=1, r=0, c=0, yLabelOn=True,
        yticks=True, exColor='white', cbar=False):
    xLabelText = '$w_I$ (nS)'
    if (yLabelOn):
        yLabelText = '$w_E$ (nS)'
    else:
        yLabelText = ''

    if (ax is None):
        ax = plt.gca()

    varList = ['bump_e', 'sigma']
    G = plotBumpSigmaTrial(dataSpace, varList, iterList,
            trialNumList=range(NTrials),
            xlabel=xLabelText,
            ylabel=yLabelText,
            colorBar=False,
            clBarLabel = "Gridness score",
            clbarNTicks=3,
            yticks=yticks,
            vmin=0,
            vmax=10)
    plt.set_cmap('jet_r')
    cax, kw = make_axes(ax, orientation='vertical', shrink=0.8,
            pad=0.05)
    globalAxesSettings(cax)
    cb = plt.colorbar(ax=ax, cax=cax, ticks=MultipleLocator(5), **kw)
    cb.set_label('Bump $\sigma$ (neurons)')
    if (cbar == False):
        cax.set_visible(False)
    ax.set_title('$\sigma$ = {0} pA'.format(int(noise_sigma)))

    return ax, cax


def plotBumpExample(exLeft, exBottom, w, h, fileName, exIdx, sweep_ax,
        sweepDataSpace, iterList, exGsCoords, **kw):
    #keyword
    wspace = kw.pop('wspace', 0)
    hspace = kw.pop('hspace', 0)
    figSize = kw.pop('figSize', (1.8, 1))
    rectColor = kw.pop('rectColor', 'black')

    # Create the example plot
    fig = plt.figure(figsize=figSize)

    exRect = [exLeft, exBottom, exLeft+w-1, exBottom+h-1]
    gs = drawBumpExamples(sweepDataSpace, exRect, iterList,
            gsCoords=exGsCoords, xlabel=False, ylabel=False, xlabel2=False,
            ylabel2=False, fontsize='xx-small', rateYPos=1.05, rateXPos=0.98,
            **kw)
    gs.update(wspace=wspace, hspace=hspace)
    plt.savefig(fileName, dpi=300, transparent=False)

    # Draw the selection into the EI plot
    if (sweep_ax is not None):
        exRow, exCol = exIdx
        Y, X = computeYX(sweepDataSpace, iterList, r=exRow, c=exCol)
        drawEIRectSelection(sweep_ax, exRect, X, Y, color=rectColor)




###############################################################################
bumpRoots = getNoiseRoots(bumpDataRoot, noise_sigmas)
bumpDataSpace0   = JobTrialSpace2D(bumpShape, bumpRoots[0])
bumpDataSpace150 = JobTrialSpace2D(bumpShape, bumpRoots[1])
bumpDataSpace300 = JobTrialSpace2D(bumpShape, bumpRoots[2])

#velRoots = getNoiseRoots(velDataRoot, noise_sigmas)
#velDataSpace0   = JobTrialSpace2D(shape, velRoots[0])
#velDataSpace150 = JobTrialSpace2D(shape, velRoots[1])
#velDataSpace300 = JobTrialSpace2D(shape, velRoots[2])

exW = 4
exH = 2
exMargin = 0.075
exGsCoords = 0.02, 0, 0.98, 1.0-exMargin
exWspace=0.2
exHspace=0.15

sweepFigSize = (3.4, 2.1)
sweepLeft   = 0.15
sweepBottom = 0.2
sweepRight  = 0.9
sweepTop    = 0.85

if (sweep0):
    # noise_sigma = 0 pA
    fig = figure("sweeps0", figsize=sweepFigSize)
    exRows = [28, 15]
    exCols = [3, 15]
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    ax, cax = drawSweeps(ax, bumpDataSpace0, iterList,
            noise_sigma=noise_sigmas[0], NTrials=NTrials, cbar=False)
    if (examples):
        exLeft = 1
        exBottom = 24
        fname = outputDir + "/figure1_examples_0pA_0.png"
        plotBumpExample(exLeft, exBottom, exW, exH, fname, exampleIdx[0],
                ax, bumpDataSpace0, iterList, exGsCoords, wspace=exWspace,
                hspace=exHspace, rectColor='red')

        exLeft = 25
        exBottom = 15
        fname = outputDir + "/figure1_examples_0pA_1.png"
        plotBumpExample(exLeft, exBottom, exW, exH, fname, exampleIdx[0],
                ax, bumpDataSpace0, iterList, exGsCoords, wspace=exWspace,
                hspace=exHspace, rectColor='red')

    fname = outputDir + "/figure1_sweeps0.png"
    fig.savefig(fname, dpi=300, transparent=True)



if (sweep150):
    # noise_sigma = 150 pA
    fig = figure("sweeps150", figsize=sweepFigSize)
    exRows = [8, 2]
    exCols = [10, 9]
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    ax, cax = drawSweeps(ax, bumpDataSpace150, iterList,
            noise_sigma=noise_sigmas[1],  NTrials=NTrials, yLabelOn=False,
            yticks=False) 
    if (examples):
        exLeft = 1
        exBottom = 24
        fname = outputDir + "/figure1_examples_150pA_0.png"
        plotBumpExample(exLeft, exBottom, exW, exH, fname, exampleIdx[1],
                ax, bumpDataSpace150, iterList, exGsCoords, wspace=exWspace,
                hspace=exHspace, rectColor='red')

        exLeft = 25
        exBottom = 15
        fname = outputDir + "/figure1_examples_150pA_1.png"
        plotBumpExample(exLeft, exBottom, exW, exH, fname, exampleIdx[1],
                ax, bumpDataSpace150, iterList, exGsCoords, wspace=exWspace,
                hspace=exHspace, rectColor='black')


    fname = outputDir + "/figure1_sweeps150.png"
    fig.savefig(fname, dpi=300, transparent=True)



if (sweep300):
    # noise_sigma = 300 pA
    fig = figure("sweeps300", figsize=sweepFigSize)
    exRows = [16, 15]
    exCols = [6, 23]
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    ax.set_clip_on(False)
    _, cax = drawSweeps(ax, bumpDataSpace300, iterList,
            noise_sigma=noise_sigmas[2],  NTrials=NTrials, yLabelOn=False,
            yticks=False, exColor='black', cbar=True)
    if (examples):
        exLeft = 1
        exBottom = 24
        fname = outputDir + "/figure1_examples_300pA_0.png"
        plotBumpExample(exLeft, exBottom, exW, exH, fname, exampleIdx[2],
                ax, bumpDataSpace300, iterList, exGsCoords, wspace=exWspace,
                hspace=exHspace, rectColor='red')

        exLeft = 25
        exBottom = 15
        fname = outputDir + "/figure1_examples_300pA_1.png"
        plotBumpExample(exLeft, exBottom, exW, exH, fname, exampleIdx[2],
                ax, bumpDataSpace300, iterList, exGsCoords, wspace=exWspace,
                hspace=exHspace, rectColor='black')


    fname = outputDir + "/figure1_sweeps300.png"
    fig.savefig(fname, dpi=300, transparent=True)


