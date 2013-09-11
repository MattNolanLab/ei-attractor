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
from matplotlib.pyplot   import figure, subplot, plot, savefig, close, \
        errorbar
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator, AutoMinorLocator, LinearLocator, MaxNLocator, \
        ScalarFormatter
from matplotlib.colorbar import make_axes
from matplotlib.transforms import Bbox

from parameters  import JobTrialSpace2D
from EI_plotting import plotBumpSigmaTrial, computeYX, aggregate2DTrial, \
        aggregate2D, drawEIRectSelection, drawBumpExamples, plotVelTrial
from plotting.grids import plotGridRateMap, plotAutoCorrelation, plotSpikes2D
from plotting.global_defs import globalAxesSettings, createColorbar
from figures_shared import plotOneHist

import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)

from matplotlib import rc
rc('pdf', fonttype=42)
rc('mathtext', default='regular')

plt.rcParams['font.size'] = 10

outputDir = "."

NTrials=10
iterList  = ['g_AMPA_total', 'g_GABA_total']

noise_sigmas = [0, 150, 300]
exampleIdx   = [(0, 0), (0, 0), (0, 0)] # (row, col)
bumpDataRoot= 'output_local/one_to_one'
velDataRoot = 'output_local/velocity'
bumpShape = (40, 40)
velShape  = (30, 30)

bumpExamples = 1
bumpSweep0   = 0
bumpSweep150 = 0
bumpSweep300 = 0
velExamples  = 1
velSweep0    = 0
velSweep150  = 0
velSweep300  = 0
hists        = 1 
velLines     = 1

##############################################################################

def getNoiseRootDir(prefix, noise_sigma):
    return  "{0}/EI_param_sweep_{1}pA".format(prefix, int(noise_sigma))


def getNoiseRoots(prefix, noise_sigmas):
    roots = []
    for s in noise_sigmas:
        roots.append(getNoiseRootDir(prefix, s))
    return roots



def drawBumpSweeps(ax, dataSpace, iterList, noise_sigma, NTrials=1, r=0, c=0, yLabelOn=True,
        yticks=True, cbar=False):
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
    cax.yaxis.set_minor_locator(AutoMinorLocator(2))

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
    plt.close()

    # Draw the selection into the EI plot
    if (sweep_ax is not None):
        exRow, exCol = exIdx
        Y, X = computeYX(sweepDataSpace, iterList, r=exRow, c=exCol)
        drawEIRectSelection(sweep_ax, exRect, X, Y, color=rectColor)


###############################################################################
def drawVelSweeps(ax, dataSpace, iterList, noise_sigma, r=0, c=0, yLabelOn=True,
        yticks=True, cbar=False):
    xLabelText = '$w_I$ (nS)'
    if (yLabelOn):
        yLabelText = '$w_E$ (nS)'
    else:
        yLabelText = ''

    if (ax is None):
        ax = plt.gca()

    varList = ['lineFitErr']
    G = plotVelTrial(dataSpace, varList, iterList,
            xlabel=xLabelText,
            ylabel=yLabelText,
            colorBar=False,
            yticks=yticks,
            vmin=0,
            vmax=60)
    plt.set_cmap('jet')
    cax, kw = make_axes(ax, orientation='vertical', shrink=0.8,
            pad=0.05)
    globalAxesSettings(cax)
    cb = plt.colorbar(ax=ax, cax=cax, ticks=MultipleLocator(20), **kw)
    cb.set_label('Fit error (neurons/s)')
    if (cbar == False):
        cax.set_visible(False)
    ax.set_title('$\sigma$ = {0} pA'.format(int(noise_sigma)))

    cax.yaxis.set_minor_locator(AutoMinorLocator(2))

    return ax, cax


def plotVelHistogram(spList, varList, xlabel="", ylabel="", **kw):
    noise_sigma = [0, 150, 300]
    colors = ['red', 'green', 'blue']
    range = kw.get('range')
    plotLegend = kw.pop('plotLegend', False)

    ax = plt.gca()
    plt.hold('on')
    globalAxesSettings(ax)

    for idx, sp in enumerate(spList):
        var = np.abs(aggregate2D(sp, varList, funReduce=None))
        filtIdx = np.logical_not(np.isnan(var))
        if (range is not None):
            var[var < range[0]] = range[0]
            var[var > range[1]] = range[1]
        plotOneHist(var[filtIdx], normed=True, **kw)

    if (plotLegend):
        leg = []
        for s in noise_sigma:
            leg.append("{0}".format(int(s)))
        l = ax.legend(leg, loc=(0.75, 0.5), title='$\sigma$ (pA)',
                frameon=False, fontsize='x-small', ncol=1)
        plt.setp(l.get_title(), fontsize='x-small')

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    f = ScalarFormatter(useMathText=True)
    f.set_scientific(True)
    f.set_powerlimits([0, 3])
    ax.yaxis.set_major_formatter(f)
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    return ax

def plotErrHistogram(spList, varList, **kw):
    ax = plotVelHistogram(spList, varList, range=[0, 60], **kw)

    ax.xaxis.set_major_locator(MultipleLocator(20))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(MaxNLocator(3))
    ax.set_ylim([-0.0025, 1.01*0.2])
    ax.margins(0.01)
    
def plotSlopeHistogram(spList, varList, **kw):
    ax = plotVelHistogram(spList, varList, **kw)

    #ax.xaxis.set_major_locator(MultipleLocator(20))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(MaxNLocator(3))
    ax.set_ylim([-0.0025, 10])
    ax.margins(0.01)
    

def plotSlopes(ax, dataSpace, pos, **kw):
    # kwargs
    trialNum = kw.pop('trialNum', 0)
    markersize = kw.pop('markersize', 4)
    color = kw.pop('color', 'blue')

    r = pos[0]
    c = pos[1]
    d = dataSpace[r][c].getAllTrialsAsDataSet().data
    a = d['analysis']
    IvelVec = dataSpace[r][c][trialNum].data['IvelVec']
    slopes = a['bumpVelAll']
    lineFit = a['lineFitLine']
    lineFitRange = a['fitRange']

    nTrials = slopes.shape[0]
    avgSlope = np.mean(slopes, axis=0)
    stdErrSlope = np.std(slopes, axis=0) / np.sqrt(nTrials)

    if (ax is None):
        ax = plt.gca()
    plt.hold('on')
    globalAxesSettings(ax)

    r = lineFitRange
    errorbar(IvelVec, -avgSlope, stdErrSlope, fmt='o-', markersize=markersize,
            color=color, alpha=0.5, **kw)
    plot(IvelVec[0:r], -lineFit, '-', linewidth=2, color=color, **kw)

def plotAllSlopes(ax, spList, positions, **kw):
    colors = kw.pop('colors', ('blue', 'green', 'red'))

    for idx, dataSpace in enumerate(spList):
        kw['color'] = colors[idx]
        plotSlopes(ax, dataSpace, positions[idx], **kw)

    ax.set_xlabel('Velocity current (pA)')
    ax.set_ylabel('$v_{bump}$ (neurons/s)')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_major_locator(MultipleLocator(50))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_major_locator(MultipleLocator(20))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.margins(0.05)
    


###############################################################################
bumpRoots = getNoiseRoots(bumpDataRoot, noise_sigmas)
bumpDataSpace0   = JobTrialSpace2D(bumpShape, bumpRoots[0])
bumpDataSpace150 = JobTrialSpace2D(bumpShape, bumpRoots[1])
bumpDataSpace300 = JobTrialSpace2D(bumpShape, bumpRoots[2])

velRoots = getNoiseRoots(velDataRoot, noise_sigmas)
velDataSpace0   = JobTrialSpace2D(velShape, velRoots[0])
velDataSpace150 = JobTrialSpace2D(velShape, velRoots[1])
velDataSpace300 = JobTrialSpace2D(velShape, velRoots[2])

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

histFigsize =(2.6, 1.7)
histLeft    = 0.22
histBottom  = 0.3
histRight   = 0.95
histTop     = 0.86

if (bumpSweep0):
    # noise_sigma = 0 pA
    fig = figure("sweeps0", figsize=sweepFigSize)
    exRows = [28, 15]
    exCols = [3, 15]
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    ax, cax = drawBumpSweeps(ax, bumpDataSpace0, iterList,
            noise_sigma=noise_sigmas[0], NTrials=NTrials, cbar=False)
    if (bumpExamples):
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



if (bumpSweep150):
    # noise_sigma = 150 pA
    fig = figure("sweeps150", figsize=sweepFigSize)
    exRows = [8, 2]
    exCols = [10, 9]
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    ax, cax = drawBumpSweeps(ax, bumpDataSpace150, iterList,
            noise_sigma=noise_sigmas[1],  NTrials=NTrials, yLabelOn=False,
            yticks=False) 
    if (bumpExamples):
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



if (bumpSweep300):
    # noise_sigma = 300 pA
    fig = figure("sweeps300", figsize=sweepFigSize)
    exRows = [16, 15]
    exCols = [6, 23]
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    ax.set_clip_on(False)
    _, cax = drawBumpSweeps(ax, bumpDataSpace300, iterList,
            noise_sigma=noise_sigmas[2],  NTrials=NTrials, yLabelOn=False,
            yticks=False, cbar=True)
    if (bumpExamples):
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

###############################################################################

velSpList = [velDataSpace0, velDataSpace150, velDataSpace300]

if (velSweep0):
    # noise_sigma = 0 pA
    fig = figure("bumpSweeps0", figsize=sweepFigSize)
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    ax, cax = drawVelSweeps(ax, velDataSpace0, iterList,
            noise_sigma=noise_sigmas[0], cbar=False)
    fname = outputDir + "/figure1_err_sweeps0.png"
    fig.savefig(fname, dpi=300, transparent=True)


if (velSweep150):
    # noise_sigma = 150 pA
    fig = figure("bumpSweeps150", figsize=sweepFigSize)
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    ax, cax = drawVelSweeps(ax, velDataSpace150, iterList, yLabelOn=False,
            yticks=False, noise_sigma=noise_sigmas[1], cbar=False)
    fname = outputDir + "/figure1_err_sweeps150.png"
    fig.savefig(fname, dpi=300, transparent=True)


if (velSweep300):
    # noise_sigma = 300 pA
    fig = figure("bumpSweeps300", figsize=sweepFigSize)
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    ax, cax = drawVelSweeps(ax, velDataSpace300, iterList, yLabelOn=False,
            yticks=False, noise_sigma=noise_sigmas[2], cbar=True)
    fname = outputDir + "/figure1_err_sweeps300.png"
    fig.savefig(fname, dpi=300, transparent=True)

# Stats
if (hists):
    fig = figure(figsize=histFigsize)
    ax = fig.add_axes(Bbox.from_extents(histLeft, histBottom, histRight,
        histTop))
    plotErrHistogram(velSpList, ['lineFitErr'], xlabel='Fit error (neurons/s)',
            ylabel='p(error)')
    fname = outputDir + "/figure1_err_histograms.pdf"
    savefig(fname, dpi=300, transparent=True)

    fig = figure(figsize=histFigsize)
    ax = fig.add_axes(Bbox.from_extents(histLeft, histBottom, histRight,
        histTop))
    plotSlopeHistogram(velSpList, ['lineFitSlope'], xlabel='Slope (neurons/s/pA)',
            ylabel='p(slope)', plotLegend=True)
    fname = outputDir + "/figure1_slope_histograms.pdf"
    savefig(fname, dpi=300, transparent=True)


if (velLines):
    positions = ((18, 2), (15, 5), (10,2))
    fig = figure(figsize=(2.5, histFigsize[1]))
    ax = fig.add_axes(Bbox.from_extents(0.3, histBottom, histRight,
        histTop))
    plotAllSlopes(ax, velSpList, positions)
    fname = outputDir + "/figure1_slope_examples.pdf"
    savefig(fname, dpi=300, transparent=True)

