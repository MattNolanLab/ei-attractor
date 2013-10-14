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
import matplotlib.ticker as ti
from matplotlib.gridspec   import GridSpec
from matplotlib.colorbar   import make_axes
from matplotlib.transforms import Bbox

import EI_plotting as EI
from plotting.global_defs import globalAxesSettings
from figures_shared       import NoiseDataSpaces

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
gridsDataRoot= 'output_local/even_spacing/grids'
bumpDataRoot= 'output_local/even_spacing/gamma_bump'
velDataRoot = 'output_local/even_spacing/velocity'
shape = (31, 31)

bumpExamples      = 0
bumpSweep0        = 0
bumpSweep150      = 0
bumpSweep300      = 0
velExamples       = 0
velSweep0         = 0
velSweep150       = 0
velSweep300       = 0
velLines          = 1
gridness_vs_error = 0

##############################################################################


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
    gs = EI.drawBumpExamples(sweepDataSpace, exRect, iterList,
            gsCoords=exGsCoords, xlabel=False, ylabel=False, xlabel2=False,
            ylabel2=False, fontsize='xx-small', rateYPos=1.05, rateXPos=0.98,
            **kw)
    gs.update(wspace=wspace, hspace=hspace)
    plt.savefig(fileName, dpi=300, transparent=False)
    plt.close()

    # Draw the selection into the EI plot
    if (sweep_ax is not None):
        exRow, exCol = exIdx
        Y, X = EI.computeYX(sweepDataSpace, iterList, r=exRow, c=exCol)
        EI.drawEIRectSelection(sweep_ax, exRect, X, Y, color=rectColor)



def plotSlopes(ax, dataSpace, pos, noise_sigma, **kw):
    # kwargs
    trialNum   = kw.pop('trialNum', 0)
    markersize = kw.pop('markersize', 4)
    color      = kw.pop('color', 'blue')
    xlabel     = kw.pop('xlabel', 'Velocity current (pA)')
    ylabel     = kw.pop('ylabel', 'Bump speed (neurons/s)')
    xticks     = kw.pop('xticks', True)
    yticks     = kw.pop('yticks', True)
    g_ann      = kw.pop('g_ann', True)
    sigma_ann  = kw.pop('sigma_ann', True)
    kw['markeredgecolor'] = color

    r = pos[0]
    c = pos[1]
    d = dataSpace[r][c].getAllTrialsAsDataSet().data
    a = d['analysis']
    IvelVec = dataSpace[r][c][trialNum].data['IvelVec']
    slopes = a['bumpVelAll']
    lineFit = a['lineFitLine']
    fitIvelVec = a['fitIvelVec']

    nTrials = slopes.shape[0]
    avgSlope = np.mean(slopes, axis=0)
    stdSlope = np.std(slopes, axis=0)

    if (ax is None):
        ax = plt.gca()
    plt.hold('on')
    globalAxesSettings(ax)

    ax.plot(IvelVec, slopes.T, 'o', color='none', markersize=markersize, **kw)
    #ax.errorbar(IvelVec, avgSlope, stdSlope, fmt='o-',
    #        markersize=markersize, color=color, alpha=0.5, **kw)
    ax.plot(fitIvelVec, lineFit, '-', linewidth=1, color=color, **kw)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_major_locator(ti.MultipleLocator(50))
    ax.xaxis.set_minor_locator(ti.AutoMinorLocator(5))
    ax.yaxis.set_major_locator(ti.MultipleLocator(20))
    ax.yaxis.set_minor_locator(ti.AutoMinorLocator(2))
    ax.margins(0.05)

    if (not xticks):
        ax.xaxis.set_ticklabels([])
    if (not yticks):
        ax.yaxis.set_ticklabels([])

    # Annotations
    if (sigma_ann):
        sigma_txt = '$\sigma$ = {0} pA'.format(noise_sigma)
    else:
        sigma_txt = ''

    if (g_ann):
        Y, X = EI.computeVelYX(dataSpace, iterList, r, c)
        gE = Y[r, c]
        gI = X[r, c]
        g_txt = '$g_E$ = {0}, $g_I$ = {1} nS'.format(gE, gI)
    else:
        g_txt = ''

    txt = '{0}\n{1}'.format(g_txt, sigma_txt)
    ax.text(0.05, 0.85, txt, transform=ax.transAxes, va='bottom',
            ha='left', size='x-small')




def plotGridnessVsFitErr(spListGrids, spListVelocity, trialNumList,
        ylabelPos=-0.2, maxErr=None):
    GVars = ['gridnessScore']
    errVars = ['lineFitErr']
    slopeVars = ['lineFitSlope']
    noise_sigma = [0, 150, 300]
    markers = ['o', '^', '*']
    colors = ['blue', 'green', 'red']
    errMins = []
    errMaxs = []

    ax = plt.gca()
    plt.hold('on')
    globalAxesSettings(ax)
    ax.set_yscale('log')

    for idx, (spGrids, spVel) in enumerate(zip(spListGrids, spListVelocity)):
        G = EI.aggregate2DTrial(spGrids, GVars, trialNumList).flatten()
        errs = EI.aggregate2D(spVel, errVars, funReduce=np.sum).flatten()
        #slopes = np.abs(EI.aggregate2D(spVel, slopeVars,
        #    funReduce=None).flatten())
        i = np.logical_not(np.logical_and(np.isnan(G), np.isnan(errs)))
        ax.scatter(G[i], errs[i],  s=5, marker=markers[idx], 
                color=colors[idx], edgecolors='None')
        errMins.append(np.min(errs[i]))
        errMaxs.append(np.max(errs[i]))

    if (maxErr is not None):
        ax.set_ylim([0, maxErr])
    else:
        ax.set_ylim([0, None])

    leg = []
    for s in noise_sigma:
        leg.append("{0}".format(int(s)))
    l = ax.legend(leg, loc=(0.3, 0.85), title='$\sigma$ (pA)', frameon=False,
            fontsize='small', ncol=3, columnspacing=1.5)
    plt.setp(l.get_title(), fontsize='small')

    ax.set_xlabel("Gridness score")
    ax.set_ylabel('Error (nrns/s/data point)')
    #ax.text(ylabelPos, 0.5, 'Error (nrns/s/data point)', rotation=90, transform=ax.transAxes,
    #        va='center', ha='right')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_major_locator(ti.MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(ti.AutoMinorLocator(2))
    ax.set_xmargin(0.05)
    ax.autoscale_view(tight=True)
    ax.set_ylim(np.min(errMins), np.max(errMaxs)*1.5)


###############################################################################
roots = NoiseDataSpaces.Roots(bumpDataRoot, velDataRoot, gridsDataRoot)
ps    = NoiseDataSpaces(roots, shape, noise_sigmas)


exW = 4
exH = 2
exMargin = 0.075
exGsCoords = 0.02, 0, 0.98, 1.0-exMargin
exWspace=0.2
exHspace=0.15

sweepFigSize = (2.5, 2.1)
sweepLeft   = 0.15
sweepBottom = 0.2
sweepRight  = 0.87
sweepTop    = 0.85

velFigsize =(2.6, 2)
velLeft    = 0.3
velBottom  = 0.35
velRight   = 0.95
velTop     = 0.65


##############################################################################
sigmaVarList = ['bump_e', 'sigma']
bumpTrialNumList = np.arange(5)
bump_vmin = 0
bump_vmax = 10
bump_cbar_kwargs = dict(
        orientation='vertical',
        shrink=0.8, pad=0.05,
        ticks=ti.MultipleLocator(5),
        label='Bump $\sigma$ (neurons)',
        extend='max', extendfrac=0.1)

if (bumpSweep0):
    # noise_sigma = 0 pA
    fig = plt.figure("sweeps0", figsize=sweepFigSize)
    exRows = [28, 15]
    exCols = [3, 15]
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    EI.plotBumpSigmaTrial(ps.bumpGamma[0], sigmaVarList, iterList,
            noise_sigma=ps.noise_sigmas[0],
            ax=ax,
            trialNumList=bumpTrialNumList,
            cbar=False, cbar_kwargs=bump_cbar_kwargs,
            vmin=bump_vmin, vmax=bump_vmax)
    if (bumpExamples):
        exLeft = 1
        exBottom = 24
        fname = outputDir + "/figure2_examples_0pA_0.png"
        plotBumpExample(exLeft, exBottom, exW, exH, fname, exampleIdx[0],
                ax, ps.bumpGamma[0], iterList, exGsCoords, wspace=exWspace,
                hspace=exHspace, rectColor='red')

        exLeft = 25
        exBottom = 15
        fname = outputDir + "/figure2_examples_0pA_1.png"
        plotBumpExample(exLeft, exBottom, exW, exH, fname, exampleIdx[0],
                ax, ps.bumpGamma[0], iterList, exGsCoords, wspace=exWspace,
                hspace=exHspace, rectColor='red')

    fname = outputDir + "/figure2_sweeps0.png"
    fig.savefig(fname, dpi=300, transparent=True)



if (bumpSweep150):
    # noise_sigma = 150 pA
    fig = plt.figure("sweeps150", figsize=sweepFigSize)
    exRows = [8, 2]
    exCols = [10, 9]
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    EI.plotBumpSigmaTrial(ps.bumpGamma[1], sigmaVarList, iterList,
            noise_sigma=noise_sigmas[1],
            ax=ax,
            trialNumList=bumpTrialNumList,
            cbar=False, cbar_kwargs=bump_cbar_kwargs,
            ylabel='', yticks=False,
            vmin=bump_vmin, vmax=bump_vmax)
    if (bumpExamples):
        exLeft = 1
        exBottom = 24
        fname = outputDir + "/figure2_examples_150pA_0.png"
        plotBumpExample(exLeft, exBottom, exW, exH, fname, exampleIdx[1],
                ax, ps.bumpGamma[1], iterList, exGsCoords, wspace=exWspace,
                hspace=exHspace, rectColor='red')

        exLeft = 25
        exBottom = 15
        fname = outputDir + "/figure2_examples_150pA_1.png"
        plotBumpExample(exLeft, exBottom, exW, exH, fname, exampleIdx[1],
                ax, ps.bumpGamma[1], iterList, exGsCoords, wspace=exWspace,
                hspace=exHspace, rectColor='black')


    fname = outputDir + "/figure2_sweeps150.png"
    fig.savefig(fname, dpi=300, transparent=True)



if (bumpSweep300):
    # noise_sigma = 300 pA
    fig = plt.figure("sweeps300", figsize=sweepFigSize)
    exRows = [16, 15]
    exCols = [6, 23]
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    ax.set_clip_on(False)
    EI.plotBumpSigmaTrial(ps.bumpGamma[2], sigmaVarList, iterList,
            noise_sigma=noise_sigmas[2],
            ax=ax,
            trialNumList=bumpTrialNumList,
            cbar=True, cbar_kwargs=bump_cbar_kwargs,
            ylabel='', yticks=False,
            vmin=bump_vmin, vmax=bump_vmax)
    if (bumpExamples):
        exLeft = 1
        exBottom = 24
        fname = outputDir + "/figure2_examples_300pA_0.png"
        plotBumpExample(exLeft, exBottom, exW, exH, fname, exampleIdx[2],
                ax, ps.bumpGamma[2], iterList, exGsCoords, wspace=exWspace,
                hspace=exHspace, rectColor='red')

        exLeft = 25
        exBottom = 15
        fname = outputDir + "/figure2_examples_300pA_1.png"
        plotBumpExample(exLeft, exBottom, exW, exH, fname, exampleIdx[2],
                ax, ps.bumpGamma[2], iterList, exGsCoords, wspace=exWspace,
                hspace=exHspace, rectColor='black')



    fname = outputDir + "/figure2_sweeps300.png"
    fig.savefig(fname, dpi=300, transparent=True)

###############################################################################

slopeVarList = ['lineFitSlope']
slope_vmin = 0
slope_vmax = 1.6
std_vmin = 0
std_vmax = 14

slope_cbar_kwargs = dict(
        orientation='vertical',
        label='Slope (neurons/s/pA)',
        ticks=ti.MultipleLocator(0.5),
        extend='max', extendfrac=0.1)

std_cbar_kwargs = dict(
        orientation='vertical',
        label='Mean $\sigma_{speed}$ (neurons/s)',
        ticks=ti.MultipleLocator(4),
        extend='max', extendfrac=0.1)


def createSweepFig(name=None):
    fig = plt.figure(name, figsize=sweepFigSize)
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    return fig, ax

if (velSweep0):
    # noise_sigma = 0 pA
    fig, ax = createSweepFig("velSlopeSweep0")
    _, ax, cax = EI.plotVelTrial(ps.v[0], slopeVarList, iterList,
            noise_sigmas[0],
            ax=ax,
            xlabel='', xticks=False,
            cbar=False, cbar_kwargs=slope_cbar_kwargs,
            vmin=slope_vmin, vmax=slope_vmax)
    fname = outputDir + "/figure2_slope_sweeps0.png"
    fig.savefig(fname, dpi=300, transparent=True)

    fig, ax = createSweepFig()
    _, ax, cax = EI.plotVelStdSweep(ps.v[0], iterList,
            noise_sigmas[0],
            ax=ax,
            sigmaTitle=False,
            cbar=False, cbar_kwargs=std_cbar_kwargs,
            vmin=std_vmin, vmax=std_vmax)
    fname = outputDir + "/figure2_vel_std_sweeps0.png"
    fig.savefig(fname, dpi=300, transparent=True)



if (velSweep150):
    # noise_sigma = 150 pA
    fig, ax = createSweepFig("velSlopeSweep150")
    _, ax, cax = EI.plotVelTrial(ps.v[1], slopeVarList, iterList,
            noise_sigma=noise_sigmas[1],
            ax=ax,
            xlabel='', xticks=False,
            ylabel='', yticks=False,
            cbar=False, cbar_kwargs=slope_cbar_kwargs,
            vmin=slope_vmin, vmax=slope_vmax)
    fname = outputDir + "/figure2_slope_sweeps150.png"
    fig.savefig(fname, dpi=300, transparent=True)

    fig, ax = createSweepFig()
    _, ax, cax = EI.plotVelStdSweep(ps.v[1], iterList,
            noise_sigmas[1],
            ax=ax,
            ylabel='', yticks=False,
            sigmaTitle=False,
            cbar=False, cbar_kwargs=std_cbar_kwargs,
            vmin=std_vmin, vmax=std_vmax)
    fname = outputDir + "/figure2_vel_std_sweeps150.png"
    fig.savefig(fname, dpi=300, transparent=True)



if (velSweep300):
    # noise_sigma = 300 pA
    fig, ax = createSweepFig("velSlopeSweep300")
    _, ax, cax = EI.plotVelTrial(ps.v[2], slopeVarList, iterList,
            noise_sigma=noise_sigmas[2],
            ax=ax,
            xlabel='', xticks=False,
            ylabel='', yticks=False,
            cbar=True, cbar_kwargs=slope_cbar_kwargs,
            vmin=slope_vmin, vmax=slope_vmax)
    fname = outputDir + "/figure2_slope_sweeps300.png"
    fig.savefig(fname, dpi=300, transparent=True)

    fig, ax = createSweepFig()
    _, ax, cax = EI.plotVelStdSweep(ps.v[2], iterList,
            noise_sigmas[2],
            ax=ax,
            ylabel='', yticks=False,
            sigmaTitle=False,
            cbar=True, cbar_kwargs=std_cbar_kwargs,
            vmin=std_vmin, vmax=std_vmax)
    fname = outputDir + "/figure2_vel_std_sweeps300.png"
    fig.savefig(fname, dpi=300, transparent=True)


if (velLines):
    positions = ((4, 27), (4, 27), (4, 27))
    #positions = ((15, 15), (15, 15), (15, 15))
    fig = plt.figure(figsize=(2.5, velFigsize[1]))
    ax = fig.add_axes(Bbox.from_extents(0.3, velBottom, velRight,
        velTop))
    plotSlopes(ax, ps.v[0], positions[0], noise_sigma=ps.noise_sigmas[0],
            xlabel='', xticks=False,
            ylabel='',
            color='blue')
    fname = outputDir + "/figure2_slope_examples_0.pdf"
    plt.savefig(fname, dpi=300, transparent=True)

    fig = plt.figure(figsize=(2.5, velFigsize[1]))
    ax = fig.add_axes(Bbox.from_extents(0.3, velBottom, velRight,
        velTop))
    plotSlopes(ax, ps.v[1], positions[1], noise_sigma=ps.noise_sigmas[1],
            xlabel='', xticks=False,
            g_ann=False,
            color='green')
    fname = outputDir + "/figure2_slope_examples_1.pdf"
    plt.savefig(fname, dpi=300, transparent=True)

    fig = plt.figure(figsize=(2.5, velFigsize[1]))
    ax = fig.add_axes(Bbox.from_extents(0.3, velBottom, velRight,
        velTop))
    plotSlopes(ax, ps.v[2], positions[2], noise_sigma=ps.noise_sigmas[2],
            ylabel='',
            g_ann=False,
            color='red')
    fname = outputDir + "/figure2_slope_examples_2.pdf"
    plt.savefig(fname, dpi=300, transparent=True)


if (gridness_vs_error):
    fig = plt.figure(figsize=(3.4, 2.5))
    plotGridnessVsFitErr(ps.grids, ps.v, range(NTrials), maxErr=None)
    fig.tight_layout()
    fname = outputDir + "/figure2_gridness_vs_error.pdf"
    plt.savefig(fname, dpi=300, transparent=True)

