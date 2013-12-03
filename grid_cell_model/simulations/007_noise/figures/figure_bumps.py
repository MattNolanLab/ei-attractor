#!/usr/bin/env python
#
#   figure3.py
#
#   Noise publication Figure 3.
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

from EI_plotting                import sweeps, details, examples, rasters
from plotting.global_defs       import globalAxesSettings
from figures_shared             import NoiseDataSpaces
from parameters                 import JobTrialSpace2D

import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)

from matplotlib import rc
rc('pdf', fonttype=42)
rc('mathtext', default='regular')

plt.rcParams['font.size'] = 11

outputDir = "panels/"

NTrials=10
iterList  = ['g_AMPA_total', 'g_GABA_total']

noise_sigmas = [0, 150, 300]
exampleIdx   = [(0, 0), (0, 0), (0, 0)] # (row, col)
gridsDataRoot= 'output_local/even_spacing/grids'
bumpDataRoot= 'output_local/even_spacing/gamma_bump'
velDataRoot = 'output_local/even_spacing/velocity'
shape = (31, 31)

bumpSweep         = 0
bumpExamples      = 1
velExamples       = 0
velSweep          = 0
gridness_vs_error = 0
detailed_noise    = 1
rastersFlag       = 0
rates             = 0

##############################################################################

def plotBumpExample(sp, rc, iterList, EIType, **kw):
    #keyword
    wspace    = kw.pop('wspace', 0)
    hspace    = kw.pop('hspace', 0)
    gsCoords  = kw.pop('exGsCoords', (0, 0, 1, 1))

    r, c = rc[0], rc[1] 
    spaceRect = [c, r, c, r]
    return examples.drawBumpExamples(sp, spaceRect, iterList, gsCoords, EIType,
            xlabel=False, ylabel=False,
            xlabel2=False, ylabel2=False,
            fontsize='x-small',
            rateYPos=1.05, rateXPos=0.98,
            **kw)


###############################################################################
roots = NoiseDataSpaces.Roots(bumpDataRoot, velDataRoot, gridsDataRoot)
ps    = NoiseDataSpaces(roots, shape, noise_sigmas)

exW = 4
exH = 2
exMargin = 0.075
exWspace=0.2
exHspace=0.15

sweepFigSize = (2.4, 2.9)
sweepLeft   = 0.2
sweepBottom = 0.1
sweepRight  = 0.87
sweepTop    = 0.9


##############################################################################
sigmaBumpText = '$\sigma_{bump}^{-1}\ (neurons^{-1})$'
sigmaVarList = ['bump_e', 'sigma']
bumpTrialNumList = np.arange(5)
bump_vmin = 0
bump_vmax = 0.58
bump_cbar_kw = dict(
        label       = sigmaBumpText,
        location    = 'bottom',
        shrink      = 0.8,
        pad         = 0.18,
        ticks       = ti.MultipleLocator(0.25),
        extend      = 'max',
        extendfrac  = 0.1,
        rasterized  = True)

exampleRC = ( (5, 15), (15, 5) )

ann0 = dict(
        txt='b',
        rc=exampleRC[0],
        xytext_offset=(1.5, 0.5),
        color='white')
ann1 = dict(
        txt='a',
        rc=exampleRC[1],
        xytext_offset=(1.2, 1.1),
        color='black')
ann = [ann0, ann1]

if (bumpSweep):
    # noise_sigma = 0 pA
    fig = plt.figure("sweeps0", figsize=sweepFigSize)
    exRows = [28, 15]
    exCols = [3, 15]
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    sweeps.plotBumpSigmaTrial(ps.bumpGamma[0], sigmaVarList, iterList,
            noise_sigma=ps.noise_sigmas[0],
            ax=ax,
            trialNumList=bumpTrialNumList,
            cbar=False, cbar_kw=bump_cbar_kw,
            vmin=bump_vmin, vmax=bump_vmax,
            annotations=ann)
    fname = outputDir + "/bumps_sweeps0.pdf"
    fig.savefig(fname, dpi=300, transparent=True)

    # noise_sigma = 150 pA
    for a in ann:
        a['color'] = 'black'
    fig = plt.figure("sweeps150", figsize=sweepFigSize)
    exRows = [8, 2]
    exCols = [10, 9]
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    sweeps.plotBumpSigmaTrial(ps.bumpGamma[1], sigmaVarList, iterList,
            noise_sigma=noise_sigmas[1],
            ax=ax,
            ylabel='', yticks=False,
            trialNumList=bumpTrialNumList,
            cbar=False, cbar_kw=bump_cbar_kw,
            vmin=bump_vmin, vmax=bump_vmax,
            annotations=ann)
    fname = outputDir + "/bumps_sweeps150.pdf"
    fig.savefig(fname, dpi=300, transparent=True)

    # noise_sigma = 300 pA
    fig = plt.figure("sweeps300", figsize=sweepFigSize)
    exRows = [16, 15]
    exCols = [6, 23]
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    ax.set_clip_on(False)
    sweeps.plotBumpSigmaTrial(ps.bumpGamma[2], sigmaVarList, iterList,
            noise_sigma=noise_sigmas[2],
            ax=ax,
            trialNumList=bumpTrialNumList,
            cbar=True, cbar_kw=bump_cbar_kw,
            vmin=bump_vmin, vmax=bump_vmax,
            annotations=ann)
    fname = outputDir + "/bumps_sweeps300.pdf"
    fig.savefig(fname, dpi=300, transparent=True)


##############################################################################
# Bump examples
exampleEFName = outputDir + "/bumps_examples_E_{0}pA_{1}.pdf"
exampleIFName = outputDir + "/bumps_examples_I_{0}pA_{1}.pdf"
bumpTrialNum = 0
exTransparent = True
exampleFigSize = (0.8, 0.8)
exampleLeft   = 0.01
exampleBottom = 0.01
exampleRight  = 0.99
exampleTop    = 0.82

if (bumpExamples):
    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        for idx, rc in enumerate(exampleRC):
            for EIType in ['E', 'I']:
                if EIType == 'E':
                    fnameTemplate =exampleEFName
                else:
                    fnameTemplate =exampleIFName
                fname = fnameTemplate.format(noise_sigma, idx)
                plt.figure(figsize=exampleFigSize)
                gs = plotBumpExample(ps.bumpGamma[ns_idx], rc, iterList,
                        EIType,
                        exIdx=exampleIdx[ns_idx],
                        trialNum=bumpTrialNum)
                gs.update(left=exampleLeft, bottom=exampleBottom,
                        right=exampleRight, top=exampleTop)
                plt.savefig(fname, dpi=300, transparent=exTransparent)
                plt.close()



###############################################################################

std_vmin = 0
std_vmax = 14

std_cbar_kw= dict(
        orientation='vertical',
        label='Mean $\sigma_{speed}$ (neurons/s)',
        shrink = 0.8,
        pad = 0.05,
        ticks=ti.MultipleLocator(4),
        extend='max', extendfrac=0.1)


def createSweepFig(name=None):
    sweepFigSize = (2.6, 1.9)
    sweepLeft   = 0.15
    sweepBottom = 0.2
    sweepRight  = 0.87
    sweepTop    = 0.85
    fig = plt.figure(name, figsize=sweepFigSize)
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    return fig, ax

if (velSweep):
    # noise_sigma = 0 pA
    fig, ax = createSweepFig()
    _, ax, cax = sweeps.plotVelStdSweep(ps.v[0], iterList,
            noise_sigmas[0],
            ax=ax,
            sigmaTitle=False,
            cbar=False, cbar_kw=std_cbar_kw,
            vmin=std_vmin, vmax=std_vmax)
    fname = outputDir + "/bumps_vel_std_sweeps0.pdf"
    fig.savefig(fname, dpi=300, transparent=True)


    # noise_sigma = 150 pA
    fig, ax = createSweepFig()
    _, ax, cax = sweeps.plotVelStdSweep(ps.v[1], iterList,
            noise_sigmas[1],
            ax=ax,
            ylabel='', yticks=False,
            sigmaTitle=False,
            cbar=False, cbar_kw=std_cbar_kw,
            vmin=std_vmin, vmax=std_vmax)
    fname = outputDir + "/bumps_vel_std_sweeps150.pdf"
    fig.savefig(fname, dpi=300, transparent=True)


    # noise_sigma = 300 pA
    fig, ax = createSweepFig()
    _, ax, cax = sweeps.plotVelStdSweep(ps.v[2], iterList,
            noise_sigmas[2],
            ax=ax,
            ylabel='', yticks=False,
            sigmaTitle=False,
            cbar=True, cbar_kw=std_cbar_kw,
            vmin=std_vmin, vmax=std_vmax)
    fname = outputDir + "/bumps_vel_std_sweeps300.pdf"
    fig.savefig(fname, dpi=300, transparent=True)


#if (gridness_vs_error):
#    fig = plt.figure(figsize=(3.4, 2.5))
#    plotGridnessVsFitErr(ps.grids, ps.v, range(NTrials), maxErr=None)
#    fig.tight_layout()
#    fname = outputDir + "/bumps_gridness_vs_error.pdf"
#    plt.savefig(fname, dpi=300, transparent=True)


##############################################################################
# Detailed noise plots
EI13Root  = 'output_local/detailed_noise/gamma_bump/EI-1_3'
EI31Root  = 'output_local/detailed_noise/gamma_bump/EI-3_1'
detailedShape = (31, 9)

EI13PS = JobTrialSpace2D(detailedShape, EI13Root)
EI31PS = JobTrialSpace2D(detailedShape, EI31Root)
detailedNTrials = 5


detailFigSize = (3.8, 2.6)
detailLeft   = 0.18
detailBottom = 0.26
detailRight  = 0.95
detailTop    = 0.95
if (detailed_noise):
    ylabelPos = -0.17

    types = ('bump', 'sigma')
    fig = plt.figure(figsize=detailFigSize)
    ax = fig.add_axes(Bbox.from_extents(detailLeft, detailBottom, detailRight,
        detailTop))
    _, p13, l13 = details.plotDetailedNoise(EI13PS, detailedNTrials, types, ax=ax,
            ylabel=sigmaBumpText, ylabelPos=ylabelPos,
            color='red', markerfacecolor='red', zorder=10)
    _, p31, l31 = details.plotDetailedNoise(EI31PS, detailedNTrials, types, ax=ax,
            ylabelPos=ylabelPos,
            color='#505050')
    #ax.set_yscale("log")
    #ax.set_ylim([1.5, 300])
    leg = ['a',  'b']
    l = ax.legend([p31, p13], leg, loc=(0.85, 0.1), fontsize='small', frameon=False,
            numpoints=1, handletextpad=0.05)
    plt.setp(l.get_title(), fontsize='small')

    fname = outputDir + "/bumps_detailed_noise_sigma.pdf"
    plt.savefig(fname, dpi=300, transparent=True)
    plt.close()


##############################################################################
rasterRC      = [(5, 15), (5, 15), (5, 15)] # (row, col)
tLimits = [2e3, 2.25e3] # ms

rasterFigSize = (3, 1.9)
transparent   = True
rasterLeft    = 0.28
rasterBottom  = 0.1
rasterRight   = 0.95
rasterTop     = 0.8

ylabelPos   = -0.35


if (rastersFlag):
    # noise_sigma = 0 pA
    fig = plt.figure("rasters0", figsize=rasterFigSize)
    ax = fig.add_axes(Bbox.from_extents(rasterLeft, rasterBottom, rasterRight,
        rasterTop))
    rasters.EIRaster(ps.bumpGamma[0], 
            noise_sigma=ps.noise_sigmas[0],
            spaceType='bump',
            r=rasterRC[0][0], c=rasterRC[0][1],
            ylabelPos=ylabelPos,
            tLimits=tLimits,
            ann_EI=True)
    fname = outputDir + "/bumps_raster0.png"
    fig.savefig(fname, dpi=300, transparent=transparent)
    plt.close()
        

    # noise_sigma = 150 pA
    fig = plt.figure("rasters150", figsize=rasterFigSize)
    ax = fig.add_axes(Bbox.from_extents(rasterLeft, rasterBottom, rasterRight,
        rasterTop))
    rasters.EIRaster(ps.bumpGamma[1], 
            noise_sigma=ps.noise_sigmas[1],
            spaceType='bump',
            r=rasterRC[1][0], c=rasterRC[1][1],
            ylabelPos=ylabelPos,
            tLimits=tLimits,
            ylabel='', yticks=False)
    fname = outputDir + "/bumps_raster150.png"
    fig.savefig(fname, dpi=300, transparent=transparent)
    plt.close()
        

    # noise_sigma = 300 pA
    fig = plt.figure("rasters300", figsize=rasterFigSize)
    ax = fig.add_axes(Bbox.from_extents(rasterLeft, rasterBottom, rasterRight,
        rasterTop))
    rasters.EIRaster(ps.bumpGamma[2], 
            noise_sigma=ps.noise_sigmas[2],
            spaceType='bump',
            r=rasterRC[2][0], c=rasterRC[2][1],
            ylabelPos=ylabelPos,
            tLimits=tLimits,
            ylabel='', yticks=False)
    fname = outputDir + "/bumps_raster300.png"
    fig.savefig(fname, dpi=300, transparent=transparent)
    plt.close()
        

##############################################################################
rateFigSize   = (rasterFigSize[0], 0.65)
rateLeft    = rasterLeft
rateBottom  = 0.2
rateRight   = rasterRight
rateTop     = 0.9


if (rates):
    for idx, noise_sigma in enumerate(ps.noise_sigmas):
        # E cells
        fig = plt.figure(figsize=rateFigSize)
        ax = fig.add_axes(Bbox.from_extents(rateLeft, rateBottom, rateRight,
            rateTop))
        kw = {}
        if (idx != 0):
            kw['ylabel'] = ''

        rasters.plotAvgFiringRate(ps.bumpGamma[idx],
                spaceType='bump',
                noise_sigma=ps.noise_sigmas[idx],
                popType='E',
                r=rasterRC[idx][0], c=rasterRC[idx][1],
                ylabelPos=ylabelPos,
                color='red',
                tLimits=tLimits,
                ax=ax, **kw)
        fname = outputDir + "/bumps_rate_e{0}.pdf".format(noise_sigma)
        fig.savefig(fname, dpi=300, transparent=transparent)
        plt.close()

        # I cells
        fig = plt.figure(figsize=rateFigSize)
        ax = fig.add_axes(Bbox.from_extents(rateLeft, rateBottom, rateRight,
            rateTop))
        kw = {}
        if (idx != 0):
            kw['ylabel'] = ''

        rasters.plotAvgFiringRate(ps.bumpGamma[idx],
                spaceType='bump',
                noise_sigma=ps.noise_sigmas[idx],
                popType='I', 
                r=rasterRC[idx][0], c=rasterRC[idx][1],
                ylabelPos=ylabelPos,
                color='blue',
                tLimits=tLimits,
                ax=ax, **kw)
        fname = outputDir + "/bumps_rate_i{0}.pdf".format(noise_sigma)
        fig.savefig(fname, dpi=300, transparent=transparent)
        plt.close()



