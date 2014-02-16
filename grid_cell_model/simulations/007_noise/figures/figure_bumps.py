#!/usr/bin/env python
#
'''
Noise publication Figure 3.
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ti
from matplotlib.gridspec   import GridSpec
from matplotlib.colorbar   import make_axes
from matplotlib.transforms import Bbox

from EI_plotting          import sweeps, details, examples, rasters
from EI_plotting          import aggregate as aggr
from plotting.global_defs import globalAxesSettings, prepareLims
from EI_plotting.base     import NoiseDataSpaces, getNoiseDataSpaces
from parameters           import JobTrialSpace2D
from analysis             import clustering
from submitting import flagparse

from matplotlib import rc
rc('pdf', fonttype=42)
rc('mathtext', default='regular')

plt.rcParams['font.size'] = 11

outputDir = "panels/"

NTrials=10
iterList  = ['g_AMPA_total', 'g_GABA_total']

noise_sigmas = [0, 150, 300]
exampleIdx   = [(0, 0), (0, 0), (0, 0)] # (row, col)
gridsDataRoot    = 'output_local/even_spacing/grids'
bumpDataRoot     = 'output_local/even_spacing/gamma_bump'
velDataRoot      = 'output_local/even_spacing/velocity'
constPosDataRoot = 'output_local/even_spacing/const_position'
shape = (31, 31)

parser = flagparse.FlagParser()
parser.add_flag('--bumpSweep')
parser.add_flag('--bumpDriftSweep')
parser.add_flag('--bumpDiffAtTSweep')
parser.add_flag('--bumpDiffResetSweep')
parser.add_flag('--bumpExamples')
parser.add_flag('--velExamples')
parser.add_flag('--velSweep')
parser.add_flag('--gridness_vs_error')
parser.add_flag('--detailed_noise')
parser.add_flag('--rastersFlag')
parser.add_flag('--rates')
args = parser.parse_args()

###############################################################################
roots = NoiseDataSpaces.Roots(bumpDataRoot, velDataRoot, gridsDataRoot,
        constPos=constPosDataRoot)
ps    = NoiseDataSpaces(roots, shape, noise_sigmas)
constPosPS = getNoiseDataSpaces(roots.constPos, noise_sigmas, shape)

exW = 4
exH = 2
exMargin = 0.075
exWspace=0.2
exHspace=0.15

sweepFigSize = (3.7, 2.6)
sweepLeft   = 0.08
sweepBottom = 0.2
sweepRight  = 0.8
sweepTop    = 0.85
transparent  = True


##############################################################################
# Bump sigma sweeps
sigmaBumpText = '$\sigma_{bump}^{-1}\ (neurons^{-1})$'
bumpTStart = 500.0
bumpNTrials = 5
bump_vmin = 0
bump_vmax = 0.421
bump_cbar_kw = dict(
        label       = sigmaBumpText,
        location    = 'right',
        shrink      = 0.8,
        pad         = -0.05,
        ticks       = ti.MultipleLocator(0.2),
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

if args.bumpSweep or args.all:
    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        fig = plt.figure(figsize=sweepFigSize)
        ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
            sweepTop))
        kw = dict(cbar=False)
        if ns_idx != 0:
            kw['ylabel'] = ''
            kw['yticks'] = False
        if ns_idx == 2:
            kw['cbar'] = True
        data = aggr.AggregateBumpReciprocal(ps.bumpGamma[ns_idx], iterList,
                bumpNTrials, tStart=bumpTStart)
        _, _, cax = sweeps.plotSweep(data,
                noise_sigma=noise_sigma,
                ax=ax,
                cbar_kw=bump_cbar_kw,
                vmin=bump_vmin, vmax=bump_vmax,
                annotations=ann, **kw)
        fname = outputDir + "/bumps_sweeps{0}.pdf"
        fig.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)
        plt.close()

##############################################################################
# Bump drift at a specified time
bumpDriftText = 'Average bump drift\n(neurons)'
bumpDriftT = 9e3 # ms
drift_vmin = 0
drift_vmax = 24
bump_drift_cbar_kw = dict(
        label       = bumpDriftText,
        location    = 'right',
        shrink      = 0.8,
        pad         = -0.05,
        ticks       = ti.MultipleLocator(10),
        rasterized  = True)


if args.bumpDriftSweep or args.all:
    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        fig = plt.figure(figsize=sweepFigSize)
        ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
            sweepTop))
        kw = dict(cbar=False)
        if ns_idx != 0:
            kw['ylabel'] = ''
            kw['yticks'] = False
        if ns_idx == 2:
            kw['cbar'] = True
        data = aggr.BumpDriftAtTime(bumpDriftT, 
                ps.bumpGamma[ns_idx],
                iterList,
                bumpNTrials,
                tStart=bumpTStart)
        _, _, cax = sweeps.plotSweep(data, noise_sigma=noise_sigma,
                ax=ax,
                cbar_kw=bump_drift_cbar_kw,
                vmin=drift_vmin, vmax=drift_vmax,
                **kw)
        fname = outputDir + "/bumps_drift_at_time_sweeps{0}.pdf"
        fig.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)
        plt.close()

##############################################################################
# Distance from init position
bumpDiffText = 'Distance from init position\n(neurons)'
bumpDiffT = 0.75e3 # ms
bumpDiff_vmin = 0
bumpDiff_vmax = 20
diffStartPos = [17, 15]
bumpDiff_cbar_kw = dict(
        label       = bumpDiffText,
        location    = 'right',
        shrink      = 0.8,
        pad         = -0.05,
        ticks       = ti.MultipleLocator(10),
        rasterized  = True)


if args.bumpDiffAtTSweep or args.all:
    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        fig = plt.figure(figsize=sweepFigSize)
        ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
            sweepTop))
        kw = dict(cbar=False)
        if ns_idx != 0:
            kw['ylabel'] = ''
            kw['yticks'] = False
        if ns_idx == 2:
            kw['cbar'] = True
        data = aggr.BumpDifferenceAtTime(diffStartPos, bumpDiffT,
                ps.bumpGamma[ns_idx],
                iterList,
                bumpNTrials)
        _, _, cax = sweeps.plotSweep(data, noise_sigma=noise_sigma,
                ax=ax,
                cbar_kw=bumpDiff_cbar_kw,
                vmin=bumpDiff_vmin, vmax=bumpDiff_vmax,
                **kw)
        fname = outputDir + "/bumps_difference_at_time_sweeps{0}.pdf"
        fig.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)
        plt.close()

##############################################################################
# Average distance from position enforced by place cells, from theta_start_t
# until the end of the simultion
bumpResetText = 'Distance from reset\nposition (neurons)'
bumpResetTStart = 0.5e3 # ms
bumpReset_vmin = 0
bumpReset_vmax = 20
bumpResetStartPos = [17, 15]
constPosNTrials = 5
bumpReset_cbar_kw = dict(
        label       = bumpResetText,
        location    = 'right',
        shrink      = 0.8,
        pad         = -0.05,
        ticks       = ti.MultipleLocator(5),
        rasterized  = True)


if args.bumpDiffResetSweep or args.all:
    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        fig = plt.figure(figsize=sweepFigSize)
        ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
            sweepTop))
        kw = dict(cbar=False)
        if ns_idx != 0:
            kw['ylabel'] = ''
            kw['yticks'] = False
        if ns_idx == 2:
            kw['cbar'] = True
        data = aggr.BumpAvgDifferenceFromPos(bumpResetStartPos,
                constPosPS[ns_idx],
                iterList,
                constPosNTrials,
                tstart=bumpResetTStart)
        _, _, cax = sweeps.plotSweep(data, noise_sigma=noise_sigma,
                ax=ax,
                cbar_kw=bumpReset_cbar_kw,
                vmin=bumpReset_vmin, vmax=bumpReset_vmax,
                **kw)
        fname = outputDir + "/bumps_avg_difference_reset_sweeps{0}.pdf"
        fig.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)
        plt.close()


##############################################################################
# Bump examples
exampleEFName = outputDir + "/bumps_examples_E_{0}pA_{1}.pdf"
exampleIFName = outputDir + "/bumps_examples_I_{0}pA_{1}.pdf"
bumpExampleTypes = ['bump_full']
bumpTrialNum = 0
exTransparent = True
exampleFigSize = (0.8, 0.8)
exampleLeft   = 0.01
exampleBottom = 0.01
exampleRight  = 0.99
exampleTop    = 0.82

if args.bumpExamples or args.all:
    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        for idx, rc in enumerate(exampleRC):
            for EIType in ['E', 'I']:
                if EIType == 'E':
                    fnameTemplate =exampleEFName
                    types = bumpExampleTypes + ['rateMap_e']
                else:
                    fnameTemplate =exampleIFName
                    types = bumpExampleTypes + ['rateMap_i']
                fname = fnameTemplate.format(noise_sigma, idx)
                plt.figure(figsize=exampleFigSize)
                gs = examples.plotOneBumpExample(ps.bumpGamma[ns_idx], rc, iterList,
                        types,
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

if args.velSweep or args.all:
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


#if args.gridness_vs_error or args.all:
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
if args.detailed_noise or args.all:
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


if args.rastersFlag or args.all:
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


if args.rates or args.all:
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



