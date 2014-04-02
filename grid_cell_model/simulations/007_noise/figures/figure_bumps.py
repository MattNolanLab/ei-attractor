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

import default_settings as ds
from EI_plotting          import sweeps, details, examples, scatter
from EI_plotting          import aggregate as aggr
from plotting.global_defs import prepareLims
from plotting.low_level   import zeroLines
from parameters           import JobTrialSpace2D
from submitting import flagparse

outputDir = ds.figOutputDir
exampleIdx   = [(0, 0), (0, 0), (0, 0)] # (row, col)

parser = flagparse.FlagParser()
parser.add_flag('--bumpSweep')
parser.add_flag('--bumpDriftSweep')
parser.add_flag('--bumpDiffAtInitSweep')
parser.add_flag('--bumpDiffResetSweep')
parser.add_flag('--bumpExamples')
parser.add_flag('--velExamples')
parser.add_flag('--velSweep')
parser.add_flag('--gridness_vs_error')
parser.add_flag('--detailed_noise')
parser.add_flag('--rastersFlag')
parser.add_flag('--rates')
parser.add_flag('--scatter_diff_fracTotal_grids')
parser.add_flag('--scatter_grids_fracTotal')
parser.add_flag('--fracTotalSweepAnn')
parser.add_flag('--isBump')
parser.add_argument('--expScatter', action="store_true")
args = parser.parse_args()


###############################################################################
ps = ds.getDefaultParamSpaces()

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
        data = aggr.AggregateBumpReciprocal(ps.bumpGamma[ns_idx], ds.iterList,
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
bumpDriftTStart = 1e3 #ms
bumpDriftT = 9e3 # ms
drift_vmin = 0
drift_vmax = 20
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
                ds.iterList,
                bumpNTrials,
                tStart=bumpDriftTStart)
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
bumpDiffText = 'Distance from init\nposition (neurons)'
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


if args.bumpDiffAtInitSweep or args.all:
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
                ds.iterList,
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
                ds.iterList,
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
                gs = examples.plotOneBumpExample(ps.bumpGamma[ns_idx], rc, ds.iterList,
                        types,
                        exIdx=exampleIdx[ns_idx],
                        trialNum=bumpTrialNum)
                gs.update(left=exampleLeft, bottom=exampleBottom,
                        right=exampleRight, top=exampleTop)
                plt.savefig(fname, dpi=300, transparent=exTransparent)
                plt.close()



###############################################################################

std_vmin = 0
std_vmax = 10.421

std_cbar_kw= dict(
        orientation='vertical',
        label='Mean $\sigma_{speed}$ (neurons/s)',
        shrink = 0.8,
        pad = 0.05,
        ticks=ti.MultipleLocator(5),
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
    _, ax, cax = sweeps.plotVelStdSweep(ps.v[0], ds.iterList,
            ps.noise_sigmas[0],
            ax=ax,
            sigmaTitle=False,
            cbar=False, cbar_kw=std_cbar_kw,
            vmin=std_vmin, vmax=std_vmax)
    fname = outputDir + "/bumps_vel_std_sweeps0.pdf"
    fig.savefig(fname, dpi=300, transparent=True)


    # noise_sigma = 150 pA
    fig, ax = createSweepFig()
    _, ax, cax = sweeps.plotVelStdSweep(ps.v[1], ds.iterList,
            ps.noise_sigmas[1],
            ax=ax,
            ylabel='', yticks=False,
            sigmaTitle=False,
            cbar=False, cbar_kw=std_cbar_kw,
            vmin=std_vmin, vmax=std_vmax)
    fname = outputDir + "/bumps_vel_std_sweeps150.pdf"
    fig.savefig(fname, dpi=300, transparent=True)


    # noise_sigma = 300 pA
    fig, ax = createSweepFig()
    _, ax, cax = sweeps.plotVelStdSweep(ps.v[2], ds.iterList,
            ps.noise_sigmas[2],
            ax=ax,
            ylabel='', yticks=False,
            sigmaTitle=False,
            cbar=True, cbar_kw=std_cbar_kw,
            vmin=std_vmin, vmax=std_vmax)
    fname = outputDir + "/bumps_vel_std_sweeps300.pdf"
    fig.savefig(fname, dpi=300, transparent=True)


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
# Correlate (difference between) isBump and gridness score.
corrDiffFigsize = (4, 5)
corrDiffXLabel = "$\Delta$ P(bumps)"
corrDiffYLabel = '$\Delta$ Gridness score'

def setCorrAxes(ax):
    ax.set_xlim(prepareLims([-1, 1]))
    ax.set_ylim(prepareLims([-1.5, 1.5]))
    ax.xaxis.set_major_locator(ti.MultipleLocator(0.5))
    ax.yaxis.set_major_locator(ti.MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(ti.MultipleLocator(0.25))
    ax.yaxis.set_minor_locator(ti.MultipleLocator(0.25))
    zeroLines(ax)


if args.scatter_diff_fracTotal_grids or args.all:
    fig = plt.Figure(corrDiffFigsize)
    ax = fig.add_subplot(111)

    isBumpData = []
    gridData = []
    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        isBumpData.append(aggr.IsBump(ps.bumpGamma[ns_idx], ds.iterList,
            ignoreNaNs=True))
        gridData.append(aggr.GridnessScore(ps.grids[ns_idx], ds.iterList,
            ignoreNaNs=True))

    which = 0
    scatterPlot = scatter.DiffScatterPlot(
            isBumpData, gridData, None, None, None, None, None, which,
            s=15,
            linewidth=0.3,
            edgecolor='white',
            xlabel = corrDiffXLabel,
            ylabel = corrDiffYLabel,
            sigmaTitle=False,
            ignoreNaNs=True,
            cmap='Set1',
            ax=ax)

    scatterPlot.plot()
    #ax.set_title('Difference\n $\sigma_{noise} = 150\ -\ \sigma_{noise} = 0$ pA')
    setCorrAxes(ax)

    fig.tight_layout()
    fname = outputDir + "/bumps_scatter_diff_bumpFracTotal_gscore.pdf"
    fig.savefig(fname, dpi=300, transparent=transparent)
    plt.close()
    

##############################################################################
# Correlate P(bump) vs gridness score
xlabel = 'P(bumps)'
ylabel = 'Gridness score'

scatterAllFigSize = (5.8, 3.2)
scatterAllLeft   = 0.05
scatterAllBottom = 0.05
scatterAllRight  = 0.95
scatterAllTop    = 0.9

if args.scatter_grids_fracTotal or args.all:
    fig = plt.figure(figsize=scatterAllFigSize)
    ax = fig.gca()

    scatterColors = ['blue', 'green', 'red']

    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        isBumpData = aggr.IsBump(ps.bumpGamma[ns_idx], ds.iterList,
                                 ignoreNaNs=True)
        gridData = aggr.GridnessScore(ps.grids[ns_idx], ds.iterList,
                                      ignoreNaNs=True)

        scatterPlot = scatter.ScatterPlot(
                isBumpData, gridData, None, None, None, None, None,
                c=scatterColors[ns_idx],
                s=15,
                linewidth=0.3,
                xlabel=xlabel,
                ylabel=ylabel,
                ax=ax)
        scatterPlot.plot()

    ax.xaxis.set_major_locator(ti.MultipleLocator(0.2))
    ax.yaxis.set_major_locator(ti.MultipleLocator(0.5))
    leg = ['0', '150', '300']
    l = ax.legend(leg, loc=(-0.28, 0.8), fontsize='small', frameon=False,
            numpoints=1, title='$\sigma$ (pA)')
    plt.setp(l.get_title(), size='small')
    ax.set_ylabel(ax.get_ylabel(), y=0., ha='left')

    fig.tight_layout(rect=[scatterAllLeft, scatterAllBottom, scatterAllRight,
                           scatterAllTop])
    fname = outputDir + "/bumps_scatter_grids_vs_bumpFracTotal.pdf"
    fig.savefig(fname, dpi=300, transparent=True)



##############################################################################
# Bump formation sweeps
bump_vmin = 0
bump_vmax = 1.
fracTotalText = 'P(bumps)'

ann0_0 = dict(
        txt='a',
        rc=(15, 5),
        xytext_offset=(.5, 2),
        color='white')
ann0_1 = dict(
        txt='b',
        rc=(5, 15),
        xytext_offset=(1.5, 1),
        color='white')

ann150_0 = dict(
        txt='c',
        rc=(5, 15),
        xytext_offset=(1.5, 1),
        color='white')

ann300_0 = dict(
        txt='d',
        rc=(15, 5),
        xytext_offset=(1.5, 1),
        color='white')
ann300_1 = dict(
        txt='e',
        rc=(5, 15),
        xytext_offset=(1.5, 1),
        color='white')


ann0   = [  ann0_0,   ann0_1]
ann150 = [ann150_0]
ann300 = [ann300_0, ann300_1]
ann = [ann0, ann150, ann300]

sweepAnn_cbar_kw = dict(
        label       = fracTotalText,
        location    = 'left',
        shrink      = 0.8,
        pad         = 0.25,
        ticks       = ti.MultipleLocator(0.5),
        rasterized  = True)
if args.fracTotalSweepAnn or args.all:
    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        fig, ax = ds.getDefaultSweepFig(scale=0.8, colorBarPos='left')
        kw = dict(cbar=False)
        if ns_idx != 0:
            kw['ylabel'] = ''
            kw['yticks'] = False
        if ns_idx == 0:
            kw['cbar'] = True
        data = aggr.IsBump(ps.bumpGamma[ns_idx], ds.iterList, ignoreNaNs=True)
        _, _, cax = sweeps.plotSweep(data,
                noise_sigma=noise_sigma,
                xlabel='', xticks=False,
                ax=ax,
                cbar_kw=sweepAnn_cbar_kw,
                vmin=bump_vmin, vmax=bump_vmax,
                annotations=ann[ns_idx],
                **kw)
        fname = outputDir + "/bumps_mainFig_isBumpFracTotal_sweeps_annotated{0}.pdf"
        fig.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)
        plt.close()


##############################################################################
# Bump formation, thresholded
bumpThreshold = 0.95
if args.isBump or args.all:
    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        fig, ax = ds.getDefaultSweepFig(scale=.8, colorBarPos='left' )
        kw = dict(cbar=False, sigmaTitle=False)
        if ns_idx != 0:
            kw['ylabel'] = ''
            kw['yticks'] = False
        data = aggr.BumpFormationFilter(bumpThreshold, ps.bumpGamma[ns_idx],
                ds.iterList, ignoreNaNs=True)
        _, _, cax = sweeps.plotSweep(data,
                noise_sigma=noise_sigma,
                ax=ax,
                cbar_kw=sweepAnn_cbar_kw,
                vmin=bump_vmin, vmax=bump_vmax,
                annotations=ann[ns_idx],
                **kw)
        fname = outputDir + "/bumps_mainFig_isBump_sweeps{0}.pdf"
        fig.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)
        plt.close()

