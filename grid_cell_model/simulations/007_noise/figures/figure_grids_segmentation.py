#!/usr/bin/env python
'''
Gridness score segmentation plots.
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ti
from matplotlib.transforms import Bbox

from EI_plotting      import sweeps, examples, details, segmentation
from EI_plotting      import aggregate as aggr
from EI_plotting.base import NoiseDataSpaces, getOption, plotStateSignal
from parameters       import JobTrialSpace2D
from data_storage     import DataStorage
from data_storage.sim_models.ei import extractSummedSignals
import plotting.low_level
from plotting.global_defs import prepareLims
from grid_cell_model.analysis import clustering
from submitting import flagparse

from matplotlib import rc
rc('pdf', fonttype=42)
rc('mathtext', default='regular')

plt.rcParams['font.size'] = 11

outputDir = "panels/"

NTrials=1
gridTrialNumList = np.arange(NTrials)
iterList  = ['g_AMPA_total', 'g_GABA_total']

noise_sigmas   = [0, 150, 300]
exampleIdx     = [(5, 15), (5, 15), (5, 15)] # (row, col)
bumpDataRoot   = None
velDataRoot    = None
gridsDataRoot  = 'output_local/even_spacing/grids_vertical'
singleDataRoot = 'output_local/single_neuron'
shape = (31, 31)


parser = flagparse.FlagParser()
parser.add_flag('-c', '--collapsed_sweeps')
parser.add_flag('-d', '--diff_all')
parser.add_flag('-s', '--sweep_segmentation')
parser.add_flag('--diff_sweep')
args = parser.parse_args()


##############################################################################
roots = NoiseDataSpaces.Roots(bumpDataRoot, velDataRoot, gridsDataRoot)
ps    = NoiseDataSpaces(roots, shape, noise_sigmas)


##############################################################################
# Collapsed sweep plots
collapsedFigSize = (3.5, 2.8)
collapsedLeft   = 0.1
collapsedBottom = 0.2
collapsedRight  = 0.95
collapsedTop    = 0.8

gridTypes = ['grids', 'gridnessScore']
NGridTrials = 3

if args.collapsed_sweeps or args.all:
    fig = plt.figure(figsize=collapsedFigSize)
    ax = fig.add_axes(Bbox.from_extents(collapsedLeft, collapsedBottom,
        collapsedRight, collapsedTop))

    data = []
    for ns_idx, _ in enumerate(ps.noise_sigmas):
        d, _, _ = aggr.aggregateType(ps.grids[ns_idx], iterList, gridTypes,
                NTrials, ignoreNaNs=True)
        data.append(d)
    sweeps.plotCollapsedSweeps(ps.noise_sigmas, data,
        ylabel='', yticks=False)
    ax.set_ylim([-0.6, 1.2])
    ax.yaxis.set_major_locator(ti.MultipleLocator(0.4))
    ax.yaxis.set_minor_locator(ti.MultipleLocator(0.2))
    
    fname = outputDir + "/grids_collapsed_sweeps.pdf"
    plt.savefig(fname, dpi=300, transparent=True)
    plt.close()



##############################################################################
# Scatter plot with histograms all in one plot
diffXlim = [-1.5, 1.5]

diffAllFigSize = (6, 6)  # Must be square
diffAllLeft = 0.2
diffAllBottom = 0.15
diffAllScatterSize = 0.55
diffAllScatterRight = diffAllLeft + diffAllScatterSize
diffAllScatterTop   = diffAllBottom + diffAllScatterSize
diffHistSize = 0.2
diffHistOffset = 0.0
diffAllXlim   = prepareLims([-1.5, 1.5])
diffAllXLabel = '$G_{150}$ -  $G_{0}$'
diffAllYLabel = '$G_{300}$ -  $G_{150}$'
histRatio     = 0.2
diffAllBins   = 40

segThresholds = [
        [-np.infty, 0, np.infty],
        [-np.infty, 0, np.infty]]
#segMergeInfo = [
#        [0, 1, 3],
#        [6, 7]]

segMergeInfo = None
gridnessThreshold = 0


if args.diff_all or args.all:
    fig = plt.figure(figsize=diffAllFigSize)

    data = []
    for ns_idx, _ in enumerate(ps.noise_sigmas):
        d, _, _ = aggr.aggregateType(ps.grids[ns_idx], iterList, gridTypes,
                NTrials, ignoreNaNs=True)
        data.append(d)

    # Scatter plot
    ax_scatter = fig.add_axes(Bbox.from_extents(diffAllLeft, diffAllBottom,
        diffAllScatterRight, diffAllScatterTop))
    segmentation.plotDiffScatter(data, ps.noise_sigmas,
            ax=ax_scatter,
            s=15,
            linewidth=0.3,
            xlabel=diffAllXLabel,
            ylabel=diffAllYLabel,
            cmap='Set1',
            filterThreshold=gridnessThreshold, 
            doSegmentation=True, thresholds=segThresholds, mergeInfo=segMergeInfo)
    ax_scatter.set_xlim(diffAllXlim)
    ax_scatter.set_ylim(diffAllXlim)
    ax_scatter.xaxis.set_major_locator(ti.MultipleLocator(1.5))
    ax_scatter.yaxis.set_major_locator(ti.MultipleLocator(1.5))
    ax_scatter.xaxis.set_minor_locator(ti.MultipleLocator(0.1))
    ax_scatter.yaxis.set_minor_locator(ti.MultipleLocator(0.1))
    ax_scatter.xaxis.set_ticks_position('both')
    ax_scatter.yaxis.set_ticks_position('both')

    # 150-0 pA histogram
    diffXHistBottom = diffAllScatterTop + diffHistOffset
    diffXHistTop = diffXHistBottom + diffHistSize

    ax_XHist = fig.add_axes(Bbox.from_extents(diffAllLeft, diffXHistBottom,
        diffAllScatterRight, diffXHistTop))
    segmentation.plotDiffHistograms(data, ps.noise_sigmas,
            ax=ax_XHist, which=0,
            bins=diffAllBins,
            range=diffXlim,
            filterThreshold=gridnessThreshold, 
            xlabel='', ylabel='')
    ax_XHist.set_xlim(diffAllXlim)
    ax_XHist.xaxis.set_major_locator(ti.MultipleLocator(1.5))
    ax_XHist.axis('off')

    # 300 - 150 pA histogram
    diffYHistLeft = diffAllScatterRight + diffHistOffset
    diffYHistRight = diffYHistLeft + diffHistSize
    ax_YHist = fig.add_axes(Bbox.from_extents(diffYHistLeft, diffAllBottom,
        diffYHistRight, diffAllScatterTop))
    segmentation.plotDiffHistograms(data, ps.noise_sigmas,
            ax=ax_YHist, which=1,
            bins=diffAllBins,
            range=diffXlim,
            orientation='horizontal',
            filterThreshold=gridnessThreshold, 
            xlabel='', ylabel='')
    ax_YHist.set_ylim(diffAllXlim)
    ax_YHist.yaxis.set_major_locator(ti.MultipleLocator(1.5))
    ax_YHist.axis('off')


    # Save
    #gs.tight_layout(fig, pad=0.01)
    fname = outputDir + "/grids_diff_all_plots.pdf"
    plt.savefig(fname, dpi=300, transparent=True)
    plt.close()




##############################################################################
# Clustering
sweepFigSize = (3.7, 2.6)
sweepLeft   = 0.08
sweepBottom = 0.2
sweepRight  = 0.8
sweepTop    = 0.85
transparent  = True

if args.sweep_segmentation or args.all:
    fig = plt.figure(figsize=sweepFigSize)
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))

    segmentation.plotSweepSegments(ps.grids, ps.noise_sigmas, iterList,
            segThresholds, gridTypes, NTrials,
            cmap='Set1',
            ax=ax,
            filterThreshold=gridnessThreshold, 
            mergeInfo=segMergeInfo)

    fname = outputDir + "/grids_diff_sweep_segments.pdf"
    plt.savefig(fname, dpi=300, transparent=True)
    plt.close()



##############################################################################
# Parameter sweep of the difference between noise_150 and noise_0
sweepFigSize = (3.7, 2.6)
sweepLeft   = 0.08
sweepBottom = 0.2
sweepRight  = 0.8
sweepTop    = 0.85
transparent  = True
gridDiffText = '$\Delta_{150 - 0}$(Gridness score)'
gridDiff_cbar_kw = dict(
        label       = gridDiffText,
        location    = 'right',
        shrink      = 0.8,
        pad         = -0.05,
        ticks       = ti.MultipleLocator(0.5),
        rasterized  = True)

if args.diff_sweep or args.all:
    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas[0:-1]):
        fig = plt.figure(figsize=sweepFigSize)
        ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
            sweepTop))

        which = ns_idx
        sweeps.plotDiffTrial(ps.grids, iterList, which, NGridTrials, gridTypes,
                ax=ax,
                ignoreNaNs=True,
                r=exampleIdx[0][0], c=exampleIdx[0][1],
                cbar=True, cbar_kw=gridDiff_cbar_kw,
                symmetricLimits=True,
                cmap='RdBu_r')

        fname = outputDir + "/grids_diff_sweep{0}.pdf"
        plt.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)
        plt.close()


