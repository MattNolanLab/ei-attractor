#!/usr/bin/env python
#
#   figure_grids_segmentation.py
#
#   Gridness score segmentation plots
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
from matplotlib.transforms import Bbox

from EI_plotting      import sweeps, examples, details, segmentation
from EI_plotting      import aggregate as aggr
from EI_plotting.base import NoiseDataSpaces, getOption, plotStateSignal
from parameters       import JobTrialSpace2D
from data_storage     import DataStorage
from data_storage.sim_models.ei import extractSummedSignals
import plotting.low_level
from plotting.global_defs import prepareLims
from analysis         import clustering

import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)

from matplotlib import rc
rc('pdf', fonttype=42)
rc('mathtext', default='regular')

plt.rcParams['font.size'] = 11

outputDir = "panels/"

NTrials=3
gridTrialNumList = np.arange(NTrials)
iterList  = ['g_AMPA_total', 'g_GABA_total']

noise_sigmas   = [0, 150, 300]
exampleIdx     = [(1, 22), (1, 22), (1, 22)] # (row, col)
bumpDataRoot   = None
velDataRoot    = None
gridsDataRoot  = 'output_local/even_spacing/grids'
singleDataRoot = 'output_local/single_neuron'
shape = (31, 31)

collapsed_sweeps   = 1
diff_distrib       = 1
diff_scatter       = 1
diff_all           = 1
sweep_segmentation = 1

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

if collapsed_sweeps:
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
# Distribution of differences
diffFigSize = (5, 3)
diffLeft   = 0.15
diffBottom = 0.2
diffRight  = 0.9
diffTop    = 0.9
diffXlim = [-1.5, 1.5]
diffLegend = ['$G_{150}$ -  $G_{0}$', '$G_{300}$ -  $G_{150}$']


if (diff_distrib):

    data = []
    for ns_idx, _ in enumerate(ps.noise_sigmas):
        d, _, _ = aggr.aggregateType(ps.grids[ns_idx], iterList, gridTypes,
                NTrials, ignoreNaNs=True)
        data.append(d)

    for ns_idx in xrange(len(ps.noise_sigmas) - 1):
        fig = plt.figure(figsize=diffFigSize)
        ax = fig.add_axes(Bbox.from_extents(diffLeft, diffBottom, diffRight,
            diffTop))

        segmentation.plotDiffHistograms(data, ps.noise_sigmas,
                ax=ax, which=ns_idx,
                bins=80,
                range=diffXlim,
                xlabel='Gridness score difference')
        ax.set_xlim(diffXlim)
        l = ax.legend([diffLegend[ns_idx]], loc='best', fontsize='small',
                frameon=False)

        fname = outputDir + "/grids_diff_histograms_{0}pA.pdf"
        plt.savefig(fname.format(int(ps.noise_sigmas[ns_idx])), dpi=300, transparent=True)
        plt.close()

    # Plot them all in one plot
    fig = plt.figure(figsize=diffFigSize)
    ax = fig.add_axes(Bbox.from_extents(diffLeft, diffBottom, diffRight,
        diffTop))

    segmentation.plotDiffHistograms(data, ps.noise_sigmas,
            ax=ax, which=None,
            bins=80,
            range=diffXlim,
            xlabel='Gridness score difference')
    ax.set_xlim(diffXlim)
    l = ax.legend(diffLegend, loc='best', fontsize='small', frameon=False)

    fname = outputDir + "/grids_diff_histograms_all_in_one.pdf"
    plt.savefig(fname, dpi=300, transparent=True)
    plt.close()


##############################################################################
# Scatter plot of differences
diffScatterFigSize = (10, 10)
diffScatterLeft   = 0.1
diffScatterBottom = 0.1
diffScatterRight  = 0.95
diffScatterTop    = 0.95
diffScatterXlim   = [-1.5, 1.5]
diffScatterXLabel = '$G_{150}$ -  $G_{0}$'
diffScatterYLabel = '$G_{300}$ -  $G_{150}$'


if (diff_scatter):

    data = []
    for ns_idx, _ in enumerate(ps.noise_sigmas):
        d, _, _ = aggr.aggregateType(ps.grids[ns_idx], iterList, gridTypes,
                NTrials, ignoreNaNs=True)
        data.append(d)

    fig = plt.figure(figsize=diffScatterFigSize)
    ax = fig.add_axes(Bbox.from_extents(diffScatterLeft, diffScatterBottom,
        diffScatterRight, diffScatterTop))

    segmentation.plotDiffScatter(data, ps.noise_sigmas,
            ax=ax,
            s=25,
            linewidth=0.3,
            xlabel=diffScatterXLabel,
            ylabel=diffScatterYLabel)
    ax.set_xlim(prepareLims(diffScatterXlim))
    ax.set_ylim(prepareLims(diffScatterXlim))
    ax.xaxis.set_major_locator(ti.MultipleLocator(0.5))
    ax.yaxis.set_major_locator(ti.MultipleLocator(0.5))

    fname = outputDir + "/grids_diff_scatter.pdf"
    plt.savefig(fname, dpi=300, transparent=True)
    plt.close()



##############################################################################
# Scatter plot with histograms all in one plot
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
        [-np.infty, -0.2, 0.2, np.infty],
        [-np.infty, -0.1, 0.1, np.infty]]
#segMergeInfo = [
#        [0, 1, 3],
#        [6, 7, 8]]
segMergeInfo = None

if diff_all:
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

if sweep_segmentation:
    fig = plt.figure(figsize=sweepFigSize)
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))

    segmentation.plotSweepSegments(ps, iterList, segThresholds, gridTypes,
            NTrials,
            cmap='Set1',
            ax=ax,
            mergeInfo=segMergeInfo)

    fname = outputDir + "/grids_diff_sweep_segments.pdf"
    plt.savefig(fname, dpi=300, transparent=True)
    plt.close()
