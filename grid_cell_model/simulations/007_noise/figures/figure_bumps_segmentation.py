#!/usr/bin/env python
'''
Bump width segmentation plots.
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ti
from matplotlib.transforms import Bbox

import default_settings as ds
from EI_plotting      import sweeps, segmentation, scatter
from EI_plotting      import aggregate as aggr
from EI_plotting.base import NoiseDataSpaces
from grid_cell_model.plotting.low_level   import zeroLines
from grid_cell_model.plotting.global_defs import prepareLims
from grid_cell_model.analysis import clustering
from submitting import flagparse

outputDir = ds.figOutputDir

gridTypes = ['grids', 'gridnessScore']
bumpTypes = ['bump', 'sigma']
gridNTrials = 3
bumpNTrials = 5

exampleIdx     = [(1, 22), (1, 22), (1, 22)] # (row, col)

parser = flagparse.FlagParser()
parser.add_flag('-a', '--diff_all')
parser.add_flag('-s', '--diff_sweep')
parser.add_flag('--scatter_diff_bump_grids')
parser.add_flag('--scatter_diff_bump_grids_seg')
args = parser.parse_args()

ps = ds.getDefaultParamSpaces()


##############################################################################
# Scatter plot with histograms all in one plot
sigmaBumpText = '$\sigma_{bump}^{-1}\ (neurons^{-1})$'

diffAllFigSize = (6, 6)  # Must be square
diffAllLeft = 0.2
diffAllBottom = 0.15
diffAllScatterSize = 0.55
diffAllScatterRight = diffAllLeft + diffAllScatterSize
diffAllScatterTop   = diffAllBottom + diffAllScatterSize
diffHistSize = 0.2
diffHistOffset = 0.0
diffAllXlim   = prepareLims([-0.5, 0.5])
diffAllXLabel = '$\Delta^{{150 - 0\ pA}}$({0})'.format(sigmaBumpText)
diffAllYLabel = '$\Delta^{{300 - 150\ pA}}$({0})'.format(sigmaBumpText)
diffAllBins   = 40

#segThresholds = [
#        [-np.infty, 0, np.infty],
#        [-np.infty, 0, np.infty]]
segThresholds = None

#segMergeInfo = [
#        [0, 1, 3],
#        [6, 7]]
segMergeInfo = None


if args.diff_all or args.all:
    fig = plt.figure(figsize=diffAllFigSize)

    data = []
    for ns_idx, _ in enumerate(ps.noise_sigmas):
        d, _, _ = aggr.aggregateType(ps.bumpGamma[ns_idx], ds.iterList, bumpTypes,
                bumpNTrials, ignoreNaNs=False)
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
            doSegmentation=False, thresholds=segThresholds, mergeInfo=segMergeInfo)
    ax_scatter.set_xlim(diffAllXlim)
    ax_scatter.set_ylim(diffAllXlim)
    ax_scatter.xaxis.set_major_locator(ti.MultipleLocator(0.5))
    ax_scatter.yaxis.set_major_locator(ti.MultipleLocator(0.5))
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
            range=diffAllXlim,
            xlabel='', ylabel='')
    ax_XHist.set_xlim(diffAllXlim)
    ax_XHist.xaxis.set_major_locator(ti.MultipleLocator(0.5))
    ax_XHist.axis('off')

    # 300 - 150 pA histogram
    diffYHistLeft = diffAllScatterRight + diffHistOffset
    diffYHistRight = diffYHistLeft + diffHistSize
    ax_YHist = fig.add_axes(Bbox.from_extents(diffYHistLeft, diffAllBottom,
        diffYHistRight, diffAllScatterTop))
    segmentation.plotDiffHistograms(data, ps.noise_sigmas,
            ax=ax_YHist, which=1,
            bins=diffAllBins,
            range=diffAllXlim,
            orientation='horizontal',
            xlabel='', ylabel='')
    ax_YHist.set_ylim(diffAllXlim)
    ax_YHist.yaxis.set_major_locator(ti.MultipleLocator(0.5))
    ax_YHist.axis('off')


    # Save
    #gs.tight_layout(fig, pad=0.01)
    fname = outputDir + "/bumps_diff_all_plots.pdf"
    plt.savefig(fname, dpi=300, transparent=True)
    #plt.close()



##############################################################################
# Parameter sweep of the difference between noise_150 and noise_0
exampleRC = ( (5, 15), (15, 5) )

ann0 = dict(
        txt='b',
        rc=exampleRC[0],
        xytext_offset=(1.5, 1),
        color='black')
ann1 = dict(
        txt='a',
        rc=exampleRC[1],
        xytext_offset=(0.5, 1.5),
        color='black')
ann = [ann0, ann1]


sigmaBumpText = '$\Delta_{150 - 0}[P(bumps)]$'
bumpDiff_cbar_kw = dict(
        label       = sigmaBumpText,
        location    = 'right',
        shrink      = 0.8,
        pad         = -0.05,
        ticks       = ti.MultipleLocator(0.3),
        rasterized  = True)

if args.diff_sweep or args.all:
    dataList = []
    for ns_idx, _ in enumerate(ps.noise_sigmas):
        data = aggr.IsBump(ps.bumpGamma[ns_idx], ds.iterList, ignoreNaNs=True)
        dataList.append(data)

    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas[0:-1]):
        fig, ax = ds.getDefaultSweepFig()

        which = ns_idx
        _, ax, cax = sweeps.plotDiffTrial(dataList, None, which, None, None,
                ax=ax,
                cbar=True, cbar_kw=bumpDiff_cbar_kw,
                symmetricLimits=True,
                annotations=ann,
                cmap='RdBu_r')
        #cax.yaxis.set_major_locator(ti.MultipleLocator(.9))
        #cax.yaxis.set_minor_locator(ti.MultipleLocator(.45))

        fname = outputDir + "/bumps_isBumpFracTotal_diff_sweeps{0}.pdf"
        plt.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)
        plt.close()


##############################################################################
# Correlate (difference between sigma_{bump}) and gridness score.
corrDiffFigsize = (4, 5)
corrDiffXLabel = sigmaBumpText
corrDiffYLabel = '$\Delta$ Gridness score'

gridTypes = ['grids', 'gridnessScore']
bumpTypes = ['bump', 'sigma']
gridNTrials = 3
bumpNTrials = 5
includeGridSegmentation = True
segGridThresholds = [ # TODO: merge this with figure_grids_segmentation.py
        [-np.infty, 0, np.infty],
        [-np.infty, 0, np.infty]]

def setCorrAxes(ax):
    ax.set_xlim(prepareLims([-0.2, 0.4]))
    ax.set_ylim(prepareLims([-1.5, 1.5]))
    ax.xaxis.set_major_locator(ti.MultipleLocator(0.2))
    ax.yaxis.set_major_locator(ti.MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(ti.MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(ti.MultipleLocator(0.25))
    zeroLines(ax)

# Highlight segmented points
if args.scatter_diff_bump_grids_seg or args.all:
    fig = plt.Figure(corrDiffFigsize)
    ax = fig.add_subplot(111)

    which = 0
    scatterPlot = scatter.DiffScatterPlot(
            ps.bumpGamma, ps.grids, bumpTypes, gridTypes, ds.iterList,
            bumpNTrials, gridNTrials, which,
            s=15,
            linewidth=0.3,
            edgecolor='white',
            xlabel = corrDiffXLabel,
            ylabel = corrDiffYLabel,
            sigmaTitle=False,
            ignoreNaNs=True,
            cmap='Set1',
            ax=ax)

    if includeGridSegmentation:
        gridData = scatterPlot.diffData2
        differences = [
                gridData[0, :],
                gridData[1, :]]
        clusters = clustering.ThresholdClusters(differences, segGridThresholds)
        scatterPlot.setColors(clusters.assignClusters())

    scatterPlot.plot()
    ax.set_title('Difference\n $\sigma_{noise} = 150\ -\ \sigma_{noise} = 0$ pA')
    setCorrAxes(ax)

    fig.tight_layout()
    fname = outputDir + "/bumps_scatter_diff_bump_grids_segments.pdf"
    fig.savefig(fname, dpi=300, transparent=True)
    plt.close()


# Color coding according to E/I position
scatterColorFigSize = (1, 1)
scatterTransparent = True

if args.scatter_diff_bump_grids or args.all:
    fig = plt.Figure(corrDiffFigsize)
    ax = fig.add_subplot(111)

    which = 0
    scatterPlot = scatter.DiffScatterPlot(
            ps.bumpGamma, ps.grids, bumpTypes, gridTypes, ds.iterList,
            bumpNTrials, gridNTrials, which,
            s=15,
            linewidth=0.3,
            edgecolor='white',
            color2D=True,
            xlabel = corrDiffXLabel,
            ylabel = corrDiffYLabel,
            sigmaTitle=False,
            ignoreNaNs=True,
            ax=ax)

    scatterPlot.plot()
    setCorrAxes(ax)

    fig.tight_layout()
    fname = outputDir + "/bumps_scatter_diff_bump_grids.pdf"
    fig.savefig(fname, dpi=300, transparent=True)
    plt.close()

    fig = plt.figure(figsize=scatterColorFigSize)
    ax = fig.gca()
    scatterPlot.plotColorbar(ax)
    fig.tight_layout(pad=0)
    fname = outputDir + "/bumps_scatter_diff_bump_grids_colorbar.pdf"
    fig.savefig(fname, dpi=300, transparent=scatterTransparent)

