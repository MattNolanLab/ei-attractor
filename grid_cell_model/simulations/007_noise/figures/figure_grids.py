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
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.ticker as ti
from matplotlib.transforms import Bbox
import matplotlib.gridspec as gridspec

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

grids              = 0
examplesFlag       = 0
detailed_noise     = 0
Vm_examples        = 0
collapsed_sweeps   = 0
diff_distrib       = 0
diff_scatter       = 0
diff_all           = 1
sweep_segmentation = 1

##############################################################################
roots = NoiseDataSpaces.Roots(bumpDataRoot, velDataRoot, gridsDataRoot)
ps    = NoiseDataSpaces(roots, shape, noise_sigmas)


sweepFigSize = (3.7, 2.6)
sweepLeft   = 0.08
sweepBottom = 0.2
sweepRight  = 0.8
sweepTop    = 0.85
transparent  = True

cbar_kw= {
    'label'      : 'Gridness score',
    'location'   : 'right',
    'shrink'     : 0.8,
    'pad'        : -0.05,
    'ticks'      : ti.MultipleLocator(0.5),
    'rasterized' : True}

vmin = -0.5
vmax = 1.1

##############################################################################
exampleRC = ( (5, 15), (15, 5) )
sliceAnn = None

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


varList = ['gridnessScore']

if (grids):
    # noise_sigma = 0 pA
    fig = plt.figure("sweeps0", figsize=sweepFigSize)
    exRows = [28, 15]
    exCols = [3, 15]
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    sweeps.plotGridTrial(ps.grids[0], varList, iterList, ps.noise_sigmas[0],
            trialNumList=gridTrialNumList,
            ax=ax,
            r=exampleIdx[0][0], c=exampleIdx[0][1],
            cbar=False, cbar_kw=cbar_kw,
            vmin=vmin, vmax=vmax,
            ignoreNaNs=True,
            annotations=ann,
            sliceAnn=sliceAnn)
    fname = outputDir + "/grids_sweeps0.pdf"
    fig.savefig(fname, dpi=300, transparent=True)
    plt.close()


    # noise_sigma = 150 pA
    #for a in sliceAnn:
    #    a['letterColor'] = 'black'
    fig = plt.figure("sweeps150", figsize=sweepFigSize)
    exRows = [8, 2]
    exCols = [10, 9]
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    sweeps.plotGridTrial(ps.grids[1], varList, iterList, ps.noise_sigmas[1],
            trialNumList=gridTrialNumList,
            ax=ax,
            r=exampleIdx[1][0], c=exampleIdx[1][1],
            cbar=False, cbar_kw=cbar_kw,
            vmin=vmin, vmax=vmax,
            ignoreNaNs=True,
            ylabel='', yticks=False,
            annotations=ann,
            sliceAnn=sliceAnn)
    fname = outputDir + "/grids_sweeps150.pdf"
    fig.savefig(fname, dpi=300, transparent=True)
    plt.close()


    # noise_sigma = 300 pA
    fig = plt.figure("sweeps300", figsize=sweepFigSize)
    exRows = [16, 15]
    exCols = [6, 23]
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    _, _, cax = sweeps.plotGridTrial(ps.grids[2], varList, iterList, ps.noise_sigmas[2],
            trialNumList=gridTrialNumList,
            ax=ax,
            r=exampleIdx[2][0], c=exampleIdx[2][1],
            cbar_kw=cbar_kw,
            vmin=vmin, vmax=vmax,
            ignoreNaNs=True,
            ylabel='', yticks=False,
            annotations=ann,
            sliceAnn=sliceAnn)
    #for label in cax.yaxis.get_ticklabels():
    #    label.set_ha('right')
    #cax.tick_params(pad=30)
    fname = outputDir + "/grids_sweeps300.pdf"
    fig.savefig(fname, dpi=300, transparent=True)
    plt.close()


##############################################################################
# Grid field examples
exampleGridFName = outputDir + "/grids_examples_{0}pA_{1}.pdf"
exampleACFName = outputDir + "/grids_examples_acorr_{0}pA_{1}.pdf"
exTransparent = True
exampleFigSize = (1, 1.2)
exampleLeft   = 0.01
exampleBottom = 0.01
exampleRight  = 0.99
exampleTop    = 0.85

if (examplesFlag):
    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        for idx, rc in enumerate(exampleRC):
            # Grid field
            fname = exampleGridFName.format(noise_sigma, idx)
            fig = plt.figure(figsize=exampleFigSize)
            gs = examples.plotOneGridExample(ps.grids[ns_idx], rc, iterList,
                    exIdx=exampleIdx[idx],
                    xlabel=False, ylabel=False,
                    xlabel2=False, ylabel2=False, 
                    maxRate=True, plotGScore=False,
                    fig=fig)
            gs.update(left=exampleLeft, bottom=exampleBottom, right=exampleRight,
                    top=exampleTop)
            plt.savefig(fname, dpi=300, transparent=exTransparent)
            plt.close()

            # Autocorrelation
            fname = exampleACFName.format(noise_sigma, idx)
            fig= plt.figure(figsize=exampleFigSize)
            ax = fig.add_axes(Bbox.from_extents(exampleLeft, exampleBottom,
                exampleRight, exampleTop))
            gs = examples.plotOneGridACorrExample(ps.grids[ns_idx], rc, ax=ax)
            plt.savefig(fname, dpi=300, transparent=exTransparent)
            plt.close()
    

##############################################################################
# Membrane potential examples

def openJob(rootDir, noise_sigma):
    fileTemplate = "noise_sigma{0}_output.h5"
    fileName = rootDir + '/' + fileTemplate.format(int(noise_sigma))
    return DataStorage.open(fileName, 'r')

def drawVm(data, noise_sigma, xScaleBar=None, yScaleBar=None,
        ax=plt.gca(), sigmaTitle=True):
    yScaleX = 0.5
    yScaleY = 1.1
    yScaleXOffset = 0.06
    scaleTextSize = 'x-small'

    plotTStart = 5e3
    plotTEnd   = 5.25e3

    stateYlim = [-80, -40]

    theta_start_t = getOption(data, 'theta_start_t')
    #theta_start_t = 1e3
    simTime = getOption(data, 'time')

    mon_e = data['stateMon_e']

    # E cell Vm
    t, VmMiddle = extractSummedSignals(mon_e, ['V_m'], plotTStart,
            plotTEnd)
    plotStateSignal(ax, t, VmMiddle, labely='', color='red',
            scaleBar=xScaleBar, scaleText="ms", scaleX=0.5, scaleY=1.1,
            scaleTextSize=scaleTextSize)
    if (yScaleBar is not None):
        plotting.low_level.yScaleBar(yScaleBar, yScaleX, yScaleY,
                ax=ax,
                unitsText='mV', textXOffset=yScaleXOffset,
                size='x-small')
    ax.yaxis.set_visible(False)
    ax.spines['left'].set_visible(False)
    plt.ylim(stateYlim)

    if (sigmaTitle):
        ax.set_title('$\sigma$ = {0} pA'.format(int(noise_sigma)), loc='left',
                size=scaleTextSize, y=0.9)


VmExampleFigSize = (2.5, 1.25)
VmExampleLeft   = 0.01
VmExampleBottom = 0.01
VmExampleRight  = 0.999
VmExampleTop    = 0.6
VmExampleXScalebar = 50 # ms
VmExampleYScalebar = 10 # mV

if (Vm_examples):
    for ns_idx, noise_sigma in enumerate(noise_sigmas):
        fig = plt.figure(figsize=VmExampleFigSize)
        ax = fig.add_axes(Bbox.from_extents(VmExampleLeft, VmExampleBottom,
            VmExampleRight, VmExampleTop))
        ds = openJob(singleDataRoot, noise_sigma)
        kw = {}
        if (ns_idx == 2):
            kw['xScaleBar'] = VmExampleXScalebar
            kw['yScaleBar'] = VmExampleYScalebar
        drawVm(ds, noise_sigma=noise_sigma, ax=ax, **kw)

        fname = outputDir + "/grids_Vm_example_{0}.pdf"
        plt.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)
        plt.close()


##############################################################################
# Detailed noise plots
EI13Root  = 'output_local/detailed_noise/grids/EI-1_3'
EI31Root  = 'output_local/detailed_noise/grids/EI-3_1'
detailedShape = (31, 9)

EI13PS = JobTrialSpace2D(detailedShape, EI13Root)
EI31PS = JobTrialSpace2D(detailedShape, EI31Root)
detailedNTrials = 1


detailFigSize = (4.5, 2.8)
detailLeft   = 0.15
detailBottom = 0.2
detailRight  = 0.95
detailTop    = 0.8
if (detailed_noise):
    ylabelPos = -0.15

    types = ('grids', 'gridnessScore')
    fig = plt.figure(figsize=detailFigSize)
    ax = fig.add_axes(Bbox.from_extents(detailLeft, detailBottom, detailRight,
        detailTop))
    _, p13, l13 = details.plotDetailedNoise(EI13PS, detailedNTrials, types, ax=ax,
            ylabelPos=ylabelPos,
            color='red', markerfacecolor='red', zorder=10)
    _, p31, l31 = details.plotDetailedNoise(EI31PS, detailedNTrials, types, ax=ax,
            ylabel='Gridness score', ylabelPos=ylabelPos,
            color='#505050')
    ax.yaxis.set_major_locator(ti.MultipleLocator(0.4))
    ax.yaxis.set_minor_locator(ti.MultipleLocator(0.2))
    ax.set_ylim([-0.6, 1.2])
    leg = ['a',  'b']
    l = ax.legend([p31, p13], leg, loc=(0.8, 1), fontsize='small', frameon=False,
            numpoints=1, handletextpad=0.05)

    fname = outputDir + "/grids_detailed_noise_gscore.pdf"
    plt.savefig(fname, dpi=300, transparent=True)
    plt.close()



##############################################################################
# Collapsed sweep plots
collapsedFigSize = (3.5, detailFigSize[1])
collapsedLeft   = 0.1
collapsedBottom = detailBottom
collapsedRight  = 0.95
collapsedTop    = detailTop

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
