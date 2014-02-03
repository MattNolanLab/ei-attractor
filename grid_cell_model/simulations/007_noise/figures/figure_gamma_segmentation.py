#!/usr/bin/env python
#
#   figure_gamma_segmentation.py
#
#   Theta/gamma analysis - segmentation.
#
#       Copyright (C) 2012  Lukas Solanka <l.solanka@sms.ed.ac.uk>
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


from EI_plotting          import segmentation
from EI_plotting          import aggregate as aggr
from plotting.global_defs import prepareLims
from EI_plotting.base     import NoiseDataSpaces
import flagparse

import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)


# Other
from matplotlib import rc
rc('pdf', fonttype=42)
rc('mathtext', default='regular')

plt.rcParams['font.size'] = 11

###############################################################################
cFreq = 'blue'
cAC = 'green'
cCount = 'red'

outputDir = "panels"
NTrials = 5
iterList  = ['g_AMPA_total', 'g_GABA_total']

noise_sigmas  = [0, 150, 300]
exampleIdx    = [(0, 0), (0, 0), (0, 0)] # (row, col)
bumpDataRoot  = 'output_local/even_spacing/gamma_bump'
velDataRoot   = None
gridsDataRoot = None
shape    = (31, 31)

###############################################################################
roots = NoiseDataSpaces.Roots(bumpDataRoot, velDataRoot, gridsDataRoot)
ps    = NoiseDataSpaces(roots, shape, noise_sigmas)

parser = flagparse.FlagParser()
parser.add_flag('--diff_all')
parser.add_flag('--sweep_segmentation')
args = parser.parse_args()

##############################################################################
# Scatter plot with histograms all in one plot
diffXlim = [-0.6, 0.6] # Data analysis

diffAllFigSize = (6, 6)  # Must be square
diffAllLeft = 0.2
diffAllBottom = 0.15
diffAllScatterSize = 0.55
diffAllScatterRight = diffAllLeft + diffAllScatterSize
diffAllScatterTop   = diffAllBottom + diffAllScatterSize
diffHistSize = 0.2
diffHistOffset = 0.0
diffAllXlim   = prepareLims(diffXlim) # plotting boundaries
diffAllXLabel = '$\gamma_{150}$ -  $\gamma_{0}$'
diffAllYLabel = '$\gamma_{300}$ -  $\gamma_{150}$'
histRatio     = 0.2
diffAllBins   = 40

segThresholds = [
        [-np.infty, -0.025, np.infty],
        [-np.infty, 0.15, np.infty]]
#segMergeInfo = [
#        [3, 4]]
segMergeInfo = None
gammaThreshold = -np.infty

gammaTypes = ['gamma', 'acVal']


if args.diff_all or args.all:
    fig = plt.figure(figsize=diffAllFigSize)

    data = []
    for ns_idx, _ in enumerate(ps.noise_sigmas):
        d, _, _ = aggr.aggregateType(ps.bumpGamma[ns_idx], iterList, gammaTypes,
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
            filterThreshold=gammaThreshold,
            doSegmentation=True, thresholds=segThresholds, mergeInfo=segMergeInfo)
    ax_scatter.set_xlim(diffAllXlim)
    ax_scatter.set_ylim(diffAllXlim)
    ax_scatter.xaxis.set_major_locator(ti.MultipleLocator(0.6))
    ax_scatter.yaxis.set_major_locator(ti.MultipleLocator(0.6))
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
            filterThreshold=gammaThreshold,
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
            filterThreshold=gammaThreshold,
            xlabel='', ylabel='')
    ax_YHist.set_ylim(diffAllXlim)
    ax_YHist.yaxis.set_major_locator(ti.MultipleLocator(1.5))
    ax_YHist.axis('off')


    # Save
    #gs.tight_layout(fig, pad=0.01)
    fname = outputDir + "/gamma_diff_all_plots.pdf"
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

    segmentation.plotSweepSegments(ps.bumpGamma, ps.noise_sigmas, iterList,
            segThresholds, gammaTypes, NTrials,
            cmap='Set1',
            ax=ax,
            filterThreshold=gammaThreshold,
            mergeInfo=segMergeInfo)

    fname = outputDir + "/gamma_diff_sweep_segments.pdf"
    plt.savefig(fname, dpi=300, transparent=True)
    plt.close()

