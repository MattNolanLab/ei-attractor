#!/usr/bin/env python
#
#   suppFigure_grid_slices.py
#
#   Supplementary figure for the gridness score parameter sweeps.
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

from EI_plotting          import details
from EI_plotting          import aggregate as aggr
from plotting.global_defs import globalAxesSettings
from EI_plotting.base     import plotOneHist, NoiseDataSpaces
from parameters           import JobTrialSpace2D
import flagparse

from matplotlib import rc
rc('pdf', fonttype=42)
rc('mathtext', default='regular')

plt.rcParams['font.size'] = 11

outputDir = "panels"

NTrials=3
gridTrialNumList = np.arange(NTrials)
iterList  = ['g_AMPA_total', 'g_GABA_total']

noise_sigmas  = [0, 150, 300]
exampleIdx    = [(1, 22), (1, 22), (1, 22)] # (row, col)
bumpDataRoot  = None
velDataRoot   = None
gridsDataRoot = 'output_local/even_spacing/grids'
shape = (31, 31)

parser = flagparse.FlagParser()
parser.add_flag('--hists')
parser.add_flag('--slices')
args = parser.parse_args()

##############################################################################



def plotGridnessThresholdComparison(spList, trialNumList, thrList, r=0, c=0,
        ylabelPos=-0.2):
    varList = ['gridnessScore']
    #G = []
    noise_sigma = [0, 150, 300]

    ax = plt.gca()
    plt.hold('on')
    globalAxesSettings(ax)

    for sp in spList:
        G = aggr.aggregate2DTrial(sp, varList, trialNumList).flatten()
        counts = []
        for thr in thrList:
            thrIdx = np.logical_and(G > thr, np.logical_not(np.isnan(G)))
            counts.append(float(len(G[thrIdx])) / len(G))
        plt.plot(thrList, counts, 'o-', markersize=4)


    plt.plot([0], [1], linestyle='None', marker='None')
    ax.set_xlabel('Gridness score threshold')
    ax.text(ylabelPos, 0.5, 'Count', rotation=90, transform=ax.transAxes,
            va='center', ha='right')
    leg = []
    for s in noise_sigma:
        leg.append("{0}".format(int(s)))
    l = ax.legend(leg, loc=(0.8, 0.5), title='$\sigma$ (pA)', frameon=False,
            fontsize='x-small', ncol=1)
    plt.setp(l.get_title(), fontsize='x-small')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_major_locator(ti.MultipleLocator(0.3))
    ax.yaxis.set_major_locator(ti.MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(ti.AutoMinorLocator(3))
    ax.yaxis.set_minor_locator(ti.AutoMinorLocator(2))
    ax.margins(0.02, 0.07)


def plotGridnessHistogram(spList, trialNumList, ylabelPos=-0.2):
    varList = ['gridnessScore']
    #G = []
    noise_sigma = [0, 150, 300]
    colors = ['red', 'green', 'blue']

    ax = plt.gca()
    plt.hold('on')
    globalAxesSettings(ax)

    for idx, sp in enumerate(spList):
        G = aggr.aggregate2DTrial(sp, varList, trialNumList).flatten()
        filtIdx = np.logical_not(np.isnan(G))
        plotOneHist(G[filtIdx], normed=True)
    leg = []
    for s in noise_sigma:
        leg.append("{0}".format(int(s)))
    l = ax.legend(leg, loc=(0.8, 0.5), title='$\sigma$ (pA)', frameon=False,
            fontsize='x-small', ncol=1)
    plt.setp(l.get_title(), fontsize='x-small')

    ax.set_xlabel("Gridness score")
    ax.text(ylabelPos, 0.5, 'p(G)', rotation=90, transform=ax.transAxes,
            va='center', ha='right')
    #ax.set_ylabel("p(G)")

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_major_locator(ti.MultipleLocator(0.5))
    ax.yaxis.set_major_locator(ti.MultipleLocator(2))
    ax.xaxis.set_minor_locator(ti.AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(ti.AutoMinorLocator(2))
    ax.margins(0.05, 0.025)
    

def computeMarginal(G, type, X, Y, ignoreNaNs):
    trials = []
    for trialIdx in xrange(G.shape[2]):
        trials.append(G[:, :, trialIdx])

    res    = None
    x      = None
    xlabel = ''
    if (type == 'horizontal'):
        res = np.vstack(trials)
        if (ignoreNaNs):
            nans = np.isnan(res)
            res = ma.MaskedArray(res, mask=nans)
        res = res.T
        x   = X[0, :]
        xlabel = EI.xlabelText
    elif (type == 'vertical'):
        res = np.hstack(trials)
        if (ignoreNaNs):
            nans = np.isnan(res)
            res = ma.MaskedArray(res, mask=nans)
        x   = Y[:, 0]
        xlabel = EI.ylabelText
    else:
        raise ValueError("Marginal type must be 'horizontal' or 'vertical'")
    #import pdb; pdb.set_trace()
    return res, x, xlabel


def plotGridnessMarginal(paramSpaces, type, NTrials=1, **kw):
    # kwargs
    title        = kw.pop('title', True)
    rcG          = kw.pop('rowsCols', [(1, 22), (1, 22), (1, 22)]) # (row, col)
    ax           = kw.pop('ax', plt.gca())
    iterList     = kw.pop('iterList', ['g_AMPA_total', 'g_GABA_total'])
    kw['ylabel'] = kw.get('ylabel', 'Gridness score')
    ignoreNaNs   = kw.get('ignoreNaNs', True)

    GVars = ['gridnessScore']
    trialNumList = range(NTrials)
    sp = paramSpaces

    # Gridness score
    for idx, noise_sigma in enumerate(sp.noise_sigmas):
        space = sp.grids[idx]
        G = space.aggregateData(GVars, trialNumList, output_dtype='array',
                loadData=True, saveData=False)
        Y, X = aggr.computeYX(space, iterList, r=rcG[idx][0], c=rcG[idx][1])
        marginal, x, kw['xlabel'] = computeMarginal(G, type, X, Y, ignoreNaNs)
        details.plotOneSlice(ax, x, marginal, **kw)
    ax.yaxis.set_major_locator(ti.MaxNLocator(4))

    return ax
        


##############################################################################
roots = NoiseDataSpaces.Roots(bumpDataRoot, velDataRoot, gridsDataRoot)
ps    = NoiseDataSpaces(roots, shape, noise_sigmas)

slice_horizontal = slice(13, 18)
slice_vertical = slice(13, 18)

sweepFigSize = (3.5, 2.5)
sweepLeft   = 0.15
sweepBottom = 0.2
sweepRight  = 0.87
sweepTop    = 0.85

cbar_kw= {'label' : 'Gridness score',
    'orientation': 'vertical',
    'shrink': 0.8,
    'pad' : -0.05,
    'ticks' : ti.MultipleLocator(0.5),
    'rasterized' : True}

vmin = -0.5
vmax = 1.1

# Stats
sliceFigSize = (3.7, 2)
sliceLeft   = 0.2
sliceBottom = 0.3
sliceRight  = 0.99
sliceTop    = 0.85
if args.hists or args.all:
    ylabelPos = -0.16
    fig = plt.figure(figsize=sliceFigSize)
    ax = fig.add_axes(Bbox.from_extents(sliceLeft, sliceBottom, sliceRight,
        sliceTop))
    plotGridnessThresholdComparison(ps.grids, range(NTrials),
            thrList=np.arange(-0.4, 1.2, 0.05), ylabelPos=ylabelPos)
    fname = outputDir + "/suppFigure_grid_slices_threshold_comparison.pdf"
    plt.savefig(fname, dpi=300, transparent=True)


if args.slices or args.all:
    ylabelPos = -0.16
    fig = plt.figure(figsize=sliceFigSize)
    ax = fig.add_axes(Bbox.from_extents(sliceLeft, sliceBottom, sliceRight,
        sliceTop))
    details.plotGridnessSlice(ps, slice_horizontal, slice(None), type='horizontal',
            ax=ax)
    ax.yaxis.set_major_locator(ti.MultipleLocator(0.4))
    ax.yaxis.set_minor_locator(ti.AutoMinorLocator(2))
    ax.set_ylim([-0.5, 1.21])
    fname = outputDir + "/suppFigure_grid_slices_slice_horizontal.pdf"
    plt.savefig(fname, dpi=300, transparent=True)
    plt.close()

    fig = plt.figure(figsize=sliceFigSize)
    ax = fig.add_axes(Bbox.from_extents(sliceLeft, sliceBottom, sliceRight,
        sliceTop))
    details.plotGridnessSlice(ps, slice(None), slice_vertical, type='vertical',
            ax=ax)
    ax.yaxis.set_major_locator(ti.MultipleLocator(0.4))
    ax.yaxis.set_minor_locator(ti.AutoMinorLocator(2))
    ax.set_ylim([-0.5, 1.21])
    fname = outputDir + "/suppFigure_grid_slices_slice_vertical.pdf"
    plt.savefig(fname, dpi=300, transparent=True)
    plt.close()


