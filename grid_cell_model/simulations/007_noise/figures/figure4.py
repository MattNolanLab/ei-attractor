#!/usr/bin/env python
#
#   figure4.py
#
#   Final figure: network mechanisms
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

import EI_plotting as EI
from figures_shared import NoiseDataSpaces
from plotting.global_defs import globalAxesSettings
from analysis.visitors import SpikeTrainXCVisitor
from parameters import JobTrialSpace2D

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

outputDir = "."
NTrials = 5
iterList  = ['g_AMPA_total', 'g_GABA_total']

noise_sigmas  = [0, 150, 300]
rasterRC      = [(5, 15), (5, 15), (5, 15)] # (row, col)
bumpDataRoot  = 'output_local/even_spacing/gamma_bump'
velDataRoot   = 'output_local/even_spacing/velocity'
gridsDataRoot = 'output_local/even_spacing/grids'
shape = (31, 31)

ccExamples      = 0
spikeCCExamples = 0
velLines        = 1
detailed_noise  = 1

###############################################################################

def plotCCExamples(dataSpaces, ns_idx, r, c, trialNum=0, **kw):
    nPairs = kw.pop('npairs', 5)
    ax     = kw.pop('ax', plt.gca())

    kw['xlabel'] = kw.get('xlabel', 'Lag (ms)')
    kw['ylabel'] = kw.get('ylabel', 'C')
    kw['nticks'] = 2
    ylabel = kw.get('ylabel', 'C')

    space = dataSpaces.bumpGamma[ns_idx]
    noise_sigma = dataSpaces.noise_sigmas[ns_idx]
    data = space[r][c][trialNum].data
    xcIn = data['analysis']['x-corr'] 
    CCs =  xcIn['correlations']
    lags = xcIn['lags']
    for idx in xrange(nPairs):
        cc = CCs[idx][2*idx+1]
        signalPlot(lags, cc, ax, zeroLine=False, xmargin=0.02, nThetaTicks=2, **kw)
    ax.set_yticks([0, 1])
    ax.set_ylim([-0.02, 1.02])
    ax.grid(False)
    ax.xaxis.labelpad = 0.8

    
def plotSpikeXCExamples(dataSpaces, ns_idx, r, c, neuronIdx, nOthers, popType,
        trialNum=0, **kw):
    ax = kw.pop('ax', plt.gca())
    xlabel = kw.pop('xlabel', 'Lag (ms)')

    if (popType == 'E'):
        XCName = 'XCorrelation_e'
    elif (popType == 'I'):
        XCName = 'XCorrelation_i'
    else:
        raise ValueError()

    space = dataSpaces.grids[ns_idx]
    noise_sigma = dataSpaces.noise_sigmas[ns_idx]
    data = space[r][c][trialNum].data['analysis'][XCName]
    bin_centers = data['bin_centers']
    CCs = data['correlations']

    for other in xrange(nOthers):
        C = CCs[neuronIdx][1 + other]
        ax.plot(bin_centers, C, '-', **kw)

    globalAxesSettings(ax)
    ax.set_xticks([bin_centers[0], 0, bin_centers[-1]])
    ax.yaxis.set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    

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

    ax.plot(IvelVec, slopes.T, 'o', markerfacecolor='none', markersize=markersize, **kw)
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






###############################################################################
roots = NoiseDataSpaces.Roots(bumpDataRoot, velDataRoot, gridsDataRoot)
ps    = NoiseDataSpaces(roots, shape, noise_sigmas)


##############################################################################
velFigsize =(2.6, 2)
velLeft    = 0.3
velBottom  = 0.35
velRight   = 0.95
velTop     = 0.65

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
    fname = outputDir + "/figure4_slope_examples_0.pdf"
    plt.savefig(fname, dpi=300, transparent=True)

    fig = plt.figure(figsize=(2.5, velFigsize[1]))
    ax = fig.add_axes(Bbox.from_extents(0.3, velBottom, velRight,
        velTop))
    plotSlopes(ax, ps.v[1], positions[1], noise_sigma=ps.noise_sigmas[1],
            xlabel='', xticks=False,
            g_ann=False,
            color='green')
    fname = outputDir + "/figure4_slope_examples_1.pdf"
    plt.savefig(fname, dpi=300, transparent=True)

    fig = plt.figure(figsize=(2.5, velFigsize[1]))
    ax = fig.add_axes(Bbox.from_extents(0.3, velBottom, velRight,
        velTop))
    plotSlopes(ax, ps.v[2], positions[2], noise_sigma=ps.noise_sigmas[2],
            ylabel='',
            g_ann=False,
            color='red')
    fname = outputDir + "/figure4_slope_examples_2.pdf"
    plt.savefig(fname, dpi=300, transparent=True)


##############################################################################
EI13Root  = 'output_local/detailed_noise/velocity/EI-1_3'
EI31Root  = 'output_local/detailed_noise/velocity/EI-3_1'
detailedShape = (31, 9)

EI13PS = JobTrialSpace2D(detailedShape, EI13Root)
EI31PS = JobTrialSpace2D(detailedShape, EI31Root)
detailedNTrials = None


sliceFigSize = (4.3, 2.5)
sliceLeft   = 0.2
sliceBottom = 0.3
sliceRight  = 0.99
sliceTop    = 0.8
if (detailed_noise):
    ylabelPos = -0.2

    types = ('velocity', 'slope')
    fig = plt.figure(figsize=sliceFigSize)
    ax = fig.add_axes(Bbox.from_extents(sliceLeft, sliceBottom, sliceRight,
        sliceTop))
    EI.plotDetailedNoise(EI13PS, detailedNTrials, types, ax=ax,
            ylabelPos=ylabelPos,
            xlabel='', xticks=False,
            color='black')
    EI.plotDetailedNoise(EI31PS, detailedNTrials, types, ax=ax,
            xlabel='', xticks=False,
            ylabel='Slope\n(neurons/s/pA)', ylabelPos=ylabelPos,
            color='red')
    ax.yaxis.set_major_locator(ti.MultipleLocator(0.4))
    ax.yaxis.set_minor_locator(ti.MultipleLocator(0.2))

    fname = "figure4_detailed_noise_slope.pdf"
    plt.savefig(fname, dpi=300, transparent=True)
    plt.close()


    types = ('velocity', 'fitErr')
    fig = plt.figure(figsize=sliceFigSize)
    ax = fig.add_axes(Bbox.from_extents(sliceLeft, sliceBottom, sliceRight,
        sliceTop))
    _, p13, l13 = EI.plotDetailedNoise(EI13PS, detailedNTrials, types, ax=ax,
            ylabelPos=ylabelPos, color='black')
    _, p31, l31 = EI.plotDetailedNoise(EI31PS, detailedNTrials, types, ax=ax,
            ylabel='Fit error\n(neurons/s/trial)', ylabelPos=ylabelPos,
            color='red')
    ax.yaxis.set_major_locator(ti.MultipleLocator(4))
    ax.yaxis.set_minor_locator(ti.MultipleLocator(2))
    leg = ['($g_E,\ g_I$) = (1, 3) nS',  '($g_E,\ g_I$) = (3, 1) nS',
            'Mean', 'Mean']
    ax.legend([p13, p31, l13, l31], leg, loc=(0.55, 0.3), fontsize='small', frameon=False,
            numpoints=1)

    fname = "figure4_detailed_noise_error.pdf"
    plt.savefig(fname, dpi=300, transparent=True)
    plt.close()





##############################################################################
##############################################################################
boxFigSize = (1.75, 1.5)
boxLeft    = 0.25
boxBottom  = 0.32
boxRight   = 0.9
boxTop     = 0.93

if (ccExamples):
    for idx, noise_sigma in enumerate(ps.noise_sigmas):
        fig = plt.figure(figsize=boxFigSize)
        ax = fig.add_axes(Bbox.from_extents(boxLeft, boxBottom, boxRight,
            boxTop))
        kw = {}
        if (idx != 0):
            kw['ylabel'] = ''

        plotCCExamples(ps, ns_idx=idx,
                r=rasterRC[idx][0], c=rasterRC[idx][1],
                ax=ax, **kw)
        fname = outputDir + "/figure4_cc_examples{0}.pdf".format(noise_sigma)
        fig.savefig(fname, dpi=300, transparent=transparent)
        plt.close()


if (spikeCCExamples):
    lagRange    = (-125, 125)
    bins        = 251
    neuronIdx   = np.arange(4)
    plotIdx     = 0
    nOthers     = 3
    forceUpdate = False
    for idx, noise_sigma in enumerate(ps.noise_sigmas):
        r, c = rasterRC[idx][0], rasterRC[idx][1]
        trialNum = 0

        # E cells
        fig = plt.figure(figsize=boxFigSize)
        ax = fig.add_axes(Bbox.from_extents(boxLeft, boxBottom, boxRight,
            boxTop))

        visitor_e = SpikeTrainXCVisitor("spikeMon_e", bins,
                lagRange=lagRange, 
                neuronIdx=neuronIdx,
                forceUpdate=forceUpdate)
        dataSet_e = ps.grids[idx][r][c][trialNum]
        visitor_e.visitDictDataSet(dataSet_e)

        plotSpikeXCExamples(ps, idx, r, c, plotIdx, nOthers, 'E')

        fname = outputDir + "/figure4_spike_xc_E_examples{0}.pdf".format(noise_sigma)
        fig.savefig(fname, dpi=300, transparent=transparent)
        plt.close()

        # Icells
        fig = plt.figure(figsize=boxFigSize)
        ax = fig.add_axes(Bbox.from_extents(boxLeft, boxBottom, boxRight,
            boxTop))
        visitor_i = SpikeTrainXCVisitor("spikeMon_i", bins,
                lagRange=lagRange, 
                neuronIdx=neuronIdx,
                forceUpdate=forceUpdate)
        dataSet_i = ps.grids[idx][r][c][trialNum]
        visitor_i.visitDictDataSet(dataSet_i)

        plotSpikeXCExamples(ps, idx, r, c, plotIdx, nOthers, 'I')

        fname = outputDir + "/figure4_spike_xc_I_examples{0}.pdf".format(noise_sigma)
        fig.savefig(fname, dpi=300, transparent=transparent)
        plt.close()

