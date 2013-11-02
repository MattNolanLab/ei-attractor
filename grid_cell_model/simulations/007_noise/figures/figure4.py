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
from data_storage.sim_models.ei import MonitoredSpikes
from plotting.global_defs import globalAxesSettings
from plotting.signal import signalPlot
from analysis.visitors import SpikeTrainXCVisitor

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
velDataRoot   = None
gridsDataRoot = 'output_local/even_spacing/grids'
shape = (31, 31)

rasters         = 0
rates           = 0
ccExamples      = 0
spikeCCExamples = 1

###############################################################################

def getSpikes(space, noise_sigma, r, c, trialNum):
    data = space[r][c][trialNum].data
    ESpikes = MonitoredSpikes(data, 'spikeMon_e', 'net_Ne')
    ISpikes = MonitoredSpikes(data, 'spikeMon_i', 'net_Ni')
    return ESpikes, ISpikes


def rasterPlot(dataSpaces, noise_sigma_idx, r, c, tLimits, trialNum=0, **kw):
    title   = kw.pop('title', True)
    ann     = kw.pop('ann', False)
    ann_EI  = kw.pop('ann_EI', False)
    sigmaTitle = kw.pop('sigmaTitle', False)

    space = dataSpaces.bumpGamma[noise_sigma_idx]
    noise_sigma = dataSpaces.noise_sigmas[noise_sigma_idx]
    ESpikes, ISpikes = getSpikes(space, noise_sigma, r, c, trialNum)
    ax = EI.plotEIRaster(ESpikes, ISpikes, tLimits, **kw) 

    if (sigmaTitle):
        ax.set_title('$\sigma$ = {0} pA'.format(int(noise_sigma)), x=0.02, y=1.02,
                va='bottom', ha='left')
    if (ann):
        Y, X = EI.computeYX(space, iterList, r=r, c=c)
        gE = Y[r, c]
        gI = X[r, c]
        txt = '$g_E$ = {0} nS\n$g_I$ = {1} nS'.format(gE, gI)
        ax.text(0.99, 1.02, txt, va='bottom', ha='right', size='small',
                transform=ax.transAxes)

    if (ann_EI):
        ax.text(-0.05, 0.75, 'E', va='center', ha='center', size='small',
                transform=ax.transAxes, color='red', weight='bold')
        ax.text(-0.05, 0.25, 'I', va='center', ha='center', size='small',
                transform=ax.transAxes, color='blue', weight='bold')


    return ax


def plotAvgFiringRate(dataSpaces, noise_sigma_idx, r, c, trialNum=0, **kw):
    # keyword arguments
    ax = kw.pop('ax', plt.gca())
    kw['xlabel'] = False
    kw['ylabel'] = kw.get('ylabel', 'r (Hz)')

    space = dataSpaces.bumpGamma[noise_sigma_idx]
    noise_sigma = dataSpaces.noise_sigmas[noise_sigma_idx]
    ESpikes, ISpikes = getSpikes(space, noise_sigma, r, c, trialNum)

    tStart = tLimits[0]
    tEnd   = tLimits[1]
    dt     = 0.5  # ms
    winLen = 2.0 # ms

    ERate, eTimes = ESpikes.slidingFiringRate(tStart, tEnd, dt, winLen)
    IRate, iTimes = ISpikes.slidingFiringRate(tStart, tEnd, dt, winLen)
    meanERate = np.mean(ERate, axis=0)
    meanIRate = np.mean(IRate, axis=0)

    signalPlot(eTimes, meanERate, ax, color='red', **kw)
    signalPlot(iTimes, meanIRate, ax, color='blue', **kw)

    ax.set_ylim([0, None])
    #ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.xaxis.set_visible(False)
    #ax.yaxis.set_visible(False)
    ax.yaxis.set_major_locator(ti.LinearLocator(2))
    print eTimes, iTimes

    return ax



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

    


###############################################################################
roots = NoiseDataSpaces.Roots(bumpDataRoot, velDataRoot, gridsDataRoot)
ps    = NoiseDataSpaces(roots, shape, noise_sigmas)

tLimits = [2e3, 2.25e3] # ms

rasterFigSize = (3.5, 1.9)
transparent   = True
rasterLeft    = 0.2
rasterBottom  = 0.1
rasterRight   = 0.95
rasterTop     = 0.8

ylabelPos   = -0.2


if (rasters):
    # noise_sigma = 0 pA
    fig = plt.figure("rasters0", figsize=rasterFigSize)
    ax = fig.add_axes(Bbox.from_extents(rasterLeft, rasterBottom, rasterRight,
        rasterTop))
    rasterPlot(ps, 
            noise_sigma_idx=0,
            r=rasterRC[0][0], c=rasterRC[0][1],
            tLimits=tLimits,
            ann_EI=True)
    fname = outputDir + "/figure4_raster0.png"
    fig.savefig(fname, dpi=300, transparent=transparent)
    plt.close()
        

    # noise_sigma = 150 pA
    fig = plt.figure("rasters150", figsize=rasterFigSize)
    ax = fig.add_axes(Bbox.from_extents(rasterLeft, rasterBottom, rasterRight,
        rasterTop))
    rasterPlot(ps, 
            noise_sigma_idx=1,
            r=rasterRC[1][0], c=rasterRC[1][1],
            tLimits=tLimits,
            ylabel='', yticks=False)
    fname = outputDir + "/figure4_raster150.png"
    fig.savefig(fname, dpi=300, transparent=transparent)
    plt.close()
        

    # noise_sigma = 300 pA
    fig = plt.figure("rasters300", figsize=rasterFigSize)
    ax = fig.add_axes(Bbox.from_extents(rasterLeft, rasterBottom, rasterRight,
        rasterTop))
    rasterPlot(ps, 
            noise_sigma_idx=2,
            r=rasterRC[2][0], c=rasterRC[2][1],
            tLimits=tLimits,
            ylabel='', yticks=False)
    fname = outputDir + "/figure4_raster300.png"
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
        fig = plt.figure(figsize=rateFigSize)
        ax = fig.add_axes(Bbox.from_extents(rateLeft, rateBottom, rateRight,
            rateTop))
        kw = {}
        if (idx != 0):
            kw['ylabel'] = ''

        plotAvgFiringRate(ps, noise_sigma_idx=idx,
                r=rasterRC[idx][0], c=rasterRC[idx][1],
                ax=ax, **kw)
        fname = outputDir + "/figure4_rate{0}.pdf".format(noise_sigma)
        fig.savefig(fname, dpi=300, transparent=transparent)
        plt.close()


##############################################################################
boxFigSize = (rasterFigSize[0]/2.0, 1.5)
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

