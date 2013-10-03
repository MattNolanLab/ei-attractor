#!/usr/bin/env python
#
#   figure2.py
#
#   Noise publication: slices through the parameter region.
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
from matplotlib.gridspec   import GridSpec
from matplotlib.ticker     import MultipleLocator, AutoMinorLocator, \
        MaxNLocator
from matplotlib.transforms import Bbox

import EI_plotting as EI
from parameters           import JobTrialSpace2D
from EI_plotting          import aggregate2DTrial, aggregate2D, computeYX, \
        computeVelYX
from plotting.grids       import plotGridRateMap, plotAutoCorrelation, plotSpikes2D
from plotting.global_defs import globalAxesSettings
from figures_shared       import NoiseDataSpaces

import logging as lg
lg.basicConfig(level=lg.WARN)
#lg.basicConfig(level=lg.INFO)

from matplotlib import rc
rc('pdf', fonttype=42)
rc('mathtext', default='regular')

plt.rcParams['font.size'] = 10

outputDir = "."

NTrials=10
iterList  = ['g_AMPA_total', 'g_GABA_total']

noise_sigmas = [0, 150, 300]
rcIdx   = [(0, 0), (0, 0), (0, 0)] # (row, col)
gridsDataRoot= 'output_local/even_spacing/grids'
bumpDataRoot= 'output_local/even_spacing/gamma_bump'
velDataRoot = 'output_local/even_spacing/velocity'
shape = (31, 31)


##############################################################################

def plotOneSlice(ax, x, y, **kw):
    xlabel    = kw.pop('xlabel', '')
    ylabel    = kw.pop('ylabel', '')
    ylabelPos = kw.pop('ylabelPos', -0.2)
    fmt       = kw.pop('fmt', 'o-')
    xticks    = kw.pop('xticks', True)
    yticks    = kw.pop('yticks', True)

    globalAxesSettings(ax)
    ndim = len(y.shape)
    if (ndim == 2):
        mean = np.mean(y, axis=1) # axis 1: trials
        std  = np.std(y, axis=1)
        ax.errorbar(x, mean, std, fmt=fmt, **kw)
    elif (ndim == 1):
        ax.plot(x, y, fmt, **kw)
    ax.set_xlabel(xlabel)
    ax.text(ylabelPos, 0.5, ylabel, rotation=90, transform=ax.transAxes,
            va='center', ha='center')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(3))
    w = x[-1] - x[0]
    margin = 0.025
    ax.set_xlim([-margin*w, x[-1]+margin*w])

    if (not xticks):
        ax.xaxis.set_ticklabels([])
    if (not yticks):
        ax.yaxis.set_ticklabels([])


def extractSliceX(X, Y, rowSlice, colSlice):
    if (isinstance(rowSlice, slice)):
        return Y[rowSlice, colSlice]
    else:
        return X[rowSlice, colSlice]


def plotSlice(paramSpaces, rowSlice, colSlice, **kw):
    # kwargs
    title  = kw.pop('title', True)

    if (isinstance(rowSlice, slice)):
        xlabel = EI.ylabelText
        titleText = "$g_I$"
    else:
        xlabel = EI.xlabelText
        titleText = "$g_E$"

    GVars         = ['gridnessScore']
    gammaVars     = ['acVal']
    bumpSigmaVars = ['sigma']
    errVars       = ['lineFitErr']
    slopeVars     = ['lineFitSlope']
    gridsTrialNumList = range(1)
    gammaTrialNumList = range(5)

    rcG     = [(1, 22), (1, 22), (1, 22)] # (row, col)
    rcGamma = [(0, 0), (0, 0), (0, 0)] # (row, col)

    sp = paramSpaces

    spRows = 5
    spCols = 1

    # Gridness score
    ax_G = plt.subplot(spRows, spCols, 1)
    for idx, noise_sigma in enumerate(sp.noise_sigmas):
        space = sp.grids[idx]
        G = space.aggregateData(GVars, gridsTrialNumList,
                output_dtype='array', loadData=True, saveData=False)
        Y, X = computeYX(space, iterList, r=rcG[idx][0], c=rcG[idx][1])
        x = extractSliceX(X, Y, rowSlice, colSlice)
        y = G[rowSlice, colSlice, :]
        kw['xlabel'] = ''
        kw['ylabel'] = 'Gridness score'
        plotOneSlice(ax_G, x, y, xticks=False, **kw)
    ax_G.yaxis.set_major_locator(MaxNLocator(4))
        
    
    # Gamma
    ax_gamma = plt.subplot(spRows, spCols, 2)
    for idx, noise_sigma in enumerate(sp.noise_sigmas):
        space = sp.bumpGamma[idx]
        P = space.aggregateData(gammaVars, gammaTrialNumList,
                output_dtype='array', loadData=True, saveData=False)
        Y, X = computeYX(space, iterList, r=rcGamma[idx][0], c=rcGamma[idx][1])
        x = extractSliceX(X, Y, rowSlice, colSlice)
        y = P[rowSlice, colSlice, :]
        kw['xlabel'] = ''
        kw['ylabel'] = 'Correlation ($\gamma$)'
        plotOneSlice(ax_gamma, x, y, xticks=False, **kw)
    ax_gamma.yaxis.set_major_locator(MaxNLocator(4))
    

    # Bump sigma
    ax_sigma = plt.subplot(spRows, spCols, 3)
    for idx, noise_sigma in enumerate(sp.noise_sigmas):
        space = sp.bumpGamma[idx]
        sigma = space.aggregateData(bumpSigmaVars, gammaTrialNumList,
                output_dtype='array', loadData=True, saveData=False)
        Y, X = computeYX(space, iterList, r=rcGamma[idx][0], c=rcGamma[idx][1])
        x = extractSliceX(X, Y, rowSlice, colSlice)
        y = sigma[rowSlice, colSlice, :]
        kw['xlabel'] = ''
        kw['ylabel'] = 'Bump $\sigma$ (neurons)'
        plotOneSlice(ax_sigma, x, y, xticks=False, **kw)
    ax_sigma.set_yscale('log')


    # Vel error of fit
    ax_velErr = plt.subplot(spRows, spCols, 4)
    for idx, noise_sigma in enumerate(sp.noise_sigmas):
        space = sp.v[idx]
        V = space.aggregateData(errVars, 'all-at-once', output_dtype='array',
                loadData=True, saveData=False)
        Y, X = computeVelYX(space, iterList)
        x = extractSliceX(X, Y, rowSlice, colSlice)
        y = V[rowSlice, colSlice]
        kw['xlabel'] = ''
        kw['ylabel'] = 'Velocity error of fit\n(neurons/s)'
        plotOneSlice(ax_velErr, x, y, xticks=False, **kw)
    ax_velErr.yaxis.set_major_locator(MaxNLocator(4))

    # Vel slope
    ax_velSlope = plt.subplot(spRows, spCols, 5)
    for idx, noise_sigma in enumerate(sp.noise_sigmas):
        space = sp.v[idx]
        S = space.aggregateData(slopeVars, 'all-at-once', output_dtype='array',
                loadData=True, saveData=False)
        Y, X = computeVelYX(space, iterList)
        x = extractSliceX(X, Y, rowSlice, colSlice)
        y = S[rowSlice, colSlice]
        kw['xlabel'] = xlabel
        kw['ylabel'] = 'Bump speed slope\n(neurons/s/pA)'
        plotOneSlice(ax_velSlope, x, y, **kw)
    ax_velSlope.yaxis.set_major_locator(MaxNLocator(4))

    l_ax = ax_G
    leg = []
    for s in sp.noise_sigmas:
        leg.append("{0}".format(int(s)))
    l = l_ax.legend(leg, loc=(0.8, 0.9), title='$\sigma$ (pA)',
            frameon=False, fontsize='small', ncol=1)
    plt.setp(l.get_title(), fontsize='small')

    # Title
    if (title):
        if (isinstance(rowSlice, slice)):
            sliceG = X[0, colSlice]
        else:
            sliceG = Y[rowSlice, 0]
        ax_G.set_title('{0} = {1} nS'.format(titleText, sliceG), y=1.2,
                va='bottom', fontsize='large')
        

        


###############################################################################
roots = NoiseDataSpaces.Roots(bumpDataRoot, velDataRoot, gridsDataRoot)
ps    = NoiseDataSpaces(roots, shape, noise_sigmas)


figSize = (4, 10)
left   = 0.15
bottom = 0.2
right  = 0.87
top    = 0.85

outputDir = 'slices'

for idx in xrange(shape[0]):
    print "Horizontal slice no. {0}".format(idx)
    fig = plt.figure("slices_horizontal", figsize=figSize)
    plotSlice(ps, idx, slice(None))
    fig.tight_layout(rect=[0.1, 0.01, 0.99, 0.99], pad=0.2, h_pad=1)
    fname = "{0}/figure_slices_horizontal_{1}.png".format(outputDir, idx)
    plt.savefig(fname, dpi=300, transparent=False)
    plt.close()

for idx in xrange(shape[1]):
    print "Vertical slice no. {0}".format(idx)
    fig = plt.figure("slices_vertical", figsize=figSize)
    plotSlice(ps, slice(None), idx)
    fig.tight_layout(rect=[0.1, 0.01, 0.99, 0.99], pad=0.2, h_pad=1)
    fname = "{0}/figure_slices_vertical_{1}.png".format(outputDir, idx)
    plt.savefig(fname, dpi=300, transparent=False)
    plt.close()
