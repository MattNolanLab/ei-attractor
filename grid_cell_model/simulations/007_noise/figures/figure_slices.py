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
from grid_cell_model.parameters           import JobTrialSpace2D
from EI_plotting          import aggregate2DTrial, aggregate2D, computeYX, \
        computeVelYX
from plotting.grids       import plotGridRateMap, plotAutoCorrelation, plotSpikes2D
from plotting.global_defs import globalAxesSettings
from EI_plotting.base     import NoiseDataSpaces

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
velDataRoot = 'output_local/even_spacing/velocity_vertical'
shape = (31, 31)


##############################################################################


def plotSlice(paramSpaces, rowSlice, colSlice, **kw):
    # kwargs
    title  = kw.pop('title', True)

    gammaVars     = ['acVal']
    bumpSigmaVars = ['sigma']
    errVars       = ['lineFitErr']
    slopeVars     = ['lineFitSlope']
    gammaTrialNumList = range(5)

    labels = EI.decideLabels(rowSlice, colSlice)

    rcGamma = [(0, 0), (0, 0), (0, 0)] # (row, col)

    sp = paramSpaces

    spRows = 5
    spCols = 1

    # Gridness score
    ax_G = plt.subplot(spRows, spCols, 1)
    EI.plotGridnessSlice(paramSpaces, rowSlice, colSlice, ax=ax_G, title=False)
        
    
    # Gamma
    ax_gamma = plt.subplot(spRows, spCols, 2)
    for idx, noise_sigma in enumerate(sp.noise_sigmas):
        space = sp.bumpGamma[idx]
        P = space.aggregateData(gammaVars, gammaTrialNumList,
                output_dtype='array', loadData=True, saveData=False)
        Y, X = computeYX(space, iterList, r=rcGamma[idx][0], c=rcGamma[idx][1])
        x = EI.extractSliceX(X, Y, rowSlice, colSlice)
        y = P[rowSlice, colSlice, :]
        kw['xlabel'] = ''
        kw['ylabel'] = 'Correlation ($\gamma$)'
        EI.plotOneSlice(ax_gamma, x, y, xticks=False, **kw)
    ax_gamma.yaxis.set_major_locator(MaxNLocator(4))
    

    # Bump sigma
    ax_sigma = plt.subplot(spRows, spCols, 3)
    for idx, noise_sigma in enumerate(sp.noise_sigmas):
        space = sp.bumpGamma[idx]
        sigma = space.aggregateData(bumpSigmaVars, gammaTrialNumList,
                output_dtype='array', loadData=True, saveData=False)
        Y, X = computeYX(space, iterList, r=rcGamma[idx][0], c=rcGamma[idx][1])
        x = EI.extractSliceX(X, Y, rowSlice, colSlice)
        y = sigma[rowSlice, colSlice, :]
        kw['xlabel'] = ''
        kw['ylabel'] = 'Bump $\sigma$ (neurons)'
        EI.plotOneSlice(ax_sigma, x, y, xticks=False, **kw)
    ax_sigma.set_yscale('log')


    # Vel error of fit
    ax_velErr = plt.subplot(spRows, spCols, 4)
    for idx, noise_sigma in enumerate(sp.noise_sigmas):
        space = sp.v[idx]
        V = space.aggregateData(errVars, 'all-at-once', output_dtype='array',
                loadData=True, saveData=False)
        Y, X = computeVelYX(space, iterList)
        x = EI.extractSliceX(X, Y, rowSlice, colSlice)
        y = V[rowSlice, colSlice]
        kw['xlabel'] = ''
        kw['ylabel'] = 'Velocity error of fit\n(neurons/s)'
        EI.plotOneSlice(ax_velErr, x, y, xticks=False, **kw)
    ax_velErr.yaxis.set_major_locator(MaxNLocator(4))

    # Vel slope
    ax_velSlope = plt.subplot(spRows, spCols, 5)
    for idx, noise_sigma in enumerate(sp.noise_sigmas):
        space = sp.v[idx]
        S = space.aggregateData(slopeVars, 'all-at-once', output_dtype='array',
                loadData=True, saveData=False)
        Y, X = computeVelYX(space, iterList)
        x = EI.extractSliceX(X, Y, rowSlice, colSlice)
        y = S[rowSlice, colSlice]
        kw['xlabel'] = labels['xlabel']
        kw['ylabel'] = 'Bump speed slope\n(neurons/s/pA)'
        EI.plotOneSlice(ax_velSlope, x, y, **kw)
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
        ax_G.set_title('{0} = {1} nS'.format(labels['titleText'], sliceG),
                y=1.2, va='bottom', fontsize='large')
        

        


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
