#!/usr/bin/env python
#
#   figure1.py
#
#   Noise publication Figure 2.
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
from matplotlib.pyplot   import figure, subplot, plot, savefig
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator, AutoMinorLocator, LinearLocator, MaxNLocator, \
        ScalarFormatter
from matplotlib.colorbar import make_axes
from matplotlib.transforms import Bbox

from parameters  import JobTrialSpace2D
from EI_plotting import plotBumpSigmaTrial, computeYX, aggregate2DTrial, \
        aggregate2D, drawEIRectSelection
from plotting.grids import plotGridRateMap, plotAutoCorrelation, plotSpikes2D
from plotting.global_defs import globalAxesSettings, createColorbar
from figures_shared import plotOneHist

import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)

from matplotlib import rc
rc('pdf', fonttype=42)
rc('mathtext', default='regular')

plt.rcParams['font.size'] = 11

outputDir = "."

NTrials=10
iterList  = ['g_AMPA_total', 'g_GABA_total']

noise_sigmas = [0, 150, 300]
exampleIdx   = [(1, 2), (11, 10), (0, 5)] # (row, col)
bumpDataRoot= 'output_local/one_to_one'
velDataRoot = 'output_local/velocity'
bumpShape = (40, 40)

sweep0        = 1
sweep150      = 1
sweep300      = 1

##############################################################################

def getNoiseRootDir(prefix, noise_sigma):
    return  "{0}/EI_param_sweep_{1}pA".format(prefix, int(noise_sigma))


def getNoiseRoots(prefix, noise_sigmas):
    roots = []
    for s in noise_sigmas:
        roots.append(getNoiseRootDir(prefix, s))
    return roots



def drawSweeps(ax, dataSpace, iterList, NTrials=1, r=0, c=0, yLabelOn=True,
        yticks=True, exColor='white',
        cbar=False):
    xLabelText = '$w_I$ (nS)'
    if (yLabelOn):
        yLabelText = '$w_E$ (nS)'
    else:
        yLabelText = ''

    if (ax is None):
        ax = plt.gca()

    varList = ['bump_e', 'sigma']
    G = plotBumpSigmaTrial(dataSpace, varList, iterList,
            trialNumList=range(NTrials),
            xlabel=xLabelText,
            ylabel=yLabelText,
            colorBar=False,
            clBarLabel = "Gridness score",
            clbarNTicks=3,
            yticks=yticks,
            vmin=0,
            vmax=10)
    plt.set_cmap('jet_r')
    cax, kw = make_axes(ax, orientation='vertical', shrink=0.8,
            pad=0.05)
    globalAxesSettings(cax)
    cb = plt.colorbar(ax=ax, cax=cax, ticks=MultipleLocator(5), **kw)
    cb.set_label('Bump $\sigma$ (neurons)')
    if (cbar == False):
        cax.set_visible(False)

    return ax, cax


###############################################################################
bumpRoots = getNoiseRoots(bumpDataRoot, noise_sigmas)
bumpDataSpace0   = JobTrialSpace2D(bumpShape, bumpRoots[0])
bumpDataSpace150 = JobTrialSpace2D(bumpShape, bumpRoots[1])
bumpDataSpace300 = JobTrialSpace2D(bumpShape, bumpRoots[2])

#velRoots = getNoiseRoots(velDataRoot, noise_sigmas)
#velDataSpace0   = JobTrialSpace2D(shape, velRoots[0])
#velDataSpace150 = JobTrialSpace2D(shape, velRoots[1])
#velDataSpace300 = JobTrialSpace2D(shape, velRoots[2])

exSz = 4
exMargin = 0.2
exGsCoords = exMargin, exMargin, 1.0, 1.0

sweepFigSize = (3.4, 1.8)
sweepLeft   = 0.15
sweepBottom = 0.2
sweepRight  = 0.9
sweepTop    = 0.99

if (sweep0):
    # noise_sigma = 0 pA
    fig = figure("sweeps0", figsize=sweepFigSize)
    exRows = [28, 15]
    exCols = [3, 15]
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    ax, cax = drawSweeps(ax, bumpDataSpace0, iterList, NTrials=NTrials,
            cbar=False)
    fname = outputDir + "/figure1_sweeps0.png"
    fig.savefig(fname, dpi=300, transparent=True)



if (sweep150):
    # noise_sigma = 150 pA
    fig = figure("sweeps150", figsize=sweepFigSize)
    exRows = [8, 2]
    exCols = [10, 9]
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    ax, cax = drawSweeps(ax, bumpDataSpace150, iterList, NTrials=NTrials,
            yLabelOn=False, yticks=False) 
    fname = outputDir + "/figure1_sweeps150.png"
    fig.savefig(fname, dpi=300, transparent=True)



if (sweep300):
    # noise_sigma = 300 pA
    fig = figure("sweeps300", figsize=sweepFigSize)
    exRows = [16, 15]
    exCols = [6, 23]
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    ax.set_clip_on(False)
    _, cax = drawSweeps(ax, bumpDataSpace300, iterList, NTrials=NTrials,
            yLabelOn=False, yticks=False, exColor='black', cbar=True)
    fname = outputDir + "/figure1_sweeps300.png"
    fig.savefig(fname, dpi=300, transparent=True)


