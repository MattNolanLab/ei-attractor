#!/usr/bin/env python
#
#   suppFigure_grids_vs_bumps.py
#
#   Supplementary figure: scatter plots of gridness score vs. bump width.
#
import sys
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.ticker as ti
from matplotlib.transforms import Bbox
from copy import deepcopy

from EI_plotting import sweeps
from EI_plotting import aggregate as aggr
from figures_shared       import NoiseDataSpaces
from parameters           import JobTrialSpace2D
from plotting.global_defs import globalAxesSettings

import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)

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
bumpDataRoot= 'output_local/even_spacing/gamma_bump'
velDataRoot   = None
gridsDataRoot = 'output_local/even_spacing/grids'
shape = (31, 31)


scatterPlot  = 1

roots = NoiseDataSpaces.Roots(bumpDataRoot, velDataRoot, gridsDataRoot)
ps    = NoiseDataSpaces(roots, shape, noise_sigmas)


##############################################################################
# Scatter plot of gridness score vs. bump sigma^-1
scatterFigSize = (6.5, 3.6)
scatterLeft   = 0.12
scatterBottom = 0.15
scatterRight  = 0.98
scatterTop    = 0.92
scatterTransparent = True

scatterColorFigSize = (1.5, 1.5)

ignoreNaNs = True

if (scatterPlot):
    fname = outputDir + "/suppFigure_grids_vs_bumps{0}.pdf"

    ylabel = '$\sigma_{bump}^{-1}\ (neurons^{-1})$'

    fig = plt.figure(figsize=scatterFigSize)
    ax = fig.add_axes(Bbox.from_extents(scatterLeft, scatterBottom, scatterRight,
        scatterTop))

    NTrialsBumps = 5
    NTrialsGrids = 3
    typesBumps = ['bump', 'sigma']
    typesGrids = ['grids', 'gridnessScore']
    ax.hold('on')

    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        if (ns_idx == 2):
            xlabel = 'Gridness score'
        else:
            xlabel = ''

        fig = plt.figure(figsize=scatterFigSize)
        ax = fig.add_axes(Bbox.from_extents(scatterLeft, scatterBottom,
            scatterRight, scatterTop))

        scatterPlot = sweeps.ScatterPlot(
                ps.grids[ns_idx], ps.bumpGamma[ns_idx], typesGrids,
                typesBumps, iterList, NTrialsBumps, NTrialsGrids,
                s=25,
                linewidth=0.3,
                color2D=True,
                xlabel=xlabel,
                ylabel=ylabel,
                sigmaTitle=True,
                noise_sigma=noise_sigma,
                ignoreNaNs=ignoreNaNs,
                ax=ax)
        scatterPlot.plot()

        #ax.xaxis.set_major_locator(ti.MultipleLocator(10))
        #ax.yaxis.set_major_locator(ti.MultipleLocator(0.5))
        #ax.yaxis.set_minor_locator(ti.MultipleLocator(0.25))
        ax.set_title(ax.get_title(), size='16')
        ax.axhline(y=0.1, linestyle=':', zorder=-1, color='black')

        fig.savefig(fname.format(int(noise_sigma)), dpi=300,
                transparent=scatterTransparent)
        #plt.show()

    fig = plt.figure(figsize=scatterColorFigSize)
    ax = fig.gca()
    scatterPlot.plotColorbar(ax)
    fig.tight_layout(pad=0)
    fname = outputDir + "/suppFigure_grids_vs_bumps_colorbar.pdf"
    fig.savefig(fname, dpi=300, transparent=scatterTransparent)


