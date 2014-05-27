#!/usr/bin/env python
#
#   suppFigure_grids_vs_bumps.py
#
#   Supplementary figure: scatter plots of gridness score vs. bump width.
#
import math
import sys
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.ticker as ti
from matplotlib.transforms import Bbox
from copy import deepcopy

import default_settings as ds
from EI_plotting          import sweeps, scatter
from EI_plotting          import aggregate as aggr
from EI_plotting.base     import NoiseDataSpaces
from EI_plotting import scaling
from grid_cell_model.parameters           import JobTrialSpace2D
from grid_cell_model.plotting.global_defs import globalAxesSettings
from grid_cell_model.submitting import flagparse

parser = flagparse.FlagParser()
parser.add_flag('--scatterPlot')
parser.add_argument('--expScatter', action='store_true')
args = parser.parse_args()

outputDir = "output_figures"
ps = ds.getDefaultParamSpaces()


##############################################################################
# Scatter plot of gridness score vs. bump sigma^-1
scatterFigSize = (8.27, 11.69)
scatterLeft   = 0.12
scatterBottom = 0.17
scatterRight  = 0.98
scatterTop    = 0.92
scatterTransparent = True

scatterColorFigSize = (1.5, 1.5)

ignoreNaNs = True


##############################################################################
xlabel = 'P(bumps)'
ylabel = 'Gridness score'

if args.scatterPlot or args.all:
    fig = plt.figure(figsize=scatterFigSize)
    ax = fig.add_axes(Bbox.from_extents(scatterLeft, scatterBottom, scatterRight,
        scatterTop))

    isBumpData = []
    gridData = []
    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        isBumpData.append(aggr.IsBump(ps.bumpGamma[ns_idx], ds.iterList,
            ignoreNaNs=True))
        gridData.append(aggr.GridnessScore(ps.grids[ns_idx], ds.iterList,
            ignoreNaNs=True))

    fig = plt.figure(figsize=scatterFigSize)
    scatterPlot = scatter.FullScatterPlot(
            isBumpData, gridData, None, None, ds.iterList, None, None,
            s=25,
            linewidth=0.3,
            color2D=True,
            xlabel=xlabel,
            ylabel=ylabel,
            sigmaTitle=True,
            noise_sigmas=ps.noise_sigmas,
            ignoreNaNs=ignoreNaNs,
            captionLetters=('A', 'B', 'C'),
            fig=fig)
    scatterPlot.plot(captionLeft=-0.1, plotcolorbar=False)
    l = 0.14
    w = 0.165
    scatterPlot.plotColorbar(left=l, bottom=.85, right=l+w, top=.95)
    scatterPlot.set_titleSizes(16)
    if args.expScatter:
        for ns_idx, _ in enumerate(ps.noise_sigmas):
            ax = scatterPlot.axes[ns_idx]
            ax.set_xscale('custom')
            ax.xaxis.set_major_locator(ti.MultipleLocator(.5))
            ax.xaxis.set_minor_locator(ti.MultipleLocator(.1))
            ax.set_xlim([-0.3, 1.002])
        fname = outputDir + "/suppFigure_grids_vs_bumps_exp.pdf"
    else:
        fname = outputDir + "/suppFigure_grids_vs_bumps.pdf"

    fig.savefig(fname, dpi=300)

