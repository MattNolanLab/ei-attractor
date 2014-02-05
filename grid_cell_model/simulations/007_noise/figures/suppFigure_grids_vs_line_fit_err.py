#!/usr/bin/env python
#
#   suppFigure_grids_vs_line_fit_err.py
#
#   Supplementary figure: scatter plots of gridness score vs. line fit error
#
import sys
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.ticker as ti
from copy import deepcopy

from EI_plotting          import scatter
from EI_plotting.base     import NoiseDataSpaces
from submitting import flagparse

parser = flagparse.FlagParser()
parser.add_flag('--scatterPlot')
args = parser.parse_args()


from matplotlib import rc
rc('pdf', fonttype=42)
rc('mathtext', default='regular')
rc('font', size=11)

outputDir = "output_figures"

NTrials=3
iterList  = ['g_AMPA_total', 'g_GABA_total']

noise_sigmas  = [0, 150, 300]
exampleIdx    = [(1, 22), (1, 22), (1, 22)] # (row, col)
bumpDataRoot  = None
velDataRoot   = 'output_local/even_spacing/velocity'
gridsDataRoot = 'output_local/even_spacing/grids'
shape = (31, 31)


roots = NoiseDataSpaces.Roots(bumpDataRoot, velDataRoot, gridsDataRoot)
ps    = NoiseDataSpaces(roots, shape, noise_sigmas)


##############################################################################
# Scatter plot of gridness score vs. bump sigma^-1
scatterFigSize = (8.27, 11.69)
ignoreNaNs = True

if args.scatterPlot or args.all:
    fname = outputDir + "/suppFigure_grids_vs_line_fit_err.pdf"



    gridNTrials = 3
    velNTrials = None
    typesSlope = ['velocity', 'fitErr']
    typesGrids = ['grids', 'gridnessScore']

    xlabel = 'Bump speed fit error (neurons/s/trial)'
    ylabel = 'Gridness score'

    fig = plt.figure(figsize=scatterFigSize)
    scatterPlot = scatter.FullScatterPlot(
            ps.v, ps.grids, typesSlope, typesGrids, iterList, velNTrials,
            gridNTrials,
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
    scatterPlot.plot(captionLeft=-0.1)
    scatterPlot.set_titleSizes(16)

    fig.savefig(fname, dpi=300)

