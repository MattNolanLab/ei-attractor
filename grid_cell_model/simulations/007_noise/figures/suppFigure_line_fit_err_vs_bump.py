#!/usr/bin/env python
#
#   suppFigure_line_fit_error_vs_bump.py
#
#   Supplementary figure: scatter plots of line fit error vs. bump sigma
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
bumpDataRoot  = 'output_local/even_spacing/gamma_bump'
velDataRoot   = 'output_local/even_spacing/velocity'
gridsDataRoot = None
shape = (31, 31)


roots = NoiseDataSpaces.Roots(bumpDataRoot, velDataRoot, gridsDataRoot)
ps    = NoiseDataSpaces(roots, shape, noise_sigmas)


##############################################################################
# Scatter plot of gridness score vs. bump sigma^-1
scatterFigSize = (8.27, 11.69)
ignoreNaNs = True

if args.scatterPlot or args.all:
    fname = outputDir + "/suppFigure_line_fit_error_vs_bump.pdf"



    bumpNTrials = 5
    velNTrials = None
    typesBump = ['bump', 'sigma']
    typesFitErr = ['velocity', 'fitErr']

    xlabel = '$\sigma_{bump}^{-1}\ (neurons^{-1})$'
    ylabel = 'Bump speed fit error (neurons/s/trial)'

    fig = plt.figure(figsize=scatterFigSize)
    scatterPlot = scatter.FullScatterPlot(
            ps.bumpGamma, ps.v, typesBump, typesFitErr, iterList, bumpNTrials,
            velNTrials,
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
    scatterPlot.plot()
    scatterPlot.set_titleSizes(16)

    #ax.xaxis.set_major_locator(ti.MultipleLocator(10))
    #ax.yaxis.set_major_locator(ti.MultipleLocator(0.5))
    #ax.yaxis.set_minor_locator(ti.MultipleLocator(0.25))

    fig.savefig(fname, dpi=300)

