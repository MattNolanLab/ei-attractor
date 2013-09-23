#!/usr/bin/env python
#
#   figure4.py
#
#   Putting it all together: parameterscape plots.
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
from matplotlib.pyplot   import figure, plot, savefig
from matplotlib.transforms import Bbox
from matplotlib.ticker import MaxNLocator, MultipleLocator

import numpy.ma as ma

from parameters  import JobTrialSpace2D, DataSpace
from EI_plotting import aggregate2DTrial, computeYX
from plotting.global_defs import globalAxesSettings
from plotting.parameterscape import parameterScape
from figures_shared import getNoiseRoots, plotOneHist
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
exampleIdx    = [(0, 0), (0, 0), (0, 0)] # (row, col)
gammaDataRoot = 'output_local/even_spacing/gamma_bump'
velDataRoot = 'output_local/velocity'
gridsDataRoot= 'output_local/grids_50pA_Ivel'
shape = (31, 31)

figSize = (7.5, 9)


###############################################################################
noise0   = 1
noise150 = 1
noise300 = 1

transparent = False
dpi = 300


###############################################################################
def drawParamScape(fig, dataSpace, varList, trialNumList, noise_sigma):
    # Bump sigmas
    sigmaVarList = ['bump_e', 'sigma']
    CVarList = ['acVal']
    FVarList = ['freq']
    sigma = aggregate2DTrial(dataSpace, sigmaVarList, trialNumList)
    C = aggregate2DTrial(dataSpace, CVarList, trialNumList)
    F = aggregate2DTrial(dataSpace, FVarList, trialNumList)
    sigma = ma.MaskedArray(sigma, mask=np.isnan(sigma))
    C = ma.MaskedArray(C, mask=np.isnan(C))
    F = ma.MaskedArray(F, mask=np.isnan(F))
    Y, X = computeYX(dataSpace, iterList)
    Y = Y[:, 0]
    X = X[0, :]
    legend = ['Bump sigma (neurons)', 'Correlatin', 'Frequency (Hz)']
    legend_kw = dict(x=-0.02, y=1.07)
    tickLocators = [MultipleLocator(10), MultipleLocator(0.2),
            MultipleLocator(30)]
    colorbars = dict(
            labels=['Bump $\sigma$ (neurons)',
                'Correlation',
                'Oscillation frequency (Hz)'],
            orientation='horizontal',
            labelsize='small',
            x = 0.32,
            y = 0.85,
            div = 0.06,
            width = 0.3,
            tickLocators=tickLocators)

    ax, cax = parameterScape(fig, X, Y, [sigma, C, F],
            fmts=['o', 's', 's'],
            sizes=[25, 12, 6],
            cmaps=['jet_r', 'jet', 'jet'],
            vranges=[(0, 10), (-0.09, 0.675), (30, 90)],
            legend=(legend, legend_kw),
            xmargin=1.2,
            ymargin=1.2,
            colorbars=colorbars)
    ax.set_xlabel('$g_I$ (nS)')
    ax.set_ylabel('$g_E$ (nS)')

    main_top = 0.8
    nsX = 0.93
    nsY = main_top + 0.02
    ax.set_position(Bbox.from_extents([0.025, 0.07, 0.99, main_top]))
    noise_sigma_txt = "$\sigma_{{noise}}$ = {0} pA".format(int(noise_sigma))
    fig.text(nsX, nsY, noise_sigma_txt, va='center', ha='right', fontsize=19)




###############################################################################

gammaRoots = getNoiseRoots(gammaDataRoot, noise_sigmas)
gammaDataSpace0   = JobTrialSpace2D(shape, gammaRoots[0])
gammaDataSpace150 = JobTrialSpace2D(shape, gammaRoots[1])
gammaDataSpace300 = JobTrialSpace2D(shape, gammaRoots[2])

velRoots = getNoiseRoots(velDataRoot, noise_sigmas)
velDataSpace0   = JobTrialSpace2D(shape, velRoots[0])
velDataSpace150 = JobTrialSpace2D(shape, velRoots[1])
velDataSpace300 = JobTrialSpace2D(shape, velRoots[2])

if (noise0):
    fig = plt.figure(figsize=figSize)
    drawParamScape(fig, gammaDataSpace0, varList=None,
            trialNumList=xrange(NTrials), noise_sigma=noise_sigmas[0])
    fname = outputDir + "/figure4_parameterscape_0pA.png"
    savefig(fname, dpi=dpi, transparent=transparent)

if (noise150):
    fig = plt.figure(figsize=figSize)
    drawParamScape(fig, gammaDataSpace150, varList=None,
            trialNumList=xrange(NTrials), noise_sigma=noise_sigmas[1])
    fname = outputDir + "/figure4_parameterscape_150pA.png"
    savefig(fname, dpi=dpi, transparent=transparent)

if (noise300):
    fig = plt.figure(figsize=figSize)
    drawParamScape(fig, gammaDataSpace300, varList=None,
            trialNumList=xrange(NTrials), noise_sigma=noise_sigmas[2])
    fname = outputDir + "/figure4_parameterscape_300pA.png"
    savefig(fname, dpi=dpi, transparent=transparent)

