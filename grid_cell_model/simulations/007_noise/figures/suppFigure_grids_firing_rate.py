#!/usr/bin/env python
#
#   suppFigure_grids_firing_rate.py
#
#   Supplementary figure: average population firing rates.
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
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.ticker as ti
from matplotlib.transforms import Bbox
from copy import deepcopy

from EI_plotting import sweeps
from EI_plotting import aggregate as aggr
from figures_shared       import NoiseDataSpaces
from parameters           import JobTrialSpace2D

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
bumpDataRoot  = None
velDataRoot   = None
gridsDataRoot = 'output_local/even_spacing/grids'
shape = (31, 31)


FRSweeps    = 0
scatterPlot = 1

roots = NoiseDataSpaces.Roots(bumpDataRoot, velDataRoot, gridsDataRoot)
ps    = NoiseDataSpaces(roots, shape, noise_sigmas)


##############################################################################
# Parameter sweeps of E and I firing rates
if FRSweeps:
    varList_e = ['FR_e',  'avg']
    varList_i = ['FR_i',  'all']

    sweepFigSize = (3.5, 2.5)
    sweepLeft   = 0.15
    sweepBottom = 0.2
    sweepRight  = 0.87
    sweepTop    = 0.85
    transparent  = True

    FR_e_cbar_kw= {
        'label'      : 'E Firing rate (Hz)',
        'location'   : 'right',
        'shrink'     : 0.8,
        'pad'        : -0.05,
        'ticks'      : ti.MultipleLocator(2),
        'rasterized' : True,
        'extend'     : 'max',
        'extendfrac' : 0.1}

    FR_i_cbar_kw = deepcopy(FR_e_cbar_kw)
    FR_i_cbar_kw.update(dict(
        label='I Firing rate (Hz)',
        ticks=ti.MultipleLocator(50)))

    FR_e_vmin = 0
    FR_e_vmax = 5
    FR_i_vmin = 0
    FR_i_vmax = 100


    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        if (ns_idx == 2):
            cbar = True
        else:
            cbar = False
        kw = {}
        if (ns_idx != 0):
            kw['ylabel'] = ''
            kw['yticks'] = False

        # E cells
        fig = plt.figure(figsize=sweepFigSize)
        ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
            sweepTop))
        sweeps.plotFRTrial(ps.grids[ns_idx], varList_e, iterList,
                noise_sigma=noise_sigma,
                trialNumList=gridTrialNumList,
                ax=ax,
                xlabel='', xticks=False,
                r=exampleIdx[ns_idx][0], c=exampleIdx[ns_idx][1],
                cbar=cbar, cbar_kw=FR_e_cbar_kw,
                vmin=FR_e_vmin, vmax=FR_e_vmax,
                ignoreNaNs=False,
                **kw)
        fname = outputDir + "/suppFigure_grids_FR_E{0}.pdf"
        fig.savefig(fname.format(int(noise_sigma)), dpi=300,
                transparent=transparent)
        plt.close()

        # I cells
        fig = plt.figure(figsize=sweepFigSize)
        ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
            sweepTop))
        sweeps.plotFRTrial(ps.grids[ns_idx], varList_i, iterList,
                noise_sigma=noise_sigma,
                trialNumList=gridTrialNumList,
                ax=ax,
                r=exampleIdx[ns_idx][0], c=exampleIdx[ns_idx][1],
                cbar=cbar, cbar_kw=FR_i_cbar_kw,
                vmin=FR_i_vmin, vmax=FR_i_vmax,
                sigmaTitle=False,
                ignoreNaNs=False,
                **kw)
        fname = outputDir + "/suppFigure_grids_FR_I{0}.pdf"
        fig.savefig(fname.format(int(noise_sigma)), dpi=300,
                transparent=transparent)
        plt.close()


##############################################################################
# Scatter plot of gridness score vs. gamma power
scatterFigSize = (6.5, 3.6)
scatterLeft   = 0.15
scatterBottom = 0.2
scatterRight  = 0.98
scatterTop    = 0.95
scatterTransparent = True

ignoreNaNs = True

if (scatterPlot):
    for EIType in ['E', 'I_10']:
        if EIType == 'E':
            fname = outputDir + "/suppFigure_grids_FR-scatter_FRE_vs_grids.pdf"
            xlabel='Mean firing rate of E cells (Hz)'
            legLoc = (0.8, 0.6)
        else:
            fname = outputDir + "/suppFigure_grids_FR-scatter_FRI_vs_grids.pdf"
            xlabel='Mean firing rate of I cells (Hz)'
            legLoc = (0.1, 0.6)

        fig = plt.figure(figsize=scatterFigSize)
        ax = fig.add_axes(Bbox.from_extents(scatterLeft, scatterBottom, scatterRight,
            scatterTop))

        NTrialsGrids = 3
        typesFR = ['FR', EIType]
        typesGrids = ['grids', 'gridnessScore']
        ax.hold('on')
        scatterColors = ['blue', 'green', 'red']

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            color = scatterColors[ns_idx]
            scatterPlot = sweeps.ScatterPlot(
                    ps.grids[ns_idx], ps.grids[ns_idx], typesFR,
                    typesGrids, iterList, NTrialsGrids, NTrialsGrids,
                    c=color,
                    s=15,
                    linewidth=0.3,
                    xlabel=xlabel,
                    ylabel='Gridness score',
                    ignoreNaNs=ignoreNaNs)
            scatterPlot.plot()
        ax.xaxis.set_major_locator(ti.MultipleLocator(10))
        ax.yaxis.set_major_locator(ti.MultipleLocator(0.5))
        ax.yaxis.set_minor_locator(ti.MultipleLocator(0.25))
        ax.set_xscale('log')
        if EIType != 'E':
            ax.set_xlim([2, 300])
        leg = ['0', '150', '300']
        l = ax.legend(leg, loc=legLoc, fontsize='small', frameon=False,
                numpoints=1, title='$\sigma$ (pA)')
        plt.setp(l.get_title(), size='small')

        fig.savefig(fname, dpi=300, transparent=scatterTransparent)

