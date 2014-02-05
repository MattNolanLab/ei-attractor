#!/usr/bin/env python
#
#   suppFigure_grid_sweeps.py
#
#   Supplementary figure: grid field examples
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
import string
import matplotlib.pyplot as plt
import matplotlib.ticker as ti
from matplotlib.backends.backend_pdf import PdfPages

from EI_plotting      import sweeps, examples
from EI_plotting      import aggregate as aggr
from EI_plotting.base import getNoiseRoots
from parameters       import JobTrialSpace2D
from submitting import flagparse

parser = flagparse.FlagParser()
args = parser.parse_args()

from matplotlib import rc
rc('pdf', fonttype=42)
rc('mathtext', default='regular')

plt.rcParams['font.size'] = 12

outputDir = "output_figures"

NTrials=10
YXRC   = [(1, 22), (1, 22), (1, 22)] # (row, col)
iterList  = ['g_AMPA_total', 'g_GABA_total']

noise_sigmas = [0, 150, 300]
gridsDataRoot= 'output_local/even_spacing/grids'
shape = (31, 31)


##############################################################################

def drawSweep(ax, dataSpace, iterList, spaceRect, exIdx=(0, 0)):
    cbar_kw= {'label' : 'Gridness score',
        'orientation': 'vertical',
        'shrink': 0.8,
        'pad' : -0.05,
        'ticks' : ti.MultipleLocator(0.5),
        'rasterized' : True}

    exRow, exCol = exIdx
    Y, X = aggr.computeYX(dataSpace, iterList, r=exRow, c=exCol)


    # Replot gridness score param. sweep
    ax0 = plt.gca()
    G = sweeps.plotGridTrial(dataSpace, ['gridnessScore'], iterList,
            noise_sigma=None, sigmaTitle=False,
            trialNumList=range(NTrials),
            ax=ax0,
            r=exRow,
            c=exCol,
            xlabel="$w_I$ (nS)",
            ylabel="$w_E$ (nS)",
            cbar=True, cbar_kw=cbar_kw,
            vmin=None,
            vmax=None)

    examples.drawEIRectSelection(ax0, spaceRect, X, Y)



def drawA4RectExamples(dataSpace, noise_sigma, iterList, exRect, exIdx,
        letter=''):
    fig = plt.figure(figsize=(8.27, 11.69))
    margin    = 0.1
    sw_left   = margin
    sw_bottom = 0.82
    sw_right  = 0.4
    sw_top    = 0.96
    div       = 0.075

    letter_left = 0.03
    letter_top_off = 0.01
    letter_va='bottom'
    letter_ha='left'

    nsX = 0.9
    nsY = sw_top

    sweepsRect = sw_left, sw_bottom, sw_right-sw_left, sw_top-sw_bottom
    ax_sweeps = fig.add_axes(sweepsRect)
    drawSweep(ax_sweeps, dataSpace, iterList, exRect, exIdx)
    fig.text(letter_left, sw_top+letter_top_off, letter, va=letter_va, ha=letter_ha,
            fontsize=19, fontweight='bold')

    gsCoords = 0.12, 0.075, 0.92, sw_bottom-div
    #gsCoords = margin, 0.46, 0.5, sw_bottom-div
    gs = examples.drawGridExamples(dataSpace, exRect, iterList, gsCoords=gsCoords,
            exIdx=exIdx, fig=fig)
    #fig.text(letter_left, sw_bottom-div+letter_top_off, "B", va=letter_va,
    #        ha=letter_ha, fontsize=19, fontweight='bold')
    noise_sigma_txt = "$\sigma_{{noise}}$ = {0} pA".format(int(noise_sigma))
    fig.text(nsX, nsY, noise_sigma_txt, va='center', ha='right', fontsize=19)


gridRoots = getNoiseRoots(gridsDataRoot, noise_sigmas)
gridSpace0   = JobTrialSpace2D(shape, gridRoots[0])
gridSpace150 = JobTrialSpace2D(shape, gridRoots[1])
gridSpace300 = JobTrialSpace2D(shape, gridRoots[2])
gridSpaces = [gridSpace0, gridSpace150, gridSpace300]

exWidth = 6
exHeight = 8

exampleRC = [
        [[1, 20], [20, 10], [3, 6]],
        [[5, 2],  [15, 15], [1, 15]],
        [[3, 7], [15, 15], [1, 20]]]
enable = [
        [1, 1, 1],
        [1, 1, 1],
        [1, 1, 1]]




fname = outputDir + "/suppFigure_grid_sweeps.pdf"
outputPDF = PdfPages(fname)
strIdx = 0
for noise_idx, noise_sigma in enumerate(noise_sigmas):
    for exampleIdx, RC in enumerate(exampleRC[noise_idx]):
        print noise_idx, exampleIdx, RC

        exLeft   = RC[0]
        exBottom = RC[1]
        exRect = [exLeft, exBottom, exLeft+exWidth-1, exBottom+exHeight-1]
        drawA4RectExamples(gridSpaces[noise_idx], noise_sigma, iterList,
                exRect, YXRC[noise_idx], letter=string.ascii_uppercase[strIdx])
        
        outputPDF.savefig(dpi=300, transparent=False)
        plt.close()
        strIdx += 1

outputPDF.close()

