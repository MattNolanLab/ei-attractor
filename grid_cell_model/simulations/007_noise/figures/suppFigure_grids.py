#!/usr/bin/env python
#
#   suppFigure_grids.py
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
import matplotlib.pyplot as plt
import matplotlib.ticker as ti

import EI_plotting as EI
from parameters     import JobTrialSpace2D
from figures_shared import getNoiseRoots

import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)

from matplotlib import rc
rc('pdf', fonttype=42)
rc('mathtext', default='regular')

plt.rcParams['font.size'] = 12

outputDir = "suppFigure_grids"

NTrials=10
exampleIdx   = [(1, 22), (1, 22), (1, 22)] # (row, col)
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
    Y, X = EI.computeYX(dataSpace, iterList, r=exRow, c=exCol)


    # Replot gridness score param. sweep
    ax0 = plt.gca()
    G = EI.plotGridTrial(dataSpace, ['gridnessScore'], iterList,
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

    EI.drawEIRectSelection(ax0, spaceRect, X, Y)



def drawA4RectExamples(dataSpace, noise_sigma, iterList, exRect, exIdx):
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
    fig.text(letter_left, sw_top+letter_top_off, "A", va=letter_va, ha=letter_ha,
            fontsize=19, fontweight='bold')

    gsCoords = 0.12, 0.075, 0.92, sw_bottom-div
    #gsCoords = margin, 0.46, 0.5, sw_bottom-div
    gs = EI.drawGridExamples(dataSpace, exRect, iterList, gsCoords=gsCoords,
            exIdx=exIdx)
    fig.text(letter_left, sw_bottom-div+letter_top_off, "B", va=letter_va,
            ha=letter_ha, fontsize=19, fontweight='bold')
    noise_sigma_txt = "$\sigma_{{noise}}$ = {0} pA".format(int(noise_sigma))
    fig.text(nsX, nsY, noise_sigma_txt, va='center', ha='right', fontsize=19)


gridRoots = getNoiseRoots(gridsDataRoot, noise_sigmas)
gridDataSpace0   = JobTrialSpace2D(shape, gridRoots[0])
gridDataSpace150 = JobTrialSpace2D(shape, gridRoots[1])
gridDataSpace300 = JobTrialSpace2D(shape, gridRoots[2])

exWidth = 6
exHeight = 8

f0_0   = 1
f0_1   = 1
f0_2   = 1
f150_0 = 1
f150_1 = 1
f150_2 = 1
f300_0 = 1
f300_1 = 1
f300_2 = 1


if (f0_0):
    # High gridness Score - 0 pA
    exLeft = 1
    exBottom = 20
    exRect = [exLeft, exBottom, exLeft+exWidth-1, exBottom+exHeight-1]
    drawA4RectExamples(gridDataSpace0, noise_sigmas[0], iterList, exRect,
            exampleIdx[0])

    fname = outputDir + "/suppFigure_grids0.png"
    plt.savefig(fname, dpi=300, transparent=False)
    plt.close()



if (f0_1):
    # Low gridness Score - 0 pA
    exLeft = 20
    exBottom = 10
    exRect = [exLeft, exBottom, exLeft+exWidth-1, exBottom+exHeight-1]
    drawA4RectExamples(gridDataSpace0, noise_sigmas[0], iterList, exRect,
            exampleIdx[0])

    fname = outputDir + "/suppFigure_grids1.png"
    plt.savefig(fname, dpi=300, transparent=False)
    plt.close()

if (f0_1):
    # Transition, low g_E - 0 pA
    exLeft = 3
    exBottom = 6
    exRect = [exLeft, exBottom, exLeft+exWidth-1, exBottom+exHeight-1]
    drawA4RectExamples(gridDataSpace0, noise_sigmas[0], iterList, exRect,
            exampleIdx[0])

    fname = outputDir + "/suppFigure_grids2.png"
    plt.savefig(fname, dpi=300, transparent=False)
    plt.close()


if (f150_0):
    # Transition - 150 pA
    exLeft = 5
    exBottom = 2
    exRect = [exLeft, exBottom, exLeft+exWidth-1, exBottom+exHeight-1]
    drawA4RectExamples(gridDataSpace150, noise_sigmas[1], iterList, exRect,
            exampleIdx[1])

    fname = outputDir + "/suppFigure_grids3.png"
    plt.savefig(fname, dpi=300, transparent=False)
    plt.close()


if (f150_1):
    # Low GS region - 150 pA
    exLeft = 15
    exBottom = 15
    exRect = [exLeft, exBottom, exLeft+exWidth-1, exBottom+exHeight-1]
    drawA4RectExamples(gridDataSpace150, noise_sigmas[1], iterList, exRect,
            exampleIdx[1])

    fname = outputDir + "/suppFigure_grids4.png"
    plt.savefig(fname, dpi=300, transparent=False)
    plt.close()


if (f150_2):
    # High g_E - 150 pA
    exLeft = 1
    exBottom = 15
    exRect = [exLeft, exBottom, exLeft+exWidth-1, exBottom+exHeight-1]
    drawA4RectExamples(gridDataSpace150, noise_sigmas[1], iterList, exRect,
            exampleIdx[1])

    fname = outputDir + "/suppFigure_grids5.png"
    plt.savefig(fname, dpi=300, transparent=False)
    plt.close()


if (f300_0):
    # Transition - 300 pA
    exLeft = 3
    exBottom = 7
    exRect = [exLeft, exBottom, exLeft+exWidth-1, exBottom+exHeight-1]
    drawA4RectExamples(gridDataSpace300, noise_sigmas[2], iterList, exRect,
            exampleIdx[2])

    fname = outputDir + "/suppFigure_grids6.png"
    plt.savefig(fname, dpi=300, transparent=False)
    plt.close()


if (f300_1):
    # Low score region - 300 pA
    exLeft = 15
    exBottom = 15
    exRect = [exLeft, exBottom, exLeft+exWidth-1, exBottom+exHeight-1]
    drawA4RectExamples(gridDataSpace300, noise_sigmas[2], iterList, exRect,
            exampleIdx[2])

    fname = outputDir + "/suppFigure_grids7.png"
    plt.savefig(fname, dpi=300, transparent=False)
    plt.close()

if (f300_2):
    # Transition, high g_E - 300 pA
    exLeft = 1
    exBottom = 20
    exRect = [exLeft, exBottom, exLeft+exWidth-1, exBottom+exHeight-1]
    drawA4RectExamples(gridDataSpace300, noise_sigmas[2], iterList, exRect,
            exampleIdx[2])

    fname = outputDir + "/suppFigure_grids8.png"
    plt.savefig(fname, dpi=300, transparent=False)
    plt.close()



