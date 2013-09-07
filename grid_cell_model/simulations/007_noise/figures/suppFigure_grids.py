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
from matplotlib.pyplot   import figure, subplot, plot, savefig
from matplotlib.colorbar import make_axes

from parameters  import JobTrialSpace2D
from EI_plotting import plotGridTrial, computeYX, aggregate2D, \
        drawGridExamples,  drawEIRectSelection
from plotting.global_defs import globalAxesSettings, createColorbar

import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)

from matplotlib import rc
rc('pdf', fonttype=42)
rc('mathtext', default='regular')

plt.rcParams['font.size'] = 12

outputDir = "."

NTrials=10
exampleIdx = [(1, 2), (11, 10), (0, 5)] # (row, col)
iterList  = ['g_AMPA_total', 'g_GABA_total']

noise_sigmas = [0, 150, 300]
gridsDataRoot= 'output_local/grids'
velDataRoot = 'output_local/velocity'
shape = (30, 30)


##############################################################################

def getNoiseRootDir(prefix, noise_sigma):
    return  "{0}/EI_param_sweep_{1}pA".format(prefix, int(noise_sigma))


def getNoiseRoots(prefix, noise_sigmas):
    roots = []
    for s in noise_sigmas:
        roots.append(getNoiseRootDir(prefix, s))
    return roots


def drawSweep(ax, dataSpace, iterList, spaceRect, exIdx=(0, 0)):
    exRow, exCol = exIdx
    Y, X = computeYX(dataSpace, iterList, r=exRow, c=exCol)

    # Replot gridness score param. sweep
    ax0 = plt.gca()
    G = plotGridTrial(dataSpace, ['gridnessScore'], iterList,
            trialNumList=range(NTrials),
            r=exRow,
            c=exCol,
            xlabel="$w_I$ (nS)",
            ylabel="$w_E$ (nS)",
            colorBar=False,
            clBarLabel = "Gridness score",
            clbarNTicks=3,
            vmin=None,
            vmax=None)

    cax, kw = make_axes(ax0, orientation='vertical',
            nticks=4, shrink=0.9)
    globalAxesSettings(cax)
    cb = createColorbar(ax, None, "Gridness score", cax=cax, **kw)

    drawEIRectSelection(ax0, spaceRect, X, Y)



def drawA4RectExamples(dataSpace, noise_sigma, iterList, exRect, exIdx):
    fig = figure(figsize=(8.27, 11.69))
    margin    = 0.1
    sw_left   = margin
    sw_bottom = 0.8
    sw_right  = 0.5
    sw_top    = 0.95
    div       = 0.1

    letter_left = 0.03
    letter_top_off = 0.02
    letter_va='bottom'
    letter_ha='left'

    nsX = 0.9
    nsY = sw_top

    sweepsRect = sw_left, sw_bottom, sw_right-sw_left, sw_top-sw_bottom
    ax_sweeps = fig.add_axes(sweepsRect)
    drawSweep(ax_sweeps, dataSpace, iterList, exRect, exIdx)
    fig.text(letter_left, sw_top+letter_top_off, "A", va=letter_va, ha=letter_ha,
            fontsize=19, fontweight='bold')

    gsCoords = margin, 0.075, 1.0 - margin, sw_bottom-div
    #gsCoords = margin, 0.46, 0.5, sw_bottom-div
    gs = drawGridExamples(dataSpace, exRect, iterList, gsCoords=gsCoords,
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


# High gridness Score - 0 pA
exLeft = 2
exBottom = 21
exRect = [exLeft, exBottom, exLeft+exWidth-1, exBottom+exHeight-1]
drawA4RectExamples(gridDataSpace0, noise_sigmas[0], iterList, exRect,
        exampleIdx[0])

fname = outputDir + "/suppFigure_grids0.png"
savefig(fname, dpi=300, transparent=False)


# Low gridness Score - 0 pA
exLeft = 18
exBottom = 13
exRect = [exLeft, exBottom, exLeft+exWidth-1, exBottom+exHeight-1]
drawA4RectExamples(gridDataSpace0, noise_sigmas[0], iterList, exRect,
        exampleIdx[0])

fname = outputDir + "/suppFigure_grids1.png"
savefig(fname, dpi=300, transparent=False)


# Transition - 150 pA
exLeft = 6
exBottom = 1
exRect = [exLeft, exBottom, exLeft+exWidth-1, exBottom+exHeight-1]
drawA4RectExamples(gridDataSpace150, noise_sigmas[1], iterList, exRect,
        exampleIdx[1])

fname = outputDir + "/suppFigure_grids2.png"
savefig(fname, dpi=300, transparent=False)


# Low GS region - 150 pA
exLeft = 20
exBottom = 16
exRect = [exLeft, exBottom, exLeft+exWidth-1, exBottom+exHeight-1]
drawA4RectExamples(gridDataSpace150, noise_sigmas[1], iterList, exRect,
        exampleIdx[1])

fname = outputDir + "/suppFigure_grids3.png"
savefig(fname, dpi=300, transparent=False)


# Stripe region - 150 pA
exLeft = 1
exBottom = 18
exRect = [exLeft, exBottom, exLeft+exWidth-1, exBottom+exHeight-1]
drawA4RectExamples(gridDataSpace150, noise_sigmas[1], iterList, exRect,
        exampleIdx[1])

fname = outputDir + "/suppFigure_grids4.png"
savefig(fname, dpi=300, transparent=False)


# Stripe region - 300 pA
exLeft = 4
exBottom = 9
exRect = [exLeft, exBottom, exLeft+exWidth-1, exBottom+exHeight-1]
drawA4RectExamples(gridDataSpace300, noise_sigmas[2], iterList, exRect,
        exampleIdx[2])

fname = outputDir + "/suppFigure_grids5.png"
savefig(fname, dpi=300, transparent=False)


# Low score region - 300 pA
exLeft = 16
exBottom = 9
exRect = [exLeft, exBottom, exLeft+exWidth-1, exBottom+exHeight-1]
drawA4RectExamples(gridDataSpace300, noise_sigmas[2], iterList, exRect,
        exampleIdx[2])

fname = outputDir + "/suppFigure_grids6.png"
savefig(fname, dpi=300, transparent=False)



