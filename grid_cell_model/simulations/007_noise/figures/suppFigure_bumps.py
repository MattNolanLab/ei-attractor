#!/usr/bin/env python
#
#   suppFigure_bumps.py
#
#   Supplementary figure: bump examples
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
from EI_plotting import plotBumpSigmaTrial, computeYX, aggregate2D, \
        drawBumpExamples,  drawEIRectSelection
from plotting.global_defs import globalAxesSettings, createColorbar
from figures_shared import getNoiseRoots

import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)

from matplotlib import rc
rc('pdf', fonttype=42)
rc('mathtext', default='regular')

plt.rcParams['font.size'] = 12

outputDir = "."

NTrials=10
exampleIdx = [(0, 0), (0, 0), (0, 0)] # (row, col)
iterList  = ['g_AMPA_total', 'g_GABA_total']

noise_sigmas = [0, 150, 300]
bumpsDataRoot= 'output_local/even_spacing/gamma_bump'
shape = (31, 31)


##############################################################################


def drawSweep(ax, dataSpace, iterList, spaceRect, exIdx=(0, 0), cmap='jet_r',
        rectColor='black'):
    exRow, exCol = exIdx
    Y, X = computeYX(dataSpace, iterList, r=exRow, c=exCol)

    # Replot gridness score param. sweep
    ax0 = plt.gca()
    varList = ['bump_e', 'sigma']
    G = plotBumpSigmaTrial(dataSpace, varList, iterList,
            trialNumList=range(NTrials),
            colorBar=False,
            clbarNTicks=3,
            vmin=0,
            vmax=10,
            cmap=cmap)
    plt.set_cmap('jet_r')

    cax, kw = make_axes(ax0, orientation='vertical',
            nticks=4, shrink=0.9)
    globalAxesSettings(cax)
    cb = createColorbar(ax, None, "Bump $\sigma$ (neurons)", cax=cax, **kw)

    drawEIRectSelection(ax0, spaceRect, X, Y, color=rectColor)



def drawA4RectExamples(dataSpace, noise_sigma, iterList, exRect, exIdx,
        rectColor='black'):
    fig = figure(figsize=(8.27, 11.69))
    margin    = 0.1
    sw_left   = margin
    sw_bottom = 0.82
    sw_right  = 0.4
    sw_top    = 0.97
    div       = 0.075

    letter_left = 0.03
    letter_top_off = 0.01
    letter_va='bottom'
    letter_ha='left'

    nsX = 0.9
    nsY = sw_top

    sweepsRect = sw_left, sw_bottom, sw_right-sw_left, sw_top-sw_bottom
    ax_sweeps = fig.add_axes(sweepsRect)
    drawSweep(ax_sweeps, dataSpace, iterList, exRect, exIdx,
            rectColor=rectColor)
    fig.text(letter_left, sw_top+letter_top_off, "A", va=letter_va, ha=letter_ha,
            fontsize=19, fontweight='bold')

    gsCoords = 0.12, 0.075, 0.95, sw_bottom-div
    #gsCoords = margin, 0.46, 0.5, sw_bottom-div
    gs = drawBumpExamples(dataSpace, exRect, iterList, gsCoords=gsCoords,
            exIdx=exIdx, cmap='jet')
    fig.text(letter_left, sw_bottom-div+letter_top_off, "B", va=letter_va,
            ha=letter_ha, fontsize=19, fontweight='bold')
    noise_sigma_txt = "$\sigma_{{noise}}$ = {0} pA".format(int(noise_sigma))
    fig.text(nsX, nsY, noise_sigma_txt, va='center', ha='right', fontsize=19)


bumpRoots = getNoiseRoots(bumpsDataRoot, noise_sigmas)
bumpDataSpace0   = JobTrialSpace2D(shape, bumpRoots[0])
bumpDataSpace150 = JobTrialSpace2D(shape, bumpRoots[1])
bumpDataSpace300 = JobTrialSpace2D(shape, bumpRoots[2])

exWidth = 6
exHeight = 8

# Flags
f0_0   = 0
f0_1   = 0
f0_2   = 0
f150_0 = 0
f150_1 = 0
f150_2 = 0
f300_0 = 0
f300_1 = 0
f300_2 = 1


if (f0_0):
    # Bumps - 0 pA
    exLeft = 3
    exBottom = 19
    exRect = [exLeft, exBottom, exLeft+exWidth-1, exBottom+exHeight-1]
    drawA4RectExamples(bumpDataSpace0, noise_sigmas[0], iterList, exRect,
            exampleIdx[0])

    fname = outputDir + "/suppFigure_bumps0.png"
    savefig(fname, dpi=300, transparent=False)


if (f0_1):
    # Transition - 0 pA
    exLeft = 0
    exBottom = 3
    exRect = [exLeft, exBottom, exLeft+exWidth-1, exBottom+exHeight-1]
    drawA4RectExamples(bumpDataSpace0, noise_sigmas[0], iterList, exRect,
            exampleIdx[0], rectColor='red')

    fname = outputDir + "/suppFigure_bumps1.png"
    savefig(fname, dpi=300, transparent=False)


if (f0_2):
    # No bump region at all - 0 pA
    exLeft = 18
    exBottom = 3
    exRect = [exLeft, exBottom, exLeft+exWidth-1, exBottom+exHeight-1]
    drawA4RectExamples(bumpDataSpace0, noise_sigmas[0], iterList, exRect,
            exampleIdx[0], rectColor='red')

    fname = outputDir + "/suppFigure_bumps2.png"
    savefig(fname, dpi=300, transparent=False)


if (f150_0):
    # Transition - 150 pA
    exLeft = 5
    exBottom = 2
    exRect = [exLeft, exBottom, exLeft+exWidth-1, exBottom+exHeight-1]
    drawA4RectExamples(bumpDataSpace150, noise_sigmas[1], iterList, exRect,
            exampleIdx[1], rectColor='red')

    fname = outputDir + "/suppFigure_bumps3.png"
    savefig(fname, dpi=300, transparent=False)


if (f150_1):
    # Good bumps - 150 pA
    exLeft = 15
    exBottom = 10
    exRect = [exLeft, exBottom, exLeft+exWidth-1, exBottom+exHeight-1]
    drawA4RectExamples(bumpDataSpace150, noise_sigmas[1], iterList, exRect,
            exampleIdx[1])

    fname = outputDir + "/suppFigure_bumps4.png"
    savefig(fname, dpi=300, transparent=False)


if (f150_2):
    # Another transition - 150 pA
    exLeft = 1
    exBottom = 20
    exRect = [exLeft, exBottom, exLeft+exWidth-1, exBottom+exHeight-1]
    drawA4RectExamples(bumpDataSpace150, noise_sigmas[1], iterList, exRect,
            exampleIdx[1], rectColor='red')

    fname = outputDir + "/suppFigure_bumps5.png"
    savefig(fname, dpi=300, transparent=False)


if (f300_0):
    # Transition - 150 pA
    exLeft = 5
    exBottom = 3
    exRect = [exLeft, exBottom, exLeft+exWidth-1, exBottom+exHeight-1]
    drawA4RectExamples(bumpDataSpace300, noise_sigmas[2], iterList, exRect,
            exampleIdx[2], rectColor='red')

    fname = outputDir + "/suppFigure_bumps6.png"
    savefig(fname, dpi=300, transparent=False)


if (f300_1):
    # Small bumps - 150 pA
    exLeft = 15
    exBottom = 10
    exRect = [exLeft, exBottom, exLeft+exWidth-1, exBottom+exHeight-1]
    drawA4RectExamples(bumpDataSpace300, noise_sigmas[2], iterList, exRect,
            exampleIdx[2])

    fname = outputDir + "/suppFigure_bumps7.png"
    savefig(fname, dpi=300, transparent=False)


if (f300_2):
    # Small bumps - 150 pA
    exLeft = 1
    exBottom = 11
    exRect = [exLeft, exBottom, exLeft+exWidth-1, exBottom+exHeight-1]
    drawA4RectExamples(bumpDataSpace300, noise_sigmas[2], iterList, exRect,
            exampleIdx[2])

    fname = outputDir + "/suppFigure_bumps8.png"
    savefig(fname, dpi=300, transparent=False)


