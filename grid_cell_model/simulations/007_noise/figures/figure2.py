#!/usr/bin/env python
#
#   figure2.py
#
#   Noise publication Figure 2.
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
import matplotlib.pyplot as plt
from matplotlib.pyplot   import figure, subplot, plot, savefig
from matplotlib.gridspec import GridSpec

from parameters  import JobTrialSpace2D
from EI_plotting import plotGridTrial, computeYX
from plotting.grids import plotGridRateMap, plotAutoCorrelation, plotSpikes2D

import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)

from matplotlib import rc
rc('pdf', fonttype=42)
rc('mathtext', default='regular')

plt.rcParams['font.size'] = 12

outputDir = "."
figSize = (12, 8)

NTrials=10
iterList  = ['g_AMPA_total', 'g_GABA_total']

trialNum = 0
jobNum = 573
dataRootDir = 'output_local/grids'
root0   = "{0}/EI_param_sweep_0pA".format(dataRootDir)
root150 = "{0}/EI_param_sweep_150pA".format(dataRootDir)
root300 = "{0}/EI_param_sweep_300pA".format(dataRootDir)
shape = (30, 30)

##############################################################################



def drawGridSweeps(gs, dataSpace, iterList, NTrials=1, r=0, c=0, yLabelOn=True,
        yticks=True, exRows=[], exCols=[], exLetters=[], exDir=[],
        exColor='white'):
    if (yLabelOn):
        yLabelText = '$w_E$ (nS)'
    else:
        yLabelText = ''

    ax0 = subplot(gs[0, 0:2])
    G = plotGridTrial(dataSpace, ['gridnessScore'], iterList,
            trialNumList=range(NTrials),
            r=r,
            c=c,
            xlabel="$w_I$ (nS)",
            ylabel=yLabelText,
            colorBar=False,
            clBarLabel = "Gridness score",
            clbarNTicks=3,
            yticks=yticks)

    Y, X = computeYX(dataSpace, iterList, r=r, c=c)
    for idx in range(len(exRows)):
        row = exRows[idx]
        col = exCols[idx]
        x = X[row, col]
        y = Y[row, col]
        arrowDir = exDir[idx]
        #ax0.plot(x, y, 'o', color=clr)
        ax0.annotate(exLetters[idx],
            xy=(x, y), xycoords='data',
            xytext=(x+1.0*arrowDir[0], y+0.75*arrowDir[1]), textcoords='data',
            color=exColor, fontweight='bold',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3",
                            color=exColor, linewidth=2),
            )



def drawGridExample(gs, dataSpace, dsRows, dsCols, trialNum=0, colStart=0,
        rowStart=0, scaleBarFlag=False):

    for idx in range(len(dsRows)):
        if (idx == len(dsRows) - 1 and scaleBarFlag==True):
            scaleBar = 50
        else:
            scaleBar = None
        r = dsRows[idx]
        c = dsCols[idx]
        d = dataSpace[r][c][trialNum].data
        a = d['analysis']

        arenaDiam = d['options']['arenaSize']

        ax0 = subplot(gs[rowStart, colStart+idx]) 
        plotSpikes2D(a['spikes_e'], a['rat_pos_x'], a['rat_pos_y'],
                a['rat_dt'], diam=arenaDiam, spikeDotSize=3, scaleBar=scaleBar,
                scaleText=None)

        ax1 = subplot(gs[rowStart+1, colStart+idx]) 
        rateMap = a['rateMap_e']
        X       = a['rateMap_e_X']
        Y       = a['rateMap_e_Y']
        plotGridRateMap(rateMap, X, Y, diam=arenaDiam, scaleBar=scaleBar,
                scaleText=False)

        ax2 = subplot(gs[rowStart+2, colStart+idx]) 
        X = a['corr_X']
        Y = a['corr_Y']
        ac = a['corr']
        plotAutoCorrelation(ac, X, Y, diam=arenaDiam, scaleBar=scaleBar)


fig = figure(figsize=figSize)

th = 1 # top plot height
hr = 0.75
space = 0.5
top = 0.93
bottom = 0.05
margin = 0.05
div = 0.08
width = 0.25
hspace = 0
wspace = 0.4
gsHeightRatios = [th, space, hr, hr, hr]

letter_top=0.95
letter_left_offset=0.04
letter_va='bottom'
letter_ha='left'
letter_ex_yoffset=0.4
letter_ex_xoffset=0.01
letter_ex2_mult = 0.55

gsRows = 5
gsCols = 2

# noise_sigm = 0 pA
left = margin
right = left + width
gs = GridSpec(gsRows, gsCols, height_ratios=gsHeightRatios, hspace=hspace)
gs.update(left=left, right=right, bottom=bottom, top=top, wspace=wspace)
dataSpace = JobTrialSpace2D(shape, root0)
exRows = [28, 15]
exCols = [3, 15]
exDir  = [(1., -1), (1., 1.)]
drawGridSweeps(gs, dataSpace, iterList, NTrials=NTrials, r=1, c=2,
        exRows=exRows, exCols=exCols, exLetters=['B', 'C'], exDir=exDir)
drawGridExample(gs, dataSpace, dsRows=exRows, dsCols=exCols, trialNum=0,
        colStart=0, rowStart=2)
fig.text(left-letter_left_offset, letter_top, "A", va=letter_va, ha=letter_ha,
        fontsize=19, fontweight='bold')
fig.text(left-letter_ex_xoffset, letter_top-letter_ex_yoffset, "B", va=letter_va,
        ha=letter_ha, fontsize=19, fontweight='bold')
fig.text(left+letter_ex2_mult*width, letter_top-letter_ex_yoffset, "C",
        va=letter_va, ha=letter_ha, fontsize=19, fontweight='bold')


# noise_sigma = 150 pA
gs = GridSpec(gsRows, gsCols, height_ratios=gsHeightRatios, hspace=hspace)
left = right + div
right = left + width
gs.update(left=left, right=right, bottom=bottom, top=top, wspace=wspace)
dataSpace = JobTrialSpace2D(shape, root150)
exRows = [8, 2]
exCols = [10, 9]
exDir  = [(-1., 1), (1., 1.)]
drawGridSweeps(gs, dataSpace, iterList, NTrials=NTrials, r=11, c=10,
        yLabelOn=False, yticks=False, exRows=exRows, exCols=exCols,
        exLetters=['E', 'F'], exDir=exDir)
drawGridExample(gs, dataSpace, dsRows=exRows, dsCols=exCols, trialNum=0,
        colStart=0, rowStart=2)
fig.text(left-letter_left_offset, letter_top, "D", va=letter_va, ha=letter_ha,
        fontsize=19, fontweight='bold')
fig.text(left-letter_ex_xoffset, letter_top-letter_ex_yoffset, "E", va=letter_va,
        ha=letter_ha, fontsize=19, fontweight='bold')
fig.text(left+letter_ex2_mult*width, letter_top-letter_ex_yoffset, "F",
        va=letter_va, ha=letter_ha, fontsize=19, fontweight='bold')



# noise_sigma = 300 pA
gs = GridSpec(gsRows, gsCols, height_ratios=gsHeightRatios, hspace=hspace)
left = right + div
right = left + width
gs.update(left=left, right=right, bottom=bottom, top=top, wspace=wspace)
dataSpace = JobTrialSpace2D(shape, root300)
exRows = [16, 15]
exCols = [6, 23]
exDir  = [(1, 1), (1., 1.)]
drawGridSweeps(gs, dataSpace, iterList, NTrials=NTrials, r=0, c=5,
        yLabelOn=False, yticks=False, exRows=exRows, exCols=exCols,
        exLetters=['H', 'I'], exDir=exDir, exColor='black')
drawGridExample(gs, dataSpace, dsRows=exRows, dsCols=exCols, trialNum=0,
        colStart=0, rowStart=2, scaleBarFlag=True)
fig.text(left-letter_left_offset, letter_top, "G", va=letter_va, ha=letter_ha,
        fontsize=19, fontweight='bold')
fig.text(left-letter_ex_xoffset, letter_top-letter_ex_yoffset, "H", va=letter_va,
        ha=letter_ha, fontsize=19, fontweight='bold')
fig.text(left+letter_ex2_mult*width, letter_top-letter_ex_yoffset, "I",
        va=letter_va, ha=letter_ha, fontsize=19, fontweight='bold')

fname = outputDir + "/figure2.png"
savefig(fname, dpi=300)

