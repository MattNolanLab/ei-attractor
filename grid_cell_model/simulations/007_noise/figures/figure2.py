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
from matplotlib.pyplot   import figure, subplot, savefig
from matplotlib.gridspec import GridSpec

from parameters  import JobTrialSpace2D
from EI_plotting import plotGridTrial
from plotting.grids import plotGridRateMap, plotAutoCorrelation, plotSpikes2D

import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)
rc('pdf', fonttype=42)

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
        yticks=True):
    if (yLabelOn):
        yLabelText = 'E (nS)'
    else:
        yLabelText = ''

    ax0 = subplot(gs[0, 0:2])
    G = plotGridTrial(dataSpace, ['gridnessScore'], iterList,
            trialNumList=range(NTrials),
            r=r,
            c=c,
            xlabel="I (nS)",
            ylabel=yLabelText,
            colorBar=False,
            clBarLabel = "Gridness score",
            clbarNTicks=3,
            yticks=yticks)


    #if (yLabelOn):
    #    labelPos = -0.32
    #    ax1.text(labelPos, -0.15, "Excitatory cell",
    #            verticalalignment='center', horizontalalignment='right',
    #            transform=ax1.transAxes,
    #            rotation=90,
    #            fontsize='large')
    #    ax3.text(labelPos, -0.15, "Inhibitory cell",
    #            verticalalignment='center', horizontalalignment='right',
    #            transform=ax3.transAxes,
    #            rotation=90,
    #            fontsize='large')


def drawGridExample(gs, dataSpace, dsRows, dsCols, trialNum=0, colStart=0,
        rowStart=0):

    for idx in range(len(dsRows)):
        if (idx == len(dsRows) - 1):
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
                a['rat_dt'], spikeDotSize=3, scaleBar=scaleBar, scaleText=None)

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
wspace = 0.1
gsHeightRatios = [th, space, hr, hr, hr]

letter_top=0.97
letter_div = 0.02
letter_left=0.01
letter_va='bottom'
letter_ha='left'

gsRows = 5
gsCols = 2

# noise_sigm = 0 pA
left = margin
right = left + width
gs = GridSpec(gsRows, gsCols, height_ratios=gsHeightRatios, hspace=hspace)
gs.update(left=left, right=right, bottom=bottom, top=top)
dataSpace = JobTrialSpace2D(shape, root0)
drawGridSweeps(gs, dataSpace, iterList, NTrials=NTrials, r=1, c=2)
drawGridExample(gs, dataSpace, dsRows=[28, 15], dsCols=[3, 15], trialNum=0,
        colStart=0, rowStart=2)


# noise_sigma = 150 pA
gs = GridSpec(gsRows, gsCols, height_ratios=gsHeightRatios, hspace=hspace)
left = right + div
right = left + width
gs.update(left=left, right=right, bottom=bottom, top=top)
dataSpace = JobTrialSpace2D(shape, root150)
drawGridSweeps(gs, dataSpace, iterList, NTrials=NTrials, r=11, c=10,
        yLabelOn=False, yticks=False)
drawGridExample(gs, dataSpace, dsRows=[10, 2], dsCols=[10, 9], trialNum=0,
        colStart=0, rowStart=2)



# noise_sigma = 300 pA
gs = GridSpec(gsRows, gsCols, height_ratios=gsHeightRatios, hspace=hspace)
left = right + div
right = left + width
gs.update(left=left, right=right, bottom=bottom, top=top)
dataSpace = JobTrialSpace2D(shape, root300)
drawGridSweeps(gs, dataSpace, iterList, NTrials=NTrials, r=0, c=5,
        yLabelOn=False, yticks=False)
drawGridExample(gs, dataSpace, dsRows=[16, 15], dsCols=[7, 23], trialNum=0,
        colStart=0, rowStart=2)

fname = outputDir + "/figure2.png"
savefig(fname, dpi=300)

