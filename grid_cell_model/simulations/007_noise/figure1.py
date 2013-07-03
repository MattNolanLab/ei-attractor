#!/usr/bin/env python
#
#   figure1.py
#
#   Noise publication Figure 1.
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
from matplotlib.pyplot import *
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import AutoMinorLocator, ScalarFormatter

from figures.fig_conn_func import plotWeights
from data_storage import DataStorage
from figures_shared import plotStateSignal, plotThetaSignal, extractStateVars,\
        getOption, thetaLim

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)
rc('pdf', fonttype=42)

outputDir = "."
figSize = (13, 8)

trialNum = 0
jobNum = 573
dataRootDir = 'output_local'
root0   = "{0}/single_neuron".format(dataRootDir)
root150 = "{0}/single_neuron".format(dataRootDir)
root300 = "{0}/single_neuron".format(dataRootDir)
fileTemplate = "noise_sigma{0}_output.h5"

##############################################################################

def openJob(rootDir, noise_sigma):
    fileName = rootDir + '/' + fileTemplate.format(noise_sigma)
    return DataStorage.open(fileName, 'r')


def plotHistogram(ax, sig, color='black', labelx="", labely="",
        labelyPos=-0.5, powerLimits=(0, 3)):
    hist(sig, bins=100, normed=True, histtype='step', align='mid', color=color)

    # y label manually
    if (labely is None):
        labely = 'Count'
    ax.text(labelyPos, 0.5, labely,
        verticalalignment='center', horizontalalignment='right',
        transform=ax.transAxes,
        rotation=90)
    xlabel(labelx)
    
    ax.minorticks_on()
    ax.xaxis.set_major_locator(LinearLocator(2))
    ax.yaxis.set_major_locator(LinearLocator(2))
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax.yaxis.set_minor_locator(AutoMinorLocator(4))
    f = ScalarFormatter(useMathText=True)
    f.set_scientific(True)
    f.set_powerlimits(powerLimits)
    ax.yaxis.set_major_formatter(f)
    ax.tick_params(
            which='both',
            direction='out'
    )
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')



def drawSignals(gs, data, colStart, noise_sigma, yLabelOn=True, letter='',
        letterPos=None):
    if (yLabelOn):
        VmText = "V (mV)"
        IsynText = "I (nA)"
        countText = None
    else:
        VmText = ""
        IsynText = ""
        countText = ""
    histLabelX = "V (mV)"

    ncols = 4
    plotTStart = 5e3
    plotTEnd   = 5.25e3

    theta_start_t = getOption(data, 'theta_start_t')
    #theta_start_t = 1e3
    simTime = getOption(data, 'time')

    mon_e = data['stateMon_e']
    mon_i = data['stateMon_i']

    ax0 = subplot(gs[0, colStart:colStart+ncols])
    t, IStim = extractStateVars(mon_e, ['I_stim'], plotTStart,
            plotTEnd)
    plotThetaSignal(ax0, t, IStim, noise_sigma, yLabelOn, thetaLim)

    # E cell Vm
    ax1 = subplot(gs[1, colStart:colStart+ncols])
    t, VmMiddle = extractStateVars(mon_e, ['V_m'], plotTStart,
            plotTEnd)
    plotStateSignal(ax1, t, VmMiddle, labely=VmText, color='red')

    # I cell Vm
    ax2 = subplot(gs[2, colStart:colStart+ncols])
    t, VmMiddle = extractStateVars(mon_i, ['V_m'], plotTStart,
            plotTEnd)
    plotStateSignal(ax2, t, VmMiddle, labely=VmText, color='blue')

    # E cell Vm histogram
    ax3 = subplot(gs[3, colStart:colStart+2])
    t, VmMiddle = extractStateVars(mon_e, ['V_m'], theta_start_t,
            simTime)
    plotHistogram(ax3, VmMiddle, labelx = histLabelX, labely=countText, color='red')

    # I cell Vm histogram
    ax4 = subplot(gs[3, colStart+2:colStart+4])
    t, VmMiddle = extractStateVars(mon_i, ['V_m'], theta_start_t,
            simTime)
    plotHistogram(ax4, VmMiddle, labelx = histLabelX, labely="", color='blue')


    if (yLabelOn):
        ax1.text(-0.32, 0.5, "E cell",
                verticalalignment='center', horizontalalignment='right',
                transform=ax1.transAxes,
                rotation=90,
                fontsize=16)
        ax2.text(-0.32, 0.5, "I cell",
                verticalalignment='center', horizontalalignment='right',
                transform=ax2.transAxes,
                rotation=90,
                fontsize=16)



fig = figure(figsize=figSize)

hr = 1
th = 0.75 # top plot height
top = 0.6
bottom = 0.02
margin = 0.1
div = 0.085
width = 0.23
hspace = 0.3
wspace = 1.2

letter_top=0.97
letter_div = 0.05
letter_left=0.01
letter_va='bottom'
letter_ha='left'

# Model schematic
gs = GridSpec(1, 4)
top_margin = 0.15
top_top = 0.92
top_letter_pos = 1.5
fig.text(letter_left, letter_top, "A", va=letter_va, ha=letter_ha, fontsize=19,
        fontweight='bold')


# noise_sigm = 0 pA
left = margin
right = left + width
ds = openJob(root0, noise_sigma=0)
gs = GridSpec(4, 4, height_ratios=[th, hr, hr, hr, hr], hspace=hspace,
        wspace=wspace)
# do not update left and right
gs.update(left=left, right=right, bottom=bottom, top=top)
drawSignals(gs, ds, colStart=0, noise_sigma=0, letter="C")
fig.text(letter_left, top+letter_div, "C", va=letter_va, ha=letter_ha,
        fontsize=19, fontweight='bold')


# noise_sigma = 150 pA
ds = openJob(root150, noise_sigma=150)
gs = GridSpec(4, 4, height_ratios=[th, hr, hr, hr, hr], hspace=hspace,
        wspace=wspace)
left = right + div
right = left + width
gs.update(left=left, right=right, bottom=bottom, top=top)
drawSignals(gs, ds, colStart=0, yLabelOn=False, noise_sigma=150, letter="D",
        letterPos=-0.2)
fig.text(letter_left+margin+width+0.5*div, top+letter_div, "D", va=letter_va,
        ha=letter_ha, fontsize=19, fontweight='bold')


# Connection weights
gs = GridSpec(1, 4)
w_shift = 0.4*width + div
w_left  = left - w_shift
w_right = w_left + 0.15
gs.update(left=w_left, right=w_right, bottom=top+top_margin, top=top_top)
ax_sch = subplot(gs[0, :])
plotWeights(ax_sch)
fig.text(w_left-0.5*div, letter_top, "B", va=letter_va, ha=letter_ha, fontsize=19,
        fontweight='bold')
#gs.tight_layout(fig, rect=[w_left-0.5*div, top+top_margin, w_right, top_top])


# noise_sigma = 300 pA
ds = openJob(root300, noise_sigma=300)
gs = GridSpec(4, 4, height_ratios=[th, hr, hr, hr, hr], hspace=hspace,
        wspace=wspace)
left = right + div
right = left + width
gs.update(left=left, right=right, bottom=bottom, top=top)
drawSignals(gs, ds, colStart=0, yLabelOn=False, noise_sigma=300, letter="E",
        letterPos=-0.2)
fig.text(letter_left+margin+2*width+1.5*div, top+letter_div, "E", va=letter_va,
        ha=letter_ha, fontsize=19, fontweight='bold')

fname = outputDir + "/figure1.pdf"
savefig(fname)

## Paste-in model schematic
#fig_final = sg.fromfile(fname)
#fig_sch = sg.fromfile('figures/model_schematic.svg')
#
#plot = fig_sch.getroot()
#fig_final.append(plot)
#fig_final.save(fname)
