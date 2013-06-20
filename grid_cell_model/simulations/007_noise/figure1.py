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
from matplotlib.ticker import MaxNLocator
import svgutils.transform as sg

from plotting.global_defs import globalAxesSettings
from figures.fig_conn_func import plotWeights
from analysis.visitors.interface import extractStateVariable, sumAllVariables
from data_storage import DataStorage

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)
rc('pdf', fonttype=42)

outputDir = "."
figSize = (13, 8)

theta_T = 250.0 # ms
theta_dt = 1.0 # ms
theta_f = 8.0  # Hz
theta_DC = 300.0 # pA
theta_A  = 375.0 # pA
theta_max = 1500
theta_min = -500

trialNum = 0
jobNum = 573
dataRootDir = 'output_local'
root0   = "{0}/2013-04-24T15-27-30_EI_param_sweep_0pA_big".format(dataRootDir)
root150 = "{0}/2013-04-24T21-37-47_EI_param_sweep_150pA_big".format(dataRootDir)
root300 = "{0}/2013-04-24T21-43-32_EI_param_sweep_300pA_big".format(dataRootDir)
fileTemplate = "job{0:05}_output.h5"

##############################################################################

def openJob(rootDir, jobNum, trialNum):
    fileTemplate = "job{0:05}_output.h5"
    fileName = rootDir + '/' + fileTemplate.format(jobNum)
    return DataStorage.open(fileName, 'r')['trials'][trialNum]


def getThetaSignal(noise_sigma):
    t = np.arange(0, theta_T, theta_dt) * 1e-3
    phase = np.pi
    normSig = 0.5*(1. + np.cos(2*np.pi*theta_f*t + phase))
    noise = noise_sigma * np.random.randn(len(t))
    return t, theta_DC + theta_A*normSig + noise


def setSignalAxes(ax, leftSpineOn):
    globalAxesSettings(ax)
    ax.xaxis.set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if (leftSpineOn == False):
        ax.spines['left'].set_visible(False)
        ax.yaxis.set_visible(False)

    ax.yaxis.set_major_locator(MaxNLocator(2))

def plotStateSignal(ax, t, sig, leftSpineOn=True, labely="", color='black'):
    setSignalAxes(ax, leftSpineOn)

    if (sig is not None):
        ax.plot(t, sig, color=color)

    #ax.set_ylabel(labely)
    ax.text(-0.22, 0.5, labely,
        verticalalignment='center', horizontalalignment='right',
        transform=ax.transAxes,
        rotation=90)


def plotThetaSignal(ax, noise_sigma, yLabelOn):
    setSignalAxes(ax, leftSpineOn=False)
    t, theta = getThetaSignal(noise_sigma)
    ax.plot(t, theta, color="grey")
    ax.set_ylim([theta_min, theta_max])
    txt = '$\sigma = ' + str(noise_sigma) + '\ \mathrm{pA}$'
    ax.text(0.5, 1.1, txt,
            verticalalignment='bottom', horizontalalignment='center',
            transform=ax.transAxes,
            fontsize=17, fontweight='bold')
    ax.axhline(0.0, color='grey', linestyle=':', linewidth=0.5)
    if (yLabelOn):
        ax.text(theta_T*1e-3 - 0.01, -50, "0 pA", ha="right", va='top', fontsize='small')
    

def plotBump(ax, rateMap):
    rateMap = np.zeros((10, 10))
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.pcolormesh(rateMap)
    axis("scaled")

def plotSpikes(ax, t, trajectory, spikeTimes):
    pass



def sliceSignal(t, sig, tStart, tEnd):
    idx = np.logical_and(t >= tStart, t <= tEnd)
    return t[idx], sig[idx], idx

def extractStateVars(mon, varName, plotTStart, plotTEnd):
    '''
    Extract state variables from a pair of monitors. One in the centre, the
    other one at the edge of the neural sheet.
    '''
    nIdxMiddle = 0
    nIdxEdge = 1

    t, dt = extractStateVariable(mon, nIdxMiddle, 'times')
    sig, dt = sumAllVariables(mon, nIdxMiddle, varName)
    t, sigMiddle, idx = sliceSignal(t, sig, plotTStart, plotTEnd)
    sig, dt = sumAllVariables(mon, nIdxEdge, varName)
    sigEdge = sig[idx]
    return t, sigMiddle, sigEdge


def drawSignals(gs, data, colStart, noise_sigma, yLabelOn=True, letter='',
        letterPos=None):
    if (yLabelOn):
        VmText = "V (mV)"
        IsynText = "I (pA)"
    else:
        VmText = ""
        IsynText = ""

    ncols = 4
    plotTStart = 5e3
    plotTEnd   = 5.25e3

    ax0 = subplot(gs[0, colStart:colStart+ncols])
    plotThetaSignal(ax0, noise_sigma, yLabelOn)

    mon_e = data['stateMon_e']
    mon_i = data['stateMon_i']

    # E cell Vm
    ax1 = subplot(gs[1, colStart:colStart+ncols])
    t, VmMiddle, VmEdge = extractStateVars(mon_e, ['V_m'], plotTStart,
            plotTEnd)
    plotStateSignal(ax1, t, VmMiddle, labely=VmText, color='red')

    # E cell Isyn
    ax2 = subplot(gs[2, colStart:colStart+ncols])
    t, IsynMiddle, IsynEdge = extractStateVars(mon_e, ['I_clamp_GABA_A'],
            plotTStart, plotTEnd)
    plotStateSignal(ax2, t, IsynMiddle, labely=IsynText, color='blue')

    # I cell Vm
    ax3 = subplot(gs[3, colStart:colStart+ncols])
    t, VmMiddle, VmEdge = extractStateVars(mon_i, ['V_m'], plotTStart,
            plotTEnd)
    plotStateSignal(ax3, t, VmMiddle, labely=VmText, color='red')

    # I cell Isyn
    ax4 = subplot(gs[4, colStart:colStart+ncols])
    t, IsynMiddle, IsynEdge = extractStateVars(mon_i, ['I_clamp_AMPA',
        'I_clamp_NMDA'], plotTStart, plotTEnd)
    plotStateSignal(ax4, t, IsynMiddle, labely=IsynText, color='blue')

    #ax5 = subplot(gs[5, colStart])
    #plotBump(ax5, None)

    #ax6 = subplot(gs[5, colStart+1])
    #plotBump(ax6, None)

    #ax7 = subplot(gs[5, colStart+2])
    #plotBump(ax7, None)

    #ax8 = subplot(gs[5, colStart+3])
    #plotBump(ax8, None)

    if (yLabelOn):
        ax1.text(-0.32, -0.15, "Excitatory layer",
                verticalalignment='center', horizontalalignment='right',
                transform=ax1.transAxes,
                rotation=90,
                fontsize=16)
        ax3.text(-0.32, -0.15, "Inhibitory layer",
                verticalalignment='center', horizontalalignment='right',
                transform=ax3.transAxes,
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

letter_top=0.97
letter_div = 0.05
letter_left=0.01
letter_va='top'
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
ds = openJob(root0, jobNum, trialNum)
gs = GridSpec(5, 4, height_ratios=[th, hr, hr, hr, hr], hspace=hspace)
# do not update left and right
gs.update(left=left, right=right, bottom=bottom, top=top)
drawSignals(gs, ds, colStart=0, noise_sigma=0, letter="C")
fig.text(letter_left, top+letter_div, "C", va=letter_va, ha=letter_ha,
        fontsize=19, fontweight='bold')


# noise_sigma = 150 pA
ds = openJob(root150, jobNum, trialNum)
gs = GridSpec(5, 4, height_ratios=[th, hr, hr, hr, hr], hspace=hspace)
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
ds = openJob(root300, jobNum, trialNum)
gs = GridSpec(5, 4, height_ratios=[th, hr, hr, hr, hr], hspace=hspace)
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
