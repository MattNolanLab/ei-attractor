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

from plotting.global_defs import globalAxesSettings

#from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Bitstream Vera Sans']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)

outputDir = "."
figSize = (12, 8)

theta_T = 250.0 # ms
theta_dt = 1.0 # ms
theta_f = 8.0  # Hz
theta_DC = 300.0 # pA
theta_A  = 300.0 # pA

##############################################################################


def getThetaSignal(noise_sigma):
    t = np.arange(0, theta_T, theta_dt) * 1e-3
    phase = np.pi
    normSig = 0.5*(1. + np.cos(2*np.pi*theta_f*t + phase))
    return t, theta_DC + theta_A*normSig


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

def plotStateSignal(ax, t, sig, leftSpineOn=True, labely=""):
    setSignalAxes(ax, leftSpineOn)

    if (sig is not None):
        ax.plot(t, sig)

    ax.set_ylabel(labely)


def plotThetaSignal(ax, noise_sigma):
    setSignalAxes(ax, leftSpineOn=False)
    t, theta = getThetaSignal(noise_sigma)
    ax.fill_between(t, theta, color="grey")
    ax.set_ylim([0, np.max(theta)])
    txt = '$\sigma = ' + str(noise_sigma) + '\ \mathrm{pA}$'
    ax.text(0.5, 1.5, txt,
            verticalalignment='bottom', horizontalalignment='center',
            transform=ax.transAxes,
            fontsize=17, fontweight='bold')
    

def plotBump(ax, rateMap):
    rateMap = np.zeros((10, 10))
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.pcolormesh(rateMap)
    axis("scaled")

def plotSpikes(ax, t, trajectory, spikeTimes):
    pass


def drawSignals(gs, data, colStart, noise_sigma, yLabelOn=True):
    if (yLabelOn):
        VmText = "V (mV)"
        IsynText = "I (pA)"
    else:
        VmText = ""
        IsynText = ""

    ncols = 4

    ax0 = subplot(gs[0, colStart:colStart+ncols])
    plotThetaSignal(ax0, noise_sigma)

    ax1 = subplot(gs[1, colStart:colStart+ncols])
    plotStateSignal(ax1, None, None, labely=VmText)

    ax2 = subplot(gs[2, colStart:colStart+ncols])
    plotStateSignal(ax2, None, None, labely=IsynText)

    ax3 = subplot(gs[3, colStart:colStart+ncols])
    plotStateSignal(ax3, None, None, labely=VmText)

    ax4 = subplot(gs[4, colStart:colStart+ncols])
    plotStateSignal(ax4, None, None, labely=IsynText)

    ax5 = subplot(gs[5, colStart])
    plotBump(ax5, None)

    ax6 = subplot(gs[5, colStart+1])
    plotBump(ax6, None)

    ax7 = subplot(gs[5, colStart+2])
    plotBump(ax7, None)

    ax8 = subplot(gs[5, colStart+3])
    plotBump(ax8, None)

    if (yLabelOn):
        ax1.text(-0.3, -0.05, "Excitatory layer",
                verticalalignment='center', horizontalalignment='right',
                transform=ax1.transAxes,
                rotation=90,
                fontsize=16)

        ax3.text(-0.3, -0.05, "Inhibitory layer",
                verticalalignment='center', horizontalalignment='right',
                transform=ax3.transAxes,
                rotation=90,
                fontsize=16)



fig = figure(figsize=figSize)

hr = 1
top = 0.92
bottom = 0.05
margin = 0.05
div = 0.06
width = 0.25
gs = GridSpec(6, 4, height_ratios=[0.3, hr, hr, hr, hr, 0.75])
left = div+margin
right = left + width
gs.update(left=left, right=right, bottom=bottom, top=top)
drawSignals(gs, None, colStart=0, noise_sigma=0)

gs = GridSpec(6, 4, height_ratios=[0.3, hr, hr, hr, hr, 0.75])
left = right + div
right = left + width
gs.update(left=left, right=right, bottom=bottom, top=top)
drawSignals(gs, None, colStart=0, yLabelOn=False, noise_sigma=150)

gs = GridSpec(6, 4, height_ratios=[0.3, hr, hr, hr, hr, 0.75])
left = right + div
right = left + width
gs.update(left=left, right=right, bottom=bottom, top=top)
drawSignals(gs, None, colStart=0, yLabelOn=False, noise_sigma=300)
#gs.tight_layout(fig, h_pad=1.0)

savefig(outputDir + "/figure1.png")

