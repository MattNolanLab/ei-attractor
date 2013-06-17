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
    



def drawSignals(gs, data, colStart, yLabelOn=True, noise_sigma=0):
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


fig = figure(figsize=figSize)

hr = 3
gs = GridSpec(6, 12, height_ratios=[1, hr, hr, hr, hr, hr])
drawSignals(gs, None, colStart=0)
gs.tight_layout(fig, h_pad=1.0)

drawSignals(gs, None, colStart=4, yLabelOn=False)
gs.tight_layout(fig, h_pad=1.0)

drawSignals(gs, None, colStart=8, yLabelOn=False)
gs.tight_layout(fig, h_pad=1.0)

savefig(outputDir + "/figure1.pdf")

