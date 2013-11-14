#!/usr/bin/env python
#
#   figure_background_current.py
#
#   Figure showing the background input current: constant component and theta
#   signal.
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
from matplotlib import rcParams as rcp

import settings as se
from data_storage import DataStorage
from data_storage.sim_models.ei import extractSummedSignals
from plotting.signal import signalPlot
from plotting.global_defs import globalAxesSettings

import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)

dataRootDir = 'data/model/E_surround'

plotTStart = 5e3
plotTEnd   = 5.25e3
##############################################################################

def openJob(rootDir, jobNum=0, trialNum=0):
    fileName = "{0}/job{1:05}_output.h5".format(rootDir, jobNum)
    return DataStorage.open(fileName, 'r')['trials'][trialNum]


##############################################################################
transparent = True
nThetaTicks = 3
xmargin = 0.02

letter_top  = 0.95
letter_div  = 0.02
letter_left = 0.01
letter_va   = 'bottom'
letter_ha   = 'left'

fig = plt.figure(figsize=(5.5, 2))

# Plot actual theta signals
ax_theta = plt.subplot(1, 2, 2)
globalAxesSettings(ax_theta)
data = openJob(dataRootDir)

mon_e = data['stateMon_e']
mon_i = data['stateMon_i']

t, IStim_e = extractSummedSignals(mon_e, ['I_stim'], plotTStart, plotTEnd)
signalPlot(t, IStim_e, ax_theta, nThetaTicks=nThetaTicks, color='red',
        xmargin=xmargin, zeroLine=False)
t, IStim_i = extractSummedSignals(mon_i, ['I_stim'], plotTStart, plotTEnd)
signalPlot(t, IStim_i, ax_theta, nThetaTicks=nThetaTicks, color='blue',
        xmargin=xmargin)
ax_theta.set_ylim([-300, 1200])
ax_theta.yaxis.set_major_locator(ti.MultipleLocator(400))
ax_theta.yaxis.set_minor_locator(ti.AutoMinorLocator(2))

leg = ['E cell', 'I cell']
l = ax_theta.legend(leg, loc=(0, 1), ncol=2, frameon=False, fontsize='small')

ax_theta.set_title('B', weight='bold', size=se.letterSize, x=-0.3, y=1.1, va='bottom')



# Plot the theta signal schematic
ax_sch = plt.subplot(1, 2, 1)
globalAxesSettings(ax_sch)

t      = t - t[0]
Iconst = 300.0 # pA
A      = 360.0 # pA
phi    = np.pi
f      = 8.0   # Hz
sig    = Iconst + 0.5*A*(1 + np.cos(2.0*np.pi*f*t*1e-3 - phi))
signalPlot(t, sig, ax_sch, nThetaTicks=nThetaTicks,
        ylabel='I (pA)',
        color='blue',
        xmargin=xmargin)

color  = rcp['grid.color']
ls     = rcp['grid.linestyle']
lw     = rcp['grid.linewidth']
alphsa = rcp['grid.alpha']
ax_sch.axhline(Iconst, ls=ls, lw=lw, color=color)
ax_sch.annotate('', xy=(30, 0), xycoords='data',
        xytext=(30, Iconst), textcoords='data',
        arrowprops=dict(arrowstyle='<->', shrinkA=0, shrinkB=0))
ax_sch.text(40, Iconst/2.0, "Constant amplitude", transform=ax_sch.transData,
        ha='left', va='center', size='small')

ax_sch.annotate('', xy=(62.5, Iconst), xycoords='data',
        xytext=(62.5, Iconst+A), textcoords='data',
        arrowprops=dict(arrowstyle='<->', shrinkA=0, shrinkB=0))
ax_sch.annotate('$\\theta$ amplitude', xy=(62.5, Iconst+0.5*A), xycoords='data',
        xytext=(0.4, 1.1), textcoords='axes fraction', ha='left', va='center',
        arrowprops=dict(
            arrowstyle='->',
            linewidth=0.5,
            relpos=(0, 0.5),
            shrinkA=0))
ax_sch.set_title('A', weight='bold', size=se.letterSize, x=-0.36, y=1.1, va='bottom')


# Save
fig.tight_layout(rect=[0.01, 0, 0.99, 1], w_pad=1, pad=0)
fname = se.setFName("background_current.eps")
fig.savefig(fname, transparent=transparent)

