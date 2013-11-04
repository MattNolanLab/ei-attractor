#
#   rasters.py
#
#   Raster plot for the E/I parameter sweeps
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

from plotting.global_defs import globalAxesSettings

##############################################################################
# Raster plots
def plotEIRaster(ESpikes, ISpikes, tLimits, **kw):
    # kw arguments 
    ax               = kw.pop('ax', plt.gca())
    ylabel           = kw.pop('ylabel', 'Neuron #')
    yticks           = kw.pop('yticks', True)
    yticks_style     = kw.pop('yticks_style', 'separate')
    ylabelPos        = kw.pop('ylabelPos', -0.22)
    EColor           = kw.pop('ecolor', 'red')
    IColor           = kw.pop('icolor', 'blue')
    kw['markersize'] = kw.get('markersize', 1.0)

    ESpikes = ESpikes.windowed(tLimits)
    ISpikes = ISpikes.windowed(tLimits)

    ESenders, ETimes = ESpikes.rasterData()
    ISenders, ITimes = ISpikes.rasterData()
    ISenders += ESpikes.N

    globalAxesSettings(ax)
    ax.minorticks_on()
    ax.xaxis.set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_major_locator(ti.LinearLocator(2))
    ax.yaxis.set_minor_locator(ti.NullLocator())

    ax.plot(ETimes, ESenders+1, '.', color='red',  mec='none', **kw)
    ax.plot(ITimes, ISenders+1, '.', color='blue', mec='none', **kw)

    ax.set_xlim(tLimits)
    ax.set_ylim([1, ESpikes.N+ISpikes.N])
    if (yticks_style == 'separate'):
        ax.set_yticks([1, ESpikes.N, ESpikes.N+ISpikes.N])
    ax.invert_yaxis()
    ax.text(ylabelPos, 0.5, ylabel, va='center', ha='center',
            transform=ax.transAxes, rotation=90)
    if (not yticks):
        ax.yaxis.set_ticklabels([])
    
    return ax
