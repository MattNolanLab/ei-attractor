#
#   plotting.py
#
#   Shared plotting functions.
#
#       Copyright (C) 2012  Lukas Solanka <l.solanka@sms.ed.ac.uk>
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
from matplotlib.ticker  import MaxNLocator, LinearLocator, AutoMinorLocator

from plotting.global_defs import globalAxesSettings, createColorbar

def plot2DTrial(X, Y, C, xlabel="", ylabel="",
        colorBar=True, clBarLabel="", vmin=None, vmax=None, title="",
        clbarNTicks=2, xticks=True, yticks=True):

    ax = plt.gca()
    globalAxesSettings(ax)
    ax.minorticks_on()
    plt.pcolor(X, Y, C, vmin=vmin, vmax=vmax)
    if (colorBar):
        if (clbarNTicks == None):
            createColorbar(ax, None, clBarLabel, orientation='horizontal',
                    pad=0.2)
        else:
            createColorbar(ax, C, clBarLabel, nticks=clbarNTicks,
                    orientation='horizontal', pad=0.2)
    if (xlabel != ""):
        plt.xlabel(xlabel, va='top')
        ax.xaxis.set_label_coords(0.5, -0.15)
    if (ylabel != ""):
        plt.ylabel(ylabel, ha='right')
        ax.yaxis.set_label_coords(-0.125, 0.5)
    ax.xaxis.set_ticks([0, 6])
    ax.yaxis.set_ticks([0, 3])
    ax.xaxis.set_minor_locator(AutoMinorLocator(6))
    ax.yaxis.set_minor_locator(AutoMinorLocator(3))
    plt.axis('scaled')
    if (not xticks):
        ax.xaxis.set_ticklabels([])
    if (not yticks):
        ax.yaxis.set_ticklabels([])

    return C


