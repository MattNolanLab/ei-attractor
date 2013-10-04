#!/usr/bin/env python
#
#   connections.py
#
#   Functions/methods to plot connectivity matrices.
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

from global_defs import globalAxesSettings


def plot2DWeightMatrix(C, **kw):
    ax     = kw.pop('ax', plt.gca())
    xlabel = kw.pop('xlabel', 'X neuron #')
    ylabel = kw.pop('ylabel', 'Y neuron #')
    title  = kw.pop('title', '')

    X = np.arange(C.shape[1] + 1)
    Y = np.arange(C.shape[0] + 1)
    globalAxesSettings(ax)
    ax.pcolormesh(X, Y, C, **kw)
    ax.set_xticks([0.5, X[-1]-0.5])
    ax.set_yticks([0.5, Y[-1]-0.5])
    ax.set_xticklabels([1, X[-1]])
    ax.set_yticklabels([1, Y[-1]])
    plt.axis('scaled')

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    return ax


