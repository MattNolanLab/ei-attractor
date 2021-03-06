'''Functions for plotting attractor bump related data.

Functions
---------

.. autosummary::

    bumpPosition
    torusFiringRate
    flatFiringRate
    plotBump
'''
from __future__ import absolute_import

import matplotlib.pyplot as mpl
import numpy as np
from matplotlib.pyplot  import plot, xlabel, ylabel, legend, ylim, \
        tight_layout, axis, title, pcolormesh, colorbar, hold, subplot
from matplotlib.ticker  import MaxNLocator, LinearLocator

from ..analysis.spikes import torusPopulationVector
from .global_defs      import globalAxesSettings, createColorbar


def bumpPosition(spikes, sheetSize, tstart, tend, dt, winLen, units="s"):
    '''
    Print a plot of bump position
    '''
    Ne_x = sheetSize[0]
    Ne_y = sheetSize[1]
    (pos, bumpPos_times) = torusPopulationVector(
            spikes,
            sheetSize,
            tstart,
            tend,
            dt,
            winLen)
    plot(bumpPos_times, pos)
    xlabel('Time (%s)' % units)
    ylabel('Bump position (neurons)')
    legend(['X', 'Y'], loc='upper left')
    ylim([-Ne_x/2 - 5, Ne_x/2 + 5])
    tight_layout()


def torusFiringRate(rateMap, labelx, labely=None, titleStr="", clbar=True,
        vMin=None, vMax=None, xTicks=True, yTicks=True, xTickLabels=True,
        yTickLabels=True):
    '''Firing rate of cells on the twisted torus'''
    if (labely is None):
        labely = labelx
    ax = mpl.gca()
    ax.minorticks_on()
    pcolormesh(rateMap, vmin=vMin, vmax=vMax)
    globalAxesSettings(ax)
    locx = LinearLocator(2)
    if (xTicks):
        ax.xaxis.set_ticks([1, len(rateMap[0])])
    else:
        ax.xaxis.set_ticks([])
    if (xTicks):
        ax.yaxis.set_ticks([1, len(rateMap)])
    else:
        ax.yaxis.set_ticks([])
    if (not xTickLabels):
        ax.xaxis.set_ticklabels([])
    if (not yTickLabels):
        ax.yaxis.set_ticklabels([])

    xlabel(labelx)
    ylabel(labely)
    if (clbar):
        cb = mpl.colorbar(ticks=MaxNLocator(4))
        cb.set_label('Firing rate (Hz)')
    axis('scaled')
    title(titleStr, va='bottom')


def flatFiringRate(FR, times, labely=None, labelx=None, units="ms",
        titleStr="", nticks=2):
    globalAxesSettings(mpl.gca())
    N = len(FR) # no. of rows == no. of neurons
    T, N_id = np.meshgrid(times, np.arange(N))
    pcolormesh(T, N_id,  FR)
    loc = LinearLocator(nticks)
    loc.tick_values(0, N-1)
    mpl.gca().yaxis.set_major_locator(loc)
    if (labelx is None):
        xlabel("Time (%s)" % units)
    elif (labelx != ""):
        xlabel(labelx)

    if (labely is None):
        ylabel("Neuron #")
    elif (labely != ""):
        ylabel(labely)

    if (titleStr != ""):
        title(titleStr)

    axis('tight')
    createColorbar(mpl.gca(), data=FR, label='Firing rate (Hz)')
    tight_layout()


def plotBump(ax, rateMap, cmap='jet', minRate=False, maxRate=True, **kw):
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    fs = kw.pop('fontsize', 'small')
    rx = kw.pop('rateXPos', 0.95)
    ry = kw.pop('rateYPos', 1.025)
    ax.pcolormesh(rateMap, cmap=cmap, **kw)
    ax.axis("scaled")
    ax.axis('off')
    if (minRate):
        minrStr = '{0:.1f} Hz'.format(np.min(rateMap.flatten()))
        ax.text(0.1, ry, minrStr, ha="left", va='bottom', fontsize=fs,
                transform=ax.transAxes)
    if (maxRate):
        maxrStr = '{0:.1f} Hz'.format(np.max(rateMap.flatten()))
        ax.text(rx, ry, maxrStr, ha="right", va='bottom', fontsize=fs,
                transform=ax.transAxes)


