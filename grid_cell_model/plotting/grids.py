'''Functions for grid cell related data.'''
from __future__ import absolute_import
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as mpl
from matplotlib.pyplot  import plot, xlabel, ylabel, legend, xlim, ylim, \
        tight_layout, axis, title, pcolor, colorbar, hold, subplot, gca
from matplotlib.ticker  import MaxNLocator, LinearLocator

from .global_defs         import globalAxesSettings, createColorbar
from .low_level  import xScaleBar
from ..analysis.grid_cells import extractSpikePositions2D

lim_factor = 1.1

def gridScaleBar(scaleLen, scaleText, ax):
    if (scaleLen is not None):
        if (scaleText):
            unitsText = 'cm'
        else:
            unitsText = None
        xScaleBar(scaleLen, x=0.7, y=-0.00, ax=ax, height=0.015,
                unitsText=unitsText, size='small')

def plotSpikes2D(spikeTimes, rat_pos_x, rat_pos_y, dt, diam=np.inf, ax=None,
        titleStr='', scaleBar=None, scaleText=True, spikeDotSize=5):
    '''
    Plot spike positions into the figure. Both positions and spikes must be aligned!
    '''
    if (ax is None):
        ax = gca()
    neuronPos_x, neuronPos_y, m_i = extractSpikePositions2D(spikeTimes, rat_pos_x, rat_pos_y, dt)

    plot(rat_pos_x[0:m_i], rat_pos_y[0:m_i])
    hold('on')
    plot(neuronPos_x, neuronPos_y, 'or', markersize=spikeDotSize)
    axis('off')
    axis('scaled')
    title(titleStr, va='bottom')
    if (diam != np.inf):
        xlim([-lim_factor*diam/2.0, lim_factor*diam/2.0])
        ylim([-lim_factor*diam/2.0, lim_factor*diam/2.0])
    gridScaleBar(scaleBar, scaleText, ax)



def plotGridRateMap(rateMap, X, Y, diam, ax=None, titleStr="", scaleBar=None,
        scaleText=True, maxRate=True, G=None, rateStr='Hz', ann_div=.1, **kw):
    '''
    Plot the grid-like rate map into the current axis
    '''
    vmin = kw.pop('vmin', 0)
    if (ax is None):
        ax = gca()
    rateMap = ma.masked_array(rateMap, mask = np.sqrt(X**2 + Y**2) > diam/2.0)
    globalAxesSettings(ax)
    ax.imshow(rateMap, interpolation='none',
              extent=(X[0, 0], X[0, -1], Y[0, 0], Y[-1, 0]),
              origin='lower',
              rasterized=False)
    ax.axis('scaled')
    ax.axis('off')
    ax.set_title(titleStr, va='bottom')
    if (diam != np.inf):
        ax.set_xlim([-lim_factor*diam/2.0, lim_factor*diam/2.0])
        ax.set_ylim([-lim_factor*diam/2.0, lim_factor*diam/2.0])
    gridScaleBar(scaleBar, scaleText, ax)
    if (maxRate):
        if rateStr != '':
            rateStr = ' ' + rateStr
        rStr = '{0:.1f}{1}'.format(np.max(rateMap.flatten()), rateStr)
        ax.text(1. - ann_div, 1.025, rStr, ha="right", va='bottom',
                fontsize='xx-small', transform=ax.transAxes)
    if (G is not None):
        if (int(G*100)/100.0 == int(G)):
            gStr = '{0}'.format(int(G))
        else:
            gStr = '{0:.2f}'.format(G)
        ax.text(ann_div, 1.025, gStr, ha="left", va='bottom',
                fontsize='xx-small', transform=ax.transAxes)



def plotAutoCorrelation(ac, X, Y, diam=np.inf, ax=None, titleStr="",
        scaleBar=None, scaleText=True, **kw):
    if (ax is None):
        ax = gca()
    ac = ma.masked_array(ac, mask = np.sqrt(X**2 + Y**2) > diam)
    globalAxesSettings(ax)
    ax.pcolormesh(X, Y, ac, **kw)
    ax.axis('scaled')
    ax.axis('off')
    ax.set_title(titleStr, va='bottom')
    if (diam != np.inf):
        ax.set_xlim([-lim_factor*diam, lim_factor*diam])
        ax.set_ylim([-lim_factor*diam, lim_factor*diam])
    gridScaleBar(scaleBar, scaleText, ax)


