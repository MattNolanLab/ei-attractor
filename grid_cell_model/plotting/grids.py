#
#   grids.py
#
#   Functions for grid cell related data.
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
import numpy.ma as ma
import matplotlib.pyplot as mpl
from matplotlib.pyplot  import plot, xlabel, ylabel, legend, xlim, ylim, \
        tight_layout, axis, title, pcolor, colorbar, hold, subplot, gca
from matplotlib.ticker  import MaxNLocator, LinearLocator

from global_defs         import globalAxesSettings, createColorbar
from plotting.low_level  import xScaleBar
from analysis.grid_cells import extractSpikePositions2D

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
        scaleText=True):
    '''
    Plot the grid-like rate map into the current axis
    '''
    if (ax is None):
        ax = gca()
    rateMap = ma.masked_array(rateMap, mask = np.sqrt(X**2 + Y**2) > diam/2.0)
    globalAxesSettings(ax)
    pcolor(X, Y, rateMap)
    axis('scaled')
    axis('off')
    title(titleStr, va='bottom')
    if (diam != np.inf):
        xlim([-lim_factor*diam/2.0, lim_factor*diam/2.0])
        ylim([-lim_factor*diam/2.0, lim_factor*diam/2.0])
    gridScaleBar(scaleBar, scaleText, ax)

 
def plotAutoCorrelation(ac, X, Y, diam=np.inf, ax=None, titleStr="",
        scaleBar=None, scaleText=True):
    if (ax is None):
        ax = gca()
    ac = ma.masked_array(ac, mask = np.sqrt(X**2 + Y**2) > diam)
    globalAxesSettings(ax)
    pcolor(X, Y, ac)
    axis('scaled')
    axis('off')
    title(titleStr, va='bottom')
    if (diam != np.inf):
        xlim([-lim_factor*diam, lim_factor*diam])
        ylim([-lim_factor*diam, lim_factor*diam])
    gridScaleBar(scaleBar, scaleText, ax)


