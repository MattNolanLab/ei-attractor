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
from matplotlib.pyplot  import plot, xlabel, ylabel, legend, ylim, \
        tight_layout, axis, title, pcolormesh, colorbar, hold, subplot, gca
from matplotlib.ticker  import MaxNLocator, LinearLocator

from global_defs        import globalAxesSettings, createColorbar
from plotting.low_level import xScaleBar



def plotGridRateMap(rateMap, X, Y, diam, ax=None, titleStr="", scaleBar=None):
    '''
    Plot the grid-like rate map into the current axis
    '''
    if (ax is None):
        ax = gca()
    rateMap = ma.masked_array(rateMap, mask = np.sqrt(X**2 + Y**2) > diam/2.0)
    globalAxesSettings(ax)
    pcolormesh(X, Y, rateMap)
    axis('scaled')
    axis('off')
    title(titleStr, va='bottom')
    if (scaleBar is not None):
        xScaleBar(scaleBar, ax, bottom=0, right=diam/2.0, unitsText='cm', size='small')

 
def plotAutoCorrelation(ac, X, Y, diam=np.inf, ax=None, titleStr="",
        scaleBar=None):
    if (ax is None):
        ax = gca()
    ac = ma.masked_array(ac, mask = np.sqrt(X**2 + Y**2) > diam)
    globalAxesSettings(ax)
    pcolormesh(X, Y, ac)
    axis('scaled')
    axis('off')
    title(titleStr, va='bottom')
    if (scaleBar is not None):
        xScaleBar(scaleBar, ax, bottom=0, right=diam, unitsText='cm', size='small')


