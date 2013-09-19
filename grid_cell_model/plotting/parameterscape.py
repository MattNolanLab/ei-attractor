#
#   parameterscape.py
#
#   Parameterscape plots.
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
from matplotlib.cm import ScalarMappable

from global_defs import globalAxesSettings

def parameterScape(ax, X, Y, C, fmts, sizes, cmaps, vranges, **kwargs):
    """
      pcolor(C, **kwargs)
      pcolor(X, Y, C, **kwargs)

    """
    #alpha = kwargs.pop('alpha', None)
    #norm = kwargs.pop('norm', None)
    #cmap = kwargs.pop('cmap', None)
    #vmin = kwargs.pop('vmin', None)
    #vmax = kwargs.pop('vmax', None)
    #shading = kwargs.pop('shading', 'flat')
    sm = []
    colors = []
    for cIdx, c in enumerate(C):
        m = ScalarMappable(cmap=cmaps[cIdx])
        m.set_clim(vranges[cIdx])
        m.set_array(c)
        colors.append(m.to_rgba(c))
        sm.append(m)

    print("X : {}".format(X))
    print("Y : {}".format(Y))
    print("colors : {}".format(colors))

    for yIdx, y in enumerate(Y):
        for xIdx, x in enumerate(X):
            for cIdx, c in enumerate(C):
                if (np.isnan(c[yIdx, xIdx])):
                    continue
                ax.plot(x, y, fmts[cIdx], color=colors[cIdx][yIdx, xIdx, :],
                        markersize=sizes[cIdx], **kwargs)

    globalAxesSettings(ax)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    #ax.spines['bottom'].set_visible(False)
    #ax.spines['left'].set_visible(False)
    ax.margins(0.05, 0.05)

    #minx = np.amin(x)
    #maxx = np.amax(x)
    #miny = np.amin(y)
    #maxy = np.amax(y)

    #corners = (minx, miny), (maxx, maxy)
    #self.update_datalim( corners)
    #self.autoscale_view()
    #self.add_collection(collection)
    #return collection
