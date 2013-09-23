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
import matplotlib
import matplotlib.cm as cm
from matplotlib.cm import ScalarMappable
from matplotlib.transforms import Bbox
from matplotlib.ticker import MaxNLocator

from global_defs import globalAxesSettings

def parameterScape(fig, X, Y, C, fmts, sizes, cmaps, vranges, **kwargs):
    """
    """
    ax = fig.gca()

    defaultMargin = 1.0
    xmargin = kwargs.pop('xmargin', defaultMargin)
    ymargin = kwargs.pop('ymargin', defaultMargin)
    legend  = kwargs.pop('legend', None)
    colorbars = kwargs.pop('colorbars', None)

    sm = []
    colors = []
    for cIdx, c in enumerate(C):
        m = ScalarMappable(cmap=cmaps[cIdx])
        m.set_clim(vranges[cIdx])
        m.set_array(c)
        colors.append(m.to_rgba(c))
        sm.append(m)


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
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(bottom='off', top='off', left='off', right='off', length=0,
            pad=-5)
    ax.axis('scaled')
    dx = X[1] - X[0]
    dy = Y[1] - Y[0]
    ax.set_xlim([X[0] - xmargin*dx, X[-1] + xmargin*dx])
    ax.set_ylim([Y[0] - ymargin*dy, Y[-1] + ymargin*dy])

    if (legend is not None):
        legend_strings = legend[0]
        legend_kw = legend[1]
        if (legend_kw is None):
            legend_kw = {}
        parameterScapeLegend(ax, legend_strings, fmts, sizes, cmaps, kwargs,
                **legend_kw)

    if (colorbars is not None):
        cax = plotColorbars(fig, sm, **colorbars)
    else:
        cax = None

    return ax, cax


def parameterScapeLegend(ax, leg, fmts, sizes, cmaps, plot_kwargs, **kwargs):
    x = kwargs.pop('x', -0.1)
    y = kwargs.pop('y', 1.1)
    sizeScaleFactor = kwargs.pop('sizeScaleFactor', 1.5)
    plot_kwargs['clip_on'] = False
    colors = kwargs.pop('colors', np.linspace(0, 1, len(leg)))
    print sizes

    for idx, l in enumerate(leg):
        color = cm.get_cmap(cmaps[idx])(colors[idx])
        print color
        ax.plot(x, y, fmts[idx], markersize=sizes[idx]*sizeScaleFactor,
                color=color, transform=ax.transAxes, **plot_kwargs)


def plotColorbars(fig, mappables, **kw):
    labels = kw.pop('labels')
    height = kw.pop('height', 0.01)
    width  = kw.pop('width', 0.5)
    x      = kw.pop('x', 0.5)
    y      = kw.pop('y', 0.9)
    div    = kw.pop('div', 0.05)
    orientation = kw.get('orientation', 'horizontal')
    kw['orientation'] = orientation
    labelsize = kw.pop('labelsize', None)
    tickLocators = kw.pop('tickLocators', None)

    plotX = x
    plotY = y
    cbarAxes = []
    for idx, m in enumerate(mappables):
        cax = fig.add_axes(Bbox.from_bounds(plotX, plotY, width, height))
        #globalAxesSettings(cax)
        cbarAxes.append(cax)
        cb = fig.colorbar(m, cax=cax, ticks=tickLocators[idx], **kw)
        cb.set_label(labels[idx], size=labelsize)
        if (orientation == 'horizontal'):
            plotY += div
        elif (orientation == 'vertical'):
            plotX += div
        else:
            raise ValueError('Orientation must be either `horizontal` or `vertical`')


    return cbarAxes


