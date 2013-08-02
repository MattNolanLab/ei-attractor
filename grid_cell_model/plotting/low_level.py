#
#   low_level.py
#
#   Low level plotting functions.
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
import matplotlib.pyplot     as m
import matplotlib.transforms as transforms

def xScaleBar(scaleLen, ax=m.gca(), height=0.02, bottom=0.05, right=None,
        color='black', unitsText='ms', size='medium'):
    '''
    Plot a horizontal (X) scale bar into the axes.

    Parameters
    ----------
    scaleLen : float
        Size of the scale (X) in data coordinates.
    ax : mpl.axes.Axes
        Axes object. If unspecified, use the current axes.
    height : float, optional
        Height of the scale bar, in relative axis units.
    bottom : float, optional
        Bottom position of the scale bar, in axis units.
    right  : float, optional
        Right end of the scale bar, in data coordinates. If unspecified, use
        the right most axis limit at the time of drawing this bar.
    color
        Color of the bar.
    unitesText : string
        Units drawn below the scale bar.
    size 
        Size of the text below the scale bar.
    '''
    if (right is None):
        (left, right) = ax.get_xlim()
    scaleCenter = right - 0.5*scaleLen
    ax.axvspan(
            xmin = right - scaleLen,
            xmax = right,
            ymin = bottom,
            ymax = bottom + height,
            color = color)
    trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
    textTemplate = '{0}'
    if (unitsText != ''):
        textTemplate += ' {1}'
    ax.text(scaleCenter, bottom - 0.05,
            textTemplate.format(scaleLen, unitsText),
            va='top', ha='center', transform=trans, size=size)


def yScaleBar(scaleLen, ax=m.gca(), width=0.01, bottom=None, right=0.8,
        color='black', unitsText='ms', size='medium'):
    '''
    Plot a vertical (Y) scale bar into the axes.
    '''
    if (bottom is None):
        (bottom, top) = ax.get_ylim()
    scaleCenter = bottom + 0.5*scaleLen
    ax.axhspan(
            xmin = right - width,
            xmax = right,
            ymin = bottom,
            ymax = bottom + scaleLen,
            color = color)
    trans = transforms.blended_transform_factory(ax.transAxes, ax.transData)
    textTemplate = '{0}'
    if (unitsText != ''):
        textTemplate += ' {1}'
    ax.text(right + 0.05, scaleCenter,
            textTemplate.format(scaleLen, unitsText),
            va='center', ha='left', transform=trans, size=size)

