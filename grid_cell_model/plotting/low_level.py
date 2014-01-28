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
import matplotlib.patches    as patches

def xScaleBar(scaleLen, x, y, ax=m.gca(), height=0.02, color='black',
        unitsText='ms', size='medium', textYOffset=0.075):
    '''
    Plot a horizontal (X) scale bar into the axes.

    **Parameters:**
    scaleLen : float
        Size of the scale (X) in data coordinates.
    ax : mpl.axes.Axes
        Axes object. If unspecified, use the current axes.
    x  : float, optional
        Left end of the scale bar, in axis coordinates.
    y : float, optional
        Bottom position of the scale bar, in axis units. This excludes the scale
        text
    height : float, optional
        Height of the scale bar, in relative axis units.
    color
        Color of the bar.
    unitsText : string
        Units drawn below the scale bar.
    size 
        Size of the text below the scale bar.
    textYOffset : float
        Offset of the text from the scale bar. Positive value is a downward
        offset.
    '''
    (left, right) = ax.get_xlim()
    axisLen = scaleLen / (right - left)
    scaleCenter = x + 0.5*axisLen
    rect = patches.Rectangle((x,y), width=axisLen, height=height,
            transform=ax.transAxes, color=color)
    rect.set_clip_on(False)
    ax.add_patch(rect)
    if (unitsText is not None):
        textTemplate = '{0}'
        if (unitsText != ''):
            textTemplate += ' {1}'
        ax.text(scaleCenter, y - textYOffset,
                textTemplate.format(scaleLen, unitsText),
                va='top', ha='center', transform=ax.transAxes, size=size)


def yScaleBar(scaleLen, x, y, ax=m.gca(), width=0.0075, color='black',
        unitsText='ms', size='medium', textXOffset=0.075):
    '''
    Plot a vertical (Y) scale bar into the axes.

    **Parameters:**
    scaleLen : float
        Size of the scale (X) in data coordinates.
    ax : mpl.axes.Axes
        Axes object. If unspecified, use the current axes.
    x  : float, optional
        Left position of the scale bar, in axis coordinates. This excludes the
        scale text
    y : float, optional
        Bottom end of the scale bar, in axis units.
    width : float, optional
        Width of the scale bar, in relative axis units.
    color
        Color of the bar.
    unitsText : string
        Units drawn below the scale bar.
    size 
        Size of the text below the scale bar.
    textXOffset : float
        Offset of the text from the scale bar. Positive value is a leftward
        offset.
    '''
    (bottom, top) = ax.get_ylim()
    axisHeight = scaleLen / (top - bottom)
    scaleCenter = y + 0.5*axisHeight
    rect = patches.Rectangle((x,y), width=width, height=axisHeight,
            transform=ax.transAxes, color=color)
    rect.set_clip_on(False)
    ax.add_patch(rect)
    if (unitsText is not None):
        textTemplate = '{0}'
        if (unitsText != ''):
            textTemplate += ' {1}'
        ax.text(x - textXOffset, scaleCenter,
                textTemplate.format(scaleLen, unitsText),
                va='center', ha='left', transform=ax.transAxes, size=size,
                rotation=90)


def removeAllSpines(ax):
    for spine in ax.spines.itervalues():
        spine.set_visible(False)
