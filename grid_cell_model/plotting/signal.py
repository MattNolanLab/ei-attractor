#
#   signal.py
#
#   Signal plotting functions
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
import matplotlib.pyplot as plt
from matplotlib.ticker  import MaxNLocator

from global_defs import globalAxesSettings


def signalPlot(t, sig, ax, labelx=None, labely="", timeUnits="ms",
        leg=['middle', 'edge'], nticks=3):
    '''
    Plot a customized signal plot into an axis.

    Parameters
    ----------
    t : numpy array
        Time array
    sig : numpy array
        Signal values
    ax : matplotlib.axes object
        Axes to plot into
    labelx : string
        X label
    labely : string
        Y label
    timeUnits : string
        Time units to be appended to xlabe
    leg : list of strings or None
        Legend to print. If None, do not print any legend
    nticks : int
        Number of ticks on the yaxis
    '''
    if (labelx is None):
        labelx = 'Time (%s)' % timeUnits
    plt.hold('on')
    globalAxesSettings(ax)
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(nticks-1))
    plt.plot(t, sig)
    plt.xlabel(labelx)
    plt.ylabel(labely)
    if (leg is not None):
        plt.legend(leg)


def EIPlot(E, I, labelx=None, labely="", holdVal='on', timeUnits="ms",
        leg=['middle', 'edge'], nticks=3):
    '''
    Plot a double plot of E and I variable.

    *E*
        A pair (sig, times)
    *I*
        A pair (sig, times)
    *holdVal*
        value to pass on to the hold() function
    '''

    plt.hold('on')
    ax = plt.subplot(211)
    globalAxesSettings(plt.gca())
    plt.gca().xaxis.set_major_locator(MaxNLocator(4))
    plt.gca().yaxis.set_major_locator(MaxNLocator(nticks-1))
    plt.plot(E[1], E[0])
    if (labely != ""):
        plt.ylabel("E cell " + labely)
    plt.legend(leg)

    plt.subplot(212, sharex=ax)
    globalAxesSettings(plt.gca())
    plt.gca().xaxis.set_major_locator(MaxNLocator(4))
    plt.gca().yaxis.set_major_locator(MaxNLocator(nticks-1))
    plt.plot(I[1], I[0])
    if (labely != ""):
        plt.ylabel("I cell " + labely)
    if (labelx is None):
        plt.xlabel('Time (%s)' % timeUnits)
    elif (labelx != ""):
        plt.xlabel(labelx)
    plt.legend(leg)

    plt.tight_layout(w_pad=1.0)


