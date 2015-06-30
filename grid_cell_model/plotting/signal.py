'''Signal plotting functions.'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ti
from matplotlib import rcParams as rcp

from global_defs import globalAxesSettings


def signalPlot(t, sig, ax, timeUnits="ms", nticks=3, nThetaTicks=None, **kw):
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
    xlabel : string
        X label
    ylabel : string
        Y label
    timeUnits : string
        Time units to be appended to xlabe
    leg : list of strings or None
        Legend to print. If None, do not print any legend
    nticks : int
        Number of ticks on the yaxis
    '''
    xlabel    = kw.pop('xlabel', 'Time (%s)' % timeUnits)
    ylabel    = kw.pop('ylabel', '')
    ylabelPos = kw.pop('ylabelPos', -0.22)
    xmargin   = kw.pop('xmargin', 0)
    zeroLine  = kw.pop('zeroLine', True)

    plt.hold('on')
    globalAxesSettings(ax)
    if nThetaTicks is not None:
        xticks = np.linspace(t[0], t[-1], nThetaTicks)
        ax.set_xticks(xticks)
        plt.grid(b=True, which='major', axis='x')
    else:
        ax.xaxis.set_major_locator(ti.MaxNLocator(4))
    ax.yaxis.set_major_locator(ti.MaxNLocator(nticks-1))
    ax.xaxis.set_minor_locator(ti.AutoMinorLocator(2))
    ax.plot(t, sig, **kw)
    ax.set_xlabel(xlabel)
    ax.text(ylabelPos, 0.5, ylabel, va='center', ha='center',
            transform=ax.transAxes, rotation=90)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    w = t[-1] - t[0]
    ax.set_xlim([t[0] - xmargin*w, t[-1] + xmargin*w])

    if (zeroLine):
        color  = rcp['grid.color']
        ls     = rcp['grid.linestyle']
        lw     = rcp['grid.linewidth']
        alphsa = rcp['grid.alpha']
        ax.axhline(0, ls=ls, lw=lw, color=color)


def EIPlot(E, I, labelx=None, labely="", holdVal='on', timeUnits="ms",
        leg=['middle', 'edge'], nticks=3):
    '''
    Plot a double plot of E and I variable.

    Parameters
    ----------
    E : A pair (sig, times)
        Excitatory variable.
    I : A pair (sig, times)
        Inhibitory variable.
    holdVal : str
        Value to pass on to the hold() function
    '''

    plt.hold('on')
    ax = plt.subplot(211)
    globalAxesSettings(plt.gca())
    plt.gca().xaxis.set_major_locator(ti.MaxNLocator(4))
    plt.gca().yaxis.set_major_locator(ti.MaxNLocator(nticks-1))
    plt.plot(E[1], E[0])
    if (labely != ""):
        plt.ylabel("E cell " + labely)
    plt.legend(leg)

    plt.subplot(212, sharex=ax)
    globalAxesSettings(plt.gca())
    plt.gca().xaxis.set_major_locator(ti.MaxNLocator(4))
    plt.gca().yaxis.set_major_locator(ti.MaxNLocator(nticks-1))
    plt.plot(I[1], I[0])
    if (labely != ""):
        plt.ylabel("I cell " + labely)
    if (labelx is None):
        plt.xlabel('Time (%s)' % timeUnits)
    elif (labelx != ""):
        plt.xlabel(labelx)
    plt.legend(leg)

    plt.tight_layout(w_pad=1.0)



###############################################################################
## Theta-gamma
#def plotThetaSignal(t, theta, **kw):
#    # kw arguments
#    ax = kw.pop('ax', plt.gca())
#    ylabel = kw.pop('ylabel', 'I (pA)')
#    ylabelPos   = kw.pop('ylabelPos', None)
#
#    globalAxesSettings(ax)
#    ax.plot(t, theta, **kw)
#    if (ylabelPos is None):
#        ax.set_ylabel(ylabel)
#    else:
#        ax.text(ylabelPos, 0.5, ylabel, rotation=90, transform=ax.transAxes,
#                va='center', ha='center')
#
