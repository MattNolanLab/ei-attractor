'''Global plotting definitions.

Functions
---------

.. autosummary::

    globalAxesSettings
    createColorbar
    prepareLims
'''
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, LinearLocator


def globalAxesSettings(ax, setTickPos=True):
    if (setTickPos):
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')


def createColorbar(ax, data=None, label="", nticks=2, **kw):
    #if (data is not None):
    #    mn = np.min(data.flatten())
    #    mx = np.max(data.flatten())
    #    ticks = np.linspace(mn, mx, nticks)
    #else:
    #    ticks = None
    cb = plt.colorbar(ax=ax, ticks=MaxNLocator(nticks+1), **kw)
    if (label != ""):
        cb.set_label(label)
    return cb


def prepareLims(lim, margin=0.01):
    w = lim[1] - lim[0]
    return (lim[0] - margin*w, lim[1] + margin*w)
