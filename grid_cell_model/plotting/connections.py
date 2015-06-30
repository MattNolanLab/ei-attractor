'''Functions/methods to plot connectivity matrices.

.. currentmodule:: grid_cell_model.plotting.connections

Classes / Functions
-------------------

.. autosummary::

    plotConnHistogram
    plot2DWeightMatrix
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ti

from global_defs import globalAxesSettings


def plotConnHistogram(val, **kw):
    '''
    Plot a histogram of connection weights with some predetermined formatting.

    Parameters
    ----------
    val : numpy array
        A 1D array of connetions.
    **kw
        Keyword arguments to the hist function. See code for some additional
        keyword arguments.
    '''
    # keyword arguments
    kw['bins']      = kw.get('bins', 20)
    #kw['edgecolor'] = kw.get('edgecolor', 'none')
    ax              = kw.pop('ax', plt.gca())
    xlabel          = kw.pop('xlabel', 'g (nS)')
    ylabel          = kw.pop('ylabel', 'Count')
    title           = kw.pop('title', '')
    locators        = kw.pop('locators', {})
    ylabelPos       = kw.pop('ylabelPos', -0.2)

    globalAxesSettings(ax)
    n, _, _ = ax.hist(val, **kw)
    print(np.max(n))
    ax.set_xlabel(xlabel)
    ax.text(ylabelPos, 0.5, ylabel, rotation=90, transform=ax.transAxes,
            va='center', ha='right')
    ax.set_title(title)

    # tick formatting
    x_major = locators.get('x_major', ti.MaxNLocator(3))
    x_minor = locators.get('x_minor', ti.AutoMinorLocator(2))
    y_major = locators.get('y_major', ti.LinearLocator(2))
    y_minor = locators.get('y_minor', ti.AutoMinorLocator(4))
    ax.xaxis.set_major_locator(x_major)
    ax.yaxis.set_major_locator(y_major)
    ax.xaxis.set_minor_locator(x_minor)
    ax.yaxis.set_minor_locator(y_minor)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    return ax


def plot2DWeightMatrix(C, **kw):
    ax       = kw.pop('ax', plt.gca())
    xlabel   = kw.pop('xlabel', 'Neuron #')
    ylabel   = kw.pop('ylabel', 'Neuron #')
    labelpad = kw.pop('labelpad', None)
    title    = kw.pop('title', '')
    kw['rasterized'] = kw.get('rasterized', True)

    X = np.arange(C.shape[1] + 1)
    Y = np.arange(C.shape[0] + 1)
    globalAxesSettings(ax)
    ax.pcolormesh(X, Y, C, **kw)
    ax.set_xticks([0.5, X[-1]-0.5])
    ax.set_yticks([0.5, Y[-1]-0.5])
    ax.set_xticklabels([1, X[-1]])
    ax.set_yticklabels([1, Y[-1]])
    plt.axis('scaled')

    ax.set_xlabel(xlabel, labelpad=labelpad)
    ax.set_ylabel(ylabel, labelpad=labelpad)
    ax.set_title(title, size='small')

    return ax


