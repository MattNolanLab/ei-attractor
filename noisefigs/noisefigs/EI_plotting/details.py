'''Plotting parameter sweeps of the detailed noise levels.

.. currentmodule:: noisefigs.EI_plotting.details

Functions used in [SOLANKA2015]_
--------------------------------

.. autosummary::

    plotDetailedNoise
'''
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.ticker as ti

from . import xlabelText, ylabelText
from . import aggregate as aggr
from grid_cell_model.plotting.global_defs import globalAxesSettings


def decideLabels(type):
    if (type == 'horizontal'):
        xlabel = xlabelText
        titleText = "$g_E$"
    elif (type == 'vertical'):
        xlabel = ylabelText
        titleText = "$g_I$"
    else:
        raise ValueError("type must be 'horizontal' or 'vertical'")

    titles = dict(xlabel=xlabel, titleText=titleText)
    return titles


def plotOneSlice(ax, x, y, **kw):
    xlabel           = kw.pop('xlabel', '')
    ylabel           = kw.pop('ylabel', '')
    ylabelPos        = kw.pop('ylabelPos', -0.2)
    fmt              = kw.pop('fmt', 'o-')
    xticks           = kw.pop('xticks', True)
    yticks           = kw.pop('yticks', True)
    errors           = kw.pop('errors', True)
    kw['markersize'] = kw.get('markersize', 4)

    globalAxesSettings(ax)
    ndim = len(y.shape)
    if (ndim == 2):
        mean = np.mean(y, axis=1) # axis 1: trials
        std  = np.std(y, axis=1)
        if (errors):
            ax.errorbar(x, mean, std, fmt=fmt, **kw)
        else:
            ax.plot(x, mean, fmt, **kw)
    elif (ndim == 1):
        ax.plot(x, y, fmt, **kw)
    ax.set_xlabel(xlabel)
    ax.text(ylabelPos, 0.5, ylabel, rotation=90, transform=ax.transAxes,
            va='center', ha='center')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_major_locator(ti.MultipleLocator(1))
    ax.xaxis.set_minor_locator(ti.AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(ti.AutoMinorLocator(3))
    w = x[-1] - x[0]
    margin = 0.025
    ax.set_xlim([-margin*w, x[-1]+margin*w])

    if (not xticks):
        ax.xaxis.set_ticklabels([])
    if (not yticks):
        ax.yaxis.set_ticklabels([])


def extractSliceX(X, Y, rowSlice, colSlice, type):
    if (type == 'horizontal'):
        return X[0, colSlice]
    else:
        return Y[rowSlice, 0]


def collapseTrials(var, type, ignoreNaNs):
    if (len(var.shape) != 3):
        return var

    trials = []
    for trialIdx in xrange(var.shape[2]):
        trials.append(var[:, :, trialIdx])

    res    = None
    if (type == 'horizontal'):
        res = np.vstack(trials).T
    elif (type == 'vertical'):
        res = np.hstack(trials)
    else:
        raise ValueError("type must be 'horizontal' or 'vertical'")

    if (ignoreNaNs):
        nans = np.isnan(res)
        res = ma.MaskedArray(res, mask=nans)

    return res


def plotGridnessSlice(paramSpaces, rowSlice, colSlice, type, NTrials=1, **kw):
    # kwargs
    title        = kw.pop('title', True)
    rcG          = kw.pop('rowsCols', [(1, 22), (1, 22), (1, 22)]) # (row, col)
    ax           = kw.pop('ax', plt.gca())
    iterList     = kw.pop('iterList', ['g_AMPA_total', 'g_GABA_total'])
    kw['ylabel'] = kw.get('ylabel', 'Gridness score')
    labels = decideLabels(type)
    kw['xlabel'] = kw.get('xlabel', labels['xlabel'])
    ignoreNaNs   = kw.get('ignoreNaNs', True)

    GVars = ['analysis', 'gridnessScore']
    trialNumList = range(NTrials)
    sp = paramSpaces

    # Gridness score
    for idx, noise_sigma in enumerate(sp.noise_sigmas):
        space = sp.grids[idx]
        G = space.aggregateData(GVars, trialNumList, output_dtype='array',
                loadData=True, saveData=False)
        if (ignoreNaNs):
            nans = np.isnan(G)
            G = ma.MaskedArray(G, mask=nans)
        Y, X = aggr.computeYX(space, iterList, r=rcG[idx][0], c=rcG[idx][1])
        x = extractSliceX(X, Y, rowSlice, colSlice, type)
        y = collapseTrials(G[rowSlice, colSlice, :], type, ignoreNaNs)
        #import pdb; pdb.set_trace()
        plotOneSlice(ax, x, y, **kw)
    ax.yaxis.set_major_locator(ti.MaxNLocator(4))

    # Title
    if (title):
        if (type == 'horizontal'):
            sliceG = Y[rowSlice, 0]
        else:
            sliceG = X[0, colSlice]
        if (isinstance(sliceG, np.ndarray)):
            txt = '{0} = {1} - {2} nS'.format(labels['titleText'], sliceG[0],
                    sliceG[-1])
        else:
            txt = '{0} = {1} nS'.format(labels['titleText'], sliceG)
        ax.set_title(txt, x=0.99, y=0.92, va='bottom', ha='right',
                fontsize='small')

    return ax


def plotSliceAnnotation(ax, X, Y, sliceSpan, type, **kw):
    # kw arguments
    kw['color']     = kw.get('color', 'white')
    kw['linewidth'] = kw.get('linewidth', 1.0)
    letter      = kw.pop('letter', None)
    letterColor = kw.pop('letterColor', 'white')

    if (type == 'horizontal'):
        pos0  = Y[sliceSpan.start,0]
        pos1  = Y[sliceSpan.stop, 0]
        ax.axhline(y=pos0, **kw)
        ax.axhline(y=pos1, **kw)
        if (letter is not None):
            width = X[0, -1] - X[0, 0]
            ax.text(0.9*width, (pos1 + pos0)/2.0, letter, color=letterColor,
                    ha='center', va='center', weight='bold', zorder=2)
    elif (type == 'vertical'):
        pos0 = X[0, sliceSpan.start]
        pos1 = X[0, sliceSpan.stop]
        ax.axvline(x=pos0, **kw)
        ax.axvline(x=pos1, **kw)
        if (letter is not None):
            height = Y[-1, 0] - Y[0, 0]
            ax.text((pos1+pos0)/2.0, 0.9*height, letter, color=letterColor,
                    ha='center', va='center', weight='bold', zorder=2)

    else:
        raise ValueError()


def plotDetailedNoise(sp, NTrials, types, **kw):
    xlabel                = kw.pop('xlabel', '$\sigma_{noise}$ (pA)')
    ylabel                = kw.pop('ylabel', '')
    ylabelPos             = kw.pop('ylabelPos', -0.4)
    xticks                = kw.pop('xticks', True)
    yticks                = kw.pop('yticks', True)
    ax                    = kw.pop('ax', plt.gca())
    markersize            = kw.pop('markersize', 4)
    kw['markerfacecolor'] = kw.get('markerfacecolor', 'none')
    color                 = kw.pop('color', 'blue')
    ignoreNaNs            = kw.pop('ignoreNaNs', True)

    iterList = ['noise_sigma', 'g_AMPA_total']

    ax.hold('on')
    if isinstance(sp, aggr.AggregateData):
        data, X, Y = sp.getData()
    else:
        data, X, Y = aggr.aggregateType(sp, iterList, types, NTrials, **kw)
    if (ignoreNaNs):
        nans = np.isnan(data)
        data = ma.MaskedArray(data, mask=nans)
    print("    max(data): {0}".format(np.max(data)))
    print("    min(data): {0}".format(np.min(data)))
    noise_sigma = Y[:, 0]
    mean = np.mean(data, axis=1)

    globalAxesSettings(ax)
    noise_sigma_all = np.repeat(np.reshape(noise_sigma, (1, len(noise_sigma))),
            data.shape[1], axis=0).T
    p1, = ax.plot(noise_sigma_all.ravel(), data.ravel(), 'o', markeredgecolor=color,
            markersize=markersize, **kw)
    l1, = ax.plot(noise_sigma, mean, '-', color=color, markersize=markersize, **kw)
    ax.set_xlabel(xlabel)
    ax.text(ylabelPos, 0.5, ylabel, va='center', ha='center',
            transform=ax.transAxes, rotation=90)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if (not yticks):
        ax.yaxis.set_ticklabels([])
    if (not xticks):
        ax.xaxis.set_ticklabels([])
    ax.margins(0.03, 0.03)

    return ax, p1, l1



