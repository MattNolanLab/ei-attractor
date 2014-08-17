#
#   sweeps.py
#
#   Functions/classes related to parameter sweeps.
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
import logging 

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.ticker as ti

from grid_cell_model.plotting.global_defs import globalAxesSettings
from grid_cell_model.plotting.low_level   import symmetricDataLimits
from base                 import createColorbar

from . import xlabelText, ylabelText
from . import aggregate as aggr
from .base import filterData

logger = logging.getLogger(__name__)
plotSweepLogger = logging.getLogger('{0}.{1}'.format(__name__, 'plotSweep'))

##############################################################################
# Parameter sweeps
def plotACTrial(sp, varList, iterList, noise_sigma, trialNumList=[0], **kw):
    #kw arguments
    r           = kw.pop('r', 0)
    c           = kw.pop('c', 0)
    cbar        = kw.pop('cbar', True)
    sigmaTitle  = kw.pop('sigmaTitle', True)
    annotations = kw.pop('annotations', None)

    if isinstance(sp, aggr.AggregateData):
        C, X, Y  = sp.getData()
    else:
        C = aggr.aggregate2DTrial(sp, varList, trialNumList)
        Y, X = aggr.computeYX(sp, iterList, r=r, c=c)
    C = ma.MaskedArray(C, mask=np.isnan(C))
    C, ax, cax = plot2DTrial(X, Y, C, colorBar=cbar, **kw)

    print("    max(C): {0}".format(np.max(C)))
    print("    min(C): {0}".format(np.min(C)))

    if (sigmaTitle):
        ax.set_title('$\sigma$ = {0} pA'.format(int(noise_sigma)))

    if (annotations is not None):
        for a in annotations:
            plotSweepAnnotation(X=X, Y=Y, ax=ax, **a)

    return C, ax, cax


def plotBumpSigmaTrial(aggregateData, noise_sigma, **kw):
    #kw arguments
    r           = kw.pop('r', 0)
    c           = kw.pop('c', 0)
    cbar        = kw.pop('cbar', True)
    kw['cmap']  = kw.get('cmap', 'jet')
    ignoreNaNs  = kw.pop('ignoreNaNs', False)
    sigmaTitle  = kw.pop('sigmaTitle', True)
    annotations = kw.pop('annotations', None)

    
    C, X, Y = aggregateData.getData()
    print("min(C): {0}".format(np.min(C)))
    print("max(C): {0}".format(np.max(C)))
    C, ax, cax =  plot2DTrial(X, Y, C, colorBar=cbar, **kw)

    if (sigmaTitle):
        ax.set_title('$\sigma$ = {0} pA'.format(int(noise_sigma)))

    if (annotations is not None):
        for a in annotations:
            plotSweepAnnotation(X=X, Y=Y, ax=ax, **a)

    return C, ax, cax



def plotBumpErrTrial(sp, varList, iterList, thr=np.infty, r=0, c=0, mask=None,
        trialNumList=[0], xlabel="", ylabel="", colorBar=True, clBarLabel="",
        vmin=None, vmax=None, title="", clbarNTicks=2, xticks=True,
        yticks=True):
    C = np.sqrt(aggr.aggregate2DTrial(sp, varList, trialNumList))
    if mask is None:
        mask = False
    C = ma.MaskedArray(C, mask=np.logical_or(np.logical_or(np.isnan(C), C >
        thr), mask))
    Y, X = aggr.computeYX(sp, iterList, r=r, c=c)
    return plot2DTrial(X, Y, C, xlabel, ylabel, colorBar, clBarLabel, vmin,
            vmax, title, clbarNTicks, xticks, yticks)

def plotFRTrial(sp, varList, iterList, noise_sigma, trialNumList=[0],
        thr=np.infty, **kw):
    r           = kw.pop('r', 0)
    c           = kw.pop('c', 0)
    ignoreNaNs  = kw.pop('ignoreNaNs', False)
    sigmaTitle  = kw.pop('sigmaTitle', True)
    cbar        = kw.pop('cbar', True)
    annotations = kw.pop('annotations', None)

    C = aggr.aggregate2DTrial(sp, varList, trialNumList,
            ignoreNaNs=ignoreNaNs)
    C = ma.MaskedArray(C, mask=np.logical_or(np.isnan(C), C > thr))
    Y, X = aggr.computeYX(sp, iterList, r=r, c=c)
    C, ax, cax = plot2DTrial(X, Y, C, colorBar=cbar, **kw)

    print("    max(C): {0}".format(np.max(C)))
    print("    min(C): {0}".format(np.min(C)))

    if (sigmaTitle):
        ax.set_title('$\sigma$ = {0} pA'.format(int(noise_sigma)))

    if (annotations is not None):
        for a in annotations:
            plotSweepAnnotation(X=X, Y=Y, ax=ax, **a)

    return C, ax, cax
            

def plotGridTrial(sp, varList, iterList, noise_sigma, trialNumList=[0], **kw):
    #kw arguments
    r           = kw.pop('r', 0)
    c           = kw.pop('c', 0)
    nansAs0     = kw.pop('nansAs0', False)
    ignoreNaNs  = kw.pop('ignoreNaNs', False)
    sigmaTitle  = kw.pop('sigmaTitle', True)
    cbar        = kw.pop('cbar', True)
    annotations = kw.pop('annotations', None)

    if isinstance(sp, aggr.AggregateData):
        G, X, Y = sp.getData()
    else:
        G = aggr.aggregate2DTrial(sp, varList, trialNumList, ignoreNaNs=ignoreNaNs)
        Y, X = aggr.computeYX(sp, iterList, r=r, c=c)

    nans = np.isnan(G)
    if (nansAs0):
        G[nans] = 0
    else:
        G = ma.MaskedArray(G, mask=nans)
    G, ax, cax = plot2DTrial(X, Y, G, colorBar=cbar, **kw)

    print("    max(G): {0}".format(np.max(G)))
    print("    min(G): {0}".format(np.min(G)))

    if (sigmaTitle):
        ax.set_title('$\sigma$ = {0} pA'.format(int(noise_sigma)))

    if (annotations is not None):
        for a in annotations:
            plotSweepAnnotation(X=X, Y=Y, ax=ax, **a)

    return G, ax, cax


def plotVelTrial(sp, varList, iterList, noise_sigma, **kw):
    # process kwargs
    r           = kw.pop('r', 0)
    c           = kw.pop('c', 0)
    sigmaTitle  = kw.pop('sigmaTitle', True)
    cbar        = kw.pop('cbar', True)
    annotations = kw.pop('annotations', None)

    C    = aggr.aggregate2D(sp, varList, funReduce = np.sum)
    C    = ma.MaskedArray(C, mask = np.isnan(C))
    Y, X = aggr.computeVelYX(sp, iterList, r=r, c=c)
    C, ax, cax = plot2DTrial(X, Y, C, colorBar=cbar, **kw)

    print("plotVelTrial: max(C): {0}".format(np.max(C.ravel())))
    print("plotVelTrial: min(C): {0}".format(np.min(C.ravel())))

    if (sigmaTitle):
        ax.set_title('$\sigma$ = {0} pA'.format(int(noise_sigma)))

    if (annotations is not None):
        for a in annotations:
            plotSweepAnnotation(X=X, Y=Y, ax=ax, **a)

    return C, ax, cax


def plotVelStdSweep(sp, iterList, noise_sigma, **kw):
    # process kwargs
    r          = kw.pop('r', 0)
    c          = kw.pop('c', 0)
    sigmaTitle = kw.pop('sigmaTitle', True)
    cbar       = kw.pop('cbar', True)


    varList = ['analysis', 'bumpVelAll']
    all = sp.aggregateData(varList, funReduce=None,
        trialNumList='all-at-once', output_dtype='list')
    C = np.ndarray((len(all), len(all[0])))
    for row in xrange(len(all)):
        for col in xrange(len(all[0])):
            slopes = all[row][col]
            if (isinstance(slopes, np.ndarray)):
                std = np.std(slopes, axis=0)
                C[row, col] = np.mean(std)
            else:
                C[row, col] = np.nan

    Y, X = aggr.computeVelYX(sp, iterList, r=r, c=c)
    C = ma.MaskedArray(C, mask = np.isnan(C))
    C, ax, cax = plot2DTrial(X, Y, C, colorBar=cbar, **kw)

    print("plotVelStdSweep: max(C): {0}".format(np.max(C.ravel())))
    print("plotVelStdSweep: min(C): {0}".format(np.min(C.ravel())))

    if (sigmaTitle):
        ax.set_title('$\sigma$ = {0} pA'.format(int(noise_sigma)))
    
    return C, ax, cax



def plotDiffTrial(spList, iterList, which, NTrials, types, **kw):
    r           = kw.pop('r', 0)
    c           = kw.pop('c', 0)
    cbar        = kw.pop('cbar', True)
    ignoreNaNs  = kw.pop('ignoreNaNs', False)
    annotations = kw.pop('annotations', None)
    symmetricLimits = kw.pop('symmetricLimits', True)
    filterThreshold = kw.pop('filterThreshold', -np.infty)


    if isinstance(spList[0], aggr.AggregateData):
        dataList = spList
        stackedData, X, Y = aggr.collapseNoiseAggregated(dataList)
    else:
        stackedData, X, Y = aggr.collapseNoise(spList, iterList, types,
                NTrials, ignoreNaNs=ignoreNaNs, normalizeTicks=True, r=r, c=c)
    stackedData, _ = filterData(stackedData, filterThreshold)
    diffData = np.diff(stackedData, axis=0)[which, :]
    if symmetricLimits:
        if 'vmin' in kw.keys() or 'vmax' in kw.keys():
            logger.warn('Overriding vmin and vmax by making data limits symmetric')
        kw['vmin'], kw['vmax'] = symmetricDataLimits(diffData)
        info_msg = 'vmin: {}, vmax: {}'.format(kw['vmin'], kw['vmax'])
        logger.info(info_msg)

    space0 = spList[0]
    nY, nX = Y.shape
    diffData, ax, cax = plot2DTrial(X, Y, np.reshape(diffData, (nY, nX)),
            colorBar=cbar, **kw)

    print("plotDiffTrial: max(diffData): {0}".format(np.max(diffData.ravel())))
    print("plotDiffTrial: min(diffData): {0}".format(np.min(diffData.ravel())))

    if (annotations is not None):
        for a in annotations:
            plotSweepAnnotation(X=X, Y=Y, ax=ax, **a)

    return diffData, ax, cax



def plotSweep(aggregateData, noise_sigma, **kw):
    cbar            = kw.pop('cbar', True)
    annotations     = kw.pop('annotations', None)
    #symmetricLimits = kw.pop('symmetricLimits', True)
    filterThreshold = kw.pop('filterThreshold', -np.infty)
    sigmaTitle      = kw.pop('sigmaTitle', True)

    data, X, Y = aggregateData.getData()
    plotSweepLogger.info("min(data): {0}".format(np.min(data)))
    plotSweepLogger.info("max(data): {0}".format(np.max(data)))
    plotSweepLogger.info("median(data): {0}".format(np.median(data)))

    data, ax, cax = plot2DTrial(X, Y, data, colorBar=cbar, **kw)

    if (sigmaTitle):
        ax.set_title('$\sigma$ = {0} pA'.format(int(noise_sigma)))

    if (annotations is not None):
        for a in annotations:
            plotSweepAnnotation(X=X, Y=Y, ax=ax, **a)

    return data, ax, cax



def plot2DTrial(X, Y, C, ax=plt.gca(), xlabel=None, ylabel=None,
        colorBar=False, clBarLabel="", vmin=None, vmax=None, title="",
        clbarNTicks=2, xticks=True, yticks=True, cmap=None, cbar_kw={},
        sliceAnn=None, **kw):
    # kw arguments
    kw['rasterized'] = kw.get('rasterized', True)
    # kw arguments (cbar)
    cbar_kw['label']       = cbar_kw.get('label', '')
    cbar_kw['shrink']      = cbar_kw.get('shrink', 0.8)
    #cbar_kw['orientation'] = cbar_kw.get('orientation', 'vertical')
    cbar_kw['pad']         = cbar_kw.get('pad', 0.05)
    cbar_kw['ticks']       = cbar_kw.get('ticks', ti.MultipleLocator(5))
    cbar_kw['rasterized']  = cbar_kw.get('rasterized', True)

    if xlabel is None:
        xlabel = xlabelText
    if ylabel is None:
        ylabel = ylabelText

    globalAxesSettings(ax)
    ax.minorticks_on()
    mappable = ax.pcolor(X, Y, C, vmin=vmin, vmax=vmax, cmap=cmap, **kw)
    cax = createColorbar(ax, mappable=mappable, **cbar_kw)
    if (colorBar == False):
        cax.set_visible(False)
    if (xlabel != ""):
        ax.set_xlabel(xlabel, va='top')
        ax.xaxis.set_label_coords(0.5, -0.125)
    if (ylabel != ""):
        ax.set_ylabel(ylabel, ha='center')
        ax.yaxis.set_label_coords(-0.125, 0.5)
    ax.xaxis.set_ticks([0, 6])
    ax.yaxis.set_ticks([0, 6])
    ax.xaxis.set_minor_locator(ti.AutoMinorLocator(6))
    ax.yaxis.set_minor_locator(ti.AutoMinorLocator(6))
    ax.axis('scaled')
    if (not xticks):
        ax.xaxis.set_ticklabels([])
    if (not yticks):
        ax.yaxis.set_ticklabels([])

    # slice annotations
    if (sliceAnn is not None):
        for args in sliceAnn:
            args.update(ax=ax, X=X, Y=Y)
            plotSliceAnnotation(**args)
        

    return C, ax, cax



def plotSweepAnnotation(txt, X, Y, rc, xytext_offset, **kw):
    # keyword args
    ax = kw.pop('ax', plt.gca())
    kw['color']      = kw.get('color', 'black')
    kw['arrowprops'] = kw.get('arrowprops',
            dict(arrowstyle         = "->",
                    linewidth       = plt.rcParams['lines.linewidth']*1.5,
                    color           = kw['color'],
                    connectionstyle = "arc3,rad=0"))

    r, c = rc[0], rc[1]
    xy   = (X[r, c], Y[r, c])
    x = xy[0]
    y = xy[1]
    xo = xytext_offset[0]
    yo = xytext_offset[1]
    
    ax.annotate(txt, (x, y), xytext=(x+xo, y+yo), xycoords='data',
        textcoords='data', va='center', fontweight='bold',
        zorder=10, **kw)


##############################################################################
# Parameter sweeps collapsed into noise-dependent lines

def plotCollapsedSweeps(noise_sigmas, data, **kw):
    ax          = kw.pop('ax', plt.gca())
    xlabel      = kw.pop('xlabel', "$\sigma_{noise}$")
    ylabel      = kw.pop('ylabel', '')
    xticks      = kw.pop('xticks', True)
    yticks      = kw.pop('yticks', True)
    #kw['color'] = kw.get('color', 'black')
    
    if (len(noise_sigmas) != len(data)):
        raise ValueError("len(noise_sigmas) != len(data)")

    stackedData = aggr.collapseSweeps(data)

    globalAxesSettings(ax)
    ax.plot(noise_sigmas, stackedData, "o-", **kw)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.margins(0.03, 0.03)
    ax.set_xticks(noise_sigmas)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if (not xticks):
        ax.xaxis.set_ticklabels([])
    if (not yticks):
        ax.yaxis.set_ticklabels([])

    return ax


class Contours(object):
    '''Contour object on a 2D plot. Usually a paramter sweep plot'''
    def __init__(self, data, V):
        '''Create the contour object from ``data``, with contour values
        ``V``.'''
        self.data = data
        self.V = V
    
    def plot(self, ax, V=None, **kwargs):
        '''Plot the contours into matplotlib axis.

        Parameters
        ----------
        ax : matplotlib.Axes
            Axes to plot into
        V : array-like
            A list of contour values to plot. If not None, the internal contour
            values will be overriden during plotting, but not inside the
            object.
        kwargs : dict
            Keyword arguments to pass on to the ax.contour() method.
        '''
        if V is None:
            V = self.V
        d, X, Y = self.data.getData()
        ax.contour(X, Y, d, V, **kwargs)
