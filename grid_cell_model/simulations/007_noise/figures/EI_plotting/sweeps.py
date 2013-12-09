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
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.ticker as ti

from plotting.global_defs import globalAxesSettings
from figures_shared       import createColorbar
from plotting.colors      import Colormap2D

from . import xlabelText, ylabelText, xlabelTextShort, ylabelTextShort
from . import aggregate as aggr

##############################################################################
# Parameter sweeps
def plotACTrial(sp, varList, iterList, noise_sigma, trialNumList=[0], **kw):
    #kw arguments
    r           = kw.pop('r', 0)
    c           = kw.pop('c', 0)
    cbar        = kw.pop('cbar', True)
    sigmaTitle  = kw.pop('sigmaTitle', True)
    annotations = kw.pop('annotations', None)

    C = aggr.aggregate2DTrial(sp, varList, trialNumList)
    C = ma.MaskedArray(C, mask=np.isnan(C))
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


def plotBumpSigmaTrial(sp, varList, iterList, noise_sigma, trialNumList=[0],
        thr=np.infty, **kw):
    #kw arguments
    r           = kw.pop('r', 0)
    c           = kw.pop('c', 0)
    cbar        = kw.pop('cbar', True)
    kw['cmap']  = kw.get('cmap', 'jet')
    ignoreNaNs  = kw.pop('ignoreNaNs', False)
    sigmaTitle  = kw.pop('sigmaTitle', True)
    annotations = kw.pop('annotations', None)

    C = aggr.aggregate2DTrial(sp, varList, trialNumList, ignoreNaNs=ignoreNaNs)
    C = ma.MaskedArray(C, mask=np.logical_or(np.isnan(C), C > thr))
    C = 1. / C
    print("min(C): {0}".format(np.min(C)))
    print("max(C): {0}".format(np.max(C)))
    Y, X = aggr.computeYX(sp, iterList, r=r, c=c)
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

    G = aggr.aggregate2DTrial(sp, varList, trialNumList, ignoreNaNs=ignoreNaNs)
    nans = np.isnan(G)
    if (nansAs0):
        G[nans] = 0
    else:
        G = ma.MaskedArray(G, mask=nans)
    Y, X = aggr.computeYX(sp, iterList, r=r, c=c)
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

    C    = np.abs(aggr.aggregate2D(sp, varList, funReduce = np.sum))
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



def plot2DTrial(X, Y, C, ax=plt.gca(), xlabel=xlabelText, ylabel=ylabelText,
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
                    linewidth       = 1.5,
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
# Scatter plots
class ScatterPlot(object):

    def __init__(self, space1, space2, types1, types2, iterList, NTrials1,
            NTrials2, **kw):
        if (not np.all(space1.shape == space2.shape)):
            raise ValueError('space1.shape != space2.shape!')

        self.space1   = space1
        self.space2   = space2
        self.types1   = types1
        self.types2   = types2
        self.iterList = iterList
        self.NTrials1 = NTrials1
        self.NTrials2 = NTrials2
        self.kw       = kw

        self.xSize = self.space1.shape[1] # sp.shape is (r, c)
        self.ySize = self.space1.shape[0]

        self.color2D     = kw.pop('color2D', False)
        self.colormap2D  = kw.pop('colormap2d', None)
        self.resolveColors()


    def resolveColors(self):
        CX, CY = None, None
        if self.color2D:
            CX, CY = np.meshgrid(np.linspace(0, self.xSize-1, self.xSize),
                    np.linspace(0, self.ySize-1, self.ySize))
            if (self.colormap2D is None):
                xColor = [1, 0, 0]
                yColor = [0, 1, 0]
                cmap = Colormap2D((self.xSize, self.ySize), xColor, yColor)
            else:
                cmap = colormap2D
            colors = (CX + CY*self.xSize).flatten()
        else:
            cmap = None
            colors = None

        self.cmap = cmap
        self.colors = colors
        self.CX, self.CY = CX, CY

    def plot(self):
        kw = self.kw
        ax          = kw.pop('ax', plt.gca())
        ignoreNaNs  = kw.pop('ignoreNaNs', False)
        xlabel      = kw.pop('xlabel', '')
        ylabel      = kw.pop('ylabel', '')
        sigmaTitle  = kw.pop('sigmaTitle', False)
        noise_sigma = kw.pop('noise_sigma', None)

        X, _, _ = aggr.aggregateType(self.space1, self.iterList, self.types1,
                self.NTrials1, ignoreNaNs=ignoreNaNs, **kw)
        Y, _, _  = aggr.aggregateType(self.space2, self.iterList, self.types2,
                self.NTrials2, ignoreNaNs=ignoreNaNs, **kw)

        globalAxesSettings(ax)

                
        kw['edgecolor'] = 'white'
        if (self.colors is not None):
            kw['c'] = self.colors
        kw['cmap'] = self.cmap
        kw['picker'] = True
        kw['vmin'] = 0
        kw['vmax'] = self.xSize * self.ySize - 1

        collection = ax.scatter(X.flatten(), Y.flatten(), **kw)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.margins(0.05)
        ax.autoscale_view(tight=True)

        if sigmaTitle and noise_sigma is not None:
            ax.set_title('$\sigma$ = {0} pA'.format(int(noise_sigma)))

        # Picking - for interactive mode
        ax.figure.canvas.mpl_connect('pick_event', self.onpick)

        return ax

    def plotColorbar(self, ax):
        globalAxesSettings(ax)
        collection = ax.pcolor(self.CX + self.CY*self.xSize, cmap=self.cmap,
                rasterized=True)
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        ax.set_xlabel(xlabelTextShort, size='small')
        ax.set_ylabel(ylabelTextShort, size='small')
        ax.axis('scaled')
        return collection

    def onpick(self, event):
        idx = event.ind
        colors = event.artist.get_array()
        r = colors[idx] // self.xSize
        c = colors[idx] % self.xSize
        print idx
        print("r: {0}, c: {1}".format(r, c))

