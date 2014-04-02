#
#   scatter.py
#
#   Scatter plots for GridCells.
#
#       Copyright (C) 2013  Lukas Solanka <l.solanka@sms.ed.ac.uk>
#
from abc import ABCMeta, abstractmethod

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
from copy import deepcopy

from plotting.global_defs import globalAxesSettings
from plotting.colors      import Colormap2D

from . import xlabelTextShort, ylabelTextShort
from . import aggregate as aggr



class ScatterPlotBase(object):
    __metaclass = ABCMeta

    def __init__(self, spaceSizes, **kw):
        self.kw       = kw
        self.xSize, self.ySize = spaceSizes

        self.color2D     = self.kw.pop('color2D', False)
        self.colormap2D  = self.kw.pop('colormap2D', None)
        self.resolveColors()


    @staticmethod
    def determineSpaceSizes(space):
        xSize = space.shape[1] # sp.shape is (r, c)
        ySize = space.shape[0]
        return (xSize, ySize)


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
            cmap = self.kw.get('cmap', None)
            colors = self.kw.get('colors', None)

        self.cmap = cmap
        self.colors = colors
        self.CX, self.CY = CX, CY


    def setColors(self, c):
        self.colors = c

    def _plot(self, X, Y, **kw):
        ax          = kw.pop('ax', plt.gca())
        xlabel      = kw.pop('xlabel', '')
        ylabel      = kw.pop('ylabel', '')
        sigmaTitle  = kw.pop('sigmaTitle', False)
        noise_sigma = kw.pop('noise_sigma', None)

        globalAxesSettings(ax)
                
        kw['edgecolor'] = 'white'
        if (self.colors is not None):
            kw['c'] = self.colors
        kw['cmap'] = self.cmap
        kw['picker'] = True
        if (self.colormap2D):
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
        #ax.figure.canvas.mpl_connect('pick_event', self.onpick)

        self.ax = ax
        return ax


    @abstractmethod
    def plot(self):
        '''
        Plot the scatter plot with parameters defined in the keyword arguments
        in the constructor.
        '''
        raise NotImplementedError()


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



class RawScatterPlot(ScatterPlotBase):
    def __init__(self, templateSpace, X, Y, **kw):
        spaceSizes = ScatterPlotBase.determineSpaceSizes(templateSpace)
        ScatterPlotBase.__init__(self, spaceSizes, **kw)

        self.X = X
        self.Y = Y

    def plot(self):
        return self._plot(self.X, self.Y, **self.kw)


class ScatterPlot(RawScatterPlot):
    def __init__(self, space1, space2, types1, types2, iterList, NTrials1,
            NTrials2, **kw):

        self.space1   = space1
        self.space2   = space2
        self.types1   = types1
        self.types2   = types2
        self.iterList = iterList
        self.NTrials1 = NTrials1
        self.NTrials2 = NTrials2

        ignoreNaNs  = kw.pop('ignoreNaNs', False)
        if isinstance(space1, aggr.AggregateData):
            X, _, _ = space1.getData()
            Y, _, _ = space2.getData()
            RawScatterPlot.__init__(self, X, X, Y, **kw)
        else:
            if (not np.all(space1.shape == space2.shape)):
                raise ValueError('space1.shape != space2.shape!')
            X, _, _ = aggr.aggregateType(self.space1, self.iterList,
                                         self.types1, self.NTrials1,
                                         ignoreNaNs=ignoreNaNs, **kw)
            Y, _, _ = aggr.aggregateType(self.space2, self.iterList,
                                         self.types2, self.NTrials2,
                                         ignoreNaNs=ignoreNaNs, **kw)
            RawScatterPlot.__init__(self, self.space1, X, Y, **kw)



class FullScatterPlot(object):
    '''
    Create a figure with scatter plots for all noise levels, plus a 2D colorbar
    '''
    def __init__(self, spaces1, spaces2, types1, types2, iterList, NTrials1,
            NTrials2, **kw):
        self.noise_sigmas   = kw.pop('noise_sigmas', None)
        self.fig            = kw.pop('fig', plt.gcf())
        self.captionLetters = kw.pop('captionLetters', None)

        self.scatterPlots = []
        self.axes = []
        if (len(spaces1) != len(spaces2)):
            raise RuntimeError("len(spaces1) != len(spaces2)")
        for spIdx, (sp1, sp2) in enumerate(zip(spaces1, spaces2)):
            currentAx = self.fig.add_subplot(len(spaces1), 1, spIdx+1)
            self.axes.append(currentAx)
            currentKw = kw.copy()
            currentKw['ax'] = currentAx
            if (self.noise_sigmas is not None):
                currentKw['noise_sigma'] = self.noise_sigmas[spIdx]
            if (spIdx != len(spaces1) - 1):
                currentKw['xlabel'] = ''
           
            if isinstance(sp1, aggr.AggregateData):
                data1, _, _ = sp1.getData()
                data2, _, _ = sp2.getData()
                currentKw.pop('ignoreNaNs')
                scatterPlot = RawScatterPlot(data1, data1.flatten(),
                        data2.flatten(), **currentKw)
            else:
                scatterPlot = ScatterPlot(sp1, sp2, types1, types2, iterList,
                        NTrials1, NTrials2, **currentKw)
            self.scatterPlots.append(scatterPlot)

    def plot(self, plotcolorbar=True, captionLeft=-0.075):
        for idx, scatterPlot in enumerate(self.scatterPlots):
            scatterPlot.plot()
            if (self.captionLetters is not None):
                ax = self.axes[idx]
                ax.text(captionLeft, 1.0, self.captionLetters[idx], va='bottom',
                        ha='center', transform=ax.transAxes, size=14,
                        fontweight='bold')

        self.fig.tight_layout(h_pad=2.5)
        if (plotcolorbar):
            self.plotColorbar()


    def set_titleSizes(self, sz):
        for ax in self.axes:
            ax.set_title(ax.get_title(), size=sz)

    def plotColorbar(self, left=0.825, bottom=0.85, right=0.99, top=0.95):
        self.colorBarAx = self.fig.add_axes(Bbox.from_extents(left, bottom,
            right, top))
        self.scatterPlots[0].plotColorbar(self.colorBarAx)



class DiffScatterPlot(RawScatterPlot):
    '''
    Plot a scatter plot of differences between the spaces, i.e. noise levels.
    '''
    def __init__(self, spaces1, spaces2, types1, types2, iterList, NTrials1,
            NTrials2, which, **kw):
        '''
        **Parameters:**

        which : int
            Index into the difference array. Value of `0` means ``diffArray[1]
            - diffArray[0]``.
        '''
        if (len(spaces1) != len(spaces2)):
            raise RuntimeError("len(spaces1) != len(spaces2)")

        self.spaces1  = spaces1
        self.spaces2  = spaces2
        self.types1   = types1
        self.types2   = types2
        self.iterList = iterList
        self.NTrials1 = NTrials1
        self.NTrials2 = NTrials2
        self.which    = which

        ignoreNaNs  = kw.pop('ignoreNaNs', False)

        if isinstance(spaces1[0], aggr.AggregateData):
            data1, _, _ = aggr.collapseNoiseAggregated(spaces1)
            data2, _, _ = aggr.collapseNoiseAggregated(spaces2)
            templSpace, _, _ = spaces1[0].getData()
        else:
            data1, _, _ = aggr.collapseNoise(self.spaces1, iterList,
                    self.types1, self.NTrials1, ignoreNaNs=ignoreNaNs)
            data2, _, _ = aggr.collapseNoise(self.spaces2, iterList,
                    self.types2, self.NTrials2, ignoreNaNs=ignoreNaNs)
            templSpace = spaces1[0]

        self.diffData1 = np.diff(data1, axis=0)
        self.diffData2 = np.diff(data2, axis=0)

        X = self.diffData1[self.which, :]
        Y = self.diffData2[self.which, :]
        RawScatterPlot.__init__(self, templSpace, X, Y, **kw)

