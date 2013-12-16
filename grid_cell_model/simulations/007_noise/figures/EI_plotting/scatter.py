#
#   scatter.py
#
#   Scatter plots for GridCells.
#
#       Copyright (C) 2013  Lukas Solanka <l.solanka@sms.ed.ac.uk>
#
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
from copy import deepcopy

from plotting.global_defs import globalAxesSettings
from plotting.colors      import Colormap2D

from . import xlabelTextShort, ylabelTextShort
from . import aggregate as aggr



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

        self.color2D     = self.kw.pop('color2D', False)
        self.colormap2D  = self.kw.pop('colormap2d', None)
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

        self.ax = ax
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
           
            self.scatterPlots.append(ScatterPlot(
                sp1, sp2, types1, types2, iterList, NTrials1,
                NTrials2, **currentKw))

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
