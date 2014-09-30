'''Firing rate plotters.'''
from __future__ import absolute_import, print_function

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.ticker as ti
from matplotlib.transforms import Bbox
from copy import deepcopy

from grid_cell_model.plotting.global_defs import globalAxesSettings
from ..EI_plotting import sweeps, scatter
from ..EI_plotting import aggregate as aggr
from .base import SweepPlotter


__all__ = [
    'FRSweepPlotter',
]


##############################################################################
def plotThresholdHist(var1, threshold1, var2, **kw):
    ax = kw.pop('ax', plt.gca())
    xlabel      = kw.pop('xlabel', '')
    ylabel      = kw.pop('ylabel', 'p($\cdot$)')
    sigmaTitle  = kw.pop('sigmaTitle', False)
    noise_sigma = kw.pop('noise_sigma', None)

    grp1Idx = var1 < threshold1
    grp2Idx = np.logical_not(grp1Idx)

    globalAxesSettings(ax)
    ax.hold('on')
    #import pdb; pdb.set_trace()
    grp1 = var2[grp1Idx]
    grp2 = var2[grp2Idx]
    h1 = ax.hist(grp1, **kw)
    h2 = ax.hist(grp2, **kw)

    print("mean(grp1): {0}".format(np.mean(grp1)))
    print("mean(grp2): {0}".format(np.mean(grp2)))

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if sigmaTitle and noise_sigma is not None:
        ax.set_title('$\sigma$ = {0} pA'.format(int(noise_sigma)))
    
    return h1, h2


def plotFRGridThresholded(dataSpace, threshold, FRTypes, iterList, NTrials, **kw):
    ignoreNaNs  = kw.pop('ignoreNaNs', False)

    typesGrids = ['grids', 'gridnessScore']
     
    GS, _, _ = aggr.aggregateType(dataSpace, iterList, typesGrids, NTrials,
            ignoreNaNs=ignoreNaNs, **kw)
    FR, _, _  = aggr.aggregateType(dataSpace, iterList, FRTypes, NTrials,
            ignoreNaNs=ignoreNaNs, **kw)

    return plotThresholdHist(GS.flatten(), threshold, FR.flatten(), **kw)

    

##############################################################################
# Parameter sweeps of E and I firing rates
class FRSweepPlotter(SweepPlotter):
    def __init__(self, *args, **kwargs):
        super(FRSweepPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        myc= self._get_class_config()
        sweepc = self._get_sweep_config()
        ps = self.env.ps
        iter_list = self.config['iter_list']

        gridTrialNumList = np.arange(self.config['grids']['ntrials'])
        varList_e = ['FR_e',  'avg']
        varList_i = ['FR_i',  'all']
        exampleIdx = [(1, 22), (1, 22), (1, 22)] # (row, col)

        FR_e_vmin = 0
        FR_e_vmax = 5
        FR_i_vmin = 0
        FR_i_vmax = 100

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            # E cells
            fname_e = self.get_fname("suppFigure_grids_FR_E_{ns}.pdf",
                                     ns=noise_sigma)
            with self.figure_and_axes(fname_e, sweepc) as (fig, ax):
                if (ns_idx == 2):
                    cbar = True
                else:
                    cbar = False
                kw = {}
                if (ns_idx != 0):
                    kw['ylabel'] = ''
                    kw['yticks'] = False

                sweeps.plotFRTrial(ps.grids[ns_idx], varList_e, iter_list,
                        noise_sigma=noise_sigma,
                        trialNumList=gridTrialNumList,
                        ax=ax,
                        xlabel='', xticks=False,
                        r=exampleIdx[ns_idx][0], c=exampleIdx[ns_idx][1],
                        cbar=cbar,
                        cbar_kw=self.myc['cbar_kw_e'],
                        vmin=FR_e_vmin, vmax=FR_e_vmax,
                        ignoreNaNs=True,
                        **kw)

            # I cells
            fname_i = self.get_fname("suppFigure_grids_FR_I_{ns}.pdf",
                                     ns=noise_sigma)
            with self.figure_and_axes(fname_i, sweepc) as (fig, ax):
                sweeps.plotFRTrial(ps.grids[ns_idx], varList_i, iter_list,
                        noise_sigma=noise_sigma,
                        trialNumList=gridTrialNumList,
                        ax=ax,
                        r=exampleIdx[ns_idx][0], c=exampleIdx[ns_idx][1],
                        cbar=cbar,
                        cbar_kw=self.myc['cbar_kw_i'],
                        vmin=FR_i_vmin, vmax=FR_i_vmax,
                        sigmaTitle=False,
                        ignoreNaNs=True,
                        **kw)


###############################################################################
## Scatter plot of gridness score vs. gamma power
#scatterFigSize = (6.5, 3.6)
#scatterLeft   = 0.15
#scatterBottom = 0.2
#scatterRight  = 0.98
#scatterTop    = 0.95
#scatterTransparent = True
#
#ignoreNaNs = True
#
#if args.scatterPlot or args.all:
#    for EIType in ['E', 'I_10']:
#        if EIType == 'E':
#            fname = outputDir + "/suppFigure_grids_FR-scatter_FRE_vs_grids.pdf"
#            xlabel='Mean firing rate of E cells (Hz)'
#            legLoc = (0.8, 0.6)
#        else:
#            fname = outputDir + "/suppFigure_grids_FR-scatter_FRI_vs_grids.pdf"
#            xlabel='Mean firing rate of I cells (Hz)'
#            legLoc = (0.1, 0.6)
#
#        fig = plt.figure(figsize=scatterFigSize)
#        ax = fig.add_axes(Bbox.from_extents(scatterLeft, scatterBottom, scatterRight,
#            scatterTop))
#
#        NTrialsGrids = 3
#        typesFR = ['FR', EIType]
#        typesGrids = ['grids', 'gridnessScore']
#        ax.hold('on')
#        scatterColors = ['blue', 'green', 'red']
#
#        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
#            color = scatterColors[ns_idx]
#            scatterPlot = scatter.ScatterPlot(
#                    ps.grids[ns_idx], ps.grids[ns_idx], typesFR,
#                    typesGrids, iterList, NTrialsGrids, NTrialsGrids,
#                    c=color,
#                    s=15,
#                    linewidth=0.3,
#                    xlabel=xlabel,
#                    ylabel='Gridness score',
#                    ignoreNaNs=ignoreNaNs)
#            scatterPlot.plot()
#        ax.xaxis.set_major_locator(ti.MultipleLocator(10))
#        ax.yaxis.set_major_locator(ti.MultipleLocator(0.5))
#        ax.yaxis.set_minor_locator(ti.MultipleLocator(0.25))
#        ax.set_xscale('log')
#        if EIType == 'E':
#            ax.set_xlim([0.1, None])
#        else:
#            ax.set_xlim([2, 300])
#        leg = ['0', '150', '300']
#        l = ax.legend(leg, loc=legLoc, fontsize='small', frameon=False,
#                numpoints=1, title='$\sigma$ (pA)')
#        plt.setp(l.get_title(), size='small')
#
#        fig.savefig(fname, dpi=300, transparent=scatterTransparent)
#
#
###############################################################################
#histFigSize     = (3.75, 2.2)
#histLeft    = 0.25
#histBottom  = 0.3
#histRight   = 0.95
#histTop     = 0.8
#histTransparent = True
#
#if args.FRHistograms or args.all:
#    for EIType in ['E', 'I_10']:
#        if EIType == 'E':
#            fname = outputDir + "/suppFigure_grids_FR-histogram_FRE{0}.pdf"
#            xlabel='Firing rate of E cells (Hz)'
#            legLoc = (0.8, 0.6)
#            dataRange = (0, 10)
#            sigmaTitle = True
#        else:
#            fname = outputDir + "/suppFigure_grids_FR-histogram_FRI{0}.pdf"
#            xlabel='Firing rate of I cells (Hz)'
#            legLoc = (0.1, 0.6)
#            dataRange = (0, 200)
#            sigmaTitle = False
#
#        NTrialsGrids = 3
#        typesFR = ['FR', EIType]
#        threshold = 0.20 # Gridness score threshold
#        kw = {}
#        kw['range'] = dataRange
#        kw['bins'] = 20
#
#        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
#            fig = plt.figure(figsize=histFigSize)
#            ax = fig.add_axes(Bbox.from_extents(histLeft, histBottom, histRight,
#                histTop))
#
#            plotFRGridThresholded(
#                    ps.grids[ns_idx], threshold, typesFR,
#                    iterList, NTrialsGrids,
#                    xlabel=xlabel,
#                    ignoreNaNs=ignoreNaNs, 
#                    sigmaTitle=sigmaTitle,
#                    noise_sigma=noise_sigma,
#                    alpha=0.5,
#                    normed=True,
#                    **kw)
#            #ax.xaxis.set_major_locator(ti.MultipleLocator(10))
#            #ax.yaxis.set_major_locator(ti.MultipleLocator(0.5))
#            #ax.yaxis.set_minor_locator(ti.MultipleLocator(0.25))
#
#            if (ns_idx == 0 and EIType == 'E'):
#                leg = ['$gridness < {0}$'.format(threshold),
#                        '$gridness \geq {0}$'.format(threshold)]
#                l = ax.legend(leg, loc=(0.5, 0.5), fontsize='small',
#                        frameon=False)
#
#            fig.savefig(fname.format(int(noise_sigma)), dpi=300,
#                    transparent=histTransparent)
#            plt.close()
#
