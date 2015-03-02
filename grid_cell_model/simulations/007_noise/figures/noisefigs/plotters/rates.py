'''Firing rate plotters.'''
from __future__ import absolute_import, print_function

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.ticker as ti
from matplotlib.transforms import Bbox
from matplotlib.colors import LogNorm
from copy import deepcopy

from grid_cell_model.plotting.global_defs import globalAxesSettings
from simtools.plotting.plotters import FigurePlotter
from ..EI_plotting import sweeps, scatter
from ..EI_plotting import aggregate as aggr
from .base import SweepPlotter


__all__ = [
    'FRSweepPlotter',
    'ScatterGridsFRAllPlotter',
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

        example_idx = self.config['grids']['example_idx']

        FR_e_vmin = 0.12
        FR_e_vmax = 54.53
        FR_i_vmin = 7.1
        FR_i_vmax = 295

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            # E cells
            fname_e = self.get_fname("suppFigure_grids_FR_E_{ns}.pdf",
                                     ns=noise_sigma)
            data_e = aggr.AvgPopulationFR('FR_e/avg',
                                          ps.grids[ns_idx], iter_list,
                                          ignoreNaNs=True, normalizeTicks=True,
                                          r=example_idx[ns_idx][0],
                                          c=example_idx[ns_idx][1])
            data_e_nozero = aggr.NoZeroCouplingFilter(data_e)
            gridness_score = aggr.GridnessScore(
                ps.grids[ns_idx], iter_list,
                ignoreNaNs=True, normalizeTicks=True,
                r=example_idx[ns_idx][0],
                c=example_idx[ns_idx][1])


            with self.figure_and_axes(fname_e, sweepc) as (fig, ax):
                if (ns_idx == 2):
                    cbar = True
                else:
                    cbar = False
                kw = {}
                if (ns_idx != 0):
                    kw['ylabel'] = ''
                    kw['yticks'] = False

                sweeps.plotSweep(
                        data_e_nozero,
                        noise_sigma,
                        ax=ax,
                        #xlabel='', xticks=False,
                        cbar=cbar,
                        cbar_kw=self.myc['cbar_kw_e'],
                        #vmin=FR_e_vmin, vmax=FR_e_vmax,
                        norm=LogNorm(vmin=FR_e_vmin, vmax=FR_e_vmax),
                        **kw)

                if self.myc['plot_grid_contours'][ns_idx]:
                    contours = sweeps.Contours(gridness_score,
                            self.config['sweeps']['grid_contours'])
                    contours.plot(
                            ax,
                            **self.config['sweeps']['contours_kwargs'])

            gridness_threshold = .2
            grids_filter    = aggr.GTFilter(gridness_score, gridness_threshold)
            no_grids_filter = aggr.LEQFilter(gridness_score, gridness_threshold)
            FR_grids, _, _    = data_e_nozero.filter_data(grids_filter).getData()
            FR_no_grids, _, _ = data_e_nozero.filter_data(no_grids_filter).getData()
            print("mean E FR (grids):    {0}".format(np.mean(FR_grids)))
            print("mean E FR (no grids): {0}".format(np.mean(FR_no_grids)))

            # I cells
            fname_i = self.get_fname("suppFigure_grids_FR_I_{ns}.pdf",
                                     ns=noise_sigma)
            with self.figure_and_axes(fname_i, sweepc) as (fig, ax):
                data_i = aggr.AvgPopulationFR(
                        'FR_i/all',
                        ps.grids[ns_idx], iter_list,
                        ignoreNaNs=True, normalizeTicks=True,
                        r=example_idx[ns_idx][0],
                        c=example_idx[ns_idx][1])
                flt_data_i = aggr.NoZeroCouplingFilter(data_i)

                sweeps.plotSweep(
                        flt_data_i,
                        noise_sigma,
                        ax=ax,
                        cbar=cbar,
                        cbar_kw=self.myc['cbar_kw_i'],
                        norm=LogNorm(vmin=FR_i_vmin, vmax=FR_i_vmax),
                        **kw)


##############################################################################
# Scatter plot of gridness score vs. firing rates
class ScatterGridsFRAllPlotter(FigurePlotter):
    def __init__(self, *args, **kwargs):
        super(ScatterGridsFRAllPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ps = self.env.ps
        iter_list = self.config['iter_list']
        l, b, r, t = self.myc['bbox_rect']

        scatterColors = ['green', 'red', 'blue']
        scatterOrders = [2, 3, 1]

        # E cells
        fig = self._get_final_fig(self.myc['fig_size'])
        ax = fig.add_axes(Bbox.from_extents(l, b, r, t))
        ax.hold('on')
        xlabel='Mean E firing rate (Hz)'

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            grid_data = aggr.GridnessScore(ps.grids[ns_idx], iter_list,
                                               ignoreNaNs=True,
                                               normalizeTicks=False)
            fre_data = aggr.AvgPopulationFR('FR_e/avg',
                                            ps.grids[ns_idx], iter_list,
                                            ignoreNaNs=True,
                                            normalizeTicks=False)
            color = scatterColors[ns_idx]
            scatterPlot = scatter.ScatterPlot(
                    fre_data, grid_data,
                    None, None, None, None, None,
                    c=color,
                    s=self.myc['dot_size']*self.config['scale_factor'],
                    linewidth=0.3,
                    xlabel=xlabel,
                    ylabel=self.myc['ylabel'],
                    yticks=self.myc['yticks'],
                    zorder=scatterOrders[ns_idx])
            scatterPlot.plot()
        ax.xaxis.set_major_locator(ti.MultipleLocator(10))
        ax.yaxis.set_major_locator(ti.MultipleLocator(0.5))
        ax.yaxis.set_minor_locator(ti.MultipleLocator(0.25))

        leg = ['0', '150', '300']
        l = ax.legend(leg, **self.myc['legend_kwargs'])
        plt.setp(l.get_title(), size=self.myc['legend_kwargs']['fontsize'])

        fname = self.get_fname("/scatter_FRE_vs_grids.pdf")
        fig.savefig(fname, dpi=300, transparent=True)


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
