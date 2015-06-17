from __future__ import absolute_import, print_function

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ti
from matplotlib.gridspec   import GridSpec
from matplotlib.colorbar   import make_axes
from matplotlib.transforms import Bbox
from grid_cell_model.plotting.global_defs import prepareLims
from grid_cell_model.plotting.low_level   import zeroLines
from grid_cell_model.parameters           import JobTrialSpace2D
from grid_cell_model.parameters.metadata import (GEProfileWidthExtractor,
                                                 EISweepExtractor,
                                                 GenericExtractor,
                                                 Extractor1D)
from simtools.plotting.plotters import FigurePlotter

from ..EI_plotting import sweeps, details, examples, scatter
from ..EI_plotting import aggregate as aggr
from .base import SweepPlotter, SweepPlotter1D

__all__ = [
    'BumpSigmaSweepPlotter',
    'BumpDriftAtTimePlotter',
    'BumpDiffAtInitPlotter',
    'BumpDiffResetPlotter',
    'BumpExamplePlotter',
    'BumpExampleColorbarPlotter',
    'BumpSigmaDetailedNoisePlotter',
    'MainBumpFormationPlotter',
    'GEProfileWidthBumpPlotter',
    'Generic2DPBumpPlotter',
    'MainScatterGridsBumpsPlotter',
    'MainIsBumpPlotter',
    'Generic1DPBumpPlotter',
]

###############################################################################

bumpTStart  = 500.0
exampleRC   = ( (5, 15), (15, 5) )
exampleIdx  = [(0, 0), (0, 0), (0, 0)] # (row, col)

class BumpSigmaSweepPlotter(SweepPlotter):
    ''''Bump sigma sweeps.'''
    bump_vmin = 0
    bump_vmax = 0.421

    ann0 = dict(
            txt='b',
            rc=exampleRC[0],
            xytext_offset=(1.5, 0.5),
            color='white')
    ann1 = dict(
            txt='a',
            rc=exampleRC[1],
            xytext_offset=(1.2, 1.1),
            color='black')

    def __init__(self, *args, **kwargs):
        super(BumpSigmaSweepPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        myc= self._get_class_config()
        sweepc = self._get_sweep_config()
        ps = self.env.ps
        iter_list = self.config['iter_list']
        ann = [self.ann0, self.ann1]
        n_trials = self.config['bumps']['n_trials']

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            fname = (self.config['output_dir'] +
                     "/bumps_sweeps{0}.pdf".format(int(noise_sigma)))
            with self.figure_and_axes(fname, sweepc) as (fig, ax):
                kw = dict()
                if ns_idx != 0:
                    kw['ylabel'] = ''
                    kw['yticks'] = False
                data = aggr.AggregateBumpReciprocal(
                        ps.bumpGamma[ns_idx],
                        iter_list,
                        n_trials, tStart=bumpTStart)
                _, _, cax = sweeps.plotSweep(data,
                        noise_sigma=noise_sigma,
                        ax=ax,
                        cbar_kw=myc['cbar_kw'],
                        cbar=self.myc['cbar'][ns_idx],
                        vmin=self.bump_vmin, vmax=self.bump_vmax,
                        annotations=ann, **kw)


class BumpDriftAtTimePlotter(SweepPlotter):
    '''Bump drift at a specified time.'''
    def __init__(self, *args, **kwargs):
        super(BumpDriftAtTimePlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        myc= self._get_class_config()
        sweepc = self._get_sweep_config()
        ps = self.env.ps
        iter_list = self.config['iter_list']
        n_trials = self.config['bumps']['n_trials']
        grids_example_idx = self.config['grids']['example_idx']

        bumpDriftTStart = 1e3 #ms
        bumpDriftT = 9e3 # ms
        drift_vmin = 0
        drift_vmax = 20

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            fname = (self.config['output_dir'] +
                     "/bumps_drift_at_time_sweeps{0}.pdf".format(int(noise_sigma)))
            with self.figure_and_axes(fname, sweepc) as (fig, ax):
                kw = dict(cbar=False)
                if ns_idx != 0:
                    kw['ylabel'] = ''
                    kw['yticks'] = False
                if ns_idx == 2:
                    kw['cbar'] = True
                data = aggr.BumpDriftAtTime(bumpDriftT,
                        ps.bumpGamma[ns_idx],
                        iter_list,
                        n_trials,
                        tStart=bumpDriftTStart)
                _, _, cax = sweeps.plotSweep(
                    data, noise_sigma=noise_sigma,
                    ax=ax,
                    cbar_kw=myc['cbar_kw'],
                    vmin=drift_vmin, vmax=drift_vmax,
                    axis_setting=self.myc.get('axis_setting', 'scaled'),
                    **kw)

                if self.myc['plot_grid_contours'][ns_idx]:
                    gridData = aggr.GridnessScore(
                        ps.grids[ns_idx], iter_list,
                        ignoreNaNs=True,
                        normalizeTicks=True,
                        r=grids_example_idx[ns_idx][0],
                        c=grids_example_idx[ns_idx][1])
                    contours = sweeps.Contours(gridData,
                            self.config['sweeps']['grid_contours'])
                    contours.plot(
                            ax,
                            **self.config['sweeps']['contours_kwargs'])


class BumpDiffAtInitPlotter(SweepPlotter):
    '''Distance from init position.'''
    def __init__(self, *args, **kwargs):
        super(BumpDiffAtInitPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        myc= self._get_class_config()
        sweepc = self._get_sweep_config()
        ps = self.env.ps
        iter_list = self.config['iter_list']
        n_trials = self.config['bumps']['n_trials']

        bumpDiffT = 0.75e3 # ms
        bumpDiff_vmin = 0
        bumpDiff_vmax = 20
        diffStartPos = [17, 15]

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            fname = (self.config['output_dir'] +
                     "/bumps_difference_at_time_sweeps{0}.pdf".format(int(noise_sigma)))
            with self.figure_and_axes(fname, sweepc) as (fig, ax):
                kw = dict(cbar=False)
                if ns_idx != 0:
                    kw['ylabel'] = ''
                    kw['yticks'] = False
                if ns_idx == 2:
                    kw['cbar'] = True
                data = aggr.BumpDifferenceAtTime(diffStartPos, bumpDiffT,
                        ps.bumpGamma[ns_idx],
                        iter_list,
                        n_trials)
                _, _, cax = sweeps.plotSweep(data, noise_sigma=noise_sigma,
                        ax=ax,
                        cbar_kw=myc['cbar_kw'],
                        vmin=bumpDiff_vmin, vmax=bumpDiff_vmax,
                        **kw)


class BumpDiffResetPlotter(SweepPlotter):
    '''Average distance from position enforced by place cells.

    From theta_start_t until the end of the simultion.
    '''
    def __init__(self, *args, **kwargs):
        super(BumpDiffResetPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        myc= self._get_class_config()
        sweepc = self._get_sweep_config()
        ps = self.env.ps
        iter_list = self.config['iter_list']

        bumpResetTStart = 0.5e3 # ms
        bumpReset_vmin = 0
        bumpReset_vmax = 20
        bumpResetStartPos = [17, 15]
        constPosNTrials = 5

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            fname = (self.config['output_dir'] +
                     "/bumps_avg_difference_reset_sweeps{0}.pdf".format(int(noise_sigma)))
            with self.figure_and_axes(fname, sweepc) as (fig, ax):
                kw = dict(cbar=False)
                if ns_idx != 0:
                    kw['ylabel'] = ''
                    kw['yticks'] = False
                if ns_idx == 2:
                    kw['cbar'] = True
                data = aggr.BumpAvgDifferenceFromPos(bumpResetStartPos,
                        ps.constPos[ns_idx],
                        iter_list,
                        constPosNTrials,
                        tstart=bumpResetTStart)
                _, _, cax = sweeps.plotSweep(data, noise_sigma=noise_sigma,
                        ax=ax,
                        cbar_kw=myc['cbar_kw'],
                        vmin=bumpReset_vmin, vmax=bumpReset_vmax,
                        **kw)
                self.plot_grid_contours(ns_idx, ax, ps.grids)


class BumpExamplePlotter(FigurePlotter):
    '''Bump examples.'''
    def __init__(self, *args, **kwargs):
        super(BumpExamplePlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ps = self.env.ps
        output_dir = self.config['output_dir']
        iter_list = self.config['iter_list']

        exampleEFName = output_dir + "/bumps_examples_E_{0}pA_{1}.pdf"
        exampleIFName = output_dir + "/bumps_examples_I_{0}pA_{1}.pdf"
        bumpExampleTypes = ['bump_full']
        bumpTrialNum = 0
        exTransparent = True
        exampleFigSize = (0.8, 0.8)
        l, b, r, t = self.myc['bbox']

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            for idx, rc in enumerate(exampleRC):
                for EIType in ['E', 'I']:
                    if EIType == 'E':
                        fnameTemplate =exampleEFName
                        types = bumpExampleTypes + ['rateMap_e']
                    else:
                        fnameTemplate =exampleIFName
                        types = bumpExampleTypes + ['rateMap_i']
                    fname = fnameTemplate.format(noise_sigma, idx)
                    fig = self._get_final_fig(exampleFigSize)
                    gs = examples.plotOneBumpExample(ps.bumpGamma[ns_idx], rc, iter_list,
                            types,
                            exIdx=exampleIdx[ns_idx],
                            trialNum=bumpTrialNum,
                            fig=fig)
                    gs.update(left=l, bottom=b, right=r, top=t)
                    fig.savefig(fname, dpi=300, transparent=exTransparent)
                    plt.close()


class BumpExampleColorbarPlotter(FigurePlotter):
    def __init__(self, *args, **kwargs):
        super(BumpExampleColorbarPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        fig = self._get_final_fig(self.myc['fig_size'])
        ax_cbar = fig.add_axes([0.05, 0.07, 0.8, 0.1])
        cbar = mpl.colorbar.ColorbarBase(ax_cbar, cmap=mpl.cm.jet,
                                         norm=mpl.colors.Normalize(vmin=0,
                                                                   vmax=1),
                                         ticks=[],
                                         orientation='horizontal')
        #ax_cbar.yaxis.set_ticklabels(['0', 'Max'])
        fname_cbar = self.get_fname("bumps_examples_colorbar.pdf")
        plt.savefig(fname_cbar, dpi=300, transparent=True)
        plt.close()


class BumpSigmaDetailedNoisePlotter(FigurePlotter):
    '''Detailed noise plots.'''
    def __init__(self, *args, **kwargs):
        super(BumpSigmaDetailedNoisePlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        iter_list = self.config['iter_list']
        sigma_bump_text = self.config['bump_sigma']['sigma_bump_text']

        EI13Root  = 'simulation_data/main_network/detailed_noise/gamma_bump/EI-1_3'
        EI31Root  = 'simulation_data/main_network/detailed_noise/gamma_bump/EI-3_1'
        detailedShape = (31, 9)

        EI13PS = JobTrialSpace2D(detailedShape, EI13Root)
        EI31PS = JobTrialSpace2D(detailedShape, EI31Root)
        detailedNTrials = 5

        ylabelPos = -0.17
        detailFigSize = (3.8, 2.6)
        detailLeft   = 0.18
        detailBottom = 0.26
        detailRight  = 0.95
        detailTop    = 0.95

        iter_list = ['noise_sigma', 'g_AMPA_total']

        types = ('bump', 'sigma')
        fig = plt.figure(figsize=detailFigSize)
        ax = fig.add_axes(Bbox.from_extents(detailLeft, detailBottom, detailRight,
            detailTop))

        data13 = aggr.AggregateBumpReciprocal(
                EI13PS,
                iter_list,
                detailedNTrials, tStart=bumpTStart,
                normalizeTicks=False)
        _, p13, l13 = details.plotDetailedNoise(data13, None, None, ax=ax,
                ylabel=sigma_bump_text, ylabelPos=ylabelPos,
                color='red', markerfacecolor='red', zorder=10)

        data31 = aggr.AggregateBumpReciprocal(
                EI31PS,
                iter_list,
                detailedNTrials, tStart=bumpTStart,
                normalizeTicks=False)
        _, p31, l31 = details.plotDetailedNoise(data31, None, None, ax=ax,
                ylabelPos=ylabelPos,
                color='#505050')

        #ax.set_yscale("log")
        #ax.set_ylim([1.5, 300])
        leg = ['a',  'b']
        l = ax.legend([p31, p13], leg, loc=(0.85, 0.1), fontsize='small', frameon=False,
                numpoints=1, handletextpad=0.05)
        plt.setp(l.get_title(), fontsize='small')

        fname = self.config['output_dir'] + "/bumps_detailed_noise_sigma.pdf"
        plt.savefig(fname, dpi=300, transparent=True)
        plt.close()


class MainScatterGridsBumpsPlotter(FigurePlotter):
    '''Correlate P(bump) vs gridness score.'''
    def __init__(self, *args, **kwargs):
        super(MainScatterGridsBumpsPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ps = self.env.ps
        myc = self._get_class_config()

        iter_list = self.config['iter_list']
        l, b, r, t = self.myc['bbox_rect']
        legend = self.myc.get('legend', True)
        legend_kwargs = myc['legend_kwargs']

        xlabel = self.myc.get('xlabel', 'P(bumps)')
        ylabel = 'Gridness score'

        fig = self._get_final_fig(myc['fig_size'])
        ax = fig.add_axes(Bbox.from_extents(l, b, r, t))

        scatterColors = ['green', 'red', 'blue']
        scatterOrders = [2, 3, 1]

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            isBumpData = aggr.IsBump(ps.bumpGamma[ns_idx], iter_list,
                                     ignoreNaNs=True)
            gridData = aggr.GridnessScore(ps.grids[ns_idx], iter_list,
                                          ignoreNaNs=True,
                                          normalizeTicks=False)

            scatterPlot = scatter.ScatterPlot(
                    isBumpData, gridData, None, None, None, None, None,
                    c=scatterColors[ns_idx],
                    s=6*self.config['scale_factor'],
                    linewidth=0.3,
                    xlabel=xlabel,
                    ylabel=ylabel,
                    ax=ax,
                    zorder=scatterOrders[ns_idx])
            scatterPlot.plot()

        ax.xaxis.set_major_locator(ti.MultipleLocator(0.2))
        ax.yaxis.set_major_locator(ti.MultipleLocator(0.5))
        if legend:
            leg = ['0', '150', '300']
            l = ax.legend(leg, **legend_kwargs)
            plt.setp(l.get_title(), size='small')
        #ax.set_ylabel(ax.get_ylabel(), y=0., ha='left')

        # Normal scale
        fname = self.config['output_dir'] + "/bumps_scatter_grids_vs_bumpFracTotal.pdf"
        fig.savefig(fname, dpi=300, transparent=True)

        # Exponential scale
        ax.set_xscale('exponential')
        ax.xaxis.set_major_locator(ti.MultipleLocator(.5))
        ax.xaxis.set_minor_locator(ti.MultipleLocator(.1))
        ax.set_xlim([-0.3, 1.002])
        fname = self.config['output_dir'] + "/bumps_scatter_grids_vs_bumpFracTotal_exp.pdf"
        fig.savefig(fname, dpi=300, transparent=True)


class BumpFormationBase(SweepPlotter):
    '''Bump formation sweeps.'''
    bump_vmin = 0
    bump_vmax = 1.

    ann0_0 = dict(
            txt='a',
            rc=(15, 5),
            xytext_offset=(.5, 2),
            color='white')
    ann0_1 = dict(
            txt='b',
            rc=(5, 15),
            xytext_offset=(1.5, 1),
            color='white')

    ann150_0 = dict(
            txt='c',
            rc=(5, 15),
            xytext_offset=(1.5, 1),
            color='white')

    ann300_0 = dict(
            txt='d',
            rc=(15, 5),
            xytext_offset=(1.5, 1),
            color='white')
    ann300_1 = dict(
            txt='e',
            rc=(5, 15),
            xytext_offset=(1.5, 1),
            color='white')

    def __init__(self, *args, **kwargs):
        super(BumpFormationBase, self).__init__(*args, **kwargs)

    def get_ann(self):
        ann_config = self.myc.get('ann', None)
        if ann_config is not None:
            return ann_config

        ann0   = [  self.ann0_0,   self.ann0_1]
        ann150 = [self.ann150_0]
        ann300 = [self.ann300_0, self.ann300_1]
        return [ann0, ann150, ann300]


class MainBumpFormationPlotter(BumpFormationBase):
    '''This one goes to the main figure'''
    def __init__(self, *args, **kwargs):
        super(MainBumpFormationPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        myc= self._get_class_config()
        sweepc = self._get_sweep_config()
        ps = self.env.ps
        iter_list = self.config['iter_list']
        xlabel = self.myc.get('xlabel', None)
        xticks = myc['xticks']
        vmin = self.myc.get('vmin', self.bump_vmin)
        vmax = self.myc.get('vmax', self.bump_vmax)

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            fname = self.get_fname(
                "bumps_mainFig_isBumpFracTotal_sweeps_annotated{ns}.pdf",
                ns=noise_sigma)
            with self.figure_and_axes(fname, sweepc) as (_, ax):
                #fig, ax = ds.getDefaultSweepFig(scale=0.8, colorBarPos='left')
                kw = dict()
                if ns_idx != 0:
                    kw['ylabel'] = ''
                    kw['yticks'] = False
                metadata = EISweepExtractor(ps.bumpGamma[ns_idx])
                data = aggr.IsBump(
                    ps.bumpGamma[ns_idx],
                    iter_list,
                    ignoreNaNs=True,
                    metadata_extractor=metadata)
                sweeps.plotSweep(
                    data,
                    noise_sigma=noise_sigma,
                    sigmaTitle=self.myc.get('sigmaTitle', True),
                    xlabel='' if xticks[ns_idx] == False else xlabel,
                    xticks=xticks[ns_idx],
                    ax=ax,
                    cbar=self.myc['cbar'][ns_idx],
                    cbar_kw=myc['cbar_kw'],
                    vmin=vmin, vmax=vmax,
                    annotations=self.get_ann()[ns_idx],
                    axis_setting=self.myc.get('axis_setting', 'scaled'),
                    **kw)

                self.plot_grid_contours(ns_idx, ax, ps.grids)


class GEProfileWidthBumpPlotter(BumpFormationBase):
    '''Bump formation plotter for gE vs. pAMPA_sigma.'''
    def __init__(self, *args, **kwargs):
        super(GEProfileWidthBumpPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        myc= self._get_class_config()
        sweepc = self._get_sweep_config()
        ps = self.env.ps
        iter_list = self.config['iter_list']
        xlabel = self.myc.get('xlabel', None)
        xticks = myc['xticks']
        l, b, r, t = self.myc['bbox']
        fname = self.myc.get('fname', "bumps_Pbumps_gE_pAMPA_sigma{ns}.pdf")

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            fname = self.get_fname(fname, ns=noise_sigma)
            fig = self._get_final_fig(self.config['sweeps']['fig_size'])
            ax = fig.add_axes(Bbox.from_extents(l, b, r, t))
            kw = dict()
            if ns_idx != 0:
                kw['ylabel'] = ''
                kw['yticks'] = False
            metadata = GEProfileWidthExtractor(ps.bumpGamma[ns_idx])
            data = aggr.IsBump(
                ps.bumpGamma[ns_idx],
                iter_list,
                ignoreNaNs=True,
                metadata_extractor=metadata)
            sweeps.plotSweep(
                data,
                noise_sigma=noise_sigma,
                sigmaTitle=self.myc.get('sigmaTitle', True),
                xlabel='' if xticks[ns_idx] == False else xlabel,
                xticks=xticks[ns_idx],
                ax=ax,
                cbar=self.myc['cbar'][ns_idx],
                cbar_kw=myc['cbar_kw'],
                vmin=self.bump_vmin, vmax=self.bump_vmax,
                annotations=self.get_ann()[ns_idx],
                **kw)
            ax.axis('tight')
            fig.savefig(fname, dpi=300, transparent=True)


class Generic2DPBumpPlotter(BumpFormationBase):
    '''Bump formation plotter for a generic set of two variables.'''
    def __init__(self, *args, **kwargs):
        super(Generic2DPBumpPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        myc= self._get_class_config()
        sweepc = self._get_sweep_config()
        ps = self.env.ps
        iter_list = self.config['iter_list']
        xlabel = self.myc.get('xlabel', None)
        ylabel = self.myc.get('ylabel', None)
        xticks = myc['xticks']
        yticks = myc['yticks']
        vmin = self.myc.get('vmin', self.bump_vmin)
        vmax = self.myc.get('vmax', self.bump_vmax)
        normalize_type = self.myc.get('normalize_type', (None, None))
        l, b, r, t = self.myc['bbox']
        fname = self.myc.get('fname', "bumps_pbumps_generic_{ns}.pdf")

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            file_name = self.get_fname(fname, ns=noise_sigma)
            fig = self._get_final_fig(self.config['sweeps']['fig_size'])
            ax = fig.add_axes(Bbox.from_extents(l, b, r, t))
            metadata = GenericExtractor(ps.bumpGamma[ns_idx],
                                        normalize=self.myc['normalize_ticks'],
                                        normalize_type=normalize_type)
            data = aggr.IsBump(ps.bumpGamma[ns_idx],
                               iter_list,
                               ignoreNaNs=True,
                               metadata_extractor=metadata)
            sweeps.plotSweep(
                data,
                noise_sigma=noise_sigma,
                sigmaTitle=self.myc.get('sigmaTitle', True),
                xlabel='' if xticks[ns_idx] == False else xlabel,
                xticks=xticks[ns_idx],
                ylabel='' if yticks[ns_idx] == False else ylabel,
                yticks=yticks[ns_idx],
                ax=ax,
                cbar=self.myc['cbar'][ns_idx],
                cbar_kw=myc['cbar_kw'],
                vmin=vmin, vmax=vmax,
                annotations=self.get_ann()[ns_idx],
                axis_setting=self.myc.get('axis_setting', 'scaled'))
            ax.axis('tight')
            fig.savefig(file_name, dpi=300, transparent=True)
            plt.close(fig)


class MainIsBumpPlotter(BumpFormationBase):
    '''Bump formation, thresholded.'''
    bumpThreshold = 0.95

    def __init__(self, *args, **kwargs):
        super(MainIsBumpPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        myc= self._get_class_config()
        sweepc = self._get_sweep_config()
        ps = self.env.ps
        iter_list = self.config['iter_list']

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            fname = (self.config['output_dir'] +
                     "/{0}bumps_mainFig_isBump_sweeps{1}.pdf".format(
                         self.config.get('fname_prefix', ''),
                         int(noise_sigma)))
            with self.figure_and_axes(fname, sweepc) as (fig, ax):
                kw = dict(cbar=False, sigmaTitle=False)
                if ns_idx != 0:
                    kw['ylabel'] = ''
                    kw['yticks'] = False
                data = aggr.BumpFormationFilter(self.bumpThreshold, ps.bumpGamma[ns_idx],
                        iter_list, ignoreNaNs=True)
                _, _, cax = sweeps.plotSweep(data,
                        noise_sigma=noise_sigma,
                        ax=ax,
                        cbar_kw=myc['cbar_kw'],
                        vmin=self.bump_vmin, vmax=self.bump_vmax,
                        annotations=self.get_ann()[ns_idx],
                        **kw)


class Generic1DPBumpPlotter(SweepPlotter1D):
    '''A 1D sweep plotter.

    Plots a 1D representation of data into a figure.
    '''
    def __init__(self, *args, **kwargs):
        super(Generic1DPBumpPlotter, self).__init__(*args, **kwargs)

    def get_data(self, ns_idx):
        normalize_ticks = self.myc.get('normalize_ticks', False)
        normalize_type = self.myc.get('normalize_type', None)

        metadata = Extractor1D(self.env.ps.bumpGamma[ns_idx],
                               normalize=normalize_ticks,
                               normalize_type=normalize_type)
        return aggr.IsBump(self.env.ps.bumpGamma[ns_idx],
                           None,
                           ignoreNaNs=True,
                           metadata_extractor=metadata)
