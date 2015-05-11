from __future__ import absolute_import, print_function

import string

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ti
from matplotlib.transforms import Bbox

from grid_cell_model.parameters import JobTrialSpace2D
from grid_cell_model.data_storage import DataStorage
from grid_cell_model.data_storage.sim_models.ei import extractSummedSignals
from grid_cell_model.plotting.grids import (plotSpikes2D,
                                            plotGridRateMap,
                                            plotAutoCorrelation)
import grid_cell_model.plotting.low_level as low_level
from simtools.plotting.plotters import Computation, FigurePlotter

from ..EI_plotting import sweeps, examples, details, scatter
from ..EI_plotting import aggregate as aggr
from ..EI_plotting.base import getOption, plotStateSignal
from ..EI_plotting import scaling
from .base import SweepPlotter, ProbabilityPlotter, DummyPlotter

__all__ = [
    'GridSweepsPlotter',
    'GridExamplesPlotter',
    'GridExampleRectPlotter',
    'GridExampleColorbarPlotter',
    'VmExamplesPlotter',
    'GridDetailedNoisePlotter',
    'GridsDiffSweep',
    'GridBumpScatterPlotter',
    'GridsPBumpsProbabilityPlotter',
    'GridSimpleExamplePlotter',
    'HighGridScoreFraction',
]


class GridSweepsPlotter(SweepPlotter):
    '''Parameter sweeps of gridness score.'''
    cmap = 'jet'
    varList = ['gridnessScore']

    def __init__(self, *args, **kwargs):
        super(GridSweepsPlotter, self).__init__(*args, **kwargs)
        self.fig = None
        self.ax = None

    def plot(self, *args, **kwargs):
        sweepc = self._get_sweep_config()
        ps = self.env.ps
        example_idx = self.config['grids']['example_idx']
        iter_list = self.config['iter_list']

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            fname = (self.config['output_dir'] +
                     "/grids_sweeps{0}.pdf".format(int(noise_sigma)))
            self.data = aggr.GridnessScore(ps.grids[ns_idx], iter_list,
                                           ignoreNaNs=True, normalizeTicks=True,
                                           r=example_idx[ns_idx][0],
                                           c=example_idx[ns_idx][1])
            with self.figure_and_axes(fname, sweepc) as (self.fig, self.ax):
                # Sweep itself
                kw = dict()
                sweeps.plotGridTrial(
                    self.data,
                    None,
                    None,
                    ps.noise_sigmas[ns_idx],
                    trialNumList=None,
                    r=None, c=None,
                    ax=self.ax,
                    cbar=self.myc['cbar'][ns_idx],
                    cbar_kw=self.myc['cbar_kw'],
                    cmap=self.cmap,
                    vmin=self.myc['vmin'], vmax=self.myc['vmax'],
                    xlabel=self.myc['xlabel'][ns_idx],
                    xticks=self.myc['xticks'][ns_idx],
                    ylabel=self.myc['ylabel'][ns_idx],
                    yticks=self.myc['yticks'][ns_idx],
                    ignoreNaNs=True,
                    annotations=self.myc['ann'],
                    sliceAnn=None,
                    sigmaTitle=self.myc['sigma_title'],
                    **kw)

                # Contours
                if self.myc['plot_contours'][ns_idx]:
                    contours = sweeps.Contours(self.data,
                            self.config['sweeps']['grid_contours'])
                    contours.plot(
                            self.ax,
                            **self.config['sweeps']['contours_kwargs'])


class ContourGridSweepsPlotter(GridSweepsPlotter):
    '''Parameter sweeps of gridness score with a contour plot.'''
    def __init__(self, *args, **kwargs):
        super(ContourGridSweepsPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        super(ContourGridSweepsPlotter, self).plot(
                self, *args, plotContours=True, **kwargs)


class GridExamplesPlotter(FigurePlotter):
    '''Grid field examples for the main figure.'''
    exampleGridFName = "/grids_examples_{0}pA_{1}.pdf"
    exampleACFName = "/grids_examples_acorr_{0}pA_{1}.pdf"

    def __init__(self, *args, **kwargs):
        super(GridExamplesPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        myc= self._get_class_config()
        ps = self.env.ps
        exampleLeft, exampleBottom, exampleRight, exampleTop = myc['ax_box']
        example_idx = self.config['grids']['example_idx']

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            for idx, rc in enumerate(self.config['grids']['example_rc']):
                # Grid field
                fname = (self.config['output_dir'] +
                            self.exampleGridFName.format(noise_sigma, idx))
                fig = self._get_final_fig(myc['fig_size'])
                gs = examples.plotOneGridExample(
                        ps.grids[ns_idx],
                        rc,
                        self.config['iter_list'],
                        exIdx=example_idx[idx],
                        xlabel=False, ylabel=False,
                        xlabel2=False, ylabel2=False,
                        maxRate=True, plotGScore=False,
                        fig=fig)
                gs.update(left=exampleLeft, bottom=exampleBottom, right=exampleRight,
                        top=exampleTop)
                plt.savefig(fname, dpi=300, transparent=myc['transparent'])
                plt.close()

                # Autocorrelation
                fname = (self.config['output_dir'] +
                            self.exampleACFName.format(noise_sigma, idx))
                fig= self._get_final_fig(myc['fig_size'])
                ax = fig.add_axes(Bbox.from_extents(exampleLeft, exampleBottom,
                    exampleRight, exampleTop))
                gs = examples.plotOneGridACorrExample(ps.grids[ns_idx], rc,
                                                      ax=ax, vmin=0, vmax=None)
                plt.savefig(fname, dpi=300, transparent=myc['transparent'])
                plt.close()


class GridExampleRectPlotter(FigurePlotter):
    '''Grid field examples - plotted on the sheet.'''
    cmap = 'jet'

    def __init__(self, *args, **kwargs):
        super(GridExampleRectPlotter, self).__init__(*args, **kwargs)

    def drawSweep(self, ax, data, spaceRect):
        sweeps.plotSweep(
                data,
                None,
                ax=ax,
                cbar=True,
                cbar_kw=self.myc['cbar_kw'],
                cmap=self.cmap,
                vmin=self.myc['vmin'], vmax=self.myc['vmax'],
                sliceAnn=None,
                sigmaTitle=False)

        _, X, Y = data.getData()
        examples.drawEIRectSelection(ax, spaceRect, X, Y)

    def drawA4RectExamples(self, ns_idx, exRect, exIdx, letter=''):
        ps = self.env.ps
        dataSpace = ps.grids[ns_idx]
        noise_sigma = ps.noise_sigmas[ns_idx]
        iter_list = self.config['iter_list']
        example_idx = self.config['grids']['example_idx']

        fig = plt.figure(figsize=(8.27, 11.69))
        margin    = 0.05
        sw_width  = 0.25
        sw_left   = margin
        sw_bottom = 0.85
        sw_right  = sw_left + sw_width
        sw_top    = 0.96
        div       = 0.05

        letter_left = 0.02
        letter_top_off = 0.01
        letter_va='bottom'
        letter_ha='left'

        nsX = 0.9
        nsY = sw_top

        sweepsRect = sw_left, sw_bottom, sw_right-sw_left, sw_top-sw_bottom
        ax_sweeps = fig.add_axes(sweepsRect)
        data = aggr.GridnessScore(ps.grids[ns_idx], iter_list,
                                  ignoreNaNs=True, normalizeTicks=True,
                                  r=example_idx[ns_idx][0],
                                  c=example_idx[ns_idx][1])
        self.drawSweep(ax_sweeps, data, exRect)
        fig.text(letter_left, sw_top+letter_top_off, letter, va=letter_va,
                 ha=letter_ha, fontsize=19, fontweight='bold')

        gsCoords = 0.06, 0.05, 0.97, sw_bottom-div
        #gsCoords = margin, 0.46, 0.5, sw_bottom-div
        gs = examples.drawGridExamples(dataSpace, exRect, iter_list,
                                       gsCoords=gsCoords, exIdx=exIdx, fig=fig,
                                       maxRate=True, rateStr='',
                                       rasterized=True)
        gs.update(wspace=0.05)
        noise_sigma_txt = "$\sigma_{{noise}}$ = {0} pA".format(int(noise_sigma))
        fig.text(nsX, nsY, noise_sigma_txt, va='center', ha='right', fontsize=19)
        return fig

    def plot(self, *args, **kwargs):
        ps = self.env.ps
        iter_list = self.config['iter_list']
        YXRC = [(1, 22), (1, 22), (1, 22)] # (row, col)

        exWidth = 16
        exHeight = 16

        exampleRC = [
                [[0, 0], [16, 0], [0, 16], [16, 16]],
                [[0, 0], [16, 0], [0, 16], [16, 16]],
                [[0, 0], [16, 0], [0, 16], [16, 16]],
        ]

        saver = self.myc['fig_saver']
        saver.set_file_name(self.config['output_dir'] +
                            "/suppFigure_grid_examples")
        saver.ext = "pdf"
        saver.set_backend_params(dpi=150, transparent=True)

        strIdx = 0
        for noise_idx, noise_sigma in enumerate(ps.noise_sigmas):
            for exampleIdx, RC in enumerate(exampleRC[noise_idx]):
                print(noise_idx, exampleIdx, RC)

                exLeft   = RC[0]
                exBottom = RC[1]
                exRect = [exLeft, exBottom, exLeft+exWidth-1,
                          exBottom + exHeight - 1]
                fig = self.drawA4RectExamples(
                        noise_idx, exRect, YXRC[noise_idx],
                        letter=string.ascii_uppercase[strIdx])

                saver.savefig(fig)
                plt.close()
                strIdx += 1

        saver.close()


class GridExampleColorbarPlotter(FigurePlotter):
    def __init__(self, *args, **kwargs):
        super(GridExampleColorbarPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        fig = self._get_final_fig(self.myc['fig_size'])
        ax_cbar = fig.add_axes([0.05, 0.07, 0.1, 0.8])
        cbar = mpl.colorbar.ColorbarBase(ax_cbar, cmap=mpl.cm.jet,
                                         norm=mpl.colors.Normalize(vmin=0,
                                                                   vmax=1),
                                         ticks=[0, 1],
                                         orientation='vertical')
        ax_cbar.yaxis.set_ticklabels(['0', 'Max'])
        fname_cbar = self.get_fname("grids_examples_colorbar.pdf")
        plt.savefig(fname_cbar, dpi=300, transparent=True)
        plt.close()


def openJob(rootDir, noise_sigma):
    fileTemplate = "noise_sigma{0}_output.h5"
    fileName = rootDir + '/' + fileTemplate.format(int(noise_sigma))
    return DataStorage.open(fileName, 'r')


def drawVm(data, noise_sigma, xScaleBar=None, yScaleBar=None,
        ax=plt.gca(), sigmaTitle=True):
    '''Membrane potential examples.'''
    yScaleX = 0.5
    yScaleY = 1.1
    yScaleXOffset = 0.06
    scaleTextSize = 'x-small'

    plotTStart = 5e3
    plotTEnd   = 5.25e3

    stateYlim = [-80, -40]

    theta_start_t = getOption(data, 'theta_start_t')
    #theta_start_t = 1e3
    simTime = getOption(data, 'time')

    mon_e = data['stateMon_e']

    # E cell Vm
    t, VmMiddle = extractSummedSignals(mon_e, ['V_m'], plotTStart,
            plotTEnd)
    plotStateSignal(ax, t, VmMiddle, labely='', color='red',
            scaleBar=xScaleBar, scaleText="ms", scaleX=0.5, scaleY=1.1,
            scaleTextSize=scaleTextSize)
    if (yScaleBar is not None):
        low_level.yScaleBar(yScaleBar, yScaleX, yScaleY,
                ax=ax,
                unitsText='mV', textXOffset=yScaleXOffset,
                size='x-small')
    ax.yaxis.set_visible(False)
    ax.spines['left'].set_visible(False)
    plt.ylim(stateYlim)

    if (sigmaTitle):
        ax.set_title('$\sigma$ = {0} pA'.format(int(noise_sigma)), loc='left',
                size=scaleTextSize, y=0.95)


class VmExamplesPlotter(FigurePlotter):
    def __init__(self, *args, **kwargs):
        super(VmExamplesPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        myc = self._get_class_config()
        ps = self.env.ps

        l, b, r, t = myc['ax_rect']
        VmExampleXScalebar = 50 # ms
        VmExampleYScalebar = 10 # mV

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            fig = self._get_final_fig(myc['fig_size'])
            ax = fig.add_axes(Bbox.from_extents(l, b, r, t))
            ds = openJob(self.config['singleDataRoot'], noise_sigma)
            kw = {}
            if (ns_idx == 1):
                kw['xScaleBar'] = VmExampleXScalebar
                kw['yScaleBar'] = VmExampleYScalebar
            drawVm(ds, noise_sigma=noise_sigma, ax=ax, **kw)

            fname = self.config['output_dir'] + "/grids_Vm_example_{0}.pdf"
            plt.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)
            plt.close()


##############################################################################
EI13Root  = 'simulation_data/submission/detailed_noise/grids/EI-1_3'
EI31Root  = 'simulation_data/submission/detailed_noise/grids/EI-3_1'
detailedShape = (31, 9)

EI13PS = JobTrialSpace2D(detailedShape, EI13Root)
EI31PS = JobTrialSpace2D(detailedShape, EI31Root)
detailedNTrials = 1


detailFigSize = (4.5, 2.8)
detailLeft   = 0.15
detailBottom = 0.2
detailRight  = 0.95
detailTop    = 0.8

class GridDetailedNoisePlotter(FigurePlotter):
    '''Detailed noise plots.'''
    def __init__(self, *args, **kwargs):
        super(GridDetailedNoisePlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ylabelPos = -0.15

        types = ('grids', 'gridnessScore')
        fig = self._get_final_fig(detailFigSize)
        ax = fig.add_axes(Bbox.from_extents(detailLeft, detailBottom, detailRight,
            detailTop))
        _, p13, l13 = details.plotDetailedNoise(EI13PS, detailedNTrials, types, ax=ax,
                ylabelPos=ylabelPos,
                color='red', markerfacecolor='red', zorder=10)
        _, p31, l31 = details.plotDetailedNoise(EI31PS, detailedNTrials, types, ax=ax,
                ylabel='Gridness score', ylabelPos=ylabelPos,
                color='#505050')
        ax.yaxis.set_major_locator(ti.MultipleLocator(0.4))
        ax.yaxis.set_minor_locator(ti.MultipleLocator(0.2))
        ax.set_ylim([-0.6, 1.2])
        leg = self.myc['legend']
        l = ax.legend([p31, p13], leg, **self.myc['legend_kwargs'])

        fname = self.config['output_dir'] + "/grids_detailed_noise_gscore.pdf"
        plt.savefig(fname, dpi=300, transparent=True)
        plt.close()


gridTypes = ['grids', 'gridnessScore']
class GridsDiffSweep(SweepPlotter):
    '''Parameter sweep of the difference between noise_150 and noise_0.'''
    def __init__(self, *args, **kwargs):
        super(GridsDiffSweep, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        myc= self._get_class_config()
        sweepc = self._get_sweep_config()
        ps = self.env.ps
        example_idx = self.config['grids']['example_idx']

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas[0:-1]):
            fname = (self.config['output_dir'] +
                     "/grids_diff_sweep{0}.pdf".format(int(noise_sigma)))
            with self.figure_and_axes(fname, sweepc) as (fig, ax):
                which = ns_idx
                sweeps.plotDiffTrial(
                    ps.grids, self.config['iter_list'], which,
                    self.config['grids']['ntrials'], gridTypes,
                    ax=ax,
                    ignoreNaNs=True,
                    r=example_idx[0][0], c=example_idx[0][1],
                    cbar=True, cbar_kw=myc['cbar_kw'],
                    symmetricLimits=True,
                    cmap='RdBu_r'
                )


class GridBumpScatterPlotter(FigurePlotter):
    '''Scatter plot of gridness score vs. P(bump).'''
    def __init__(self, *args, **kwargs):
        super(GridBumpScatterPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ps = self.env.ps
        iter_list = self.config['iter_list']
        output_dir = self.config['output_dir']

        xlabel = 'P(bumps)'
        ylabel = 'Gridness score'

        isBumpData = []
        gridData = []
        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            isBumpData.append(aggr.IsBump(ps.bumpGamma[ns_idx], iter_list,
                ignoreNaNs=True, normalizeTicks=False))
            gridData.append(aggr.GridnessScore(ps.grids[ns_idx], iter_list,
                ignoreNaNs=True, normalizeTicks=False))

        fig = self._get_final_fig(self.myc['fig_size'])
        scatterPlot = scatter.FullScatterPlot(
                isBumpData, gridData, None, None, iter_list, None, None,
                s=8,
                linewidth=0.3,
                color2D=True,
                xlabel=xlabel,
                ylabel=ylabel,
                sigmaTitle=True,
                noise_sigmas=ps.noise_sigmas,
                ignoreNaNs=True,
                captionLetters=('A', 'B', 'C'),
                fig=fig)
        scatterPlot.plot(captionLeft=-0.1, plotcolorbar=False)
        scatterPlot.plotColorbar(**self.myc['color_box_coords'])
        scatterPlot.set_titleSizes(16)

        # Normal scale
        fname = output_dir + "/suppFigure_grids_vs_bumps.pdf"
        fig.savefig(fname, dpi=300)

        # Exponential scale
        for ns_idx, _ in enumerate(ps.noise_sigmas):
            ax = scatterPlot.axes[ns_idx]
            ax.set_xscale('exponential')
            ax.xaxis.set_major_locator(ti.MultipleLocator(.5))
            ax.xaxis.set_minor_locator(ti.MultipleLocator(.1))
            ax.set_xlim([-0.3, 1.002])
        fname = output_dir + "/suppFigure_grids_vs_bumps_exp.pdf"
        fig.savefig(fname, dpi=300)


class GridsPBumpsProbabilityPlotter(ProbabilityPlotter):
    '''Probability plots of gridness score vs gamma power.'''
    def __init__(self, *args, **kwargs):
        super(GridsPBumpsProbabilityPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ps = self.env.ps
        myc = self._get_class_config()
        iter_list = self.config['iter_list']
        l, b, r, t = myc['bbox_rect']
        drange = [[0, 1], [-.5, 1.2]]

        pbumps_all = np.empty(0)
        gridness_all = np.empty(0)

        # Separate noise sigmas
        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            if ns_idx == 0:
                mi_title = 'gridness score vs. P_bumps'
            else:
                mi_title = None

            pbumpsData = aggr.IsBump(ps.bumpGamma[ns_idx], iter_list,
                                     ignoreNaNs=True, collapseTrials=True)
            gridnessData = aggr.GridnessScore(ps.grids[ns_idx], iter_list,
                    normalizeTicks=False, collapseTrials=True)


            pbumpsData, _, _ = pbumpsData.getData()
            gridnessData, _, _ = gridnessData.getData()
            pbumps_all = np.hstack((pbumps_all, pbumpsData.flatten()))
            gridness_all = np.hstack((gridness_all, gridnessData.flatten()))

            # Gamma power vs. gridness score
            fig = self._get_final_fig(myc['fig_size'])
            ax = fig.add_axes(Bbox.from_extents(l, b, r, t))
            self.plotDistribution(pbumpsData, gridnessData, ax,
                                  noise_sigma=noise_sigma,
                                  range=drange,
                                  xlabel='P(bumps)',
                                  ylabel='', yticks=False,
                                  title_size=self.myc['title_size'])
            ax.axis('tight')
            ax.xaxis.set_major_locator(ti.MultipleLocator(0.5))
            if 'strip_axis' in self.myc and self.myc['strip_axis']:
                ax.axis('off')
            fname = self.config['output_dir'] + "/grids_pbumps_probability_{0}.pdf"
            fig.savefig(fname.format(int(noise_sigma)), dpi=300,
                             transparent=True)
            plt.close(fig)

            self.mutual_information(pbumpsData, gridnessData,
                                    noise_sigma=noise_sigma,
                                    title=mi_title)

        # All together
        fig = self._get_final_fig(myc['fig_size'])
        ax = fig.add_axes(Bbox.from_extents(l, b, r, t))
        self.plotDistribution(pbumps_all, gridness_all, ax,
                xlabel='P(bumps)',
                ylabel='Gridness score',
                range=drange,
                title_size=self.myc['title_size'])
        ax.axis('tight')
        ax.xaxis.set_major_locator(ti.MultipleLocator(0.5))
        if 'strip_axis' in self.myc and self.myc['strip_axis']:
            ax.axis('off')
        fname = self.config['output_dir'] + "/grids_pbumps_probability_all.pdf"
        fig.savefig(fname, dpi=300, transparent=True)
        plt.close(fig)

        self.mutual_information(pbumps_all, gridness_all)


class GridSimpleExamplePlotter(FigurePlotter):
    arenaDiam = 180.0 # cm

    def __init__(self, *args, **kwargs):
        super(GridSimpleExamplePlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        myc = self._get_class_config()
        ps  = self.env.ps

        ns_idx = self.myc['ns_idx']
        r, c = self.myc['rc']
        trial_no = self.myc['trial_no']

        fig = self._get_final_fig(myc['fig_size'])

        # Spikes
        ax_spikes = fig.add_subplot(1, 3, 1)
        data = ps.grids[ns_idx][r][c][trial_no].data['analysis']
        plotSpikes2D(
            data['spikes_e'],
            data['rat_pos_x'],
            data['rat_pos_y'],
            data['rat_dt'],
            ax=ax_spikes,
            spikeDotSize=2)

        # Rate map
        ax_field = fig.add_subplot(1, 3, 2)
        plotGridRateMap(
            data['rateMap_e'],
            data['rateMap_e_X'],
            data['rateMap_e_Y'],
            self.arenaDiam,
            rasterized=True,
            ax=ax_field)

        # Rate map dimmed
        ax_field = fig.add_subplot(1, 3, 3)
        plotGridRateMap(
            data['rateMap_e'],
            data['rateMap_e_X'],
            data['rateMap_e_Y'],
            self.arenaDiam,
            rasterized=True,
            maxRate=False,
            ax=ax_field,
            scaleBar=50,
            alpha=0.6)

        # Auto correlation
        #ax_corr = fig.add_subplot(1, 3, 3)
        #plotAutoCorrelation(
        #    data['corr'],
        #    data['corr_X'],
        #    data['corr_Y'],
        #    self.arenaDiam,
        #    rasterized=True,
        #    ax=ax_corr,
        #    scaleBar=50)

        fname = (self.config['output_dir'] +
                 'grid_spiking_and_field_example.pdf')
        fig.savefig(fname, dpi=300, transparent=myc['transparent'])
        plt.close()


class HighGridScoreFraction(Computation):
    '''Generate statistics of how many simulations in the Sweep are above
    threshold.'''
    def __init__(self, *args, **kwargs):
        super(HighGridScoreFraction, self).__init__(*args, **kwargs)

    def run_all(self, *args, **kwargs):
        ps = self.env.ps
        example_idx = self.config['grids']['example_idx']
        iter_list = self.config['iter_list']
        threshold = self.myc['threshold']

        print("Gridness score > %f" % threshold)
        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            data = aggr.GridnessScore(ps.grids[ns_idx], iter_list,
                                      ignoreNaNs=True, normalizeTicks=True,
                                      r=example_idx[ns_idx][0],
                                      c=example_idx[ns_idx][1])
            gscore, _, _ = data.getData()
            above = np.sum(gscore > threshold)
            print("\tsigma: %s pA, %d/%d" % (str(noise_sigma),
                                             above,
                                             gscore.size))

