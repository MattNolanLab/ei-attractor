from __future__ import absolute_import, print_function

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.ticker as ti
from matplotlib.transforms import Bbox

from grid_cell_model.parameters           import JobTrialSpace2D, DataSpace
from grid_cell_model.plotting.global_defs import globalAxesSettings, prepareLims
from simtools.plotting.plotters import FigurePlotter

from ..EI_plotting      import sweeps, examples, details, scatter
from ..EI_plotting      import aggregate as aggr, scaling
from ..EI_plotting.base import plotOneHist, NoiseDataSpaces
from .base import SweepPlotter, ProbabilityPlotter, DummyPlotter

__all__ = [
    'GammaSweepsPlotter',
    'GammaDetailedNoisePlotter',
    'GammaExamplePlotter',
    'GammaScatterAllPlotter',
    'GammaFreqGridsScatterAllPlotter',
    'ScatterGammaGridsSeparatePlotter',
    'ScatterGammaFGridsSeparatePlotter',
    'GammaScatterPBumpsAllPlotter',
    'GammaPBumpsProbabilityPlotter',
    'GammaFreqPBumpsProbabilityPlotter',
    'GammaGridsProbabilityPlotter',
    'GammaFreqGridsProbabilityPlotter',
]


NTrials = 5

###############################################################################

def extractACExample(sp, r, c, trialNum):
    data = sp[r][c][trialNum].data
    ac = data['analysis']['acVec'][0]
    dt = data['stateMonF_e'][0]['interval']
    freq = data['analysis']['freq'][0]
    acVal = data['analysis']['acVal'][0]
    noise_sigma = data['options']['noise_sigma']
    return ac, dt, freq, acVal, noise_sigma


def aggregateBar2(spList, varLists, trialNumList, func=(None, None)):
    vars = ([], [])
    noise_sigma = []
    for idx in xrange(len(spList)):
        for varIdx in range(len(varLists)):
            f = func[varIdx]
            if f is None:
                f = lambda x: x
            vars[varIdx].append(f(aggr.aggregate2DTrial(spList[idx], varLists[varIdx],
                trialNumList).flatten()))
        noise_sigma.append(spList[idx][0][0][0].data['options']['noise_sigma'])

    noise_sigma = np.array(noise_sigma, dtype=int)
    return vars, noise_sigma



def getACFreqThreshold(spList, trialNumList, ACThr):
    varLists = [['acVal'], ['freq']]
    vars, noise_sigma = aggregateBar2(spList, varLists, trialNumList)
    AC = vars[0]
    freq = vars[1]
    ACMean   = []
    freqMean = []
    ACStd    = []
    freqStd  = []
    thrCount = []
    for spIdx in range(len(spList)):
        thrIdx = np.logical_and(AC[spIdx] >= ACThr,
                np.logical_not(np.isnan(AC[spIdx])))
        ac_filt = AC[spIdx][thrIdx]
        ACMean.append(np.mean(ac_filt))
        ACStd.append(np.std(ac_filt))
        freq_filt = freq[spIdx][thrIdx]
        freqMean.append(np.mean(freq_filt))
        freqStd.append(np.std(freq_filt))
        thrCount.append(float(len(AC[spIdx][thrIdx])) / len (AC[spIdx]))

    return (ACMean, ACStd), (freqMean, freqStd), thrCount, noise_sigma


def plotThresholdComparison(spList, trialNumList, ACThrList):
    counts = []
    noise_sigma = None
    for ACThr in ACThrList:
        _, _, thrCount, noise_sigma = getACFreqThreshold(spList, trialNumList,
                ACThr)
        counts.append(thrCount)
    counts = np.array(counts)

    print(ACThrList, counts)

    ax = plt.gca()
    globalAxesSettings(ax)
    plt.plot(ACThrList, counts, 'o-')
    plt.plot([0], [1], linestyle='None', marker='None')
    ax.set_xlabel('Correlation threshold', labelpad=5)
    ax.set_ylabel('Count', labelpad=5)
    leg = []
    for s in noise_sigma:
        leg.append("{0}".format(int(s)))
    ax.legend(leg, loc=(0.8, 0.5), title='$\sigma$ (pA)', frameon=False,
            fontsize='small')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_major_locator(ti.MultipleLocator(0.3))
    ax.yaxis.set_major_locator(ti.MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(ti.AutoMinorLocator(3))
    ax.yaxis.set_minor_locator(ti.AutoMinorLocator(2))
    ax.margins(0.025)


def plotFreqHistogram(spList, trialNumList, ylabelPos=-0.2, CThreshold=0.1):
    FVarList = ['freq']
    CVarList = ['acVal']
    noise_sigma = [0, 150, 300]
    colors = ['red', 'green', 'blue']

    ax = plt.gca()
    plt.hold('on')
    globalAxesSettings(ax)

    for idx, sp in enumerate(spList):
        F = aggr.aggregate2DTrial(sp, FVarList, trialNumList).flatten()
        C = aggr.aggregate2DTrial(sp, CVarList, trialNumList).flatten()
        filtIdx = np.logical_and(np.logical_not(np.isnan(F)), C > CThreshold)
        plotOneHist(F[filtIdx], bins=20, normed=True)
    leg = []
    for s in noise_sigma:
        leg.append("{0}".format(int(s)))
    l = ax.legend(leg, loc=(0.8, 0.5), title='$\sigma$ (pA)', frameon=False,
            fontsize='x-small', ncol=1)
    plt.setp(l.get_title(), fontsize='x-small')

    ax.set_xlabel("Oscillation frequency (Hz)")
    #ax.text(ylabelPos, 0.5, 'p(F)', rotation=90, transform=ax.transAxes,
    #        va='center', ha='right')
    ax.set_ylabel('p(Frequency)')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_major_locator(ti.MultipleLocator(20))
    ax.yaxis.set_major_locator(ti.MaxNLocator(4))
    ax.xaxis.set_minor_locator(ti.AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(ti.AutoMinorLocator(2))
    f = ti.ScalarFormatter(useMathText=True)
    f.set_scientific(True)
    f.set_powerlimits([0, 3])
    ax.yaxis.set_major_formatter(f)
    #ax.margins(0.01, 0.00)

    thStr = 'Frequencies with C > {0}'.format(CThreshold)
    ax.text(0.99, 1.1, thStr, transform=ax.transAxes, va='bottom',
            ha='right')


###############################################################################


# gamma example rows and columns

ACVarList = ['acVal']
FVarList  = ['freq']

class GammaSweepsPlotter(SweepPlotter):
    def __init__(self, *args, **kwargs):
        super(GammaSweepsPlotter, self).__init__(*args, **kwargs)

    def get_fig(self):
        if 'fig_size' in self.myc:
            return self._get_final_fig(self.myc['fig_size'])
        else:
            return super(GammaSweepsPlotter, self).get_fig()

    def plot(self, *args, **kwargs):
        sweepc = self._get_sweep_config()
        ps = self.env.ps
        ac_xticks = self.myc['AC_xticks']
        ac_yticks = self.myc['AC_yticks']
        f_xticks = self.myc['F_xticks']
        f_yticks = self.myc['F_yticks']
        iter_list = self.config['iter_list']
        grids_example_idx = self.config['grids']['example_idx']
        fname_prefix = self.config.get('fname_prefix', '')
        filter_with_gridness = self.myc.get('filter_with_gridness', False)

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            ACData = aggr.GammaAggregateData('acVal', ps.bumpGamma[ns_idx],
                                             iter_list, normalizeTicks=True,
                                             ignoreNaNs=True)
            gammaFData = aggr.GammaAggregateData('freq', ps.bumpGamma[ns_idx],
                                                 iter_list,
                                                 normalizeTicks=True)
            gridData = aggr.GridnessScore(ps.grids[ns_idx], iter_list,
                                              normalizeTicks=True,
                                              collapseTrials=True,
                                              ignoreNaNs=True,
                                              r=grids_example_idx[ns_idx][0],
                                              c=grids_example_idx[ns_idx][1])

            if filter_with_gridness:
                gridFilter = aggr.GTFilter(gridData,
                                           self.myc['gridness_threshold'])
                ACData = ACData.filter_data(gridFilter)
                gammaFData = gammaFData.filter_data(gridFilter)

            # Gamma power
            fname = self.get_fname("gamma_sweeps{ns}.pdf", ns=noise_sigma)
            with self.figure_and_axes(fname, sweepc) as (fig, ax):
                sweeps.plotSweep(
                        ACData,
                        noise_sigma=ps.noise_sigmas[ns_idx],
                        ax=ax,
                        xlabel='' if ac_xticks[ns_idx] == False else None,
                        xticks=ac_xticks[ns_idx],
                        ylabel='' if ac_yticks[ns_idx] == False else None,
                        yticks=ac_yticks[ns_idx],
                        sigmaTitle=self.myc['AC_sigma_title'],
                        cbar=self.myc['cbar'][ns_idx],
                        cbar_kw=self.myc['AC_cbar_kw'],
                        vmin=self.myc['AC_vmin'],
                        vmax=self.myc['AC_vmax'],
                        annotations=self.myc['ann'])
                if self.myc['plot_grid_contours'][ns_idx]:
                    gridData = aggr.GridnessScore(
                            ps.grids[ns_idx],
                            iter_list,
                            ignoreNaNs=True,
                            normalizeTicks=True,
                            r=grids_example_idx[ns_idx][0],
                            c=grids_example_idx[ns_idx][1])
                    contours = sweeps.Contours(gridData,
                            self.config['sweeps']['grid_contours'])
                    contours.plot(
                            ax,
                            **self.config['sweeps']['contours_kwargs'])

            # Gamma frequency
            fname = self.get_fname("gamma_freq_sweeps{ns}.pdf", ns=noise_sigma)
            with self.figure_and_axes(fname, sweepc) as (fig, ax):
                sweeps.plotACTrial(
                        gammaFData,
                        None,
                        None,
                        noise_sigma=ps.noise_sigmas[ns_idx],
                        ax=ax,
                        xlabel='' if f_xticks[ns_idx] == False else None,
                        xticks=f_xticks[ns_idx],
                        ylabel='' if f_yticks[ns_idx] == False else None,
                        yticks=f_yticks[ns_idx],
                        trialNumList=xrange(NTrials),
                        sigmaTitle=self.myc['F_sigma_title'],
                        cbar=self.myc['cbar'][ns_idx],
                        cbar_kw=self.myc['F_cbar_kw'],
                        vmin=self.myc['F_vmin'],
                        vmax=self.myc['F_vmax'],
                        annotations=self.myc['annF'])
                if self.myc['plot_grid_contours'][ns_idx]:
                    contours = sweeps.Contours(gridData,
                            self.config['sweeps']['grid_contours'])
                    contours.plot(
                            ax,
                            **self.config['sweeps']['contours_kwargs'])




#if args.threshold or args.all:
#    ###############################################################################
#    plt.figure(figsize=(3.5, 2))
#    plotThresholdComparison(ps.bumpGamma,
#            trialNumList=range(NTrials),
#            ACThrList=np.arange(0, 0.65, 0.05))
#    plt.tight_layout()
#    fname = outputDir + '/gamma_AC_threshold_comparison.pdf'
#    plt.savefig(fname, transparent=True, dpi=300)
#
#
#if args.freqHist or args.all:
#    ylabelPos = -0.16
#    fig = plt.figure(figsize=(3.7, 2.5))
#    plotFreqHistogram(ps.bumpGamma, range(NTrials), ylabelPos=ylabelPos)
#    plt.tight_layout()
#    fname = outputDir + "/gamma_freq_histograms.pdf"
#    plt.savefig(fname, dpi=300, transparent=True)
#


##############################################################################
EI13Root  = 'simulation_data/submission/detailed_noise/gamma_bump/EI-1_3'
EI31Root  = 'simulation_data/submission/detailed_noise/gamma_bump/EI-3_1'
detailedShape = (31, 9)

EI13PS = JobTrialSpace2D(detailedShape, EI13Root)
EI31PS = JobTrialSpace2D(detailedShape, EI31Root)
detailedNTrials = 5

detailFigSize = (4.6, 1.8)
detailLeft   = 0.2
detailBottom = 0.3
detailRight  = 0.98
detailTop    = 0.9

class GammaDetailedNoisePlotter(FigurePlotter):
    ylabelPos = -0.17

    def __init__(self, *args, **kwargs):
        super(GammaDetailedNoisePlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        # 1st autocorrelation peak (gamma power)
        types = ('gamma', 'acVal')
        fig = self._get_final_fig(detailFigSize)
        ax = fig.add_axes(Bbox.from_extents(detailLeft, detailBottom, detailRight,
            detailTop))
        _, p13, l13 = details.plotDetailedNoise(EI13PS, detailedNTrials, types, ax=ax,
                ylabelPos=self.ylabelPos,
                xlabel='', xticks=False,
                color='red', markerfacecolor='red', zorder=10)
        _, p31, l31 = details.plotDetailedNoise(EI31PS, detailedNTrials, types, ax=ax,
                xlabel='', xticks=False,
                ylabel='$1^{st}$\nautocorrelation\npeak', ylabelPos=self.ylabelPos,
                color='#505050')
        ax.xaxis.set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.yaxis.set_major_locator(ti.MultipleLocator(0.6))
        ax.yaxis.set_minor_locator(ti.AutoMinorLocator(6))
        ax.set_ylim(prepareLims((0, 0.6), margin=0.03))
        leg = self.myc['legend']
        l = ax.legend([p31, p13], leg, **self.myc['legend_kwargs'])
        plt.setp(l.get_title(), fontsize='small')


        fname = self.config['output_dir'] + "/gamma_detailed_noise_power.pdf"
        plt.savefig(fname, dpi=300, transparent=True)
        plt.close()

        # Gamma frequency
        types = ('gamma', 'freq')
        fig = self._get_final_fig(detailFigSize)
        ax = fig.add_axes(Bbox.from_extents(detailLeft, detailBottom, detailRight,
            detailTop))
        _, p13, l13 = details.plotDetailedNoise(EI13PS, detailedNTrials, types, ax=ax,
                ylabelPos=self.ylabelPos,
                xlabel='',
                color='red', markerfacecolor='red', zorder=10)
        _, p31, l31 = details.plotDetailedNoise(EI31PS, detailedNTrials, types, ax=ax,
                ylabel='Oscillation\nfrequency (Hz)', ylabelPos=self.ylabelPos,
                color='#505050')
        ax.yaxis.set_major_locator(ti.MultipleLocator(30))
        ax.yaxis.set_minor_locator(ti.AutoMinorLocator(3))
        ax.set_ylim(prepareLims((30, 90), margin=0.03))

        fname = self.config['output_dir'] + "/gamma_detailed_noise_freq.pdf"
        plt.savefig(fname, dpi=300, transparent=True)
        plt.close()



##############################################################################
class GammaExamplePlotter(FigurePlotter):
    def __init__(self, *args, **kwargs):
        super(GammaExamplePlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ps = self.env.ps
        example_rc = self.config['gamma']['example_rc']

        exampleTrialNum = 0
        exampleFigSize = (2, 1.1)
        exampleLeft   = 0.08
        exampleBottom = 0.2
        exampleRight  = 0.99
        exampleTop    = 0.85

        for nsIdx, ns in enumerate(ps.noise_sigmas):
            for idx, rc in enumerate(example_rc):
                exampleFName = (self.config['output_dir'] +
                                "/{0}gamma_example{1}_{2}.pdf")
                fname = exampleFName.format(self.config.get('fname_prefix', ''),
                                            ns, idx)
                fig = self._get_final_fig(exampleFigSize)
                ax = fig.add_axes(Bbox.from_extents(exampleLeft, exampleBottom,
                    exampleRight, exampleTop))
                nsAnn = None
                xscale_kw = None
                if self.myc['xscales'][idx][nsIdx]:
                    xscale_kw = self.myc['xscale_kw']
                if self.myc['sigma_titles'][idx][nsIdx]:
                    nsAnn = ns
                examples.plotGammaExample(ps.bumpGamma[nsIdx], ax=ax,
                        r=example_rc[idx][0], c=example_rc[idx][1],
                        trialNum=exampleTrialNum,
                        tStart = 2e3, tEnd=2.25e3,
                        noise_sigma=nsAnn, noise_sigma_xy=(0.95, 1),
                        xscale_kw=xscale_kw,
                        yscale_kw=self.myc['yscale_kw'][idx][nsIdx])
                plt.savefig(fname, dpi=300, transparent=True)
                plt.close()


##############################################################################
# Separate scatter plot of gridness score vs. gamma power
class ScatterGammaGridsSeparatePlotter(FigurePlotter):
    def __init__(self, *args, **kwargs):
        super(ScatterGammaGridsSeparatePlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ps = self.env.ps
        iter_list = self.config['iter_list']
        output_dir = self.config['output_dir']
        ignoreNaNs = True

        xlabel = '$1^{st}$ autocorrelation peak'
        ylabel = 'Gridness score'

        fig = self._get_final_fig(self.myc['fig_size'])

        gammaData = []
        gridData = []
        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            gammaData.append(aggr.GammaAggregateData('acVal',
                                                     ps.bumpGamma[ns_idx],
                                                     iter_list,
                                                     normalizeTicks=False,
                                                     collapseTrials=True))
            gridData.append(aggr.GridnessScore(ps.grids[ns_idx], iter_list,
                                               ignoreNaNs=True,
                                               normalizeTicks=False))

        scatterPlot = scatter.FullScatterPlot(
                gammaData, gridData, None, None, iter_list, None, None,
                s=12,
                linewidth=0.3,
                color2D=True,
                xlabel=xlabel,
                ylabel=ylabel,
                sigmaTitle=True,
                noise_sigmas=ps.noise_sigmas,
                ignoreNaNs=ignoreNaNs,
                captionLetters=('A', 'B', 'C'),
                fig=fig)
        scatterPlot.plot(captionLeft=-0.1, plotcolorbar=False)
        l = 0.87
        w = 0.12
        scatterPlot.plotColorbar(left=l, bottom=.9, right=l+w, top=1)
        scatterPlot.set_titleSizes(16)
        for ax in scatterPlot.axes:
            #ax.xaxis.set_major_locator(ti.MultipleLocator(0.2))
            ax.yaxis.set_major_locator(ti.MultipleLocator(0.3))
            #ax.set_ylim(prepareLims((-0.5, 1.2), margin=0.02))

        fname = output_dir + "/suppFigure_gamma.pdf"
        fig.savefig(fname, dpi=300, transparent=True)


##############################################################################
# Separate scatter plot of gridness score vs. gamma frequency
class ScatterGammaFGridsSeparatePlotter(FigurePlotter):
    def __init__(self, *args, **kwargs):
        super(ScatterGammaFGridsSeparatePlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ps = self.env.ps
        iter_list = self.config['iter_list']
        output_dir = self.config['output_dir']
        ignoreNaNs = True

        xlabel = 'Oscillation frequency (Hz)'
        ylabel = 'Gridness score'

        fig = self._get_final_fig(self.myc['fig_size'])

        gammaFData = []
        gridData = []
        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            gammaFData.append(aggr.GammaAggregateData('freq',
                                                     ps.bumpGamma[ns_idx],
                                                     iter_list,
                                                     normalizeTicks=False,
                                                     collapseTrials=True))
            gridData.append(aggr.GridnessScore(ps.grids[ns_idx], iter_list,
                                               ignoreNaNs=True,
                                               normalizeTicks=False))

        scatterPlot = scatter.FullScatterPlot(
                gammaFData, gridData, None, None, iter_list, None, None,
                s=12,
                linewidth=0.3,
                color2D=True,
                xlabel=xlabel,
                ylabel=ylabel,
                sigmaTitle=True,
                noise_sigmas=ps.noise_sigmas,
                ignoreNaNs=ignoreNaNs,
                captionLetters=('A', 'B', 'C'),
                fig=fig)
        scatterPlot.plot(captionLeft=-0.1, plotcolorbar=False)
        for ax in scatterPlot.axes:
            ax.yaxis.set_major_locator(ti.MultipleLocator(0.3))
        l = 0.87
        w = 0.12
        scatterPlot.plotColorbar(left=l, bottom=.22, right=l+w, top=.32)
        scatterPlot.set_titleSizes(16)

        fname = output_dir + "/suppFigure_gammaF_grids_scatter.pdf"
        fig.savefig(fname, dpi=300, transparent=True)


##############################################################################
# Scatter plot of gridness score vs. gamma power
# All in one plot
class GammaScatterAllPlotter(FigurePlotter):
    def __init__(self, *args, **kwargs):
        super(GammaScatterAllPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ps = self.env.ps
        myc = self._get_class_config()
        plot_legend = self.myc.get('plot_legend', False)
        legend_kwargs = myc['legend_kwargs']
        l, b, r, t = self.myc['bbox_rect']

        self.fig = self._get_final_fig(myc['fig_size'])
        self.ax = self.fig.add_axes(Bbox.from_extents(l, b, r, t))

        NTrialsGamma = 5
        NTrialsGrids = 3
        typesGamma = ['gamma', 'acVal']
        typesGrids = ['grids', 'gridnessScore']
        self.ax.hold('on')
        scatterColors = ['green', 'red', 'blue']
        scatterOrders = [2, 3, 1]

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            color = scatterColors[ns_idx]
            scatterPlot = scatter.ScatterPlot(
                    ps.bumpGamma[ns_idx], ps.grids[ns_idx], typesGamma,
                    typesGrids, self.config['iter_list'], NTrialsGamma, NTrialsGrids,
                    c=color,
                    s=self.myc['dot_size']*self.config['scale_factor'],
                    linewidth=0.3,
                    xlabel='$1^{st}$ autocorrelation peak',
                    ylabel=self.myc['ylabel'],
                    zorder=scatterOrders[ns_idx])
            scatterPlot.plot()
        self.ax.xaxis.set_major_locator(ti.MultipleLocator(0.2))
        self.ax.yaxis.set_major_locator(ti.MultipleLocator(0.5))
        if plot_legend:
            leg = ['0', '150', '300']
            l = self.ax.legend(leg, **legend_kwargs)
            plt.setp(l.get_title(), size=legend_kwargs['fontsize'])

    def save(self, *args, **kwargs):
        fname = self.config['output_dir'] + "/gamma_scatter_gamma_grids_all.pdf"
        self.fig.savefig(fname, dpi=300, transparent=True)


##############################################################################
# Scatter plot of gridness score vs. gamma frequency
# All in one plot
class GammaFreqGridsScatterAllPlotter(FigurePlotter):
    def __init__(self, *args, **kwargs):
        super(GammaFreqGridsScatterAllPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ps = self.env.ps
        myc = self._get_class_config()
        legend_kwargs = myc['legend_kwargs']
        iter_list = self.config['iter_list']
        l, b, r, t = self.myc['bbox_rect']

        self.fig = self._get_final_fig(myc['fig_size'])
        self.ax = self.fig.add_axes(Bbox.from_extents(l, b, r, t))

        scatterColors = ['green', 'red', 'blue']
        scatterOrders = [2, 3, 1]

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            gammaFData = aggr.GammaAggregateData('freq', ps.bumpGamma[ns_idx],
                                                 iter_list,
                                                 normalizeTicks=False,
                                                 collapseTrials=True)
            gridnessData = aggr.GridnessScore(ps.grids[ns_idx], iter_list,
                                              normalizeTicks=False,
                                              collapseTrials=True)

            color = scatterColors[ns_idx]
            scatterPlot = scatter.ScatterPlot(
                    gammaFData, gridnessData,
                    None, None, None, None, None,
                    c=color,
                    s=self.myc['dot_size']*self.config['scale_factor'],
                    linewidth=0.3,
                    xlabel='Oscillation frequency (Hz)',
                    ylabel=self.myc['ylabel'],
                    yticks=self.myc['yticks'],
                    zorder=scatterOrders[ns_idx])
            scatterPlot.plot()
        self.ax.xaxis.set_major_locator(ti.MultipleLocator(20))
        self.ax.yaxis.set_major_locator(ti.MultipleLocator(0.5))
        leg = ['0', '150', '300']
        l = self.ax.legend(leg, **legend_kwargs)
        plt.setp(l.get_title(), size=legend_kwargs['fontsize'])

    def save(self, *args, **kwargs):
        fname = self.config['output_dir'] + "/gamma_scatter_gammaF_grids_all.pdf"
        self.fig.savefig(fname, dpi=300, transparent=True)


##############################################################################
# Probability plots of gridness score vs gamma power
class GammaGridsProbabilityPlotter(ProbabilityPlotter):
    def __init__(self, *args, **kwargs):
        super(GammaGridsProbabilityPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ps = self.env.ps
        myc = self._get_class_config()
        iter_list = self.config['iter_list']
        l, b, r, t = myc['bbox_rect']
        drange = [[-.2, .8], [-.5, 1.2]]

        gamma_all = np.empty(0)
        gridness_all = np.empty(0)

        # Separate noise sigmas
        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            if ns_idx == 0:
                mi_title = 'Gamma power vs. gridness score'
            else:
                mi_title = None

            gammaData = aggr.GammaAggregateData('acVal', ps.bumpGamma[ns_idx],
                                                iter_list,
                                                normalizeTicks=False,
                                                collapseTrials=True)
            gridnessData = aggr.GridnessScore(ps.grids[ns_idx], iter_list,
                    normalizeTicks=False, collapseTrials=True)


            gammaData, _, _ = gammaData.getData()
            gridnessData, _, _ = gridnessData.getData()
            gamma_all = np.hstack((gamma_all, gammaData.flatten()))
            gridness_all = np.hstack((gridness_all, gridnessData.flatten()))

            # Gamma power vs. gridness score
            fig = self._get_final_fig(myc['fig_size'])
            ax = fig.add_axes(Bbox.from_extents(l, b, r, t))
            self.plotDistribution(gammaData, gridnessData, ax,
                                  noise_sigma=noise_sigma,
                                  range=drange,
                                  xlabel='$Power_\gamma$',
                                  ylabel='', yticks=False,
                                  title_size=self.myc['title_size'])
            ax.axis('tight')
            ax.xaxis.set_major_locator(ti.MultipleLocator(0.5))
            if 'strip_axis' in self.myc and self.myc['strip_axis']:
                ax.axis('off')
            fname = self.config['output_dir'] + "/gamma_gridness_probability_{0}.pdf"
            fig.savefig(fname.format(int(noise_sigma)), dpi=300,
                             transparent=True)
            plt.close(fig)

            self.mutual_information(gammaData, gridnessData,
                                    noise_sigma=noise_sigma,
                                    title=mi_title)

        # All together
        fig = self._get_final_fig(myc['fig_size'])
        ax = fig.add_axes(Bbox.from_extents(l, b, r, t))
        self.plotDistribution(gamma_all, gridness_all, ax,
                xlabel='$Power_\gamma$',
                ylabel='Gridness score',
                range=drange,
                title_size=self.myc['title_size'])
        ax.axis('tight')
        ax.xaxis.set_major_locator(ti.MultipleLocator(0.5))
        if 'strip_axis' in self.myc and self.myc['strip_axis']:
            ax.axis('off')
        fname = self.config['output_dir'] + "/gamma_gridness_probability_all.pdf"
        fig.savefig(fname, dpi=300, transparent=True)
        plt.close(fig)

        self.mutual_information(gamma_all, gridness_all)


##############################################################################
# Probability plots of gridness score vs gamma frequency
class GammaFreqGridsProbabilityPlotter(ProbabilityPlotter):
    def __init__(self, *args, **kwargs):
        super(GammaFreqGridsProbabilityPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ps = self.env.ps
        myc = self._get_class_config()
        iter_list = self.config['iter_list']
        l, b, r, t = myc['bbox_rect']
        drange = [[20, 160], [-.5, 1.2]]

        gammaF_all = np.empty(0)
        gridness_all = np.empty(0)

        # Separate noise sigmas
        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            if ns_idx == 0:
                mi_title = 'Gamma frequency vs. gridness score'
            else:
                mi_title = None

            gammaFData = aggr.GammaAggregateData('freq', ps.bumpGamma[ns_idx],
                                                 iter_list,
                                                 normalizeTicks=False,
                                                 collapseTrials=True)
            gridnessData = aggr.GridnessScore(ps.grids[ns_idx], iter_list,
                    normalizeTicks=False, collapseTrials=True)


            gammaFData, _, _ = gammaFData.getData()
            gridnessData, _, _ = gridnessData.getData()
            gammaF_all = np.hstack((gammaF_all, gammaFData.flatten()))
            gridness_all = np.hstack((gridness_all, gridnessData.flatten()))

            # Gamma frequency vs. gridness score
            fig = self._get_final_fig(myc['fig_size'])
            ax = fig.add_axes(Bbox.from_extents(l, b, r, t))
            self.plotDistribution(gammaFData, gridnessData, ax,
                                  noise_sigma=noise_sigma,
                                  range=drange,
                                  xlabel='Gridness score',
                                  ylabel='', yticks=False,
                                  title_size=self.myc['title_size'])
            ax.axis('tight')
            ax.xaxis.set_major_locator(ti.MultipleLocator(0.5))
            if 'strip_axis' in self.myc and self.myc['strip_axis']:
                ax.axis('off')
            fname = self.config['output_dir'] + "/gammaFreq_gridness_probability_{0}.pdf"
            fig.savefig(fname.format(int(noise_sigma)), dpi=300,
                             transparent=True)
            plt.close(fig)

            self.mutual_information(gridnessData, gammaFData,
                    noise_sigma=noise_sigma, title=mi_title)

        # All together
        fig = self._get_final_fig(myc['fig_size'])
        ax = fig.add_axes(Bbox.from_extents(l, b, r, t))
        self.plotDistribution(gammaF_all, gridness_all, ax,
                xlabel='Gridness score',
                ylabel='Oscillation frequency (Hz)',
                range=drange,
                title_size=self.myc['title_size'])
        ax.axis('tight')
        ax.xaxis.set_major_locator(ti.MultipleLocator(0.5))
        if 'strip_axis' in self.myc and self.myc['strip_axis']:
            ax.axis('off')
        fname = self.config['output_dir'] + "/gammaFreq_gridness_probability_all.pdf"
        fig.savefig(fname, dpi=300, transparent=True)
        plt.close(fig)

        self.mutual_information(gammaF_all, gridness_all)


##############################################################################
# Scatter plot of gamma power vs P_{bumps}
# All in one plot
class GammaScatterPBumpsAllPlotter(FigurePlotter):
    def __init__(self, *args, **kwargs):
        super(GammaScatterPBumpsAllPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ps = self.env.ps
        myc = self._get_class_config()
        iter_list = self.config['iter_list']
        l, b, r, t = self.myc['bbox_rect']
        legend_kwargs = myc['legend_kwargs']
        xlabel = self.myc.get('xlabel', 'P(bumps)')

        self.fig = self._get_final_fig(myc['fig_size'])
        self.ax = self.fig.add_axes(Bbox.from_extents(l, b, r, t))

        self.ax.hold('on')
        scatterColors = ['green', 'red', 'blue']
        scatterOrders = [2, 3, 1]

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            gammaData = aggr.GammaAggregateData('acVal', ps.bumpGamma[ns_idx],
                                                iter_list,
                                                normalizeTicks=False,
                                                collapseTrials=False)
            pbumpsData = aggr.IsBump(ps.bumpGamma[ns_idx], iter_list,
                                     ignoreNaNs=True, collapseTrials=False)
            color = scatterColors[ns_idx]
            scatterPlot = scatter.ScatterPlot(
                    pbumpsData, gammaData, None, None, None, None, None,
                    c=color,
                    s=6*self.config['scale_factor'],
                    linewidth=0.3,
                    xlabel=xlabel,
                    ylabel='$1^{st}$ autocorrelation peak',
                    zorder=scatterOrders[ns_idx])
            scatterPlot.plot()
        self.ax.xaxis.set_major_locator(ti.MultipleLocator(0.2))
        self.ax.yaxis.set_major_locator(ti.MultipleLocator(0.2))
        leg = ['0', '150', '300']
        l = self.ax.legend(leg, **legend_kwargs)
        plt.setp(l.get_title(), size=legend_kwargs['fontsize'])

    def save(self, *args, **kwargs):
        # Linear scale
        fname = self.config['output_dir'] + "/gamma_scatter_gamma_pbumps_all.pdf"
        self.fig.savefig(fname, dpi=300, transparent=True)

        # Exponential scale
        fname = self.config['output_dir'] + "/gamma_scatter_gamma_pbumps_all_exp.pdf"
        self.ax.set_xscale('exponential')
        self.ax.xaxis.set_major_locator(ti.MultipleLocator(.5))
        self.ax.xaxis.set_minor_locator(ti.MultipleLocator(.1))
        self.ax.set_xlim([-0.3, 1.002])
        self.fig.savefig(fname, dpi=300, transparent=True)


class GammaPBumpsProbabilityPlotter(ProbabilityPlotter):
    def __init__(self, *args, **kwargs):
        super(GammaPBumpsProbabilityPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ps = self.env.ps
        myc = self._get_class_config()
        iter_list = self.config['iter_list']
        l, b, r, t = myc['bbox_rect']

        gamma_all = np.empty(0)
        pbumps_all = np.empty(0)

        # Separate noise sigmas
        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            if ns_idx == 0:
                mi_title = 'Gamma power vs. P_Bumps'
            else:
                mi_title = None

            gammaData = aggr.GammaAggregateData('acVal', ps.bumpGamma[ns_idx],
                                                iter_list,
                                                normalizeTicks=False,
                                                collapseTrials=False)
            pbumpsData = aggr.IsBump(ps.bumpGamma[ns_idx], iter_list,
                                     ignoreNaNs=True, collapseTrials=False)
            gammaData, _, _ = gammaData.getData()
            pbumpsData, _, _ = pbumpsData.getData()
            gamma_all = np.hstack((gamma_all, gammaData.flatten()))
            pbumps_all = np.hstack((pbumps_all, pbumpsData.flatten()))

            fig = self._get_final_fig(myc['fig_size'])
            ax = fig.add_axes(Bbox.from_extents(l, b, r, t))
            self.plotDistribution(pbumpsData, gammaData, ax,
                                  noise_sigma=noise_sigma,
                                  ylabel='', yticks=False)
            fname = self.config['output_dir'] + "/gamma_pbumps_probability_{0}.pdf"
            fig.savefig(fname.format(int(noise_sigma)), dpi=300,
                             transparent=True)
            plt.close(fig)

            self.mutual_information(gammaData, pbumpsData,
                    noise_sigma=noise_sigma,
                    title=mi_title)

        # All together
        fig = self._get_final_fig(myc['fig_size'])
        ax = fig.add_axes(Bbox.from_extents(l, b, r, t))
        self.plotDistribution(pbumps_all, gamma_all, ax)
        #fig.tight_layout(**myc['tight_layout_kwargs'])
        fname = self.config['output_dir'] + "/gamma_pbumps_probability_all.pdf"
        fig.savefig(fname, dpi=300, transparent=True)
        plt.close(fig)

        self.mutual_information(gamma_all, pbumps_all)


##############################################################################
# Gamma frequency vs. P_Bumps
class GammaFreqPBumpsProbabilityPlotter(ProbabilityPlotter):
    def __init__(self, *args, **kwargs):
        super(GammaFreqPBumpsProbabilityPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ps = self.env.ps
        myc = self._get_class_config()
        iter_list = self.config['iter_list']
        l, b, r, t = myc['bbox_rect']
        drange = [[0, 1], [0, 140]]

        gammaF_all = np.empty(0)
        pbumps_all = np.empty(0)

        # Separate noise sigmas
        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            if ns_idx == 0:
                mi_title = 'Gamma freq vs. P_Bumps'
            else:
                mi_title = None

            gammaFData = aggr.GammaAggregateData('freq', ps.bumpGamma[ns_idx],
                                                 iter_list,
                                                 normalizeTicks=False,
                                                 collapseTrials=False)
            pbumpsData = aggr.IsBump(ps.bumpGamma[ns_idx], iter_list,
                                     ignoreNaNs=True, collapseTrials=False)
            gammaFData, _, _ = gammaFData.getData()
            pbumpsData, _, _ = pbumpsData.getData()
            gammaF_all = np.hstack((gammaF_all, gammaFData.flatten()))
            pbumps_all = np.hstack((pbumps_all, pbumpsData.flatten()))

            fig = self._get_final_fig(myc['fig_size'])
            ax = fig.add_axes(Bbox.from_extents(l, b, r, t))
            self.plotDistribution(pbumpsData, gammaFData, ax,
                                  noise_sigma=noise_sigma,
                                  range=drange,
                                  ylabel='', yticks=False)
            fname = self.config['output_dir'] + "/gammaFreq_pbumps_probability_{0}.pdf"
            fig.savefig(fname.format(int(noise_sigma)), dpi=300,
                             transparent=True)
            plt.close(fig)

            self.mutual_information(gammaFData, pbumpsData,
                    noise_sigma=noise_sigma,
                    title=mi_title)

        # All together
        fig = self._get_final_fig(myc['fig_size'])
        ax = fig.add_axes(Bbox.from_extents(l, b, r, t))
        self.plotDistribution(pbumps_all, gammaF_all, ax,
                range=drange,
                ylabel='Oscillation frequency (Hz)')
        #fig.tight_layout(**myc['tight_layout_kwargs'])
        fname = self.config['output_dir'] + "/gammaFreq_pbumps_probability_all.pdf"
        fig.savefig(fname, dpi=300, transparent=True)
        plt.close(fig)

        self.mutual_information(gammaF_all, pbumps_all)
