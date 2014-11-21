'''
Figure illustrating seizures.
'''
from __future__ import absolute_import, print_function

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ti
from matplotlib.colors import LogNorm
from matplotlib.transforms import Bbox
import matplotlib.gridspec as gridspec

from ..EI_plotting import sweeps, rasters, base, scatter
from ..EI_plotting import aggregate as aggr
from .base import FigurePlotter, SweepPlotter, ProbabilityPlotter

__all__ = [
    'EIRasterPlotter',
    'EIRatePlotter',
    'MaxPopulationFRSweepsPlotter',
    'MaxMeanThetaFRSweepPlotter',
    'MaxStdThetaFRSweepPlotter',
    'MaxMedianThetaFRSweepPlotter',
    'MaxThetaFRHistPlotter',
    'PSeizureSweepPlotter',
    'MaxFRGridsProbabilityPlotter',
    'MaxFRGridsScatterAllPlotter',
    'PSeizureGridsProbabilityPlotter',
    'PSeizureGridsScatterAllPlotter',
    'RasterExamplePlotter'
]

##############################################################################
rasterRC      = [(5, 15), (5, 15), (5, 15)] # (row, col)
tLimits = [2e3, 2.25e3] # ms

transparent   = True
rasterLeft    = 0.28
rasterBottom  = 0.15
rasterRight   = 0.95
rasterTop     = 0.8


class EIRasterPlotter(FigurePlotter):
    def __init__(self, *args, **kwargs):
        super(EIRasterPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ps = self.env.ps

        output_dir = self.config['output_dir']

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            scaleBar = 25 if ns_idx == 2 else None
            fig = self._get_final_fig(self.myc['fig_size'])
            ax = fig.add_axes(Bbox.from_extents(rasterLeft, rasterBottom, rasterRight,
                rasterTop))
            rasters.EIRaster(ps.bumpGamma[ns_idx],
                    noise_sigma=noise_sigma,
                    spaceType='bump',
                    r=rasterRC[ns_idx][0], c=rasterRC[ns_idx][1],
                    ylabelPos=self.myc['ylabelPos'],
                    tLimits=tLimits,
                    markersize=2*self.config['scale_factor'],
                    ylabel='' if self.myc['yticks'][ns_idx] == False else None,
                    yticks=self.myc['yticks'][ns_idx],
                    scaleBar=scaleBar, scaleX=.85, scaleY=-.1,
                    scaleTextYOffset=.03, scaleHeight=.005,
                    ann_EI=True)

            fname = "%s/bumps_raster%d.%s" % (output_dir, int(noise_sigma),
                                              self.myc['fig_ext'])
            fig.savefig(fname, dpi=300, transparent=transparent)
            plt.close()


##############################################################################

class EIRatePlotter(FigurePlotter):
    def __init__(self, *args, **kwargs):
        super(EIRatePlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ps = self.env.ps

        output_dir = self.config['output_dir']

        rateLeft    = rasterLeft
        rateBottom  = 0.2
        rateRight   = rasterRight
        rateTop     = self.myc['rateTop']

        for idx, noise_sigma in enumerate(ps.noise_sigmas):
            # E cells
            fig = self._get_final_fig(self.myc['fig_size'])
            ax = fig.add_axes(Bbox.from_extents(rateLeft, rateBottom, rateRight,
                rateTop))
            kw = {}
            if (idx != 0):
                kw['ylabel'] = ''

            rasters.plotAvgFiringRate(ps.bumpGamma[idx],
                    spaceType='bump',
                    noise_sigma=ps.noise_sigmas[idx],
                    popType='E',
                    r=rasterRC[idx][0], c=rasterRC[idx][1],
                    ylabelPos=self.myc['ylabelPos'],
                    color='red',
                    tLimits=tLimits,
                    ax=ax, **kw)
            fname = output_dir + "/bumps_rate_e{0}.pdf".format(noise_sigma)
            fig.savefig(fname, dpi=300, transparent=transparent)
            plt.close()

            # I cells
            fig = self._get_final_fig(self.myc['fig_size'])
            ax = fig.add_axes(Bbox.from_extents(rateLeft, rateBottom, rateRight,
                rateTop))
            kw = {}
            if (idx != 0):
                kw['ylabel'] = ''

            rasters.plotAvgFiringRate(ps.bumpGamma[idx],
                    spaceType='bump',
                    noise_sigma=ps.noise_sigmas[idx],
                    popType='I',
                    r=rasterRC[idx][0], c=rasterRC[idx][1],
                    ylabelPos=self.myc['ylabelPos'],
                    color='blue',
                    tLimits=tLimits,
                    ax=ax, **kw)
            fname = output_dir + "/bumps_rate_i{0}.pdf".format(noise_sigma)
            fig.savefig(fname, dpi=300, transparent=transparent)
            plt.close()




##############################################################################
# Seizure measure - max firing rate for the whole simulation
maxFR_vmin = 0
maxFR_vmax = 500.

class MaxPopulationFRSweepsPlotter(SweepPlotter):
    def __init__(self, *args, **kwargs):
        super(MaxPopulationFRSweepsPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        myc= self._get_class_config()
        sweepc = self._get_sweep_config()
        ps = self.env.ps
        iter_list = self.config['iter_list']

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            fname = (self.config['output_dir'] +
                     "/bumps_popMaxFR_sweep{0}.pdf".format(int(noise_sigma)))
            with self.figure_and_axes(fname, sweepc) as (fig, ax):
                kw = dict()
                if ns_idx != 0:
                    kw['ylabel'] = ''
                    kw['yticks'] = False
                data = aggr.MaxPopulationFR(ps.bumpGamma[ns_idx], iter_list,
                        ignoreNaNs=True, normalizeTicks=True)
                _, _, cax = sweeps.plotSweep(data,
                        noise_sigma=noise_sigma,
                        ax=ax,
                        cbar=self.myc['cbar'][ns_idx],
                        cbar_kw=myc['cbar_kw'],
                        vmin=maxFR_vmin, vmax=maxFR_vmax,
                        **kw)

                # Contours
                if self.myc['plot_grid_contours'][ns_idx]:
                    grids_example_idx = self.config['grids']['example_idx']
                    gridData = aggr.GridnessScore(ps.grids[ns_idx], iter_list,
                                                  ignoreNaNs=True, normalizeTicks=True,
                                                  r=grids_example_idx[ns_idx][0],
                                                  c=grids_example_idx[ns_idx][1])
                    contours = sweeps.Contours(gridData,
                            self.myc['grid_contours'])
                    contours.plot(
                            ax,
                            **self.config['sweeps']['contours_kwargs'])

##############################################################################
# Seizure measure - max firing rate per theta cycle
# mean
class MaxMeanThetaFRSweepPlotter(SweepPlotter):
    def __init__(self, *args, **kwargs):
        super(MaxMeanThetaFRSweepPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        myc= self._get_class_config()
        sweepc = self._get_sweep_config()
        ps = self.env.ps
        iter_list = self.config['iter_list']

        maxThetaFR_vmin = 2.
        maxThetaFR_vmax = 500.

        thetaT = self.config['seizures']['thetaT']
        sig_dt = self.config['seizures']['sig_dt']

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            fname = (self.config['output_dir'] +
                    "/bumps_popMaxThetaFR_sweep{0}.pdf".format(int(noise_sigma)))
            with self.figure_and_axes(fname, sweepc) as (fig, ax):
                kw = dict(cbar=False)
                if ns_idx != 0:
                    kw['ylabel'] = ''
                    kw['yticks'] = False
                if ns_idx == 0:
                    kw['cbar'] = True
                data = aggr.MaxThetaPopulationFR(
                        thetaT, sig_dt, np.mean,
                        ps.bumpGamma[ns_idx], iter_list,
                        ignoreNaNs=True, normalizeTicks=True)
                _, _, cax = sweeps.plotSweep(data,
                        noise_sigma=noise_sigma,
                        ax=ax,
                        cbar_kw=myc['cbar_kw'],
                        vmin=maxThetaFR_vmin, vmax=maxThetaFR_vmax,
                        #norm=LogNorm(maxThetaFR_vmin, maxThetaFR_vmax),
                        sigmaTitle=False,
                        **kw)


# median
class MaxMedianThetaFRSweepPlotter(SweepPlotter):
    def __init__(self, *args, **kwargs):
        super(MaxMedianThetaFRSweepPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        myc = self._get_class_config()
        sweepc = self._get_sweep_config()
        ps = self.env.ps
        iter_list = self.config['iter_list']

        thetaT = self.config['seizures']['thetaT']
        sig_dt = self.config['seizures']['sig_dt']

        vmin = 2.
        vmax = 500.

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            fname = (self.config['output_dir'] +
                     "/bumps_popMaxThetaFR_median_sweep{0}.pdf".format(int(noise_sigma)))
            with self.figure_and_axes(fname, sweepc) as (fig, ax):
                kw = dict(cbar=False)
                if ns_idx != 0:
                    kw['ylabel'] = ''
                    kw['yticks'] = False
                if ns_idx == 0:
                    kw['cbar'] = True
                data = aggr.MaxThetaPopulationFR(
                        thetaT, sig_dt, np.median,
                        ps.bumpGamma[ns_idx], iter_list,
                        ignoreNaNs=True, normalizeTicks=True)
                _, _, cax = sweeps.plotSweep(data,
                        noise_sigma=noise_sigma,
                        ax=ax,
                        cbar_kw=myc['cbar_kw'],
                        vmin=vmin, vmax=vmax,
                        sigmaTitle=False,
                        **kw)


# std
class MaxStdThetaFRSweepPlotter(SweepPlotter):
    def __init__(self, *args, **kwargs):
        super(MaxStdThetaFRSweepPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        myc = self._get_class_config()
        sweepc = self._get_sweep_config()
        ps = self.env.ps
        iter_list = self.config['iter_list']

        thetaT = self.config['seizures']['thetaT']
        sig_dt = self.config['seizures']['sig_dt']

        maxThetaFR_std_vmin = 0.
        maxThetaFR_std_vmax = None

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            fname = (self.config['output_dir'] +
                     "/bumps_popMaxThetaFR_std_sweep{0}.pdf".format(int(noise_sigma)))
            with self.figure_and_axes(fname, sweepc) as (fig, ax):
                kw = dict(cbar=True)
                data = aggr.MaxThetaPopulationFR(
                        thetaT, sig_dt, np.std,
                        ps.bumpGamma[ns_idx], iter_list,
                        ignoreNaNs=True, normalizeTicks=True)
                _, _, cax = sweeps.plotSweep(data,
                        noise_sigma=noise_sigma,
                        ax=ax,
                        cbar_kw=myc['cbar_kw'],
                        vmin=maxThetaFR_std_vmin, vmax=maxThetaFR_std_vmax,
                        sigmaTitle=False,
                        **kw)


##############################################################################
# Proportion of cycles with max firing rate larger than threshold, i.e. number
# of seizures during the simulation
class thresholdReduction(object):
    def __init__(self, threshold):
        self.threshold = threshold

    def __call__(self, data):
        return float(np.count_nonzero(data >= self.threshold)) / len(data)

class PSeizureSweepPlotter(SweepPlotter):
    def __init__(self, *args, **kwargs):
        super(PSeizureSweepPlotter, self).__init__(*args, **kwargs)

    def get_fig(self):
        if 'fig_size' in self.myc:
            return self._get_final_fig(self.myc['fig_size'])
        else:
            return super(PSeizureSweepPlotter, self).get_fig()

    def plot(self, *args, **kwargs):
        myc = self._get_class_config()
        sweepc = self._get_sweep_config()
        ps = self.env.ps
        iter_list = self.config['iter_list']

        FRThreshold = myc['FRThreshold']
        thetaT = self.config['seizures']['thetaT']
        sig_dt = self.config['seizures']['sig_dt']

        vmin = 0.
        vmax = 1.

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            fname = (self.config['output_dir'] +
                     "/bumps_seizureProportion_sweep{0}.pdf".format(int(noise_sigma)))
            with self.figure_and_axes(fname, sweepc) as (fig, ax):
                kw = dict(cbar=False)
                if ns_idx != 0:
                    kw['ylabel'] = ''
                    kw['yticks'] = False
                if ns_idx == 0:
                    kw['cbar'] = True
                data = aggr.MaxThetaPopulationFR(
                        thetaT, sig_dt, thresholdReduction(FRThreshold),
                        ps.bumpGamma[ns_idx], iter_list,
                        ignoreNaNs=True, normalizeTicks=True)
                _, _, cax = sweeps.plotSweep(data,
                        noise_sigma=noise_sigma,
                        ax=ax,
                        cbar_kw=myc['cbar_kw'],
                        vmin=vmin, vmax=vmax,
                        sigmaTitle=False,
                        **kw)
                # Contours
                if self.myc['plot_grid_contours'][ns_idx]:
                    grids_example_idx = self.config['grids']['example_idx']
                    gridData = aggr.GridnessScore(ps.grids[ns_idx], iter_list,
                                                  ignoreNaNs=True, normalizeTicks=True,
                                                  r=grids_example_idx[ns_idx][0],
                                                  c=grids_example_idx[ns_idx][1])
                    contours = sweeps.Contours(gridData,
                            self.myc['grid_contours'])
                    contours.plot(
                            ax,
                            **self.config['sweeps']['contours_kwargs'])


##############################################################################
#       Histograms of maxima of firing rates during theta cycles
class MaxThetaFRHistPlotter(FigurePlotter):
    def __init__(self, *args, **kwargs):
        super(MaxThetaFRHistPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ps = self.env.ps
        output_dir = self.config['output_dir']
        iter_list = self.config['iter_list']

        thetaT = self.config['seizures']['thetaT']
        sig_dt = self.config['seizures']['sig_dt']

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            fig = self._get_final_fig(self.config['sweeps']['fig_size'])
            ax = plt.subplot(111)
            kw = dict(cbar=True)
            data = aggr.MaxThetaPopulationFR(
                    thetaT, sig_dt, None,
                    ps.bumpGamma[ns_idx], iter_list,
                    ignoreNaNs=True, normalizeTicks=True)
            flatData = np.hstack(data.getNonReducedData().flatten().tolist())
            flatData = flatData[np.logical_not(np.isnan(flatData))]
            base.plotOneHist(flatData, bins=80, ax=ax, rwidth=.8, normed=True)
            ax.set_xlabel('Max rate in $\\theta$ cycle (Hz)')
            ax.set_ylabel('p(rate)')
            fig.tight_layout()
            fname = output_dir + "/bumps_popMaxThetaFR_hist{0}.pdf"
            fig.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)
            plt.close()



##############################################################################
# Probability plots of gridness score vs max. firing rate
class MaxFRGridsProbabilityPlotter(ProbabilityPlotter):
    def __init__(self, *args, **kwargs):
        super(MaxFRGridsProbabilityPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ps = self.env.ps
        myc = self._get_class_config()
        iter_list = self.config['iter_list']
        l, b, r, t = myc['bbox_rect']
        drange = [[0, 500], [-.4, .8]]

        maxFR_all = np.empty(0)
        gridness_all = np.empty(0)

        # Separate noise sigmas
        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            if ns_idx == 0:
                mi_title = '$E-rate_{max}$ vs. gridness score'
            else:
                mi_title = None

            maxFRData = aggr.MaxPopulationFR(ps.bumpGamma[ns_idx], iter_list,
                    ignoreNaNs=True, normalizeTicks=False, collapseTrials=True)
            gridnessData = aggr.GridnessScore(ps.grids[ns_idx], iter_list,
                    normalizeTicks=False, collapseTrials=True)


            maxFRData, _, _ = maxFRData.getData()
            gridnessData, _, _ = gridnessData.getData()
            maxFR_all = np.hstack((maxFR_all, maxFRData.flatten()))
            gridness_all = np.hstack((gridness_all, gridnessData.flatten()))

            # Gamma power vs. gridness score
            fig = self._get_final_fig(myc['fig_size'])
            ax = fig.add_axes(Bbox.from_extents(l, b, r, t))
            self.plotDistribution(maxFRData, gridnessData, ax,
                                  noise_sigma=noise_sigma,
                                  range=drange,
                                  xlabel='$E-rate_{max}$',
                                  ylabel='Gridness score')
            ax.axis('tight')
            ax.xaxis.set_major_locator(ti.MultipleLocator(250))
            fname = self.config['output_dir'] + "/maxFR_gridness_probability_{0}.pdf"
            fig.savefig(fname.format(int(noise_sigma)), dpi=300,
                             transparent=True)
            plt.close(fig)

            self.mutual_information(maxFRData, gridnessData,
                                    noise_sigma=noise_sigma,
                                    title=mi_title)

        # All together
        fig = self._get_final_fig(myc['fig_size'])
        ax = fig.add_axes(Bbox.from_extents(l, b, r, t))
        self.plotDistribution(maxFR_all, gridness_all, ax,
                xlabel='$E-rate_{max}$',
                ylabel='Gridness score',
                range=drange)
        ax.axis('tight')
        ax.xaxis.set_major_locator(ti.MultipleLocator(100))
        fname = self.config['output_dir'] + "/maxFR_gridness_probability_all.pdf"
        fig.savefig(fname, dpi=300, transparent=True)
        plt.close(fig)

        self.mutual_information(maxFR_all, gridness_all)

##############################################################################
# Scatter plot of gridness score vs. max firing rate
class MaxFRGridsScatterAllPlotter(FigurePlotter):
    def __init__(self, *args, **kwargs):
        super(MaxFRGridsScatterAllPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ps = self.env.ps
        myc = self._get_class_config()
        iter_list = self.config['iter_list']
        l, b, r, t = myc['bbox_rect']
        legend_kwargs = myc['legend_kwargs']

        self.fig = self._get_final_fig(myc['fig_size'])
        self.ax = self.fig.add_axes(Bbox.from_extents(l, b, r, t))

        self.ax.hold('on')
        scatterColors = ['green', 'red', 'blue']
        scatterOrders = [3, 2, 1]

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            maxFRData = aggr.MaxPopulationFR(ps.bumpGamma[ns_idx], iter_list,
                    ignoreNaNs=True, normalizeTicks=False, collapseTrials=True)
            gridnessData = aggr.GridnessScore(ps.grids[ns_idx], iter_list,
                    normalizeTicks=False, collapseTrials=True)

            color = scatterColors[ns_idx]
            scatterPlot = scatter.ScatterPlot(
                    maxFRData, gridnessData, None, None, None, None, None,
                    c=color,
                    s=6*self.config['scale_factor'],
                    linewidth=0.3,
                    xlabel='$E-rate_{max}$',
                    ylabel='Gridness score',
                    zorder=scatterOrders[ns_idx])
            scatterPlot.plot()
        self.ax.xaxis.set_major_locator(ti.MultipleLocator(250))
        self.ax.yaxis.set_major_locator(ti.MultipleLocator(0.4))
        if self.myc['plot_legend']:
            leg = ['0', '150', '300']
            l = self.ax.legend(leg, **legend_kwargs)
            plt.setp(l.get_title(), size=legend_kwargs['fontsize'])
        #self.fig.tight_layout(**myc['tight_layout_kwargs'])

    def save(self, *args, **kwargs):
        # Linear scale
        fname = self.config['output_dir'] + "/maxFR_gridness_scatter_all.pdf"
        self.fig.savefig(fname, dpi=300, transparent=True)


##############################################################################
# Probability plots of gridness score vs seizure proportion
class PSeizureGridsProbabilityPlotter(ProbabilityPlotter):
    def __init__(self, *args, **kwargs):
        super(PSeizureGridsProbabilityPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ps = self.env.ps
        myc = self._get_class_config()
        iter_list = self.config['iter_list']
        l, b, r, t = myc['bbox_rect']
        drange = [[0, 1], [-.5, .8]]
        FRThreshold = myc['FRThreshold']
        thetaT = self.config['seizures']['thetaT']
        sig_dt = self.config['seizures']['sig_dt']

        PSeizure_all = np.empty(0)
        gridness_all = np.empty(0)

        # Separate noise sigmas
        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            if ns_idx == 0:
                mi_title = '$P(E-rate_{max} > 300)$ vs. gridness score'
            else:
                mi_title = None

            PSeizureData = aggr.MaxThetaPopulationFR(
                    thetaT, sig_dt, thresholdReduction(FRThreshold),
                    ps.bumpGamma[ns_idx], iter_list,
                    ignoreNaNs=True, normalizeTicks=True)
            gridnessData = aggr.GridnessScore(ps.grids[ns_idx], iter_list,
                    normalizeTicks=False, collapseTrials=True)


            PSeizureData, _, _ = PSeizureData.getData()
            gridnessData, _, _ = gridnessData.getData()
            PSeizure_all = np.hstack((PSeizure_all, PSeizureData.flatten()))
            gridness_all = np.hstack((gridness_all, gridnessData.flatten()))

            # Gamma power vs. gridness score
            fig = self._get_final_fig(myc['fig_size'])
            ax = fig.add_axes(Bbox.from_extents(l, b, r, t))
            self.plotDistribution(PSeizureData, gridnessData, ax,
                                  noise_sigma=noise_sigma,
                                  range=drange,
                                  xlabel='$P(E-rate_{max} > 300$',
                                  ylabel='', yticks=False)
            ax.axis('tight')
            ax.xaxis.set_major_locator(ti.MultipleLocator(.5))
            fname = self.config['output_dir'] + "/PSeizure_gridness_probability_{0}.pdf"
            fig.savefig(fname.format(int(noise_sigma)), dpi=300,
                             transparent=True)
            plt.close(fig)

            self.mutual_information(gridnessData, PSeizureData,
                                    noise_sigma=noise_sigma,
                                    title=mi_title)

        # All together
        fig = self._get_final_fig(myc['fig_size'])
        ax = fig.add_axes(Bbox.from_extents(l, b, r, t))
        self.plotDistribution(PSeizure_all, gridness_all, ax,
                xlabel='$P(E-rate_{max} > 300)$',
                ylabel='Gridness score',
                range=drange)
        ax.axis('tight')
        ax.xaxis.set_major_locator(ti.MultipleLocator(.5))
        fname = self.config['output_dir'] + "/PSeizure_gridness_probability_all.pdf"
        fig.savefig(fname, dpi=300, transparent=True)
        plt.close(fig)

        self.mutual_information(gridness_all, PSeizure_all)


##############################################################################
# P(E-rate_max) vs gridness score scatter plot
class PSeizureGridsScatterAllPlotter(FigurePlotter):
    def __init__(self, *args, **kwargs):
        super(PSeizureGridsScatterAllPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ps = self.env.ps
        myc = self._get_class_config()
        iter_list = self.config['iter_list']
        l, b, r, t = myc['bbox_rect']
        legend_kwargs = myc['legend_kwargs']
        FRThreshold = myc['FRThreshold']
        thetaT = self.config['seizures']['thetaT']
        sig_dt = self.config['seizures']['sig_dt']

        self.fig = self._get_final_fig(myc['fig_size'])
        self.ax = self.fig.add_axes(Bbox.from_extents(l, b, r, t))

        self.ax.hold('on')
        scatterColors = ['green', 'red', 'blue']
        scatterOrders = [3, 2, 1]

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            PSeizureData = aggr.MaxThetaPopulationFR(
                    thetaT, sig_dt, thresholdReduction(FRThreshold),
                    ps.bumpGamma[ns_idx], iter_list,
                    ignoreNaNs=True, normalizeTicks=True)
            gridnessData = aggr.GridnessScore(ps.grids[ns_idx], iter_list,
                    normalizeTicks=False, collapseTrials=True)

            color = scatterColors[ns_idx]
            scatterPlot = scatter.ScatterPlot(
                    PSeizureData, gridnessData, None, None, None, None, None,
                    c=color,
                    s=6*self.config['scale_factor'],
                    linewidth=0.3,
                    xlabel='$P(E-rate_{max} > 300)$',
                    ylabel='',
                    zorder=scatterOrders[ns_idx])
            scatterPlot.plot()
        self.ax.xaxis.set_major_locator(ti.MultipleLocator(.5))
        self.ax.yaxis.set_major_locator(ti.MultipleLocator(0.4))
        self.ax.yaxis.set_ticklabels([])
        leg = ['0', '150', '300']
        l = self.ax.legend(leg, **legend_kwargs)
        plt.setp(l.get_title(), size=legend_kwargs['fontsize'])
        #self.fig.tight_layout(**myc['tight_layout_kwargs'])

    def save(self, *args, **kwargs):
        # Linear scale
        fname = self.config['output_dir'] + "/PSeizure_gridness_scatter_all.pdf"
        self.fig.savefig(fname, dpi=300, transparent=True)



##############################################################################
# Seizure measure - max firing rate for the whole simulation
class RasterExamplePlotter(SweepPlotter):
    dt = .1  # ms
    freq = 8. # Hz
    const = .4 # Fraction of max. theta

    def __init__(self, *args, **kwargs):
        super(RasterExamplePlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        myc= self._get_class_config()
        sweepc = self._get_sweep_config()
        ps = self.env.ps
        iter_list = self.config['iter_list']
        sl, sb, sr, st = self.myc['sweep_rect']
        tLimits = [2e3, 4e3]

        ann150_0 = dict(
                txt='',
                rc=(4, 4),
                xytext_offset=(1, 1.5),
                color='white')
        ann0_0 = dict(
                txt='',
                rc=(5, 15),
                xytext_offset=(1.5, 1),
                color='white')
        ann0_1 = dict(
                txt='',
                rc=(20, 15),
                xytext_offset=(1.5, 1),
                color='white')
        ann300_0 = dict(
                txt='',
                rc=(5, 15),
                xytext_offset=(1.5, 1),
                color='white')
        ann0_2 = dict(
                txt='',
                rc=(15, 5),
                xytext_offset=(.5, 2),
                color='black')
        ann150_1 = dict(
                txt='',
                rc=(5, 15),
                xytext_offset=(1.5, 1),
                color='white')
        ann150_2 = dict(
                txt='',
                rc=(15, 5),
                xytext_offset=(1.5, 1),
                color='white')
        ann300_1 = dict(
                txt='',
                rc=(15, 5),
                xytext_offset=(1.5, -1),
                color='white')


        ann0   = [  ann0_0]
        ann150 = [ann150_1]
        ann300 = [ann300_0]
        ann = [ann0, ann150, ann300]
        labels = ['A', 'B', 'C']

        FRThreshold = myc['FRThreshold']
        thetaT = self.config['seizures']['thetaT']
        sig_dt = self.config['seizures']['sig_dt']

        vmin = 0.
        vmax = 500.

        saver = self.myc['fig_saver']
        saver.set_file_name(self.get_fname('raster_examples'))
        saver.ext = "pdf"
        saver.set_backend_params(dpi=300, transparent=True)

        label_it = 0
        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            ann_noise = ann[ns_idx]
            #data = aggr.MaxThetaPopulationFR(
            #        thetaT, sig_dt, thresholdReduction(FRThreshold),
            #        ps.bumpGamma[ns_idx], iter_list,
            #        ignoreNaNs=True, normalizeTicks=True)
            data = aggr.MaxThetaPopulationFR(
                    thetaT, sig_dt, np.mean,
                    ps.bumpGamma[ns_idx], iter_list,
                    ignoreNaNs=True, normalizeTicks=True)
            for annotation in ann_noise:
                r, c = annotation['rc']
                if not self.myc['plot_ann_txt']:
                    annotation['txt'] = ''
                fig = self._get_final_fig(self.myc['fig_size'])

                # Sweep
                ax_sweep = fig.add_axes(Bbox.from_extents(sl, sb, sr, st))
                _, _, cax = sweeps.plotSweep(data,
                        noise_sigma=noise_sigma,
                        ax=ax_sweep,
                        cbar=True, cbar_kw=myc['cbar_kw'],
                        vmin=vmin, vmax=vmax,
                        annotations=[annotation])
                fig.text(0.01, st, labels[label_it], size=16, weight='bold',
                         va='bottom', ha='left')

                gs = gridspec.GridSpec(3, 1, height_ratios=(2.5, 1, 1))

                # EI Raster
                ax_raster = fig.add_subplot(gs[0, 0])
                rasters.EIRaster(ps.bumpGamma[ns_idx],
                        noise_sigma=noise_sigma,
                        spaceType='bump',
                        r=r, c=c,
                        ylabelPos=self.myc['ylabelPos'],
                        tLimits=tLimits,
                        markersize=self.config['scale_factor']*self.myc['markersize'],
                        yticks=True,
                        sigmaTitle=False,
                        ann_EI=True,
                        scaleBar=125, scaleX=.85, scaleY=-.05,
                        scaleTextYOffset=.03, scaleHeight=.01,
                        rasterized=False)



                # EI rates
                ax_erates = fig.add_subplot(gs[1, 0])
                ax_irates = fig.add_subplot(gs[2, 0])
                rasters.plotAvgFiringRate(ps.bumpGamma[ns_idx],
                        spaceType='bump',
                        noise_sigma=noise_sigma,
                        popType='E',
                        r=r, c=c,
                        ylabelPos=self.myc['ylabelPos'],
                        color='red',
                        tLimits=tLimits,
                        ax=ax_erates)
                rasters.plotAvgFiringRate(ps.bumpGamma[ns_idx],
                        spaceType='bump',
                        noise_sigma=noise_sigma,
                        popType='I',
                        r=r, c=c,
                        ylabelPos=self.myc['ylabelPos'],
                        color='blue',
                        tLimits=tLimits,
                        ax=ax_irates)

                gsl = .12
                gsb = .05
                gsr = .95
                gst = .65
                #fig.text(0.01, gst, 'B', size=16, weight='bold',
                #         va='bottom', ha='left')
                gs.update(left=gsl, bottom=gsb, right=gsr, top=gst, hspace=.2)

                ax_theta = fig.add_axes(Bbox.from_extents(gsl, gst - .015,
                                                          gsr, gst + .01))
                t = np.arange(tLimits[0], tLimits[1]+self.dt, self.dt)
                theta = self.const + .5 * (1. +
                        np.cos(2*np.pi*self.freq*1e-3*t - np.pi)) * (1 - self.const)
                ax_theta.fill_between(t, theta, edgecolor='None',
                                      color=self.myc['theta_color'])
                ax_theta.set_xlim([tLimits[0], tLimits[1]])
                ax_theta.set_ylim(-.02, 1.02)
                ax_theta.axis('off')

                saver.savefig(fig)
                plt.close(fig)
                label_it += 1
        saver.close()
