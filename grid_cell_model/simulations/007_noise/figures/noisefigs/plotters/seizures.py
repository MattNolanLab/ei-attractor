'''
Figure illustrating seizures.
'''
from __future__ import absolute_import, print_function

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ti
from matplotlib.colors import LogNorm
from matplotlib.transforms import Bbox

from ..EI_plotting import sweeps, rasters, base
from ..EI_plotting import aggregate as aggr
from .base import FigurePlotter, SweepPlotter

__all__ = [
    'EIRasterPlotter',
    'EIRatePlotter',
    'MaxPopulationFRSweepsPlotter',
    'MaxMeanThetaFRSweepPlotter',
    'MaxStdThetaFRSweepPlotter',
    'MaxMedianThetaFRSweepPlotter',
    'MaxThetaFRHistPlotter',
    'PSeizureSweepPlotter',
]

##############################################################################
rasterRC      = [(5, 15), (5, 15), (5, 15)] # (row, col)
tLimits = [2e3, 2.25e3] # ms

transparent   = True
rasterLeft    = 0.28
rasterBottom  = 0.1
rasterRight   = 0.95
rasterTop     = 0.8
        

class EIRasterPlotter(FigurePlotter):
    def __init__(self, *args, **kwargs):
        super(EIRasterPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ps = self.env.ps

        output_dir = self.config['output_dir']
        
        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            fig = self._get_final_fig(self.myc['fig_size'])
            ax = fig.add_axes(Bbox.from_extents(rasterLeft, rasterBottom, rasterRight,
                rasterTop))
            rasters.EIRaster(ps.bumpGamma[ns_idx], 
                    noise_sigma=noise_sigma,
                    spaceType='bump',
                    r=rasterRC[ns_idx][0], c=rasterRC[ns_idx][1],
                    ylabelPos=self.myc['ylabelPos'],
                    tLimits=tLimits,
                    markersize=self.config['scale_factor'],
                    ylabel='' if self.myc['yticks'][ns_idx] == False else None,
                    yticks=self.myc['yticks'][ns_idx],
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
                kw = dict(cbar=False)
                if ns_idx != 0:
                    kw['ylabel'] = ''
                    kw['yticks'] = False
                if ns_idx == 0:
                    kw['cbar'] = True
                data = aggr.MaxPopulationFR(ps.bumpGamma[ns_idx], iter_list,
                        ignoreNaNs=True, normalizeTicks=True)
                _, _, cax = sweeps.plotSweep(data,
                        noise_sigma=noise_sigma,
                        ax=ax,
                        cbar_kw=myc['cbar_kw'],
                        vmin=maxFR_vmin, vmax=maxFR_vmax,
                        **kw)


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
