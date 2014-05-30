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
]

##############################################################################
rasterRC      = [(5, 15), (5, 15), (5, 15)] # (row, col)
tLimits = [2e3, 2.25e3] # ms

rasterFigSize = (3, 1.9)
transparent   = True
rasterLeft    = 0.28
rasterBottom  = 0.1
rasterRight   = 0.95
rasterTop     = 0.8

ylabelPos   = -0.35


class EIRasterPlotter(FigurePlotter):
    def __init__(self, *args, **kwargs):
        super(EIRasterPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        output_dir = self.config['output_dir']
        ps = self.env.ps
        
        # noise_sigma = 0 pA
        fig = self._get_final_fig(rasterFigSize)
        ax = fig.add_axes(Bbox.from_extents(rasterLeft, rasterBottom, rasterRight,
            rasterTop))
        rasters.EIRaster(ps.bumpGamma[0], 
                noise_sigma=ps.noise_sigmas[0],
                spaceType='bump',
                r=rasterRC[0][0], c=rasterRC[0][1],
                ylabelPos=ylabelPos,
                tLimits=tLimits,
                ann_EI=True)
        fname = output_dir + "/bumps_raster0.png"
        fig.savefig(fname, dpi=300, transparent=transparent)
        plt.close()
            

        # noise_sigma = 150 pA
        fig = self._get_final_fig(rasterFigSize)
        ax = fig.add_axes(Bbox.from_extents(rasterLeft, rasterBottom, rasterRight,
            rasterTop))
        rasters.EIRaster(ps.bumpGamma[1], 
                noise_sigma=ps.noise_sigmas[1],
                spaceType='bump',
                r=rasterRC[1][0], c=rasterRC[1][1],
                ylabelPos=ylabelPos,
                tLimits=tLimits,
                ylabel='', yticks=False)
        fname = output_dir + "/bumps_raster150.png"
        fig.savefig(fname, dpi=300, transparent=transparent)
        plt.close()
            

        # noise_sigma = 300 pA
        fig = self._get_final_fig(rasterFigSize)
        ax = fig.add_axes(Bbox.from_extents(rasterLeft, rasterBottom, rasterRight,
            rasterTop))
        rasters.EIRaster(ps.bumpGamma[2], 
                noise_sigma=ps.noise_sigmas[2],
                spaceType='bump',
                r=rasterRC[2][0], c=rasterRC[2][1],
                ylabelPos=ylabelPos,
                tLimits=tLimits,
                ylabel='', yticks=False)
        fname = output_dir + "/bumps_raster300.png"
        fig.savefig(fname, dpi=300, transparent=transparent)
        plt.close()
        

##############################################################################
rateFigSize   = (rasterFigSize[0], 0.65)
rateLeft    = rasterLeft
rateBottom  = 0.2
rateRight   = rasterRight
rateTop     = 0.9


class EIRatePlotter(FigurePlotter):
    def __init__(self, *args, **kwargs):
        super(EIRatePlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ps = self.env.ps
        output_dir = self.config['output_dir']
        for idx, noise_sigma in enumerate(ps.noise_sigmas):
            # E cells
            fig = plt.figure(figsize=rateFigSize)
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
                    ylabelPos=ylabelPos,
                    color='red',
                    tLimits=tLimits,
                    ax=ax, **kw)
            fname = output_dir + "/bumps_rate_e{0}.pdf".format(noise_sigma)
            fig.savefig(fname, dpi=300, transparent=transparent)
            plt.close()

            # I cells
            fig = plt.figure(figsize=rateFigSize)
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
                    ylabelPos=ylabelPos,
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


###############################################################################
## Seizure measure - max firing rate per theta cycle
#maxThetaFR_cbar_kw = dict(
#        label       = "max(E rate)/$\\theta$ cycle (Hz)",
#        location    = 'left',
#        shrink      = 0.8,
#        pad         = 0.25,
#        ticks       = ti.MultipleLocator(100),
#        #ticks       = ti.LogLocator(base=4),
#        #format      = ti.LogFormatter(4),
#        rasterized  = True)
#maxThetaFR_vmin = 2.
#maxThetaFR_vmax = 500.
#
#thetaT = 125.  # ms
#sig_dt = .5    # ms
#
## mean
#if args.maxThetaFRSweeps or args.all:
#    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
#        fig, ax = ds.getDefaultSweepFig(scale=1., colorBarPos='left')
#        kw = dict(cbar=False)
#        if ns_idx != 0:
#            kw['ylabel'] = ''
#            kw['yticks'] = False
#        if ns_idx == 0:
#            kw['cbar'] = True
#        data = aggr.MaxThetaPopulationFR(
#                thetaT, sig_dt, np.mean,
#                ps.bumpGamma[ns_idx], ds.iterList,
#                ignoreNaNs=True, normalizeTicks=True)
#        _, _, cax = sweeps.plotSweep(data,
#                noise_sigma=noise_sigma,
#                ax=ax,
#                cbar_kw=maxThetaFR_cbar_kw,
#                vmin=maxThetaFR_vmin, vmax=maxThetaFR_vmax,
#                #norm=LogNorm(maxThetaFR_vmin, maxThetaFR_vmax),
#                sigmaTitle=False,
#                **kw)
#        fname = outputDir + "/bumps_popMaxThetaFR_sweep{0}.pdf"
#        fig.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)
#        #plt.close()
#
#
## median
#if args.maxThetaFRSweeps_median or args.all:
#    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
#        fig, ax = ds.getDefaultSweepFig(scale=1., colorBarPos='left')
#        kw = dict(cbar=False)
#        if ns_idx != 0:
#            kw['ylabel'] = ''
#            kw['yticks'] = False
#        if ns_idx == 0:
#            kw['cbar'] = True
#        data = aggr.MaxThetaPopulationFR(
#                thetaT, sig_dt, np.median,
#                ps.bumpGamma[ns_idx], ds.iterList,
#                ignoreNaNs=True, normalizeTicks=True)
#        _, _, cax = sweeps.plotSweep(data,
#                noise_sigma=noise_sigma,
#                ax=ax,
#                cbar_kw=maxThetaFR_cbar_kw,
#                vmin=maxThetaFR_vmin, vmax=maxThetaFR_vmax,
#                sigmaTitle=False,
#                **kw)
#        fname = outputDir + "/bumps_popMaxThetaFR_median_sweep{0}.pdf"
#        fig.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)
#        plt.close()
#
#
## std
#maxThetaFR_std_cbar_kw = dict(
#        label       = "max(E rate)/$\\theta$ cycle (Hz)",
#        location    = 'left',
#        shrink      = 0.8,
#        pad         = 0.25,
#        ticks       = ti.MaxNLocator(4),
#        rasterized  = True)
#maxThetaFR_std_vmin = 0.
#maxThetaFR_std_vmax = None
#if args.maxThetaFRSweeps_std or args.all:
#    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
#        fig, ax = ds.getDefaultSweepFig(scale=1., colorBarPos='left')
#        kw = dict(cbar=True)
#        data = aggr.MaxThetaPopulationFR(
#                thetaT, sig_dt, np.std,
#                ps.bumpGamma[ns_idx], ds.iterList,
#                ignoreNaNs=True, normalizeTicks=True)
#        _, _, cax = sweeps.plotSweep(data,
#                noise_sigma=noise_sigma,
#                ax=ax,
#                cbar_kw=maxThetaFR_std_cbar_kw,
#                vmin=maxThetaFR_std_vmin, vmax=maxThetaFR_std_vmax,
#                sigmaTitle=False,
#                **kw)
#        fname = outputDir + "/bumps_popMaxThetaFR_std_sweep{0}.pdf"
#        fig.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)
#        plt.close()
#
#
## Proportion of cycles with max firing rate larger than threshold, i.e. number
## of seizures during the simulation
#FRThreshold = 300
#
#maxThetaFR_std_cbar_kw = dict(
#        label       = "P[max(rate during $\\theta$) > {0}]".format(FRThreshold),
#        location    = 'left',
#        shrink      = 0.8,
#        pad         = 0.25,
#        ticks       = ti.MultipleLocator(0.5),
#        rasterized  = True)
#seizureThreshold_vmin = 0.
#seizureThreshold_vmax = 1.
#
#class thresholdReduction(object):
#    def __init__(self, threshold):
#        self.threshold = threshold
#
#    def __call__(self, data):
#        return float(np.count_nonzero(data >= self.threshold)) / len(data)
#
#if args.seizureProportion or args.all:
#    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
#        fig, ax = ds.getDefaultSweepFig(scale=1., colorBarPos='left')
#        kw = dict(cbar=False)
#        if ns_idx != 0:
#            kw['ylabel'] = ''
#            kw['yticks'] = False
#        if ns_idx == 0:
#            kw['cbar'] = True
#        data = aggr.MaxThetaPopulationFR(
#                thetaT, sig_dt, thresholdReduction(FRThreshold),
#                ps.bumpGamma[ns_idx], ds.iterList,
#                ignoreNaNs=True, normalizeTicks=True)
#        _, _, cax = sweeps.plotSweep(data,
#                noise_sigma=noise_sigma,
#                ax=ax,
#                cbar_kw=maxThetaFR_std_cbar_kw,
#                vmin=maxThetaFR_std_vmin, vmax=maxThetaFR_std_vmax,
#                sigmaTitle=False,
#                **kw)
#        fname = outputDir + "/bumps_seizureProportion_sweep{0}.pdf"
#        fig.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)
#        plt.close()
#
#
#
#
###############################################################################
##       Histograms of maxima of firing rates during theta cycles
#if args.maxThetaFRHist or args.all:
#    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
#        fig, _ = ds.getDefaultSweepFig(scale=1., colorBarPos='left')
#        ax = plt.subplot(111)
#        kw = dict(cbar=True)
#        data = aggr.MaxThetaPopulationFR(
#                thetaT, sig_dt, None,
#                ps.bumpGamma[ns_idx], ds.iterList,
#                ignoreNaNs=True, normalizeTicks=True)
#        flatData = np.hstack(data.getNonReducedData().flatten().tolist())
#        flatData = flatData[np.logical_not(np.isnan(flatData))]
#        base.plotOneHist(flatData, bins=80, ax=ax, rwidth=.8, normed=True)
#        ax.set_xlabel('Max rate in $\\theta$ cycle (Hz)')
#        ax.set_ylabel('p(rate)')
#        fig.tight_layout()
#        fname = outputDir + "/bumps_popMaxThetaFR_hist{0}.pdf"
#        fig.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)
#        plt.close()
