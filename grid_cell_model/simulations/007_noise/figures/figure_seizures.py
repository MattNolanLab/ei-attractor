#!/usr/bin/env python
#
'''
Figure illustrating seizures.
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ti
from matplotlib.colors import LogNorm
from matplotlib.transforms import Bbox

import default_settings as ds
from EI_plotting import sweeps, rasters, base
from EI_plotting import aggregate as aggr
from submitting  import flagparse


outputDir = ds.figOutputDir

parser = flagparse.FlagParser()
parser.add_flag('--rastersFlag')
parser.add_flag('--rates')
parser.add_flag('--maxFRSweeps')
parser.add_flag('--maxThetaFRSweeps')
parser.add_flag('--maxThetaFRSweeps_median')
parser.add_flag('--maxThetaFRSweeps_std')
parser.add_flag('--maxThetaFRHist')
args = parser.parse_args()

###############################################################################
ps = ds.getDefaultParamSpaces()

sweepFigSize = (3.7, 2.6)
sweepLeft   = 0.08
sweepBottom = 0.2
sweepRight  = 0.8
sweepTop    = 0.85
transparent  = True



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


if args.rastersFlag or args.all:
    # noise_sigma = 0 pA
    fig = plt.figure("rasters0", figsize=rasterFigSize)
    ax = fig.add_axes(Bbox.from_extents(rasterLeft, rasterBottom, rasterRight,
        rasterTop))
    rasters.EIRaster(ps.bumpGamma[0], 
            noise_sigma=ps.noise_sigmas[0],
            spaceType='bump',
            r=rasterRC[0][0], c=rasterRC[0][1],
            ylabelPos=ylabelPos,
            tLimits=tLimits,
            ann_EI=True)
    fname = outputDir + "/bumps_raster0.png"
    fig.savefig(fname, dpi=300, transparent=transparent)
    plt.close()
        

    # noise_sigma = 150 pA
    fig = plt.figure("rasters150", figsize=rasterFigSize)
    ax = fig.add_axes(Bbox.from_extents(rasterLeft, rasterBottom, rasterRight,
        rasterTop))
    rasters.EIRaster(ps.bumpGamma[1], 
            noise_sigma=ps.noise_sigmas[1],
            spaceType='bump',
            r=rasterRC[1][0], c=rasterRC[1][1],
            ylabelPos=ylabelPos,
            tLimits=tLimits,
            ylabel='', yticks=False)
    fname = outputDir + "/bumps_raster150.png"
    fig.savefig(fname, dpi=300, transparent=transparent)
    plt.close()
        

    # noise_sigma = 300 pA
    fig = plt.figure("rasters300", figsize=rasterFigSize)
    ax = fig.add_axes(Bbox.from_extents(rasterLeft, rasterBottom, rasterRight,
        rasterTop))
    rasters.EIRaster(ps.bumpGamma[2], 
            noise_sigma=ps.noise_sigmas[2],
            spaceType='bump',
            r=rasterRC[2][0], c=rasterRC[2][1],
            ylabelPos=ylabelPos,
            tLimits=tLimits,
            ylabel='', yticks=False)
    fname = outputDir + "/bumps_raster300.png"
    fig.savefig(fname, dpi=300, transparent=transparent)
    plt.close()
        

##############################################################################
rateFigSize   = (rasterFigSize[0], 0.65)
rateLeft    = rasterLeft
rateBottom  = 0.2
rateRight   = rasterRight
rateTop     = 0.9


if args.rates or args.all:
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
        fname = outputDir + "/bumps_rate_e{0}.pdf".format(noise_sigma)
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
        fname = outputDir + "/bumps_rate_i{0}.pdf".format(noise_sigma)
        fig.savefig(fname, dpi=300, transparent=transparent)
        plt.close()




##############################################################################
# Seizure measure - max firing rate for the whole simulation
maxFR_cbar_kw = dict(
        label       = "max(E rate) (Hz)",
        location    = 'left',
        shrink      = 0.8,
        pad         = 0.25,
        ticks       = ti.MultipleLocator(100),
        rasterized  = True)
maxFR_vmin = 0
maxFR_vmax = 500.

if args.maxFRSweeps or args.all:
    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        fig, ax = ds.getDefaultSweepFig(scale=1., colorBarPos='left')
        kw = dict(cbar=False)
        if ns_idx != 0:
            kw['ylabel'] = ''
            kw['yticks'] = False
        if ns_idx == 0:
            kw['cbar'] = True
        data = aggr.MaxPopulationFR(ps.bumpGamma[ns_idx], ds.iterList,
                ignoreNaNs=True, normalizeTicks=True)
        _, _, cax = sweeps.plotSweep(data,
                noise_sigma=noise_sigma,
                xlabel='', xticks=False,
                ax=ax,
                cbar_kw=maxFR_cbar_kw,
                vmin=maxFR_vmin, vmax=maxFR_vmax,
                **kw)
        fname = outputDir + "/bumps_popMaxFR_sweep{0}.pdf"
        fig.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)
        plt.close()



##############################################################################
# Seizure measure - max firing rate per theta cycle (mean)
maxThetaFR_cbar_kw = dict(
        label       = "max(E rate)/$\\theta$ cycle (Hz)",
        location    = 'left',
        shrink      = 0.8,
        pad         = 0.25,
        ticks       = ti.MultipleLocator(100),
        #ticks       = ti.LogLocator(base=4),
        #format      = ti.LogFormatter(4),
        rasterized  = True)
maxThetaFR_vmin = 2.
maxThetaFR_vmax = 500.

thetaT = 125.  # ms
sig_dt = .5    # ms

if args.maxThetaFRSweeps or args.all:
    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        fig, ax = ds.getDefaultSweepFig(scale=1., colorBarPos='left')
        kw = dict(cbar=False)
        if ns_idx != 0:
            kw['ylabel'] = ''
            kw['yticks'] = False
        if ns_idx == 0:
            kw['cbar'] = True
        data = aggr.MaxThetaPopulationFR(
                thetaT, sig_dt, np.mean,
                ps.bumpGamma[ns_idx], ds.iterList,
                ignoreNaNs=True, normalizeTicks=True)
        _, _, cax = sweeps.plotSweep(data,
                noise_sigma=noise_sigma,
                ax=ax,
                cbar_kw=maxThetaFR_cbar_kw,
                vmin=maxThetaFR_vmin, vmax=maxThetaFR_vmax,
                #norm=LogNorm(maxThetaFR_vmin, maxThetaFR_vmax),
                sigmaTitle=False,
                **kw)
        fname = outputDir + "/bumps_popMaxThetaFR_sweep{0}.pdf"
        fig.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)
        #plt.close()


if args.maxThetaFRSweeps_median or args.all:
    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        fig, ax = ds.getDefaultSweepFig(scale=1., colorBarPos='left')
        kw = dict(cbar=False)
        if ns_idx != 0:
            kw['ylabel'] = ''
            kw['yticks'] = False
        if ns_idx == 0:
            kw['cbar'] = True
        data = aggr.MaxThetaPopulationFR(
                thetaT, sig_dt, np.median,
                ps.bumpGamma[ns_idx], ds.iterList,
                ignoreNaNs=True, normalizeTicks=True)
        _, _, cax = sweeps.plotSweep(data,
                noise_sigma=noise_sigma,
                ax=ax,
                cbar_kw=maxThetaFR_cbar_kw,
                vmin=maxThetaFR_vmin, vmax=maxThetaFR_vmax,
                sigmaTitle=False,
                **kw)
        fname = outputDir + "/bumps_popMaxThetaFR_median_sweep{0}.pdf"
        fig.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)
        plt.close()


maxThetaFR_std_cbar_kw = dict(
        label       = "max(E rate)/$\\theta$ cycle (Hz)",
        location    = 'left',
        shrink      = 0.8,
        pad         = 0.25,
        ticks       = ti.MaxNLocator(4),
        rasterized  = True)
maxThetaFR_std_vmin = 0.
maxThetaFR_std_vmax = None
if args.maxThetaFRSweeps_std or args.all:
    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        fig, ax = ds.getDefaultSweepFig(scale=1., colorBarPos='left')
        kw = dict(cbar=True)
        data = aggr.MaxThetaPopulationFR(
                thetaT, sig_dt, np.std,
                ps.bumpGamma[ns_idx], ds.iterList,
                ignoreNaNs=True, normalizeTicks=True)
        _, _, cax = sweeps.plotSweep(data,
                noise_sigma=noise_sigma,
                ax=ax,
                cbar_kw=maxThetaFR_std_cbar_kw,
                vmin=maxThetaFR_std_vmin, vmax=maxThetaFR_std_vmax,
                sigmaTitle=False,
                **kw)
        fname = outputDir + "/bumps_popMaxThetaFR_std_sweep{0}.pdf"
        fig.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)
        plt.close()

# Histograms of maxima of firing rates during theta cycles
if args.maxThetaFRHist or args.all:
    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        fig, _ = ds.getDefaultSweepFig(scale=1., colorBarPos='left')
        ax = plt.subplot(111)
        kw = dict(cbar=True)
        data = aggr.MaxThetaPopulationFR(
                thetaT, sig_dt, None,
                ps.bumpGamma[ns_idx], ds.iterList,
                ignoreNaNs=True, normalizeTicks=True)
        flatData = np.hstack(data.getNonReducedData().flatten().tolist())
        flatData = flatData[np.logical_not(np.isnan(flatData))]
        base.plotOneHist(flatData, bins=80, ax=ax, rwidth=.8, normed=True)
        ax.set_xlabel('Max rate in $\\theta$ cycle (Hz)')
        ax.set_ylabel('p(rate)')
        fig.tight_layout()
        fname = outputDir + "/bumps_popMaxThetaFR_hist{0}.pdf"
        fig.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)
        plt.close()
