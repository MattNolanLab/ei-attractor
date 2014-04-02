#!/usr/bin/env python
#
'''
Figure illustrating seizures.
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ti
from matplotlib.transforms import Bbox

import default_settings as ds
from EI_plotting          import rasters
from submitting import flagparse


outputDir = ds.figOutputDir

parser = flagparse.FlagParser()
parser.add_flag('--rastersFlag')
parser.add_flag('--rates')
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




