#!/usr/bin/env python
#
'''
Everything related to velocity: sweeps, lines, etc.
'''
import logging

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ti
from matplotlib.transforms import Bbox

from grid_cell_model.plotting.global_defs import globalAxesSettings
from grid_cell_model.parameters import JobTrialSpace2D
from grid_cell_model.submitting import flagparse

from ..EI_plotting import sweeps, details, rasters, scatter
from ..EI_plotting import aggregate as aggr
from .base import FigurePlotter, SweepPlotter

logger = logging.getLogger(__name__)

__all__ = [
    'GridsLineErrScatterPlotter',
    'GridsLineSlopeScatterPlotter',
    'LineErrSlopeScatterPlotter',
]

###############################################################################

#rasterRC      = [(5, 15), (5, 15), (5, 15)] # (row, col)
#
#parser = flagparse.FlagParser()
#parser.add_flag('--velLines')
#parser.add_flag('--detailed_noise')
#parser.add_flag('--slope_sweeps')
#parser.add_flag('--rasters')
#parser.add_flag('--rastersZoom')
#parser.add_flag('--rates')
#args = parser.parse_args()
#
#ps = ds.getDefaultParamSpaces()
#
################################################################################
#
#def plotSlopes(ax, dataSpace, pos, noise_sigma, **kw):
#    # kwargs
#    trialNum   = kw.pop('trialNum', 0)
#    markersize = kw.pop('markersize', 4)
#    color      = kw.pop('color', 'blue')
#    xlabel     = kw.pop('xlabel', 'Velocity current (pA)')
#    ylabel     = kw.pop('ylabel', 'Bump speed\n(neurons/s)')
#    xticks     = kw.pop('xticks', True)
#    yticks     = kw.pop('yticks', True)
#    g_ann      = kw.pop('g_ann', True)
#    sigma_ann  = kw.pop('sigma_ann', True)
#    kw['markeredgecolor'] = color
#
#    r = pos[0]
#    c = pos[1]
#    d = dataSpace[r][c].getAllTrialsAsDataSet().data
#    a = d['analysis']
#    IvelVec = dataSpace[r][c][trialNum].data['IvelVec']
#    slopes = a['bumpVelAll']
#    lineFit = a['lineFitLine']
#    fitIvelVec = a['fitIvelVec']
#
#    nTrials = slopes.shape[0]
#    avgSlope = np.mean(slopes, axis=0)
#    stdSlope = np.std(slopes, axis=0)
#
#    if (ax is None):
#        ax = plt.gca()
#    plt.hold('on')
#    globalAxesSettings(ax)
#
#    ax.plot(IvelVec, slopes.T, 'o', markerfacecolor='none', markersize=markersize, **kw)
#    #ax.errorbar(IvelVec, avgSlope, stdSlope, fmt='o-',
#    #        markersize=markersize, color=color, alpha=0.5, **kw)
#    ax.plot(fitIvelVec, lineFit, '-', linewidth=1, color=color, **kw)
#
#    ax.set_xlabel(xlabel)
#    ax.set_ylabel(ylabel)
#    ax.spines['top'].set_visible(False)
#    ax.spines['right'].set_visible(False)
#    ax.xaxis.set_major_locator(ti.MultipleLocator(50))
#    ax.xaxis.set_minor_locator(ti.AutoMinorLocator(5))
#    ax.yaxis.set_major_locator(ti.MultipleLocator(20))
#    ax.yaxis.set_minor_locator(ti.MultipleLocator(10))
#    ax.margins(0.05)
#
#    if (not xticks):
#        ax.xaxis.set_ticklabels([])
#    if (not yticks):
#        ax.yaxis.set_ticklabels([])
#
#    # Annotations
#    if (sigma_ann):
#        sigma_txt = '$\sigma$ = {0} pA'.format(noise_sigma)
#        ax.set_title(sigma_txt, y=1.1, va='bottom')
#
#    if (g_ann):
#        Y, X = aggr.computeVelYX(dataSpace, ds.iterList, r, c)
#        gE = Y[r, c]
#        gI = X[r, c]
#        g_txt = '$g_E$ = {0}\n$g_I$ = {1} nS'.format(gE, gI)
#    else:
#        g_txt = ''
#
#    txt = '{0}'.format(g_txt)
#    ax.text(0.05, 1.1, txt, transform=ax.transAxes, va='top',
#            ha='left', size='x-small')
#
#
###############################################################################
#velFigsize =(2.6, 2)
#velLeft    = 0.3
#velBottom  = 0.35
#velRight   = 0.95
#velTop     = 0.65
#
#if args.velLines or args.all:
#    #positions = ((4, 27), (4, 27), (4, 27))
#    positions = ((5, 15), (5, 15), (5, 15))
#    fig = plt.figure(figsize=(2.5, velFigsize[1]))
#    ax = fig.add_axes(Bbox.from_extents(0.3, velBottom, velRight,
#        velTop))
#    plotSlopes(ax, ps.v[0], positions[0], noise_sigma=ps.noise_sigmas[0],
#            xlabel='',
#            color='blue')
#    fname = outputDir + "/velocity_slope_examples_0.pdf"
#    plt.savefig(fname, dpi=300, transparent=True)
#
#    fig = plt.figure(figsize=(2.5, velFigsize[1]))
#    ax = fig.add_axes(Bbox.from_extents(0.3, velBottom, velRight,
#        velTop))
#    plotSlopes(ax, ps.v[1], positions[1], noise_sigma=ps.noise_sigmas[1],
#            ylabel='',
#            g_ann=False,
#            color='green')
#    fname = outputDir + "/velocity_slope_examples_1.pdf"
#    plt.savefig(fname, dpi=300, transparent=True)
#
#    fig = plt.figure(figsize=(2.5, velFigsize[1]))
#    ax = fig.add_axes(Bbox.from_extents(0.3, velBottom, velRight,
#        velTop))
#    plotSlopes(ax, ps.v[2], positions[2], noise_sigma=ps.noise_sigmas[2],
#            xlabel='',
#            ylabel='',
#            g_ann=False,
#            color='red')
#    fname = outputDir + "/velocity_slope_examples_2.pdf"
#    plt.savefig(fname, dpi=300, transparent=True)
#
#
###############################################################################
#EI13Root  = 'simulation_data/submission/detailed_noise/velocity/EI-1_3'
#EI31Root  = 'simulation_data/submission/detailed_noise/velocity/EI-3_1'
#detailedShape = (31, 9)
#
#EI13PS = JobTrialSpace2D(detailedShape, EI13Root)
#EI31PS = JobTrialSpace2D(detailedShape, EI31Root)
#detailedNTrials = None
#
#
#sliceFigSize = (3.5, 2.5)
#sliceLeft   = 0.25
#sliceBottom = 0.3
#sliceRight  = 0.95
#sliceTop    = 0.8
#if args.detailed_noise or args.all:
#    ylabelPos = -0.27
#
#    types = ('velocity', 'slope')
#    fig = plt.figure(figsize=sliceFigSize)
#    ax = fig.add_axes(Bbox.from_extents(sliceLeft, sliceBottom, sliceRight,
#        sliceTop))
#    details.plotDetailedNoise(EI13PS, detailedNTrials, types, ax=ax,
#            ylabelPos=ylabelPos,
#            xlabel='', xticks=False,
#            color='red', markerfacecolor='red', zorder=10)
#    details.plotDetailedNoise(EI31PS, detailedNTrials, types, ax=ax,
#            xlabel='', xticks=False,
#            ylabel='Slope\n(neurons/s/pA)', ylabelPos=ylabelPos,
#            color='#505050')
#    ax.xaxis.set_visible(False)
#    ax.spines['bottom'].set_visible(False)
#    ax.yaxis.set_major_locator(ti.MultipleLocator(0.4))
#    ax.yaxis.set_minor_locator(ti.MultipleLocator(0.2))
#    ax.spines['bottom'].set_visible(False)
#
#    fname = outputDir + "/velocity_detailed_noise_slope.pdf"
#    plt.savefig(fname, dpi=300, transparent=True)
#    plt.close()
#
#
#    types = ('velocity', 'fitErr')
#    fig = plt.figure(figsize=sliceFigSize)
#    ax = fig.add_axes(Bbox.from_extents(sliceLeft, sliceBottom, sliceRight,
#        sliceTop))
#    _, p13, l13 = details.plotDetailedNoise(EI13PS, detailedNTrials, types, ax=ax,
#            ylabelPos=ylabelPos,
#            color='red', markerfacecolor='red', zorder=10)
#    _, p31, l31 = details.plotDetailedNoise(EI31PS, detailedNTrials, types, ax=ax,
#            ylabel='Fit error\n(neurons/s/trial)', ylabelPos=ylabelPos,
#            color='#505050')
#    ax.yaxis.set_major_locator(ti.MultipleLocator(4))
#    ax.yaxis.set_minor_locator(ti.MultipleLocator(2))
#    leg = ['(3, 1)',  '(1, 3)']
#    l = ax.legend([p31, p13], leg, loc=(0.55, 0.3), fontsize='small', frameon=False,
#            numpoints=1, title='($g_E,\ g_I$) (nS)')
#    plt.setp(l.get_title(), fontsize='x-small')
#
#    fname = outputDir + "/velocity_detailed_noise_error.pdf"
#    plt.savefig(fname, dpi=300, transparent=True)
#    plt.close()
#
#
#
#
#
###############################################################################
##                           Velocity (slope) sweeps
#
## This should be corresponding to the velLine examples as well !!
#rasterRC      = [(5, 15), (5, 15), (5, 15)] # (row, col)
#slopeVarList = ['lineFitSlope']
#slope_vmin = None
#slope_vmax = None
#velSweep_cmap = 'jet'
#
#slope_cbar_kw= dict(
#        location='right',
#        shrink = 0.8,
#        pad = -0.05,
#        label='Slope\n(neurons/s/pA)',
#        ticks=ti.MultipleLocator(0.4),
#        extend='max', extendfrac=0.1)
#
#ann0 = dict(
#        txt='A',
#        rc=rasterRC[0],
#        xytext_offset=(1.5, 0.5),
#        color='white')
#
#ann = [ann0]
#
#def createSweepFig(name=None):
#    sweepFigSize = (3.7, 2.6)
#    sweepLeft   = 0.08
#    sweepBottom = 0.2
#    sweepRight  = 0.8
#    sweepTop    = 0.85
#    transparent  = True
#    fig = plt.figure(name, figsize=sweepFigSize)
#    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
#        sweepTop))
#    return fig, ax
#
#if args.slope_sweeps or args.all:
#    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
#        fig, ax = createSweepFig(None)
#        kw = dict(cbar=True)
#        #if (ns_idx == 2):
#        #    kw['cbar'] = True
#        if (ns_idx != 0):
#            kw['ylabel'] = ''
#            kw['yticks'] = False
#        _, ax, cax = sweeps.plotVelTrial(ps.v[ns_idx], slopeVarList, ds.iterList,
#                noise_sigma, sigmaTitle=True,
#                ax=ax,
#                cbar_kw=slope_cbar_kw,
#                cmap=velSweep_cmap,
#                vmin=slope_vmin, vmax=slope_vmax,
#                annotations=ann,
#                **kw)
#        fname = outputDir + "/velocity_slope_sweeps{}.pdf"
#        fig.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)
#
#
#
###############################################################################
##                           Raster and rate plots
###############################################################################
#tLimits  = [2e3, 3e3] # ms
#trialNum = 0
#
#rasterFigSize = (3.75, 2.2)
#transparent   = True
#rasterLeft    = 0.2
#rasterBottom  = 0.2
#rasterRight   = 0.99
#rasterTop     = 0.8
#
#ylabelPos   = -0.22
#
#
#if args.rasters or args.all:
#    logger.info("Plotting rasters")
#    for idx, noise_sigma in enumerate(ps.noise_sigmas):
#        logger.info("   Rasters: %d pA", noise_sigma)
#        fig = plt.figure(figsize=rasterFigSize)
#        ax = fig.add_axes(Bbox.from_extents(rasterLeft, rasterBottom, rasterRight,
#            rasterTop))
#        kw = dict(scaleBar=None)
#        if (idx != 0):
#            kw['ylabel'] = ''
#            kw['yticks'] = False
#        if idx == 2:
#            kw['scaleBar'] = 125
#        rasters.EIRaster(ps.v[idx], 
#                noise_sigma=noise_sigma,
#                spaceType='velocity',
#                r=rasterRC[idx][0], c=rasterRC[idx][1],
#                ylabelPos=ylabelPos,
#                tLimits=tLimits,
#                trialNum=trialNum,
#                ann_EI=True,
#                scaleX=0.85,
#                scaleY=-0.15,
#                **kw)
#        fname = outputDir + "/velocity_raster{0}.png"
#        fig.savefig(fname.format(int(noise_sigma)), dpi=300,
#                transparent=transparent)
#        plt.close()
#        
#        
#
###############################################################################
#rateFigSize   = (rasterFigSize[0], 1)
#rateLeft    = rasterLeft
#rateBottom  = 0.2
#rateRight   = rasterRight
#rateTop     = 0.7
#
#
#if args.rates or args.all:
#    for idx, noise_sigma in enumerate(ps.noise_sigmas):
#        # E cells
#        fig = plt.figure(figsize=rateFigSize)
#        ax = fig.add_axes(Bbox.from_extents(rateLeft, rateBottom, rateRight,
#            rateTop))
#        kw = {}
#        if (idx != 0):
#            kw['ylabel'] = ''
#
#        rasters.plotAvgFiringRate(ps.v[idx],
#                spaceType='velocity',
#                noise_sigma=noise_sigma,
#                popType='E',
#                r=rasterRC[idx][0], c=rasterRC[idx][1],
#                ylabelPos=ylabelPos,
#                color='red',
#                tLimits=tLimits,
#                trialNum=trialNum,
#                sigmaTitle=False,
#                ax=ax, **kw)
#        fname = outputDir + "/velocity_rate_e{0}.pdf".format(noise_sigma)
#        fig.savefig(fname, dpi=300, transparent=transparent)
#        plt.close()
#
#        # I cells
#        fig = plt.figure(figsize=rateFigSize)
#        ax = fig.add_axes(Bbox.from_extents(rateLeft, rateBottom, rateRight,
#            rateTop))
#        kw = {}
#        if (idx != 0):
#            kw['ylabel'] = ''
#
#        rasters.plotAvgFiringRate(ps.v[idx],
#                spaceType='velocity',
#                noise_sigma=noise_sigma,
#                popType='I', 
#                r=rasterRC[idx][0], c=rasterRC[idx][1],
#                ylabelPos=ylabelPos,
#                color='blue',
#                tLimits=tLimits,
#                trialNum=trialNum,
#                sigmaTitle=False,
#                ax=ax, **kw)
#        fname = outputDir + "/velocity_rate_i{0}.pdf".format(noise_sigma)
#        fig.savefig(fname, dpi=300, transparent=transparent)
#        plt.close()
#
#
#
#
###############################################################################
##                           Raster and rate zoom-ins
###############################################################################
#tLimits  = [2.75e3, 2.875e3] # ms
#trialNum = 0
#
#c = 0.75
#rasterZoomFigSize = c*rasterFigSize[0], c*rasterFigSize[1]
#ylabelPos   = -0.22
#
#
#if args.rastersZoom or args.all:
#    for idx, noise_sigma in enumerate(ps.noise_sigmas):
#        fig = plt.figure(figsize=rasterZoomFigSize)
#        ax = fig.add_axes(Bbox.from_extents(rasterLeft, rasterBottom, rasterRight,
#            rasterTop))
#        kw = dict(scaleBar=None)
#        if idx == 2:
#            kw['scaleBar'] = 25
#        rasters.EIRaster(ps.v[idx], 
#                noise_sigma=noise_sigma,
#                spaceType='velocity',
#                r=rasterRC[idx][0], c=rasterRC[idx][1],
#                ylabelPos=ylabelPos,
#                tLimits=tLimits,
#                trialNum=trialNum,
#                sigmaTitle=False,
#                ann_EI=True,
#                scaleX=0.75,
#                scaleY=-0.15,
#                ylabel='', yticks=False,
#                **kw)
#        fname = outputDir + "/velocity_raster_zooms{0}.png"
#        fig.savefig(fname.format(int(noise_sigma)), dpi=300,
#                transparent=transparent)
#        plt.close()
        

##############################################################################
# Scatter plot of gridness score vs. bump speed line fit error
class GridsLineErrScatterPlotter(FigurePlotter):
    def __init__(self, *args, **kwargs):
        super(GridsLineErrScatterPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ps = self.env.ps
        output_dir = self.config['output_dir']
        iter_list = self.config['iter_list']

        scatterFigSize = (8.27, 11.69)
        ignoreNaNs = True
        fname = output_dir + "/suppFigure_grids_vs_line_fit_err.pdf"

        gridNTrials = 3
        velNTrials = None
        typesSlope = ['velocity', 'fitErr']
        typesGrids = ['grids', 'gridnessScore']

        xlabel = 'Bump speed fit error (neurons/s/trial)'
        ylabel = 'Gridness score'

        fig = plt.figure(figsize=scatterFigSize)
        scatterPlot = scatter.FullScatterPlot(
                ps.v, ps.grids, typesSlope, typesGrids, iter_list, velNTrials,
                gridNTrials,
                s=25,
                linewidth=0.3,
                color2D=True,
                xlabel=xlabel,
                ylabel=ylabel,
                sigmaTitle=True,
                noise_sigmas=ps.noise_sigmas,
                ignoreNaNs=ignoreNaNs,
                captionLetters=('A', 'B', 'C'),
                fig=fig)
        scatterPlot.plot(captionLeft=-0.1)
        scatterPlot.set_titleSizes(16)

        fig.savefig(fname, dpi=300)


##############################################################################
# Scatter plot of gridness score vs. bump speed line slope
class GridsLineSlopeScatterPlotter(FigurePlotter):
    def __init__(self, *args, **kwargs):
        super(GridsLineSlopeScatterPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ps = self.env.ps
        output_dir = self.config['output_dir']
        iter_list = self.config['iter_list']

        scatterFigSize = (8.27, 11.69)
        ignoreNaNs = True
        fname = output_dir + "/suppFigure_grids_vs_line_slope.pdf"

        gridNTrials = 3
        velNTrials = None
        typesSlope = ['velocity', 'slope']
        typesGrids = ['grids', 'gridnessScore']

        xlabel = 'Bump slope (neurons/s/pA)'
        ylabel = 'Gridness score'

        fig = plt.figure(figsize=scatterFigSize)
        scatterPlot = scatter.FullScatterPlot(
                ps.v, ps.grids, typesSlope, typesGrids, iter_list, velNTrials,
                gridNTrials,
                s=25,
                linewidth=0.3,
                color2D=True,
                xlabel=xlabel,
                ylabel=ylabel,
                sigmaTitle=True,
                noise_sigmas=ps.noise_sigmas,
                ignoreNaNs=ignoreNaNs,
                captionLetters=('A', 'B', 'C'),
                fig=fig)
        scatterPlot.plot()
        scatterPlot.set_titleSizes(16)

        fig.savefig(fname, dpi=300)


##############################################################################
# Scatter plot of line fit error vs line slope
class LineErrSlopeScatterPlotter(FigurePlotter):
    def __init__(self, *args, **kwargs):
        super(LineErrSlopeScatterPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ps = self.env.ps
        output_dir = self.config['output_dir']
        iter_list = self.config['iter_list']

        scatterFigSize = (8.27, 11.69)
        ignoreNaNs = True
        fname = output_dir + "/suppFigure_line_fit_error_vs_slope.pdf"

        NTrials = None
        typesSlope = ['velocity', 'slope']
        typesFitErr = ['velocity', 'fitErr']

        xlabel = 'Bump slope (neurons/s/pA)'
        ylabel = 'Fit error (neurons/s/trial)'

        fig = plt.figure(figsize=scatterFigSize)
        scatterPlot = scatter.FullScatterPlot(
                ps.v, ps.v, typesSlope, typesFitErr, iter_list, NTrials, NTrials,
                s=25,
                linewidth=0.3,
                color2D=True,
                xlabel=xlabel,
                ylabel=ylabel,
                sigmaTitle=True,
                noise_sigmas=ps.noise_sigmas,
                ignoreNaNs=ignoreNaNs,
                captionLetters=('A', 'B', 'C'),
                fig=fig)
        scatterPlot.plot()
        scatterPlot.set_titleSizes(16)

        #ax.xaxis.set_major_locator(ti.MultipleLocator(10))
        #ax.yaxis.set_major_locator(ti.MultipleLocator(0.5))
        #ax.yaxis.set_minor_locator(ti.MultipleLocator(0.25))

        fig.savefig(fname, dpi=300)

