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
    'VelocityRasterPlotter',
    'VelocityRatePlotter',
    'VelocityRasterZoomPlotter',
    'GridsLineErrScatterPlotter',
    'GridsLineSlopeScatterPlotter',
    'LineErrSlopeScatterPlotter',
    'VelFitErrSweepPlotter',
    'VelFitStdSweepPlotter',
    'VelSlopeSweepPlotter',
    'VelLinesPlotter',
]

###############################################################################

rasterRC      = [(5, 15), (5, 15), (5, 15)] # (row, col)


################################################################################
# Examples of velocity fitting lines
class VelLinesPlotter(FigurePlotter):
    def __init__(self, *args, **kwargs):
        super(VelLinesPlotter, self).__init__(*args, **kwargs)

    def plotSlopes(self, ax, dataSpace, pos, noise_sigma, iterList, **kw):
        # kwargs
        trialNum   = kw.pop('trialNum', 0)
        markersize = kw.pop('markersize', 4)
        color      = kw.pop('color', 'blue')
        lineColor  = kw.pop('lineColor', 'black')
        xlabel     = kw.pop('xlabel', 'Velocity current (pA)')
        ylabel     = kw.pop('ylabel', 'Bump speed\n(neurons/s)')
        xticks     = kw.pop('xticks', True)
        yticks     = kw.pop('yticks', True)
        g_ann      = kw.pop('g_ann', True)
        sigma_ann  = kw.pop('sigma_ann', True)
        kw['markeredgecolor'] = color

        r = pos[0]
        c = pos[1]
        d = dataSpace[r][c].getAllTrialsAsDataSet().data
        a = d['analysis']
        IvelVec = dataSpace[r][c][trialNum].data['IvelVec']
        slopes = a['bumpVelAll']
        lineFit = a['lineFitLine']
        fitIvelVec = a['fitIvelVec']

        nTrials = slopes.shape[0]
        avgSlope = np.mean(slopes, axis=0)
        stdSlope = np.std(slopes, axis=0)

        if (ax is None):
            ax = plt.gca()
        plt.hold('on')
        globalAxesSettings(ax)

        ax.plot(IvelVec, slopes.T, 'o', markerfacecolor='none', markersize=markersize, **kw)
        #ax.errorbar(IvelVec, avgSlope, stdSlope, fmt='o-',
        #        markersize=markersize, color=color, alpha=0.5, **kw)
        ax.plot(fitIvelVec, lineFit, '-', linewidth=1, color=lineColor, **kw)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.xaxis.set_major_locator(ti.MultipleLocator(50))
        #ax.xaxis.set_minor_locator(ti.AutoMinorLocator(5))
        ax.yaxis.set_major_locator(ti.MultipleLocator(25))
        #ax.yaxis.set_minor_locator(ti.MultipleLocator(10))
        ax.margins(0.05)

        if (not xticks):
            ax.xaxis.set_ticklabels([])
        if (not yticks):
            ax.yaxis.set_ticklabels([])

        # Annotations
        if (sigma_ann):
            sigma_txt = '$\sigma$ = {0} pA'.format(noise_sigma)
            ax.set_title(sigma_txt, y=1.3, va='bottom')

        if (g_ann):
            Y, X = aggr.computeVelYX(dataSpace, iterList, r, c)
            gE = Y[r, c]
            gI = X[r, c]
            g_txt = '$g_E$ = {0}\n$g_I$ = {1} nS'.format(gE, gI)
        else:
            g_txt = ''

        txt = '{0}'.format(g_txt)
        ax.text(0.05, 1.3, txt, transform=ax.transAxes, va='top',
                ha='left', size='x-small')

    def plot(self, *args, **kwargs):
        ps = self.env.ps
        output_dir = self.config['output_dir']
        iter_list = self.config['iter_list']
        l, b, r, t = self.myc['bbox_rect']

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            fig = self._get_final_fig(self.myc['fig_size'])
            ax = fig.add_axes(Bbox.from_extents(l, b, r, t))
            kwargs = dict()
            if ns_idx != 1:
                kwargs['xlabel'] = ''
            if ns_idx != 0:
                kwargs['ylabel'] = ''
            self.plotSlopes(
                ax,
                ps.v[ns_idx],
                self.myc['positions'][ns_idx],
                noise_sigma=noise_sigma,
                iterList=iter_list,
                color='blue',
                **kwargs)

            fname = (self.config['output_dir'] +
                     "/velocity_slope_examples_{0}.pdf".format(int(noise_sigma)))
            plt.savefig(fname, dpi=300, transparent=True)
        

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


##############################################################################
#                           Velocity (slope) sweeps
##############################################################################
class VelSlopeSweepPlotter(SweepPlotter):
    def __init__(self, *args, **kwargs):
        super(VelSlopeSweepPlotter, self).__init__(*args, **kwargs)
        self.fig = None
        self.ax = None

    def plot(self, *args, **kwargs):
        sweepc = self._get_sweep_config()
        ps = self.env.ps
        iter_list = self.config['iter_list']
        grid_example_idx = self.config['grids']['example_idx']

        # This should be corresponding to the velLine examples as well !!
        slopeVarList = ['lineFitSlope']
        slope_vmin = -.472
        slope_vmax = 1.353

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            fname = (self.config['output_dir'] +
                     "/velocity_slope_sweeps{0}.pdf".format(int(noise_sigma)))
            self.data = aggr.GridnessScore(ps.grids[ns_idx], iter_list,
                                           ignoreNaNs=True, normalizeTicks=True,
                                           r=grid_example_idx[ns_idx][0],
                                           c=grid_example_idx[ns_idx][1])
            with self.figure_and_axes(fname, sweepc) as (self.fig, self.ax):
                kw = dict()
                if (ns_idx != 0):
                    kw['ylabel'] = ''
                    kw['yticks'] = False
                sweeps.plotVelTrial(
                    ps.v[ns_idx],
                    slopeVarList,
                    iter_list,
                    noise_sigma, sigmaTitle=False,
                    ax=self.ax,
                    cbar=self.myc['cbar'][ns_idx],
                    cbar_kw=self.myc['cbar_kw'],
                    vmin=slope_vmin, vmax=slope_vmax,
                    **kw)

                # Contours
                if self.myc['plot_contours'][ns_idx]:
                    contours = sweeps.Contours(self.data,
                            self.config['sweeps']['grid_contours'])
                    contours.plot(
                            self.ax,
                            **self.config['sweeps']['contours_kwargs'])


##############################################################################
#                           Raster and rate plots
##############################################################################

class VelocityRasterPlotter(FigurePlotter):
    def __init__(self, *args, **kwargs):
        super(VelocityRasterPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ps = self.env.ps
        output_dir = self.config['output_dir']

        logger.info("Plotting rasters")
        for idx, noise_sigma in enumerate(ps.noise_sigmas):
            logger.info("   Rasters: %d pA", noise_sigma)
            fig = self._get_final_fig(self.myc['fig_size'])
            l, b, r, t = self.myc['bbox']
            ax = fig.add_axes(Bbox.from_extents(l, b, r, t))
            kw = dict(scaleBar=None)
            if (idx != 0):
                kw['ylabel'] = ''
                kw['yticks'] = False
            if idx == 2:
                kw['scaleBar'] = 125
            rasters.EIRaster(ps.v[idx], 
                    noise_sigma=noise_sigma,
                    spaceType='velocity',
                    r=rasterRC[idx][0], c=rasterRC[idx][1],
                    ylabelPos=self.config['vel_rasters']['ylabelPos'],
                    tLimits=self.config['vel_rasters']['tLimits'],
                    trialNum=self.config['vel_rasters']['trialNum'],
                    ann_EI=True,
                    scaleX=0.85,
                    scaleY=-0.15,
                    **kw)
            fname = output_dir + "/velocity_raster{0}.png"
            fig.savefig(fname.format(int(noise_sigma)), dpi=300,
                    transparent=self.myc['transparent'])
            plt.close()
        
        

##############################################################################

class VelocityRatePlotter(FigurePlotter):
    def __init__(self, *args, **kwargs):
        super(VelocityRatePlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ps = self.env.ps
        output_dir = self.config['output_dir']
        transparent = self.myc['transparent']

        for idx, noise_sigma in enumerate(ps.noise_sigmas):
            # E cells
            fig = self._get_final_fig(self.myc['fig_size'])
            l, b, r, t = self.myc['bbox']
            ax = fig.add_axes(Bbox.from_extents(l, b, r, t))
            kw = {}
            if (idx != 0):
                kw['ylabel'] = ''

            rasters.plotAvgFiringRate(ps.v[idx],
                    spaceType='velocity',
                    noise_sigma=noise_sigma,
                    popType='E',
                    r=rasterRC[idx][0], c=rasterRC[idx][1],
                    color='red',
                    ylabelPos=self.config['vel_rasters']['ylabelPos'],
                    tLimits=self.config['vel_rasters']['tLimits'],
                    trialNum=self.config['vel_rasters']['trialNum'],
                    sigmaTitle=False,
                    ax=ax, **kw)
            fname = output_dir + "/velocity_rate_e{0}.pdf".format(noise_sigma)
            fig.savefig(fname, dpi=300, transparent=transparent)
            plt.close()

            # I cells
            fig = self._get_final_fig(self.myc['fig_size'])
            ax = fig.add_axes(Bbox.from_extents(l, b, r, t))
            kw = {}
            if (idx != 0):
                kw['ylabel'] = ''

            rasters.plotAvgFiringRate(ps.v[idx],
                    spaceType='velocity',
                    noise_sigma=noise_sigma,
                    popType='I', 
                    r=rasterRC[idx][0], c=rasterRC[idx][1],
                    color='blue',
                    ylabelPos=self.config['vel_rasters']['ylabelPos'],
                    tLimits=self.config['vel_rasters']['tLimits'],
                    trialNum=self.config['vel_rasters']['trialNum'],
                    sigmaTitle=False,
                    ax=ax, **kw)
            fname = output_dir + "/velocity_rate_i{0}.pdf".format(noise_sigma)
            fig.savefig(fname, dpi=300, transparent=transparent)
            plt.close()


##############################################################################
#                           Raster and rate zoom-ins
##############################################################################
class VelocityRasterZoomPlotter(FigurePlotter):
    def __init__(self, *args, **kwargs):
        super(VelocityRasterZoomPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ps = self.env.ps
        output_dir = self.config['output_dir']
        l, b, r, t = self.myc['bbox']
        transparent = self.myc['transparent']

        tLimits  = [2.75e3, 2.875e3] # ms
        trialNum = 0

        for idx, noise_sigma in enumerate(ps.noise_sigmas):
            fig = self._get_final_fig(self.myc['fig_size'])
            ax = fig.add_axes(Bbox.from_extents(l, b, r, t))
            kw = dict(scaleBar=None)
            if idx == 2:
                kw['scaleBar'] = 25
            rasters.EIRaster(ps.v[idx], 
                    noise_sigma=noise_sigma,
                    spaceType='velocity',
                    r=rasterRC[idx][0], c=rasterRC[idx][1],
                    ylabelPos=self.myc['ylabelPos'],
                    tLimits=tLimits,
                    trialNum=trialNum,
                    sigmaTitle=False,
                    ann_EI=True,
                    scaleX=0.75,
                    scaleY=-0.15,
                    ylabel='', yticks=False,
                    **kw)
            fname = output_dir + "/velocity_raster_zooms{0}.png"
            fig.savefig(fname.format(int(noise_sigma)), dpi=300,
                    transparent=transparent)
            plt.close()
       

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



##############################################################################
# Bump velocity line fit error
#exampleRC = ( (5, 15), (15, 5) )
class VelFitErrSweepPlotter(SweepPlotter):
    def __init__(self, *args, **kwargs):
        super(VelFitErrSweepPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        myc= self._get_class_config()
        sweepc = self._get_sweep_config()
        ps = self.env.ps
        iter_list = self.config['iter_list']
        grid_example_idx = self.config['grids']['example_idx']

        errVarList = ['lineFitErr']
        err_vmin = 0
        err_vmax = 11.2

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            fname = (self.config['output_dir'] +
                     "/suppFigure_velocity_err_sweeps{}.pdf".format(int(noise_sigma)))
            self.data = aggr.GridnessScore(ps.grids[ns_idx], iter_list,
                                           ignoreNaNs=True, normalizeTicks=True,
                                           r=grid_example_idx[ns_idx][0],
                                           c=grid_example_idx[ns_idx][1])
            with self.figure_and_axes(fname, sweepc) as (fig, ax):
                kw = {'cbar': False}
                if (ns_idx != 0):
                    kw['ylabel'] = ''
                    kw['yticks'] = False
                if (ns_idx == 2):
                    kw['cbar'] = True
                _, ax, cax = sweeps.plotVelTrial(ps.v[ns_idx], errVarList,
                        iter_list,
                        noise_sigma,
                        ax=ax,
                        cbar_kw=myc['cbar_kw'],
                        vmin=err_vmin, vmax=err_vmax,
                        sigmaTitle=False,
                        **kw)

                # Contours
                if self.myc['plot_contours'][ns_idx]:
                    contours = sweeps.Contours(self.data,
                            self.config['sweeps']['grid_contours'])
                    contours.plot(
                            ax,
                            **self.config['sweeps']['contours_kwargs'])



## Stats
#if args.hists or args.all:
#    fig = plt.figure(figsize=histFigsize)
#    ax = fig.add_axes(Bbox.from_extents(histLeft, histBottom, histRight,
#        histTop))
#    plotErrHistogram(ps.v, ['lineFitErr'], xlabel='Fit error (neurons/s)',
#            ylabel='p(error)')
#    fname = outputDir + "/suppFigure_velocity_err_histograms.pdf"
#    plt.savefig(fname, dpi=300, transparent=True)
#
#    fig = plt.figure(figsize=histFigsize)
#    ax = fig.add_axes(Bbox.from_extents(histLeft, histBottom, histRight,
#        histTop))
#    plotSlopeHistogram(ps.v, ['lineFitSlope'], xlabel='Slope (neurons/s/pA)',
#            ylabel='p(slope)', plotLegend=True)
#    fname = outputDir + "/suppFigure_velocity_slope_histograms.pdf"
#    plt.savefig(fname, dpi=300, transparent=True)


################################################################################
# Velocity line fit std. deviation
class VelFitStdSweepPlotter(SweepPlotter):
    def __init__(self, *args, **kwargs):
        super(VelFitStdSweepPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        myc= self._get_class_config()
        sweepc = self._get_sweep_config()
        ps = self.env.ps

        iter_list = self.config['iter_list']
        output_dir = self.config['output_dir']

        std_vmin = 0
        std_vmax = 10.421

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            fname = (output_dir +
                    "/bumps_vel_std_sweeps{0}.pdf".format(int(noise_sigma)))
            with self.figure_and_axes(fname, sweepc) as (fig, ax):
                kw = {'cbar': False}
                if (ns_idx != 0):
                    kw['ylabel'] = ''
                    kw['yticks'] = False
                if (ns_idx == 2):
                    kw['cbar'] = True
                _, ax, cax = sweeps.plotVelStdSweep(ps.v[ns_idx], iter_list,
                        ps.noise_sigmas[ns_idx],
                        ax=ax,
                        sigmaTitle=False,
                        cbar_kw=myc['cbar_kw'],
                        vmin=std_vmin, vmax=std_vmax,
                        **kw)

