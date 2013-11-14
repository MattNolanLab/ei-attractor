#!/usr/bin/env python
#
#   suppFigure_velocity.py
#
#   Noise publication: supplementary figure: bump speed/velocity estimations.
#
#       Copyright (C) 2013  Lukas Solanka <l.solanka@sms.ed.ac.uk>
#       
#       This program is free software: you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation, either version 3 of the License, or
#       (at your option) any later version.
#       
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#       
#       You should have received a copy of the GNU General Public License
#       along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ti
from matplotlib.transforms import Bbox

import EI_plotting as EI
from parameters.param_space import JobTrialSpace2D
from plotting.global_defs import globalAxesSettings
from figures_shared       import plotOneHist, NoiseDataSpaces

import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)

from matplotlib import rc
rc('pdf', fonttype=42)
rc('mathtext', default='regular')

plt.rcParams['font.size'] = 10

outputDir = "."

NTrials=10
iterList  = ['g_AMPA_total', 'g_GABA_total']

noise_sigmas = [0, 150, 300]
exampleIdx   = [(0, 0), (0, 0), (0, 0)] # (row, col)
gridsDataRoot= None
bumpDataRoot= None
velDataRoot = 'output_local/even_spacing/velocity'
shape = (31, 31)

velSweep       = 1
hists          = 1

##############################################################################



def plotVelHistogram(spList, varList, xlabel="", ylabel="", **kw):
    noise_sigma = [0, 150, 300]
    colors = ['red', 'green', 'blue']
    range = kw.get('range')
    plotLegend = kw.pop('plotLegend', False)

    ax = plt.gca()
    plt.hold('on')
    globalAxesSettings(ax)

    for idx, sp in enumerate(spList):
        var = np.abs(EI.aggregate2D(sp, varList, funReduce=None))
        filtIdx = np.logical_not(np.isnan(var))
        if (range is not None):
            var[var < range[0]] = range[0]
            var[var > range[1]] = range[1]
        plotOneHist(var[filtIdx], normed=True, **kw)

    if (plotLegend):
        leg = []
        for s in noise_sigma:
            leg.append("{0}".format(int(s)))
        l = ax.legend(leg, loc=(0.75, 0.5), title='$\sigma$ (pA)',
                frameon=False, fontsize='x-small', ncol=1)
        plt.setp(l.get_title(), fontsize='x-small')

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    f = ti.ScalarFormatter(useMathText=True)
    f.set_scientific(True)
    f.set_powerlimits([0, 3])
    ax.yaxis.set_major_formatter(f)
    ax.yaxis.set_minor_locator(ti.AutoMinorLocator(2))
    return ax

def plotErrHistogram(spList, varList, **kw):
    ax = plotVelHistogram(spList, varList, range=[0, 10], **kw)

    ax.xaxis.set_major_locator(ti.MultipleLocator(2))
    ax.xaxis.set_minor_locator(ti.AutoMinorLocator(2))
    ax.yaxis.set_major_locator(ti.MaxNLocator(4))
    ax.set_ylim([-0.0025, 2])
    #ax.margins(0.01)
    
def plotSlopeHistogram(spList, varList, **kw):
    ax = plotVelHistogram(spList, varList, range=[0, 1.5], **kw)

    ax.xaxis.set_major_locator(ti.MultipleLocator(0.4))
    ax.xaxis.set_minor_locator(ti.AutoMinorLocator(2))
    ax.yaxis.set_major_locator(ti.MaxNLocator(4))
    ax.set_ylim([-0.0025, 9])
    #ax.margins(0.01)
    

###############################################################################
roots = NoiseDataSpaces.Roots(bumpDataRoot, velDataRoot, gridsDataRoot)
ps    = NoiseDataSpaces(roots, shape, noise_sigmas)

sweepFigSize = (2.6, 1.9)
sweepLeft   = 0.15
sweepBottom = 0.2
sweepRight  = 0.87
sweepTop    = 0.85

histFigsize =(2.6, 1.7)
histLeft    = 0.22
histBottom  = 0.3
histRight   = 0.95
histTop     = 0.86


###############################################################################

errVarList = ['lineFitErr']
err_vmin = 0
err_vmax = 10


def createSweepFig(name):
    fig = plt.figure(name, figsize=sweepFigSize)
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    return fig, ax


exampleRC = ( (5, 15), (15, 5) )
ann0 = dict(
        txt='D,E',
        rc=exampleRC[0],
        xytext_offset=(1.5, 1),
        ha='center',
        color='black')
ann1 = dict(
        txt='D,E',
        rc=exampleRC[1],
        xytext_offset=(0.5, 1.5),
        ha='center',
        color='red')
ann = [ann0, ann1]
cbar_kw = dict(
    label = 'Fit error (neurons/s)',
    orientation = 'vertical',
    shrink = 0.8,
    pad = 0.05,
    ticks = ti.MultipleLocator(4),
    rasterized = True)

if (velSweep):
    # noise_sigma = 0 pA
    fig, ax = createSweepFig("velErrSweeps0")
    _, ax, cax = EI.plotVelTrial(ps.v[0], errVarList, iterList,
            noise_sigmas[0],
            ax=ax,
            cbar=False, cbar_kw=cbar_kw,
            xlabel='', xticks=False,
            vmin=err_vmin, vmax=err_vmax,
            annotations=ann)
    fname = outputDir + "/suppFigure_velocity_err_sweeps0.pdf"
    fig.savefig(fname, dpi=300, transparent=True)


    # noise_sigma = 150 pA
    fig, ax = createSweepFig("velErrSweeps150")
    _, ax, cax = EI.plotVelTrial(ps.v[1], errVarList, iterList,
            noise_sigma=noise_sigmas[1],
            ax=ax,
            ylabel='', yticks=False,
            cbar=False, cbar_kw=cbar_kw,
            xlabel='', xticks=False,
            vmin=err_vmin, vmax=err_vmax,
            annotations=ann)
    fname = outputDir + "/suppFigure_velocity_err_sweeps150.pdf"
    fig.savefig(fname, dpi=300, transparent=True)


    # noise_sigma = 300 pA
    fig, ax = createSweepFig("velErrSweeps300")
    _, ax, cax = EI.plotVelTrial(ps.v[2], errVarList, iterList,
            noise_sigma=noise_sigmas[2],
            ax=ax,
            ylabel='', yticks=False,
            cbar=True,
            xlabel='', xticks=False,
            vmin=err_vmin, vmax=err_vmax,
            cbar_kw = cbar_kw,
            annotations=ann)
    fname = outputDir + "/suppFigure_velocity_err_sweeps300.pdf"
    fig.savefig(fname, dpi=300, transparent=True)



# Stats
if (hists):
    fig = plt.figure(figsize=histFigsize)
    ax = fig.add_axes(Bbox.from_extents(histLeft, histBottom, histRight,
        histTop))
    plotErrHistogram(ps.v, ['lineFitErr'], xlabel='Fit error (neurons/s)',
            ylabel='p(error)')
    fname = outputDir + "/suppFigure_velocity_err_histograms.pdf"
    plt.savefig(fname, dpi=300, transparent=True)

    fig = plt.figure(figsize=histFigsize)
    ax = fig.add_axes(Bbox.from_extents(histLeft, histBottom, histRight,
        histTop))
    plotSlopeHistogram(ps.v, ['lineFitSlope'], xlabel='Slope (neurons/s/pA)',
            ylabel='p(slope)', plotLegend=True)
    fname = outputDir + "/suppFigure_velocity_slope_histograms.pdf"
    plt.savefig(fname, dpi=300, transparent=True)



