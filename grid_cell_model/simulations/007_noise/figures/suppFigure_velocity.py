#!/usr/bin/env python
#
'''
Noise publication: supplementary figure: bump speed/velocity estimations.
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ti
from matplotlib.transforms import Bbox

from EI_plotting            import sweeps
from EI_plotting            import aggregate as aggr
from EI_plotting.base       import plotOneHist, NoiseDataSpaces
from grid_cell_model.parameters.param_space import JobTrialSpace2D
from plotting.global_defs   import globalAxesSettings
from submitting import flagparse

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
velDataRoot = 'output_local/even_spacing/velocity_vertical'
shape = (31, 31)


parser = flagparse.FlagParser()
parser.add_flag('--velSweep')
parser.add_flag('--hists')
args = parser.parse_args()

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
        var = np.abs(aggr.aggregate2D(sp, varList, funReduce=None))
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
err_vmax = 3


def createSweepFig(name):
    fig = plt.figure(name, figsize=sweepFigSize)
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    return fig, ax


exampleRC = ( (5, 15), (15, 5) )
cbar_kw = dict(
    label       = 'Fit error (neurons/s)',
    orientation = 'vertical',
    shrink      = 0.8,
    pad         = 0.05,
    ticks       = ti.MultipleLocator(1),
    extend      = 'max',
    extendfrac  = 0.1,
    rasterized  = True)

if args.velSweep or args.all:
    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        fig, ax = createSweepFig(None)
        kw = {'cbar': False}
        if (ns_idx != 0):
            kw['ylabel'] = ''
            kw['yticks'] = False
        if (ns_idx == 2):
            kw['cbar'] = True
        _, ax, cax = sweeps.plotVelTrial(ps.v[ns_idx], errVarList, iterList,
                noise_sigma,
                ax=ax,
                cbar_kw=cbar_kw,
                xlabel='', xticks=False,
                vmin=err_vmin, vmax=err_vmax,
                **kw)
        fname = outputDir + "/suppFigure_velocity_err_sweeps{}.pdf"
        fig.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)


# Stats
if args.hists or args.all:
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



