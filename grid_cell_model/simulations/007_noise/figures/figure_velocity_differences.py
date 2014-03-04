#!/usr/bin/env python
#
'''
Everything related to velocity difference plots: mainly sweeps.
'''
import copy

from submitting import flagparse
parser = flagparse.FlagParser()
parser.add_flag('--slope_sweeps')
parser.add_flag('--fitErr_sweeps')
parser.add_flag('--scatter_diff_slope_grids')
parser.add_flag('--scatter_diff_fitErr_grids')
parser.add_flag('--slope_diff_horiz_vert')
args = parser.parse_args()


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ti
from matplotlib.transforms import Bbox

import default_settings as ds
from EI_plotting        import sweeps, scatter
from EI_plotting.base   import NoiseDataSpaces, getNoiseDataSpaces
from plotting.low_level import zeroLines


###############################################################################

outputDir = "panels"
NTrials = 5

noise_sigmas  = [0, 150, 300]
bumpDataRoot  = 'output_local/even_spacing/gamma_bump'
velDataRoot   = 'output_local/even_spacing/velocity_vertical'
gridsDataRoot = 'output_local/even_spacing/grids'
shape = (31, 31)

###############################################################################
roots = NoiseDataSpaces.Roots(bumpDataRoot, velDataRoot, gridsDataRoot)
ps    = NoiseDataSpaces(roots, shape, noise_sigmas)



slopeTypes = ['velocity', 'slope']
fitErrTypes = ['velocity', 'fitErr']
velSweep_cmap = 'RdBu_r'

##############################################################################
#                           Slope difference plots

# Slope difference sweeps
slopeDiffLabel = '$\Delta$ Slope\n(neurons/s/pA)'
slope_vmin = -1
slope_vmax = 1
cbarBase = dict(
        location='right',
        shrink = 0.8,
        pad = -0.05)


slope_cbar_kw = copy.deepcopy(cbarBase)
slope_cbar_kw.update(
        label=slopeDiffLabel,
        ticks=ti.MultipleLocator(0.5),
        extend='min', extendfrac=0.1)

def createSweepFig(name=None):
    sweepFigSize = (3.7, 2.6)
    sweepLeft   = 0.05
    sweepBottom = 0.2
    sweepRight  = 0.8
    sweepTop    = 0.85
    fig = plt.figure(name, figsize=sweepFigSize)
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    return fig, ax


if args.slope_sweeps or args.all:
    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas[0:-1]):
        fig, ax = createSweepFig(None)
        which = ns_idx
        _, ax, cax = sweeps.plotDiffTrial(ps.v, ds.iterList, which, None,
                slopeTypes,
                ax=ax,
                cbar=True, cbar_kw=slope_cbar_kw,
                vmin=slope_vmin, vmax=slope_vmax,
                cmap=velSweep_cmap,
                ignoreNaNs=True,
                symmetricLimits=False)
        fname = outputDir + "/velocity_slope_diff_sweeps{}.pdf"
        fig.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)
        plt.close()


# Correlate (difference between slope) and gridness score.
corrDiffFigsize = (4, 5)
corrDiffXLabel = slopeDiffLabel
corrDiffYLabel = '$\Delta$ Gridness score'

gridTypes = ['grids', 'gridnessScore']
gridNTrials = 3
def setCorrAxes(ax):
    #ax.set_xlim(prepareLims([-0.2, 0.4]))
    #ax.set_ylim(prepareLims([-1.5, 1.5]))
    #ax.xaxis.set_major_locator(ti.MultipleLocator(0.2))
    #ax.yaxis.set_major_locator(ti.MultipleLocator(0.5))
    #ax.xaxis.set_minor_locator(ti.MultipleLocator(0.1))
    #ax.yaxis.set_minor_locator(ti.MultipleLocator(0.25))
    zeroLines(ax)

# Highlight segmented points
if args.scatter_diff_slope_grids or args.all:
    fig = plt.Figure(corrDiffFigsize)
    ax = fig.add_subplot(111)

    which = 0
    scatterPlot = scatter.DiffScatterPlot(
            ps.v, ps.grids, slopeTypes, gridTypes, ds.iterList,
            None, gridNTrials, which,
            s=15,
            linewidth=0.3,
            edgecolor='white',
            color2D=True,
            xlabel = corrDiffXLabel,
            ylabel = corrDiffYLabel,
            sigmaTitle=False,
            ignoreNaNs=True,
            ax=ax)

    scatterPlot.plot()
    setCorrAxes(ax)

    fig.tight_layout()
    fname = outputDir + "/velocity_scatter_diff_slope_grids.pdf"
    fig.savefig(fname, dpi=300, transparent=True)
    plt.close()



##############################################################################
#                      Fit error difference plots
fitErrDiffLabel = '$\Delta$ Fit error\n(neurons/s)'

# Fit error sweeps
fitErr_cbar_kw = copy.deepcopy(cbarBase)
fitErr_cbar_kw.update(
        ticks = ti.MultipleLocator(5),
        label = fitErrDiffLabel)
fitErr_vmin = -5.404
fitErr_vmax = 5.404

if args.fitErr_sweeps or args.all:
    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas[0:-1]):
        fig, ax = createSweepFig(None)
        which = ns_idx
        _, ax, cax = sweeps.plotDiffTrial(ps.v, ds.iterList, which, None,
                fitErrTypes,
                ax=ax,
                cbar=True, cbar_kw=fitErr_cbar_kw,
                vmin=fitErr_vmin, vmax=fitErr_vmax,
                cmap=velSweep_cmap,
                ignoreNaNs=True,
                symmetricLimits=False)
        fname = outputDir + "/velocity_fitErr_diff_sweeps{}.pdf"
        fig.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)
        plt.close()



# Correlate (difference between line fit error) and gridness score.
corrDiffFigsize = (4, 5)
corrDiffXLabel = fitErrDiffLabel
corrDiffYLabel = '$\Delta$ Gridness score'

fitErrTypes = ['velocity', 'fitErr']


# Highlight segmented points
if args.scatter_diff_fitErr_grids or args.all:
    fig = plt.Figure(corrDiffFigsize)
    ax = fig.add_subplot(111)

    which = 0
    scatterPlot = scatter.DiffScatterPlot(
            ps.v, ps.grids, fitErrTypes, gridTypes, ds.iterList,
            None, gridNTrials, which,
            s=15,
            linewidth=0.3,
            edgecolor='white',
            color2D=True,
            xlabel = corrDiffXLabel,
            ylabel = corrDiffYLabel,
            sigmaTitle=False,
            ignoreNaNs=True,
            ax=ax)

    scatterPlot.plot()
    setCorrAxes(ax)

    fig.tight_layout()
    fname = outputDir + "/velocity_scatter_diff_fitErr_grids.pdf"
    fig.savefig(fname, dpi=300, transparent=True)
    plt.close()



##############################################################################
# Difference between vertical and horizontal slope of bump speed estimation.
horizVelDataRoot   = 'output_local/even_spacing/velocity'
horizPS = getNoiseDataSpaces(horizVelDataRoot, ps.noise_sigmas, shape)
horizVert_vmin = -1.05
horizVert_vmax =  1.05

if args.slope_diff_horiz_vert or args.all:
    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        kw = dict(cbar=False)
        if ns_idx == 2:
            kw['cbar'] = True
        if ns_idx != 0:
            kw['ylabel'] = ''
            kw['yticks'] = False
        fig, ax = createSweepFig(None)
        which = 0
        spaces = [horizPS[ns_idx], ps.v[ns_idx]]
        _, ax, cax = sweeps.plotDiffTrial(spaces, ds.iterList, which, None,
                slopeTypes,
                ax=ax,
                cbar_kw=slope_cbar_kw,
                vmin=horizVert_vmin, vmax=horizVert_vmax,
                cmap=velSweep_cmap,
                ignoreNaNs=True,
                symmetricLimits=False,
                **kw)
        fname = outputDir + "/velocity_slope_diff_sweeps_horiz_vert{}.pdf"
        fig.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)
        plt.close()
