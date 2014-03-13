#!/usr/bin/env python
#
'''
Figure showing where bumps are and where they are not.
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ti
from matplotlib.gridspec   import GridSpec
from matplotlib.colorbar   import make_axes
from matplotlib.transforms import Bbox

from EI_plotting          import sweeps, details, examples, rasters
from EI_plotting          import aggregate as aggr
from plotting.global_defs import globalAxesSettings, prepareLims
from EI_plotting.base     import NoiseDataSpaces, getNoiseDataSpaces
from parameters           import JobTrialSpace2D
from analysis             import clustering
from submitting import flagparse

from matplotlib import rc
rc('pdf', fonttype=42)
rc('mathtext', default='regular')

plt.rcParams['font.size'] = 11

outputDir = "panels/"

NTrials=10
iterList  = ['g_AMPA_total', 'g_GABA_total']

noise_sigmas = [0, 150, 300]
exampleIdx   = [(0, 0), (0, 0), (0, 0)] # (row, col)
gridsDataRoot    = 'output_local/even_spacing/grids'
bumpDataRoot     = 'output_local/even_spacing/gamma_bump'
velDataRoot      = 'output_local/even_spacing/velocity_vertical'
constPosDataRoot = 'output_local/even_spacing/const_position'
shape = (31, 31)

parser = flagparse.FlagParser()
parser.add_flag('--fracTotalSweep')
parser.add_flag('--isBump')
parser.add_flag('--fracTotalHist')
args = parser.parse_args()

###############################################################################
roots = NoiseDataSpaces.Roots(bumpDataRoot, velDataRoot, gridsDataRoot,
        constPos=constPosDataRoot)
ps    = NoiseDataSpaces(roots, shape, noise_sigmas)

sweepFigSize = (3.7, 2.6)
sweepLeft   = 0.08
sweepBottom = 0.2
sweepRight  = 0.8
sweepTop    = 0.85
transparent  = True


##############################################################################
# Bump sigma sweeps
fracTotalText = 'Fraction/simulation'
bumpTStart = 500.0
bumpNTrials = 5
bump_vmin = 0
bump_vmax = 1.
bump_cbar_kw = dict(
        label       = fracTotalText,
        location    = 'right',
        shrink      = 0.8,
        pad         = -0.05,
        ticks       = ti.MultipleLocator(0.5),
        rasterized  = True)

if args.fracTotalSweep or args.all:
    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        fig = plt.figure(figsize=sweepFigSize)
        ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
            sweepTop))
        kw = dict(cbar=False, xlabel='', xticks=False)
        if ns_idx != 0:
            kw['ylabel'] = ''
            kw['yticks'] = False
        if ns_idx == 2:
            kw['cbar'] = True
        data = aggr.IsBump(ps.bumpGamma[ns_idx], iterList, bumpNTrials,
                ignoreNaNs=True)
        _, _, cax = sweeps.plotSweep(data,
                noise_sigma=noise_sigma,
                ax=ax,
                cbar_kw=bump_cbar_kw,
                vmin=bump_vmin, vmax=bump_vmax,
                **kw)
        fname = outputDir + "/bumps_isBumpFracTotal_sweeps{0}.pdf"
        fig.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)
        plt.close()



##############################################################################
# Bump formation, thresholded
bumpThreshold = 0.9
if args.isBump or args.all:
    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        fig = plt.figure(figsize=sweepFigSize)
        ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
            sweepTop))
        kw = dict(cbar=False, sigmaTitle=False)
        if ns_idx != 0:
            kw['ylabel'] = ''
            kw['yticks'] = False
        data = aggr.BumpFormationFilter(bumpThreshold,
                ps.bumpGamma[ns_idx], iterList, bumpNTrials,
                ignoreNaNs=True)
        _, _, cax = sweeps.plotSweep(data,
                noise_sigma=noise_sigma,
                ax=ax,
                cbar_kw=bump_cbar_kw,
                vmin=bump_vmin, vmax=bump_vmax,
                **kw)
        fname = outputDir + "/bumps_isBump_sweeps{0}.pdf"
        fig.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)
        plt.close()


histFigSize = (5, sweepFigSize[1])
histBins = 40
if args.fracTotalHist or args.all:
    fig = plt.figure(figsize=sweepFigSize)
    ax = fig.add_subplot(111)
    globalAxesSettings(ax)
    data = aggr.IsBumpCollapsed(ps.bumpGamma, iterList, bumpNTrials,
            ignoreNaNs=True)
    fractions, X, Y = data.getData()
    ax.hist(fractions, normed=True, rwidth=0.8, linewidth=0, bins=histBins)
    ax.set_xlabel('Fraction of bump formed/simulation')
    ax.set_ylabel('$p(\cdot)$')
    ax.margins(0.01, 0.01)
    fig.tight_layout()
    fname = outputDir + "/bumps_isBumpFracTotal_hist.pdf"
    fig.savefig(fname, dpi=300, transparent=True)
    plt.close()


