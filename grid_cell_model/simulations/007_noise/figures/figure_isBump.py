#!/usr/bin/env python
'''
Figure showing where bumps are and where they are not.
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ti

import default_settings as ds
from EI_plotting          import sweeps, examples
from EI_plotting          import aggregate as aggr
from plotting.global_defs import globalAxesSettings, prepareLims
from submitting import flagparse


outputDir = ds.figOutputDir
ps = ds.getDefaultParamSpaces()

exampleIdx   = [(0, 0), (0, 0), (0, 0)] # (row, col)

parser = flagparse.FlagParser()
parser.add_flag('--fracTotalSweep')
parser.add_flag('--isBump')
parser.add_flag('--fracTotalHist')
args = parser.parse_args()


##############################################################################
# Bump sigma sweeps
fracTotalText = 'Bump quality'
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
        fig, ax = ds.getDefaultSweepFig()
        kw = dict(cbar=False, xlabel='', xticks=False)
        if ns_idx != 0:
            kw['ylabel'] = ''
            kw['yticks'] = False
        if ns_idx == 2:
            kw['cbar'] = True
        data = aggr.IsBump(ps.bumpGamma[ns_idx], ds.iterList, bumpNTrials,
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
        fig, ax = ds.getDefaultSweepFig()
        kw = dict(cbar=False, sigmaTitle=False)
        if ns_idx != 0:
            kw['ylabel'] = ''
            kw['yticks'] = False
        data = aggr.BumpFormationFilter(bumpThreshold,
                ps.bumpGamma[ns_idx], ds.iterList, bumpNTrials,
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


#histFigSize = (5, sweepFigSize[1])
histBins = 40
if args.fracTotalHist or args.all:
    fig = plt.figure(figsize=ds.sweepFigSize)
    ax = fig.add_subplot(111)
    globalAxesSettings(ax)
    data = aggr.IsBumpCollapsed(ps.bumpGamma, ds.iterList, bumpNTrials,
            ignoreNaNs=True)
    fractions, X, Y = data.getData()
    ax.hist(fractions, normed=True, rwidth=0.8, linewidth=0, bins=histBins)
    ax.set_xlabel('Bump quality')
    ax.set_ylabel('p(Bump quality)')
    ax.margins(0.01, None)
    ax.set_ylim(prepareLims((0, 20), 0.01))
    fig.tight_layout()
    fname = outputDir + "/bumps_isBumpFracTotal_hist.pdf"
    fig.savefig(fname, dpi=300, transparent=True)
    plt.close()


