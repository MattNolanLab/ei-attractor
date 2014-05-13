#!/usr/bin/env python
#
'''
Bump formation: long refractory periods
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ti

import default_settings as ds
from EI_plotting          import sweeps
from EI_plotting          import aggregate as aggr
from plotting.global_defs import prepareLims
from grid_cell_model.parameters           import JobTrialSpace2D
from submitting import flagparse

outputDir = ds.figOutputDir + '/long_refractory'

parser = flagparse.FlagParser()
parser.add_flag('--fracTotalSweep')
parser.add_flag('--isBump')
args = parser.parse_args()


###############################################################################
bumpDataRoot     = 'output_local/even_spacing/gamma_bump_long_refr/0pA'
bumpPS = JobTrialSpace2D(ds.evenShape, bumpDataRoot)

exW = 4
exH = 2
exMargin = 0.075
exWspace=0.2
exHspace=0.15

sweepFigSize = (3.7, 2.6)
sweepLeft   = 0.08
sweepBottom = 0.2
sweepRight  = 0.8
sweepTop    = 0.85
transparent  = True


##############################################################################
# Bump formation sweeps
bump_vmin = 0
bump_vmax = 1.
fracTotalText = 'P(bumps)'

sweep_cbar_kw = dict(
        label       = fracTotalText,
        location    = 'left',
        shrink      = 0.8,
        pad         = 0.25,
        ticks       = ti.MultipleLocator(0.5),
        rasterized  = True)
if args.fracTotalSweep or args.all:
    for ns_idx, noise_sigma in enumerate([0]):
        fig, ax = ds.getDefaultSweepFig(scale=0.8, colorBarPos='left')
        kw = dict(cbar=False)
        if ns_idx != 0:
            kw['ylabel'] = ''
            kw['yticks'] = False
        if ns_idx == 0:
            kw['cbar'] = True
        data = aggr.IsBump(bumpPS, ds.iterList, ignoreNaNs=True)
        _, _, cax = sweeps.plotSweep(data,
                noise_sigma=noise_sigma,
                xlabel='', xticks=False,
                ax=ax,
                cbar_kw=sweep_cbar_kw,
                vmin=bump_vmin, vmax=bump_vmax,
                **kw)
        fname = outputDir + "/bumps_isBumpFracTotal_sweeps{0}.pdf"
        fig.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)
        plt.close()


##############################################################################
# Bump formation, thresholded
bumpThreshold = 0.95
if args.isBump or args.all:
    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        fig, ax = ds.getDefaultSweepFig(scale=.8, colorBarPos='left' )
        kw = dict(cbar=False, sigmaTitle=False)
        if ns_idx != 0:
            kw['ylabel'] = ''
            kw['yticks'] = False
        data = aggr.BumpFormationFilter(bumpThreshold, ps.bumpGamma[ns_idx],
                ds.iterList, ignoreNaNs=True)
        _, _, cax = sweeps.plotSweep(data,
                noise_sigma=noise_sigma,
                ax=ax,
                cbar_kw=sweep_cbar_kw,
                vmin=bump_vmin, vmax=bump_vmax,
                annotations=ann[ns_idx],
                **kw)
        fname = outputDir + "/bumps_mainFig_isBump_sweeps{0}.pdf"
        fig.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)
        plt.close()

