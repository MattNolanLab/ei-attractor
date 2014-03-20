#!/usr/bin/env python
'''
Figure showing where bumps are and where they are not.
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ti

import default_settings as ds
from EI_plotting          import sweeps, examples, base
from EI_plotting          import aggregate as aggr
from plotting.global_defs import globalAxesSettings, prepareLims
from submitting import flagparse


outputDir = ds.figOutputDir
ps = ds.getDefaultParamSpaces()

exampleIdx   = [(0, 0), (0, 0), (0, 0)] # (row, col)

parser = flagparse.FlagParser()
parser.add_flag('--isBump')
parser.add_flag('--fracTotalHist')
parser.add_flag('--fracTotalSweepAnn')
args = parser.parse_args()


bump_vmin = 0
bump_vmax = 1.
fracTotalText = 'P(bumps)'

#histFigSize = (5, sweepFigSize[1])
histBins = 20
if args.fracTotalHist or args.all:
    fig = plt.figure(figsize=(2.2, ds.sweepFigSize[1]*.85))
    ax = fig.add_subplot(111)
    globalAxesSettings(ax)
    data = aggr.IsBumpCollapsed(ps.bumpGamma, ds.iterList, ignoreNaNs=True)
    fractions, X, Y = data.getData()
    ax.hist(fractions, normed=True, rwidth=0.8, linewidth=0, bins=histBins)
    ax.set_xlabel(fracTotalText)
    ax.set_ylabel('p[%s]' % fracTotalText)
    ax.margins(0.02, None)
    #ax.yaxis.tick_right()
    #ax.yaxis.set_label_position('right')
    ax.xaxis.set_major_locator(ti.MultipleLocator(.5))
    ax.xaxis.set_minor_locator(ti.MultipleLocator(.1))
    ax.yaxis.set_major_locator(ti.MultipleLocator(4))
    ax.yaxis.set_minor_locator(ti.MultipleLocator(2))
    f = ti.ScalarFormatter()
    f.set_powerlimits((0, 4))
    f.set_scientific(True)
    ax.xaxis.set_major_formatter(f)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    #ax.yaxis.set_tick_params(which='both', direction='in')
    ax.set_ylim(prepareLims((0, 12), 0.01))
    fig.tight_layout()
    fname = outputDir + "/bumps_isBumpFracTotal_hist.pdf"
    fig.savefig(fname, dpi=300, transparent=True)
    plt.close()


ann150_0 = dict(
        txt='a',
        rc=(4, 4),
        xytext_offset=(1, 1.5),
        color='white')
ann0_0 = dict(
        txt='b',
        rc=(5, 15),
        xytext_offset=(1.5, 1),
        color='white')
ann0_1 = dict(
        txt='c',
        rc=(20, 15),
        xytext_offset=(1.5, 1),
        color='black')
ann300_0 = dict(
        txt='d',
        rc=(20, 25),
        xytext_offset=(-1.5, 1),
        color='white')
ann0_2 = dict(
        txt='e',
        rc=(15, 5),
        xytext_offset=(.5, 2),
        color='white')
ann150_1 = dict(
        txt='f',
        rc=(5, 15),
        xytext_offset=(1.5, 1),
        color='white')
ann150_2 = dict(
        txt='g',
        rc=(15, 5),
        xytext_offset=(1.5, 1),
        color='white')
ann300_1 = dict(
        txt='h',
        rc=(15, 5),
        xytext_offset=(1.5, -1),
        color='white')


ann0   = [  ann0_0,   ann0_1,   ann0_2]
ann150 = [ann150_0, ann150_1, ann150_2]
ann300 = [ann300_0, ann300_1]
ann = [ann0, ann150, ann300]

sweepAnn_cbar_kw = dict(
        label       = fracTotalText,
        location    = 'left',
        shrink      = 0.8,
        pad         = 0.25,
        ticks       = ti.MultipleLocator(0.5),
        rasterized  = True)
if args.fracTotalSweepAnn or args.all:
    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        fig, ax = ds.getDefaultSweepFig(scale=0.8, colorBarPos='left')
        kw = dict(cbar=False)
        if ns_idx != 0:
            kw['ylabel'] = ''
            kw['yticks'] = False
        if ns_idx == 0:
            kw['cbar'] = True
        data = aggr.IsBump(ps.bumpGamma[ns_idx], ds.iterList, ignoreNaNs=True)
        _, _, cax = sweeps.plotSweep(data,
                noise_sigma=noise_sigma,
                xlabel='', xticks=False,
                ax=ax,
                cbar_kw=sweepAnn_cbar_kw,
                vmin=bump_vmin, vmax=bump_vmax,
                annotations=ann[ns_idx],
                **kw)
        fname = outputDir + "/bumps_isBumpFracTotal_sweeps_annotated{0}.pdf"
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
                cbar_kw=sweepAnn_cbar_kw,
                vmin=bump_vmin, vmax=bump_vmax,
                **kw)
        fname = outputDir + "/bumps_isBump_sweeps{0}.pdf"
        fig.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)
        plt.close()



