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
from data_storage.sim_models import ei
from analysis import spikes as aspikes
from submitting import flagparse


outputDir = ds.figOutputDir
ps = ds.getDefaultParamSpaces()

exampleIdx   = [(0, 0), (0, 0), (0, 0)] # (row, col)

parser = flagparse.FlagParser()
parser.add_flag('--fracTotalSweep')
parser.add_flag('--isBump')
parser.add_flag('--fracTotalHist')
parser.add_flag('--bumpExamples')
parser.add_flag('--fracTotalSweepAnn')
args = parser.parse_args()


bumpNTrials = 5
bump_vmin = 0
bump_vmax = 1.
fracTotalText = '$P_{time}$(Bump formed)'

#histFigSize = (5, sweepFigSize[1])
histBins = 20
if args.fracTotalHist or args.all:
    fig = plt.figure(figsize=(2.2, ds.sweepFigSize[1]*.85))
    ax = fig.add_subplot(111)
    globalAxesSettings(ax)
    data = aggr.IsBumpCollapsed(ps.bumpGamma, ds.iterList, bumpNTrials,
            ignoreNaNs=True)
    fractions, X, Y = data.getData()
    ax.hist(fractions, normed=True, rwidth=0.8, linewidth=0, bins=histBins)
    ax.set_xlabel(fracTotalText)
    ax.set_ylabel('p(%s)' % fracTotalText)
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



##############################################################################
# Plot examples of bump snapshots for different settings of row/col and
# annotated sweeps corresponding to the examples.
trialNum = 0
nSnapshots = 8
exampleFigSize = (8, 1.9)

class ExampleSetting(object):
    def __init__(self, r, c, trialNum, ps, noise_sigma):
        self.r = r
        self.c = c
        self.trialNum = trialNum
        self.ps = ps
        self.noise_sigma = noise_sigma

exampleSettings = [
        ExampleSetting( 4,  4, trialNum, ps.bumpGamma[1], ps.noise_sigmas[1]),
        ExampleSetting( 5, 15, trialNum, ps.bumpGamma[0], ps.noise_sigmas[0]),
        ExampleSetting(20, 15, trialNum, ps.bumpGamma[0], ps.noise_sigmas[0]),
        ExampleSetting(20, 25, trialNum, ps.bumpGamma[2], ps.noise_sigmas[2]),
        ExampleSetting(15,  5, trialNum, ps.bumpGamma[0], ps.noise_sigmas[0]),
        ExampleSetting( 5, 15, trialNum, ps.bumpGamma[1], ps.noise_sigmas[1]),
        ExampleSetting(15,  5, trialNum, ps.bumpGamma[1], ps.noise_sigmas[1]),
        ExampleSetting(15,  5, trialNum, ps.bumpGamma[2], ps.noise_sigmas[2]),
]


def extractRateMaps(ps, r, c, trialNum):
    data = ps[r][c][trialNum].data
    tStart = 0.
    win_dt = 125.      # ms
    tEnd = ei.getOption(data, 'time') - win_dt
    winLen = 250.0 # ms
    Ne_x = ei.getNetParam(data, 'Ne_x')
    Ne_y = ei.getNetParam(data, 'Ne_y')
    mon = data['spikeMon_e']
    senders, times = ei.extractSpikes(mon)
    pop = aspikes.TorusPopulationSpikes(senders, times, (Ne_x, Ne_y))
    return pop.slidingFiringRate(tStart, tEnd, win_dt, winLen)


if args.bumpExamples or args.all:
    for setIdx, set in enumerate(exampleSettings):
        fig = plt.figure(figsize=exampleFigSize)
        FR, FRt = extractRateMaps(set.ps, set.r, set.c, set.trialNum)
        isBump = aggr.IsBump(set.ps, ds.iterList, bumpNTrials,
                ignoreNaNs=True)
        bumpQuality, _, _ = isBump._getRawData()
        timeTitles = setIdx == 0
        examples.plotBumpSnapshots(FR, FRt, nSnapshots,
                fig=fig, bumpQuality=bumpQuality[set.r, set.c, set.trialNum],
                timeTitles=timeTitles, maxRate=False,
                bumpQualityText="$P_{time}$\nBump\nis\nformed")
        fname = outputDir + "/bumps_isBumpSnapshotExamples_{0}pA_{1}_{2}.pdf"
        fig.savefig(fname.format(set.noise_sigma, set.r, set.c), dpi=300, transparent=True)
        plt.close()


es = exampleSettings
ann0_0 = dict(
        txt='b',
        rc=(es[1].r, es[1].c),
        xytext_offset=(1.5, 1),
        color='white')
ann0_1 = dict(
        txt='c',
        rc=(es[2].r, es[2].c),
        xytext_offset=(1.5, 1),
        color='black')

ann0_2 = dict(
        txt='e',
        rc=(es[4].r, es[4].c),
        xytext_offset=(.5, 2),
        color='white')

ann150_0 = dict(
        txt='a',
        rc=(es[0].r, es[0].c),
        xytext_offset=(1, 1.5),
        color='white')
ann150_1 = dict(
        txt='f',
        rc=(es[5].r, es[5].c),
        xytext_offset=(1.5, 1),
        color='white')
ann150_2 = dict(
        txt='g',
        rc=(es[6].r, es[6].c),
        xytext_offset=(1.5, 1),
        color='white')

ann300_0 = dict(
        txt='d',
        rc=(es[3].r, es[3].c),
        xytext_offset=(-1.5, 1),
        color='white')
ann300_1 = dict(
        txt='h',
        rc=(es[7].r, es[7].c),
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
        data = aggr.IsBump(ps.bumpGamma[ns_idx], ds.iterList, bumpNTrials,
                ignoreNaNs=True)
        _, _, cax = sweeps.plotSweep(data,
                noise_sigma=noise_sigma,
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
bumpThreshold = 0.9
if args.isBump or args.all:
    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        fig, ax = ds.getDefaultSweepFig(scale=.8, colorBarPos='left' )
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
                cbar_kw=sweepAnn_cbar_kw,
                vmin=bump_vmin, vmax=bump_vmax,
                **kw)
        fname = outputDir + "/bumps_isBump_sweeps{0}.pdf"
        fig.savefig(fname.format(int(noise_sigma)), dpi=300, transparent=True)
        plt.close()



