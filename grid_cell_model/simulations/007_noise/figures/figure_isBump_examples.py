#!/usr/bin/env python
'''
Figure showing examples of consecutive snapshots of population firing rates,
together with the fraction of snaphots classified as those with bumps.
'''
import numpy as np
import matplotlib.pyplot as plt

import default_settings as ds
from EI_plotting          import examples, base
from EI_plotting          import aggregate as aggr
from submitting import flagparse


outputDir = ds.figOutputDir
ps = ds.getDefaultParamSpaces()


parser = flagparse.FlagParser()
parser.add_flag('--bumpExamples')
args = parser.parse_args()

##############################################################################
# Plot examples of bump snapshots for different settings of row/col and
# annotated sweeps corresponding to the examples.
trialNum = 0
nSnapshots = 8
exampleFigSize = (8, 1.2)

class ExampleSetting(object):
    def __init__(self, r, c, trialNum, ps, noise_sigma):
        self.r = r
        self.c = c
        self.trialNum = trialNum
        self.ps = ps
        self.noise_sigma = noise_sigma

exampleSettings = [
        ExampleSetting( 5, 15, trialNum, ps.bumpGamma[0], ps.noise_sigmas[0]),
        ExampleSetting(20, 15, trialNum, ps.bumpGamma[0], ps.noise_sigmas[0]),
        ExampleSetting(15,  5, trialNum, ps.bumpGamma[0], ps.noise_sigmas[0]),
        ExampleSetting(14,  1, trialNum, ps.bumpGamma[0], ps.noise_sigmas[0]),
        ExampleSetting( 4,  4, trialNum, ps.bumpGamma[1], ps.noise_sigmas[1]),
        ExampleSetting( 5, 15, trialNum, ps.bumpGamma[1], ps.noise_sigmas[1]),
        ExampleSetting(15,  5, trialNum, ps.bumpGamma[1], ps.noise_sigmas[1]),
        ExampleSetting(15,  5, trialNum, ps.bumpGamma[2], ps.noise_sigmas[2]),
        ExampleSetting(20, 25, trialNum, ps.bumpGamma[2], ps.noise_sigmas[2]),
        ExampleSetting( 5, 15,        3, ps.bumpGamma[2], ps.noise_sigmas[2]),
]


if args.bumpExamples or args.all:
    for setIdx, set in enumerate(exampleSettings):
        for timeTitles in [False, True]:
            fig = plt.figure(figsize=exampleFigSize)
            FR, FRt = base.extractRateMaps(set.ps, set.r, set.c, set.trialNum)
            isBump = aggr.IsBump(set.ps, ds.iterList, ignoreNaNs=True)
            bumpQuality, _, _ = isBump._getRawData()
            examples.plotBumpSnapshots(FR, FRt, nSnapshots,
                    fig=fig, bumpQuality=bumpQuality[set.r, set.c, set.trialNum],
                    timeTitles=timeTitles, maxRate=False,
                    bumpQualityText="P\n(bumps)")
            fname = outputDir + "/bumps_isBumpSnapshotExamples_{0}pA_{1}_{2}{3}.pdf"
            timeStr = '_times' if timeTitles else ''
            fig.savefig(fname.format(set.noise_sigma, set.r, set.c, timeStr),
                    dpi=300, transparent=True)
            plt.close()


