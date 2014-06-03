'''
Figure showing examples of consecutive snapshots of population firing rates,
together with the fraction of snaphots classified as those with bumps.
'''
from __future__ import absolute_import, print_function

import numpy as np
import matplotlib.pyplot as plt

from ..EI_plotting import examples, base
from ..EI_plotting import aggregate as aggr
from .base import FigurePlotter

__all__ = [
    'IsBumpExamplePlotter',
]

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

class IsBumpExamplePlotter(FigurePlotter):
    def __init__(self, *args, **kwargs):
        super(IsBumpExamplePlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ps = self.env.ps
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

        for setIdx, set in enumerate(exampleSettings):
            for timeTitles in [False, True]:
                iter_list = self.config['iter_list']
                fig = self._get_final_fig(exampleFigSize)

                FR, FRt = base.extractRateMaps(set.ps, set.r, set.c, set.trialNum)
                isBump = aggr.IsBump(set.ps, iter_list, ignoreNaNs=True)
                bumpQuality, _, _ = isBump._getRawData()
                examples.plotBumpSnapshots(FR, FRt, nSnapshots,
                        fig=fig, bumpQuality=bumpQuality[set.r, set.c, set.trialNum],
                        timeTitles=timeTitles, maxRate=False,
                        bumpQualityText="P\n(bumps)")
                fname = self.config['output_dir'] + "/bumps_isBumpSnapshotExamples_{0}pA_{1}_{2}{3}.pdf"
                timeStr = '_times' if timeTitles else ''
                fig.savefig(fname.format(set.noise_sigma, set.r, set.c, timeStr),
                        dpi=300, transparent=True)
                plt.close()


