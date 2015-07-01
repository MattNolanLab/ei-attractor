'''Examples of consecutive snapshots of network activity.

.. currentmodule:: noisefigs.plotters.isbump_examples

Figure showing examples of consecutive snapshots of population firing rates,
together with the fraction of snaphots classified as those with bumps.

Classes
-------

.. autosummary::

    IsBumpExamplePlotter
'''
from __future__ import absolute_import, print_function

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from simtools.plotting.plotters import FigurePlotter

from ..EI_plotting import examples, base
from ..EI_plotting import aggregate as aggr
from .base import ExampleSetting

__all__ = [
    'IsBumpExamplePlotter',
]

##############################################################################
trialNum = 0
snapshot_tstep = 10
exampleFigSize = (8, 1.2)

class IsBumpExamplePlotter(FigurePlotter):
    '''Plot examples of successive bump snapshots.

    These snapshots can be plotted for different settings of row/col and
    annotated sweeps corresponding to the examples.
    '''
    def __init__(self, *args, **kwargs):
        super(IsBumpExamplePlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        ps = self.env.ps
        rateColors = self.myc['rateColors']
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
            maxRate = None
            for timeTitles in [False, True]:
                iter_list = self.config['iter_list']
                fig = self._get_final_fig(exampleFigSize)

                FR, FRt = base.extractRateMaps(set.ps, set.r, set.c, set.trialNum)
                isBump = aggr.IsBump(set.ps, iter_list, ignoreNaNs=True)
                bumpQuality, _, _ = isBump._getRawData()
                maxRate = examples.plotBumpSnapshots(FR, FRt, snapshot_tstep,
                        fig=fig, bumpQuality=bumpQuality[set.r, set.c, set.trialNum],
                        timeTitles=timeTitles, maxRate=False,
                        bumpQualityText="P\n(bumps)",
                        bumpQualityX=self.myc['bumpQualityX'],
                        maxRateColor=rateColors[setIdx])
                timeStr = '_times' if timeTitles else ''
                fname = self.get_fname(
                    "/bumps_isBumpSnapshotExamples_{ns}pA_{r}_{c}{timeStr}.pdf",
                    ns=set.noise_sigma,
                    r=set.r,
                    c=set.c,
                    timeStr=timeStr)
                fig.savefig(fname, dpi=300, transparent=True)
                plt.close()

            cbar_fig = self._get_final_fig(self.myc['cbar_fig_size'])
            ax_cbar = cbar_fig.add_axes([0.05, 0.07, 0.1, 0.8])
            cbar = mpl.colorbar.ColorbarBase(ax_cbar, cmap=mpl.cm.jet,
                                             norm=mpl.colors.Normalize(vmin=0,
                                                                       vmax=1),
                                             ticks=[0, 1],
                                             orientation='vertical')
            ax_cbar.yaxis.set_ticklabels(['0', "%.0f" % maxRate])
            ax_cbar.set_ylabel('R (Hz)')
            cbar_fname = self.get_fname(
                "/bumps_isBumpSnapshotExamples_{ns}pA_{r}_{c}_colorbar.pdf",
                ns=set.noise_sigma,
                r=set.r,
                c=set.c)
            plt.savefig(cbar_fname, dpi=300, transparent=True)
            plt.close()



