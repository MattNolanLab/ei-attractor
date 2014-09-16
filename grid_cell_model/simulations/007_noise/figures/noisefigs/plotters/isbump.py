#!/usr/bin/env python
'''
Figure showing where bumps are and where they are not.
'''
from __future__ import absolute_import, print_function

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ti

from grid_cell_model.plotting.global_defs import globalAxesSettings, prepareLims

from ..EI_plotting import sweeps, examples, base
from ..EI_plotting import aggregate as aggr
from .base import FigurePlotter, SweepPlotter

__all__ = [
    'FracTotalHistPlotter',
    'FracTotalSweepAnnPlotter',
    'IsBumpPlotter',
]

bump_vmin = 0
bump_vmax = 1.

class FracTotalHistPlotter(FigurePlotter):
    histBins = 20

    def __init__(self, *args, **kwargs):
        super(FracTotalHistPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        iter_list = self.config['iter_list']
        ps = self.env.ps
        frac_total_text = self.config['p_bumps']['frac_total_text']

        fig = self._get_final_fig(
            (2.2, self.config['sweeps']['fig_size'][1]*.85)
        )
        ax = fig.add_subplot(111)
        globalAxesSettings(ax)
        data = aggr.IsBumpCollapsed(ps.bumpGamma, iter_list, ignoreNaNs=True)
        fractions, X, Y = data.getData()
        ax.hist(fractions, normed=True, rwidth=0.8, linewidth=0, bins=self.histBins)
        ax.set_xlabel(frac_total_text)
        ax.set_ylabel('p[%s]' % frac_total_text)
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
        fname = self.config['output_dir'] + "/bumps_isBumpFracTotal_hist.pdf"
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

class FracTotalSweepAnnPlotter(SweepPlotter):
    def __init__(self, *args, **kwargs):
        super(FracTotalSweepAnnPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        myc= self._get_class_config()
        sweepc = self._get_sweep_config()
        ps = self.env.ps

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            fname = (self.config['output_dir'] +
                     "/bumps_isBumpFracTotal_sweeps_annotated{0}.pdf".format(int(noise_sigma)))
            with self.figure_and_axes(fname, sweepc) as (fig, ax):
                #fig, ax = ds.getDefaultSweepFig(scale=0.8, colorBarPos='left')
                kw = dict()
                if ns_idx != 0:
                    kw['ylabel'] = ''
                    kw['yticks'] = False
                data = aggr.IsBump(ps.bumpGamma[ns_idx],
                                   self.config['iter_list'],
                                   ignoreNaNs=True)
                _, _, cax = sweeps.plotSweep(data,
                        noise_sigma=noise_sigma,
                        ax=ax,
                        cbar=self.myc['cbar'][ns_idx],
                        cbar_kw=myc['cbar_kw'],
                        vmin=bump_vmin, vmax=bump_vmax,
                        annotations=ann[ns_idx],
                        **kw)


##############################################################################
# Bump formation, thresholded
class IsBumpPlotter(SweepPlotter):
    bumpThreshold = 0.95

    def __init__(self, *args, **kwargs):
        super(IsBumpPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        myc= self._get_class_config()
        sweepc = self._get_sweep_config()
        ps = self.env.ps
        iter_list = self.config['iter_list']

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            fname = (self.config['output_dir'] +
                     "bumps_isBump_sweeps{0}.pdf".format(int(noise_sigma)))
            with self.figure_and_axes(fname, sweepc) as (fig, ax):
                kw = dict(cbar=False, sigmaTitle=False)
                if ns_idx != 0:
                    kw['ylabel'] = ''
                    kw['yticks'] = False
                data = aggr.BumpFormationFilter(self.bumpThreshold, ps.bumpGamma[ns_idx],
                        iter_list, ignoreNaNs=True)
                _, _, cax = sweeps.plotSweep(data,
                        noise_sigma=noise_sigma,
                        ax=ax,
                        cbar_kw=myc['cbar_kw'],
                        vmin=bump_vmin, vmax=bump_vmax,
                        annotations=ann[ns_idx],
                        **kw)

