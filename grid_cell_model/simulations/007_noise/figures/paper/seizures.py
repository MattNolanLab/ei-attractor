from __future__ import absolute_import, print_function, division

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
import matplotlib.gridspec as gridspec

from noisefigs.EI_plotting import sweeps, rasters
from noisefigs.EI_plotting import aggregate as aggr
from noisefigs.plotters.base import SweepPlotter

__all__ = [
    'RasterExamplePlotter'
]

##############################################################################
# Seizure measure - max firing rate for the whole simulation
maxFR_vmin = 0
maxFR_vmax = 500.

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
        color='white')
ann300_0 = dict(
        txt='d',
        rc=(20, 25),
        xytext_offset=(-1.5, 1),
        color='white')
ann0_2 = dict(
        txt='e',
        rc=(15, 5),
        xytext_offset=(.5, 2),
        color='black')
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

class RasterExamplePlotter(SweepPlotter):
    def __init__(self, *args, **kwargs):
        super(RasterExamplePlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        myc= self._get_class_config()
        sweepc = self._get_sweep_config()
        ps = self.env.ps
        iter_list = self.config['iter_list']
        sl, sb, sr, st = (.1, .73, .4, .95)
        tLimits = [2e3, 3e3]

        for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
            ann_noise = ann[ns_idx]
            data = aggr.MaxPopulationFR(ps.bumpGamma[ns_idx], iter_list,
                    ignoreNaNs=True, normalizeTicks=True)
            for annotation in ann_noise:
                r, c = annotation['rc']
                fig = plt.figure(figsize=(8.3, 8.3))

                # Sweep
                ax_sweep = fig.add_axes(Bbox.from_extents(sl, sb, sr, st))
                _, _, cax = sweeps.plotSweep(data,
                        noise_sigma=noise_sigma,
                        ax=ax_sweep,
                        cbar=True, cbar_kw=myc['cbar_kw'],
                        vmin=maxFR_vmin, vmax=maxFR_vmax,
                        annotations=[annotation])

                gs = gridspec.GridSpec(3, 1, height_ratios=(2.5, 1, 1))

                # EI Raster
                ax_raster = fig.add_subplot(gs[0, 0])
                rasters.EIRaster(ps.bumpGamma[ns_idx], 
                        noise_sigma=noise_sigma,
                        spaceType='bump',
                        r=r, c=c,
                        ylabelPos=self.myc['ylabelPos'],
                        tLimits=tLimits,
                        markersize=self.config['scale_factor']*self.myc['markersize'],
                        yticks=True,
                        sigmaTitle=False,
                        ann_EI=True,
                        scaleBar=125, scaleX=.85, scaleY=-.1,
                        scaleTextYOffset=.03, scaleHeight=.01,
                        rasterized=True)

                # EI rates
                ax_erates = fig.add_subplot(gs[1, 0])
                ax_irates = fig.add_subplot(gs[2, 0])
                rasters.plotAvgFiringRate(ps.bumpGamma[ns_idx],
                        spaceType='bump',
                        noise_sigma=noise_sigma,
                        popType='E',
                        r=r, c=c,
                        ylabelPos=self.myc['ylabelPos'],
                        color='red',
                        tLimits=tLimits,
                        ax=ax_erates)
                rasters.plotAvgFiringRate(ps.bumpGamma[ns_idx],
                        spaceType='bump',
                        noise_sigma=noise_sigma,
                        popType='I',
                        r=r, c=c,
                        ylabelPos=self.myc['ylabelPos'],
                        color='blue',
                        tLimits=tLimits,
                        ax=ax_irates)

                # Save
                gs.update(left=.12, bottom=.05, right=.95, top=.65, hspace=.2)

                fname = (self.config['output_dir'] +
                         "/raster_examples_{0}pA_{1}_{2}.pdf".format(
                             int(noise_sigma), r, c))

                fig.savefig(fname, dpi=300, transparent=True)
                plt.close(fig)
