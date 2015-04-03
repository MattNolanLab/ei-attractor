#!/usr/bin/env python
from __future__ import absolute_import, print_function

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.transforms import Bbox
import numpy as np
from grid_cell_model.submitting import flagparse
from grid_cell_model.data_storage import DataStorage
from grid_cell_model.analysis.spikes import PopulationSpikes
import noisefigs
from noisefigs.env import NoiseEnvironment
from noisefigs.plotters.base import FigurePlotter
from noisefigs.EI_plotting import rasters

import config

class RasterExamples(FigurePlotter):
    dt = .1  # ms
    freq = 8. # Hz
    const = .4 # Fraction of max. theta

    def __init__(self, *args, **kwargs):
        super(RasterExamples, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        myc= self._get_class_config()
        tLimits = [2e3, 2.5e3]


        fig = self._get_final_fig(self.myc['fig_size'])
        gs = gridspec.GridSpec(3, 1, height_ratios=(2.5, 1, 1))

        data = DataStorage.open(self.myc['data_file'], 'r')
        trial_data = data['trials'][0]
        events_e = trial_data['spikeMon_e']['events']
        ESpikes = PopulationSpikes(trial_data['net_attr']['net_Ne'],
                                   events_e['senders'],
                                   events_e['times'] * 1e3)
        events_i = trial_data['spikeMon_i']['events']
        ISpikes = PopulationSpikes(trial_data['net_attr']['net_Ni'],
                                   events_i['senders'],
                                   events_i['times'] * 1e3)

        # EI Raster
        ax_raster = fig.add_subplot(gs[0, 0])
        rasters.plotEIRaster(
            ESpikes, ISpikes,
            ylabelPos=self.myc['ylabelPos'],
            tLimits=tLimits,
            markersize=self.config['scale_factor']*self.myc['markersize'],
            yticks=True,
            sigmaTitle=False,
            ann_EI=False,
            scaleBar=125, scaleX=.7, scaleY=-.05, scaleText='ms',
            scaleTextYOffset=.03, scaleHeight=.01,
            rasterized=False,
            reshape_senders=False)

        # EI rates
        ax_erates = fig.add_subplot(gs[1, 0])
        ax_irates = fig.add_subplot(gs[2, 0])
        rasters.plot_avg_firing_rate_spikes(ESpikes,
                                            ylabelPos=self.myc['ylabelPos'],
                                            color='red',
                                            tLimits=tLimits,
                                            ax=ax_erates,
                                            dt=.5,
                                            winLen=2.)
        rasters.plot_avg_firing_rate_spikes(ISpikes,
                                            ylabelPos=self.myc['ylabelPos'],
                                            color='blue',
                                            tLimits=tLimits,
                                            ax=ax_irates,
                                            dt=.5,
                                            winLen=2.)

        gsl = .12
        gsb = .02
        gsr = .95
        gst = .95
        #fig.text(0.01, gst, 'B', size=16, weight='bold',
        #         va='bottom', ha='left')
        gs.update(left=gsl, bottom=gsb, right=gsr, top=gst, hspace=.2)

        ax_theta = fig.add_axes(Bbox.from_extents(gsl, gst - .015,
                                                    gsr, gst + .01))
        t = np.arange(tLimits[0], tLimits[1]+self.dt, self.dt)
        theta = self.const + .5 * (1. +
                np.cos(2*np.pi*self.freq*1e-3*t - np.pi)) * (1 - self.const)
        ax_theta.fill_between(t, theta, edgecolor='None',
                                color=self.myc['theta_color'])
        ax_theta.set_xlim([tLimits[0], tLimits[1]])
        ax_theta.set_ylim(-.02, 1.02)
        ax_theta.axis('off')

        plt.savefig(self.get_fname('pastoll_et_al_rasters.pdf'), dpi=300,
                    transparent=True)
        plt.close(fig)


parser = flagparse.FlagParser()
parser.add_flag('--raster_examples')
args = parser.parse_args()

env = NoiseEnvironment(user_config=config.get_config())

if args.raster_examples or args.all:
    env.register_plotter(RasterExamples)


env.plot()
