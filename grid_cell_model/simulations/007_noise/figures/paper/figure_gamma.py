#!/usr/bin/env python
from __future__ import absolute_import, print_function

import copy

import matplotlib
from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

import config

parser = flagparse.FlagParser()
parser.add_flag('--gammaSweep')
parser.add_flag('--threshold')
parser.add_flag('--freqHist')
parser.add_flag('--detailed_noise')
parser.add_flag('--examples')
parser.add_flag('--scatter_all')
parser.add_flag('--gamma_pbumps_prob')
parser.add_flag('--gamma_grids_prob')
parser.add_flag('--gridness_stats')
args = parser.parse_args()


env = NoiseEnvironment(user_config=config.get_config())

if args.gammaSweep or args.all:
    env.register_plotter(noisefigs.plotters.GammaSweepsPlotter)

if args.detailed_noise or args.all:
    env.register_plotter(noisefigs.plotters.GammaDetailedNoisePlotter)

if args.examples or args.all:
    env.register_plotter(noisefigs.plotters.GammaExamplePlotter)

if args.scatter_all or args.all:
    env.register_plotter(noisefigs.plotters.GammaScatterAllPlotter)
    env.register_plotter(noisefigs.plotters.GammaFreqGridsScatterAllPlotter)

if args.gamma_pbumps_prob or args.all:
    env.register_plotter(noisefigs.plotters.GammaPBumpsProbabilityPlotter)
    env.register_plotter(noisefigs.plotters.GammaFreqPBumpsProbabilityPlotter)

if args.gamma_grids_prob or args.all:
    env.register_plotter(noisefigs.plotters.GammaGridsProbabilityPlotter)
    env.register_plotter(noisefigs.plotters.GammaFreqGridsProbabilityPlotter)

if args.gridness_stats or args.all:
    # Paper figures
    stats_config = {
        'GammaSweepsPlotter': {
            'filter_with_gridness': True,
            'gridness_threshold': .5,
            'fname_prefix': 'gridness_filt_',
            'ann': None,
            'annF': None,
            'plot_grid_contours': [0, 0, 0],
            'AC_vmin': 0.03,
            'AC_vmax': 0.594,
            'AC_cbar_kw': dict(
                ticks = matplotlib.ticker.MultipleLocator(0.1),
            ),
            'F_vmin': 31.2,
            'F_vmax': 103,
            'F_cbar_kw': dict(
                ticks = matplotlib.ticker.MultipleLocator(25),
                extend = 'neither',
            ),
        }
    }
    env.register_plotter(
        noisefigs.plotters.GammaSweepsPlotter,
        config=stats_config,
    )


env.plot()
