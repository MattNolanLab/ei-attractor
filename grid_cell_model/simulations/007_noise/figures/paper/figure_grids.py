#!/usr/bin/env python
from __future__ import absolute_import, print_function

import matplotlib.ticker as ti
from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

import config

parser = flagparse.FlagParser()
parser.add_flag('--grids')
parser.add_flag('--spatial_info')
parser.add_flag('--sparsity')
parser.add_flag('--spatial_info_stats')
parser.add_flag('--examplesFlag')
parser.add_flag('--examples_colorbar')
parser.add_flag('--correlation_angles')
parser.add_flag('--detailed_noise')
parser.add_flag('--diff_sweep')
parser.add_flag('--grids_pbumps_prob')
parser.add_flag('--high_gscore_frac')
args = parser.parse_args()

env = NoiseEnvironment(user_config=config.get_config())

if args.grids or args.all:
    env.register_plotter(noisefigs.plotters.GridSweepsPlotter)
    env.register_plotter(noisefigs.plotters.GridSweepsPlotter,
                         config={
                             'GridSweepsPlotter': {
                                 'population_type': 'I',
                                 'plot_contours': [1, 1, 1],
                                 'vmin': -0.32,
                                 'vmax': 0.774,
                                 'cbar_kw': {
                                     'ticks': ti.MultipleLocator(0.2),
                                 },
                             }
                         })

if args.spatial_info or args.all:
    env.register_plotter(noisefigs.plotters.SpatialInfoPlotter)
    env.register_plotter(noisefigs.plotters.SpatialInfoPlotter,
                         config={
                             'SpatialInfoPlotter': {
                                 'population_type': 'I',
                                 'cbar': [0, 0, 1],
                                 'sigma_title': False,
                                 'xlabel': [None, None, None],
                                 'xticks': [True, True, True],
                             },
                         })

if args.spatial_info_stats or args.all:
    env.register_plotter(noisefigs.plotters.SpatialInfoStats)
    env.register_plotter(noisefigs.plotters.SpatialSparsityStats)

if args.sparsity or args.all:
    env.register_plotter(noisefigs.plotters.SpatialSparsityPlotter)
    env.register_plotter(noisefigs.plotters.SpatialSparsityPlotter,
                         config={
                             'SpatialSparsityPlotter': {
                                 'population_type': 'I',
                                 'cbar': [0, 0, 1],
                                 'sigma_title': False,
                                 'xlabel': [None, None, None],
                                 'xticks': [True, True, True],
                             },
                         })

if args.examplesFlag or args.all:
    env.register_plotter(noisefigs.plotters.GridExamplesPlotter)
    env.register_plotter(noisefigs.plotters.GridExamplesPlotter,
                         config={
                            'GridExamplesPlotter': {
                                'population_type': 'I'
                            },
                         })

if args.examples_colorbar or args.all:
    env.register_plotter(noisefigs.plotters.GridExampleColorbarPlotter)

if args.correlation_angles or args.all:
    env.register_plotter(noisefigs.plotters.GridnessCorrelationPlotter)
    env.register_plotter(noisefigs.plotters.GridnessCorrelationPlotter,
                         config={
                             'GridnessCorrelationPlotter': {
                                 'population_type': 'I',
                             },
                         })

if args.detailed_noise or args.all:
    env.register_plotter(noisefigs.plotters.GridDetailedNoisePlotter)

if args.diff_sweep or args.all:
    env.register_plotter(noisefigs.plotters.GridsDiffSweep)

if args.grids_pbumps_prob or args.all:
    env.register_plotter(noisefigs.plotters.GridsPBumpsProbabilityPlotter)

if args.high_gscore_frac or args.all:
    env.register_plotter(noisefigs.plotters.HighGridScoreFraction)


env.plot()
