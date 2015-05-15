#!/usr/bin/env python
'''Grid field figures - I-surround (all other config values as in the default
E-surround configuration (main figures).'''
from __future__ import absolute_import, print_function

from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

import config_original as config

parser = flagparse.FlagParser()
parser.add_flag('--grids')
parser.add_flag('--examplesFlag')
parser.add_flag('--examples_colorbar')
parser.add_flag('--diff_sweep')
args = parser.parse_args()

env = NoiseEnvironment(user_config=config.get_config())

if args.grids or args.all:
    env.register_plotter(noisefigs.plotters.GridSweepsPlotter)

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

if args.diff_sweep or args.all:
    env.register_plotter(noisefigs.plotters.GridsDiffSweep)

env.plot()
