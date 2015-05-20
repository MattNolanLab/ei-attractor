#!/usr/bin/env python
'''Grid field figures - probabilistic connections.'''
from __future__ import absolute_import, print_function

import matplotlib.ticker as ti
from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

import config

parser = flagparse.FlagParser()
parser.add_flag('--grids')
args = parser.parse_args()

env = NoiseEnvironment(user_config=config.get_config())

if args.grids or args.all:
    env.register_plotter(noisefigs.plotters.GenericGridSweepsPlotter)
    env.register_plotter(noisefigs.plotters.GenericGridSweepsPlotter,
                         config={
                             'GenericGridSweepsPlotter': {
                                 'population_type': 'I',
                             }
                         })

env.plot()
