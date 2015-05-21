#!/usr/bin/env python
'''Grid field figures - probabilistic connections.'''
from __future__ import absolute_import, print_function

import matplotlib.ticker as ti
from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

import config_trials_weight_0 as config

parser = flagparse.FlagParser()
parser.add_flag('--grids')
parser.add_flag('--scatter')
parser.add_flag('--hist')
parser.add_flag('--examples')
parser.add_flag('--examples_cbar')
args = parser.parse_args()

env = NoiseEnvironment(user_config=config.get_config())

if args.hist or args.all:
    env.register_plotter(noisefigs.plotters.IPCHistogramsPlotter,
                         config={
                             'IPCHistogramsPlotter': {
                                 'weight_row': 0,
                             },
                         })

if args.examples or args.all:
    env.register_plotter(noisefigs.plotters.IPCExamplePlotter,
                         config={
                             'IPCExamplePlotter': {
                                 'weight_row': 0,
                             },
                         })
    env.register_plotter(noisefigs.plotters.IPCExamplePlotter,
                         config={
                             'IPCExamplePlotter': {
                                 'population_type': 'I',
                                 'weight_row': 0,
                             },
                         })

if args.examples_cbar or args.all:
    env.register_plotter(noisefigs.plotters.IPCExampleColorbarPlotter)

env.plot()
