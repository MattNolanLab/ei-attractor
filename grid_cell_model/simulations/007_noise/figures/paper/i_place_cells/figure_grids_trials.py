#!/usr/bin/env python
from __future__ import absolute_import, print_function

import matplotlib.ticker as ti
from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

import config_trials as config

parser = flagparse.FlagParser()
#parser.add_flag('--grids')
parser.add_flag('--scatter')
parser.add_flag('--hist')
parser.add_flag('--examples')
parser.add_flag('--examples_cbar')
args = parser.parse_args()

env = NoiseEnvironment(user_config=config.get_config())

#if args.grids or args.all:
#    env.register_plotter(noisefigs.plotters.IPCGridSweepsPlotter)
#    env.register_plotter(noisefigs.plotters.IPCGridSweepsPlotter,
#                         config={
#                             'IPCGridSweepsPlotter': {
#                                 'population_type': 'I',
#                             }
#                         })

if args.scatter or args.all:
    env.register_plotter(noisefigs.plotters.IPCScatterPlotter)

if args.hist or args.all:
    env.register_plotter(noisefigs.plotters.IPCHistogramsPlotter)

if args.examples or args.all:
    env.register_plotter(noisefigs.plotters.IPCExamplePlotter)
    env.register_plotter(noisefigs.plotters.IPCExamplePlotter,
                         config={
                             'IPCExamplePlotter': {
                                 'population_type': 'I'
                             },
                         })

if args.examples_cbar or args.all:
    env.register_plotter(noisefigs.plotters.IPCExampleColorbarPlotter)

env.plot()
