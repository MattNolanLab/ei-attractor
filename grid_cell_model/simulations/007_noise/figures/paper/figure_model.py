#!/usr/bin/env python
'''Figures for the model description.'''
from __future__ import absolute_import, print_function

from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

import config

parser = flagparse.FlagParser()
parser.add_flag('--Vm_examples')
parser.add_flag('--conn_func')
parser.add_flag('--example_hists')
parser.add_flag('--weight_plots')
args = parser.parse_args()

env = NoiseEnvironment(user_config=config.get_config())

if args.Vm_examples or args.all:
    env.register_plotter(noisefigs.plotters.VmExamplesPlotter)

if args.conn_func or args.all:
    env.register_plotter(noisefigs.plotters.ConnectionFunctionPlotter)

if args.example_hists or args.all:
    env.register_plotter(noisefigs.plotters.WeightExamplesHists)

if args.weight_plots or args.all:
    env.register_plotter(noisefigs.plotters.WeightOutE2IPlotter)
    env.register_plotter(noisefigs.plotters.WeightOutI2EPlotter)
    env.register_plotter(noisefigs.plotters.WeightInE2IPlotter)
    env.register_plotter(noisefigs.plotters.WeightInI2EPlotter)

env.plot()
