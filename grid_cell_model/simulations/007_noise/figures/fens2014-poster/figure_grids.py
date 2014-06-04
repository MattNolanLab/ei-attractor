#!/usr/bin/env python
from __future__ import absolute_import, print_function

from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

import config

singleDataRoot = 'simulation_data/submission/single_neuron'

parser = flagparse.FlagParser()
parser.add_flag('--grids')
parser.add_flag('--examplesFlag')
parser.add_flag('--Vm_examples')
parser.add_flag('--conn_func')
parser.add_flag('--example_hists')
args = parser.parse_args()

#ps = ds.getDefaultParamSpaces()

env = NoiseEnvironment(user_config=config.get_config())

if args.grids or args.all:
    env.register_plotter(noisefigs.plotters.GridSweepsPlotter)

if args.examplesFlag or args.all:
    env.register_plotter(noisefigs.plotters.GridExamplesPlotter)

if args.Vm_examples or args.all:
    env.register_plotter(noisefigs.plotters.VmExamplesPlotter)

if args.conn_func or args.all:
    env.register_plotter(noisefigs.plotters.ConnectionFunctionPlotter)

if args.example_hists or args.all:
    env.register_plotter(noisefigs.plotters.WeightExamplesHists)

env.plot()
