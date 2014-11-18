#!/usr/bin/env python
from __future__ import absolute_import, print_function

from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

import config

singleDataRoot = 'simulation_data/submission/single_neuron'

parser = flagparse.FlagParser()
parser.add_flag('--grids')
parser.add_flag('--Vm_examples')
parser.add_flag('--diff_sweep')
parser.add_flag('--conn_func')
args = parser.parse_args()

env = NoiseEnvironment(user_config=config.get_config())

if args.grids or args.all:
    env.register_plotter(noisefigs.plotters.GridSweepsPlotter)

if args.Vm_examples or args.all:
    env.register_plotter(noisefigs.plotters.VmExamplesPlotter)

if args.conn_func or args.all:
    env.register_plotter(noisefigs.plotters.ConnectionFunctionPlotter)

env.plot()
