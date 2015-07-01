#!/usr/bin/env python
from __future__ import absolute_import, print_function

from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

import config

singleDataRoot = 'simulation_data/single_neuron'

parser = flagparse.FlagParser()
parser.add_flag('--grids')
parser.add_flag('--examplesFlag')
parser.add_flag('--examples_colorbar')
parser.add_flag('--detailed_noise')
parser.add_flag('--diff_sweep')
parser.add_flag('--conn_func')
parser.add_flag('--example_hists')
parser.add_flag('--grids_pbumps_prob')
parser.add_flag('--high_gscore_frac')
args = parser.parse_args()

env = NoiseEnvironment(user_config=config.get_config())

if args.grids or args.all:
    env.register_plotter(noisefigs.plotters.GridSweepsPlotter)

if args.examplesFlag or args.all:
    env.register_plotter(noisefigs.plotters.GridExamplesPlotter)

if args.examples_colorbar or args.all:
    env.register_plotter(noisefigs.plotters.GridExampleColorbarPlotter)

if args.detailed_noise or args.all:
    env.register_plotter(noisefigs.plotters.GridDetailedNoisePlotter)

if args.diff_sweep or args.all:
    env.register_plotter(noisefigs.plotters.GridsDiffSweep)

#if args.conn_func or args.all:
#    env.register_plotter(noisefigs.plotters.ConnectionFunctionPlotter)

#if args.example_hists or args.all:
#    env.register_plotter(noisefigs.plotters.WeightExamplesHists)

if args.grids_pbumps_prob or args.all:
    env.register_plotter(noisefigs.plotters.GridsPBumpsProbabilityPlotter)

if args.high_gscore_frac or args.all:
    env.register_plotter(noisefigs.plotters.HighGridScoreFraction)

env.plot()
