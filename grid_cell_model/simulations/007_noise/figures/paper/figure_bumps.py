#!/usr/bin/env python
from __future__ import absolute_import, print_function

from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

import config

parser = flagparse.FlagParser()
parser.add_flag('--bumpSweep')
parser.add_flag('--bumpExamples')
parser.add_flag('--bump_running_examples')
parser.add_flag('--bump_examples_colorbar')
parser.add_flag('--detailed_noise')
parser.add_flag('--scatter_grids_fracTotal')
parser.add_flag('--scatter_gamma_pbumps_all')
parser.add_flag('--fracTotalSweepAnn')
parser.add_flag('--isBump')
parser.add_argument('--expScatter', action="store_true")
args = parser.parse_args()


env = NoiseEnvironment(user_config=config.get_config())

if args.fracTotalSweepAnn or args.all:
    env.register_plotter(noisefigs.plotters.MainBumpFormationPlotter)

if args.scatter_grids_fracTotal or args.all:
    env.register_plotter(noisefigs.plotters.MainScatterGridsBumpsPlotter)

if args.scatter_gamma_pbumps_all or args.all:
    env.register_plotter(noisefigs.plotters.GammaScatterPBumpsAllPlotter)

if args.isBump or args.all:
    env.register_plotter(noisefigs.plotters.MainIsBumpPlotter)

if args.bumpSweep or args.all:
    env.register_plotter(noisefigs.plotters.BumpSigmaSweepPlotter)

if args.bumpExamples or args.all:
    env.register_plotter(noisefigs.plotters.BumpExamplePlotter)

if args.bump_running_examples or args.all:
    env.register_plotter(noisefigs.plotters.IsBumpExamplePlotter)

if args.bump_examples_colorbar or args.all:
    env.register_plotter(noisefigs.plotters.BumpExampleColorbarPlotter)

if args.detailed_noise or args.all:
    env.register_plotter(noisefigs.plotters.BumpSigmaDetailedNoisePlotter)

env.plot()

