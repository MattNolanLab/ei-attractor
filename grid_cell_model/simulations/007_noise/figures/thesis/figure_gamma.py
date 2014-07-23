#!/usr/bin/env python
from __future__ import absolute_import, print_function

from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

import config

parser = flagparse.FlagParser()
parser.add_flag('--gammaSweep')
parser.add_flag('--threshold')
parser.add_flag('--freqHist')
parser.add_flag('--detailed_noise')
parser.add_flag('--examples')
parser.add_flag('--scatter_all')
args = parser.parse_args()


env = NoiseEnvironment(user_config=config.get_config())

if args.gammaSweep or args.all:
    env.register_plotter(noisefigs.plotters.GammaSweepsPlotter)

if args.detailed_noise or args.all:
    env.register_plotter(noisefigs.plotters.GammaDetailedNoisePlotter)

if args.examples or args.all:
    env.register_plotter(noisefigs.plotters.GammaExamplePlotter)

if args.scatter_all or args.all:
    env.register_plotter(noisefigs.plotters.GammaScatterAllPlotter)

env.plot()
