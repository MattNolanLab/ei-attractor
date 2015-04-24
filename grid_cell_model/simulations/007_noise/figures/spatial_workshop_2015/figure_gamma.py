#!/usr/bin/env python
from __future__ import absolute_import, print_function

import copy

import matplotlib
from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

import config

parser = flagparse.FlagParser()
parser.add_flag('--gammaSweep')
parser.add_flag('--examples')
args = parser.parse_args()


env = NoiseEnvironment(user_config=config.get_config())

if args.gammaSweep or args.all:
    env.register_plotter(noisefigs.plotters.GammaSweepsPlotter)

if args.examples or args.all:
    env.register_plotter(noisefigs.plotters.GammaExamplePlotter)


env.plot()
