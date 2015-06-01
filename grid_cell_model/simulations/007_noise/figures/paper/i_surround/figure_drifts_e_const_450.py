#!/usr/bin/env python
'''I-surround bump drift plots with the original E-surround settings but with
an increased constant drive to E cells.'''
from __future__ import absolute_import, print_function

from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

import config_e_const_450 as config

parser = flagparse.FlagParser()
parser.add_flag('--drift_sweep')
parser.add_flag('--theta_signal')
args = parser.parse_args()

env = NoiseEnvironment(user_config=config.get_config())

if args.drift_sweep or args.all:
    env.register_plotter(noisefigs.plotters.BumpDriftAtTimePlotter)

if args.theta_signal or args.all:
    env.register_plotter(noisefigs.plotters.ThetaSignalPlotter)

env.plot()

