#!/usr/bin/env python
from __future__ import absolute_import, print_function

from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

import config

parser = flagparse.FlagParser()
parser.add_flag('--fracTotalSweepAnn')
args = parser.parse_args()


env = NoiseEnvironment(user_config=config.get_config())

if args.fracTotalSweepAnn or args.all:
    env.register_plotter(noisefigs.plotters.GEProfileWidthBumpPlotter)

env.plot()

