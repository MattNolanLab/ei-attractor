#!/usr/bin/env python
from __future__ import absolute_import, print_function

from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import MplEnvironment

import config

parser = flagparse.FlagParser()
parser.add_flag('--rasters_and_bumps')
args = parser.parse_args()

env = MplEnvironment(config=config.get_config())

if args.rasters_and_bumps or args.all:
    env.register_class(noisefigs.plotters.PopulationActivityPlotter)

env.plot()
