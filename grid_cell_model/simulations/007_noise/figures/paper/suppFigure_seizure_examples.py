#!/usr/bin/env python
from __future__ import absolute_import, print_function

from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

import seizures
import config

singleDataRoot = 'simulation_data/submission/single_neuron'

parser = flagparse.FlagParser()
parser.add_flag('--raster-examples')
args = parser.parse_args()

env = NoiseEnvironment(user_config=config.get_config())

if args.raster_examples or args.all:
    env.register_plotter(seizures.RasterExamplePlotter)


env.plot()
