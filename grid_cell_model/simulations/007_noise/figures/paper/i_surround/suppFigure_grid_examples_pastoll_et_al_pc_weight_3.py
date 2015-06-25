#!/usr/bin/env python
'''Examples of grid firing fields in the Pastoll et al. configuration with PC
max. weight set to 3 nS.'''
from __future__ import absolute_import, print_function

from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

import config_pastoll_pc_weight_3 as config

parser = flagparse.FlagParser()
parser.add_flag('--grids_examples')
args = parser.parse_args()

env = NoiseEnvironment(user_config=config.get_config())

if args.grids_examples or args.all:
    env.register_plotter(noisefigs.plotters.GridExampleRectPlotter)

env.plot()
