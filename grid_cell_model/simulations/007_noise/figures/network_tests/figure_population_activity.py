#!/usr/bin/env python
from __future__ import absolute_import, print_function

from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import MplEnvironment

import config

parser = flagparse.FlagParser()
parser.add_flag('--rasters_and_bumps')
parser.add_flag('--rasters_and_bumps_test')
args = parser.parse_args()

env = MplEnvironment(config=config.get_config())

if args.rasters_and_bumps or args.all:
    env.register_class(noisefigs.plotters.PopulationActivityPlotter)

if args.rasters_and_bumps_test or args.all:
    env.register_class(
        noisefigs.plotters.PopulationActivityPlotter,
        config={
            'data_root'     : 'simulation_data_local/network_test/150pA',

            'PopulationActivityPlotter': {
                'fname_prefix': 'test_',
                'fig_size': (4, 6),
                't_limits': (0, 2.5e3),
            },
        })

env.plot()
