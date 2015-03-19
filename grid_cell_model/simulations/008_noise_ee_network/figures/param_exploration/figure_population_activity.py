#!/usr/bin/env python
from __future__ import absolute_import, print_function

import os

from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import MplEnvironment

import config

parser = flagparse.FlagParser()
parser.add_flag('--param_exploration')
args = parser.parse_args()

env = MplEnvironment(config=config.get_config())

data_root = 'simulation_data/network_test/150pA'

if args.param_exploration or args.all:
    for file_name in os.listdir(data_root):
        if 'job' in file_name:
            print(file_name)
            env.register_class(
                noisefigs.plotters.PopulationActivityPlotter,
                config={
                    'data_root'     : data_root,
                    'data_file_name': file_name,

                    'PopulationActivityPlotter': {
                        'fname_prefix': 'test_%s_' % file_name,
                        'fig_size': (4, 6),
                        't_limits': (0, 2.5e3),
                    },
                })

env.plot()
