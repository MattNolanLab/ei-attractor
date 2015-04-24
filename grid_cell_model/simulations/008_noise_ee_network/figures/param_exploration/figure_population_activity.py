#!/usr/bin/env python
from __future__ import absolute_import, print_function

import os

import matplotlib
matplotlib.use('Agg')
from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import MplEnvironment
from noisefigs.plotters.base import SeparateMultipageSaver

import config

parser = flagparse.FlagParser()
parser.add_flag('--param_exploration')
parser.add_argument('--data_root', type=str, help='Data root directory',
                    default='simulation_data/network_test/150pA')
args = parser.parse_args()

env = MplEnvironment(config=config.get_config())


if args.param_exploration or args.all:
    for file_name in os.listdir(args.data_root):
        if 'job' in file_name:
            print(file_name)
            env.register_class(
                noisefigs.plotters.PopulationActivityPlotter,
                config={
                    'data_root'     : args.data_root,
                    'data_file_name': file_name,

                    'PopulationActivityPlotter': {
                        'fname_prefix': 'test_%s_' % file_name,
                        'fig_size': (4, 6),
                        't_limits': (0, 2.5e3),
                        'fig_saver': SeparateMultipageSaver(None, 'pdf'),
                    },
                })

env.plot()
