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
        if 'rasters' in file_name:
            print(file_name)
            env.register_class(
                noisefigs.plotters.PopulationActivityPlotter,
                config={
                    'data_root'     : args.data_root,
                    'data_file_name': file_name,

                    'PopulationActivityPlotter': {
                        'fig_size': (6, 5),
                        't_limits': (5e3, 6e3),

                        'raster_rect': (.125, 0.35, 0.99, 0.97),
                        'snapshot_tstep': 1,
                        'e_snapshots_rect': (.125, .15, 0.99, 0.25),
                        'i_snapshots_rect': (.125, .02, 0.99, 0.12),

                        'fname_prefix': 'thesis_rasters_%s_' % file_name,
                        'fig_saver': SeparateMultipageSaver(None, 'pdf'),
                        'reshape_senders' : False,
                        'max_e_rate': False,
                        'max_i_rate': False,
                        'scale_bar': 250,
                        'scale_x': .73,
                        'ann_ei': False,
                        'y_label_pos': -.1
                    },
                })

env.plot()
