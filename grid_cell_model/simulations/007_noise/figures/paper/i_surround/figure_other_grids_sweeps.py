#!/usr/bin/env python
from __future__ import absolute_import, print_function
from copy import deepcopy
import os.path

from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

import config_original as config
import config_pastoll

parser = flagparse.FlagParser()
parser.add_flag('--grids_pastoll_pc_weight_3')
parser.add_flag('--grids_pastoll_pc_weight')
parser.add_flag('--grids_pastoll_pc_weight_no_vel')
args = parser.parse_args()

if args.grids_pastoll_pc_weight_3 or args.all:
    new_config = deepcopy(config_pastoll.get_config())
    new_config.update({
        'bump_data_root': None,
        'grids_data_root': os.path.join('simulation_data',
                                        'i_surround',
                                        'pastoll_et_al',
                                        'grids_pc_weight_3'),
        'vel_data_root':  None,
        'const_pos_data_root': None,
        'singleDataRoot': None,

        'noise_sigmas': [150],
        'even_shape': None,
        'output_dir': 'panels',
    })
    env = NoiseEnvironment(user_config=new_config)
    env.register_class(noisefigs.plotters.GridSweepsPlotter,
                       config={

                           'GridSweepsPlotter': {
                               'fname_prefix': "pastoll_pc_weight_3_",
                               'vmin': -0.4286,
                               'vmax': 0.8582,
                               'ann': None,
                           },
                       })
    env.plot()


if args.grids_pastoll_pc_weight or args.all:
    new_config = deepcopy(config.get_config())
    new_config.update({
        'bump_data_root': None,
        'grids_data_root': os.path.join('simulation_data',
                                        'i_surround',
                                        'pastoll_et_al_pc_max_rate_vs_weight',
                                        'grids'),
        'vel_data_root':  None,
        'const_pos_data_root': None,
        'singleDataRoot': None,

        'noise_sigmas': [150],
        'even_shape': None,
        'output_dir': 'panels',
    })
    env = NoiseEnvironment(user_config=new_config)
    env.register_class(
        noisefigs.plotters.GenericGridSweepsPlotter,
        config={
            'GenericGridSweepsPlotter': {
                'fname_prefix': "pastoll_pc_weight_",
                'normalize_ticks': (False, False),  # (Y, X)
                'normalize_type': (None, None),
                'xlabel': '',
                'ylabel': '',
                'xticks': [0, 0, 0],
                'yticks': [1, 0, 0],

                'bbox': (.15, .2, .85, .9),
                'cbar': [1, 1, 1],
            },
        })
    env.plot()

if args.grids_pastoll_pc_weight_no_vel or args.all:
    new_config = deepcopy(config.get_config())
    new_config.update({
        'bump_data_root': None,
        'grids_data_root': os.path.join('simulation_data',
                                        'i_surround',
                                        'pastoll_et_al_pc_max_rate_vs_weight',
                                        'grids_no_velocity'),
        'vel_data_root':  None,
        'const_pos_data_root': None,
        'singleDataRoot': None,

        'noise_sigmas': [150],
        'even_shape': None,
        'output_dir': 'panels',
    })
    env = NoiseEnvironment(user_config=new_config)
    env.register_class(
        noisefigs.plotters.GenericGridSweepsPlotter,
        config={
            'GenericGridSweepsPlotter': {
                'fname_prefix': "pastoll_pc_weight_no_velocity_",
                'normalize_ticks': (False, False),  # (Y, X)
                'normalize_type': (None, None),
                'xlabel': 'PC max. weight (nS)',
                'ylabel': '',
                'yticks': [1, 0, 0],
                'sigmaTitle': False,

                'bbox': (.15, .2, .85, .9),
                'cbar': [1, 1, 1],
            },
        })
    env.plot()
