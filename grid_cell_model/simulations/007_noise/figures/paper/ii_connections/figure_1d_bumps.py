#!/usr/bin/env python
'''Plot figures for 1D parameter explorations.'''
from __future__ import absolute_import, print_function

from copy import deepcopy
import os.path
import matplotlib; matplotlib.use('Agg')
from grid_cell_model.submitting import flagparse
from grid_cell_model.parameters.param_space import JobTrialSpace1D
import noisefigs
from noisefigs.env import NoiseEnvironment

import config

parser = flagparse.FlagParser()
parser.add_flag('--bumps_g_II')
args = parser.parse_args()


if args.bumps_g_II or args.all:
    new_config = deepcopy(config.get_config())
    new_config.update({
        'grids_data_root': None,
        'bump_data_root': os.path.join('simulation_data',
                                       'ii_connections',
                                       'g_II_total_sweep',
                                       'gamma_bump'),
        'vel_data_root':  None,
        'const_pos_data_root': None,
        'singleDataRoot': None,

        'noise_sigmas': [150],
        'even_shape': None,
    })
    env = NoiseEnvironment(user_config=new_config, space_cls=JobTrialSpace1D)
    env.register_class(
        noisefigs.plotters.Generic1DPBumpPlotter,
        config={
            'Generic1DPBumpPlotter' : {
                'fname' : 'bumps_Pbumps_g_II_{ns}.pdf',
                'normalize_ticks': False,
                'xlim': (-10, 410),
                'ylim': (-0.05, 1.05),
            },
        }
    )
    env.plot()


