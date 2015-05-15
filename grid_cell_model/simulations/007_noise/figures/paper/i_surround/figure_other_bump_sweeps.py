#!/usr/bin/env python
from __future__ import absolute_import, print_function
from copy import deepcopy

from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

import config_original as config

parser = flagparse.FlagParser()
parser.add_flag('--pbumps_Iext_e')
args = parser.parse_args()

if args.pbumps_Iext_e or args.all:
    new_config = deepcopy(config.get_config())
    new_config.update({
        'grids_data_root': None,
        'bump_data_root': ('simulation_data/submission/i_surround/'
                           'gamma_bump_Iext_e_theta_vs_Iext_e_const'),
        'vel_data_root':  None,
        'const_pos_data_root': None,
        'singleDataRoot': None,

        'noise_sigmas': [0, 150, 300],
        'even_shape': None,
        'output_dir': 'panels',
    })
    env = NoiseEnvironment(user_config=new_config)
    env.register_class(
        noisefigs.plotters.Generic2DPBumpPlotter,
        config={
            'Generic2DPBumpPlotter': {
                'fname': "bumps_pbumps_Iext_e_theta_vs_Iext_e_const_{ns}.pdf",
                'normalize_ticks': (False, False),  # (Y, X)
                'normalize_type': (None, None),
                'xlabel': 'Constant amplitude\n(pA)',
                'ylabel': '$\\theta$ amplitude (pA)',
                'bbox': (.15, .2, .85, .9),
            },
        })
    env.plot()
