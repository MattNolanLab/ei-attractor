#!/usr/bin/env python
from __future__ import absolute_import, print_function
from copy import deepcopy

from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

import config

parser = flagparse.FlagParser()
parser.add_flag('--pbumps_gEE_EE_sigma')
args = parser.parse_args()

if args.pbumps_gEE_EE_sigma or args.all:
    new_config = deepcopy(config.get_config())
    new_config.update({
        'grids_data_root': None,
        'bump_data_root': ('simulation_data/g_EE_total_vs_pEE_sigma/'
                           'gamma_bump'),
        'vel_data_root':  None,
        'const_pos_data_root': None,
        'singleDataRoot': None,

        'noise_sigmas': [150],
        'even_shape': None,
    })
    env = NoiseEnvironment(user_config=new_config)
    env.register_class(
        noisefigs.plotters.Generic2DPBumpPlotter,
        config={
            'Generic2DPBumpPlotter': {
                'fname': "bumps_Pbumps_gEE_pEE_sigma_{ns}.pdf",
                'normalize_ticks': (True, False),  # (Y, X)
                'normalize_type': ('E', None),
                'xlabel': '$\sigma_{EE}$',
                'ylabel': '$g_E$',
                'bbox': (.2, .17, .85, .9),
            },
        })
    env.plot()
