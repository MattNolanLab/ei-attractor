#!/usr/bin/env python
from __future__ import absolute_import, print_function, division
from copy import deepcopy
from os.path import join

import matplotlib.ticker as ti
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
        'bump_data_root': join('simulation_data', 'submission',
                               'ee_connections_ei_flat',
                               'g_EE_total_vs_pEE_sigma', 'gamma_bump'),
        'vel_data_root':  None,
        'const_pos_data_root': None,
        'singleDataRoot': None,

        'noise_sigmas': [0, 150, 300],
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
                'xlabel': '$\sigma_{E{\\rightarrow}E}$',
                'ylabel': '$g_{E{\\rightarrow}E}$',
                'bbox': (.2, .17, .85, .9),

                'cbar': [0, 0, 1],
            },
        })
    env.plot()
