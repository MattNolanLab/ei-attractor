#!/usr/bin/env python
'''Plot figures for 1D parameter explorations (gamma).'''
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
parser.add_flag('--gamma_g_II')
args = parser.parse_args()


if args.gamma_g_II or args.all:
    new_config = deepcopy(config.get_config())
    new_config.update({
        'grids_data_root': None,
        'bump_data_root': os.path.join('simulation_data', 'submission',
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
        noisefigs.plotters.Generic1DGammaPlotter,
        config={
            'Generic1DGammaPlotter' : {
                'what': 'acVal',
                'fname' : 'gamma_power_g_II_{ns}.pdf',
                'normalize_ticks': False,
                'ylabel': '$1^{st}$ autocorrelation peak',
                'xlim': (-10, 410),
                'ylim': (-0.04, 0.58),
            },
        }
    )
    env.register_class(
        noisefigs.plotters.Generic1DGammaPlotter,
        config={
            'Generic1DGammaPlotter' : {
                'what': 'freq',
                'fname' : 'gamma_freq_g_II_{ns}.pdf',
                'normalize_ticks': False,
                'ylabel': 'Frequency (Hz)',
                'xlim': (-10, 410),
                'ylim': (22, 175),
            },
        }
    )
    env.plot()


