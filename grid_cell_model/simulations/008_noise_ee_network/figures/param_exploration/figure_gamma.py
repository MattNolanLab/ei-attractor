#!/usr/bin/env python
from __future__ import absolute_import, print_function
from copy import deepcopy

import matplotlib.ticker as ti
from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

import config

parser = flagparse.FlagParser()
parser.add_flag('--gamma_gEE_EE_sigma')
args = parser.parse_args()

if args.gamma_gEE_EE_sigma or args.all:
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
        noisefigs.plotters.GenericGammaPlotter,
        config={
            'GenericGammaPlotter': {
                'what': 'acVal',
                'fname': "gamma_power_gEE_pEE_sigma_{ns}.pdf",
                'normalize_ticks': (True, False),  # (Y, X)
                'normalize_type': ('E', None),
                'xlabel': '$\sigma_{E{\\rightarrow}E}$',
                'ylabel': '$g_{E{\\rightarrow}E}$',
                'bbox': (.2, .17, .85, .9),
                'vmin': .0,
                'vmax': .505,
                'cbar_kw': dict(
                    ticks = ti.MultipleLocator(0.1),
                )
                },
        })
    env.register_class(
        noisefigs.plotters.GenericGammaPlotter,
        config={
            'GenericGammaPlotter': {
                'what': 'freq',
                'fname': "gamma_freq_gEE_pEE_sigma_{ns}.pdf",
                'normalize_ticks': (True, False),  # (Y, X)
                'normalize_type': ('E', None),
                'xlabel': '$\sigma_{E{\\rightarrow}E}$',
                'ylabel': '',
                'yticks': [False],
                'bbox': (.2, .17, .85, .9),
                'vmin': 30,
                'vmax': 72,
                'cbar_kw': dict(
                    label='Frequency (Hz)',
                    ticks=ti.MultipleLocator(10),
                )
                },
        })
    env.plot()
