#!/usr/bin/env python
from __future__ import absolute_import, print_function
from copy import deepcopy

import matplotlib.ticker as ti
from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

import config_other as config

parser = flagparse.FlagParser()
parser.add_flag('--gamma_Iext_e')
parser.add_flag('--gamma_Iext_e_vs_uni_GABA')
args = parser.parse_args()

if args.gamma_Iext_e or args.all:
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
        noisefigs.plotters.GenericGammaPlotter,
        config={
            'GenericGammaPlotter': {
                'what': 'acVal',
                'fname': "gamma_power_Iext_e_theta_vs_Iext_e_const_{ns}.pdf",
                'normalize_ticks': (False, False),  # (Y, X)
                'normalize_type': (None, None),
                'xlabel': '',
                'ylabel': '$\\theta$ amplitude (pA)',
                'xticks': [0, 0, 0],
                'bbox': (.15, .2, .85, .9),

                'cbar_kw': dict(
                    label      = '$1^{st}$ autocorrelation\npeak',
                    location    = 'right',
                    shrink      = 0.8,
                    pad         = .05,
                    ticks      = ti.MultipleLocator(0.1),
                    rasterized  = True
                ),

                'vmin': 0,
                'vmax': 0.43,
            },
        })
    env.register_class(
        noisefigs.plotters.GenericGammaPlotter,
        config={
            'GenericGammaPlotter': {
                'what': 'freq',
                'fname': "gamma_freq_Iext_e_theta_vs_Iext_e_const_{ns}.pdf",
                'normalize_ticks': (False, False),  # (Y, X)
                'normalize_type': (None, None),
                'xlabel': 'Constant amplitude\n(pA)',
                'ylabel': '$\\theta$ amplitude (pA)',
                'bbox': (.15, .2, .85, .9),

                'sigma_title': False,
                'cbar_kw': dict(
                    label      = 'Oscillation\nfrequency (Hz)',
                    location    = 'right',
                    shrink      = 0.8,
                    pad         = .05,
                    ticks      = ti.MultipleLocator(30),
                    rasterized  = True
                ),

                'vmin': 30,
                'vmax': 90,
            },
        })
    env.plot()


if args.gamma_Iext_e_vs_uni_GABA or args.all:
    new_config = deepcopy(config.get_config())
    new_config.update({
        'grids_data_root': None,
        'bump_data_root': ('simulation_data/submission/i_surround/'
                           'Iext_e_const_vs_uni_GABA/gamma_bump'),
        'vel_data_root':  None,
        'const_pos_data_root': None,
        'singleDataRoot': None,

        'noise_sigmas': [0, 150, 300],
        'even_shape': None,
        'output_dir': 'panels',
    })
    env = NoiseEnvironment(user_config=new_config)
    env.register_class(
        noisefigs.plotters.GenericGammaPlotter,
        config={
            'GenericGammaPlotter': {
                'what': 'acVal',
                'fname': "gamma_power_Iext_e_const_vs_uni_GABA_{ns}.pdf",
                'normalize_ticks': (False, False),  # (Y, X)
                'normalize_type': (None, None),
                'xlabel': '',
                'ylabel': 'Constant amplitude\n(pA)',
                'xticks': [0, 0, 0],
                'bbox': (.17, .2, .87, .9),
                'sigma_title': False,

                'cbar_kw': dict(
                    label      = '$1^{st}$ autocorrelation\npeak',
                    location    = 'right',
                    shrink      = 0.8,
                    pad         = .05,
                    ticks      = ti.MultipleLocator(0.1),
                    rasterized  = True
                ),

                'vmin': -0.03,
                'vmax': 0.692,
            },
        })
    env.register_class(
        noisefigs.plotters.GenericGammaPlotter,
        config={
            'GenericGammaPlotter': {
                'what': 'freq',
                'fname': "gamma_freq_Iext_e_const_vs_uni_GABA_{ns}.pdf",
                'normalize_ticks': (False, False),  # (Y, X)
                'normalize_type': (None, None),
                'xlabel': '',
                'ylabel': 'Constant amplitude\n(pA)',
                'bbox': (.17, .2, .87, .9),

                'sigma_title': False,
                'cbar_kw': dict(
                    label      = 'Oscillation\nfrequency (Hz)',
                    location    = 'right',
                    shrink      = 0.8,
                    pad         = .05,
                    ticks      = ti.MultipleLocator(10),
                    rasterized  = True
                ),

                'vmin': 28.57,
                'vmax': 61.03,
            },
        })
    env.plot()
