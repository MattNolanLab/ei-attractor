'''Network test configuration file.'''
from __future__ import absolute_import, print_function

import os.path
from configobj import ConfigObj
import matplotlib.ticker as ti


scale_factor = 1.

tick_width = 1. * scale_factor
tick_len   = 6. * scale_factor

DATA_ROOT = ['simulation_data', 'i_place_cells']


def get_config():
    '''Return the configuration object.'''
    _default_config = ConfigObj()
    _default_config.merge({
        #'grids_data_root': os.path.join(*(DATA_ROOT + ['10_trials_rate_200_field_std_40'])),
        'grids_data_root': os.path.join(*(DATA_ROOT + ['10_trials_rate_100_field_std_80'])),
        'bump_data_root': None,
        'vel_data_root':  None,
        'const_pos_data_root': None,
        'singleDataRoot': None,
        'connection_data_root': None,

        'scale_factor': scale_factor,

        'output_dir'    : 'panels_weight_sparsity_trials_weight_0/',
        'noise_sigmas': [150],
        'even_shape': None,

        # Sections
        'mpl': {
            'font.size': 11,
            'pdf.fonttype': 42,
            'mathtext.default': 'regular',
            'font.sans-serif': ['Helvetica', 'Avant Garde',
                                'Computer Modern Sans serif'],

            'xtick.major.size'  : tick_len,
            'xtick.major.width' : tick_width,
            'xtick.minor.size'  : tick_len / 2.,
            'xtick.minor.width' : tick_width,
            'xtick.direction'   : 'out',

            'ytick.major.size'  : tick_len,
            'ytick.major.width' : tick_width,
            'ytick.minor.size'  : tick_len / 2.,
            'ytick.minor.width' : tick_width,
            'ytick.direction'   : 'out',
        },

        'IPCGridSweepsPlotter': {
            'scale_factor': 1.,
            'cbar': [1, 1, 1],
            'cbar_kw': dict(
                label       = "Gridness score",
                location    = 'right',
                shrink      = 0.8,
                pad         = .05,
                ticks       = ti.MultipleLocator(0.2),
                rasterized  = True
            ),
            'xlabel': 'Random seed',
            'ylabel': 'Weight (nS)',
            'xticks': [True]*3,
            'yticks': [True, False, False],
            'ann': [None, None, None],
            'bbox': (.15, .17, .9, .9),

            'normalize_ticks': [False, False],
            'vmin': -0.216,
            'vmax': 1.006,
        },

        'IPCScatterPlotter': {
            'scale_factor': 1.,
            'fig_size': (10, 5),

            'xlabel': 'Random seed',
            'ylabel': 'Weight (nS)',

            'normalize_ticks': [False, False],
        },

        'IPCHistogramsPlotter': {
            'scale_factor': 1.,
            'fig_size': (4.2, 2.7),
            'bbox': (.17, .2, .55, .85),

            'normalize_ticks': [False, False],
        },

        'IPCExamplePlotter': {
            'fig_size': (1, 1),
            'bbox': (.02, .02, .85, .85),

            'normalize_ticks': [False, False],
            'population_type': 'E',
        },

        'IPCExampleColorbarPlotter': {
            'fig_size': (0.6, 0.75),
            'bbox': (.02, .02, .90, .90),
        },
    })

    ##########################################################################
    return _default_config
