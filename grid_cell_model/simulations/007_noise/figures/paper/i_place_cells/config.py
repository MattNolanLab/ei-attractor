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
        'grids_data_root': os.path.join(*(DATA_ROOT + ['grids_max_rate_100_field_std_80'])),
        'bump_data_root': None,
        'vel_data_root':  None,
        'const_pos_data_root': None,
        'singleDataRoot': None,
        'connection_data_root': None,

        'scale_factor': scale_factor,

        'output_dir'    : 'panels_weight_sparsity/',
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
            'xlabel': 'Weight (nS)',
            'ylabel': '# PCs connected',
            'xticks': [True]*3,
            'yticks': [True, False, False],
            'ann': [None, None, None],
            'bbox': (.15, .17, .9, .9),

            'normalize_ticks': [False, False],
            'vmin': None,
            'vmax': None,
        },
    })

    ##########################################################################
    return _default_config
