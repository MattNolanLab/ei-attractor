
'''Configuration file for the noise paper.'''
from __future__ import absolute_import, print_function

import matplotlib.ticker as ti
from configobj import ConfigObj

scale_factor = 1.

tick_width = 1. * scale_factor
tick_len   = 6. * scale_factor


def get_config():
    '''Return the configuration.'''
    return _config


_config = ConfigObj({
    'iter_list': ['g_AMPA_total', 'pAMPA_sigma'],

    'grids_data_root':      None,
    'bump_data_root':       'simulation_data/other/ee_fine_tuning/pAMPA_sigma/sigma_sweep',
    'vel_data_root':        None,
    'const_pos_data_root':  None,
    'singleDataRoot':       None,

    'even_shape': (61, 11),
    'noise_sigmas': [150],

    'scale_factor' : scale_factor,

    'mpl': {
        'font.size': 11,
        'pdf.fonttype': 42,
        'mathtext.default': 'regular',
        'font.sans-serif'    : ['Helvetica', 'Avant Garde', 'Computer Modern Sans serif'],

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

    'GEProfileWidthBumpPlotter': {
        'scale_factor': 1.,
        'cbar': [1],
        'cbar_kw': dict(
            label       = "P(bumps)",
            location    = 'right',
            shrink      = 0.8,
            pad         = .05,
            ticks       = ti.MultipleLocator(0.5),
            rasterized  = True
        ),
        'xticks': [True]*3,
        'plot_grid_contours': [0],
        'ann': [None],
        'xlabel': "$\sigma_{E{\\rightarrow}I}$",
        'bbox': (.15, .17, .9, .9),
    },
})

