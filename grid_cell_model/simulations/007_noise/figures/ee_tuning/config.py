
'''Configuration file for the noise paper.'''
from __future__ import absolute_import, print_function

import matplotlib.ticker as ti
from noisefigs.plotters.base import SeparateMultipageSaver

def get_config():
    return _config


_config = {
    'iter_list': ['g_AMPA_total', 'pAMPA_sigma'],

    'grids_data_root':      None,
    'bump_data_root':       'simulation_data/pAMPA_sigma/sigma_sweep',
    'vel_data_root':        None,
    'const_pos_data_root':  None,
    'singleDataRoot':       None,

    'even_shape': (61, 11),
    'noise_sigmas': [150],

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
}

