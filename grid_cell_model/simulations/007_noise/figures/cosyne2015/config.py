
'''Configuration file for the noise paper.'''
from __future__ import absolute_import, print_function

import matplotlib.ticker as ti

def get_config():
    return _config


_config = {
    'GridSweepsPlotter': {
        'scale_factor': 0.75,
        'cbar': [0, 1, 0],
        'cbar_kw': {
            'fraction': .25,
        },
        'ann': None,
    },

    'ConnectionFunctionPlotter': {
        'fig_size': (2, 1.5),
        'x_range': (0, .5, .001),
        'xticks': (0, .5),
        'bbox_rect': (.31, .25, .93, .75),

        'leg1_kwargs': dict(
            loc=(-.15, 1.0),
        ),
        'leg2_kwargs': dict(
            loc=(.35, 1.03),
        ),
    },

    'VmExamplesPlotter': {
        'fig_size': (1.7, .7),
    },

    'GammaSweepsPlotter': {
        'scale_factor': .75,
        'cbar': [1, 1, 1],
        'ann': None,

        'AC_yticks': [1, 1, 1],
        'F_yticks': [1, 1, 1],

        'cbar_kw' : { # This has to match cbar_kw-s below
            'location': 'right',
        },
        'AC_cbar_kw': dict(
            location   = 'right',
            pad        = -.05,
        ),
        'F_cbar_kw': dict(
            location   = 'right',
            pad        = -.05,
        ),

        'ann': None,
        'annF': None,
    },
}

