'''Configuration file for the FENS2014 poster figures.'''
from __future__ import absolute_import, print_function

def get_config():
    return _config


_config = {
    'GridSweepsPlotter': {
        'scale_factor': .75,
        'cbar': [1, 0, 0],
        'cbar_kw': {
            'location': 'left',
            'pad': .2,
            'labelpad': 0,
            'rasterized': False,
        },
    },

    'GridExamplesPlotter': {
        'scale_factor': .75,
    },

    'GammaSweepsPlotter': {
        'scale_factor': .75,
        'AC_cbar_kw': {
            'labelpad': 2,
        },
        #'AC_xticks': [True]*3,
        'F_cbar_kw': {
            'labelpad': 2,
        },
        #'F_sigma_title': True,
    },

    'GammaExamplePlotter': {
        'scale_factor': .75,
    },

    'GammaScatterAllPlotter': {
        'fig_size': (5.25, 3.2),
        'tight_layout_kwargs': {
            'pad': 0,
            'rect': (.01, 0, 0.95, .99),
        },
    },

    'MainBumpFormationPlotter': {
        'scale_factor': .75,
    },

    'IsBumpExamplePlotter': {
        'scale_factor': .75,
        'bumpQualityX': -.8,
    },

    'MainScatterGridsBumpsPlotter': {
        'fig_size': (5.25, 3.2),
        'tight_layout_kwargs': {
            'pad': 0,
            'rect': (.01, 0.02, .99, 0.85),
        },
        'legend_kwargs': {
            'frameon': False,
        },
    },

    'EIRasterPlotter': {
        'fig_size': (2.9, 2),
        'scale_factor': .75,
        'ylabelPos': -.37,
    },

    'EIRatePlotter': {
        'fig_size': (2.9, .65),
        'scale_factor': .75,
        'rateTop': .85,
        'ylabelPos': -.37,
    },

    'MaxPopulationFRSweepsPlotter': {
        'scale_factor': .75,
        'cbar_kw': {
            'labelpad': 2,
        },
    },
}
