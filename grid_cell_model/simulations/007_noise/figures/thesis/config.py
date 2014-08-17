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

    'BumpSigmaSweepPlotter': {
        'cbar': [1, 0, 0],
        'cbar_kw': {
            'location': 'left',
            'pad': .25,
            'labelpad': 3,
        },
        'scale_factor': .75,
    },

    'BumpExamplePlotter': {
        'scale_factor': .75,
        'bbox': (0.01, 0.01, 0.99, 0.8),
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

    'BumpDiffAtInitPlotter': {
        'scale_factor': .75,
        'cbar_kw': {
            'label': 'Initialisation error\n(neurons)',
            'pad': -.15,
        },
    },

    'BumpDriftAtTimePlotter': {
        'scale_factor': .75,
        'cbar_kw': {
            'pad': -.15,
        },
    },

    'BumpDiffResetPlotter': {
        'scale_factor': .75,
        'cbar_kw': {
            'pad': -.15,
        },
    },

    'PSeizureGridsScatterAllPlotter': {
        'fig_size': (2.5, 3),         # inches
        'scale_factor': .75,
        'bbox_rect': (0.35, 0.22, 0.92, 0.9),
        'legend_kwargs': dict(
            loc=(0.3, 0.6),
        ),
    },

    'MaxFRGridsScatterAllPlotter': {
        'fig_size': (2.5, 3),         # inches
        'scale_factor': .75,
        'bbox_rect': (0.35, 0.22, 0.92, 0.9),
    },

    'PSeizureSweepPlotter': {
        'fig_size': (2.6, 3.7),  # inches
        'cbar_kw': {
            'location': 'top',
            'pad': .1,
            'labelpad': 3,
        },
        'scale_factor': .75,
    },

    'RasterExamplePlotter': {
        'fig_size': (5.5, 8.3),
        'sweep_rect' : (.1, .73, .5, .95),
        'cbar_kw': dict(
            pad         = .1,
        ),
        'plot_ann_txt' : False,
    },
##############################################################################
}
