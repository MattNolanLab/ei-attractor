'''Configuration file for the FENS2014 poster figures.'''
from __future__ import absolute_import, print_function

from noisefigs.plotters.base import SeparateMultipageSaver

def get_config():
    return _config


_config = {
    'GridSweepsPlotter': {
        'cbar': [1, 1, 1],
        'scale_factor': .75,
    },

    'GridExamplesPlotter': {
        'scale_factor': .75,
    },

    'GridDetailedNoisePlotter': {
        'legend':  ['$(g_E,\ gI)$ = (3, 1) nS',  '$(g_E,\ gI)$ = (1, 3) nS'],
        'legend_kwargs': dict(
            loc=(0.6, 1),
        ),
    },

    'GammaSweepsPlotter': {
        'fig_size': (2.6, 4.5),  # inches
        'scale_factor': .67,
        'cbar': [1, 0, 1],
        'AC_cbar_kw': {
            'location': 'bottom',
            'pad': .2,
            'labelpad': 2,
        },
        'AC_xticks': [1, 0, 1],
        'AC_yticks': [1, 1, 1],
        'AC_sigma_title': False,

        'F_cbar_kw': {
            'location': 'bottom',
            'pad': .2,
            'labelpad': 2,
        },
        'F_xticks': [1, 0, 1],
        'F_yticks': [0, 0, 0],
        'F_sigma_title': False,
    },

    'GammaExamplePlotter': {
        'scale_factor': .75,
        'sigma_titles': [
            [0, 0, 0],
            [0, 0, 0],
        ],
        'yscale_kw': [[
            dict(
                scaleLen=5,
                unitsText='nA',
                x=.3, y=0,
                size='x-small'
            ),
            dict(
                scaleLen=0.5,
                unitsText='nA',
                x=.55, y=0,
                size='x-small'
            ),
            dict(
                scaleLen=0.5,
                unitsText='nA',
                x=.55, y=-.1,
                size='x-small'
            )],

            [dict(
                scaleLen=5,
                unitsText='nA',
                x=.3, y=0,
                size='x-small'
            ),
            dict(
                scaleLen=0.5,
                unitsText='nA',
                x=.55, y=0,
                size='x-small'
            ),
            dict(
                scaleLen=0.5,
                unitsText='nA',
                x=.55, y=-.1,
                size='x-small'
            )]],
        },

    'GammaDetailedNoisePlotter': {
        'legend':  ['$(g_E,\ gI)$ = (3, 1) nS',  '$(g_E,\ gI)$ = (1, 3) nS'],
        'legend_kwargs': dict(
            loc=(.6, .75),
        ),
    },

    'GridExampleRectPlotter': {
        'fig_saver': SeparateMultipageSaver(None, 'pdf')
    },

    'GammaScatterAllPlotter': {
        'fig_size': (3, 2.75),
        'bbox_rect': (.25, .35, .95, .85),
        'ylabel': 'Gridness score',
    },

    'GammaFreqGridsScatterAllPlotter': {
        'fig_size': (3, 2.75),
        'bbox_rect': (.25, .35, .95, .85),
        'legend_kwargs': dict(
            loc=(0.75, 0.5),
        ),
        'ylabel': '',
        'yticks': False,
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
        'cbar': [0, 0, 1],
        'cbar_kw': {
            'location': 'right',
            'pad': -.05,
            'labelpad': 5,
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
        'ylabelPos': -0.12,
        'plot_ann_txt' : False,
    },

    'GammaGridsProbabilityPlotter': {
        'fig_size': (.65, .65),  # inches
        'scale_factor': .8,
        'bbox_rect': (0.05, 0.05, 0.95, 0.7),
        'strip_axis': True,
        'title_size': 'x-small',
    },

    'GammaFreqGridsProbabilityPlotter': {
        'fig_size': (.65, .65),  # inches
        'scale_factor': .8,
        'bbox_rect': (0.05, 0.05, 0.95, 0.7),
        'strip_axis': True,
        'title_size': 'x-small',
    },

    'ScatterGammaGridsSeparatePlotter': {
        'fig_size': (5.6, 7.9),
    },

    'GridBumpScatterPlotter': {
        'fig_size': (5.6, 7.9),
        'color_box_coords': {
            'left': .2,  # w = 0.165
            'bottom': .87,
            'right': .2 + .165,
            'top': .97
        }
    },

    'FracTotalSweepAnnPlotter': {
        'cbar': (0, 0, 1),
        'cbar_kw': dict(
            location    = 'right',
            pad         = -.05,
        )
    },

##############################################################################
}
