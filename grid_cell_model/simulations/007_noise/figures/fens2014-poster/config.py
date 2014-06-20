'''Configuration file for the FENS2014 poster figures.'''
from __future__ import absolute_import, print_function

scale_factor = 2.5

tick_width = 1. * scale_factor
tick_len   = 6. * scale_factor

def get_config():
    return _config


_config = {
    'scale_factor': scale_factor,

    # Sections
    'mpl': {
        'font.size'        : 30,
        'lines.linewidth'  : scale_factor,
        'lines.markersize' : 6. * scale_factor,
        'axes.linewidth'   : scale_factor,

        'xtick.major.size'  : tick_len,
        'xtick.major.width' : tick_width,
        'xtick.major.pad'   : 4*scale_factor,
        'xtick.minor.size'  : tick_len / 2.,
        'xtick.minor.width' : tick_width,
        'xtick.direction'   : 'out',

        'ytick.major.size'  : tick_len,
        'ytick.major.width' : tick_width,
        'ytick.major.pad'   : 4*scale_factor,
        'ytick.minor.size'  : tick_len / 2.,
        'ytick.minor.width' : tick_width,
        'ytick.direction'   : 'out',

    },


    'GridSweepsPlotter': {
        'scale_factor' : .8,
        'cbar': [1, 0, 0],
        'cbar_kw' : {
            'label': '',
            'location': 'left',
            'pad': .2,
        },
        'sigma_title': False,

        'ann': [
            dict(
                txt='a',
                rc=(5, 15),
                xytext_offset=(0.5, 1.5),
                color='black'
            )
        ],
    },

    'GridExamplesPlotter': {
        'scale_factor': .8,
    },

    'ConnectionFunctionPlotter': {
        'fig_size': (3, 2),
        'ylabel_coords': (-.1, .5),
    },

    'VmExamplesPlotter': {
        'fig_size': (2.3, 1.25),
        'scale_factor': .9,
    },

    'GammaScatterAllPlotter': {
        'fig_size': (3.2, 3.2),
        'legend_kwargs': dict(
            loc=(0, 1.),
            fontsize='small',
            frameon=False,
            scatterpoints=1,
            ncol=3,
        ),
        'tight_layout_kwargs': {
            'pad': 3.,
            'rect': (.01, .01, .99, .85),
            'pad': 0,
        },
    },

    'MainScatterGridsBumpsPlotter': {
        'fig_size': (4., 2.5),
        'tight_layout_kwargs': {
            'rect': (0.05, 0.05, 0.95, 0.85),
        },
        'legend_kwargs': dict(
            loc=(0.05, 1.02),
            frameon=False,
        ),
    },

    'EIRasterPlotter': {
        'fig_size': (3, 1.5),
        'fig_ext': 'pdf',
    },

    'EIRatePlotter': {
        'fig_size': (3, .5),
        'rateTop': .85
    },

    'MainBumpFormationPlotter': {
        'xticks' : [True]*3,
        'ann': ([], [], []),
    },

    'GammaSweepsPlotter': {
        'AC_xticks': [True]*3,
        'ann': [
            dict(
                txt='a',
                rc=(5, 15),
                xytext_offset=(1.5, 1),
                color='white',
            ),
        ],
    },

    'GammaExamplePlotter': {
        'scale_factor': .9,
        'xscales': [
            [0, 0, 1],
            [0, 0, 0],
        ],
        'sigma_titles': [
            [1, 1, 1],
            [0, 0, 0],
        ],

        'xscale_kw': dict(
            x=0.75, y=.2,
        ),
    },
}


def get_config():
    return _config

##############################################################################

