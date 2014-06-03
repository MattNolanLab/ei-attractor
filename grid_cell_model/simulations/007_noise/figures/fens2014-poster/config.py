'''Configuration file for the FENS2014 poster figures.'''
from __future__ import absolute_import, print_function

scale_factor = 2.5

tick_width = 1. * scale_factor
tick_len   = 6. * scale_factor

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
        'xtick.minor.size'  : tick_len / 2.,
        'xtick.minor.width' : tick_width,
        'xtick.direction'   : 'out',

        'ytick.major.size'  : tick_len,
        'ytick.major.width' : tick_width,
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
    },

    'GridExamplesPlotter': {
        'scale_factor': .8,
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
    },

    'EIRatePlotter': {
        'fig_size': (3, .5),
        'rateTop': .85
    },

    'MainBumpFormationPlotter': {
        'xticks' : [True]*3,
    },

    'GammaSweepsPlotter': {
        'AC_xticks': [True]*3,
    }
}


def get_config():
    return _config

##############################################################################

