'''Configuration file for the FENS2014 poster figures.'''
from __future__ import absolute_import, print_function

import matplotlib.ticker as ti

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

    'sweeps': {
        'contours_kwargs': {
            'hold': True,
            'colors': 'k',
            'linewidths': [1.5 * scale_factor]
        },
    },

    'GridSweepsPlotter': {
        'scale_factor' : .7,
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
        'scale_factor': .7,
    },

    'ConnectionFunctionPlotter': {
        'fig_size': (2, 1.5),
        'x_range': (0, .5, .001),
        'xlabel': 'Distance',
        'xticks': (0, .5),
        'bbox_rect': (.31, .25, .93, .75),
        'uniform_random': False,

        'leg1_kwargs': dict(
            loc=(.2, .9),
        ),
        'leg2_kwargs': dict(
            loc=(.35, 1.03),
        ),
    },

    'VmExamplesPlotter': {
        'fig_size': (2.3, 1.25),
        'scale_factor': .9,
    },

    'GammaScatterAllPlotter': {
        'fig_size': (3, 2.3),
        'dot_size': 10 * scale_factor,
        'legend_kwargs': dict(
            loc=(0, 1.),
            fontsize='small',
            frameon=False,
            scatterpoints=1,
            ncol=3,
        ),
        'bbox_rect': (.225, .25, .95, .85),
        'ylabel': 'Gridness score',
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
        'scale_factor': .9,
        'fig_size': (3, 1.5),
        'fig_ext': 'pdf',
        'scaleBar': [None, 25, 25],
        'scaleX': .85,
        'scaleY': -.05,

        'plot_theta' : True,
        'theta_rect' : [.28, .77, .95, .83],
        'theta_color': (0, 0, 0, .3),
    },

    'EIRatePlotter': {
        'scale_factor': .9,
        'fig_size': (3, .5),
        'rateTop': .8
    },

    'MainBumpFormationPlotter': {
        'scale_factor': .7,
        'xticks' : [True]*3,
        'ann': [
            [dict(txt='a',
                  rc=(5, 15),
                  xytext_offset=(1.5, 1),
                  color='white')],
            [dict(txt='a',
                  rc=(5, 15),
                  xytext_offset=(1.5, 1),
                  color='white')],
            [dict(txt='a',
                  rc=(5, 15),
                  xytext_offset=(1.5, 1),
                  color='white')],
        ],
        'cbar': [0, 1, 1],
    },

    'GammaSweepsPlotter': {
        'scale_factor': .8,
        'ann': [
            dict(
                txt='a',
                rc=(5, 15),
                xytext_offset=(1.5, 1),
                color='white',
            ),
        ],
        'annF': [
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
            [0, 0, 0],
            [0, 0, 0],
        ],

        'xscale_kw': dict(
            x=0.75, y=.2,
        ),
    },

    'PSeizureSweepPlotter': {
        'scale_factor': .7,
        'cbar_kw': dict(
            location    = 'right',
            shrink      = 0.8,
            pad         = -.05,
            ticks       = ti.MultipleLocator(0.5),
            rasterized  = True
        ),

        'ann': [
            dict(
                txt='a',
                rc=(5, 15),
                xytext_offset=(1.5, 1),
                color='white',
            ),
        ],
    },
}


def get_config():
    return _config

##############################################################################

