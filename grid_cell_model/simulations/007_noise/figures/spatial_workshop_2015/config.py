'''Network test configuration file.'''
from __future__ import absolute_import, print_function

from configobj import ConfigObj
import matplotlib.ticker as ti

from noisefigs.plotters.base import PdfOutputSaver


scale_factor = 1.

tick_width = 1. * scale_factor
tick_len   = 6. * scale_factor


def get_config():
    '''Return the configuration object.'''
    _default_config = ConfigObj()
    _default_config.merge({
        'scale_factor': scale_factor,

        'data_root'     : 'simulation_data/network_test/10s_4trials/150pA',
        'data_file_name': 'job00000_output.h5',
        'output_dir'    : 'panels/',

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

        'GridSweepsPlotter': {
            'cbar': [1, 1, 1],
            'sigma_title': False,
            'ylabel': [None, None, None],
            'yticks': [True, True, True],
        },

        'GammaSweepsPlotter': {
            'cbar': [1, 1, 1],
            'cbar_kw' : { # This has to match cbar_kw-s below
                'location': 'right',
            },

            'AC_cbar_kw': dict(
                location   = 'right',
                ticks      = ti.MultipleLocator(0.3),
                fraction   = 0.25,
                shrink     = 0.8,
                pad        = -.05,
                labelpad   = 8,
                label      = '$1^{st}$ autocorrelation\npeak',
                rasterized = True,
            ),
            'AC_xticks': [False]*3,
            'AC_yticks': [1, 1, 1],
            'AC_sigma_title': False,
            'AC_vmin': -0.09,
            'AC_vmax': 0.675,

            'F_cbar_kw': dict(
                location   = 'right',
                ticks      = ti.MultipleLocator(30),
                fraction   = 0.25,
                shrink     = 0.8,
                pad        = -.05,
                labelpad   = 8,
                label      = 'Oscillation\nfrequency (Hz)',
                extend     = 'max',
                extendfrac = 0.1,
                rasterized = True
            ),
            'F_xticks': [True]*3,
            'F_yticks': [1, 1, 1],
            'F_sigma_title': False,
            'F_vmin': 30,
            'F_vmax': 120,

        },

        'GammaExamplePlotter': {
            'xscales': [
                [0, 1, 0],
                [0, 0, 1],
            ],

            'sigma_titles': [
                [0, 0, 0],
                [0, 0, 0],
            ],

            'yscale_kw': [[
                dict(
                    scaleLen=5,
                    unitsText='nA',
                    x=.5, y=.1,
                    size='x-small'
                ),
                dict(
                    scaleLen=0.5,
                    unitsText='nA',
                    x=.5, y=.05,
                    size='x-small'
                ),
                dict(
                    scaleLen=0.5,
                    unitsText='nA',
                    x=.5, y=.05,
                    size='x-small'
                )],

                [dict(
                    scaleLen=5,
                    unitsText='nA',
                    x=.5, y=.1,
                    size='x-small'
                ),
                None,
                dict(
                    scaleLen=0.5,
                    unitsText='nA',
                    x=.55, y=0,
                    size='x-small'
                )]],
        },

        'PopulationActivityPlotter': {
            'raster_rect': (.075, 0.35, 0.99, 0.97),
            'fig_size': (15, 6),
            't_limits': (0, 10e3),

            'snapshot_tstep': 4,
            'e_snapshots_rect': (.075, .15, 0.99, 0.25),
            'i_snapshots_rect': (.075, .02, 0.99, 0.1),
        },

        'Generic2DPBumpPlotter': {
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
            'yticks': [True]*3,
            'plot_grid_contours': [0],
            'ann': [None],
            'bbox': (.15, .17, .9, .9),
        },

        'GenericGammaPlotter': {
            'scale_factor': 1.,
            'cbar': [1],
            'sigma_title': True,
            'cbar_kw': dict(
                label      = '$1^{st}$ autocorrelation\npeak',
                location    = 'right',
                shrink      = 0.8,
                pad         = .05,
                ticks      = ti.MultipleLocator(0.3),
                rasterized  = True
            ),
            'xticks': [True]*3,
            'yticks': [True]*3,
            'plot_grid_contours': [0],
            'ann': [None],
            'bbox': (.15, .17, .9, .9),
            'vmin': None,
            'vmax': None,
        },
    })

    ##########################################################################
    return _default_config
