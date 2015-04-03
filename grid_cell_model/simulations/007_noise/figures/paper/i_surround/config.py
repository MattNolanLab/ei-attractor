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

        'grids_data_root':      'None',
        'bump_data_root':       'simulation_data/submission/i_surround/gamma_bump',
        'vel_data_root':        None,
        'const_pos_data_root':  None,
        'singleDataRoot':       None,

        'output_dir'    : 'panels/',
        'noise_sigmas'  : [0, 150, 300],

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

        'MainBumpFormationPlotter': {
            'plot_grid_contours': [0, 0, 0],
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

        'GammaSweepsPlotter': {
            'plot_grid_contours': [0, 0, 0],
            'AC_vmin': .08,
            'AC_vmax': .65,
            'AC_cbar_kw': {
                'ticks': ti.MultipleLocator(0.1),
            },
            'F_vmin': 26,
            'F_vmax': 76.4,
            'F_cbar_kw': {
                'ticks': ti.MultipleLocator(10),
                'extend': 'neither',
            },
            },

        'GammaExamplePlotter': {
            'mon_idx_e': 1,
            'mon_idx_i': 1,
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
                dict(
                    scaleLen=0.5,
                    unitsText='nA',
                    x=.5, y=.05,
                    size='x-small'
                ),
                dict(
                    scaleLen=0.5,
                    unitsText='nA',
                    x=.55, y=0,
                    size='x-small'
                )]],
        },

        'BumpDriftAtTimePlotter': {
            'plot_grid_contours': [0, 0, 0],
        },
    })

    ##########################################################################
    return _default_config
