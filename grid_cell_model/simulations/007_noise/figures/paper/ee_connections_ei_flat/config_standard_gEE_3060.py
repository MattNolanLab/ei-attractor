
'''Configuration file for the noise paper.'''
from __future__ import absolute_import, print_function

import os.path

from configobj import ConfigObj
import matplotlib.ticker as ti
from noisefigs.plotters.base import SeparateMultipageSaver


DATA_ROOT = ['simulation_data', 'submission', 'ee_connections_ei_flat',
             'standard_sweep_g_EE_3060_pEE_sigma_0_0833']

scale_factor = 1.

tick_width = 1. * scale_factor
tick_len   = 6. * scale_factor


def get_config():
    return _config


_config = ConfigObj({
    'grids_data_root': None,
    'bump_data_root': os.path.join(*(DATA_ROOT + ['gamma_bump'])),
    'vel_data_root':  None,
    'const_pos_data_root': None,
    'singleDataRoot': None,
    'connection_data_root': os.path.join(*(DATA_ROOT + ['connections'])),

    'output_dir': 'panels_standard_gEE_3060',

    'scale_factor': scale_factor,

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

    'GridExampleRectPlotter': {
        'fig_saver': SeparateMultipageSaver(None, 'pdf')
    },


    'MainBumpFormationPlotter': {
        'plot_grid_contours': [0, 0, 0],
        'ann': [None, None, None],
    },

    'GammaSweepsPlotter': {

        'AC_vmin': -0.15,
        'AC_vmax': 0.6,

        'F_vmin': 30,
        'F_vmax': 142.9,

        'F_cbar_kw': dict(
            extend     = 'neither',
        ),
    },

    'MaxPopulationFRSweepsPlotter': {
    },

    'PSeizureSweepPlotter': {
        'FRThreshold': 300,
    },

    'BumpDriftAtTimePlotter': {
        'plot_grid_contours': [0, 0, 0],
    },

    'VelFitErrSweepPlotter': {
    },

    'VelFitStdSweepPlotter': {
    },

    'VelSlopeSweepPlotter': {
        'vmin': -.24,
        'vmax': 1.3306,
    },

    'GridSweepsPlotter': {
        'vmin': -0.5,
        'vmax': 1.089,
    },

    'IsBumpExamplePlotter': {
    }
})

