
'''Configuration file for the noise paper.'''
from __future__ import absolute_import, print_function

import os.path

import matplotlib.ticker as ti
from noisefigs.plotters.base import SeparateMultipageSaver


DATA_ROOT = ['simulation_data', 'submission', 'probabilistic_connections']


def get_config():
    return _config


_config = {
    'grids_data_root': os.path.join(*(DATA_ROOT + ['grids'])),
    'bump_data_root': os.path.join(*(DATA_ROOT + ['gamma_bump'])),
    'vel_data_root':  os.path.join(*(DATA_ROOT + ['velocity'])),
    'const_pos_data_root': None,
    'singleDataRoot': None,
    'connection_data_root': os.path.join(*(DATA_ROOT + ['connections'])),

    'GridExampleRectPlotter': {
        'fig_saver': SeparateMultipageSaver(None, 'pdf')
    },


    'MainBumpFormationPlotter': {
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
}

