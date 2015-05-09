
'''Configuration file for the noise paper.'''
from __future__ import absolute_import, print_function

import os.path

import matplotlib.ticker as ti
from noisefigs.plotters.base import SeparateMultipageSaver

def get_config():
    return _config


ROOT_DIR = ['simulation_data', 'submission', 'ii_connections', 'gE_vs_gI']

_config = {
    'grids_data_root': os.path.join(*(ROOT_DIR + ['grids'])),
    'bump_data_root': os.path.join(*(ROOT_DIR + ['gamma_bump'])),
    'vel_data_root':  os.path.join(*(ROOT_DIR + ['velocity'])),
    'const_pos_data_root': None,
    'singleDataRoot': None,

    'GridExampleRectPlotter': {
        'fig_saver': SeparateMultipageSaver(None, 'pdf')
    },


    'MainBumpFormationPlotter': {
    },

    'GammaSweepsPlotter': {

        'F_vmin': 30,
        'F_vmax': 167,
    },

    'MaxPopulationFRSweepsPlotter': {
    },

    'PSeizureSweepPlotter': {
        'FRThreshold': 300,
    },

    'BumpDriftAtTimePlotter': {
    },

    'VelFitErrSweepPlotter': {
        'vmin': 0.1,
        'vmax': 12.101,
    },

    'VelFitStdSweepPlotter': {
    },

    'VelSlopeSweepPlotter': {
        'vmin': -.6,
        'vmax': 1.64,
    },
}
