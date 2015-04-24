
'''Configuration file for the noise paper.'''
from __future__ import absolute_import, print_function

import os.path

import matplotlib.ticker as ti
from noisefigs.plotters.base import SeparateMultipageSaver

def get_config():
    return _config


_config = {
    'grids_data_root': None,
    'bump_data_root': os.path.join('simulation_data', 'submission',
                                    'probabilistic_connections',
                                    'gamma_bump'),
    'vel_data_root':  None,
    'const_pos_data_root': None,
    'singleDataRoot': None,

    'GridExampleRectPlotter': {
        'fig_saver': SeparateMultipageSaver(None, 'pdf')
    },


    'MainBumpFormationPlotter': {
        'plot_grid_contours': [0, 0, 0],
    },

    'GammaSweepsPlotter': {
        'plot_grid_contours': [0, 0, 0],

        'F_vmin': 30,
        'F_vmax': 167,
    },

    'MaxPopulationFRSweepsPlotter': {
        'plot_grid_contours': [0, 0, 0],
    },

    'PSeizureSweepPlotter': {
        'FRThreshold': 300,
        'plot_grid_contours': [0, 0, 0],
    },

    'BumpDriftAtTimePlotter': {
        'plot_grid_contours': [0, 0, 0],
    },
}

