
'''Configuration file for the noise paper.'''
from __future__ import absolute_import, print_function

import os.path

import matplotlib.ticker as ti
from noisefigs.plotters.base import SeparateMultipageSaver

def get_config():
    return _config


ROOT_DIR = ['simulation_data', 'submission', 'ee_connections']

_config = {
    'grids_data_root': None,
    'bump_data_root': os.path.join(*(ROOT_DIR + ['gamma_bump'])),
    'vel_data_root':  os.path.join(*(ROOT_DIR + ['velocity'])),
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

        'AC_vmin': -0.12,
        'AC_vmax': 0.68,

        'F_vmin': 30,
        'F_vmax': 121.5,

        'F_cbar_kw': dict(
            extend = 'neither',
        ),
    },

    'GammaExamplePlotter': {
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
                scaleLen=3,
                unitsText='nA',
                x=.5, y=.05,
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

    'MaxPopulationFRSweepsPlotter': {
        'plot_grid_contours': [0, 0, 0],
    },

    'PSeizureSweepPlotter': {
        'plot_grid_contours': [0, 0, 0],
        'FRThreshold': 300,
    },

    'BumpDriftAtTimePlotter': {
        'plot_grid_contours': [0, 0, 0],
    },

    'VelFitErrSweepPlotter': {
        'plot_contours': [0, 0, 0],
        'vmin': 0.2,
        'vmax': 13.93,
    },

    'VelFitStdSweepPlotter': {
        'plot_contours': [0, 0, 0],
    },

    'VelSlopeSweepPlotter': {
        'plot_contours': [0, 0, 0],
        'vmin': -0.26,
        'vmax': 1.26,
    },
}
