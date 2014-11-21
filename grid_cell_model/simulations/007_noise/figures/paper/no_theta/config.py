
'''Configuration file for the noise paper.'''
from __future__ import absolute_import, print_function

import matplotlib.ticker as ti

def get_config():
    return _config


_config = {
    #"fname_prefix": "no_theta_",

    'grids_data_root':      'simulation_data/submission/no_theta/grids',
    'bump_data_root':       'simulation_data/submission/no_theta/gamma_bump',
    'vel_data_root':        'simulation_data/submission/no_theta/velocity',
    'const_pos_data_root':  'simulation_data/submission/no_theta/const_position',
    'singleDataRoot':       'simulation_data/submission/no_theta/single_neuron',

    'GammaSweepsPlotter': {
        'plot_grid_contours': [0, 0, 0],
        'AC_vmin': -0.21,
        'AC_vmax': 0.97,
        'F_vmin': 25,
        'F_vmax': 142,
        'F_cbar_kw': dict(
            extend     = 'neither',
        ),
    },

    'GammaExamplePlotter' : {
        'yscale_kw': [[
            dict(
                scaleLen=1,
                unitsText='nA',
                x=.7, y=-.1,
                size='x-small'
            ),
            dict(
                scaleLen=0.25,
                unitsText='nA',
                x=.5, y=-.2,
                size='x-small'
            ),
            dict(
                scaleLen=0.5,
                unitsText='nA',
                x=.5, y=-.1,
                size='x-small'
            )],

            [dict(
                scaleLen=1,
                unitsText='nA',
                x=.6, y=-.1,
                size='x-small'
            ),
            dict(
                scaleLen=0.5,
                unitsText='nA',
                x=.6, y=-.2,
                size='x-small'
            ),
            dict(
                scaleLen=0.5,
                unitsText='nA',
                x=.55, y=-.25,
                size='x-small'
            )]],
    },

    'MainBumpFormationPlotter': {
        'plot_grid_contours': [0, 0, 0],
    },

    'BumpDriftAtTimePlotter': {
        'plot_grid_contours': [0, 0, 0],
    },

    'MaxPopulationFRSweepsPlotter': {
        'plot_grid_contours': [0, 0, 0],
    },

    'PSeizureSweepPlotter': {
        'plot_grid_contours': [0, 0, 0],
    },
}

