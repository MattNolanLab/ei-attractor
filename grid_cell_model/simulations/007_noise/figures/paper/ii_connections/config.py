
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

    'GridSweepsPlotter': {
        'scale_factor': .9,
        'cbar': [1, 0, 0],
        'cbar_kw': {
            'label': 'Gridness score',
            'fraction': 0.25,
            'location': 'left',
            'shrink': 0.8,
            'pad': .2,
            'labelpad': 8,
            'ticks': ti.MultipleLocator(0.5),
            'rasterized': True
        },
        'ann': None,
    },

    'MainBumpFormationPlotter': {
    },

    'GammaSweepsPlotter': {

        'F_vmin': 30,
        'F_vmax': 167,
    },

    'GammaExamplePlotter': {
        'yscale_kw': [[
            dict(
                scaleLen=3,
                unitsText='nA',
                x=.5, y=.1,
                size='x-small'
            ),
            dict(
                scaleLen=0.5,
                unitsText='nA',
                x=.5, y=.01,
                size='x-small'
            ),
            dict(
                scaleLen=0.5,
                unitsText='nA',
                x=.5, y=-.1,
                size='x-small'
            )],

            [dict(
                scaleLen=5,
                unitsText='nA',
                x=.5, y=.01,
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
