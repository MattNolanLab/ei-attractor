#!/usr/bin/env python
'''Ratio corrected gamma figures.'''
from __future__ import absolute_import, print_function

import copy

import matplotlib
import matplotlib.ticker as ti
from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

import config_ratio_corrected as config

parser = flagparse.FlagParser()
parser.add_flag('--gamma_sweep')
parser.add_flag('--examples')
args = parser.parse_args()


env = NoiseEnvironment(user_config=config.get_config())

if args.gamma_sweep or args.all:
    env.register_class(
        noisefigs.plotters.GenericGammaPlotter,
        config={
            'GenericGammaPlotter': {
                'scale_factor': .7,
                'what': 'acVal',
                'fname': "gamma_sweeps{ns}.pdf",
                'normalize_ticks': (True, True),  # (Y, X)
                'normalize_type': ('E', 'I'),
                'xlabel': '',
                'xticks': [False, False, False],
                'ylabel': None,
                'yticks': [True, False, False],
                'bbox': (.2, .17, .85, .87),
                'vmin': .0,
                'vmax': .505,
                'cbar': [0, 0, 1],
                'cbar_kw': dict(
                    ticks = ti.MultipleLocator(0.1),
                ),
                'plot_grid_contours': [0, 0, 0],

                'ann': [
                    dict(
                        txt='b',
                        rc=(5, 15),
                        xytext_offset=(.1, 0),
                        color='white',
                    ),
                    dict(
                        txt='a',
                        rc=(15, 5),
                        xytext_offset=(-.05, 20.),
                        color='white',
                    ),
                ],
            },
        })
    env.register_class(
        noisefigs.plotters.GenericGammaPlotter,
        config={
            'GenericGammaPlotter': {
                'scale_factor': .7,
                'what': 'freq',
                'fname': "gamma_freq_sweeps{ns}.pdf",
                'normalize_ticks': (True, True),  # (Y, X)
                'normalize_type': ('E', 'I'),
                'xlabel': None,
                'ylabel': None,
                'yticks': [True, False, False],
                'bbox': (.2, .17, .85, .87),
                'vmin': 24,
                'vmax': 145.2,
                'cbar': [0, 0, 1],
                'cbar_kw': dict(
                    label='Frequency (Hz)',
                    ticks=ti.MultipleLocator(20),
                ),
                'plot_grid_contours': [0, 0, 0],
                'ann': [
                    dict(
                        txt='b',
                        rc=(5, 15),
                        xytext_offset=(.1, 0),
                        color='white',
                    ),
                    dict(
                        txt='a',
                        rc=(15, 5),
                        xytext_offset=(-.05, 20.),
                        color='white',
                    ),
                ],
                'sigma_title': False,
            },
        })

if args.examples or args.all:
    env.register_plotter(noisefigs.plotters.GammaExamplePlotter)

env.plot()
