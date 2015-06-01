#!/usr/bin/env python
'''Ratio corrected bump figures.'''
from __future__ import absolute_import, print_function

from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

import config_ratio_corrected as config

parser = flagparse.FlagParser()
parser.add_flag('--pbumps_sweep')
parser.add_flag('--pbumps_threshold_sweep')
parser.add_flag('--bump_examples')
parser.add_flag('--bump_examples_isbump')
args = parser.parse_args()


env = NoiseEnvironment(user_config=config.get_config())

if args.pbumps_sweep or args.all:
    env.register_plotter(
        noisefigs.plotters.Generic2DPBumpPlotter,
        config={
            'Generic2DPBumpPlotter': {
                'scale_factor': 0.8,
                'cbar': [0, 0, 1],
                'fname': "bumps_mainFig_isBumpFracTotal_sweeps_annotated{ns}.pdf",
                'normalize_ticks': (True, True),  # (Y, X)
                'normalize_type': ('E', 'I'),
                'xlabel': None,
                'yticks': [True, False, False],
                'ylabel': None,
                'bbox': (.15, .2, .85, .9),
                'ann': [
                    [dict(txt='a',
                          rc=(15, 5),
                          xytext_offset=(.05, 20),
                          color='white'),
                     dict(txt='b',
                          rc=(5, 15),
                          xytext_offset=(.15, 10),
                          color='white')],

                    [dict(txt='c',
                          rc=(5, 15),
                          xytext_offset=(.15, 10),
                          color='white'),
                     dict(txt='d',
                          rc=(15, 5),
                          xytext_offset=(.15, 10),
                          color='white')],

                    [dict(
                        txt='e',
                        rc=(5, 15),
                        xytext_offset=(.15, 10),
                        color='white')],
                ],
            },
        })

if args.pbumps_threshold_sweep or args.all:
    env.register_plotter(noisefigs.plotters.MainIsBumpPlotter)

if args.bump_examples or args.all:
    env.register_plotter(noisefigs.plotters.BumpExamplePlotter)

if args.bump_examples_isbump or args.all:
    env.register_plotter(noisefigs.plotters.IsBumpExamplePlotter)

env.plot()

