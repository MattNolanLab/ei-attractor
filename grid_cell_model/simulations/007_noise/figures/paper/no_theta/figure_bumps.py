#!/usr/bin/env python
from __future__ import absolute_import, print_function

from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

import config

parser = flagparse.FlagParser()
parser.add_flag('--bumpExamples')
parser.add_flag('--scatter_grids_fracTotal')
parser.add_flag('--fracTotalSweepAnn')
parser.add_flag('--isBump')
parser.add_argument('--expScatter', action="store_true")
args = parser.parse_args()


env = NoiseEnvironment(user_config=config.get_config())

if args.fracTotalSweepAnn or args.all:
    env.register_plotter(noisefigs.plotters.MainBumpFormationPlotter)
    env.register_plotter(
        noisefigs.plotters.MainBumpFormationPlotter,
        config={
            'MainBumpFormationPlotter': {
                'fname_prefix': 'paper_',
                'sigmaTitle': False,
                'cbar_kw': dict(
                    labelpad    = 8,
                    location    = 'right',
                    fraction    = .25,
                    pad         = -.05,
                ),
                'ann': [None]*3,
            },
        })


if args.scatter_grids_fracTotal or args.all:
    env.register_plotter(noisefigs.plotters.MainScatterGridsBumpsPlotter)

if args.isBump or args.all:
    env.register_plotter(noisefigs.plotters.MainIsBumpPlotter)

if args.bumpExamples or args.all:
    env.register_plotter(noisefigs.plotters.BumpExamplePlotter)

env.plot()

