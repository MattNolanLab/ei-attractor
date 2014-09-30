#!/usr/bin/env python
from __future__ import absolute_import, print_function

from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

import config
import nmda_0p1_config

parser = flagparse.FlagParser()
parser.add_flag('--bumpExamples')
parser.add_flag('--fracTotalSweepAnn')
parser.add_flag('--isBump')
parser.add_argument('--expScatter', action="store_true")
args = parser.parse_args()


env_0mM = NoiseEnvironment(user_config=config.get_config())
env_0p1mM = NoiseEnvironment(user_config=nmda_0p1_config.get_config())

if args.fracTotalSweepAnn or args.all:
    env_0mM.register_plotter(noisefigs.plotters.MainBumpFormationPlotter)
    env_0p1mM.register_plotter(noisefigs.plotters.MainBumpFormationPlotter)

if args.isBump or args.all:
    env_0mM.register_plotter(noisefigs.plotters.MainIsBumpPlotter)

if args.bumpExamples or args.all:
    env_0mM.register_plotter(noisefigs.plotters.BumpExamplePlotter)

env_0mM.plot()
env_0p1mM.plot()
