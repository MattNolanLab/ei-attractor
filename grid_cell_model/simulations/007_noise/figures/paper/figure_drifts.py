#!/usr/bin/env python
from __future__ import absolute_import, print_function

from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

parser = flagparse.FlagParser()
parser.add_flag('--bumpDriftSweep')
parser.add_flag('--bumpDiffAtInitSweep')
parser.add_flag('--bumpDiffResetSweep')
args = parser.parse_args()

env = NoiseEnvironment()

if args.bumpDriftSweep or args.all:
    env.register_plotter(noisefigs.plotters.BumpDriftAtTimePlotter)

if args.bumpDiffAtInitSweep or args.all:
    env.register_plotter(noisefigs.plotters.BumpDiffAtInitPlotter)

if args.bumpDiffResetSweep or args.all:
    env.register_plotter(noisefigs.plotters.BumpDiffResetPlotter)

env.plot()

