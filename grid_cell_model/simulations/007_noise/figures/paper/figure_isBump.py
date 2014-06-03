#!/usr/bin/env python
from __future__ import absolute_import, print_function

from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

parser = flagparse.FlagParser()
parser.add_flag('--isBump')
parser.add_flag('--fracTotalHist')
parser.add_flag('--fracTotalSweepAnn')
args = parser.parse_args()

env = NoiseEnvironment()

if args.fracTotalSweepAnn or args.all:
    env.register_plotter(noisefigs.plotters.FracTotalSweepAnnPlotter)

if args.isBump or args.all:
    env.register_plotter(noisefigs.plotters.IsBumpPlotter)

if args.fracTotalHist or args.all:
    env.register_plotter(noisefigs.plotters.FracTotalHistPlotter)

env.plot()

