#!/usr/bin/env python
from __future__ import absolute_import, print_function

from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

parser = flagparse.FlagParser()
parser.add_flag('--rastersFlag')
parser.add_flag('--rates')
parser.add_flag('--maxFRSweeps')
args = parser.parse_args()

env = NoiseEnvironment()

if args.rastersFlag or args.all:
    env.register_plotter(noisefigs.plotters.EIRasterPlotter)

if args.rates or args.all:
    env.register_plotter(noisefigs.plotters.EIRatePlotter)

if args.maxFRSweeps or args.all:
    env.register_plotter(noisefigs.plotters.MaxPopulationFRSweepsPlotter)


env.plot()

