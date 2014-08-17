#!/usr/bin/env python
from __future__ import absolute_import, print_function

from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

import config

parser = flagparse.FlagParser()
parser.add_flag('--rastersFlag')
parser.add_flag('--rates')
parser.add_flag('--maxFRSweeps')
parser.add_flag('--maxFRGridsProbability')
parser.add_flag('--maxFRGridsScatter')
parser.add_flag('--PSeizureGridsProbability')
parser.add_flag('--PSeizureGridsScatter')
args = parser.parse_args()

env = NoiseEnvironment(user_config=config.get_config())

if args.rastersFlag or args.all:
    env.register_plotter(noisefigs.plotters.EIRasterPlotter)

if args.rates or args.all:
    env.register_plotter(noisefigs.plotters.EIRatePlotter)

if args.maxFRSweeps or args.all:
    env.register_plotter(noisefigs.plotters.MaxPopulationFRSweepsPlotter)

if args.maxFRGridsProbability or args.all:
    env.register_plotter(noisefigs.plotters.MaxFRGridsProbabilityPlotter)

if args.maxFRGridsScatter or args.all:
    env.register_plotter(noisefigs.plotters.MaxFRGridsScatterAllPlotter)

if args.PSeizureGridsProbability or args.all:
    env.register_plotter(noisefigs.plotters.PSeizureGridsProbabilityPlotter)

if args.PSeizureGridsScatter or args.all:
    env.register_plotter(noisefigs.plotters.PSeizureGridsScatterAllPlotter)

env.plot()

