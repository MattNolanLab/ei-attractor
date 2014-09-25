#!/usr/bin/env python
from __future__ import absolute_import, print_function

from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

import config

parser = flagparse.FlagParser()
parser.add_flag('--rasters')
parser.add_flag('--rastersZoom')
parser.add_flag('--rates')
parser.add_flag('--vel_slope_sweeps')
args = parser.parse_args()

env = NoiseEnvironment(user_config=config.get_config())

if args.rasters or args.all:
    env.register_plotter(noisefigs.plotters.VelocityRasterPlotter)

if args.rates or args.all:
    env.register_plotter(noisefigs.plotters.VelocityRatePlotter)

if args.rastersZoom or args.all:
    env.register_plotter(noisefigs.plotters.VelocityRasterZoomPlotter)

if args.vel_slope_sweeps or args.all:
    env.register_plotter(noisefigs.plotters.VelSlopeSweepPlotter)

env.plot()
