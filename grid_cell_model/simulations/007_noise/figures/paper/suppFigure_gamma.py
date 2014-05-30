#!/usr/bin/env python
from __future__ import absolute_import, print_function

from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

parser = flagparse.FlagParser()
parser.add_flag('--scatterPlot')
args = parser.parse_args()


env = NoiseEnvironment()

if args.scatterPlot or args.all:
    env.register_plotter(noisefigs.plotters.ScatterGammaGridsSeparatePlotter)

env.plot()
