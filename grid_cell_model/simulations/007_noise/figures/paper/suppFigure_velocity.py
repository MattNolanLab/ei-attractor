#!/usr/bin/env python
from __future__ import absolute_import, print_function

from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

singleDataRoot = 'simulation_data/submission/single_neuron'

parser = flagparse.FlagParser()
parser.add_flag('--velSweep')
parser.add_flag('--velFitStdSweep')
parser.add_flag('--velLines')
args = parser.parse_args()

#ps = ds.getDefaultParamSpaces()

env = NoiseEnvironment()

if args.velSweep or args.all:
    env.register_plotter(noisefigs.plotters.VelFitErrSweepPlotter)

if args.velFitStdSweep or args.all:
    env.register_plotter(noisefigs.plotters.VelFitStdSweepPlotter)

if args.velLines or args.all:
    env.register_plotter(noisefigs.plotters.VelLinesPlotter)

env.plot()
