#!/usr/bin/env python
from __future__ import absolute_import, print_function

from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

singleDataRoot = 'simulation_data/submission/single_neuron'

parser = flagparse.FlagParser()
parser.add_flag('--maxThetaFRSweeps')
parser.add_flag('--maxThetaFRSweeps_median')
parser.add_flag('--maxThetaFRSweeps_std')
parser.add_flag('--maxThetaFRHist')
parser.add_flag('--seizureProportion')
args = parser.parse_args()

#ps = ds.getDefaultParamSpaces()

env = NoiseEnvironment()

if args.maxThetaFRSweeps or args.all:
    env.register_plotter(noisefigs.plotters.MaxMeanThetaFRSweepPlotter)

if args.seizureProportion or args.all:
    env.register_plotter(noisefigs.plotters.PSeizureSweepPlotter)

if args.maxThetaFRSweeps_std or args.all:
    env.register_plotter(noisefigs.plotters.MaxStdThetaFRSweepPlotter)

if args.maxThetaFRSweeps_median or args.all:
    env.register_plotter(noisefigs.plotters.MaxMedianThetaFRSweepPlotter)

if args.maxThetaFRHist or args.all:
    env.register_plotter(noisefigs.plotters.MaxThetaFRHistPlotter)

env.plot()
