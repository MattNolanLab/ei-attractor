#!/usr/bin/env python
from __future__ import absolute_import, print_function

from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

parser = flagparse.FlagParser()
parser.add_flag('--bumpSweep')
parser.add_flag('--bumpDriftSweep')
parser.add_flag('--bumpDiffAtInitSweep')
parser.add_flag('--bumpDiffResetSweep')
parser.add_flag('--bumpExamples')
parser.add_flag('--velExamples')
parser.add_flag('--velSweep')
parser.add_flag('--gridness_vs_error')
parser.add_flag('--detailed_noise')
parser.add_flag('--rastersFlag')
parser.add_flag('--rates')
parser.add_flag('--scatter_diff_fracTotal_grids')
parser.add_flag('--scatter_grids_fracTotal')
parser.add_flag('--fracTotalSweepAnn')
parser.add_flag('--isBump')
parser.add_argument('--expScatter', action="store_true")
args = parser.parse_args()


env = NoiseEnvironment()

if args.fracTotalSweepAnn or args.all:
    env.register_plotter(noisefigs.plotters.MainBumpFormationPlotter)

if args.scatter_grids_fracTotal or args.all:
    env.register_plotter(noisefigs.plotters.MainScatterGridsBumpsPlotter)

if args.isBump or args.all:
    env.register_plotter(noisefigs.plotters.MainIsBumpPlotter)

if args.bumpSweep or args.all:
    env.register_plotter(noisefigs.plotters.BumpSigmaSweepPlotter)

if args.bumpExamples or args.all:
    env.register_plotter(noisefigs.plotters.BumpExamplePlotter)

if args.detailed_noise or args.all:
    env.register_plotter(noisefigs.plotters.BumpSigmaDetailedNoisePlotter)

env.plot()

