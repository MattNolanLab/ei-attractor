#!/usr/bin/env python
from __future__ import absolute_import, print_function

from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

import config
import nmda_0p1_config

parser = flagparse.FlagParser()
parser.add_flag('--gammaSweep')
#parser.add_flag('--threshold')
#parser.add_flag('--freqHist')
#parser.add_flag('--detailed_noise')
parser.add_flag('--examples')
#parser.add_flag('--scatter_all')
#parser.add_flag('--scatter_gamma_pbumps_all')
#parser.add_flag('--gamma_pbumps_prob')
#parser.add_flag('--gamma_grids_prob')
args = parser.parse_args()


env_0mM = NoiseEnvironment(user_config=config.get_config())
env_0p1mM = NoiseEnvironment(user_config=nmda_0p1_config.get_config())

if args.gammaSweep or args.all:
    env_0mM.register_plotter(noisefigs.plotters.GammaSweepsPlotter)
    env_0p1mM.register_plotter(noisefigs.plotters.GammaSweepsPlotter)

if args.examples or args.all:
    env_0mM.register_plotter(noisefigs.plotters.GammaExamplePlotter)
    env_0p1mM.register_plotter(noisefigs.plotters.GammaExamplePlotter)

#if args.scatter_all or args.all:
#    env.register_plotter(noisefigs.plotters.GammaScatterAllPlotter)
#    env.register_plotter(noisefigs.plotters.GammaFreqGridsScatterAllPlotter)

#if args.scatter_gamma_pbumps_all or args.all:
#    env.register_plotter(noisefigs.plotters.GammaScatterPBumpsAllPlotter)

#if args.gamma_pbumps_prob or args.all:
#    env.register_plotter(noisefigs.plotters.GammaPBumpsProbabilityPlotter)
#    env.register_plotter(noisefigs.plotters.GammaFreqPBumpsProbabilityPlotter)
#
#if args.gamma_grids_prob or args.all:
#    env.register_plotter(noisefigs.plotters.GammaGridsProbabilityPlotter)
#    env.register_plotter(noisefigs.plotters.GammaFreqGridsProbabilityPlotter)

env_0mM.plot()
env_0p1mM.plot()
