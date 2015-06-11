#!/usr/bin/env python
'''I-surround gamma plots with the Pastoll et al. settings.'''
from __future__ import absolute_import, print_function

import copy

import matplotlib
from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

import config_pastoll_pc_weight_3 as config

parser = flagparse.FlagParser()
parser.add_flag('--gamma_sweep')
parser.add_flag('--examples')
parser.add_flag('--scatter_gamma_pbumps_all')
args = parser.parse_args()


env = NoiseEnvironment(user_config=config.get_config())

if args.gamma_sweep or args.all:
    env.register_plotter(noisefigs.plotters.GammaSweepsPlotter)

if args.examples or args.all:
    env.register_plotter(noisefigs.plotters.GammaExamplePlotter)

if args.scatter_gamma_pbumps_all or args.all:
    env.register_plotter(noisefigs.plotters.GammaScatterPBumpsAllPlotter)

env.plot()
