#!/usr/bin/env python
from __future__ import absolute_import, print_function

from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

import config

parser = flagparse.FlagParser()
parser.add_flag('--pac_example')
parser.add_flag('--burak_conn_weights')
args = parser.parse_args()


env = NoiseEnvironment(user_config=config.get_config())

if args.pac_example or args.all:
    env.register_plotter(noisefigs.plotters.PACExamplePlotter)

if args.burak_conn_weights or args.all:
    env.register_plotter(noisefigs.plotters.Burak2009ConnectionPlotter)

env.plot()

