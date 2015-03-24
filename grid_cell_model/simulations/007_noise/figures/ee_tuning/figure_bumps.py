#!/usr/bin/env python
from __future__ import absolute_import, print_function
from copy import deepcopy

from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment

import config

parser = flagparse.FlagParser()
parser.add_flag('--pbumps')
parser.add_flag('--pbumps_narrow_sigma')
args = parser.parse_args()


env = NoiseEnvironment(user_config=config.get_config())

if args.pbumps or args.all:
    env.register_plotter(noisefigs.plotters.GEProfileWidthBumpPlotter)

if args.pbumps_narrow_sigma or args.all:
    new_config = deepcopy(config.get_config())
    new_config['bump_data_root'] = ('simulation_data/pAMPA_sigma/'
                                    'sigma_sweep_narrow_range')
    narrow_env = NoiseEnvironment(user_config=new_config)
    narrow_env.register_class(
        noisefigs.plotters.GEProfileWidthBumpPlotter,
        config={
            'GEProfileWidthBumpPlotter': {
                'fname': "bumps_Pbumps_gE_pAMPA_sigma_narrow{ns}.pdf"
            },
        })

env.plot()
narrow_env.plot()
