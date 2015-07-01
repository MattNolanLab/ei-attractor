#!/usr/bin/env python
'''Bump attractor-related plots.'''
from __future__ import absolute_import, print_function
from copy import deepcopy
from os.path import join

import matplotlib.ticker as ti
from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import NoiseEnvironment, MplEnvironment
from noisefigs.plotters.base import SeparateMultipageSaver
from grid_cell_model.parameters import JobTrialSpace2D

import config_standard_gEE_3060 as config

parser = flagparse.FlagParser()
parser.add_flag('--pbumps_sweep')
parser.add_flag('--rasters')
parser.add_flag('--bump_examples')
parser.add_flag('--bump_examples_colorbar')
parser.add_flag('--bump_running_examples')
parser.add_argument('--expScatter', action="store_true")
args = parser.parse_args()


env = NoiseEnvironment(user_config=config.get_config())

if args.pbumps_sweep or args.all:
    env.register_plotter(noisefigs.plotters.MainBumpFormationPlotter,
                         config={
                             'MainBumpFormationPlotter': {
                                'ann': [
                                    None,
                                    [dict(
                                        txt='',
                                        rc=(15, 15),
                                        xytext_offset=(1.5, 1.5),
                                        color='black'
                                    )],
                                    None,
                                ],
                                'vmin': 0,
                                'vmax': 0.98,
                                'cbar_kw': dict(
                                    ticks = ti.MultipleLocator(0.2),
                                ),
                             },
                         })
    env.plot()

if args.rasters or args.all:
    shape = (31, 31)

    output_dir_0 = join('simulation_data',
                      'ee_connections_ei_flat',
                      'standard_sweep_g_EE_3060_pEE_sigma_0_0833',
                      'gamma_bump', '0pA')
    sp_0 = JobTrialSpace2D(shape, output_dir_0)
    new_config = deepcopy(config.get_config())
    env = MplEnvironment(config=new_config)
    env.register_class(
        noisefigs.plotters.PopulationActivityPlotter,
        config={
            'data_root'     : output_dir_0,
            'data_file_name': sp_0[1][22].file_name_base,
            'output_dir'    : 'panels_standard_gEE_3060',

            'PopulationActivityPlotter': {
                'fname_prefix': '0pA_r1_c22_',
                'raster_rect': (.075, 0.35, 0.93, 0.97),
                'fig_saver': SeparateMultipageSaver(None, 'pdf'),
                'fig_size': (8, 6),
                't_limits': (0, 5e3),

                'snapshot_tstep': 4,
                'e_snapshots_rect': (.075, .15, 0.93, 0.25),
                'i_snapshots_rect': (.075, .02, 0.93, 0.12),
            },
        })
    env.plot()


    output_dir = join('simulation_data',
                      'ee_connections_ei_flat',
                      'standard_sweep_g_EE_3060_pEE_sigma_0_0833',
                      'gamma_bump', '150pA')
    sp = JobTrialSpace2D(shape, output_dir)
    new_config = deepcopy(config.get_config())
    env = MplEnvironment(config=new_config)
    env.register_class(
        noisefigs.plotters.PopulationActivityPlotter,
        config={
            'data_root'     : output_dir,
            'data_file_name': sp[20][5].file_name_base,
            'output_dir'    : 'panels_standard_gEE_3060',

            'PopulationActivityPlotter': {
                'fname_prefix': 'r15_c15_',
                'raster_rect': (.075, 0.35, 0.93, 0.97),
                'fig_saver': SeparateMultipageSaver(None, 'pdf'),
                'fig_size': (8, 6),
                't_limits': (0, 5e3),

                'snapshot_tstep': 4,
                'e_snapshots_rect': (.075, .15, 0.93, 0.25),
                'i_snapshots_rect': (.075, .02, 0.93, 0.12),
            },
        })
    env.plot()
