#!/usr/bin/env python
from __future__ import absolute_import, print_function, division
from copy import deepcopy
from os.path import join

import matplotlib.ticker as ti
from grid_cell_model.submitting import flagparse
from grid_cell_model.parameters import JobTrialSpace2D
import noisefigs
from noisefigs.env import NoiseEnvironment, MplEnvironment
from noisefigs.plotters.base import SeparateMultipageSaver

import config

parser = flagparse.FlagParser()
parser.add_flag('--pbumps_gEE_EE_sigma')
parser.add_flag('--rasters_gEE_EE_sigma')
parser.add_flag('--pbumps_gEE_EE_sigma_AMPA_3060_GABA_1020')
parser.add_flag('--rasters_gEE_EE_sigma_AMPA_3060_GABA_1020')
args = parser.parse_args()

if args.pbumps_gEE_EE_sigma or args.all:
    new_config = deepcopy(config.get_config())
    new_config.update({
        'grids_data_root': None,
        'bump_data_root': join('simulation_data', 'submission',
                               'ee_connections_ei_flat',
                               'g_EE_total_vs_pEE_sigma', 'gamma_bump'),
        'vel_data_root':  None,
        'const_pos_data_root': None,
        'singleDataRoot': None,

        'noise_sigmas': [0, 150, 300],
        'even_shape': None,
    })
    env = NoiseEnvironment(user_config=new_config)
    env.register_class(
        noisefigs.plotters.Generic2DPBumpPlotter,
        config={
            'Generic2DPBumpPlotter': {
                'fname': "bumps_Pbumps_gEE_pEE_sigma_{ns}.pdf",
                'normalize_ticks': (True, False),  # (Y, X)
                'normalize_type': ('E', None),
                'xlabel': '$\sigma_{E{\\rightarrow}E}$',
                'ylabel': '$g_{E{\\rightarrow}E}$',
                'bbox': (.2, .17, .85, .9),
                'cbar': [0, 0, 1],
                'ann': [
                    None,
                    [dict(
                        txt='',
                        rc=(20, 5),
                        xytext_offset=(0.03, 1.5),
                        color='black'
                    )],
                    None,
                ],
            },
        })
    env.plot()


if args.rasters_gEE_EE_sigma or args.all:
    shape = (31, 11)
    output_dir = join('simulation_data', 'submission',
                      'ee_connections_ei_flat', 'g_EE_total_vs_pEE_sigma',
                      'gamma_bump', '150pA')
    sp = JobTrialSpace2D(shape, output_dir)
    new_config = deepcopy(config.get_config())
    env = MplEnvironment(config=new_config)
    env.register_class(
        noisefigs.plotters.PopulationActivityPlotter,
        config={
            'data_root'     : output_dir,
            'data_file_name': sp[20][5].file_name_base,
            'output_dir'    : 'panels',

            'PopulationActivityPlotter': {
                'fname_prefix': 'r20_c5_',
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


if args.pbumps_gEE_EE_sigma_AMPA_3060_GABA_1020 or args.all:
    new_config = deepcopy(config.get_config())
    new_config.update({
        'grids_data_root': None,
        'bump_data_root': join('simulation_data', 'submission',
                               'ee_connections_ei_flat',
                               'g_EE_total_vs_pEE_sigma_AMPA_3060_GABA_1020',
                               'gamma_bump'),
        'vel_data_root':  None,
        'const_pos_data_root': None,
        'singleDataRoot': None,

        'output_dir': 'panels_AMPA_3060_GABA_1020',
        'noise_sigmas': [0, 150, 300],
        'even_shape': None,
    })
    env = NoiseEnvironment(user_config=new_config)
    env.register_class(
        noisefigs.plotters.Generic2DPBumpPlotter,
        config={
            'Generic2DPBumpPlotter': {
                'fname': "bumps_Pbumps_gEE_pEE_sigma_{ns}.pdf",
                'normalize_ticks': (True, False),  # (Y, X)
                'normalize_type': ('E', None),
                'xlabel': '$\sigma_{E{\\rightarrow}E}$',
                'ylabel': '$g_{E{\\rightarrow}E}$',
                'bbox': (.2, .17, .85, .9),

                'cbar': [0, 0, 1],
                'ann': [
                    None,
                    [dict(
                        txt='',
                        rc=(15, 5),
                        xytext_offset=(0.03, 1.5),
                        color='black'
                    )],
                    None,
                ],
            },
        })
    env.plot()


if args.rasters_gEE_EE_sigma_AMPA_3060_GABA_1020 or args.all:
    shape = (31, 11)

    output_dir_0 = join('simulation_data', 'submission',
                      'ee_connections_ei_flat',
                      'g_EE_total_vs_pEE_sigma_AMPA_3060_GABA_1020',
                      'gamma_bump', '0pA')
    sp_0 = JobTrialSpace2D(shape, output_dir_0)
    new_config = deepcopy(config.get_config())
    env = MplEnvironment(config=new_config)
    env.register_class(
        noisefigs.plotters.PopulationActivityPlotter,
        config={
            'data_root'     : output_dir_0,
            'data_file_name': sp_0[15][5].file_name_base,
            'output_dir'    : 'panels_AMPA_3060_GABA_1020',

            'PopulationActivityPlotter': {
                'fname_prefix': '0pA_r15_c5_',
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


    output_dir = join('simulation_data', 'submission',
                      'ee_connections_ei_flat',
                      'g_EE_total_vs_pEE_sigma_AMPA_3060_GABA_1020',
                      'gamma_bump', '150pA')
    sp_150 = JobTrialSpace2D(shape, output_dir)
    new_config = deepcopy(config.get_config())
    env = MplEnvironment(config=new_config)
    env.register_class(
        noisefigs.plotters.PopulationActivityPlotter,
        config={
            'data_root'     : output_dir,
            'data_file_name': sp_150[15][5].file_name_base,
            'output_dir'    : 'panels_AMPA_3060_GABA_1020',

            'PopulationActivityPlotter': {
                'fname_prefix': 'r15_c5_',
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
