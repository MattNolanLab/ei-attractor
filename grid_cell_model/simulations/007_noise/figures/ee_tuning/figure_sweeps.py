#!/usr/bin/env python
#
'''
Perform analysis on whole 2D data sets.
'''
from __future__ import absolute_import, print_function, division
import time

import matplotlib; matplotlib.use('agg')
from grid_cell_model.parameters import JobTrialSpace2D
from grid_cell_model.submitting import flagparse
import noisefigs
from noisefigs.env import MplEnvironment
from noisefigs.plotters.base import SeparateMultipageSaver

import common
import config

###############################################################################
parser = flagparse.FlagParser()
parser.add_argument('--row',          type=int, required=True)
parser.add_argument('--col',          type=int, required=True)
parser.add_argument('--shapeRows',    type=int, required=True)
parser.add_argument('--shapeCols',    type=int, required=True)
parser.add_argument("--output_dir",   type=str, required=True)
parser.add_argument("--figure_dir",   type=str, required=True)
parser.add_argument("--job_num",      type=int) # unused
parser.add_argument("--type",         type=str, choices=common.allowed_types,
                    required=True, nargs="+")
o = parser.parse_args()

###############################################################################
startT = time.time()

shape = (o.shapeRows, o.shapeCols)
sp = JobTrialSpace2D(shape, o.output_dir)

if common.pop_type in o.type:
    env = MplEnvironment(config=config.get_config())
    env.register_class(
        noisefigs.plotters.PopulationActivityPlotter,
        config={
            'data_root'     : o.output_dir,
            'data_file_name': sp[o.row][o.col].file_name_base,
            'output_dir'    : o.figure_dir,

            'PopulationActivityPlotter': {
                'fname_prefix': 'r%03d_c%03d_' % (o.row, o.col),
                'raster_rect': (.075, 0.35, 0.95, 0.97),
                'fig_saver': SeparateMultipageSaver(None, 'pdf'),
                'fig_size': (10, 6),
                't_limits': (0, 5e3),

                'snapshot_tstep': 4,
                'e_snapshots_rect': (.075, .15, 0.95, 0.25),
                'i_snapshots_rect': (.075, .02, 0.95, 0.12),

            },
        })
    env.plot()


print('Total time: %.3f s' % (time.time() - startT))
