#!/usr/bin/env python
'''Submit job(s) to the cluster/workstation: grid field parameter sweeps with
extrace place cells connected to I cells.

These simulations do the following:
    - Bump is initialized, in the beginning [0, theta_start_t], by a very
      strong place cell input.
    - The network receives velocity input for the rest of the simulation.
    - The bump moves around and (should) track the position of the animal.
    - Spikes from the whole E population and some cells from the I population
      are then exported to the output file.

.. note::
    When simulating on a machine with a run time limit, use a limit around
    08h:00:00
'''
from __future__ import absolute_import, print_function, division

import numpy as np
from grid_cell_model.submitting.noise.templates import ParameterSweep
from grid_cell_model.submitting.noise.slopes import PickedDefaultSelector
from default_params import defaultParameters as dp

sweep = ParameterSweep('simulation_grids.py', dp)
parser = sweep.parser
parser.add_argument('--g_AMPA_total', type=float)
parser.add_argument('--g_AMPA_row', type=float, help='Row index into bump slope data.')
parser.add_argument('--g_GABA_total', type=float)
parser.add_argument('--g_GABA_col', type=float, help='Column index into bump slope data.')
parser.add_argument('--ipc_field_std', type=float, help='I PC field std. dev..')
parser.add_argument('--ipc_max_rate', type=float)
parser.parse_args()
o = parser.options

sweep.set_bump_slope_selector(PickedDefaultSelector('bump_slope_data',
                                                    -np.infty,
                                                    o.g_AMPA_row,
                                                    o.g_GABA_col,
                                                    parser.dimensions))

p = {}
p['master_seed']      = 123456
p['g_AMPA_total'] = o.g_AMPA_total
p['g_GABA_total'] = o.g_GABA_total

p['velON']            = 1
p['pcON']             = 1
p['constantPosition'] = 0

p['ipc_ON'] = 1
if parser.options.ipc_field_std is not None:
    p['ipc_field_std'] = parser.options.ipc_field_std
if parser.options.ipc_max_rate is not None:
    p['ipc_max_rate'] = parser.options.ipc_max_rate

sweep.update_user_parameters(p)
sweep.run()
