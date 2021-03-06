#!/usr/bin/env python
'''Grid field parameter sweeps - I-surround configuration (all other config
values however match the original E-surround configuration (main figures)).

These simulations do the following:
    - Bump is initialized, in the beginning [0, theta_start_t], by a very
      strong place cell input.
    - The network receives velocity input for the rest of the simulation.
    - The bump moves around and (should) track the position of the animal.
    - Spikes from the whole E population and some cells from the I population
      are then exported to the output file.

.. note::
    When simulating on a machine with a run time limit, use a limit around
    08h:00:00 or more.
'''
from __future__ import absolute_import, print_function, division

import numpy as np
from grid_cell_model.submitting.noise.templates import ParameterSweep
from grid_cell_model.submitting.noise.slopes import ISurroundOrigSelector
from default_params import defaultParameters as dp

sweep = ParameterSweep('simulation_grids.py', dp)
sweep.set_bump_slope_selector(ISurroundOrigSelector('bump_slope_data',
                                                    -np.infty))

p = {}
p['master_seed']      = 123456
p['velON']            = 1
p['pcON']             = 1
p['constantPosition'] = 0

p['AMPA_gaussian'] = 1

sweep.update_user_parameters(p)
sweep.run()
