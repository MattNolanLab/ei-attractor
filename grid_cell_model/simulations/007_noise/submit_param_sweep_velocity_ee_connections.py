#!/usr/bin/env python
'''
Bump velocity gain estimation - E->E connections

In these simulations, we add E->E connections (with a Gaussian profile). All
other parameters are from the default network.

These simulations do the following:
    - Bump is initialized, in the beginning [0, theta_start_t], by a very
      strong place cell input.
    - The bump moves, with a constant speed, following **vertical** direction.
      for the duration of the simulation.
    - The speed is controlled by IvelMax and dIvel parameters (currently
      [0, 100] pA, with a step of 10 pA
    - At the end of each run, spikes from E and I population are exported to
      the output file.

.. note::
    When simulating on a machine with a run time limit, use a limit around
    08h:00:00
'''
from __future__ import absolute_import, print_function, division

from grid_cell_model.submitting.noise.templates import ParameterSweep
from default_params import defaultParameters as dp

sweep = ParameterSweep('simulation_velocity.py', dp)

p = {}
p['master_seed'] = 123456

p['EI_flat'] = 0
p['IE_flat'] = 0
p['use_EE']  = 1
p['g_EE_total'] = 510.      # nS
p['pEE_sigma'] = 0.5 / 6

p['IvelMax']                = 100   # pA
p['dIvel']                  = 10    # pA

sweep.update_user_parameters(p)
sweep.run()
