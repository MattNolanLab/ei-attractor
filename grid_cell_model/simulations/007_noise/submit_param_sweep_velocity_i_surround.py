#!/usr/bin/env python
'''
Bump velocity gain estimation - I-surround configuration with all the other
configuration settings as in the E-surround variant (e.g. theta inputs, etc.)

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
    12h:00:00
'''
from __future__ import absolute_import, print_function, division

from grid_cell_model.submitting.noise.templates import ParameterSweep
from default_params import defaultParameters as dp

sweep = ParameterSweep('simulation_velocity.py', dp)

p = {}
p['master_seed'] = 123456

p['AMPA_gaussian'] = 1

p['IvelMax']                = 100
p['dIvel']                  = 10

sweep.update_user_parameters(p)
sweep.run()
