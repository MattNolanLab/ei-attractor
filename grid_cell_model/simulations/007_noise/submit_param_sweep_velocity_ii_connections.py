#!/usr/bin/env python
'''
Bump velocity gain estimation - probabilistic connections

In these simulations, the synapse weights have a constant value, but the
probability of connection between 2 neurons is generated according to synaptic
profile.

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

p['use_II'] = 1
p['g_II_total'] = 70.               # nS

p['IvelMax']                = 100   # pA
p['dIvel']                  = 10    # pA

sweep.update_user_parameters(p)
sweep.run()
