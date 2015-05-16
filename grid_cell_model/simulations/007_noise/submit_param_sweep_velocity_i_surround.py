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
parser = sweep.parser
parser.add_argument('--Iext_e_theta', type=float, help='Cosine theta input to E cells(pA).')
parser.add_argument('--Iext_e_const', type=float, help='Constant input to E cells(pA).')
parser.add_argument('--Iext_i_theta', type=float, help='Cosine theta input to I cells (pA).')
parser.add_argument('--Iext_i_const', type=float, help='Constant input to I cells(pA).')
parser.add_argument('--g_uni_GABA_frac', type=float,
                    help='Strength or uniform I->E synapses relative to the '
                         'structured I->E.')
parser.add_argument('--g_AMPA_total', type=float)
parser.add_argument('--g_GABA_total', type=float)
parser.parse_args()

p = {}
p['master_seed'] = 123456
if parser.options.Iext_e_theta is not None:
    p['Iext_e_theta'] = parser.options.Iext_e_theta
if parser.options.Iext_e_const is not None:
    p['Iext_e_const'] = parser.options.Iext_e_const
if parser.options.Iext_i_theta is not None:
    p['Iext_i_theta'] = parser.options.Iext_i_theta
if parser.options.Iext_i_const is not None:
    p['Iext_i_const'] = parser.options.Iext_i_const
if parser.options.g_uni_GABA_frac is not None:
    p['g_uni_GABA_frac'] = parser.options.g_uni_GABA_frac
if parser.options.g_AMPA_total is not None:
    p['g_AMPA_total'] = parser.options.g_AMPA_total
if parser.options.g_GABA_total is not None:
    p['g_GABA_total'] = parser.options.g_GABA_total

p['AMPA_gaussian'] = 1

p['IvelMax']                = 100
p['dIvel']                  = 10

sweep.update_user_parameters(p)
sweep.run()
