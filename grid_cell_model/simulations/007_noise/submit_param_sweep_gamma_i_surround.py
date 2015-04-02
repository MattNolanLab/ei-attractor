#!/usr/bin/env python
'''Parameter sweeps of the I-surround network.

2D parameter sweep that simulates a stationary bump and records spiking
activity and synaptic currents from selected neurons.'''
from __future__ import absolute_import, print_function, division

from grid_cell_model.submitting.noise.templates import ParameterSweep
from default_params import defaultParameters as dp

sweep = ParameterSweep('../common/simulation_stationary.py', dp)
parser = sweep.parser
parser.add_argument('--Iext_e_theta', type=float, help='Cosine theta input to E cells(pA).')
parser.add_argument('--Iext_i_theta', type=float, help='Cosine theta input to I cells (pA).')
parser.parse_args()

p = {}
p['master_seed'] = 123456
if parser.options.Iext_e_theta is not None:
    p['Iext_e_theta'] = parser.options.Iext_e_theta
if parser.options.Iext_i_theta is not None:
    p['Iext_i_theta'] = parser.options.Iext_i_theta

p['AMPA_gaussian'] = 1

sweep.update_user_parameters(p)
sweep.run()
