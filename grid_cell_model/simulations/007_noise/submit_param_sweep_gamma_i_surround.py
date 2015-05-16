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

sweep.update_user_parameters(p)
sweep.run()
