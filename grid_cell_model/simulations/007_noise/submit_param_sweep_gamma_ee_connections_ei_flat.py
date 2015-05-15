#!/usr/bin/env python
'''Parameter sweeps of network with extra E->E connections and flat E-I
profiles.

2D parameter sweep that simulates a stationary bump and records spiking
activity and synaptic currents from selected neurons.'''
from __future__ import absolute_import, print_function, division

from grid_cell_model.submitting.noise.templates import ParameterSweep
from default_params import defaultParameters as dp

sweep = ParameterSweep('../common/simulation_stationary.py', dp)
parser = sweep.parser
parser.add_argument('--g_AMPA_total', type=float,
                    help='Strength of E->I synapses')
parser.add_argument('--g_GABA_total', type=float,
                    help='Strength of I->E synapses')
parser.add_argument('--g_EE_total', type=float,
                    help='Strength of E->E synapses')
parser.add_argument('--pEE_sigma', type=float, help='Width of E->E synapses')
sweep.parse_args()
o = sweep.options

p = {}
p['master_seed'] = 123456

p['EI_flat'] = 1
p['IE_flat'] = 1
p['use_EE']  = 1

if o.g_AMPA_total is not None:
    p['g_AMPA_total'] = o.g_AMPA_total
if o.g_GABA_total is not None:
    p['g_GABA_total'] = o.g_GABA_total
if o.g_EE_total is not None:
    p['g_EE_total'] = o.g_EE_total
if o.pEE_sigma is not None:
    p['pEE_sigma'] = o.pEE_sigma

sweep.update_user_parameters(p)
sweep.run()
