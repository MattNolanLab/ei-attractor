#!/usr/bin/env python
'''Parameter sweeps of E-surround network with extra E->E connections.

2D parameter sweep that simulates a stationary bump and records spiking
activity and synaptic currents from selected neurons.'''
from __future__ import absolute_import, print_function, division

from grid_cell_model.submitting.noise.templates import ParameterSweep
from default_params import defaultParameters as dp

sweep = ParameterSweep('../common/simulation_stationary.py', dp)

p = {}
p['master_seed'] = 123456

p['EI_flat'] = 0
p['IE_flat'] = 0
p['use_EE']  = 1
p['g_EE_total'] = 510.      # nS
p['pEE_sigma'] = 0.5 / 6

sweep.update_user_parameters(p)
sweep.run()
