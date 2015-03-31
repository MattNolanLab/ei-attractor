#!/usr/bin/env python
'''Parameter sweeps of the I-surround network.

2D parameter sweep that simulates a stationary bump and records spiking
activity and synaptic currents from selected neurons.'''
from __future__ import absolute_import, print_function, division

from grid_cell_model.submitting.noise.templates import ParameterSweep
from default_params import defaultParameters as dp

sweep = ParameterSweep('../common/simulation_stationary.py', dp)

p = {}
p['master_seed'] = 123456

p['AMPA_gaussian'] = 1
p['Iext_e_theta']  = 650.      # pA
p['Iext_i_theta']  = 50.       # pA
p['g_uni_GABA_frac'] = 0.3125  


sweep.update_user_parameters(p)
sweep.run()
