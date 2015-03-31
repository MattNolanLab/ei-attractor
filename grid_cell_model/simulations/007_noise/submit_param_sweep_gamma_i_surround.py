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

sweep.update_user_parameters(p)
sweep.run()
