#!/usr/bin/env python
'''Parameter sweeps in EI-structured networks with added EE connections.

This simulation is meant to test the width of the E->I connectivity when EE
connections are added.

Run a network for a short time and save data. These simulations are meant to be
run for interactive tuning of the network.
'''
from __future__ import absolute_import, print_function, division

from grid_cell_model.submitting.noise.templates import ParameterSweep
from default_params import defaultParameters as dp

sweep = ParameterSweep('../common/simulation_stationary.py', dp)

p = {}
p['master_seed'] = 123456

p['EI_flat'] = 0
p['IE_flat'] = 0
p['use_EE']  = 1

p['g_AMPA_total'] = 2040.    # nS
p['g_GABA_total'] = 2160.    # nS
p['g_EE_total']   = 510.     # nS
p['pEE_sigma']    = 0.05     # nS

sweep.update_user_parameters(p)
sweep.run()
