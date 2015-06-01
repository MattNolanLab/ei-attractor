#!/usr/bin/env python
'''Run a 1-variable sweep in which I->I connections are enabled.

Run a network for a short time; save data and plot activity.
'''
from __future__ import absolute_import, print_function, division

from grid_cell_model.submitting.noise.templates import ParameterSweep
from default_params import defaultParameters as dp

sweep = ParameterSweep('../common/simulation_stationary.py', dp)

parser = sweep.parser
parser.add_argument('--nthreads', type=int, default=1,
                    help='Number of simulation threads.')
parser.add_argument('--Ivel', type=float,
                    help='Velocity input (pA). Default is 50 pA.')
parser.add_argument('--g_AMPA_total', type=float,
                    help='Total E->I synapse strength (nS)')
parser.add_argument('--g_GABA_total', type=float,
                    help='Total I->E synapse strength (nS)')
o = parser.parse_args()

p = {}
p['master_seed'] = 123456

p['time']      = 10e3 if o.time is None else o.time  # ms
p['nthreads']  = o.nthreads
p['verbosity'] = o.verbosity
p['Ivel']      = 50. if o.Ivel is None else o.Ivel  # mA
p['use_II']    = 1

p['g_AMPA_total'] = 3060. if o.g_AMPA_total is None else o.g_AMPA_total
p['g_GABA_total'] = 1020. if o.g_GABA_total is None else o.g_GABA_total

sweep.update_user_parameters(p)
sweep.run()
