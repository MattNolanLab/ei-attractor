
'''Configuration file for NMDA voltage dependence.
In this config, C_Mg = 0.1 mM
'''
from __future__ import absolute_import, print_function

import matplotlib.ticker as ti

def get_config():
    config = {
        'fname_prefix': 'nmda_0p1_mM_',

        'grids_data_root': None,
        'bump_data_root':  'simulation_data/nmda_vm_dep/C_Mg_0p1_mM/gamma_bump',
        'vel_data_root':   None,

        'noise_sigmas': [0, 150, 300],


        #'0p1mM': {
        #    'grids_data_root': None,
        #    'bump_data_root':  'simulation_data/nmda_vm_dep/C_Mg_0p1_mM/gamma_bump',
        #    'vel_data_root':   None,
        #},

        'GammaSweepsPlotter': {
            'plot_grid_contours': [0, 0, 0],
        },

        'MainBumpFormationPlotter': {
            'plot_grid_contours': [0, 0, 0],
        },

    }

    return config
