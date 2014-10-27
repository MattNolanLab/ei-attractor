
'''Configuration file for NMDA voltage dependence.
In this default config, C_Mg = 0 mM
'''
from __future__ import absolute_import, print_function

import matplotlib.ticker as ti

def get_config():
    config = {
        # Defaults for 0 mM C_Mg
        'grids_data_root': None,
        'bump_data_root':  'simulation_data/nmda_vm_dep/C_Mg_0_mM/gamma_bump',
        'vel_data_root':   None,

        'noise_sigmas': [0],


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
