
'''Configuration file for the noise paper: no place cell input.'''
from __future__ import absolute_import, print_function

from matplotlib import ticker as ti

def get_config():
    return _config


_config = {
    'grids_data_root':      'simulation_data/main_network/grids_no_velocity',
    'bump_data_root':       None,
    'vel_data_root':        None,
    'const_pos_data_root':  None,
    'singleDataRoot':       None,

    'GridSweepsPlotter': {
        'vmin': -.67,
        'vmax': .42,

        'cbar_kw': {
            'ticks': ti.MultipleLocator(0.25),
        },
    },
}
