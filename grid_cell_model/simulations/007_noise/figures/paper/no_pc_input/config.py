
'''Configuration file for the noise paper: no place cell input.'''
from __future__ import absolute_import, print_function


def get_config():
    return _config


_config = {
    'grids_data_root':      'simulation_data/submission/grids_no_pc_input',
    'bump_data_root':       None,
    'vel_data_root':        None,
    'const_pos_data_root':  None,
    'singleDataRoot':       None,

    'GridSweepsPlotter': {
        'vmin': -.57,
        'vmax': .66,
    },
}
