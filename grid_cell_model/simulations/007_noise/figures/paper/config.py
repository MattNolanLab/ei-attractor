
'''Configuration file for the noise paper.'''
from __future__ import absolute_import, print_function

import matplotlib.ticker as ti
from noisefigs.plotters.base import SeparateMultipageSaver

def get_config():
    return _config


_config = {
    'GridExampleRectPlotter': {
        'fig_saver': SeparateMultipageSaver(None, 'pdf')
    },
}

