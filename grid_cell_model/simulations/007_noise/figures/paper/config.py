
'''Configuration file for the noise paper.'''
from __future__ import absolute_import, print_function

import matplotlib.ticker as ti

def get_config():
    return _config


_config = {
    'RasterExamplePlotter' : {
        'cbar_kw': dict(
            label       = "max(E rate) (Hz)",
            location    = 'right',
            shrink      = 0.8,
            pad         = -.02,
            ticks       = ti.MultipleLocator(100),
            rasterized  = True
        ),
        'ylabelPos': -0.1,
        'markersize': 1.5,
    }
}

