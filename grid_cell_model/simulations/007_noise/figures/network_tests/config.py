'''Network test configuration file.'''
from __future__ import absolute_import, print_function

from configobj import ConfigObj

from noisefigs.plotters.base import PdfOutputSaver
from noisefigs.plotters.base import SeparateMultipageSaver


scale_factor = 1.

tick_width = 1. * scale_factor
tick_len   = 6. * scale_factor


def get_config():
    '''Return the configuration object.'''
    _default_config = ConfigObj()
    _default_config.merge({
        'scale_factor': scale_factor,

        'data_root'     : 'simulation_data/network_test/10s_4trials/150pA',
        'data_file_name': 'job00000_output.h5',
        'output_dir'    : 'panels/',

        # Sections
        'mpl': {
            'font.size': 11,
            'pdf.fonttype': 42,
            'mathtext.default': 'regular',
            'font.sans-serif': ['Helvetica', 'Avant Garde',
                                'Computer Modern Sans serif'],

            'xtick.major.size'  : tick_len,
            'xtick.major.width' : tick_width,
            'xtick.minor.size'  : tick_len / 2.,
            'xtick.minor.width' : tick_width,
            'xtick.direction'   : 'out',

            'ytick.major.size'  : tick_len,
            'ytick.major.width' : tick_width,
            'ytick.minor.size'  : tick_len / 2.,
            'ytick.minor.width' : tick_width,
            'ytick.direction'   : 'out',
        },

        'PopulationActivityPlotter': {
            'raster_rect': (.075, 0.35, 0.99, 0.97),
            'fig_saver': SeparateMultipageSaver(None, 'pdf'),
            'fig_size': (15, 6),
            't_limits': (0, 10e3),

            'snapshot_tstep': 4,
            'e_snapshots_rect': (.075, .15, 0.99, 0.25),
            'i_snapshots_rect': (.075, .02, 0.99, 0.1),
        },
    })

    ##########################################################################
    return _default_config
