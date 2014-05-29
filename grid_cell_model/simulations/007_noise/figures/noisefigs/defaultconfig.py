from __future__ import absolute_import, print_function

import matplotlib.ticker as ti

exampleRC = ((5, 15), (15, 5))

_default_config = {
    'scale_factor': 1.,

    'iter_list': ['g_AMPA_total', 'g_GABA_total'],

    'output_dir': 'panels/',

    'grids_data_root':     'simulation_data/submission/grids',
    'bump_data_root':      'simulation_data/submission/gamma_bump',
    'vel_data_root':       'simulation_data/submission/velocity',
    'const_pos_data_root': 'simulation_data/submission/const_position',
    'singleDataRoot':      'simulation_data/submission/single_neuron',

    'even_shape': (31, 31),
    'noise_sigmas': [0, 150, 300],

    # Sections
    'mpl': {
        'font.size': 11,
        'pdf.fonttype': 42,
        'mathtext.default': 'regular',
    },

    'sweeps': {
        'fig_size': (3.7, 2.6),         # inches
        'bbox': (0.08, 0.2, .72, .65),  # l, b, w, h
        'transparent': True,
    },

    'grids': {
        'example_rc': ((5, 15), (15, 5)),
        'example_idx': [(5, 15), (5, 15), (5, 15)],  # (row, col)
        'ntrials': 3,
    },
}


def get_config():
    return _default_config


##############################################################################

GridSweepsPlotter_config = {
    'cbar_kw': {
        'label': 'Gridness score',
        'location': 'right',
        'shrink': 0.8,
        'pad': -0.05,
        'ticks': ti.MultipleLocator(0.5),
        'rasterized': True
    },
    'vmin': -0.505,
    'vmax': 1.141,
    'ann': [
        dict(
            txt='b',
            rc=get_config()['grids']['example_rc'][0],
            xytext_offset=(1.5, 1),
            color='black'
        ),
        dict(
            txt='a',
            rc=get_config()['grids']['example_rc'][1],
            xytext_offset=(0.5, 1.5),
            color='black'
        )
    ],
}
_default_config['GridSweepsPlotter'] = GridSweepsPlotter_config

##############################################################################
GridExamplesPlotter_config = {
    'fig_size': (1, 1.2),
    'ax_box': (0.01, 0.01, 0.99, 0.85),  # l, b, r, t
    'transparent': True,
}
_default_config['GridExamplesPlotter'] = GridExamplesPlotter_config

##############################################################################

GridsDiffSweep_config = {
    'cbar_kw': dict(
        label      = '$\Delta_{150 - 0}$(Gridness score)',
        location   = 'right',
        shrink     = 0.8,
        pad        = -0.05,
        ticks      = ti.MultipleLocator(0.5),
        rasterized = True
    )
}
_default_config['GridsDiffSweep'] = GridsDiffSweep_config

##############################################################################

GammaSweepsPlotter_config = {
    'AC_cbar_kw': dict(
        location   = 'left',
        ticks      = ti.MultipleLocator(0.3),
        fraction   = 0.25,
        shrink     = 0.8,
        pad        = .2,
        labelpad   = 8,
        label      = '$1^{st}$ autocorrelation\npeak',
        rasterized = True,
    ),
    'F_cbar_kw': dict(
        location   = 'left',
        ticks      = ti.MultipleLocator(30),
        fraction   = 0.25,
        shrink     = 0.8,
        pad        = .2,
        labelpad   = 8,
        label      = 'Oscillation\nfrequency (Hz)',
        extend     = 'max',
        extendfrac = 0.1,
        rasterized = True
    ),

    'cbar_kw' : {
        'location': 'left',
    }
}
_default_config['GammaSweepsPlotter'] = GammaSweepsPlotter_config

##############################################################################

fracTotalText = 'P(bumps)'

FracTotalSweepAnnPlotter_config = {
    'scale_factor': .8,
    'cbar_kw': dict(
        label       = fracTotalText,
        location    = 'left',
        shrink      = 0.8,
        pad         = 0.25,
        ticks       = ti.MultipleLocator(0.5),
        rasterized  = True
    )
}
_default_config['FracTotalSweepAnnPlotter'] = FracTotalSweepAnnPlotter_config

##############################################################################

_default_config['MainBumpFormationPlotter'] = FracTotalSweepAnnPlotter_config

##############################################################################

_default_config['MainIsBumpPlotter'] = FracTotalSweepAnnPlotter_config

