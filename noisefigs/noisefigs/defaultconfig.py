from __future__ import absolute_import, print_function

from copy import deepcopy

import matplotlib.ticker as ti
from configobj import ConfigObj

from .plotters.base import PdfOutputSaver
from .plotters.base import SeparateMultipageSaver


scale_factor = 1.

exampleRC = ((5, 15), (15, 5))
tick_width = 1. * scale_factor
tick_len   = 6. * scale_factor

def get_config():
    return _get_default_config()


def _get_default_config():
    _default_config = ConfigObj()
    _default_config.merge({
        'scale_factor': 1.,

        'iter_list': ['g_AMPA_total', 'g_GABA_total'],

        'output_dir': 'panels/',

        'grids_data_root':      'simulation_data/submission/grids',
        'bump_data_root':       'simulation_data/submission/gamma_bump',
        'vel_data_root':        'simulation_data/submission/velocity',
        'const_pos_data_root':  'simulation_data/submission/const_position',
        'singleDataRoot':       'simulation_data/submission/single_neuron',

        'even_shape': (31, 31),
        'noise_sigmas': [0, 150, 300],

        # Sections
        'mpl': {
            'font.size': 11,
            'pdf.fonttype': 42,
            'mathtext.default': 'regular',
            'font.sans-serif'    : ['Helvetica', 'Avant Garde', 'Computer Modern Sans serif'],

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

        'sweeps': {
            'fig_size': (3.7, 2.6),         # inches
            'bbox': (0.08, 0.2, .72, .65),  # l, b, w, h
            'transparent': True,
            'grid_contours': [.5],
            'contours_kwargs': {
                'hold': True,
                'colors': 'k',
                'linewidths': [1.5]
            },
        },

        'grids': {
            'example_rc': ((5, 15), (15, 5)),
            'example_idx': [(5, 15), (5, 15), (5, 15)],  # (row, col)
            'ntrials': 3,
        },

        'gamma': {
            'example_rc': ((5, 15), (15, 5)),
        },

        'bumps': {
            'n_trials': 5,
        },

        'p_bumps': {
            'frac_total_text' : 'P(bumps)'
        },

        'bump_sigma': {
            'sigma_bump_text': '$\sigma_{bump}^{-1}\ (neurons^{-1})$',
        },

        'seizures': {
            'thetaT': 125.,  # ms
            'sig_dt': .5     # ms
        },

        'vel_rasters': {
            'tLimits': [2e3, 3e3],  # ms
            'trialNum': 0,
            'ylabelPos': -0.22,
        },
    })


    ##############################################################################

    GridSweepsPlotter_config = {
        'cbar': [0, 0, 1],
        'cbar_kw': {
            'label': 'Gridness score',
            'location': 'right',
            'shrink': 0.8,
            'pad': -0.05,
            'ticks': ti.MultipleLocator(0.5),
            'rasterized': True
        },
        'sigma_title': True,
        'vmin': -0.5,
        'vmax': 1.111,
        'xlabel': [None, None, None],
        'xticks': [True, True, True],
        'ylabel': [None, '',    ''],
        'yticks': [True, False, False],
        'ann': [
            dict(
                txt='b',
                rc=_default_config['grids']['example_rc'][0],
                xytext_offset=(1.5, 1),
                color='black'
            ),
            dict(
                txt='a',
                rc=_default_config['grids']['example_rc'][1],
                xytext_offset=(0.5, 1.5),
                color='black'
            )
        ],

        'plot_contours': [0, 0, 0],
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

    GridExampleRectPlotter_config = {
        'cbar_kw': {
            'label': 'Gridness score',
            'location': 'right',
            'shrink': 0.8,
            'pad': -0.05,
            'ticks': ti.MultipleLocator(0.5),
            'rasterized': True
        },
        'vmin': -0.505,
        'vmax': 1.111,

        'fig_saver': PdfOutputSaver(None, 'pdf')
    }
    _default_config['GridExampleRectPlotter'] = GridExampleRectPlotter_config

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

    GridDetailedNoisePlotter_config = {
        'legend':  ['a',  'b'],
        'legend_kwargs': dict(
            loc=(0.8, 1),
            fontsize='small',
            frameon=False,
            numpoints=1,
            handletextpad=0.05,
        )
    }
    _default_config['GridDetailedNoisePlotter'] = GridDetailedNoisePlotter_config

    ##############################################################################

    GridDetailedNoisePlotter_config = {
        'legend':  ['a',  'b'],
        'legend_kwargs': dict(
            loc=(0.8, 1),
            fontsize='small',
            frameon=False,
            numpoints=1,
            handletextpad=0.05,
        )
    }
    _default_config['GridDetailedNoisePlotter'] = GridDetailedNoisePlotter_config

    ##############################################################################

    GammaDetailedNoisePlotter_config = {
        'legend':  ['a',  'b'],
        'legend_kwargs': dict(
            loc=(0.85, 0.7),
            fontsize='small',
            frameon=False,
            numpoints=1,
            handletextpad=0.05,
        )
    }
    _default_config['GammaDetailedNoisePlotter'] = GammaDetailedNoisePlotter_config

    ##############################################################################

    VmExamplesPlotter_config = {
        'fig_size': (2.5, 1.25),
        'ax_rect': (0.01, 0.01, 0.999, 0.6),  # l, b, r, t
    }
    _default_config['VmExamplesPlotter'] = VmExamplesPlotter_config

    ##############################################################################

    ConnectionFunctionPlotter_config = {
        'fig_size': (3, 1.5),
        'bbox_rect': (.2, .25, .95, .75),
        'uniform_random': False,

        'leg1_kwargs': dict(
            loc=(.6, .9),
            frameon=False,
            fontsize='x-small',
            ncol=1
        ),
        'leg2_kwargs': dict(
            loc=(0.45, 1.03),
            frameon=False,
            fontsize='x-small'
        ),
    }
    _default_config['ConnectionFunctionPlotter'] = ConnectionFunctionPlotter_config

    ##############################################################################

    GammaSweepsPlotter_config = {
        'scale_factor': .9,
        'cbar': [1, 0, 0],
        'cbar_kw' : { # This has to match cbar_kw-s below
            'location': 'left',
        },

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
        'AC_xticks': [False]*3,
        'AC_yticks': [1, 0, 0],
        'AC_sigma_title': True,
        'AC_vmin': -0.09,
        'AC_vmax': 0.675,

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
        'F_xticks': [True]*3,
        'F_yticks': [1, 0, 0],
        'F_sigma_title': False,
        'F_vmin': 30,
        'F_vmax': 120,


        'ann': [
            dict(
                txt='b',
                rc=None,
                xytext_offset=(1.5, 0),
                color='white',
            ),
            dict(
                txt='a',
                rc=None,
                xytext_offset=(-.5, 2.),
                color='white',
            ),
        ],

        'plot_grid_contours': [0, 1, 0],
    }
    _default_config['GammaSweepsPlotter'] = GammaSweepsPlotter_config
    tmp = GammaSweepsPlotter_config
    _default_config['GammaSweepsPlotter']['ann'][0]['rc'] = \
        _default_config['gamma']['example_rc'][0]
    _default_config['GammaSweepsPlotter']['ann'][1]['rc'] = \
        _default_config['gamma']['example_rc'][1]

    _default_config['GammaSweepsPlotter'].update({
        'annF': deepcopy(tmp['ann']),
    })

    ##############################################################################

    GammaExamplePlotter_config = {
        # index0: noise_sigma
        # index1: example index
        'xscales': [
            [0, 0, 0],
            [0, 0, 1],
        ],
        'sigma_titles': [
            [0, 0, 0],
            [1, 1, 1],
        ],

        'xscale_kw': dict(
            scaleLen=50,
            x=0.75, y=-0.07,
            size='x-small'
        ),

        'yscale_kw': [[
            dict(
                scaleLen=5,
                unitsText='nA',
                x=.5, y=.1,
                size='x-small'
            ),
            dict(
                scaleLen=0.5,
                unitsText='nA',
                x=.5, y=.05,
                size='x-small'
            ),
            dict(
                scaleLen=0.5,
                unitsText='nA',
                x=.5, y=.05,
                size='x-small'
            )],

            [dict(
                scaleLen=5,
                unitsText='nA',
                x=.5, y=.1,
                size='x-small'
            ),
            dict(
                scaleLen=0.5,
                unitsText='nA',
                x=.5, y=.05,
                size='x-small'
            ),
            dict(
                scaleLen=0.5,
                unitsText='nA',
                x=.55, y=0,
                size='x-small'
            )]],
    }
    _default_config['GammaExamplePlotter'] = GammaExamplePlotter_config

    ##############################################################################

    GammaScatterAllPlotter_config = {
        'fig_size': (4.2, 2),
        'dot_size': 6,
        'legend_kwargs': dict(
            loc=(0.9, 0.4),
            fontsize='small',
            frameon=False,
            numpoints=1,
            title='$\sigma$ (pA)'
        ),
        'bbox_rect': (.1, .35, .95, .85),
        'ylabel': '',
    }
    _default_config['GammaScatterAllPlotter'] = GammaScatterAllPlotter_config

    ##############################################################################

    GammaFreqGridsScatterAllPlotter_config = {
        'fig_size': (4.2, 2),
        'dot_size': 6,
        'legend_kwargs': dict(
            loc=(0.8, 0.4),
            fontsize='small',
            frameon=False,
            numpoints=1,
            title='$\sigma$ (pA)'
        ),
        'bbox_rect': (.1, .35, .95, .85),
        'ylabel': '',
        'yticks': True,
    }
    _default_config['GammaFreqGridsScatterAllPlotter'] = GammaFreqGridsScatterAllPlotter_config

    ##############################################################################

    GammaScatterPBumpsAllPlotter_config = {
        'fig_size': (4.5, 2.6),
        'bbox_rect': (0.3, 0.22, 0.82, 0.95),
        'xlabel': '',
        'legend_kwargs': dict(
            loc=(1.05, 0.5),
            fontsize='small',
            frameon=False,
            title='$\sigma$ (pA)'
        ),
    }
    _default_config['GammaScatterPBumpsAllPlotter'] = GammaScatterPBumpsAllPlotter_config

    ##############################################################################

    GammaPBumpsProbabilityPlotter_config = {
        'fig_size': (2.7, 2.7),         # inches
        'bbox_rect': (0.25, 0.2, 0.95, 0.9),
    }
    _default_config['GammaPBumpsProbabilityPlotter'] = GammaPBumpsProbabilityPlotter_config

    ##############################################################################

    GammaFreqPBumpsProbabilityPlotter_config = {
        'fig_size': (2.7, 2.7),         # inches
        'bbox_rect': (0.25, 0.2, 0.95, 0.9),
    }
    _default_config['GammaFreqPBumpsProbabilityPlotter'] = GammaFreqPBumpsProbabilityPlotter_config

    ##############################################################################

    GammaGridsProbabilityPlotter_config = {
        'fig_size': (2.7, 2.7),         # inches
        'bbox_rect': (0.25, 0.2, 0.95, 0.9),
        'title_size': 'medium',
    }
    _default_config['GammaGridsProbabilityPlotter'] = GammaGridsProbabilityPlotter_config

    ##############################################################################

    GammaFreqGridsProbabilityPlotter_config = {
        'fig_size': (2.7, 2.7),         # inches
        'bbox_rect': (0.25, 0.2, 0.95, 0.9),
        'title_size': 'x-small',
    }
    _default_config['GammaFreqGridsProbabilityPlotter'] = GammaFreqGridsProbabilityPlotter_config

    ##############################################################################

    fracTotalText = _default_config['p_bumps']['frac_total_text']

    FracTotalSweepAnnPlotter_config = {
        'scale_factor': .8,
        'cbar': (1, 0, 0),
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

    MainBumpFormationPlotter_config = {
        'scale_factor': 1.,
        'cbar': [0, 0, 1],
        'cbar_kw': dict(
            label       = "P(bumps)",
            location    = 'right',
            shrink      = 0.8,
            pad         = -.05,
            ticks       = ti.MultipleLocator(0.5),
            rasterized  = True
        ),
        'xticks': [True]*3,
        'plot_grid_contours': [1, 1, 1],
    }
    _default_config['MainBumpFormationPlotter'] = MainBumpFormationPlotter_config

    ##############################################################################

    _default_config['MainIsBumpPlotter'] = FracTotalSweepAnnPlotter_config

    ##############################################################################

    _default_config['IsBumpPlotter'] = FracTotalSweepAnnPlotter_config

    ##############################################################################

    IsBumpExamplePlotter_config = {
        'bumpQualityX': -.9,
    }
    _default_config['IsBumpExamplePlotter'] = IsBumpExamplePlotter_config

    ##############################################################################

    MainScatterGridsBumpsPlotter_config = {
        'fig_size': (4.5, 2.6),
        'bbox_rect': (0.3, 0.22, 0.82, 0.95),
        'xlabel' : '',
        'legend': False,
        'legend_kwargs': dict(
            loc=(1.05, 0.5),
            fontsize='small',
            frameon=False,
            handletextpad=0,
            title='$\sigma$ (pA)'
        ),
    }
    _default_config['MainScatterGridsBumpsPlotter'] = MainScatterGridsBumpsPlotter_config

    ##############################################################################
    BumpDriftAtTimePlotter_config = {
        'scale_factor': .8,
        'cbar_kw': dict(
            label       = 'Average bump drift\n(neurons)',
            location    = 'right',
            shrink      = 0.8,
            pad         = -0.05,
            ticks       = ti.MultipleLocator(10),
            rasterized  = True
        ),
        'plot_grid_contours': [1, 1, 1],
    }
    _default_config['BumpDriftAtTimePlotter'] = BumpDriftAtTimePlotter_config

    ##############################################################################
    BumpDiffAtInitPlotter_config = {
        'cbar_kw': dict(
            label       = 'Distance from init\nposition (neurons)',
            location    = 'right',
            shrink      = 0.8,
            pad         = -0.05,
            ticks       = ti.MultipleLocator(10),
            rasterized  = True
        )
    }
    _default_config['BumpDiffAtInitPlotter'] = BumpDiffAtInitPlotter_config

    ##############################################################################
    BumpDiffResetPlotter_config = {
        'scale_factor': .8,
        'cbar_kw': dict(
            label       = 'Distance from reset\nposition (neurons)',
            location    = 'right',
            shrink      = 0.8,
            pad         = -0.05,
            ticks       = ti.MultipleLocator(5),
            rasterized  = True
        ),
        'plot_grid_contours': [1, 1, 1],
    }
    _default_config['BumpDiffResetPlotter'] = BumpDiffResetPlotter_config

    ##############################################################################

    MaxPopulationFRSweepsPlotter_config = {
        'cbar': [1, 0, 0],
        'cbar_kw': dict(
            label       = "$E-rate_{max}$ (Hz)",
            location    = 'left',
            shrink      = 0.8,
            pad         = 0.25,
            ticks       = ti.MultipleLocator(100),
            rasterized  = True
        ),

        'plot_grid_contours': [1, 1, 1],
        'grid_contours': [.5],
    }
    _default_config['MaxPopulationFRSweepsPlotter'] = MaxPopulationFRSweepsPlotter_config

    ##############################################################################

    BumpSigmaSweepPlotter_config = {
        'cbar': [0, 0, 1],
        'cbar_kw': dict(
            label       = _default_config['bump_sigma']['sigma_bump_text'],
            location    = 'right',
            shrink      = 0.8,
            pad         = -0.05,
            ticks       = ti.MultipleLocator(0.2),
            rasterized  = True
        )
    }
    _default_config['BumpSigmaSweepPlotter'] = BumpSigmaSweepPlotter_config

    ##############################################################################

    BumpExamplePlotter_config = {
        'bbox': (0.01, 0.01, 0.99, 0.82),
    }
    _default_config['BumpExamplePlotter'] = BumpExamplePlotter_config

    ##############################################################################

    EIRasterPlotter_config = {
        'fig_size': (3, 1.9),
        'fig_ext': 'pdf',

        'yticks': [1, 0, 0],
        'ylabelPos': -0.35,

        'scaleBar': [None, None, 25],
        'scaleX': .85,
        'scaleY': -.1,
    }
    _default_config['EIRasterPlotter'] = EIRasterPlotter_config

    ##############################################################################

    EIRatePlotter_config = {
        'fig_size': (3, .65),
        'rateTop': .9,
        'ylabelPos': -0.35,
    }
    _default_config['EIRatePlotter'] = EIRatePlotter_config

    ##############################################################################

    MaxMeanThetaFRSweepPlotter_config = {
        'cbar_kw': dict(
            label       = "max(E rate)/$\\theta$ cycle (Hz)",
            location    = 'left',
            shrink      = 0.8,
            pad         = 0.25,
            ticks       = ti.MultipleLocator(100),
            #ticks       = ti.LogLocator(base=4),
            #format      = ti.LogFormatter(4),
            rasterized  = True
        )
    }
    _default_config['MaxMeanThetaFRSweepPlotter'] = MaxMeanThetaFRSweepPlotter_config

    ##############################################################################

    PSeizureSweepPlotter_config = {
        'FRThreshold': 300,

        'plot_grid_contours': [1, 1, 1],
        'grid_contours': [.5],
    }
    PSeizureSweepPlotter_config.update({
        'cbar_kw': dict(
            label       = "P($E-rate_{{max}}$ > {0})".format(
                            PSeizureSweepPlotter_config['FRThreshold']),
            location    = 'left',
            shrink      = 0.8,
            pad         = 0.25,
            ticks       = ti.MultipleLocator(0.5),
            rasterized  = True
        )
    })
    _default_config['PSeizureSweepPlotter'] = PSeizureSweepPlotter_config

    ##############################################################################

    MaxFRGridsProbabilityPlotter_config = {
        'fig_size': (2.7, 2.7),         # inches
        'scale_factor': .85,
        'bbox_rect': (0.3, 0.22, 0.92, 0.9),
    }
    _default_config['MaxFRGridsProbabilityPlotter'] = MaxFRGridsProbabilityPlotter_config

    ##############################################################################

    PSeizureGridsProbabilityPlotter_config = {
        'FRThreshold': 300,
        'fig_size': (2.7, 2.7),         # inches
        'scale_factor': .85,
        'bbox_rect': (0.3, 0.22, 0.92, 0.9),
    }
    _default_config['PSeizureGridsProbabilityPlotter'] = PSeizureGridsProbabilityPlotter_config

    ##############################################################################

    PSeizureGridsScatterAllPlotter_config = {
        'FRThreshold': 300,
        'fig_size': (2.5, 2.2),         # inches
        'bbox_rect': (0.3, 0.22, 0.92, 0.9),
        'tight_layout_kwargs': {
            'pad': .2,
        },

        'legend_kwargs': dict(
            loc=(0.5, 0.6),
            fontsize='small',
            frameon=False,
            numpoints=1,
            title='$\sigma$ (pA)'
        ),
    }
    _default_config['PSeizureGridsScatterAllPlotter'] = PSeizureGridsScatterAllPlotter_config

    ##############################################################################

    MaxFRGridsScatterAllPlotter_config = {
        'fig_size': (2.5, 2.2),         # inches
        'bbox_rect': (0.3, 0.22, 0.92, 0.9),
        'tight_layout_kwargs': {
            'pad': .2,
        },

        'plot_legend' : False,
        'legend_kwargs': dict(
            loc=(0.6, 0.5),
            fontsize='small',
            frameon=False,
            numpoints=1,
            title='$\sigma$ (pA)'
        ),
    }
    _default_config['MaxFRGridsScatterAllPlotter'] = MaxFRGridsScatterAllPlotter_config

    ##############################################################################

    MaxStdThetaFRSweepPlotter_config = {
        'cbar_kw': dict(
            label       = "max(E rate)/$\\theta$ cycle (Hz)",
            location    = 'left',
            shrink      = 0.8,
            pad         = 0.25,
            ticks       = ti.MaxNLocator(4),
            rasterized  = True
        )
    }
    _default_config['MaxStdThetaFRSweepPlotter'] = MaxStdThetaFRSweepPlotter_config

    ##############################################################################

    _default_config['MaxMedianThetaFRSweepPlotter'] = MaxStdThetaFRSweepPlotter_config

    ##############################################################################

    VelSlopeSweepPlotter_config = {
        'scale_factor': .8,
        'vmin': -.472,
        'vmax': 1.353,
        'cbar': [0, 0, 1],
        'cbar_kw': dict(
            location='right',
            shrink = 0.8,
            pad = -0.1,
            label='Slope\n(neurons/s/pA)',
            ticks=ti.MultipleLocator(0.4),
        ),
        'plot_contours' : [1, 1, 1],

    }
    _default_config['VelSlopeSweepPlotter'] = VelSlopeSweepPlotter_config

    ##############################################################################

    VelFitErrSweepPlotter_config = {
        'scale_factor': .8,
        'cbar': [0, 0, 1],
        'cbar_kw': dict(
            label       = 'Fit error (neurons/s)',
            location    = 'right',
            shrink      = 0.8,
            pad         = -0.1,
            ticks       = ti.MultipleLocator(2),
            rasterized  = True
        ),
        'ylabel': [None, '',    ''],
        'yticks': [1, 0, 0],
        'plot_contours' : [1, 1, 1],
        'vmin': 0,
        'vmax': 11.2,

    }
    _default_config['VelFitErrSweepPlotter'] = VelFitErrSweepPlotter_config

    ##############################################################################

    VelLinesPlotter_config = {
        'scale_factor': .8,
        'fig_size': (3., 2),
        'bbox_rect': (0.4, 0.35, 0.95, 0.65),
        'positions': ((5, 15), (5, 15), (5, 15)),
        'ivel_range': 11,
        'g_ann': False,
    }
    _default_config['VelLinesPlotter'] = VelLinesPlotter_config

    ##############################################################################

    VelFitStdSweepPlotter_config = {
        'scale_factor': .7,
        'cbar_kw': dict(
            location='right',
            label='Mean $\sigma_{spd}$ (neurons/s)',
            shrink = 0.8,
            pad = 0.05,
            ticks=ti.MultipleLocator(5),
            extend='max', extendfrac=0.1
        )
    }
    _default_config['VelFitStdSweepPlotter'] = VelFitStdSweepPlotter_config

    ##############################################################################

    VelocityRasterPlotter_config = {
        'fig_size': (3.75, 2.2),
        'transparent': True,
        'bbox': (0.2, 0.2, 0.99, 0.8)
    }
    _default_config['VelocityRasterPlotter'] = VelocityRasterPlotter_config

    ##############################################################################

    VelocityRatePlotter_config = {
        'fig_size': (3.75, 1),
        'bbox': (.2, .2, .99, 0.70),
        'transparent' : True,
    }
    _default_config['VelocityRatePlotter'] = VelocityRatePlotter_config

    ##############################################################################

    VelocityRasterZoomPlotter_config = {
        'fig_size': (3.75*.75, 1.2),
        'ylabelPos': -0.22,
        'bbox': (0.2, 0.25, 0.99, 0.95),
        'transparent' : True,
    }
    _default_config['VelocityRasterZoomPlotter'] = VelocityRasterZoomPlotter_config

    ##############################################################################

    ThetaSignalPlotter_config = {
        'fig_size': (3, .5),
        'T': .5e3,  # ms
        'bbox': (0, .05, 1., .95),  # l, b, r, t
        'color': (0, 0, 0, .3),
    }
    _default_config['ThetaSignalPlotter'] = ThetaSignalPlotter_config

    ##############################################################################

    PACExamplePlotter_config = {
        'fig_size': (5, 3.5),
        'bbox': (0, .05, 1., .95),  # l, b, r, t
        'letter_xy': (0, 1.),
        'theta_color': 'k',
        'gamma_color': 'b',
    }
    _default_config['PACExamplePlotter'] = PACExamplePlotter_config

    ##############################################################################

    RasterExamplePlotter_config = {
        'fig_size': (6.2, 8.3),
        'sweep_rect' : (.12, .73, .45, .95),
        'cbar_kw': dict(
            label       = "Mean $E-rate_{max}^{\\theta}$ (Hz)",
            location    = 'right',
            shrink      = 0.8,
            pad         = .05,
            ticks       = ti.MultipleLocator(250),
            rasterized  = True
        ),
        'FRThreshold': 300.,
        'ylabelPos': -0.1,
        'markersize': 1.5,
        'plot_ann_txt' : True,
        'theta_color': (0, 0, 0, .3),
        'fig_saver': SeparateMultipageSaver(None, 'pdf')
    }
    _default_config['RasterExamplePlotter'] = RasterExamplePlotter_config

    ##############################################################################

    ScatterGammaGridsSeparatePlotter_config = {
        'fig_size': (5., 6.7),
        #'bbox_rect': (0.12, 0.17, 0.98, 0.92),
    }
    _default_config['ScatterGammaGridsSeparatePlotter'] = ScatterGammaGridsSeparatePlotter_config

    ##############################################################################

    ScatterGammaFGridsSeparatePlotter_config = {
        'fig_size': (5., 6.7),
        #'bbox_rect': (0.12, 0.17, 0.98, 0.92),
    }
    _default_config['ScatterGammaFGridsSeparatePlotter'] = ScatterGammaFGridsSeparatePlotter_config

    ##############################################################################

    GridsPBumpsProbabilityPlotter_config = {
        'fig_size': (2.7, 2.7),         # inches
        'bbox_rect': (0.25, 0.2, 0.95, 0.9),
        'title_size': 'medium',
    }
    _default_config['GridsPBumpsProbabilityPlotter'] = GridsPBumpsProbabilityPlotter_config

    ##############################################################################

    GridBumpScatterPlotter_config = {
        'fig_size': (8.27, 11.69),
        'color_box_width': .165
    }
    GridBumpScatterPlotter_config.update({
        'color_box_coords': {
            'left': 0.14,  # w = 0.165
            'bottom': .85,
            'right': .14 + GridBumpScatterPlotter_config['color_box_width'],
            'top': .95
        }
        #'bbox_rect': (0.12, 0.17, 0.98, 0.92),
    })
    _default_config['GridBumpScatterPlotter'] = GridBumpScatterPlotter_config

    ##############################################################################

    GridSimpleExamplePlotter_config = {
        'fig_size' : (5.4, 2.5),
        'transparent' : True,
        'ns_idx' : 0,
        'rc' : (25, 2),
        'trial_no' : 0,
    }
    _default_config['GridSimpleExamplePlotter'] = GridSimpleExamplePlotter_config

    ##############################################################################

    Burak2009ConnectionPlotter_config = {
        'fig_size' : (2, 2),
    }
    _default_config['Burak2009ConnectionPlotter'] = Burak2009ConnectionPlotter_config

    ##############################################################################

    FRSweepPlotter_config = {
        'scale_factor': .8,
        'cbar_kw' : {
            'location' : 'right',  # This has to match cbar_kw_e and cbar_kw_i
        },

        'plot_grid_contours': [1, 1, 1],

        'cbar_kw_e': {
            'label'      : 'Mean E Firing rate (Hz)',
            'location'   : 'right',
            'shrink'     : 0.8,
            'pad'        : -0.05,
            'ticks'      : ti.LogLocator(subs=[1, 2, 4, 6, 8]),
            'rasterized' : True,
        },

        'cbar_kw_i': {
            'label'      : 'Mean I Firing rate (Hz)',
            'location'   : 'right',
            'shrink'     : 0.8,
            'pad'        : -0.05,
            'ticks'      : ti.LogLocator(subs=[1, 2, 4, 6, 8]),
            'rasterized' : True,
        },

    }
    _default_config['FRSweepPlotter'] = FRSweepPlotter_config

    ##############################################################################

    ScatterGridsFRAllPlotter_config = {
        'fig_size': (4.2, 3),
        'dot_size': 6,
        'legend_kwargs': dict(
            loc=(0.4, 0.6),
            fontsize='small',
            frameon=False,
            numpoints=1,
            title='$\sigma$ (pA)'
        ),
        'bbox_rect': (.2, .2, .95, .95),
        'ylabel': 'Gridness score',
        'yticks': True,
    }
    _default_config['ScatterGridsFRAllPlotter'] = ScatterGridsFRAllPlotter_config

    ###########################################################################
    GridsVelFitErrProbabilityPlotter_config = {
        'fig_size': (2.7, 2.7),         # inches
        'bbox_rect': (0.25, 0.2, 0.95, 0.9),
        'title_size': 'medium',
        'data_range': [[0, 11.2], [-.5, 1.2]],
    }
    _default_config['GridsVelFitErrProbabilityPlotter'] = GridsVelFitErrProbabilityPlotter_config

    ###########################################################################
    _default_config['HighGridScoreFraction'] = {
        'threshold': .5,
    }


    ##############################################################################
    return _default_config
