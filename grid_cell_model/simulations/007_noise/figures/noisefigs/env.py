from __future__ import absolute_import, print_function

import logging

import matplotlib as mpl

from .EI_plotting.base import NoiseDataSpaces
from . import defaultconfig as defc

logger = logging.getLogger(__name__)

class NoiseEnvironment(object):

    def __init__(self, user_config=None, param_spaces=None):
        self._config = self._load_default_config()
        self._merge_user_config(user_config)

        self._plotters = []

        if param_spaces is None:
            self.ps = self._get_default_param_spaces()
        else:
            self.ps = param_spaces

        # self.config must be initialised before this
        self._init_mpl()

    @property
    def config(self):
        return self._config

    def _load_default_config(self):
        return defc.get_config()

    def _get_default_param_spaces(self):
        roots = NoiseDataSpaces.Roots(
            self.config['bump_data_root'],
            self.config['vel_data_root'],
            self.config['grids_data_root'],
            constPos=self.config['const_pos_data_root'])
        return NoiseDataSpaces(roots, self.config['even_shape'],
                               self.config['noise_sigmas'])

    def _merge_user_config(self, user_config):
        '''Update user configuration on top of defaults.'''
        self._user_config = user_config
        if user_config is None:
            self._user_config = {}
        else:
            logger.info('Merging user configuration.')
            logger.debug('    User configuration: %s', self._user_config)
            logger.debug('    Final configuration: %s', self._config)
        self._config.merge(self._user_config)


    def _init_mpl(self):
        mpl.rcParams.update(self.config['mpl'])

    def register_plotter(self, plotter):
        if plotter not in self._plotters:
            self._plotters.append(plotter)

    def plot(self):
        for plot_cls in self._plotters:
            plotter = plot_cls(self.config, self)
            plotter.run_all()
