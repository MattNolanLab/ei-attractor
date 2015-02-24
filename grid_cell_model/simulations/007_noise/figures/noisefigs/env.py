from __future__ import absolute_import, print_function

import logging
import pprint

import matplotlib as mpl
from configobj import ConfigObj
from simtools.plotting.env import Environment

from .EI_plotting.base import NoiseDataSpaces
from . import defaultconfig as defc

pp = pprint.PrettyPrinter(indent=2)
logger = logging.getLogger(__name__)


class NoiseEnvironment(Environment):
    '''Plotting environment for noise data.
    Parameters
    ----------
    user_config : dict or ConfigObj
        A user-defined config that will be merged into the default config.
    param_spaces
        Do not use this really...
    '''
    def __init__(self, user_config=None, param_spaces=None):
        self._user_config = user_config
        def_config = self._load_default_config()
        final_config = self._merge_user_config(def_config, user_config)
        super(NoiseEnvironment, self).__init__(final_config)

        # This must be after configs are initialised. They are used to determine
        # the parameter spaces
        if param_spaces is None:
            self.ps = self._get_default_param_spaces()
        else:
            self.ps = param_spaces

        # self.config must be initialised before this
        self._init_mpl()

    def _load_default_config(self):
        '''Load a default configuration for the Noise project.'''
        return defc.get_config()

    def _get_default_param_spaces(self):
        roots = NoiseDataSpaces.Roots(
            self.config['bump_data_root'],
            self.config['vel_data_root'],
            self.config['grids_data_root'],
            constPos=self.config['const_pos_data_root'])
        return NoiseDataSpaces(roots, self.config['even_shape'],
                               self.config['noise_sigmas'])

    def _merge_user_config(self, def_config, user_config):
        '''Update user configuration on top of defaults.'''
        new_config = ConfigObj()
        new_config.merge(def_config)
        if user_config is None:
            self._user_config = {}
        else:
            logger.info('Merging user configuration.')
            logger.debug('    User configuration:\n %s',
                         pp.pformat(self._user_config))
        new_config.merge(self._user_config)
        logger.debug('    Final configuration:\n %s',
            pp.pformat(new_config.dict()))
        return new_config

    def _init_mpl(self):
        mpl.rcParams.update(self.config['mpl'])
