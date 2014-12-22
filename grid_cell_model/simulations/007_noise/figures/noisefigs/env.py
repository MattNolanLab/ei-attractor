from __future__ import absolute_import, print_function

import logging
import pprint

import matplotlib as mpl
from configobj import ConfigObj

from .EI_plotting.base import NoiseDataSpaces
from . import defaultconfig as defc

pp = pprint.PrettyPrinter(indent=2)
logger = logging.getLogger(__name__)

class Plotter(object):
    def __init__(self, cls, config):
        self.cls = cls
        self.config = config

class NoiseEnvironment(object):

    def __init__(self, user_config=None, param_spaces=None):
        self._config = self._load_default_config()
        self._merge_user_config(user_config)

        self._generators = []

        if param_spaces is None:
            self.ps = self._get_default_param_spaces()
        else:
            self.ps = param_spaces

        # self.config must be initialised before this
        self._init_mpl()

    @property
    def config(self):
        return self._config

    def _get_config_copy(self):
        new_config = ConfigObj()
        new_config.merge(self.config.dict())
        return new_config

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
            logger.debug('    User configuration:\n %s',
                         pp.pformat(self._user_config))
        self._config.merge(self._user_config)
        logger.debug('    Final configuration:\n %s',
            pp.pformat(self._config.dict()))

    def _init_mpl(self):
        mpl.rcParams.update(self.config['mpl'])

    def register_class(self, cls, config=None, merge_in=True):
        '''Register computation class with the environment.
        This creates a list of computation classes together with their
        configuration files. These will be then instantiated during a call to
        the :meth:`~run` method.

        Parameters
        ----------
        cls : FigurePlotter
            Plotter class to register
        config  : ConfigObj
            User-defined configuration
        merge_in : bool
            If ``True``, this will be merged with the environment config
            already present. If ``False``, ``config`` will be the only
            configuration file present for this particular object.
        '''
        obj_config = None
        if config is not None:
            if merge_in:
                obj_config = self._get_config_copy()
                obj_config.merge(config)
            else:
                obj_config = config
        else:
            obj_config = self.config

        self._generators.append(Plotter(cls, obj_config))

    def register_plotter(self, plotter_cls, config=None, merge_in=True):
        '''Register plotter class with the environment.
        This creates a list of plotting classes together with their
        configuration files. These will be then instantiated during a call to
        the :meth:`~plot` method.

        Parameters
        ----------
        plotter_cls : FigurePlotter
            Plotter class to register
        config  : ConfigObj
            User-defined configuration
        merge_in : bool
            If ``True``, this will be merged with the environment config
            already present. If ``False``, ``config`` will be the only
            configuration file present for this particular object.
        '''
        self.register_class(plotter_cls, config, merge_in)

    def run(self):
        '''Run everything.'''

        for g in self._generators:
            instance = g.cls(g.config, self)
            instance.run_all()

    def plot(self):
        '''Run everything.

        This runs all registered objects' computation, not just the plotters.
        Only for compatibility with older code.
        '''
        self.run()
