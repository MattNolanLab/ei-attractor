from __future__ import absolute_import, print_function

import logging
import pprint

from configobj import ConfigObj

pp = pprint.PrettyPrinter(indent=2)
logger = logging.getLogger(__name__)


class Generator(object):
    def __init__(self, cls, config, obj_args, obj_kwargs):
        self.cls = cls
        self.config = config
        self.obj_args = obj_args
        self.obj_kwargs = obj_kwargs


class Environment(object):
    '''Environment that runs all the computations.

    Parameters
    ----------
    config : ConfigObj
        A configuration file.
    '''
    def __init__(self, config=None):
        if config is None:
            self._config = ConfigObj()
        else:
            self._config = config
        self._generators = []

    @property
    def config(self):
        '''Environment configuration.'''
        return self._config

    def _get_config_copy(self):
        '''Get a copy of the configuration.'''
        new_config = ConfigObj()
        new_config.merge(self.config.dict())
        return new_config

    def register_class(self, cls, config=None, merge_in=True, *args, **kwargs):
        '''Register computation class with the environment.
        This creates a list of computation classes together with their
        configuration files. These will be then instantiated during a call to
        the :meth:`~run` method.

        Parameters
        ----------
        cls : FigurePlotter
            Plotter class to register
        *args : any positional arguments to pass on to the computation class.
            The class ``cls`` will receive the arguments in the form of
            ```cls(*args, *config_args, **kwargs)``` and must implement its
            ``__init__`` method accordingly.
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

        self._generators.append(Generator(cls, obj_config, args, kwargs))

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
            arg_list = g.obj_args + (g.config, self)
            instance = g.cls(*arg_list, **g.obj_kwargs)
            instance.run_all()


    def plot(self):
        '''Run everything.

        This runs all registered objects' computation, not just the plotters.
        Only for compatibility with older code.
        '''
        self.run()


