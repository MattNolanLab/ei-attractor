'''Noise related submission templates.'''
from __future__ import absolute_import, print_function, division

from grid_cell_model.submitting.factory import SubmitterFactory
from grid_cell_model.submitting.arguments import ArgumentCreator
from grid_cell_model.submitting.noise import (ParameterSweepParser,
                                              SubmissionParserBase)


class SimulationTemplate(object):
    '''A template for submitting a simulation.

    Parameters
    ----------
    app_name : str
        Path to the application being run.
    parser : ArgumentParser
        The command line argument parser.
    default_parameters : dict
        Contains the default parameters
    user_parameters : dict, optional
        User-defined parameters that will be merged into the default parameters
        (and will overwrite them).
    '''
    def __init__(self, app_name, parser, default_parameters, user_parameters=None):
        self._parser = parser
        self._o = None
        self._app_name = app_name
        self._dp = default_parameters
        self._up = {} if user_parameters is None else user_parameters

    @property
    def options(self):
        '''Command line arguments.'''
        return self._o

    @property
    def parser(self):
        '''Return the command line argument parser.'''
        return self._parser

    def parse_args(self):
        '''Parse the command line arguments.'''
        if self._o is None:
            self._o = self.parser.parse_args()

    def update_user_parameters(self, new_parameters):
        '''Change the user-defined parameters if necessary.'''
        self._up.update(new_parameters)

    def run(self):
        '''Run the simulations.'''
        raise NotImplementedError()


class BasicNoiseSimulation(SimulationTemplate):
    '''A simulation setup that runs only a single simulation per noise
    level.'''
    def __init__(self, app_name, default_parameters, user_parameters=None):
        super(BasicNoiseSimulation, self).__init__(app_name,
                                                   SubmissionParserBase(),
                                                   default_parameters,
                                                   user_parameters)

    def run(self):
        if self.parser.options is None:
            self.parser.parse_args()
        o = self.parser.options

        for noise_sigma in self.parser.noise_sigmas:
            p = self._dp.copy()
            p['noise_sigma'] = noise_sigma  # pA
            p['time']        = 10e3 if o.time is None else o.time  # ms

            p['nthreads']    = 1
            p['ntrials']     = o.ntrials
            p['verbosity']   = o.verbosity

            p.update(self._up)

            # Submitting
            ac = ArgumentCreator(p, printout=True)

            simLabel    = '{0}pA'.format(int(p['noise_sigma']))
            numRepeat   = 1
            submitter = SubmitterFactory.getSubmitter(ac,
                                                      self._app_name,
                                                      envType=o.env,
                                                      rtLimit=o.rtLimit,
                                                      output_dir=o.where,
                                                      label=simLabel,
                                                      blocking=True,
                                                      timePrefix=False,
                                                      numCPU=o.nCPU)
            ac.setOption('output_dir', submitter.outputDir())
            startJobNum = 0
            submitter.submitAll(startJobNum, numRepeat, dry_run=o.dry_run)


class ParameterSweep(SimulationTemplate):
    '''One or two parameter sweep exploration.'''

    def __init__(self, app_name, default_parameters, user_parameters=None):
        super(ParameterSweep, self).__init__(app_name, ParameterSweepParser(),
                                             default_parameters,
                                             user_parameters)
        self._slope_selector = None

    def set_bump_slope_selector(self, selector):
        '''When submitting, also include the bump slope data with the network.

        This should be used pretty much only with the grid field simulations
        because all other types of simulations (stationary, velocity, etc.)
        will ignore those options.

        Parameters
        ----------
        selector : :class:`~grid_cell_model.submitting.noise.slopes.SlopeSelector`
            Bump slope selector. This will determine the specific selector,
            e.g. for the default grid fields simulations/no
            velocity/probabilistic connections, etc.
        '''
        self._slope_selector = selector

    def run(self):
        '''Run the parameter sweep.'''
        if self.parser.options is None:
            self.parser.parse_args()
        o = self.parser.options

        for noise_sigma in self.parser.noise_sigmas:
            p = self._dp.copy()
            p['noise_sigma'] = noise_sigma  # pA
            p['time']        = 10e3 if o.time is None else o.time  # ms

            p['nthreads']    = 1
            p['ntrials']     = o.ntrials
            p['verbosity']   = o.verbosity

            p.update(self._up)

            # Submitting
            iterparams = self.parser.iter_params
            if self._slope_selector is not None:
                iterparams.update({
                    'bumpCurrentSlope' : self._slope_selector.get_slopes(
                        noise_sigma)})
            ac = ArgumentCreator(p, printout=True)
            ac.insertDict(iterparams, mult=False)

            ###################################################################
            simLabel    = '{0}pA'.format(int(p['noise_sigma']))
            numRepeat   = 1
            submitter = SubmitterFactory.getSubmitter(ac,
                                                      self._app_name,
                                                      envType=o.env,
                                                      rtLimit=o.rtLimit,
                                                      output_dir=o.where,
                                                      label=simLabel,
                                                      blocking=True,
                                                      timePrefix=False,
                                                      numCPU=o.nCPU)
            ac.setOption('output_dir', submitter.outputDir())
            startJobNum = 0
            submitter.submitAll(startJobNum, numRepeat, dry_run=o.dry_run,
                                filter=self.parser.flat_filter)
            submitter.saveIterParams(iterparams, self.parser.options.param_list,
                                     self.parser.dimensions, dry_run=o.dry_run)

