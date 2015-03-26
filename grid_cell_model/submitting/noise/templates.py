'''Noise related submission templates.'''
from __future__ import absolute_import, print_function, division

import numpy as np
from grid_cell_model.submitting.factory import SubmitterFactory
from grid_cell_model.submitting.arguments import ArgumentCreator
from grid_cell_model.submitting.noise import ParameterSweepParser
from grid_cell_model.submitting.flagparse import positive_int


class ParameterSweep(object):
    '''One or two parameter sweep exploration.'''

    def __init__(self, app_name, default_parameters, user_parameters=None):
        self._parser = ParameterSweepParser()
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
        '''Parser the command line arguments.'''
        if self._o is None:
            self._o = self.parser.parse_args()

    def update_user_parameters(self, new_parameters):
        '''Change the user-defined parameters if necessary.'''
        self._up.update(new_parameters)

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
            ac = ArgumentCreator(p, printout=True)
            ac.insertDict(iterparams, mult=False)

            ###########################################################################
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

