'''Noise-related submission parsers.'''
from __future__ import absolute_import, print_function, division
import itertools

import numpy as np
from grid_cell_model.submitting           import flagparse
from grid_cell_model.submitting.flagparse import positive_int

class SubmissionParserBase(flagparse.FlagParser):
    '''Parse arguments for parameter sweep submission process.'''
    def __init__(self, **kwargs):
        super(SubmissionParserBase, self).__init__(**kwargs)
        self.add_argument('env',     type=str,
                          choices=['workstation', 'cluster'])
        self.add_argument("where",      type=str)
        self.add_argument("--ns",       type=int, choices=[0, 150, 300])
        self.add_argument("--time",     type=float)
        self.add_argument('--ntrials',  type=positive_int, required=True)
        self.add_argument('--rtLimit',  type=str)
        self.add_argument('--printout', type=int, choices=[0, 1], default=1)
        self.add_argument('--nCPU',     type=positive_int, default=1)
        self.add_flag('--dry_run',
                      help='Do no run anything nor save any meta-data')

        self._opts = None

    def _check_opts(self):
        '''Check whether options have been parsed.'''
        if self._opts is None:
            raise RuntimeError("You need to parse the arguments first.")

    @property
    def noise_sigmas(self):
        '''Get noise standard deviations.'''
        self._check_opts()
        ns_all = [0.0, 150.0, 300.0] # pA
        return ns_all if self._opts.ns is None  else [self._opts.ns]

    @property
    def options(self):
        '''Return the parsed options.'''
        return self._opts

    def parse_args(self, args=None, namespace=None):
        '''Parse the arguments.'''
        self._opts = super(SubmissionParserBase, self).parse_args(args,
                                                                  namespace)
        return self._opts


class SubmissionParser(SubmissionParserBase):
    '''Submission parser for parameter sweeps.'''
    def __init__(self, **kwargs):
        super(SubmissionParser, self).__init__(**kwargs)

        self.add_argument('--row',      type=int)
        self.add_argument('--col',      type=int)

    @property
    def rowcol(self):
        '''Get row and column restriction or none if whole sweep.'''
        self._check_opts()
        if self._opts.row is not None:
            rc = (self._opts.row, self._opts.col)
        else:
            rc = None
        return rc

    def parse_args(self, args=None, namespace=None):
        self._opts = super(SubmissionParser, self).parse_args(args, namespace)

        if (self._opts.row is None) ^ (self._opts.col is None):
            raise ValueError("Specify either both --row and --col or None!")

        return self._opts


class ParameterSweepParser(SubmissionParserBase):
    '''A command line parser that allows for exploration of up to 2
    parameters.'''
    def __init__(self, **kwargs):
        super(ParameterSweepParser, self).__init__(**kwargs)

        self.add_argument('param_list',
                          type=str,
                          nargs="+",
                          help=('Explored parameter names. The length will '
                                'determine whether param2_range should be used.'))
        self.add_argument('--range1',
                          type=float,
                          nargs=3,
                          required=True,
                          help=('Parameter 1 range (start, stop, step). Stop '
                                'will be included.'))
        self.add_argument('--range2',
                          type=float,
                          nargs=3,
                          help=('Parameter 2 range (start, stop, step). Stop '
                                'will be included.'))
        self.add_argument('--filter',
                          type=int,
                          nargs="+",
                          help=('Restrict submission to only one value in the '
                                'parameter space.'))

        self.range1 = None
        self.range2 = None

    def parse_args(self, args=None, namespace=None):
        '''Parse the arguments.'''
        super(ParameterSweepParser, self).parse_args(args, namespace)

        # Check consistency
        param_len = len(self.options.param_list)
        if param_len > 2:
            raise ValueError('The parameter label list cannot be longer than '
                             '2 items.')
        if self.options.filter is not None:
            if len(self.options.filter) != param_len:
                raise ValueError('filter length must match the number of '
                                 'explored dimensions.')

        if param_len == 1:
            if self.options.range2 is not None:
                raise ValueError('range2 parameter cannot be used with only '
                                 'one dimensional parameter sweep.')

        if param_len == 2:
            if self.options.range2 is None:
                raise ValueError('2-dimensional sweeps need also the range2 '
                                 'parameter.')

        self._parse_ranges()
        self._check_filter()

        return self.options

    def _check_filter(self):
        '''Check whether the filter has a reasonable range.'''
        msg = 'Filter cannot be greater than the size of the parameter space.'
        if self.filter is not None:
            if self.filter[0] >= self.dimensions[0]:
                raise ValueError(msg)
            if self.n_dims == 2:
                if self.filter[1] >= self.dimensions[1]:
                    raise ValueError(msg)

    @property
    def filter(self):
        '''Return the filtered restriction submission index.'''
        self._check_opts()
        return self._opts.filter

    @property
    def flat_filter(self):
        '''Return a linear index into iterparams for filtering the
        submission.
        '''
        if self.filter is None:
            return None

        if self.n_dims == 1:
            return self.filter[0]
        else:
            return self.filter[0] * self.dimensions[1] + self.filter[1]

    @property
    def n_dims(self):
        '''Return the number of dimensions to explore.'''
        self._check_opts()
        return len(self.options.param_list)

    def _parse_ranges(self):
        '''Parse the range parameters.'''
        o = self.options
        self.range1 = np.arange(o.range1[0], o.range1[1] + o.range1[2],
                                o.range1[2])
        if self.n_dims == 2:
            self.range2 = np.arange(o.range2[0], o.range2[1] + o.range2[2],
                            o.range2[2])

    @property
    def iter_params(self):
        '''Generate the iteration parameters.'''
        iterparams = {}

        if self.n_dims == 1:
            iterparams[self.options.param_list[0]] = np.asanyarray(self.range1)
        elif self.n_dims == 2:
            data = zip(*itertools.product(self.range1, self.range2))
            iterparams[self.options.param_list[0]] = np.asanyarray(data[0])
            iterparams[self.options.param_list[1]] = np.asanyarray(data[1])
        else:
            raise ValueError('This is a bug. n_dims must be 1 or 2. Got %s' %
                             str(self.n_dims))

        return iterparams

    @property
    def dimensions(self):
        '''Return the sizes of the iteration parameters.'''
        dimensions = [len(self.range1)]
        if self.range2 is not None:
            dimensions.append(len(self.range2))
        return dimensions



class SingleParameterSweepParser(SubmissionParserBase):
    '''A submission parser that runs a single parameter sweep in a specified
    range.
    '''
    def __init__(self, **kwargs):
        super(SingleParameterSweepParser, self).__init__(**kwargs)

        self.add_argument('explored_param', type=str,   help='Explored parameter name.')
        self.add_argument('param_start',    type=float, help='Parameter start value')
        self.add_argument('param_stop',     type=float, help='Parameter stop value')
        self.add_argument('param_step',     type=float, help='Parameter step value')
        self.add_argument('--filter',       type=int,   help='Restrict submission to only one value')

    @property
    def filter(self):
        '''Return the filtered restriction submission index.'''
        self._check_opts()
        return self._opts.filter
