'''Noise-related submission code.'''
from __future__ import absolute_import, print_function, division

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
