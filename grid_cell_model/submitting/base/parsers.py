'''Basic interfaces for submission parsers.'''
from __future__ import absolute_import, print_function, division

from grid_cell_model.submitting           import flagparse
from grid_cell_model.submitting.flagparse import positive_int

class GenericSubmissionParser(flagparse.FlagParser):
    '''Parse arguments for submitting generic simulation runs.'''
    def __init__(self, **kwargs):
        super(GenericSubmissionParser, self).__init__(**kwargs)
        self.add_argument('env',     type=str,
                          choices=['workstation', 'cluster'])
        self.add_argument("where",      type=str)
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
    def options(self):
        '''Return the parsed options.'''
        return self._opts

    def parse_args(self, args=None, namespace=None):
        '''Parse the arguments.'''
        self._opts = super(GenericSubmissionParser, self).parse_args(args,
                                                                     namespace)
        return self._opts


