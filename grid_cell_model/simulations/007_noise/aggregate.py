'''Aggregation helpers.'''
from grid_cell_model.submitting import flagparse

evenSpacingType = 'even-spacing'
detailedNoiseType = 'detailed-noise'

class AggregationParser(flagparse.FlagParser):
    '''Parse arguments for aggregation scripts.'''
    # Positions
    allowedPositions = ['EI-1_3', 'EI-3_1']
    detailedShape = (31, 9)
    evenShape     = (31, 31)
    allowedTypes = [evenSpacingType, detailedNoiseType]
    ns_all = ['0pA', '150pA', '300pA']

    def __init__(self, **kwargs):
        super(AggregationParser, self).__init__(**kwargs)
        self.add_argument("type", type=str, choices=self.allowedTypes,
                          metavar='type',
                          help=('Type of the aggregation. Can be one of %s' %
                                str(self.allowedTypes)))
        self.add_argument("where",      type=str, help='Root directory')
        self.add_argument("--ns",       type=str, choices=['0pA', '150pA', '300pA'])
        self.add_argument('--ntrials',  type=int)
        self.add_argument('--position', type=str, choices=self.allowedPositions)
        self.add_argument('--noLoadData', action='store_true')
        self.add_argument('--shape', type=int, nargs=2,
                          help='Shape of the 2D parameter space (rows, columns)')

        self._opts = None

    def parse_args(self, args=None, namespace=None):
        '''Parse the arguments.'''
        self._opts = super(AggregationParser, self).parse_args(args, namespace)
        return self._opts

    def _check_opts(self):
        '''Check whether options have been parsed.'''
        if self._opts is None:
            raise RuntimeError("You need to parse the arguments first.")

    @property
    def noise_sigmas(self):
        '''Get noise standard deviations.'''
        self._check_opts()
        if self.options.ns is None:
            return self.ns_all
        else:
            return [self.options.ns]

    @property
    def detailed_positions(self):
        '''Get detailed positions.'''
        if self.options.position is None:
            return self.allowedPositions
        else:
            return [self.options.position]

    @property
    def options(self):
        '''Returned parsed options.'''
        self._check_opts()
        return self._opts

    @property
    def subdirs(self):
        '''Return subdirectories to run through.'''
        if self.options.type == evenSpacingType:
            return self.noise_sigmas
        else:
            return self.detailed_positions

    @property
    def shape(self):
        '''Return the 2D shapes.'''
        if self.options.shape is None:
            if self.options.type == evenSpacingType:
                return self.evenShape
            else:
                return self.detailedShape
        else:
            return (self.options.shape[0], self.options.shape[1])

    @property
    def load_data(self):
        '''Whether to load data or not during aggregation.'''
        return not self.options.noLoadData
