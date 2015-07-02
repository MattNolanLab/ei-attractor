'''
Simulation parsers with several pre-defined arguments/flags.
'''
import argparse
import logging

__all__ = ['SimulationParser', 'FlagParser', 'FlagRunner']


class SimulationParser(argparse.ArgumentParser):
    '''Parser that allows users to quickly define simulation parameters from
    command line.

    It has one pre-defined parameter: ``--verbosity`` or ``-v``, that specifies
    the logging level for the python logging module.
    '''
    def __init__(self):
        argparse.ArgumentParser.__init__(self)

        self.add_argument(
            '-v', '--verbosity',
            type=str,
            choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
            default='WARNING',
            help='Verbosity level. Corresponds to the python logging module. '
                 'Default is WARNING.'
        )

    def add_flag(self, *args, **kwargs):
        '''Add a flag to the parser.'''
        kwargs['action'] = 'store_true'
        self.add_argument(*args, **kwargs)

    def parse_args(self, args=None, namespace=None):
        args = super(SimulationParser, self).parse_args(args, namespace)
        try:
            level = getattr(logging, args.verbosity)
            logging.basicConfig(level=level)
            logging.getLogger().setLevel(level)
        except AttributeError:
            raise RuntimeError('Something went wrong! Cannot fetch logging '
                               'attribute')
        return args


class FlagParser(SimulationParser):
    '''A parser that parser flags for execution by scripts.

    Parameters
    ----------
    allDefault : bool, optional
        Determines the default value of the ``--all`` option. When any other
        flag is set by :meth:`~add_flag`, ``--all`` will be forced to False.
        This option thus defines the ``--all`` parameter as a reserved keyword
        and it should not be used for any other purpose.
    '''
    def __init__(self, allDefault=True):
        super(FlagParser, self).__init__()

        self.allDefault = allDefault
        self.add_argument(
            '--' + FlagAction.allArg,
            action=FlagAction,
            default=allDefault,
            help='Whether to perform all flags. This will be disabled by '
                 'setting any flag')

    def add_flag(self, *args, **kwargs):
        '''Add a flag to the parser.

        Using this method will set the value of ``--all`` parameter to False.
        '''
        kwargs['action'] = FlagAction
        self.add_argument(*args, **kwargs)

FlagRunner = FlagParser


class FlagAction(argparse.Action):
    '''Action to define a flag.'''
    allArg = 'all'

    def __init__(self,
                 option_strings,
                 dest,
                 default=False,
                 required=False,
                 help=None):
        super(FlagAction, self).__init__(
            option_strings=option_strings,
            dest=dest,
            nargs=0,
            const=True,
            default=default,
            required=required,
            help=help)

    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, self.const)
        print(self.dest)
        if self.dest == self.allArg:
            setattr(namespace, '_force_all', True)
        if not hasattr(namespace, '_force_all'):
            setattr(namespace, self.allArg, False)


def positive_int(arg):
    '''A type checker that accepts only positive integers.'''
    value = int(arg)
    if value <= 0:
        msg = '%s is not a positive integer' % arg
        raise argparse.ArgumentTypeError(msg)
    return value


def nonnegative_int(arg):
    '''A type checker that accepts only non-negative integers.'''
    value = int(arg)
    if value < 0:
        msg = '%s is not a non-negative integer' % arg
        raise argparse.ArgumentTypeError(msg)
    return value


def positive_float(arg):
    '''A type checker that accepts only positive floats.'''
    value = float(arg)
    if value <= 0.:
        msg = '%s is not a positive float' % arg
        raise argparse.ArgumentTypeError(msg)
    return value


def nonnegative_float(arg):
    '''A type checker that accepts only non-negative floats.'''
    value = float(arg)
    if value < 0.:
        msg = '%s is not a non-negative floag' % arg
        raise argparse.ArgumentTypeError(msg)
    return value
