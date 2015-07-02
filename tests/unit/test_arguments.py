'''Test argument parsing.'''
from __future__ import absolute_import, print_function

import sys
import logging

import pytest
from simtools.arguments import (SimulationParser, FlagParser, FlagRunner,
                                positive_int, nonnegative_int, positive_float,
                                nonnegative_float)


@pytest.fixture(params=['-v', '--verbosity'])
def fix_verbosity_argument(request):
    '''Generate long and short verbosity arguments.'''
    return request.param


@pytest.fixture(scope='module', params=[SimulationParser, FlagParser,
                                        FlagRunner])
def fix_parser(request):
    '''Generate classes for all parsers.'''
    return request.param


@pytest.fixture(scope='module', params=[FlagParser, FlagRunner])
def fix_flag_parser(request):
    '''Generate classes for flag parsing.'''
    return request.param


class TestVerbosity(object):
    '''Test verbosity handling.'''

    def _patch_args(self, monkeypatch, parsercls, verbosity,
                    verbosity_argument):
        '''Patch the script arguments and return the parser.'''
        argv = ['pytest_test.py', verbosity_argument, verbosity]
        monkeypatch.setattr(sys, 'argv', argv)
        parser = parsercls()
        o = parser.parse_args()
        return parser, o

    def _test_verbosity(self, monkeypatch, parsercls, wanted_verbosity,
                        verbosity_argument):
        '''Test one verbosity level'''
        _, o = self._patch_args(monkeypatch, parsercls, wanted_verbosity,
                                verbosity_argument)
        assert o.verbosity == wanted_verbosity
        logger = logging.getLogger()
        assert logging.getLevelName(logger.level) == wanted_verbosity

    def test_valid_verbosities(self, monkeypatch, fix_verbosity_argument,
                               fix_parser):
        '''Test all verbosity levels'''

        # Test supported verbosity levels
        verbosities = ['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'ERROR',
                       'CRITICAL']
        for v in verbosities:
            self._test_verbosity(monkeypatch, fix_parser, v,
                                 fix_verbosity_argument)

    def test_invalid_verbosity(self, monkeypatch, fix_verbosity_argument,
                               fix_parser):
        '''Test handling of unsupported verbosity level.'''
        with pytest.raises(SystemExit):
            self._patch_args(monkeypatch, fix_parser, 'INVALID_VERBOSITY',
                             fix_verbosity_argument)


def prepare_parser_with_flags(parsercls=SimulationParser, flagv=None):
    '''Create a parser which contains one flag.'''
    if flagv is None:
        flagv = ['--test_flag']
    parser = parsercls()
    for flag in flagv:
        parser.add_flag(flag)
    return parser


@pytest.fixture(scope='module', params=[(10.0, True),
                                        (10, True),
                                        (-10.0, False),
                                        (-10, False),
                                        (0, False),
                                        (0., False)
                                       ])
def fix_positive_float_data(request):
    '''Generate test data for positive floats.'''
    return request.param


@pytest.fixture(scope='module', params=[(10.0, True),
                                        (10, True),
                                        (-10.0, False),
                                        (-10, False),
                                        (0, True),
                                        (0., True)
                                       ])
def fix_nonnegative_float_data(request):
    '''Generate test data for non-negative floats.'''
    return request.param


class TestSimulationParser(object):
    '''Test aspects of the SimulationParser class.'''

    def test_flag(self, monkeypatch):
        '''Test adding and parsing a flag.'''
        argv = ['test_flag.py', '--test_flag']
        monkeypatch.setattr(sys, 'argv', argv)
        parser = prepare_parser_with_flags()
        o = parser.parse_args()
        assert o.test_flag is True

    def test_no_all(self, monkeypatch):
        '''Test that bare SimulationParser does not contain the --all flag.'''
        argv = ['test_flag.py']
        monkeypatch.setattr(sys, 'argv', argv)
        parser = SimulationParser()
        o = parser.parse_args()
        assert not hasattr(o, 'all')

    def test_positive_int(self, monkeypatch):
        '''Test the positive_int type checker.'''
        argv = ['test_flag.py', '--positive_int', '10']
        monkeypatch.setattr(sys, 'argv', argv)
        parser = SimulationParser()
        parser.add_argument('--positive_int', type=positive_int)
        o = parser.parse_args()
        assert o.positive_int == 10

        argv = ['test_flag.py', '--positive_int', '0']
        monkeypatch.setattr(sys, 'argv', argv)
        parser = SimulationParser()
        parser.add_argument('--positive_int', type=positive_int)
        with pytest.raises(SystemExit):
            o = parser.parse_args()

    def test_nonnegative_int(self, monkeypatch):
        '''Test the nonnegative_int type checker.'''
        argv = ['test_flag.py', '--nonnegative_int', '0']
        monkeypatch.setattr(sys, 'argv', argv)
        parser = SimulationParser()
        parser.add_argument('--nonnegative_int', type=nonnegative_int)
        o = parser.parse_args()
        assert o.nonnegative_int == 0

        argv = ['test_flag.py', '--nonnegative_int', '-1']
        monkeypatch.setattr(sys, 'argv', argv)
        parser = SimulationParser()
        parser.add_argument('--nonnegative_int', type=nonnegative_int)
        with pytest.raises(SystemExit):
            o = parser.parse_args()

    def test_positive_float(self, monkeypatch, fix_positive_float_data):
        '''Test the positive_float type checker.'''
        test_data, is_good = fix_positive_float_data
        argv = ['test_flag.py', '--positive_float', str(test_data)]
        monkeypatch.setattr(sys, 'argv', argv)
        parser = SimulationParser()
        parser.add_argument('--positive_float', type=positive_float)
        if is_good:
            o = parser.parse_args()
            assert o.positive_float == test_data
        else:
            with pytest.raises(SystemExit):
                o = parser.parse_args()

    def test_nonnegative_float(self, monkeypatch, fix_nonnegative_float_data):
        '''Test the nonnegative_float type checker.'''
        test_data, is_good = fix_nonnegative_float_data
        argv = ['test_flag.py', '--nonnegative_float', str(test_data)]
        monkeypatch.setattr(sys, 'argv', argv)
        parser = SimulationParser()
        parser.add_argument('--nonnegative_float', type=nonnegative_float)
        if is_good:
            o = parser.parse_args()
            assert o.nonnegative_float == test_data
        else:
            with pytest.raises(SystemExit):
                o = parser.parse_args()


class TestFlagRunner():
    '''Test aspects of FlagRunner/FlagParser classes.'''
    def test_all_set(self, monkeypatch, fix_flag_parser):
        '''Test whether --all is set.'''
        argv = ['test_flag.py']
        monkeypatch.setattr(sys, 'argv', argv)
        parser = prepare_parser_with_flags(
            parsercls=fix_flag_parser,
            flagv=['--test_flag_1', '--test_flag_2'])
        o = parser.parse_args()
        assert hasattr(o, 'all')
        assert o.all is True

    def test_at_least_one_flag_set(self, monkeypatch, fix_flag_parser):
        '''Test that --all is unset when setting at least one flag.'''
        argv = ['test_flag.py', '--test_flag_1']
        monkeypatch.setattr(sys, 'argv', argv)
        parser = prepare_parser_with_flags(
            parsercls=fix_flag_parser,
            flagv=['--test_flag_1', '--test_flag_2'])
        o = parser.parse_args()
        assert o.all is False

    def test_force_all(self, monkeypatch, fix_flag_parser):
        '''Test forcing --all.'''
        argv = ['test_flag.py', '--all', '--test_flag_1']
        monkeypatch.setattr(sys, 'argv', argv)
        parser = prepare_parser_with_flags(
            parsercls=fix_flag_parser,
            flagv=['--test_flag_1', '--test_flag_2'])
        o = parser.parse_args()
        assert o.all is True
