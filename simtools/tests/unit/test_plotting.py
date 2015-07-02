from __future__ import absolute_import, print_function, division

import os.path

import matplotlib
matplotlib.use('agg')   # To prevent problems when DISPLAY is non-existent

import pytest
import matplotlib.pyplot as plt
from configobj import ConfigObj
from simtools.plotting.plotters import (_FigureContextManager, FigurePlotter,
                                        Computation)
from simtools.plotting.env import Environment


class FigurePlotterForTesting(FigurePlotter):
    '''Figure plotter for testing.'''

    def get_fig(self):
        '''Get empty figure.'''
        return plt.figure()

    def get_ax(self, fig):
        '''Get new axes.'''
        return fig.add_subplot(111)


def test_FigureContextManager(tmpdir):
    '''Test the figure context manager.'''
    subdir = tmpdir.mkdir('sub')
    fname = str(subdir.join("test_figure.pdf"))
    plotter = FigurePlotterForTesting(config={}, env=None)
    with _FigureContextManager(fname, plotter, True) as (fig, ax):
        ax.plot([1], [1])
    lst = subdir.listdir()
    assert len(lst) == 1
    assert lst[0].basename == 'test_figure.pdf'


class TestComputation(object):
    '''Test the Computation class.'''

    def test_init(self):
        '''Test basic object creation.'''
        config = ConfigObj()
        env = None
        comp = Computation(config, env)
        assert comp.config is config
        assert comp.env is None

    def test_class_config(self):
        '''Test class-specific configuration.'''
        # Config section is not present
        config = ConfigObj()
        comp = Computation(config, None)
        assert comp.myc == {}

        # Config section contains 1 item
        config = ConfigObj({
            'Computation' : {
                'test': 10,
            }
        })
        comp = Computation(config, None)
        assert len(comp.myc) == 1
        assert comp.myc['test'] == 10

    def test_abstract_methods(self):
        '''Test abstract methods.'''
        config = ConfigObj()
        comp = Computation(config, None)
        with pytest.raises(NotImplementedError):
            comp.run_all()

    def test_get_fname(self):
        '''Test file name generation.'''
        # TODO: test keyword arguments
        output_dir = 'test_output_dir'
        fname_prefix = 'global_fname_prefix'
        obj_fname_prefix = 'obj_fname_prefix'
        filename = 'test_file.pdf'
        config = ConfigObj({
            'output_dir' : output_dir,
            'fname_prefix' : fname_prefix,
            'Computation' : {
                'fname_prefix' : obj_fname_prefix
            },
        })
        correct_path = os.path.join(
            output_dir,
            '{fname_prefix}{obj_prefix}{file_base}'.format(
                fname_prefix=fname_prefix,
                obj_prefix=obj_fname_prefix,
                file_base=filename
            )
        )
        comp = Computation(config, None)
        assert comp.get_fname(filename) == correct_path


class TestEnvironment(object):
    '''Test the Environment class.'''

    class NumPrinter(Computation):
        '''Simply prints a number to stdout.'''
        def run_all(self, *args, **kwargs):
            print(self.myc['number'])

    class TestArguments(Computation):
        '''Require a single positional argument before *args, and **kwargs.'''
        def __init__(self, my_param, *args, **kwargs):
            self.my_param = my_param
            self.kwarg = kwargs.pop('kwarg')
            super(TestEnvironment.TestArguments, self).__init__(*args, **kwargs)

        def run_all(self, *args, **kwargs):
            print(self.my_param)
            print(self.kwarg)

    def test_config(self):
        '''Test config correctness.'''
        config = ConfigObj({'test' : 10})
        env = Environment(config)
        assert env.config == config

    def test_running(self, capsys):
        '''Test runs with several classes.'''
        config = ConfigObj({
            'NumPrinter' : {
                'number' : 1,
            },
        })

        # No object-specific configurations
        env = Environment(config)
        for i in range(3):
            env.register_class(self.NumPrinter)
        env.run()
        out, _ = capsys.readouterr()
        assert out == "1\n1\n1\n"

        # Each object except the first has a specific configuration
        env = Environment(config)
        env.register_plotter(self.NumPrinter)
        for i in range(3):
            env.register_class(self.NumPrinter, config={
                'NumPrinter': {
                    'number' : i + 2,
                },
            })
        env.run()
        out, _ = capsys.readouterr()
        assert out == "1\n2\n3\n4\n"

    def test_extra_arguments(self, capsys):
        '''Test that extra arguments to Computation subclasses are handled
        correctly.
        '''
        env = Environment(ConfigObj())
        my_param = 3
        env.register_class(self.TestArguments, None, True,
                           my_param, kwarg=10)
        env.run()
        out, _ = capsys.readouterr()
        assert out == "3\n10\n"

