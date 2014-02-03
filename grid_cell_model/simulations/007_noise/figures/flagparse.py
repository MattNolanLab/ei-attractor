'''
.. currentmodule:: flagparser

Figure plotting flag parser
'''
import argparse


class FlagParser(argparse.ArgumentParser):

    def __init__(self, allDefault=True):
        argparse.ArgumentParser.__init__(self)

        self.allDefault = allDefault
        self.add_argument('--' + FlagAction.allArg, action='store_true',
                default=allDefault, help='Whether to perform all flags. This'+\
                        ' will be disabled by setting any flag')


    def add_flag(self, *args, **kwargs):
        kwargs['action'] = FlagAction
        self.add_argument(*args, **kwargs)



class FlagAction(argparse.Action):
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
        setattr(namespace, self.allArg, False)
        setattr(namespace, self.dest, self.const)
