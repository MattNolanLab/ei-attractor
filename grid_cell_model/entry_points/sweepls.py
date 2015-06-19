'''List important information about parameter sweeps.'''
from __future__ import absolute_import, print_function, division
from enum import Enum, IntEnum, unique
import os.path
import collections

from simtools.storage import DataStorage
from grid_cell_model.submitting.base.parsers import BaseParser


_DESCRIPTION = ("List various forms of data from parameter sweeps of grid cell "
                "networks.")
LIST_JUST_WIDTH = 6

@unique
class Errnum(IntEnum):
    '''Error numbers.'''
    success = 0
    fail = 1


@unique
class Commands(str, Enum):
    '''Command line commands to execute.'''
    print_data = "print-data"

    @classmethod
    def get_choice_list(cls):
        '''Get a list of the values of the commands.'''
        result = []
        for command in cls:
            result.append(command.value)
        return result

    @classmethod
    def validate(cls, name):
        '''Validate whether the ``name`` is a valid string representing a
        command.'''
        for command in cls:
            if name == command.value:
                return command

        raise ValueError("Command '%s' is not a valid command. Commands can "
                         "only be one of %s", name, cls.get_choice_list())


def split_data_path(path):
    '''Split data path components separated by / into a list.'''
    split_path = path.split('/')
    result = []
    for part in split_path:
        if part != '':
            result.append(part)
    return result


def print_data_type(data):
    '''Print a data structure ``data``, depending on its type.'''
    if isinstance(data, collections.Sequence):
        for idx, item in enumerate(data):
            print("%s%s" % (str(idx).ljust(LIST_JUST_WIDTH), str(item)))
    elif isinstance(data, collections.Mapping):
        max_key_len = 0
        for key in data.keys():
            if len(key) > max_key_len:
                max_key_len = len(key)

        for key, value in data.items():
            print("%s%s" % (key.ljust(max_key_len + 2), str(value)))
    else:
        print(data)

def print_data(args):
    '''Print data from a file or directory containing the parameter sweep.

    Parameters
    ----------
    args : argparse.Namespace
        Arguments obtained from ``argparse.Argumentparser``.

    Returns
    -------
    errno : int
        Error number.
    '''
    data = DataStorage.open(args.path, 'r')
    dataset_path = split_data_path(args.data)
    if len(dataset_path) == 0:
        print_data_type(data)
    else:
        print_data_type(data.get_item_chained(dataset_path))

    return Errnum.success


def perform_command(args):
    '''Run the specified command.'''
    command = Commands.validate(args.command)

    if command is Commands.print_data:
        return print_data(args)
    else:
        print("The command '%s' is not yet supported" % command.value)
        return Errnum.fail


def main():
    parser = BaseParser(description=_DESCRIPTION)
    parser.add_argument('command', type=str,
                        choices=Commands.get_choice_list(),
                        help="Command to execute.")
    parser.add_argument('path', type=str,
                        help='Path to the structure to list data from. '
                             'Depending on the command this can either be a '
                             'file or a directory.')
    parser.add_argument('-d', '--data', type=str,
                        help="Path to the data with a file or parameter "
                             "sweep.")
    args = parser.parse_args()
    return perform_command(args)


if __name__ == "__main__":
    exit(main())
