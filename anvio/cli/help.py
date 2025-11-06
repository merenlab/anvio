#!/usr/bin/env python

import sys

import anvio

from anvio.errors import ConfigError
from anvio.argparse import ArgumentParser
from anvio.programsearch import ProgramSearch

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ekiefl']
__description__ = "Search for anvi'o programs by keyword, inputs/outputs, etc"


def main():
    args = get_args()

    try:
        ProgramSearch(args).process()
    except ConfigError as e:
        print(e)
        sys.exit()


def get_args():
    parser = ArgumentParser(description=__description__)
    parser.add_argument('search-term', help='Find programs associated with this search term (optional)', nargs='?', default='ALL')
    parser.add_argument('--requires', '-r', action='store_true', help='Restrict to programs that require this search term')
    parser.add_argument('--provides', '-p', action='store_true', help='Restrict to programs that provide this search term')
    parser.add_argument('--name', '-n', action='store_true', help='Restrict to programs that contain this search term in their name')
    parser.add_argument('--report', '-R', help='Which information would you like to be in the report? Mess with this if you \
                                                are disappointed with the default. Possibles are Description, Tags, Requires, \
                                                Provides, Status, and Resources. Add multiple of them with commas (no whitespace). \
                                                For example, if you wanted Description and Resources, you would put \
                                                here Description,Resources')

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
