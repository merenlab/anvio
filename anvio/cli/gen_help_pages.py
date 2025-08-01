#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.programs as programs
import anvio.terminal as terminal

from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren', 'Jessica-Pan']
__description__ = "Generate a static web site for anvi'o help pages"


@terminal.time_program
def main():
    args = get_args()

    try:
        docs = programs.AnvioDocs(args)
        docs.generate()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    parser = ArgumentParser(description=__description__)

    parser.add_argument(*anvio.A('output-dir'), **anvio.K('output-dir'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
