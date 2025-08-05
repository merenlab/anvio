#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio

from anvio.errors import ConfigError, FilesNPathsError
from anvio.utils.files import transpose_tab_delimited_file

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__provides__ = ["view-data", "functions-txt", "misc-data-items-txt", "misc-data-layers-txt", "gene-calls-txt", "linkmers-txt"]
__requires__ = ["view-data", "functions-txt", "misc-data-items-txt", "misc-data-layers-txt", "gene-calls-txt", "linkmers-txt"]
__description__ = "Transpose a TAB-delimited file"


def main():
    args = get_args()

    try:
        transpose_tab_delimited_file(args.input_file, args.output_file)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser
    parser = ArgumentParser(description=__description__)

    parser.add_argument('input_file', help="Input matrix.", metavar='MATRIX_FILE')
    parser.add_argument(*anvio.A('output-file'), **anvio.K('output-file', {'required': True, 'metavar':'MATRIX_FILE'}))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
