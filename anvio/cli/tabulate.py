#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio

from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__description__ = ("Tabulates TAB-delmited data with headers in terminal: `cat table.txt | anvi-script-tabulate`")


def main():
    args = get_args()

    try:
        lines = []
        for line in sys.stdin:
            lines.append(line.strip('\n').split('\t'))

        anvio.TABULATE(lines[1:], lines[0], max_width=args.max_width)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser
    parser = ArgumentParser(description=__description__)
    parser.add_argument('--max-width', type=int, default=120, help="Maximum number of characters to be displayed in the output "
                                        "table. The default is %(default)d to make sure tables will fit to most displays. Set "
                                        "to 0 to see the entire table.")
    return parser.get_args(parser)


if __name__ == '__main__':
    main()
