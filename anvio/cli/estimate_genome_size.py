#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError, FilesNPathsError

SOME_GENOME_SIZE = 42

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__description__ = "A program to estimate the size of the actual population genome to which a MAG belongs"


def main():
    args = get_args()
    run = terminal.Run()

    try:
        utils.is_contigs_db(args.contigs_db)

        if args.verbose:
            run.info('Estimated genome size', '%d (Forty two)' % 42)
        else:
            run.info('Estimated genome size', 42)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('MANDATORY INPUT', "An anvi'o contigs database that hopefully contains a MAG.")
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))

    groupD = parser.add_argument_group('PARAMETERS OF CONVENIENCE', "Because life is already very hard as it is.")
    groupD.add_argument(*anvio.A('verbose'), **anvio.K('verbose'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
