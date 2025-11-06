#!/usr/bin/env python
# -*- coding: utf-8
"""Tabulate data on tRNA seeds from a set of tRNA-seq samples"""


import sys

import anvio
import anvio.trnaseq as trnaseq

from anvio.errors import ConfigError
from anvio.argparse import ArgumentParser


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"
__requires__ = ['trnaseq-contigs-db', 'trnaseq-profile-db']
__provides__ = ['trnaseq-seed-txt', 'modifications-txt']
__description__ = "A program to write standardized tab-delimited files of tRNA-seq seed coverage and modification results"


def main():
    args = get_args()

    try:
        result_tabulator = trnaseq.ResultTabulator(args)
        result_tabulator.process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)


def get_args():
    parser = ArgumentParser(description=__description__)

    group1 = parser.add_argument_group('MANDATORY')
    group1.add_argument(*anvio.A('contigs-db'),
                        **anvio.K('contigs-db',
                                  {'help': "Path to a `trnaseq`-variant contigs database, as produced by `anvi-merge-trnaseq`."}))
    group1.add_argument(*anvio.A('specific-profile-db'), **anvio.K('specific-profile-db'))

    group2 = parser.add_argument_group('OPTIONAL')
    group2.add_argument(*anvio.A('nonspecific-profile-db'), **anvio.K('nonspecific-profile-db'))
    group2.add_argument(*anvio.A('output-dir'), **anvio.K('output-dir'))
    group2.add_argument(*anvio.A('overwrite-output-destinations'), **anvio.K('overwrite-output-destinations'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
