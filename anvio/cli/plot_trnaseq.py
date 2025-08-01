#!/usr/bin/env python
# -*- coding: utf-8
"""Plot coverage and modification data from flexible groups of tRNA-seq seeds"""

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
__requires__ = ['trnaseq-contigs-db', 'trnaseq-seed-txt', 'modifications-txt']
__provides__ = ['trnaseq-plot']
__description__ = "A program to write plots of coverage and modification data from flexible groups of tRNA-seq seeds"


def main():
    args = get_args()

    try:
        result_plotter = trnaseq.ResultPlotter(args)
        result_plotter.go()
    except ConfigError as e:
        print(e)
        sys.exit(-1)


def get_args():
    parser = ArgumentParser(description=__description__)

    group1 = parser.add_argument_group('MANDATORY')
    group1.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))
    group1.add_argument(*anvio.A('seeds-specific-txt'), **anvio.K('seeds-specific-txt'))
    group1.add_argument(*anvio.A('modifications-txt'), **anvio.K('modifications-txt'))

    group2 = parser.add_argument_group('OPTIONAL')
    group2.add_argument(*anvio.A('output-dir'), **anvio.K('output-dir'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
