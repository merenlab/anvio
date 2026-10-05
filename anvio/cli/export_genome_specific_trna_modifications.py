#!/usr/bin/env python
"""Export genome-specific tRNA modification profiles linked to known modification enzymes"""

import sys

import anvio
import anvio.genomictrnaseq as genomictrnaseq
import anvio.terminal as terminal

from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['tucker4']
__resources__ = []
__tags__ = ['trnaseq']
__requires__ = ['trnaseq-contigs-db', 'modifications-txt', 'seeds-specific-txt',
                'external-genomes', 'trna-modification-enzyme-list']
__provides__ = ['genome-specific-trna-modifications-txt']
__description__ = ("Export tRNA-seq modification profiles linked to known modification enzymes "
                   "and the genomes that encode them")


run = terminal.Run()


def main():
    args = get_args()

    try:
        exporter = genomictrnaseq.GenomeSpecificModificationExporter(args)
        exporter.process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)


def get_args():
    parser = ArgumentParser(description=__description__)

    group1 = parser.add_argument_group('TRNASEQ INPUTS')
    group1.add_argument(*anvio.A('trnaseq-contigs-db'), **anvio.K('trnaseq-contigs-db'))
    group1.add_argument(*anvio.A('modifications-txt'), **anvio.K('modifications-txt'))
    group1.add_argument(*anvio.A('seeds-specific-txt'), **anvio.K('seeds-specific-txt'))

    group2 = parser.add_argument_group('GENOMIC INPUTS')
    group2.add_argument(*anvio.A('external-genomes'), **anvio.K('external-genomes'))

    group3 = parser.add_argument_group('ENZYME LIST')
    group3.add_argument(*anvio.A('modification-enzyme-list'), **anvio.K('modification-enzyme-list'))

    group4 = parser.add_argument_group('OUTPUT')
    group4.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))
    group4.add_argument(*anvio.A('enzyme-distribution-output'), **anvio.K('enzyme-distribution-output'))
    group4.add_argument(*anvio.A('overwrite-output-destinations'), **anvio.K('overwrite-output-destinations'))

    group5 = parser.add_argument_group('FILTERS')
    group5.add_argument('--min-coverage-for-detection', metavar='INT', type=int,
                        default=genomictrnaseq.GenomeSpecificModificationExporter.MIN_COVERAGE_FOR_DETECTION_DEFAULT,
                        help="Minimum coverage at a modification position required to report "
                             "modification fractions. Positions with coverage below this threshold "
                             "are reported as NA. Default: %(default)s.")

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
