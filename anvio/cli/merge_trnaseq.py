#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.trnaseq as trnaseq
import anvio.terminal as terminal

from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['semiller10']
__requires__ = ['trnaseq-db']
__provides__ = ['trnaseq-contigs-db', 'trnaseq-profile-db']
__description__ = ("This program processes one or more anvi'o tRNA-seq databases produced by `anvi-trnaseq` "
                   "and outputs anvi'o contigs and merged profile databases accessible to other tools in the anvi'o ecosystem. "
                   "Final tRNA \"seed sequences\" are determined from a set of samples. "
                   "Each sample yields a set of tRNA predictions stored in a tRNA-seq database, "
                   "and these tRNAs may be shared among the samples. "
                   "tRNA may be 3' fragments and thereby subsequences of longer tRNAs from other samples which would become seeds. "
                   "The profile database produced by this program records the coverages of seeds in each sample. "
                   "This program finalizes predicted nucleotide modification sites using tunable substitution rate parameters")


@terminal.time_program
def main():
    args = get_args()

    try:
        converter = trnaseq.DatabaseMerger(args)
        converter.process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('MANDATORY')
    groupA.add_argument('input', metavar='TRNASEQ_DB(S)', nargs='+',
                        help="Anvi'o tRNA-seq databases representing samples in an experiment")
    groupA.add_argument(*anvio.A('output-dir'), **anvio.K('output-dir'))
    groupA.add_argument(*anvio.A('project-name'), **anvio.K('project-name'))

    groupB = parser.add_argument_group('EXTRAS')
    groupB.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))
    groupB.add_argument(*anvio.A('max-reported-trna-seeds'), **anvio.K('max-reported-trna-seeds'))
    groupB.add_argument(*anvio.A('overwrite-output-destinations'), **anvio.K('overwrite-output-destinations'))
    groupB.add_argument(*anvio.A('description'), **anvio.K('description'))

    groupC = parser.add_argument_group('ADVANCED')
    groupC.add_argument(*anvio.A('feature-threshold'), **anvio.K('feature-threshold'))
    groupC.add_argument(*anvio.A('preferred-treatment'), **anvio.K('preferred-treatment'))
    groupC.add_argument(*anvio.A('nonspecific-output'), **anvio.K('nonspecific-output'))
    groupC.add_argument(*anvio.A('min-variation'), **anvio.K('min-variation'))
    groupC.add_argument(*anvio.A('min-third-fourth-nt'), **anvio.K('min-third-fourth-nt'))
    groupC.add_argument(*anvio.A('min-indel-fraction'), **anvio.K('min-indel-fraction'))
    groupC.add_argument(*anvio.A('distance'), **anvio.K('distance'))
    groupC.add_argument(*anvio.A('linkage'), **anvio.K('linkage'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
