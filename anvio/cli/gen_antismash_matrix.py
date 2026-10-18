#!/usr/bin/env python

import sys

import anvio

from anvio.antismash import AntismashMatrix
from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError
from anvio.terminal import time_program

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['jessika-fuessel']
__requires__ = ['functions', 'external-genomes']
__can_use__ = ['internal-genomes']
__provides__ = ['functions-across-genomes-txt']
__description__ = ("With this script you can generate a matrix across many contig-dbs or bins in a collection to "
                   "compare their biosynthetic capacity. You can choose between a per-gene matrix that resolves "
                   "the smCOG annotations and domains within each gene and their role within the BGC and a "
                   "cluster resolved matrix that summarizes the individual genes into BGC and their product. To "
                   "generate these matrices your contig-dbs already need to be annotated with antiSMASH. Use the "
                   "--genes flag to get the 'per-gene' output or the --regions flag to get the 'BGC-resolved' "
                   "output.")


@time_program
def main():
    args = get_args()

    try:
        AntismashMatrix(args).process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group("INPUT", "This script accepts external- and internal-genomes files as input, "
                             "simply use either the -e or -i flag (or both).")
    groupA.add_argument(*anvio.A('external-genomes'), **anvio.K('external-genomes', {'required': False}))
    groupA.add_argument(*anvio.A('internal-genomes'), **anvio.K('internal-genomes', {'required': False}))

    groupB = parser.add_argument_group("MATRIX TYPE", "Choose exactly one of the following.")
    groupB.add_argument(*anvio.A('antismash-genes-matrix'), **anvio.K('antismash-genes-matrix'))
    groupB.add_argument(*anvio.A('antismash-regions-matrix'), **anvio.K('antismash-regions-matrix'))

    groupC = parser.add_argument_group("OUTPUT", "A prefix for the two matrix files (frequency and "
                             "presence/absence) this program will write.")
    groupC.add_argument(*anvio.A('output-file-prefix'), **anvio.K('output-file-prefix', {'required': True}))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
