#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio

from anvio.contigops import GenbankToAnvioWrapper
from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren', 'blankenberg']
__requires__ = []
__provides__ = ["contigs-fasta", "functions-txt", "external-gene-calls"]
__resources__ = [("A tutorial on using this program to access NCBI genomes for 'omics analyses in Anvi'o", "http://merenlab.org/2019/03/14/ncbi-genome-download-magic/")]
__description__ = ("This script takes the 'metadata' output of the program `ncbi-genome-download` (see "
                   "[https://github.com/kblin/ncbi-genome-download](https://github.com/kblin/ncbi-genome-download) for details), and processes each "
                   "GenBank file found in the metadata file to generate a FASTA file, as well as genes "
                   "and functions files for each entry. Plus, it autmatically generates a FASTA TXT "
                   "file descriptor for anvi'o snakemake workflows. So it is a multi-talented program "
                   "like that")


def main():
    args = get_args()

    try:
        genbank_to_anvio = GenbankToAnvioWrapper(args)
        genbank_to_anvio.process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)


    groupA = parser.add_argument_group('INPUT', 'Give us the preciousss...')
    groupA.add_argument('-m', '--metadata', required=True, help='This is the file you get from the\
                            program `ncbi-genome-download` when you use the parameter `--metadata-table`.',
                            metavar='GENBANK_METADATA')

    groupB = parser.add_argument_group('OUTPUT', "Where to find your precioussesss...")
    groupB.add_argument(*anvio.A('output-dir'), **anvio.K('output-dir'))
    groupB.add_argument('--output-fasta-txt', help="This is not a FASTA file, but a TAB-delimited file with all\
                            the file names and paths processed by this program. This output can directly go into\
                            the anvi'o snakemake workflows because magic.", metavar='OUTPUT_FASTA_TXT')

    groupC = parser.add_argument_group('ADDITIONAL PARAMETERS', 'Additional things you can set.')
    groupC.add_argument('-E', '--exclude-gene-calls-from-fasta-txt', required=False, help="This flag will exclude\
                        the external gene calls and functions from the fasta.txt file. Files for external gene calls\
                        and functions according to the information stored in GenBank file, but they will simply\
                        not be included in your fasta.txt file. By doing so you will *gurantee* that when you use\
                        this file from within a workflow, anvi'o wil use its default gene caller to identify genes.", action='store_true')

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
