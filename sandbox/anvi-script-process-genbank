#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio

from anvio.contigops import GenbankToAnvio
from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['AstrobioMike']
__provides__ = ["contigs-fasta", "external-gene-calls", "functions-txt"]
__requires__ = ["genbank-file"]
__description__ = ("This script takes a GenBank file, and outputs a FASTA file, as well as two "
                   "additional TAB-delimited output files for external gene calls and gene "
                   "functions that can be used with the programs `anvi-gen-contigs-database` "
                   "and `anvi-import-functions`")


def main():
    args = get_args()

    try:
        genbank_to_anvio = GenbankToAnvio(args)
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
    groupA.add_argument('-i', '--input-genbank', action='store', required=True, help='Input GenBank file', metavar='GENBANK')

    groupB = parser.add_argument_group('OUTPUT', "You either provide a 'prefix', or provide specific output\
                                                  file names/paths. You can't mix the two (well, you can try).")
    groupB.add_argument(*anvio.A('output-file-prefix'), **anvio.K('output-file-prefix'))
    groupB.add_argument("--output-fasta", help='Output FASTA file path.', metavar='FASTA')
    groupB.add_argument("--output-gene-calls", help='Output file path for external gene calls', metavar='TAB DELIMITED FILE')
    groupB.add_argument("--output-functions", help="Output file path for anvi'o-importable gene functions file", metavar='TAB DELIMITED FILE')

    groupC = parser.add_argument_group('DETAILS', "Setting the annotation source and version data to appear in the output\
                                                   file for functional annotations file.")
    groupC.add_argument('--annotation-source', action='store', default='NCBI_PGAP', help='Annotation \
                             source (default: "NCBI_PGAP")')
    groupC.add_argument("--annotation-version", action="store", default="v4.6", help='Annotation \
                             source version to be stored in the database (default: "v4.6")')
    groupC.add_argument("--omit-aa-sequences-column", default=False, action="store_true", help="If amino acid sequences for "
                            "gene calls are present in the GFF file, anvi'o by default will include an `aa_sequences` column "
                            "in the external gene calls output file. This is good because then anvi-gen-contigs-database can "
                            "use that information directly, instead of trying to predict the translated sequences itself. "
                            "If you would like this program to not report amino acid sequences even when they present and use "
                            "anvi'o to predict translated sequences later, you can use this flag.")
    groupC.add_argument("--include-locus-tags-as-functions", default=False, action="store_true", help="When anvi'o processes "
                            "genes in a GenBank file, it will assign new gene caller ids to each entry. Using this flag will "
                            "will instruct anvi'o to add a new 'function' source in the annotation of each gene call so the "
                            "the user can trace back a given gene call to the original locus tag number in the GenBank file.")


    return parser.get_args(parser)


if __name__ == '__main__':
    main()
