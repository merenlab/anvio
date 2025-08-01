#!/usr/bin/env python
# -*- coding: utf-8

import sys

from argparse import Namespace

import anvio
import anvio.terminal as terminal

from anvio.cogs import COGs
from anvio.terminal import time_program
from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren', 'ge0rges']
__requires__ = ['cogs-data', 'contigs-db', 'fasta']
__provides__ = ['functions', 'functions-txt']
__description__ = ("This program runs NCBI's COGs to associate genes in an anvi'o contigs database with functions. "
                   "This program can also run NCBI's COGs to annotate an amino acid sequence with function. "
                   "COGs database was been designed as an attempt to classify proteins from completely "
                   "sequenced genomes on the basis of the orthology concept.")


@time_program
def main():
    args = get_args()

    try:
        cogs = COGs(args)
        aa_file = args.fasta_file or None
        cogs.process(aa_sequences_file_path=aa_file)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def get_args():
    try:
        cogs = COGs(Namespace())
        available_search_methods = cogs.available_search_methods
        default_search_method = cogs.default_search_method
    except ConfigError as e:
        print(e)
        sys.exit(-1)

    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('INPUT #1 - ANNOTATION OF CONTIGS DB', "One type of input this function can annotate is a contigs database. "
                                       "In this case anvi'o will attempt to annotate every gene present.")
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {"required": False}))

    groupB = parser.add_argument_group('INPUT #2 - ANNOTATION OF FASTA FILE', "Another type of input this function can annotate is a FASTA file "
                                       "containing one or more amino acid sequences. In this case each gene sequence will be annotated "
                                       "and the functions will be reported as TAB-delimited output. In this mode you should specify an output file.")
    groupB.add_argument(*anvio.A('fasta-file'), **anvio.K('fasta-file', {"help": "A FASTA formatted input file containing protein sequences."}))
    groupB.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))

    groupC = parser.add_argument_group('GENERAL OPTIONS', "These are general options which apply to both input modes.")
    groupC.add_argument(*anvio.A('cog-version'), **anvio.K('cog-version'))
    groupC.add_argument(*anvio.A('cog-data-dir'), **anvio.K('cog-data-dir'))
    groupC.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))
    groupC.add_argument(*anvio.A('temporary-dir-path'), **anvio.K('temporary-dir-path'))

    groupC.add_argument('--search-with', default=default_search_method, metavar="SEARCH_METHOD",
                        help="What program to use for database searching. The default search uses %(default)s.\
                              All available options include: %(serach_methods)s." \
                                        % {'serach_methods': ', '.join(available_search_methods),
                                           'default': default_search_method})

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
