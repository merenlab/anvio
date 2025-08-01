#!/usr/bin/env python
# -*- coding: utf-8
"""anvi-analyze-synteny ... let's extract some order from this chaos!

What do we say to loci that appear to have no coherent synteny patterns...?

... not today ;)

The goal of this program is to analyze the synteny patterns within a group of loci. The main idea is to have a sliding window
that will record a sequence of genes of size N to form an Ngram (similar to a kmer profile). This program will output
a tsv with counts for each Ngram per contig provided. The resulting Ngram count table can then be used
to analyze the presence and absence of groups of genes by the user.

See anvio/docs/programs/anvi-analyze-synteny.md for more information and how to create simple test cases.
"""

import sys

import anvio

from anvio.synteny import NGram
from anvio.terminal import time_program
from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['mschecht']
__requires__ = ['genomes-storage-db', 'functions', 'pan-db']
__provides__ = ['ngrams']
__description__ = ("Extract ngrams, as in 'co-occurring genes in synteny', from genomes")


@time_program
def main():
    args = get_args()

    try:
        ngram = NGram(args)
        ngram.report_ngrams_to_user()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def get_args():
    parser = ArgumentParser(description=__description__)

    # Essential input
    groupA = parser.add_argument_group('Essential INPUT')
    # FIXME: add optional input where user can just provide list of contigsDBs

    groupA.add_argument(*anvio.A('genomes-storage'), **anvio.K('genomes-storage', {'required': True}))
    groupA.add_argument(*anvio.A('ngram-window-range'), **anvio.K('ngram-window-range', {"required": False}))
    groupA.add_argument(*anvio.A('output-file'), **anvio.K('output-file', {"required": False}))

    # Annotation sources
    groupB = parser.add_argument_group('Annotation sources for Ngrams', "Choose one source of annotations for your Ngrams.")
    groupB.add_argument(*anvio.A('annotation-source'), **anvio.K('annotation-source', {"required": False}))
    groupB.add_argument(*anvio.A('pan-db'), **anvio.K('pan-db', {"required": False}))
    groupB.add_argument('-n', '--ngram-source', help='If two annotation sources are provided, please choose one annotation source'
                                                     ' that will be used to calcuate Ngrams (e.g. gene_clusters, COG_FUNCTION)')

    # Optional arguments
    groupC = parser.add_argument_group('Optional arguments')
    groupC.add_argument(*anvio.A('list-annotation-sources'), **anvio.K('list-annotation-sources', {"required": False}))
    groupC.add_argument('--analyze-unknown-functions',action="store_true",
                        help="Provide this flag if you want anvi-analyze-synteny to report "
                              'Ngrams that contain gene calls that have no annotation.')
    groupC.add_argument(*anvio.A('genomes-names'), **anvio.K('genomes-names', {"required": False}))
    groupC.add_argument('--first-functional-hit-only',action="store_true",
                        help="Use this flag if you want to use on the first functional annotation when making ngrams and assigning annotations. "
                             "In some cases, anvio reports more than one annotation when there are multiple good hits to the gene. "
                             "When this happens, all annotations will be reported in order of alignment score and delimited by '!!!' e.g. 'COG123!!!COG456!!!COG789'. "
                             "This flag will report 'COG123!!!COG456!!!COG789' as 'COG123'.")

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
