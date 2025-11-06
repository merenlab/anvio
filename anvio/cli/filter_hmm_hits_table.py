#!/usr/bin/env python
# -*- coding: utf-8

"""
This script removes the hmm_hits table and replaces it with a filtered version. Filtering of hmm_hits
is done using model and/or gene coverage
"""

import sys

import anvio
import anvio.data.hmm
import anvio.terminal as terminal

from anvio.tables.hmmhits import FilterHmmHitsTable
from anvio.errors import ConfigError, FilesNPathsError

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['mschecht']
__provides__ = ["hmm-hits"]
__requires__ = ["contigs-db","hmm-source", "hmm-hits"]
__description__ = ("Filter weak HMM hits from a given contigs database using a domain hits table "
                   "reported by `anvi-run-hmms`.")


@terminal.time_program
def main():
    args = get_args()

    try:
        p = FilterHmmHitsTable(args)
        p.process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser
    parser = ArgumentParser(description=__description__)

    parser.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))
    parser.add_argument(*anvio.A('hmm-source'), **anvio.K('hmm-source'))
    parser.add_argument(*anvio.A('list-hmm-sources'), **anvio.K('list-hmm-sources'))
    parser.add_argument(*anvio.A('hmm-profile-dir'), **anvio.K('hmm-profile-dir'))
    parser.add_argument('--domain-hits-table', metavar='PATH', help="Please provide the path to the domain-table-output. "
                        "You can get this file from running anvi-run-hmms with the flag --domain-hits-table.")
    parser.add_argument('--min-gene-coverage',type=float, help="The minimum percent of the gene that is covered by the profile HMM after hmmsearch. "
                        "This is the formula using the domtblout from hmmsearch: (ali_coord_to - ali_coord_from)/target_length")
    parser.add_argument('--min-model-coverage', type=float, help="The minimum percent of the profile HMM model that is covered by the gene after hmmsearch. "
                        "This is the formula using the domtblout from hmmsearch: (hmm_coord_to - hmm_coord_from)/hmm_length")
    parser.add_argument('--merge-partial-hits-within-X-nts', metavar='LENGTH', default=0, type=int, help="Filtering HMM hits based on "
                        "target or query coverage can be difficult when a gene is covered by multiple independent hits of a single "
                        "model. In these cases, the best way forward is to merge independent hits if they are close enough (for "
                        "instance, if model X covers a gene A between nucleotide positions 0 and 100, and then again between "
                        "nucleotide positions 105 and 200, one would like to merge those two hits into a single one before "
                        "calculating the approximate coverage of the gene by the model). If you set this parameter to any distance, "
                        "anvi'o will first merge independent hits for the same gene/model that are closer to each other than that "
                        "distance.")
    parser.add_argument('--filter-out-partial-gene-calls', action='store_true', help="Partial genes can lead to spurious branches "
                        "and/or inflate the number of observed populations or functions in a given set of genomes/metagenomes. "
                        "Using this flag you can instruct anvi'o to only keep HMM hits from open reading frames that represent complete "
                        "genes (i.e., genes that are not partial and that start with a start codon and end with a stop codon).")

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
