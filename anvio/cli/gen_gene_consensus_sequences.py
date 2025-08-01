#!/usr/bin/env python
# -*- coding: utf-8
"""Collapses variation for a given list of gene caller ids"""

import sys

import anvio
from anvio.argparse import ArgumentParser
import anvio.variabilityops as variabilityops

from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ["profile-db", "contigs-db"]
__provides__ = ["genes-fasta"]
__description__ = "Collapse variability for a set of genes across samples"


def main():
    args = get_args()

    try:
        c = variabilityops.ConsensusSequences(args)
        c.process(exit_if_data_empty=False)
        c.report()
    except ConfigError as e:
        print(e)
        sys.exit(1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(2)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupD = parser.add_argument_group('DATABASES', "Declaring relevant anvi'o databases. First things first.")
    groupD.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db'))
    groupD.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))

    groupS = parser.add_argument_group('FOCUS', 'What do we want? A consensus sequence for a gene, or a list of\
                                                 genes. From where do we want it? All samples, by default. When\
                                                 do we want it? Whenever it is convenient.')
    groupS.add_argument(*anvio.A('gene-caller-ids'), **anvio.K('gene-caller-ids'))
    groupS.add_argument(*anvio.A('genes-of-interest'), **anvio.K('genes-of-interest'))
    groupS.add_argument(*anvio.A('samples-of-interest'), **anvio.K('samples-of-interest'))

    groupO = parser.add_argument_group('OUTPUT', 'Output file and output style')
    groupO.add_argument(*anvio.A('output-file'), **anvio.K('output-file', {'default': 'genes.fa', 'help': 'The output\
                                                  file name. The boring default is "%(default)s". You can change\
                                                  the output file format to a TAB-delimited file using teh flag \
                                                  `--tab-delimited`, in which case please do not forget to change the\
                                                  file name, too.'}))
    groupO.add_argument(*anvio.A('tab-delimited'), **anvio.K('tab-delimited'))


    groupE = parser.add_argument_group('EXTRAS', 'Parameters that will help you to do a very precise analysis.\
                                                  If you declare nothing from this bunch, you will get "everything"\
                                                  to play with, which is not necessarily a good thing...')
    groupE.add_argument(*anvio.A('engine'), **anvio.K('engine'))
    groupE.add_argument(*anvio.A('contigs-mode'), **anvio.K('contigs-mode', {'help': "Use this flag to output consensus\
                                                                                      sequences of contigs, instead of the\
                                                                                      default, which is genes"}))
    groupE.add_argument(*anvio.A('quince-mode'), **anvio.K('quince-mode', {'help': "Use this flag to output consensus \
                                                                                    sequences for cases even where there is\
                                                                                    no variability"}))

    parser.add_argument('--compress-samples', action='store_true', required=False,
                        help="Normally all samples with variation will have their own consensus sequence. If this flag\
                              is provided, the coverages from each sample of interest will be summed and only a single\
                              consenus sequence for each gene/contig will be output.")


    return parser.get_args(parser)


if __name__ == '__main__':
    main()
