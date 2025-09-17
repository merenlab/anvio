#!/usr/bin/env python
# coding: utf-8

import sys

import anvio
import anvio.terminal as terminal
import anvio.pangenome_graph_preprocess as preprocess

from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ahenoch', 'meren']
__requires__ = []
__provides__ = []
__description__ = "Reorient contigs based on the longest complete genome presend in the given folder"


@terminal.time_program
def main():
    args = get_args()

    try:
        prep = preprocess.external_genomes_preprocess(args)
        prep.process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)

def get_args():
    parser = ArgumentParser(description=__description__)
    groupA = parser.add_argument_group('INPUT FILES')
    groupA.add_argument(*anvio.A('fasta-text-file'), **anvio.K('fasta-text-file', {'required': True}))
    groupA.add_argument(*anvio.A('output-dir'), **anvio.K('output-dir', {'help': ""}))

    groupB = parser.add_argument_group('REFERECE GENOME OPTIONS', description="Orienting contigs will require a reference genome that "
                "must be chosen by none other than YOU.")
    groupB.add_argument('--prioritize-genome-size', action="store_true", default=False, help="Of all the candidates, choose the genome "
                        "with largest size as 'reference' to orient contigs in other FASTA files.")
    groupB.add_argument('--prioritize-number-of-contigs', action="store_true", default=False, help="Of all the candidates, choose the genome "
                        "with the smallest number of contigs as 'reference' to orient contigs in other FASTA files.")

    return parser.get_args(parser)


if __name__ == '__main__':
    main()