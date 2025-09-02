#!/usr/bin/env python
# -*- coding: utf-8
"""Generator of genome storages"""

import sys

import anvio
import anvio.genomestorage as genomestorage
import anvio.genomedescriptions as genomedescriptions

from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['external-genomes', 'internal-genomes']
__provides__ = ['genomes-storage-db']
__description__ = "Create a genome storage from internal and/or external genomes for a pangenome analysis"
__resources__ = [("A tutorial on pangenomics", "http://merenlab.org/2016/11/08/pangenomics-v2/"),]


def main():
    args = get_args()

    try:
        genome_descriptions = genomedescriptions.GenomeDescriptions(args)
        genome_descriptions.load_genomes_descriptions()

        genome_storage = genomestorage.GenomeStorage(args.output_file, create_new=True)
        genome_storage.store_genomes(genome_descriptions)

    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('EXTERNAL GENOMES', "External genomes listed as anvi'o contigs databases. As in, you have one\
                                                    or more genomes say from NCBI you want to work with, and you created an\
                                                    anvi'o contigs database for each one of them.")
    groupA.add_argument(*anvio.A('external-genomes'), **anvio.K('external-genomes'))

    groupB = parser.add_argument_group("INTERNAL GENOMES", "Genome bins stored in an anvi'o profile databases as collections.")
    groupB.add_argument(*anvio.A('internal-genomes'), **anvio.K('internal-genomes'))

    groupC = parser.add_argument_group("PRO STUFF", "Things you may not have to change. But you never know (unless you read the help).")
    groupC.add_argument(*anvio.A('gene-caller'), **anvio.K('gene-caller'))

    groupD = parser.add_argument_group("OUTPUT", "Give it a nice name. Must end with '-GENOMES.db'. This is primarily due to the fact\
                                                  that there are other .db files used throughout anvi'o and it would be better to\
                                                  distinguish this very special file from them.")
    groupD.add_argument(*anvio.A('output-file'), **anvio.K('output-file', {'required': True, 'metavar': 'GENOMES_STORAGE'}))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
