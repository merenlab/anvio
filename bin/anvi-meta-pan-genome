#!/usr/bin/env python
# -*- coding: utf-8
"""Entering anvi'o metapangenomic workflow"""

import sys
from anvio.argparse import ArgumentParser

import anvio
import anvio.metapanops as metapanops
import anvio.terminal as terminal

from anvio.errors import ConfigError, FilesNPathsError, HDF5Error


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['internal-genomes', 'pan-db', 'genomes-storage-db',]
__provides__ = ['metapangenome',]
__description__ = "Convert a pangenome into a metapangenome"


def main():
    args = get_args()
    run = terminal.Run()
    progress = terminal.Progress()

    try:
        pan = metapanops.MetaPangenome(args, run, progress)
        pan.process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)
    except HDF5Error as e:
        print(e)
        sys.exit(-2)


def get_args():
    parser = ArgumentParser(description=__description__)
    groupA = parser.add_argument_group("PANGENOME", "Files for the pangenome.")
    groupA.add_argument(*anvio.A('pan-db'), **anvio.K('pan-db'))
    groupA.add_argument(*anvio.A('genomes-storage'), **anvio.K('genomes-storage'))

    groupB = parser.add_argument_group("METAGENOME", "Genome bins stored in an anvi'o profile databases as collections.")
    groupB.add_argument('-i', '--internal-genomes', metavar = 'FILE', default = None,
                        help = "A four-column TAB-delimited flat text file. This file should be identical to the internal\
                                genomes file you used for your pangenomics analysis. Anvi'o will use this file to find all\
                                profile and contigs databases that contain the information for each gene and genome across\
                                metagenomes.")

    groupC = parser.add_argument_group("CRITERION FOR DETECTION", "This is tricky. What we want to do is to identify genes that are\
                                        occurring uniformly across samples.")
    groupC.add_argument('--fraction-of-median-coverage', metavar="FLOAT", default=0.25, type=float, help="The value set here\
                        will be used to remove a gene if its total coverage across environments is less than the median coverage\
                        of all genes multiplied by this value. The default is 0.25, which means, if the median total coverage of\
                        all genes across all samples is 100X, then, a gene with a total coverage of less than 25X across all\
                        samples will be assumed not a part of the 'environmental core'.")
    groupC.add_argument('--min-detection', metavar="FLOAT", default=0.50, type=float, help="For this entire thing to work, the\
                        genome you are focusing on should be detected in at least one metagenome. If that is not the case, it would\
                        mean that you do not have any sample that represents the niche for this organism (or you do not have enough\
                        depth of coverage) to investigate the detection of genes in the environment. By default, this script requires\
                        at least '0.5' of the genome to be detected in at least one metagenome. This parameter allows you to change\
                        that. 0 would mean no detection test required, 1 would mean the entire genome must be detected.")

    groupD = parser.add_argument_group("PRO STUFF", "Things you may not have to change unless you are doing something extra cool and edgy.")
    groupD.add_argument(*anvio.A('gene-caller'), **anvio.K('gene-caller'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
