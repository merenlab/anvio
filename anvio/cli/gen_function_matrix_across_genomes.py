#!/usr/bin/env python
# -*- coding: utf-8
"""A program that uses AggregateFunctions class to generate human-readable output for functions across genomes."""

import os
import sys

import anvio
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.argparse import ArgumentParser
from anvio.genomedescriptions import AggregateFunctions
from anvio.errors import ConfigError, FilesNPathsError, DictIOError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['functions', 'genomes-storage-db', 'internal-genomes', 'external-genomes', 'groups-txt']
__provides__ = ['functional-enrichment-txt', 'functions-across-genomes-txt']
__description__ = "A program to generate reports for the distribution of functions across genomes"


@terminal.time_program
def main():
    args = get_args()

    try:
        output_directory_for_output_files = os.path.dirname(os.path.abspath(args.output_file_prefix))

        # just a few quick checks before passing everything to AggregateFunctions:
        filesnpaths.is_output_dir_writable(output_directory_for_output_files)

        if args.groups_txt:
            output_file_path_for_functional_enrichment = f"{os.path.abspath(args.output_file_prefix)}-FUNCTIONAL-ENRICHMENT.txt"
            args.output_file = output_file_path_for_functional_enrichment

        facc = AggregateFunctions(args)
        facc.report_functions_across_genomes(args.output_file_prefix, with_function_accession_ids=args.also_report_accession_ids)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)
    except DictIOError as e:
        print(e)
        sys.exit(-3)


def get_args():
    # setup the command line user interface
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('GENOMES', "Tell anvi'o where your genomes are.")
    groupA.add_argument(*anvio.A('internal-genomes'), **anvio.K('internal-genomes'))
    groupA.add_argument(*anvio.A('external-genomes'), **anvio.K('external-genomes'))
    groupA.add_argument(*anvio.A('genomes-storage'), **anvio.K('genomes-storage'))

    groupB = parser.add_argument_group('GROUPS', "If you want, you can also tell anvi'o how to group your genomes "
                                "so it can also compute functional enrichment between them.")
    groupB.add_argument(*anvio.A('groups-txt'), **anvio.K('groups-txt'))
    groupB.add_argument('--print-genome-names-and-quit', default=False, action='store_true', help="Sometimes, "
                                "especially when you are interested in creating a groups file for your genomes you "
                                "gather from multiple different sources, it may be difficult to know every single "
                                "genome name that will go into your analysis. If you declare this flag, after "
                                "initializing everything, anvi'o will print out every genome name it found and quit, "
                                "so you can actually put together a groups file for them.")

    groupC = parser.add_argument_group('FUNCTIONS', "Tell anvi'o which functional annotation source you like above all, and other "
                                "important details you like about your analysis.")
    groupC.add_argument(*anvio.A('annotation-source'), **anvio.K('annotation-source', {'required': True}))
    groupC.add_argument(*anvio.A('aggregate-based-on-accession'), **anvio.K('aggregate-based-on-accession'))
    groupC.add_argument(*anvio.A('aggregate-using-all-hits'), **anvio.K('aggregate-using-all-hits'))
    groupC.add_argument(*anvio.A('also-report-accession-ids'), **anvio.K('also-report-accession-ids'))
    groupC.add_argument('--min-occurrence', metavar="NUM GENOMES", default=1, help=("The minimum number of occurrence of any "
                                "given function accross genomes. If you set a value, those functions that occur in less number "
                                "of genomes will be excluded."), type=int)

    groupE = parser.add_argument_group('GENES', "By default, anvi'o will look for genes in contigs databases that are identified "
                                "by `pyrodigal-gv`. But if you have generated your contigs databse with external gene calls, or have "
                                "otherwise used another gene caller than the default, you can explicitly ask anvi'o to use that "
                                "one to recover your genes.")
    groupE.add_argument(*anvio.A('gene-caller'), **anvio.K('gene-caller'))

    groupF = parser.add_argument_group('ADVANCED', "Kinds of parameters that you should use at your own risk.")
    groupF.add_argument(*anvio.A('skip-checking-genome-hashes'), **anvio.K('skip-checking-genome-hashes'))

    groupG = parser.add_argument_group('OUTPUT', "Provide a 'prefix' for anvi'o to add to the beginning of the output files anvi'o "
                                "will generate for you.")
    groupG.add_argument(*anvio.A('output-file-prefix'), **anvio.K('output-file-prefix', {'required': True}))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
