#!/usr/bin/env python
# -*- coding: utf-8
"""A program that computes function enrichment across genomees"""

import sys

import anvio
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError
from anvio.genomedescriptions import AggregateFunctions

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ivagljiva', 'adw96']
__resources__ = [("A description of the enrichment script run by this program can be found in Shaiber et al 2020", "https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02195-w")]
__tags__ = ["functions"]
__requires__ = [ 'groups-txt', 'genomes-storage-db', 'external-genomes', 'internal-genomes', 'functions']
__provides__ = ['functional-enrichment-txt']
__description__ = ("A program that computes functional enrichment across groups of genomes.")


@terminal.time_program
def main():
    args = get_args()
    run = terminal.Run()
    progress = terminal.Progress()

    # make sure the output files is not going to overwrite anything
    if filesnpaths.is_file_exists(args.output_file, dont_raise=True):
        raise ConfigError(f"There is already a file at '{args.output_file}' :/ Anvi'o has a thing against overwriting existing files. "
                          f"Please remove the exiting file first, or give a different output file name to this program.")

    # make sure we can write to the output file
    filesnpaths.is_output_file_writable(args.output_file)

    if not args.groups_txt:
        raise ConfigError("You must have a groups file for this to work. Please see the help menu.")

    # how disturbingly simple is this? :(
    try:
        AggregateFunctions(args, r=run, p=progress)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('SOURCES OF GENOMES', "To estimate enriched functions across your geonomes, you can provide any combination "
                                        "of the following: an external genomes file, an internal genomes file, and/or a genomes storage file. "
                                        "Anvi'o will aggregate all functions in genomes found in any of these sources, and will do its magic.")
    groupA.add_argument(*anvio.A('external-genomes'), **anvio.K('external-genomes'))
    groupA.add_argument(*anvio.A('internal-genomes'), **anvio.K('internal-genomes'))
    groupA.add_argument(*anvio.A('genomes-storage'), **anvio.K('genomes-storage', {'required': False}))

    groupB = parser.add_argument_group('RELATIONSHIPS BETWEEN GENOMES', "Here you need to provide a groups file so anvi'o can understand which "
                                        "genome is in which group.")
    groupB.add_argument(*anvio.A('groups-txt'), **anvio.K('groups-txt', {'required': True}))

    groupC = parser.add_argument_group('FUNCTION ANNOTATION SOURCE', "Here you tell anvi'o which function annotation source to use. Obviously, "
                                        "any given function annotation source will have to be common to all genomes. You can see what you have "
                                        "in any of the contigs databases you can use the program `anvi-db-info`. YOU'RE ALMOST THERE.")
    groupC.add_argument(*anvio.A('annotation-source'), **anvio.K('annotation-source'))
    groupC.add_argument(*anvio.A('list-annotation-sources'), **anvio.K('list-annotation-sources'))

    groupD = parser.add_argument_group('OUTPUT', "What comes out the other end. The only thing mandatory here is the output file.")
    groupD.add_argument(*anvio.A('output-file'), **anvio.K('output-file', {'required': True}))
    groupD.add_argument(*anvio.A('include-ungrouped'), **anvio.K('include-ungrouped'))
    groupD.add_argument(*anvio.A('functional-occurrence-table-output'), **anvio.K('functional-occurrence-table-output'))

    groupE = parser.add_argument_group('OPTIONAL THINGIES', "If you want it, here it is, come and get it.")
    groupE.add_argument(*anvio.A('qlambda'), **anvio.K('qlambda'))
    groupE.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
