#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.parsers import parser_modules
from anvio.errors import ConfigError, FilesNPathsError
from anvio.tables.genefunctions import TableForGeneFunctions


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__resources__ = [("Importing functions", "http://merenlab.org/2016/06/18/importing-functions/"),
                 ("Importing GhostKOALA/KEGG annotations", "http://merenlab.org/2018/01/17/importing-ghostkoala-annotations/"),
                 ("Importing VirSorter phage annotaions", "http://merenlab.org/2018/02/08/importing-virsorter-annotations/")]
__requires__ = ['contigs-db', 'functions-txt',]
__provides__ = ['functions',]
__description__ = "Parse and store functional annotation of genes"


@terminal.time_program
def main():
    args = get_args()
    run = terminal.Run()
    progress = terminal.Progress()

    drop_previous_annotations = args.drop_previous_annotations

    try:
        if args.parser:
            functions_dict = get_functions_dict(args)
        else:
            if len(args.input_files) != 1:
                raise ConfigError("When you are using a simple gene calls matrix file to import functions, you can't use more than "
                                   "one file as an input.")

            input_matrix = args.input_files[0]
            functions_dict = utils.get_TAB_delimited_file_as_dictionary(input_matrix,
                                                                        expected_fields = ['gene_callers_id', 'source', 'accession', 'function', 'e_value'],
                                                                        column_mapping  = [int, str, str, str, float],
                                                                        only_expected_fields = True,
                                                                        indexing_field = -1)

        gene_function_calls_table = TableForGeneFunctions(args.contigs_db, run, progress)
        gene_function_calls_table.create(functions_dict, drop_previous_annotations_first = drop_previous_annotations)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_functions_dict(args):
    if args.parser not in parser_modules['functions']:
        raise ConfigError("Anvi'o does not know what to do with '%s'. You must use one of the available parsers "
                           "to make sense of genes (well, open reading frames) found in your contigs (please see "
                           "the documentation for a more detailed explanation): %s"\
                                                 % (args.parser, ', '.join(parser_modules['functions'])))
    if not args.input_files:
        raise ConfigError("You need to use '--input-files' parameter to list file(s) that is/are required the "
                           "parser you chose. Please see the documentation for details.")


    parser = parser_modules['functions'][args.parser](args.input_files)
    genes_dict = parser.get_dict()

    if not len(genes_dict):
        raise ConfigError("Your parser (%s) returned an empty dictionary for your input file. Something must have gone wrong. "
                           "Maybe you selected a wrong parser, or the input file format has changed between when this "
                           "parser was implemented and now. Either ways, if you have exhausted your ideas for troubleshooting "
                           "you should send an e-mail to anvio developers! Sorry for this!" % (args.parser))

    return genes_dict


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    parser.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))

    parser.add_argument('-p', '--parser', default = None,
                        help = "Parser to make sense of the input files (if you need one). There are currently\
                                %d parsers readily available: %s. IT IS OK if you do not select a parser if you\
                                have a standard, TAB-delimited input file for funcitonal annotation of genes. If\
                                this is not like 2018 and everything is already outdated, you should be able to\
                                go to this address and learn everything you need like a boss: http://merenlab.org/2016/06/18/importing-functions/"\
                                            % (len(parser_modules['functions']), list(parser_modules['functions'].keys())))
    parser.add_argument('-i', '--input-files', metavar = 'FILE(S)', nargs='+', default = None, required = True,
                        help = 'One or more input files should follow this parameter. The way these files will be handled\
                                will depend on which parser you selected (if you did select any).')

    parser.add_argument(*anvio.A('drop-previous-annotations'), **anvio.K('drop-previous-annotations'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

