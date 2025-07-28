#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.tables as t
import anvio.terminal as terminal
import anvio.parsers as parsers

from anvio.errors import ConfigError, FilesNPathsError
from anvio.tables.taxonomy import TablesForGeneLevelTaxonomy


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['contigs-db', 'gene-taxonomy-txt',]
__provides__ = ['gene-taxonomy',]
__description__ = "Import gene-level taxonomy into an anvi'o contigs database"


def main():
    try:
        run_program()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def run_program():
    args = get_args()
    run = terminal.Run()
    progress = terminal.Progress()

    if not args.input_files:
        raise ConfigError("You need to use '--input-files' parameter to list file(s) that is/are required the "
                           "parser you chose. Please see the documentation for details.")

    parser_obj = parsers.get_parser_obj('taxonomy_genes', args.parser)

    parser = parser_obj(args.input_files, t.splits_taxonomy_table_structure, run=run, progress=progress)
    parser.just_do_it = args.just_do_it

    genes_taxonomy_dict, taxon_names_dict = parser.process()

    if not len(genes_taxonomy_dict):
        raise ConfigError("Your parser (%s) returned an empty result for your input file. Something must have gone wrong. "
                           "Maybe you selected a wrong parser, or the input file format has changed between when this "
                           "parser was implemented and now. Either ways, if you have exhausted your ideas for troubleshooting "
                           "you should send an e-mail to anvio developers! Sorry for this!" % (args.parser))


    tables_for_taxonomy = TablesForGeneLevelTaxonomy(args.contigs_db, run, progress)
    tables_for_taxonomy.create(genes_taxonomy_dict, taxon_names_dict, args.parser)


def get_args():
    from anvio.argparse import ArgumentParser

    taxonomy_gene_parsers = parsers.get_parser_names('taxonomy_genes')

    parser = ArgumentParser(description=__description__)

    parser.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))

    parser.add_argument('-p', '--parser', default = None,
                        help = 'Parser to make sense of the input files. There are %d parsers readily available: %s.\
                                It is OK if you do not select a parser, but in that case there will be no additional\
                                contigs available except the identification of single-copy genes in your contigs\
                                for later use. Using a parser will not prevent the analysis of single-copy genes,\
                                but make anvio more powerful to help you make sense of your results. Please see the\
                                documentation, or get in touch with the developers if you have any questions\
                                regarding parsers.' % (len(taxonomy_gene_parsers), taxonomy_gene_parsers))
    parser.add_argument('-i', '--input-files', metavar = 'FILE(S)', nargs='+', default = None, required = True,
                        help = 'Input file(s) for selected parser. Each parser (except "blank") requires input files to\
                                process that you generate before running anvio. Please see the documentation for details.')
    parser.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
