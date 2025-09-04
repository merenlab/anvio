#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.terminal as terminal

from anvio.dbops import ContigsDatabase
from anvio.errors import ConfigError, FilesNPathsError
from anvio.tables.genefunctions import TableForGeneFunctions
from anvio.dbinfo import is_contigs_db


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ["contigs-db", "functions"]
__description__ = "Remove functional annotation sources from an anvi'o contigs database"


@terminal.time_program
def main():
    args = get_args()

    A = lambda x: args.__dict__[x] if x in args.__dict__ else None
    contigs_db_path = A('contigs_db')
    annotation_sources = A('annotation_sources')
    list_annotation_sources = A('list_annotation_sources')

    try:
        is_contigs_db(contigs_db_path)

        if list_annotation_sources:
            ContigsDatabase(contigs_db_path).list_function_sources()
            sys.exit()

        if not annotation_sources:
            raise ConfigError("You must provide this program with some functional annotation sources to drop :/")

        sources_to_drop = [s.strip() for s in annotation_sources.split(',')]

        functions_table = TableForGeneFunctions(contigs_db_path)

        contigs_db = ContigsDatabase(contigs_db_path)
        functions_table.drop_functions(database=contigs_db.db, sources_to_drop=sources_to_drop)
        contigs_db.disconnect()
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
    parser.add_argument(*anvio.A('annotation-sources'), **anvio.K('annotation-sources', {'help': "One or more functional annotations sources "
                                        "to drop. If you wish to remove more than one, separate them from each other using a comma character "
                                        "without a space. For example: 'SOURCE_1,SOURCE_2,SOURCE_3'."}))
    parser.add_argument(*anvio.A('list-annotation-sources'), **anvio.K('list-annotation-sources'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

