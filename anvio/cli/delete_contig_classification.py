#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError, FilesNPathsError
from anvio.tables.contigclassification import TablesForContigClassification

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = []
__requires__ = ['contigs-db', 'contig-classification']
__provides__ = []
__description__ = "Remove contig classification data from an anvi'o contigs database"


@terminal.time_program
def main():
    args = get_args()

    A = lambda x: args.__dict__[x] if x in args.__dict__ else None
    contigs_db_path = A('contigs_db')
    source = A('source')
    just_do_it = A('just_do_it')
    list_sources = A('list_contig_classification_sources')

    try:
        utils.is_contigs_db(contigs_db_path)
        tables = TablesForContigClassification(contigs_db_path)

        if list_sources:
            tables.list_sources()
            sys.exit()

        tables.delete(source=source, just_do_it=just_do_it)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser
    parser = ArgumentParser(description=__description__)

    group1 = parser.add_argument_group('REQUIRED')
    group1.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))

    group2 = parser.add_argument_group('OPTIONAL')
    group2.add_argument('--source', metavar='SOURCE', default=None,
                        help="The classification source to delete. If not provided, ALL classification "
                             "data will be deleted from the database.")
    group2.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))
    group2.add_argument(*anvio.A('list-contig-classification-sources'),
                        **anvio.K('list-contig-classification-sources'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
