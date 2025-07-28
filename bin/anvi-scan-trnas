#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.terminal as terminal

from anvio.tables.trnahits import TablesForTransferRNAs
from anvio.errors import ConfigError, FilesNPathsError
from anvio.terminal import time_program


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['contigs-db',]
__provides__ = ['hmm-hits',]
__description__ = ("Identify and store tRNA genes in a contigs database")


@time_program
def main():
    args = get_args()
    try:
        tables_for_trna_hits = TablesForTransferRNAs(args)
        tables_for_trna_hits.populate_search_tables(args.contigs_db)
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
    parser.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))
    parser.add_argument(*anvio.A('log-file'), **anvio.K('log-file'))
    parser.add_argument(*anvio.A('trna-hits-file'), **anvio.K('trna-hits-file'))
    parser.add_argument(*anvio.A('trna-cutoff-score'), **anvio.K('trna-cutoff-score'))
    parser.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
