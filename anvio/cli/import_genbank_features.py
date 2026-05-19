#!/usr/bin/env python

import sys

import anvio
import anvio.terminal as terminal

from anvio.errors import ConfigError, FilesNPathsError
from anvio.genbankimporter import GenbankFeatureImporter


__copyright__ = "Copyleft 2015-2026, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['semiller10']
__requires__ = ['contigs-db', 'genbank-file']
__provides__ = []
__description__ = ("Imports sequence features from a GenBank file into the "
                   "`contigs_sequence_features` table (and its companions) of "
                   "a contigs database. Strictly additive: the existing "
                   "`genes_in_contigs` table is left untouched.")


@terminal.time_program
def main():
    args = get_args()

    try:
        GenbankFeatureImporter(args).process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    parser.add_argument(*anvio.A('contigs-db'),   **anvio.K('contigs-db'))
    parser.add_argument(*anvio.A('input-genbank'),**anvio.K('input-genbank'))
    parser.add_argument(*anvio.A('source-name'),  **anvio.K('source-name'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
