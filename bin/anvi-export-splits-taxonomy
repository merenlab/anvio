#!/usr/bin/env python
# -*- coding: utf-8
"""Takes a contigs database, and export taxonomy information from it."""

import sys
from anvio.argparse import ArgumentParser

import anvio
import anvio.dbops as dbops

from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ["contigs-db"]
__provides__ = ["splits-taxonomy-txt"]
__description__ = "Export taxonomy for splits found in an anvi'o contigs database"


def main():
    args = get_args()

    try:
        c = dbops.ContigsSuperclass(args)
        c.gen_TAB_delimited_file_for_split_taxonomies(args.output_file)
    except ConfigError as e:
        print(e)
        sys.exit(1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(2)


def get_args():
    parser = ArgumentParser(description=__description__)

    parser.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))
    parser.add_argument(*anvio.A('output-file'), **anvio.K('output-file', {'required': True}))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
