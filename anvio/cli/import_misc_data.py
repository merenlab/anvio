#!/usr/bin/env python
# -*- coding: utf-8
"""Populate misc data tables"""

import sys
from anvio.argparse import ArgumentParser

import anvio
import anvio.terminal as terminal

from anvio.errors import ConfigError, FilesNPathsError
from anvio.tables.miscdata import MiscDataTableFactory


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren', 'ekiefl']
__resources__ = [("A primer on anvi'o misc data tables", "http://merenlab.org/2017/12/11/additional-data-tables/"),]
__requires__ = ['pan-db', 'profile-db', 'contigs-db', 'misc-data-items-txt', 'dendrogram', 'phylogeny', 'misc-data-layers-txt', 'misc-data-layer-orders-txt', 'misc-data-nucleotides-txt',  'misc-data-amino-acids-txt']
__provides__ = ['misc-data-items', 'misc-data-layers', 'misc-data-layer-orders', 'misc-data-nucleotides', 'misc-data-amino-acids']
__description__ = ("Populate additional data or order tables in pan or profile databases for "
                   "items and layers, OR additional data in contigs databases for nucleotides "
                   "and amino acids (the Swiss army knife-level serious stuff)")


@terminal.time_program
def main():
    # - you're ugly.
    # - no u :(

    args = get_args()

    d = args.__dict__

    try:
        if not d['contigs_db'] and not d['pan_or_profile_db']:
            raise ConfigError("Please provide either a contigs database (--contigs-db) or a profile/pan "
                              "database (--pan-or-profile-db)")

        if d['contigs_db'] and d['pan_or_profile_db']:
            if d['contigs_mode']:
                pass
            else:
                raise ConfigError("You provided a contigs db (--contigs-db) and a profile/pan db (--pan-or-profile-db). "
                                  "Please provide only one.")

        MiscDataTableFactory(args).populate_from_file(args.data_file)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    parser = ArgumentParser(description=__description__)

    parser.add_argument('data_file', metavar = "TAB DELIMITED FILE", help = 'The input file that describes an additional data\
                        for layers or items. The expected format of this file depends on the data table you will target. This\
                        can feel complicated, but we promise it is not (you probably have a PhD or working on one, so trust us\
                        when we say "it is not complicated"). You need to read the online documentation if this is your first\
                        time with this.')


    group1 = parser.add_argument_group("INPUT DATABASES", "The vast majority of cases will require you to provide only "
                        "one of these databases. A special case is when you wish to add `items` data to a profile-db, "
                        "but you only have contig names and not split names for your items. In that case you will want "
                        "to use the `--contigs-mode` flag, and anvi'o will ask you to also provide a contigs database "
                        "path.")
    group1.add_argument(*anvio.A('pan-or-profile-db'), **anvio.K('pan-or-profile-db', {'required': False}))
    group1.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': False}))

    group2 = parser.add_argument_group('TARGET TABLE/GROUP', "Parameters to help you to define where do you want to add "
                        "your data and whether you wish to group it or not (see the online documentation for details).")
    group2.add_argument(*anvio.A('target-data-table'), **anvio.K('target-data-table', {'required': True}))
    group2.add_argument(*anvio.A('target-data-group'), **anvio.K('target-data-group'))

    group3 = parser.add_argument_group('CONVENIENCE', "Parameters of convenience.")
    group3.add_argument(*anvio.A('transpose'), **anvio.K('transpose'))
    group3.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))

    group4 = parser.add_argument_group('CONTIGS MODE', "If you are working with contig names and if you need to  convert "
                        "those contig names to split names, these are the parameters you will need. Please note that this "
                        "is only relevant if you are adding `items` data.")
    group4.add_argument(*anvio.A('contigs-mode'), **anvio.K('contigs-mode'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

