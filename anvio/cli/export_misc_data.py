#!/usr/bin/env python
# -*- coding: utf-8
"""Export misc data from target tables into TAB-delimited files"""

import sys
from anvio.argparse import ArgumentParser

import anvio

from anvio.errors import ConfigError, FilesNPathsError
from anvio.tables.miscdata import MiscDataTableFactory


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren', 'ekiefl']
__requires__ = ['pan-db', 'profile-db', 'contigs-db', 'misc-data-items', 'misc-data-layers', 'misc-data-layer-orders', 'misc-data-nucleotides', 'misc-data-amino-acids']
__provides__ = ['misc-data-items-txt', 'misc-data-layers-txt', 'misc-data-layer-orders-txt', 'misc-data-nucleotides-txt', 'misc-data-amino-acids-txt']
__description__ = "Export additional data or order tables in pan or profile databases for items or layers"
__resources__ = [("Working with anvi'o additional data tables", "http://merenlab.org/2017/12/11/additional-data-tables/#views-items-layers-orders-some-anvio-terminology")]


def main():
    args = get_args()

    try:
        d = args.__dict__
        if not d['contigs_db'] and not d['pan_or_profile_db']:
            raise ConfigError("Please provide either a contigs database (--contigs-db) or a profile/pan "
                              "database (--pan-or-profile-db)")

        if d['contigs_db'] and d['pan_or_profile_db']:
            raise ConfigError("You provided a contigs db (--contigs-db) and a profile/pan "
                              "db (--pan-or-profile-db). Please provide only one.")

        table_for_additional_data = MiscDataTableFactory(args)
        table_for_additional_data.export(output_file_path=args.output_file)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    parser = ArgumentParser(description=__description__)

    group1 = parser.add_argument_group('Database input', "Provide 1 of these")
    group2 = parser.add_argument_group('Details', "Everything else.")

    group1.add_argument(*anvio.A('pan-or-profile-db'), **anvio.K('pan-or-profile-db', {'required': False}))
    group1.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': False}))

    group2.add_argument(*anvio.A('target-data-table'), **anvio.K('target-data-table', {'required': True}))
    group2.add_argument(*anvio.A('target-data-group'), **anvio.K('target-data-group'))
    group2.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
