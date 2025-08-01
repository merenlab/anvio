#!/usr/bin/env python
# -*- coding: utf-8
"""Gives you back items in a collection.

   Output files can directly be used by anvi-import-collection"""

import sys
from anvio.argparse import ArgumentParser

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.ccollections as ccollections

from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['profile-db', 'collection']
__provides__ = ['collection-txt']
__description__ = "Export a collection from an anvi'o database"


def main():
    args = get_args()
    run = terminal.Run()
    progress = terminal.Progress()

    try:
        utils.is_pan_or_profile_db(args.pan_or_profile_db)

        c = ccollections.Collections(r=run, p=progress)
        c.populate_collections_dict(args.pan_or_profile_db)

        if args.list_collections:
            c.list_collections()
            sys.exit()

        if not args.collection_name:
            raise ConfigError("You must declare a colection name. Use --list-collections flag "
                               "to list available ones.")

        if args.collection_name not in c.collections_dict:
            raise ConfigError("Are you sure collection name '%s' stored in this database? You know,\
                                probably it is not at all :/ You can always use --list-collections\
                                flag to see what is in there." % args.collection_name)

        # {'read_only': False, 'source_db_path': 'delete.db', 'num_splits': 1453, 'source_db_version': '13', 'bin_names': 'Bin_1', 'num_bins': 1}
        collection_info = c.collections_dict[args.collection_name]
        run.info('Collection found', '"%s" (%d items in %d bins)' % (args.collection_name, collection_info['num_splits'], collection_info['num_bins']))

        c.export_collection(args.collection_name, output_file_prefix=args.output_file_prefix, include_unbinned=args.include_unbinned)
    except ConfigError as e:
        print(e)
        sys.exit(1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(2)


def get_args():
    parser = ArgumentParser(description=__description__)

    parser.add_argument(*anvio.A('pan-or-profile-db'), **anvio.K('pan-or-profile-db'))
    parser.add_argument(*anvio.A('collection-name'), **anvio.K('collection-name'))
    parser.add_argument(*anvio.A('output-file-prefix'), **anvio.K('output-file-prefix'))
    parser.add_argument(*anvio.A('list-collections'), **anvio.K('list-collections'))
    parser.add_argument('--include-unbinned', default=False, action='store_true', help="When this flag is used, anvi'o\
                         will also store in the output file the items that do not appear in any of your bins. This new\
                         bin will be called 'UNBINNED_ITEMS_BIN'. Yes. The ugly name is intentional.")

    return parser.get_args(parser)

if __name__ == '__main__':
    main()
