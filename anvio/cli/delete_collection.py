#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.terminal as terminal
import anvio.ccollections as ccollections

from anvio.errors import ConfigError, FilesNPathsError
from anvio.tables.collections import TablesForCollections


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ["profile-db", "collection"]
__description__ = "Remove a collection from a given profile database"


@terminal.time_program
def main():
    args = get_args()
    run = terminal.Run()

    A = lambda x: args.__dict__[x] if x in args.__dict__ else None
    profile_db_path = A('profile_db')
    list_collections = A('list_collections')
    collection_name = A('collection_name')

    try:
        collections = ccollections.Collections()
        collections.populate_collections_dict(profile_db_path)

        if list_collections:
            collections.list_collections()
            sys.exit()

        if not collection_name:
            raise ConfigError("So you are asking anvi'o to delete a collection from the profile database "
                               "without actually providing a collection name. Very cute. You are free to "
                               "go ahead and try again.")

        if collection_name not in collections.collections_dict:
            raise ConfigError("The collection name '%s' is not in this database :/" % collection_name)

        collections_table = TablesForCollections(profile_db_path)
        collections_table.delete(collection_name)
        run.info('Collection delete', 'Done. The collection "%s" is no more.' % (collection_name))
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    parser.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db'))
    parser.add_argument(*anvio.A('collection-name'), **anvio.K('collection-name'))
    parser.add_argument(*anvio.A('list-collections'), **anvio.K('list-collections'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
