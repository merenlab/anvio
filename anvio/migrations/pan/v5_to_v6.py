#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio.db as db
import anvio.tables as t
import anvio.dbops as dbops
import anvio.terminal as terminal 

from anvio.errors import ConfigError


run = terminal.Run()
progress = terminal.Progress()

current_version = '5'
next_version    = '6'


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    # make sure someone is not being funny
    dbops.is_pan_db(db_path)

    # make sure the version is accurate
    pan_db = db.DB(db_path, None, ignore_version = True)
    if str(pan_db.get_version()) != current_version:
        raise ConfigError("Version of this pan database is not %s (hence, this script cannot really do anything)." % current_version)

    progress.new("Trying to upgrade the pan database")
    progress.update('...')

    try:
        pan_db.create_table(t.item_orders_table_name, t.item_orders_table_structure, t.item_orders_table_types)
    except:
        pass

    clusterings = pan_db.get_table_as_dict('clusterings')

    # move clustering data into the new table
    for clustering in clusterings:
        newick = clusterings[clustering]['newick']
        pan_db._exec('''INSERT INTO %s VALUES (?,?,?)''' % t.item_orders_table_name, tuple([clustering, 'newick', newick]))

    # update keys
    for old_key, new_key in [('available_clusterings', 'available_item_orders'),
                             ('PCs_clustered', 'PCs_ordered'),
                             ('default_clustering', 'default_item_order')]:
        try:
            pan_db.set_meta_value(new_key, pan_db.get_meta_value(old_key))
        except:
            pass

    # remove stuff that are not irrelevant
    try:
        pan_db._exec('DROP TABLE clusterings;')
        pan_db.remove_meta_key_value_pair('available_clusterings')
        pan_db.remove_meta_key_value_pair('PCs_clustered')
        pan_db.remove_meta_key_value_pair('default_clustering')
    except:
        pass

    # commit
    try:
        pan_db._exec('COMMIT')
    except:
        pass

    # cleanup
    try:
        pan_db._exec('vacuum')
    except:
        pass

    # set the version
    pan_db.remove_meta_key_value_pair('version')
    pan_db.set_version(next_version)

    # bye
    pan_db.disconnect()
    progress.end()

    run.info_single('Your pan db is now %s.' % next_version, nl_after=1, nl_before=1, mc='green')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade pan database from version %s to version %s' % (current_version, next_version))
    parser.add_argument('pan_db', metavar = 'PAN_DB', help = "An anvi'o pan database of version %s" % current_version)
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.pan_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
