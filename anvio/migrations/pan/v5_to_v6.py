#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio.db as db
import anvio.terminal as terminal

from anvio.errors import ConfigError
from anvio.dbinfo import is_pan_db


run = terminal.Run()
progress = terminal.Progress()

current_version = '5'
next_version    = '6'

item_orders_table_name               = 'item_orders'
item_orders_table_structure          = ['name', 'type', 'data']
item_orders_table_types              = [ 'str',  'str',  'str']

def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    # make sure someone is not being funny
    is_pan_db(db_path)

    # make sure the version is accurate
    pan_db = db.DB(db_path, None, ignore_version = True)
    if str(pan_db.get_version()) != current_version:
        raise ConfigError("Version of this pan database is not %s (hence, this script cannot really do anything)." % current_version)

    progress.new("Trying to upgrade the pan database")
    progress.update('...')

    try:
        pan_db.create_table(item_orders_table_name, item_orders_table_structure, item_orders_table_types)
    except:
        pass

    clusterings = pan_db.get_table_as_dict('clusterings')

    # move clustering data into the new table
    for clustering in clusterings:
        newick = clusterings[clustering]['newick']
        pan_db._exec('''INSERT INTO %s VALUES (?,?,?)''' % item_orders_table_name, tuple([clustering, 'newick', newick]))

    # update keys
    for old_key, new_key in [('pc_min_occurrence', 'gene_cluster_min_occurrence'),
                             ('num_protein_clusters', 'num_gene_clusters'),
                             ('num_genes_in_protein_clusters', 'num_genes_in_gene_clusters'),
                             ('available_clusterings', 'available_item_orders'),
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
        pan_db.remove_meta_key_value_pair('num_protein_clusters')
        pan_db.remove_meta_key_value_pair('num_genes_in_protein_clusters')
        pan_db.remove_meta_key_value_pair('pc_min_occurrence')
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
