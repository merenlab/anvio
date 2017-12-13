#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio.db as db
import anvio.dbops as dbops
import anvio.terminal as terminal

from anvio.errors import ConfigError


run = terminal.Run()
progress = terminal.Progress()

current_version = '6'
next_version    = '7'

item_additional_data_table_name      = 'item_additional_data'
item_additional_data_table_structure = ['entry_id', 'item_name', 'key', 'value', 'type']
item_additional_data_table_types     = [ 'numeric',   'text'   , 'str',  'str' ,  'str']


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    # make sure someone is not being funny
    dbops.is_pan_db(db_path)

    # make sure the version is accurate
    pan_db = db.DB(db_path, None, ignore_version = True)
    if str(pan_db.get_version()) != current_version:
        raise ConfigError("Version of this pan database is not %s (hence, this script cannot really do anything)." % current_version)

    # update keys
    for old_key, new_key in [('maxbit', 'minbit')]:
        try:
            pan_db.set_meta_value(new_key, pan_db.get_meta_value(old_key))
        except:
            pass

    # remove stuff that are not irrelevant
    try:
        pan_db.remove_meta_key_value_pair('maxbit')
    except:
        pass


    # learn additional_data_headers for later:
    additional_data_headers = pan_db.get_meta_value('additional_data_headers').split(',')

    # take care of the self table
    self_table = pan_db.get_table_as_list_of_tuples('self')
    pan_db.cursor.execute('ALTER TABLE self RENAME TO self_TEMP;')
    pan_db.cursor.execute('CREATE TABLE self (key text, value text);')

    for key, val in self_table:
        new_key = key.replace('PC', 'gene_cluster').replace('pc', 'gene_cluster')
        new_val = val.replace('PC', 'gene_cluster').replace('pc', 'gene_cluster')

        pan_db.set_meta_value(new_key, new_val)

    pan_db.cursor.execute('DROP TABLE self_TEMP;')

    # take care of the views table
    views_table = pan_db.get_table_as_list_of_tuples('views')
    pan_db.cursor.execute('ALTER TABLE views RENAME TO views_TEMP;')
    pan_db.cursor.execute('CREATE TABLE views (view_id str, target_table str);')

    values = []
    for view, target in views_table:
        new_view = view.replace('PC', 'gene_cluster').replace('pc', 'gene_cluster')
        new_target = target.replace('PC', 'gene_cluster').replace('pc', 'gene_cluster')
        values.append((new_view, new_target),)

    pan_db.insert_many('views', values)

    pan_db.cursor.execute('DROP TABLE views_TEMP;')

    # rename tables
    pan_db._exec('ALTER TABLE PC_frequencies RENAME TO gene_cluster_frequencies;')
    pan_db._exec('ALTER TABLE PC_presence_absence RENAME TO gene_cluster_presence_absence;')
    pan_db._exec('ALTER TABLE protein_clusters RENAME TO gene_clusters;')


    # protein_cluster_id -> gene_cluster_id in table gene_clusters.
    pan_db.cursor.execute('ALTER TABLE gene_clusters RENAME TO gene_clusters_TEMP;')
    pan_db.cursor.execute('CREATE TABLE gene_clusters (entry_id numeric, gene_caller_id numeric, gene_cluster_id str, genome_name str, alignment_summary str);')
    pan_db.cursor.execute('INSERT INTO gene_clusters(entry_id, gene_caller_id, gene_cluster_id, genome_name, alignment_summary) SELECT entry_id, gene_caller_id, protein_cluster_id, genome_name, alignment_summary FROM gene_clusters_TEMP;')
    pan_db.cursor.execute('DROP TABLE gene_clusters_TEMP;')

    # commit
    try:
        pan_db._exec('COMMIT')
    except:
        pass

    # we also added a totally new table to this version:
    pan_db.create_table(item_additional_data_table_name, item_additional_data_table_structure, item_additional_data_table_types)

    # set the version
    pan_db.remove_meta_key_value_pair('version')
    pan_db.set_version(next_version)

    # we have one more thing to do: getting rid of the 'additional_data' table without losing data, by carrying
    # its content into our new item_additional_data_table
    additional_data_table_dict = pan_db.get_table_as_dict('additional_data')

    # close the db temporarily
    pan_db.disconnect()

    # update the contents of the item_additional_data_table
    args = argparse.Namespace(pan_db=db_path, just_do_it=True)
    item_additional_data_table = dbops.TableForItemAdditionalData(args)
    item_additional_data_table.add(additional_data_headers, additional_data_table_dict)

    # open the database again to remove stuff
    pan_db = db.DB(db_path, None, ignore_version = True)
    pan_db.remove_meta_key_value_pair('additional_data_headers')
    pan_db._exec("DROP TABLE additional_data")

    # now bye for real!
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
