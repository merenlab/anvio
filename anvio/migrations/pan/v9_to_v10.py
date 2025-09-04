#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse
from ete3 import Tree

import anvio.db as db
import anvio.terminal as terminal

from anvio.errors import ConfigError
from anvio.dbinfo import is_pan_db


run = terminal.Run()
progress = terminal.Progress()

current_version, next_version = [x[1:] for x in __name__.split('_to_')]

item_orders_table_name               = 'item_orders'
item_orders_table_structure         = ['name', 'type', 'data']
item_orders_table_types             = ['text', 'text', 'text']

layer_orders_table_name              = 'layer_orders'
layer_orders_table_structure         = ['data_key', 'data_type', 'data_value']
layer_orders_table_types             = [  'text'  ,    'text'  ,    'text'   ]


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    # make sure someone is not being funny
    is_pan_db(db_path)

    # make sure the version is accurate
    pan_db = db.DB(db_path, None, ignore_version = True)
    if str(pan_db.get_version()) != current_version:
        raise ConfigError("Version of this pan database is not %s (hence, this script cannot really do anything)." % current_version)

    # migrate item orders
    item_orders = pan_db.get_table_as_dict(item_orders_table_name)
    for order_name in item_orders:
        if item_orders[order_name]['type'] == 'newick':
            newick = Tree(item_orders[order_name]['data'], format=1)
            newick = newick.write(format=2)
            pan_db._exec("""UPDATE %s SET "data" = ? WHERE "name" LIKE ?""" % item_orders_table_name, (newick, order_name))

    # migrate layer orders
    layer_orders = pan_db.get_table_as_dict(layer_orders_table_name)
    for order_name in layer_orders:
        if layer_orders[order_name]['data_type'] == 'newick':
            newick = Tree(layer_orders[order_name]['data_value'], format=1)
            newick = newick.write(format=2)
            pan_db._exec("""UPDATE %s SET "data_value" = ? WHERE "data_key" LIKE ?""" % layer_orders_table_name, (newick, order_name))

    # set the version
    pan_db.remove_meta_key_value_pair('version')
    pan_db.set_version(next_version)

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
