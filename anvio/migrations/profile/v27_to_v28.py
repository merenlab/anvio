#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio.db as db
import anvio.terminal as terminal

from anvio.errors import ConfigError
from anvio.dbinfo import is_profile_db

run = terminal.Run()
progress = terminal.Progress()

current_version, next_version = [x[1:] for x in __name__.split('_to_')]


def divide_stackedbar_to_multiple_entries(table_dict):
    if len(table_dict):
        key_counter = max(table_dict.keys()) + 1

        for entry_id in list(table_dict.keys()):
            data_type = table_dict[entry_id]['data_type']

            if data_type == 'stackedbar':
                item_name = table_dict[entry_id]['item_name']
                data_keys = table_dict[entry_id]['data_key'].split('!')[1].split(';')
                data_values = table_dict[entry_id]['data_value'].split(';')
                data_type = table_dict[entry_id]['data_type']
                data_group = table_dict[entry_id]['data_group']

                stackbar_name = table_dict[entry_id]['data_key'].split('!')[0]

                for i in range(len(data_keys)):
                    table_dict[key_counter] = {'item_name': item_name, 
                                                'data_key': stackbar_name + '!' + data_keys[i],
                                                'data_value': data_values[i],
                                                'data_type': 'stackedbar',
                                                'data_group': data_group}
                    key_counter += 1

                del table_dict[entry_id]

    return table_dict


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    # make sure someone is not being funny
    is_profile_db(db_path)

    # make sure the version is accurate
    profile_db = db.DB(db_path, None, ignore_version = True)
    if str(profile_db.get_version()) != current_version:
        raise ConfigError("Version of this profile database is not %s (hence, this script cannot really do anything)." % current_version)

    for table_name in ['layer_additional_data', 'item_additional_data']:
        new_table = divide_stackedbar_to_multiple_entries(profile_db.get_table_as_dict(table_name))
        profile_db._exec("DELETE FROM '%s'" % table_name)
        for entry_id in new_table:
            profile_db.insert(table_name, tuple(new_table[entry_id].values()))

    #set the version
    profile_db.remove_meta_key_value_pair('version')
    profile_db.set_version(next_version)

    # bye
    profile_db.disconnect()
    progress.end()

    run.info_single("Your profile db is now %s (and anvi'o is as surprised as you are)." % next_version, nl_after=1, nl_before=1, mc='green')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade profile database from version %s to version %s' % (current_version, next_version))
    parser.add_argument('profile_db', metavar = 'PROFILE_DB', help = "An anvi'o profile database of version %s" % current_version)
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.profile_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
