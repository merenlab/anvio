# -*- coding: utf-8
# pylint: disable=line-too-long

"""Implements the item additional data class.

This is the module where we maintain the item additional data tables in anvi'o
profile databases. Related issue is here: https://github.com/merenlab/anvio/issues/662.
"""

from collections import Counter

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.tableops import Table


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2018, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


run = terminal.Run()
progress = terminal.Progress()


class ItemAdditionalData(Table):
    def __init__(self, args, r=run, p=progress):
        self.run = r
        self.progress = p

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.db_path = A('pan_or_profile_db') or A('profile_db') or A('pan_db')
        self.just_do_it = A('just_do_it')

        if not self.db_path:
            raise ConfigError("ItemAdditionalData class is inherited with args object that did not\
                               contain any database path :/ Even though any of the following would\
                               have worked: `pan_or_profile_db`, `profile_db`, `pan_db` :(")

        dbops.is_pan_or_profile_db(self.db_path)
        self.db_type = dbops.get_db_type(self.db_path)
        self.db_version = dbops.get_required_version_for_db(self.db_path)

        database = db.DB(self.db_path, self.db_version)
        self.item_additional_data_keys = sorted(database.get_single_column_from_table(t.item_additional_data_table_name, 'key'))
        database.disconnect()

        Table.__init__(self, self.db_path, self.db_version, self.run, self.progress)



    def get(self):
        """Will return the additional data what is in the database."""

        database = db.DB(self.db_path, dbops.get_required_version_for_db(self.db_path))
        item_additional_data = database.get_table_as_dict(t.item_additional_data_table_name)
        item_additional_data_keys = sorted(database.get_single_column_from_table(t.item_additional_data_table_name, 'key'))
        item_names = database.get_single_column_from_table(t.item_additional_data_table_name, 'split_name')
        database.disconnect()

        if not len(item_names):
            return None, {}

        d = dict.fromkeys(item_names, Counter({}))

        for entry in item_additional_data.values():
            split_name = entry['split_name']
            key = entry['key']
            value = entry['value']

            if entry['type'] in ['int', 'float']:
                d[split_name][key] = eval(entry['type'])(value)
            else:
                d[split_name][key] = value

        for split_name in d:
            for key in item_additional_data_keys:
                if key not in d[split_name]:
                    d[split_name][key] = None

        return item_additional_data_keys, d


    def remove(self, keys_list):
        '''Give this guy a list of key, and watch their demise.'''

        if not isinstance(keys_list, list):
            raise ConfigError("The remove function in ItemAdditionalData class wants you to watch\
                               yourself before you wreck yourself. In other words, can you please\
                               make sure the keys you send is of type `list` thankyouverymuch.")

        for key in keys_list:
            if key not in self.item_additional_data_keys:
                # what the hell, user?
                return

            database = db.DB(self.db_path, dbops.get_required_version_for_db(self.db_path))
            database._exec('''DELETE from %s WHERE key="%s"''' % (t.item_additional_data_table_name, key))
            database.disconnect()

        self.run.warning("Data for the following keys removed from the database #SAD: %s" % (', '.join(keys_list)))


    def add(self, keys_list, data_dict):
        """Main function to add data into the item additional data table.

           * `data_dict`: a dictionary that should follow this format:

                d = {
                        'split_name_01': {'key_01': value,
                                          'key_02': value,
                                          'key_03': value
                                          },
                        'split_name_02': {'key_01': value,
                                          'key_03': value,
                                          },
                        (...)
                    }

           * `keys_list`: is a list of keys one or more of which should appear for each item
                          in `data_dict`.
        """

        if not isinstance(keys_list, list):
            raise ConfigError("List of keys must be of type `list`. Go away.")

        if not isinstance(data_dict, dict):
            raise ConfigError("Nope. Your data must be of type `dict`.")

        self.run.warning(None, 'New additional data...', lc="yellow")
        key_types = {}
        for key in keys_list:
            if '!' in key:
                predicted_key_type = "stackedbar"
            else:
                type_class = utils.get_predicted_type_of_items_in_a_dict(data_dict, key)
                predicted_key_type = type_class.__name__ if type_class else None

            key_types[key] = predicted_key_type
            self.run.info('Key "%s"' % key, 'Predicted type: %s' % (key_types[key]))

        # we be responsible here.
        keys_already_in_db = [c for c in keys_list if c in self.item_additional_data_keys]
        if len(keys_already_in_db):
            if self.just_do_it:
                self.run.warning('The following keys in your data dict will replace the ones that are already\
                                  in your %s database: %s.' % (self.db_type, ', '.join(keys_already_in_db)))

                self.remove(keys_already_in_db)
            else:
                raise ConfigError("Some of the keys in your new data appear to be in the database already. If you\
                                   want to replace those in the database with the ones in your new data use the\
                                   `--just-do-it` flag, and watch anvi'o make an exception just for you and complain\
                                   about nothin' for this once.")

        db_entries = []
        self.set_next_available_id(t.item_additional_data_table_name)
        for item_name in data_dict:
            for key in data_dict[item_name]:
                db_entries.append(tuple([self.next_id(t.item_additional_data_table_name),
                                         item_name,
                                         key,
                                         data_dict[item_name][key],
                                         key_types[key]]))

        database = db.DB(self.db_path, dbops.get_required_version_for_db(self.db_path))
        database._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?)''' % t.item_additional_data_table_name, db_entries)
        database.disconnect()

        self.run.info('New data added to the db', '%s.' % (', '.join(keys_list)))


    def populate_from_file(self, additional_data_file_path):
        filesnpaths.is_file_tab_delimited(additional_data_file_path)

        keys = utils.get_columns_of_TAB_delim_file(additional_data_file_path)
        data = utils.get_TAB_delimited_file_as_dictionary(additional_data_file_path)

        if not len(keys):
            raise ConfigError("There is something wrong with the additional data file at %s.\
                               It does not seem to have any additional keys for data :/" \
                                            % (additional_data_file_path))

        self.add(keys, data)

