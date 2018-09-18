# -*- coding: utf-8
# pylint: disable=line-too-long

"""The fancy additioanl data module"""

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.tables.tableops import Table


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class AdditionalAndOrderDataBaseClass(Table, object):
    """This is a base class for common operations between order and additional data classes."""

    def __init__(self, args):
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.db_path = A('pan_or_profile_db') or A('profile_db') or A('pan_db')
        self.just_do_it = A('just_do_it')
        self.target_data_group_set_by_user = A('target_data_group') or None
        self.target_data_group = self.target_data_group_set_by_user or 'default'

        if not self.db_path:
            raise ConfigError("The AdditionalAndOrderDataBaseClass is inherited with an args object that did not\
                               contain any database path :/ Even though any of the following would\
                               have worked: `pan_or_profile_db`, `profile_db`, `pan_db` :(")

        if not self.table_name:
            raise ConfigError("The AdditionalAndOrderDataBaseClass does not know anything about the table it should\
                               be working with.")

        utils.is_pan_or_profile_db(self.db_path)
        self.db_type = utils.get_db_type(self.db_path)
        self.db_version = utils.get_required_version_for_db(self.db_path)

        database = db.DB(self.db_path, self.db_version)
        self.additional_data_keys = database.get_single_column_from_table(self.table_name, 'data_key')
        database.disconnect()

        Table.__init__(self, self.db_path, self.db_version, self.run, self.progress)

        self.nulls_per_type = {'str': '',
                               'int': 0,
                               'float': 0,
                               'stackedbar': 0,
                               'unknown': None}


    def populate_from_file(self, additional_data_file_path, skip_check_names=None):

        if skip_check_names is None and utils.is_blank_profile(self.db_path):
            # FIXME: this BS is here because blank abvi'o profiles do not know what items they have,
            #        hence the utils.get_all_item_names_from_the_database function eventually explodes if we
            #        don't skip check names.
            skip_check_names = True

        filesnpaths.is_file_tab_delimited(additional_data_file_path)

        data_keys = utils.get_columns_of_TAB_delim_file(additional_data_file_path)
        data_dict = utils.get_TAB_delimited_file_as_dictionary(additional_data_file_path)

        if not len(data_keys):
            raise ConfigError("There is something wrong with the additional data file for %s at %s.\
                               It does not seem to have any additional keys for data :/" \
                                            % (self.target_table, additional_data_file_path))

        if self.target_table == 'layer_orders':
            OrderDataBaseClass.add(self, data_dict, skip_check_names)
        else:
            AdditionalDataBaseClass.add(self, data_dict, data_keys, skip_check_names)


    def remove(self, data_keys_list):
        '''Give this guy a list of key for additional data, and watch their demise.'''

        if not isinstance(data_keys_list, list):
            raise ConfigError("The remove function in AdditionalDataBaseClass wants you to watch\
                               yourself before you wreck yourself. In other words, can you please\
                               make sure the keys you send is of type `list` thankyouverymuch?")

        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))

        additional_data_keys = sorted(database.get_single_column_from_table(self.table_name, 'data_key', unique=True))

        if not len(additional_data_keys):
            self.run.info_single('There is nothing to remove --the %s additional data table is already empty :(' % self.target_table)
            database.disconnect()

            return

        missing_keys = [k for k in data_keys_list if k not in additional_data_keys]
        if len(missing_keys) and not self.just_do_it:
            database.disconnect()
            raise ConfigError("The following keys you wanted to remove from the items additional data table are\
                               not really in the table: '%s'. Anvi'o is confused :/" % (', '.join(missing_keys)))

        if data_keys_list:
            for key in data_keys_list:
                if key not in additional_data_keys:
                    # what the hell, user?
                    return

                if 'data_group' in database.get_table_structure(self.table_name):
                    database._exec('''DELETE from %s WHERE data_key="%s" and data_group="%s"''' % (self.table_name, key, self.target_data_group))
                else:
                    database._exec('''DELETE from %s WHERE data_key="%s"''' % (self.table_name, key))

            self.run.warning("Data from the table '%s' for the following data keys in data group '%s' \
                              removed from the database: '%s'. #SAD." % (self.target_table, self.target_data_group, ', '.join(data_keys_list)))
        else:
            if not self.just_do_it:
                raise ConfigError("You did not provide a list of data keys to remove, which means you are about to delete everything in the\
                                   %s additional data table. Just to be on the safe side, anvi'o is looking for a confirmation. If you\
                                   try again with the --just-do-it flag, anvi'o will put on its business socks, and burn this table\
                                   and everything in it to the ground." % self.target_table)

            database._exec('''DELETE from %s''' % (self.table_name))

            self.run.warning("All data from the %s additional data table is removed (ouch)." % self.target_table)

        database.disconnect()


    def export(self, output_file_path):
        filesnpaths.is_output_file_writable(output_file_path)

        if self.target_table in ['layers', 'items']:
            keys, data = AdditionalDataBaseClass.get(self)
            if keys:
                if len(self.available_group_names) - 1:
                    self.run.warning("You are exporting data from the additional data table '%s' for the\
                                      data group '%s'. Great. Just remember that there are %d more data\
                                      groups in your database, and you are not exporting anything from them\
                                      at this point (they know you're the boss, so they're not upset)." \
                                        % (self.target_table, self.target_data_group, len(self.available_group_names) - 1),
                                      header="FRIENDLY REMINDER", lc='yellow')

                self.run.info('Target data group', self.target_data_group, mc='green')

        elif self.target_table in ['layer_orders']:
            data = OrderDataBaseClass.get(self, native_form=True)
            keys = ['data_type', 'data_value']
        else:
            raise ConfigError("Your target table '%s' does not make any sense" % self.target_table)

        if not(len(data)):
            raise ConfigError("Additional data table for %s is empty. There is nothing to export :/" % self.target_table)

        utils.store_dict_as_TAB_delimited_file(data, output_file_path, headers=[self.target_table] + keys)

        self.run.info('Target data table', self.target_table)
        self.run.info('Output file', output_file_path)


    def list_data_keys(self):
        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))

        NOPE = lambda: self.run.info_single("There are no additional data for '%s' in this database :/" \
                                                    % (self.target_table), nl_before=1, nl_after=1, mc='red')

        additional_data_keys = {}
        # here is where things get tricky. if we are dealing with additional data layers or items, we will have
        # data groups that are not relevant for order data. this will affect the listing of data keys in either
        # of these table types. hence we get group names first here, and then will do a bunch of if/else checks
        # based on their availability
        if self.target_table in ['layers', 'items']:
            if not self.available_group_names:
                NOPE()
                return
            else:
                self.check_target_data_group()

            # if the user set a target data group, let's focus on that here. otherwise we will
            # use all data groups available for a messy output.
            group_names = [self.target_data_group] if self.target_data_group_set_by_user else self.available_group_names

            for group_name in group_names:
                data_keys_in_group = database.get_single_column_from_table(self.table_name, \
                                                                           'data_key', \
                                                                           unique=True, \
                                                                           where_clause="data_group='%s'" % group_name)
                additional_data_keys[group_name] = sorted(data_keys_in_group)

        elif self.target_table in ['layer_orders']:
            data_keys = sorted(database.get_single_column_from_table(self.table_name, 'data_key', unique=True))

            if not len(data_keys):
                self.run.info_single("There are no additional data for '%s' in this database :/" \
                                                    % (self.target_table), nl_before=1, nl_after=1, mc='red')
                database.disconnect()
                return

            additional_data_keys['default'] = data_keys
            group_names = ['default']

        self.run.warning('', 'DATA KEYS FOR "%s" in %d DATA GROUP(S)' % (self.target_table.upper(), len(group_names)), lc='green')

        for group_name in group_names:
            num_keys = len(additional_data_keys[group_name])

            self.run.info_single('DATA GROUP "%s" WITH %d KEYS' % (group_name, num_keys), nl_before = 1)

            if anvio.DEBUG:
                num_keys_to_display = num_keys
            else:
                num_keys_to_display = min([5, num_keys])

            for key_index in range(0, num_keys_to_display):
                data_key = additional_data_keys[group_name][key_index]
                rows = database.get_some_rows_from_table_as_dict(self.table_name, 'data_key="%s"' % data_key)

                if self.target_table == 'layer_orders':
                    self.run.info_single('%s (%s)' % (data_key, list(rows.values())[0]['data_type']),
                                         nl_after = 1 if data_key == additional_data_keys[group_name][-1] else 0, level=2)
                else:
                    self.run.info_single('%s (%s, describes %d %s)' % (data_key, list(rows.values())[0]['data_type'], len(rows), self.target_table),
                                         nl_after = 1 if data_key == additional_data_keys[group_name][-1] else 0, level=2)

            num_keys_not_displayed = num_keys - num_keys_to_display
            if num_keys_not_displayed > 0:
                self.run.info_single('(... %d more; use `--debug` to list all ...)' % \
                                                                (num_keys_not_displayed), nl_after = 1, mc='cyan', level=3)

        database.disconnect()


    def data_dict_sanity_check(self, data_dict, data_keys_list=None, treat_data_dict_as=None):
        data_dict_type = treat_data_dict_as or self.target_table

        if not isinstance(data_dict, dict):
            raise ConfigError("Nope. Your data must be of type %s, but it is a %s." % (type(dict()), type(data_dict)))

        if not len(data_dict):
            raise ConfigError("The data sent for sanity check seems to be empty.")

        if data_keys_list and not isinstance(data_keys_list, list):
            raise ConfigError("List of keys must be of type `list`. Go away (and come back).")

        utils.check_misc_data_keys_for_format(data_keys_list)

        # FIXME: we have two controls here. The first one is how we work with order data natively. The second one is how it
        #        looks like when it is read through the .get() member function of the TableForLayerOrders because rest of
        #        anvi'o does not know how to work with the native format. This should be fixed by teaching the rest of anvi'o
        #        how to work with order data dicts in the native form.
        looks_like_layer_orders = sorted(list(data_dict.values())[0].keys()) == sorted(['data_type', 'data_value']) or \
                                  sorted(list(data_dict.values())[0].keys()) == sorted(['basic', 'newick'])

        if looks_like_layer_orders and data_dict_type is not 'layer_orders':
            raise ConfigError("The data you sent here seems to describe an order, but you want anvi'o to treat it\
                               as additional data for %s. Not cool." % self.target_table)

        if not looks_like_layer_orders and data_dict_type is 'layer_orders':
            raise ConfigError("The data that claims to be a layer order data do not seem to be one.")

        if data_keys_list:
            for item_or_layer_name in data_dict:
                for key in data_keys_list:
                    if key not in data_dict[item_or_layer_name]:
                        raise ConfigError("Your data dictionary fails the sanity check since at least one item in it (i.e., %s) is\
                                           missing any data for the key '%s'." % (item_or_layer_name, key))


class OrderDataBaseClass(AdditionalAndOrderDataBaseClass, object):
    """Implements a base class to deal with tables that keep order data."""

    def __init__(self, args):
        AdditionalAndOrderDataBaseClass.__init__(self, args)


    def get(self, additional_data_keys=None, additional_data_dict=None, native_form=False):
        """Will return the layer order data dict.

           If `additional_data_keys` and `additional_data_dict` variables are provided, it will return an order dict
           with additioanl orders generated from data found in those. An example:

                >>> layers_additional_data_table = TableForLayerAdditionalData(args)
                >>> layer_orders_table = TableForLayerOrders(args)
                >>> layer_additional_data_keys, layer_additional_data_dict = layers_additional_data_table.get()
                >>> layer_orders_data_dict = layer_orders_table.get(layer_additional_data_keys, layer_additional_data_dict

           The `native_form` is tricky. Normaly in this function we turn the data read from the table into
           the legacy data structure anvi'o has been using before all these were implemented. So, by default,
           we return the data in that legacy format so everyting else would continue working. But if the
           programmer is interested in exporting the order data as an output file, we don't want to follow that
           silly legacy format. So the parameter `native_form` simply returns the data in native form. Export
           functions use the native form. At some point the rest of the anvi'o codebase should be fixed to
           work only with the native form, and there is no need to this tyranny. Here is a FIXME for that
           in case one day someone wants to address that.
        """

        self.progress.new('Recovering layer order data')
        self.progress.update('...')
        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))
        order_data = database.get_table_as_dict(self.table_name)
        database.disconnect()
        self.progress.end()

        if native_form:
            return order_data

        d = {}
        for order_name in order_data:
            data_type = order_data[order_name]['data_type']
            data_value = order_data[order_name]['data_value']
            if data_type == 'newick':
                d[order_name] = {'newick': data_value, 'basic': None}
            elif data_type == 'basic':
                d[order_name] = {'newick': None, 'basic': data_value}
            else:
                raise ConfigError("Something is wrong :( Anvi'o just found an entry in %s table with an unrecognized\
                                   type of '%s'." % (self.target_table, data_type))

        if additional_data_keys and additional_data_dict:
            return self.update_orders_dict_using_additional_data_dict(d, additional_data_keys, additional_data_dict)
        else:
            return d


    def update_orders_dict_using_additional_data_dict(self, order_data_dict, additional_data_keys, additional_data_dict):
        if order_data_dict:
            self.data_dict_sanity_check(order_data_dict, treat_data_dict_as='layer_orders')

        # FIXME: here we need to check whether the two dictionaries are in fact 'compatible' with respect to sample names
        #        they describe.
        self.data_dict_sanity_check(additional_data_dict, data_keys_list=additional_data_keys, treat_data_dict_as='layers')

        sum_stackbar_items = {}
        for data_key in additional_data_keys:
            if '!' in data_key:
                stackbar_name = data_key.split('!')[0]

                if stackbar_name not in sum_stackbar_items:
                    sum_stackbar_items[stackbar_name] = {}

                for layer in additional_data_dict:
                    if layer not in sum_stackbar_items[stackbar_name]:
                        sum_stackbar_items[stackbar_name][layer] = 0.0

                    if additional_data_dict[layer][data_key]:
                        sum_stackbar_items[stackbar_name][layer] += float(additional_data_dict[layer][data_key])

        for data_key in additional_data_keys:
            if '!' in data_key:
                predicted_key_type = "stackedbar"
                stacked_bar_name, item_name = data_key.split('!')
                data_key_name = '%s [%s]' % (stacked_bar_name, item_name)
            else:
                type_class = utils.get_predicted_type_of_items_in_a_dict(additional_data_dict, data_key)
                predicted_key_type = type_class.__name__ if type_class else 'unknown'
                data_key_name = data_key

            if predicted_key_type == "stackedbar":
                stackbar_name = data_key.split('!')[0]
                layer_name_layer_data_tuples = []
                for layer in additional_data_dict:
                    if additional_data_dict[layer][data_key]:
                        if sum_stackbar_items[stackbar_name][layer] == 0:
                            layer_name_layer_data_tuples.append((0, layer))
                        else:
                            layer_name_layer_data_tuples.append(((float(additional_data_dict[layer][data_key]) / (1.0 * sum_stackbar_items[stackbar_name][layer])), layer))
                    else:
                        layer_name_layer_data_tuples.append((self.nulls_per_type[predicted_key_type], layer))
            else:
                layer_name_layer_data_tuples = [(additional_data_dict[layer][data_key] if additional_data_dict[layer][data_key] else self.nulls_per_type[predicted_key_type], layer) for layer in additional_data_dict]

            order_data_dict['>> ' + data_key_name] = {'newick': None, 'basic': ','.join([t[1] for t in sorted(layer_name_layer_data_tuples)])}
            order_data_dict['>> ' + data_key_name + ' (reverse)'] = {'newick': None, 'basic': ','.join([t[1] for t in sorted(layer_name_layer_data_tuples, reverse=True)])}

        return order_data_dict


    def add(self, data_dict, skip_check_names=False):
        """
            The default function to add data into the orders table.

             * `data_dict`: this variable for layer orders is expected to follow this format:

                  d = {
                          'data_key_01': {'data_type': 'newick',
                                          'data_value': '(item_or_layer_name_01:0.0370199,(item_or_layer_name_02:0.0227268,item_or_layer_name_01:0.0227268)Int3:0.0370199);'
                                          },
                          'data_key_02': {'data_type': 'basic',
                                          'data_value': 'item_or_layer_name_02,item_or_layer_name_01,item_or_layer_name_03'
                                          },
                          (...)
                  }
        """

        self.data_dict_sanity_check(data_dict)

        if self.target_table not in ['layer_orders']:
            raise ConfigError("You are using an OrderDataBaseClass instance to add %s data into your %s database. This is\
                               illegal and if you are here, it means someone made a mistake somewhere. If you are a user,\
                               check your flags to make sure you are targeting the right data table. If you are a programmer,\
                               you are fired." % (self.target_table, self.db_type))

        self.run.warning(None, 'New %s data...' % self.target_table, lc="yellow")
        data_keys_list = list(data_dict.keys())
        data_key_types = {}
        for key in data_keys_list:
            predicted_key_type = data_dict[key]['data_type']

            data_key_types[key] = predicted_key_type
            self.run.info('Data key "%s"' % key, 'Type: %s' % (data_key_types[key]), \
                                            nl_after = 1 if key == data_keys_list[-1] else 0)

        # we be responsible here.
        keys_already_in_db = [c for c in data_keys_list if c in self.additional_data_keys]
        if len(keys_already_in_db):
            if self.just_do_it:
                self.run.warning('The following keys in your data dict will replace the ones that are already\
                                  in your %s database: %s.' % (self.db_type, ', '.join(keys_already_in_db)))

                self.remove(keys_already_in_db)
            else:
                run.info('Data keys already in the db', ', '.join(keys_already_in_db), nl_before=2, mc='red')

                raise ConfigError("Some of the keys in your new order data appear to be in the database already. If you\
                                   want to replace those in the database with the ones in your new data use the\
                                   `--just-do-it` flag.")

        if skip_check_names:
            self.run.warning("You (or the programmer) asked anvi'o to NOT check the consistency of the names of your %s\
                              between your additional data and the %s database you are attempting to update. So be it.\
                              Anvi'o will not check anything, but if things don't look the way you expected them to look,\
                              you will not blame anvi'o for your poorly prepared data, but choose between yourself or\
                              Obama." % (self.target_table, self.db_type))
        else:
            TableForLayerOrders.check_names(self, data_dict)

        db_entries = []
        for item_name in data_dict:
            db_entries.append(tuple([item_name,
                                     data_dict[item_name]['data_type'],
                                     data_dict[item_name]['data_value']]))

        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))
        database._exec_many('''INSERT INTO %s VALUES (?,?,?)''' % self.table_name, db_entries)
        database.disconnect()

        self.run.info('New order data added to the db for %s' % self.target_table, '%s.' % (', '.join(data_keys_list)))


    def get_layer_names(self, data_dict):
        layer_names = {}

        for data_key in data_dict:
            try:
                if data_dict[data_key]['data_type'] == 'newick':
                    layer_names[data_key] = utils.get_names_order_from_newick_tree(data_dict[data_key]['data_value'])
                else:
                    layer_names[data_key] = [s.strip() for s in data_dict[data_key]['data_value'].split(',')]
            except:
                raise ConfigError("Parsing the %s data for %s failed :/ We don't know why, because we are lazy. Please\
                                   take a loook at your input data and figure out :(" \
                                                                    % (data_dict[data_key]['data_type'], data_key))

        return layer_names


class AdditionalDataBaseClass(AdditionalAndOrderDataBaseClass, object):
    """Implements additional data ops base class.

       See TableForItemAdditionalData or TableForLayerAdditionalData for usage example.

       See AdditionalAndOrderDataBaseClass for inherited functionality.
    """

    def __init__(self, args):
        AdditionalAndOrderDataBaseClass.__init__(self, args)

        self.available_group_names = self.get_group_names()


    def check_target_data_group(self):
        """A function to check whether the data group set is among the available ones.

           The reason this function is a separate one and it is not being called by the
           init function of the base class is because the user *can* set a new data group
           name when they want to 'import' things into the database. So, while exporting
           data or displaying data will require the requested data group to be a proper one,
           it is better to check that explicitly."""

        if self.target_data_group not in self.available_group_names:
            raise ConfigError("You (or the programmer) requested to initiate the additional data table for '%s' with\
                               the data group '%s', which is not really in that table :/ If it helps at all,\
                               the target table happened to have these ones instead: %s. What to do now? If you are\
                               here becasue the last command you run was something like 'show me all the data in misc'\
                               data tables, then you may try to be more specific by explicitly defining your target\
                               data table. If you think you have already been as sepecific as you could be, then anvi'o\
                               is as frustrated as you are right now :(" %\
                                    (self.target_table, self.target_data_group, ', '.join(['"%s"' % d for d in self.available_group_names])))

    def get_available_data_keys(self):
        """Will only return the additional data keys so the client can do some controls."""

        if not self.target_data_group:
            raise ConfigError("The target data group is not set. Which should never be the case at this stage. You shall\
                               not break anvi'o and go back to where you came from, devil :(")

        self.progress.new('Recovering additional keys and data for %s' % self.target_table)
        self.progress.update('...')

        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))
        additional_data_keys_in_db = database.get_single_column_from_table(self.table_name, 'data_key', unique=True, \
                        where_clause="""data_group LIKE '%s'""" % self.target_data_group)
        database.disconnect()

        self.progress.end()

        return additional_data_keys_in_db


    def get(self, additional_data_keys_requested=[]):
        """Will return the additional data keys as well as the data dict."""

        if not self.target_data_group:
            raise ConfigError("It seems the target data group is not set, which makes zero sense and should never happen\
                               unless you are doing some hacker stuff. Are you doing hacker stuff? Awesome! Tell us about\
                               it!")

        if not isinstance(additional_data_keys_requested, list):
            raise ConfigError("The `get` function in AdditionalDataBaseClass is upset with you. You could change that\
                               by making sure you request additional data keys with a variable of type `list`.")

        additional_data_keys_in_db = self.get_available_data_keys()

        self.progress.new('Recovering additional data for %s' % self.target_table)
        self.progress.update('...')

        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))
        if not len(additional_data_keys_requested):
            additional_data_keys = additional_data_keys_in_db
            additional_data = database.get_some_rows_from_table_as_dict(self.table_name,
                                                                        where_clause = """data_group LIKE '%s'""" % self.target_data_group,
                                                                        error_if_no_data=False)
        else:
            if not len(additional_data_keys_in_db):
                raise ConfigError("The %s database at %s does not contain any additional data for its %s to return. Usually this\
                                   would not have resulted in an exception, and anvi'o would simply returned an empty data\
                                   dictionary. But since you are requested to get a specific list of keys ('%s'), we are raising\
                                   an exception and break the flow just to make sure you are aware of the fact that we don't\
                                   see the stuff you're requesting :(" \
                                        % (self.db_type, self.db_path, self.target_table, ' ,'.join(additional_data_keys_requested)))

            if set(additional_data_keys_requested) - set(additional_data_keys_in_db):
                raise ConfigError("The keys you requested does not seem to appear in the additional data table of this %s db\
                                   at '%s' :/ Here is the list of keys you requested: '%s'. And here is the list of keys that anvi'o\
                                   knows about: '%s'." % (self.db_type, self.db_path,
                                                          ', '.join(additional_data_keys_requested),
                                                          ', '.join(additional_data_keys_in_db)))

            additional_data = database.get_some_rows_from_table_as_dict(self.table_name,
                                                where_clause = """data_group LIKE '%s' and data_key IN (%s)""" % (self.target_data_group,
                                                                                                                  ",".join(['"%s"' % key for key in additional_data_keys_requested])))
            additional_data_keys = additional_data_keys_requested

        additional_data_item_names = database.get_single_column_from_table(self.table_name,
                                                                           'item_name',
                                                                           unique=True,
                                                                           where_clause="""data_group LIKE '%s'""" % self.target_data_group)

        database.disconnect()

        if not len(additional_data_item_names):
            self.progress.end()
            return [], {}

        d = {}
        for additional_data_item_name in additional_data_item_names:
            d[additional_data_item_name] = {}

        for entry in additional_data.values():
            additional_data_item_name = entry['item_name']
            key = entry['data_key']
            value = entry['data_value']

            if entry['data_type'] in ['int', 'float']:
                d[additional_data_item_name][key] = eval(entry['data_type'])(value or self.nulls_per_type[entry['data_type']])
            else:
                d[additional_data_item_name][key] = value

        for additional_data_item_name in d:
            for key in additional_data_keys:
                if key not in d[additional_data_item_name]:
                    d[additional_data_item_name][key] = None

        self.progress.end()

        return additional_data_keys, d


    def add(self, data_dict, data_keys_list, skip_check_names=False):
        """Function to add data into the item additional data table.

           * `data_dict`: a dictionary for items or layers additional should follow this format:

                d = {
                        'item_or_layer_name_01': {'data_key_01': value,
                                                  'data_key_02': value,
                                                  'data_key_03': value
                                                  },
                        'item_or_layer_name_02': {'data_key_01': value,
                                                  'data_key_03': value,
                                                  },
                        (...)
                    }

           * `data_keys_list`: is a list of keys one or more of which should appear for each item
                               in `data_dict`.
        """

        if self.target_table not in ['items', 'layers']:
            raise ConfigError("You are using an AdditionalDataBaseClass instance to add %s data into your %s database. But\
                               you know what? You can't do that :/ Someone made a mistake somewhere. If you are a user,\
                               check your flags to make sure you are targeting the right data table. If you are a programmer,\
                               you are fired." % (self.target_table, self.db_type))

        self.data_dict_sanity_check(data_dict, data_keys_list=data_keys_list)

        self.run.warning(None, "New data for '%s' in data group '%s'" % (self.target_table, self.target_data_group), lc="yellow")
        key_types = {}
        for key in data_keys_list:
            if '!' in key:
                predicted_key_type = "stackedbar"
            else:
                type_class = utils.get_predicted_type_of_items_in_a_dict(data_dict, key)
                predicted_key_type = type_class.__name__ if type_class else None

            key_types[key] = predicted_key_type
            self.run.info('Data key "%s"' % key, 'Predicted type: %s' % (key_types[key]), \
                                            nl_after = 1 if key == data_keys_list[-1] else 0)

        # we be responsible here.
        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))
        all_keys_for_group = database.get_single_column_from_table(self.table_name, 
            'data_key', unique=True, where_clause="""data_group='%s'""" % self.target_data_group)
        database.disconnect()

        keys_already_in_db = [c for c in data_keys_list if c in all_keys_for_group]
        if len(keys_already_in_db):
            if self.just_do_it:
                self.run.warning('The following keys in your data dict will replace the ones that are already\
                                  in your %s database %s table and %s data group: %s.' \
                                                % (self.db_type, self.target_table,
                                                   self.target_data_group, ', '.join(keys_already_in_db)))

                self.remove(keys_already_in_db)
            else:
                self.run.info('Database', self.db_type, nl_before=1)
                self.run.info('Data group', self.target_data_group)
                self.run.info('Data table', self.target_table)
                self.run.info('Data keys already in db', ', '.join(keys_already_in_db), mc='red')

                raise ConfigError("Some of the data keys in your new data are already in the database. If you\
                                   want to replace them with the ones in your new data input use the `--just-do-it` flag,\
                                   and watch anvi'o make an exception just for you and complain about nothin' for this\
                                   once.")

        if skip_check_names:
            self.run.warning("You (or the programmer) asked anvi'o to NOT check the consistency of the names of your %s\
                              between your additional data and the %s database you are attempting to update. So be it.\
                              Anvi'o will not check anything, but if things don't look the way you expected them to look,\
                              you will not blame anvi'o for your poorly prepared data, but choose between yourself or\
                              Obama." % (self.target_table, self.db_type))
        else:
            if self.target_table == 'layers':
                TableForLayerAdditionalData.check_names(self, data_dict)
            elif self.target_table == 'items':
                TableForItemAdditionalData.check_names(self, data_dict)
            else:
                raise ConfigError("Congratulations, you managed to hit an uncharted are in anvi'o. It is cerrtainly very\
                                   curious how you got here unless you are trying to implement a new functionality. Are\
                                   you? What *IS* it? IS IT FUN?")

        db_entries = []
        self.set_next_available_id(self.table_name)
        for item_name in data_dict:
            for key in data_keys_list:
                db_entries.append(tuple([self.next_id(self.table_name),
                                         item_name,
                                         key,
                                         data_dict[item_name][key],
                                         key_types[key],
                                         self.target_data_group]))

        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))
        database._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?)''' % self.table_name, db_entries)
        database.disconnect()


        self.run.warning('', 'NEW DATA', lc='green')
        self.run.info('Database', self.db_type)
        self.run.info('Data group', self.target_data_group)
        self.run.info('Data table', self.target_table)
        self.run.info('New data keys', '%s.' % (', '.join(data_keys_list)), nl_after=1)


    def get_all(self):
        keys_dict, data_dict = {}, {}

        for group_name in self.available_group_names:
            self.target_data_group = group_name
            keys_dict[group_name], data_dict[group_name] = self.get()

        return keys_dict, data_dict


    def get_group_names(self):
        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))
        group_names = database.get_single_column_from_table(self.table_name, 'data_group', unique=True)
        database.disconnect()

        return group_names


class TableForItemAdditionalData(AdditionalDataBaseClass):
    """
       This is the class where we maintain the item additional data table in anvi'o
       pan and profile databases. Related issue: https://github.com/merenlab/anvio/issues/662.
    """
    def __init__(self, args, r=run, p=progress):
        self.run = r
        self.progress = p

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.table_name = A('table_name') or t.item_additional_data_table_name

        self.target_table = 'items'

        AdditionalDataBaseClass.__init__(self, args)


    def check_names(self, data_dict):
        """Compares item names found in the data dict to the ones in the db"""

        items_in_db = utils.get_all_item_names_from_the_database(self.db_path)
        items_in_data = set(data_dict.keys())

        items_in_data_but_not_in_db = items_in_data.difference(items_in_db)
        if len(items_in_data_but_not_in_db):
            raise ConfigError("Well. %d of %d item names in your additional data are only in your data (which\
                               that they are not in the %s database you are working with (which means bad news)).\
                               Since there is no reason to add additional data for items that do not exist in your\
                               database, anvi'o will stop you right there. Please fix your data and come again. In\
                               case you want to see a random item that is only in your data, here is one: %s. Stuff\
                               in your db looks like this: %s." \
                                    % (len(items_in_data_but_not_in_db), len(items_in_data), self.db_type, \
                                       items_in_data_but_not_in_db.pop(), 
                                       items_in_db.pop() if items_in_db else "No entries found in databse"))

        items_in_db_but_not_in_data = items_in_db.difference(items_in_data)
        if len(items_in_db_but_not_in_data):
            self.run.warning("Your input contains additional data for only %d of %d total number of items in your %s\
                              database. Just wanted to make sure you know what's up, but we cool." \
                                % (len(items_in_db) - len(items_in_db_but_not_in_data), len(items_in_db), self.db_type))


class TableForLayerAdditionalData(AdditionalDataBaseClass):
    """
       This is the class where we maintain the layer additional data table in anvi'o
       pan and profile databases.

       Once upon a time there was something called an anvi'o samples database. This
       is one of two tables that made it irrelevant.

       Related issue: https://github.com/merenlab/anvio/issues/674.
    """

    def __init__(self, args, r=run, p=progress):
        self.run = r
        self.progress = p

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.table_name = A('table_name') or t.layer_additional_data_table_name

        self.target_table = 'layers'

        AdditionalDataBaseClass.__init__(self, args)


    def check_names(self, data_dict):
        """Compares layer names found in the data dict to the ones in the db"""

        layers_in_db = utils.get_all_sample_names_from_the_database(self.db_path)
        layers_in_data = set(data_dict.keys())

        layers_in_data_but_not_in_db = layers_in_data.difference(layers_in_db)
        if len(layers_in_data_but_not_in_db) and not self.just_do_it:
            raise ConfigError("Grande problemo. %d of %d layers in your additional data are *only* in your data (which\
                               means they are not in the %s database you are working with). Since there is no reason to\
                               add additional data for layers that do not exist in your database, anvi'o refuses to\
                               continue, and hopes that you will try again. In case you want to see a random layer name\
                               that is only in your data, here is one: %s. In comparison, here is a random layer name\
                               from your database: %s. If you don't want to deal with this, you could use the flag\
                               `--just-do-it`, and anvi'o would do something." \
                                    % (len(layers_in_data_but_not_in_db), len(layers_in_data), self.db_type, \
                                       layers_in_data_but_not_in_db.pop(), layers_in_db.pop()))
        elif len(layers_in_data_but_not_in_db) and self.just_do_it:
            self.run.warning("Listen up! %d of %d layers in your additional data were *only* in your data (which\
                              means they are not in the %s database you are working with). But since you asked anvi'o to\
                              keep its mouth shut, it removed the ones that were not in your database from your input\
                              data, hoping that the rest of your probably very dubious operation will go just fine :/" \
                                   % (len(layers_in_data_but_not_in_db), len(layers_in_data), self.db_type))
            for layer_name in layers_in_data_but_not_in_db:
                data_dict.pop(layer_name)

        layers_in_db_but_not_in_data = layers_in_db.difference(layers_in_data)
        if len(layers_in_db_but_not_in_data):
            self.run.warning("Your input contains additional data for only %d of %d total number of layers in your %s\
                              database. Just wanted to make sure you know what's up, but we cool." \
                                % (len(layers_in_db) - len(layers_in_db_but_not_in_data), len(layers_in_db), self.db_type))


class TableForLayerOrders(OrderDataBaseClass):
    """
       This is the class where we maintain the layer order data table in anvi'o pan and profile
       databases.
    """

    def __init__(self, args, r=run, p=progress):
        self.run = r
        self.progress = p

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.table_name = A('table_name') or t.layer_orders_table_name

        self.allowde_types = ['newick', 'basic']
        self.target_table = 'layer_orders'

        OrderDataBaseClass.__init__(self, args)


    def check_names(self, data_dict):
        """Compares layer names found in the data dict to the ones in the db"""

        layers_in_db = set(utils.get_all_sample_names_from_the_database(self.db_path))
        layers_in_data = self.get_layer_names(data_dict)

        for data_key in data_dict:
            layers_in_data_for_key = set(layers_in_data[data_key])

            if len(layers_in_data_for_key.difference(layers_in_db)) > 0:
                raise ConfigError("The incoming layer orders data for %s include layer names that do not match the ones in the database :/\
                                   Here they are: '%s'" % (data_key, ', '.join(list(layers_in_data_for_key.difference(layers_in_db)))))

            if len(layers_in_data_for_key) != len(layers_in_db):
                self.run.warning("The incoming layer orders data for %s of type %s include layer names your %s database does not know\
                                  about :/ Anvi'o will let you get away with a warning, but if things go south later on, you know who\
                                  to blame. For your records, here are the layer names for %s: '%s'. And in contrast, here are the\
                                  layer names your db recognizes: '%s'. :/" % (data_key,
                                                                               data_dict[data_key]['data_type'],
                                                                               self.db_type,
                                                                               data_key,
                                                                               ', '.join(sorted(layers_in_data_for_key)),
                                                                               ', '.join(layers_in_db)))


class MiscDataTableFactory(TableForItemAdditionalData, TableForLayerAdditionalData, TableForLayerOrders):
    """Gives seamless access to additional data or order tables in pan or profile databases.

       Create an instance with args.target_data_table = [items|layers|layer_orders], and you will be golden.
    """

    def __init__(self, args, r=run, p=progress):
        self.run = r
        self.progress = p

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        target_data_table = A('target_data_table')

        if not target_data_table:
            raise ConfigError("When creating an instance from the MiscDataTableFactory class, the `args` object\
                               must contain the `target_data_table` variable.")

        if target_data_table == 'items':
            TableForItemAdditionalData.__init__(self, args, r=self.run, p=self.progress)
        elif target_data_table == 'layers':
            TableForLayerAdditionalData.__init__(self, args, r=self.run, p=self.progress)
        elif target_data_table == 'layer_orders':
            TableForLayerOrders.__init__(self, args, r=self.run, p=self.progress)
        else:
            raise ConfigError("MiscDataTableFactory does not know about target data tables for '%s' :(\
                               You can go to the online documentation, or you can try either 'items', 'layers',\
                               or 'layer_orders'" % target_data_table)


