# -*- coding: utf-8
# pylint: disable=line-too-long

"""The fancy additioanl data module"""

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, GenesDBError
from anvio.tables.tableops import Table

import pandas as pd
from anvio.dbinfo import (
    is_blank_profile,
    is_contigs_db,
    is_pan_or_profile_db,
    is_profile_db_and_contigs_db_compatible
)
from anvio.utils.anviohelp import get_contig_name_to_splits_dict
from anvio.utils.database import (
    get_all_item_names_from_the_database,
    get_all_sample_names_from_the_database,
    get_db_type,
    get_genes_database_path_for_bin,
    get_required_version_for_db
)
from anvio.utils.files import get_TAB_delimited_file_as_dictionary, get_columns_of_TAB_delim_file, store_dict_as_TAB_delimited_file
from anvio.utils.misc import get_predicted_type_of_items_in_a_dict
from anvio.utils.phylogenetics import get_names_order_from_newick_tree
from anvio.utils.validation import check_misc_data_keys_for_format


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


pp = terminal.pretty_print


class AdditionalAndOrderDataBaseClass(Table, object):
    """This is a base class for common operations between order and additional data classes."""

    def __init__(self, args):
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None

        acceptable_db_inputs = ('pan_or_profile_db', 'profile_db', 'pan_db', 'contigs_db', 'genes_db')

        # FIXME This is a first-come, first-serve business. As of now, it is the programmer's
        # responsibility to ensure they are not passing more inputs than they need.
        for db_input in acceptable_db_inputs:
            self.db_path = A(db_input)

            if self.db_path is not None:
                break

        # We just set the path for the database we are going to be working with. but if we seem to
        # be in 'gene mode', then the actual database we want to work with through this module is
        # in fact the genes database, so here we will do quite a sketchy thing, and will update our
        # db path with that.
        if A('gene_mode'):
            gene_db_path = get_genes_database_path_for_bin(profile_db_path=A('profile_db'),
                                                                 collection_name=A('collection_name'),
                                                                 bin_name=A('bin_id'))

            if filesnpaths.is_file_exists(gene_db_path, dont_raise=True):
                self.db_path = gene_db_path
            else:
                raise GenesDBError("The misc data module can't find a genes databse :( It is inherited with an args "
                                   "object that wanted to initiate anvi'o operations for 'gene mode', which is a special "
                                   "mode of operation where gene-level coverage statistics per collection is "
                                   "read from a special database. In this mode anvi'o also tries to initialize "
                                   "additional data tables from the genes database, instead of the profile "
                                   "database with which it is associated. However, in the current run, it seems the "
                                   "genes database has not yet been initiated for the collection '%s' and bin '%s'. "
                                   "Probably this will be handled by a higher power, and the genes database will "
                                   "be generated automatically, but this very part of the code has no idea how to "
                                   "deal with this awkward situation, hence throwing this exception and waves its hand "
                                   "to you from a wild wild corner of the anvi'o codebase." % (A('collection_name'), A('bin_id')))

        self.just_do_it = A('just_do_it')
        self.target_data_group_set_by_user = A('target_data_group') or None
        self.skip_check_names = A('skip_check_names')
        self.target_data_group = self.target_data_group_set_by_user or 'default'

        # optional arguments to keep track of
        self.contigs_mode = A('contigs_mode')
        self.contigs_db_path = A('contigs_db')

        if not self.db_path:
            raise ConfigError("The AdditionalAndOrderDataBaseClass is inherited with an args object that did not "
                              "contain any database :/ Even though any of the following argument names would "
                              "have worked: %s :(" % ', '.join(["'%s'" % x for x in acceptable_db_inputs]))

        if not self.table_name:
            raise ConfigError("The AdditionalAndOrderDataBaseClass does not know anything about the table it should "
                              "be working with.")

        self.db_type = get_db_type(self.db_path)

        if self.db_type in ['pan', 'profile', 'genes']:
            is_pan_or_profile_db(self.db_path, genes_db_is_also_accepted=True)

            if self.target_table not in ['layers', 'items', 'layer_orders']:
                raise ConfigError("The only target data table names possible for DBs of type '%s' are 'layers' "
                                  "and 'items'" % self.db_type)

        elif self.db_type == 'contigs':
            is_contigs_db(self.db_path)

            if self.target_table not in ['nucleotides', 'amino_acids']:
                raise ConfigError("The only target data table names possible for DBs of type '%s' are 'nucleotides' "
                                  "and 'amino_acids'" % self.db_type)

        self.db_version = get_required_version_for_db(self.db_path)

        self.progress.new("Fetching existing data keys from database")
        self.progress.update("...")

        database = db.DB(self.db_path, self.db_version)
        self.additional_data_keys = database.get_single_column_from_table(self.table_name, 'data_key')
        database.disconnect()

        self.progress.end()

        Table.__init__(self, self.db_path, self.db_version, self.run, self.progress)

        self.nulls_per_type = {'str': '',
                               'int': 0,
                               'float': 0,
                               'stackedbar': 0,
                               'unknown': None}

        self.df = None


    def populate_from_file(self, additional_data_file_path, skip_check_names=None):

        if skip_check_names is None and is_blank_profile(self.db_path):
            # FIXME: this BS is here because blank anvi'o profiles do not know what items they have,
            #        hence the get_all_item_names_from_the_database function eventually explodes if we
            #        don't skip check names.
            skip_check_names = True

        self.progress.new("Processing input")

        self.progress.update("Checking integrity of input file")
        filesnpaths.is_file_tab_delimited(additional_data_file_path)

        self.progress.update("Loading input file into memory")
        data_keys = get_columns_of_TAB_delim_file(additional_data_file_path)
        data_dict = get_TAB_delimited_file_as_dictionary(additional_data_file_path)

        bad_keys = [k for k in data_keys if k.strip() != k]
        if len(bad_keys):
            raise ConfigError("Some of the keys in your input file contain white characters at the beginning "
                              "or at the end. This is your lucky day, because anvi'o caught it, and you will "
                              "fix it before continuing, and those weird data keys will not cause any headaches "
                              "in the future. Thank you and you are welcome. Here are the offending keys: %s." % \
                                        (', '.join(['"%s"' % k for k in bad_keys])))

        if not len(data_keys):
            raise ConfigError("There is something wrong with the additional data file for %s at %s. "
                              "It does not seem to have any additional keys for data :/" \
                                            % (self.target_table, additional_data_file_path))

        self.progress.end()
        self.run.info_single("%s successfully loaded" % additional_data_file_path, nl_after=1, nl_before=1)

        if self.target_table == 'layer_orders':
            OrderDataBaseClass.add(self, data_dict, skip_check_names)
        else:
            AdditionalDataBaseClass.add(self, data_dict, data_keys, skip_check_names)


    def remove(self, data_keys_list=[], data_groups_list=[]):
        """Give this guy a list of keys or groups for additional data, and watch their demise."""

        if not isinstance(data_keys_list, list) or not isinstance(data_groups_list, list):
            raise ConfigError("The remove function in AdditionalDataBaseClass wants you to watch "
                              "yourself before you wreck yourself. In other words, can you please "
                              "make sure the keys or groups you send is of type `list` thankyouverymuch?")

        if data_keys_list and data_groups_list:
            raise ConfigError("You seem to be interested in removing both data keys and data groups from "
                              "misc data tables. Using both of these parameters is quite risky, so anvi'o "
                              "would like you to use only one of them at a time. Apologies.")

        database = db.DB(self.db_path, get_required_version_for_db(self.db_path))

        # not all additional data table types have data_group concept
        if data_groups_list:
            if 'data_group' not in database.get_table_structure(self.table_name):
                raise ConfigError("You seem to be interested in removing some data groups from your profile database "
                                  "but the additional data table for %s does not have data groups :/ Perhaps you are "
                                  "looking for the parameter `--keys-to-remove`?" % self.target_table)
            else:
                additional_data_groups = sorted(database.get_single_column_from_table(self.table_name, 'data_group', unique=True))

        additional_data_keys = sorted(database.get_single_column_from_table(self.table_name, 'data_key', unique=True))

        if not len(additional_data_keys):
            self.run.info_single('There is nothing to remove--the %s additional data table is already empty :(' % self.target_table)
            database.disconnect()
            return

        if data_keys_list:
            # see if there are any missing keys
            missing_keys = [k for k in data_keys_list if k not in additional_data_keys]
            if len(missing_keys) and not self.just_do_it:
                database.disconnect()
                raise ConfigError("The following keys you wanted to remove from the %s additional data table are "
                                  "not really in the table: '%s'. Anvi'o is confused but can totally ignore the missing "
                                  "keys and just remove whatever is matching to your list if you were to use the flag "
                                  "`--just-do-it`." % (self.target_table, ', '.join(missing_keys)))

            for key in data_keys_list:
                if key not in additional_data_keys:
                    # what the hell, user?
                    continue

                if 'data_group' in database.get_table_structure(self.table_name):
                    database._exec('''DELETE from %s WHERE data_key="%s" and data_group="%s"''' % (self.table_name, key, self.target_data_group))
                else:
                    database._exec('''DELETE from %s WHERE data_key="%s"''' % (self.table_name, key))

            self.run.warning("Data from the table '%s' for the following data keys in data group '%s' "
                             "removed from the database: '%s'. #SAD." % (self.target_table, self.target_data_group, ', '.join(data_keys_list)))
        elif data_groups_list:
            # see if there are any missing groups
            missing_groups = [k for k in data_groups_list if k not in additional_data_groups]
            if len(missing_groups) and not self.just_do_it:
                database.disconnect()
                raise ConfigError("The following groups you wanted to remove from the %s additional data table are "
                                  "not really in the table: '%s'. Anvi'o is confused but can totally ignore the missing "
                                  "keys and just remove whatever is matching to your list if you were to use the flag "
                                  "`--just-do-it`." % (self.target_table, ', '.join(missing_groups)))

            for group in data_groups_list:
                if group not in additional_data_groups:
                    # user is cray.
                    continue
                else:
                    database._exec('''DELETE from %s WHERE data_group="%s"''' % (self.table_name, group))

            self.run.warning("Data from the table '%s' for the following data groups "
                             "are now removed from the database: '%s'. #SAD." % (self.target_table, ', '.join(data_groups_list)))
        else:
            if not self.just_do_it:
                raise ConfigError("You did not provide a list of data keys to remove, which means you are about to delete everything in the "
                                  "%s additional data table. Just to be on the safe side, anvi'o is looking for a confirmation. If you "
                                  "try again with the --just-do-it flag, anvi'o will put on its business socks, and burn this table "
                                  "and everything in it to the ground." % self.target_table)

            database._exec('''DELETE from %s''' % (self.table_name))

            self.run.warning("All data from the %s additional data table is removed (ouch)." % self.target_table)

        database.disconnect()


    def export(self, output_file_path):
        filesnpaths.is_output_file_writable(output_file_path)

        if self.target_table in ['layers', 'items', 'nucleotides', 'amino_acids']:
            keys, data = AdditionalDataBaseClass.get(self)
            if keys:
                if len(self.available_group_names) - 1:
                    self.run.warning("You are exporting data from the additional data table '%s' for the "
                                     "data group '%s'. Great. Just remember that there are %d more data "
                                     "groups in your database, and you are not exporting anything from them "
                                     "at this point (they know you're the boss, so they're not upset)." \
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

        store_dict_as_TAB_delimited_file(data, output_file_path, headers=[self.target_table] + keys)

        self.run.info('Target data table', self.target_table)
        self.run.info('Output file', output_file_path)


    def list_data_keys(self):
        database = db.DB(self.db_path, get_required_version_for_db(self.db_path))

        NOPE = lambda: self.run.info_single("There are no additional data for '%s' in this database :/" \
                                                % (self.target_table), nl_before=1, nl_after=1, mc='red')

        additional_data_keys = {}
        # here is where things get tricky. if we are dealing with additional data layers or items, we will have
        # data groups that are not relevant for order data. this will affect the listing of data keys in either
        # of these table types. hence we get group names first here, and then will do a bunch of if/else checks
        # based on their availability
        if self.target_table in ['layers', 'items', 'nucleotides', 'amino_acids']:
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

        check_misc_data_keys_for_format(data_keys_list)

        # FIXME: we have two controls here. The first one is how we work with order data natively. The second one is how it
        #        looks like when it is read through the .get() member function of the TableForLayerOrders because rest of
        #        anvi'o does not know how to work with the native format. This should be fixed by teaching the rest of anvi'o
        #        how to work with order data dicts in the native form.
        looks_like_layer_orders = sorted(list(data_dict.values())[0].keys()) == sorted(['data_type', 'data_value']) or \
                                  sorted(list(data_dict.values())[0].keys()) == sorted(['basic', 'newick'])

        if looks_like_layer_orders and data_dict_type != 'layer_orders':
            raise ConfigError("The data you sent here seems to describe an order, but you want anvi'o to treat it "
                              "as additional data for %s. Not cool." % self.target_table)

        if not looks_like_layer_orders and data_dict_type == 'layer_orders':
            raise ConfigError("The data that claims to be a layer order data do not seem to be one.")

        if data_keys_list:
            # data_name is an item name, layer name, nucleotide identifier, or amino acid identifier
            for data_name in data_dict:
                for key in data_keys_list:
                    if key not in data_dict[data_name]:
                        raise ConfigError("Your data dictionary fails the sanity check since at least one item in it (i.e., %s) is "
                                          "missing any data for the key '%s'." % (data_name, key))


    def init_table_as_dataframe(self):
        """Init the raw table as a dataframe"""

        database = db.DB(self.db_path, get_required_version_for_db(self.db_path))
        self.df = database.get_table_as_dataframe(self.table_name, error_if_no_data=False)
        database.disconnect()


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
             >>> layer_orders_data_dict = layer_orders_table.get(layer_additional_data_keys, layer_additional_data_dict)

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
        database = db.DB(self.db_path, get_required_version_for_db(self.db_path))
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
                raise ConfigError("Something is wrong :( Anvi'o just found an entry in %s table with an unrecognized "
                                  "type of '%s'." % (self.target_table, data_type))

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
                type_class = get_predicted_type_of_items_in_a_dict(additional_data_dict, data_key)
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
        """The default function to add data into the orders table.

        Parameters
        ==========
        data_dict : dict
            This variable for layer orders is expected to follow this format:
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
            raise ConfigError("You are using an OrderDataBaseClass instance to add %s data into your %s database. This is "
                              "illegal and if you are here, it means someone made a mistake somewhere. If you are a user, "
                              "check your flags to make sure you are targeting the right data table. If you are a programmer, "
                              "you are fired." % (self.target_table, self.db_type))

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
                self.run.warning('The following keys in your data dict will replace the ones that are already '
                                 'in your %s database: %s.' % (self.db_type, ', '.join(keys_already_in_db)))

                self.remove(keys_already_in_db)
            else:
                self.run.info('Data keys already in the db', ', '.join(keys_already_in_db), nl_before=2, mc='red')

                raise ConfigError("Some of the keys in your new order data appear to be in the database already. If you "
                                  "want to replace those in the database with the ones in your new data use the "
                                  "`--just-do-it` flag.")

        if skip_check_names:
            self.run.warning("You (or the programmer) asked anvi'o to NOT check the consistency of the names of your %s "
                             "between your additional data and the %s database you are attempting to update. So be it. "
                             "Anvi'o will not check anything, but if things don't look the way you expected them to look, "
                             "you will not blame anvi'o for your poorly prepared data, but choose between yourself or "
                             "Obama." % (self.target_table, self.db_type))
        else:
            TableForLayerOrders.check_names(self, data_dict)

        db_entries = []
        for item_name in data_dict:
            db_entries.append(tuple([item_name,
                                     data_dict[item_name]['data_type'],
                                     data_dict[item_name]['data_value']]))

        database = db.DB(self.db_path, get_required_version_for_db(self.db_path))
        database._exec_many('''INSERT INTO %s VALUES (?,?,?)''' % self.table_name, db_entries)
        database.disconnect()

        self.run.info('New order data added to the db for %s' % self.target_table, '%s.' % (', '.join(data_keys_list)))


    def get_layer_names(self, data_dict):
        layer_names = {}

        for data_key in data_dict:
            try:
                if data_dict[data_key]['data_type'] == 'newick':
                    layer_names[data_key] = get_names_order_from_newick_tree(data_dict[data_key]['data_value'])
                else:
                    layer_names[data_key] = [s.strip() for s in data_dict[data_key]['data_value'].split(',')]
            except:
                raise ConfigError("Parsing the %s data for %s failed :/ We don't know why, because we are lazy. Please "
                                  "take a loook at your input data and figure out :(" \
                                                                    % (data_dict[data_key]['data_type'], data_key))

        return layer_names


class AdditionalDataBaseClass(AdditionalAndOrderDataBaseClass, object):
    """Implements additional data ops base class.

    Notes
    =====
    - See TableForItemAdditionalData or TableForLayerAdditionalData for usage example.
    - See AdditionalAndOrderDataBaseClass for inherited functionality.
    """

    def __init__(self, args):
        AdditionalAndOrderDataBaseClass.__init__(self, args)

        self.available_group_names = self.get_group_names()

        self.storage_buffer = []


    def check_target_data_group(self):
        """A function to check whether the data group set is among the available ones.

        The reason this function is a separate one and it is not being called by the
        init function of the base class is because the user *can* set a new data group
        name when they want to 'import' things into the database. So, while exporting
        data or displaying data will require the requested data group to be a proper one,
        it is better to check that explicitly.
        """

        if self.target_data_group not in self.available_group_names:
            raise ConfigError("You (or the programmer) requested to initiate the additional data table for '%s' with "
                              "the data group '%s', which is not really in that table :/ If it helps at all, "
                              "the target table happened to have these ones instead: %s. What to do now? If you are "
                              "here because the last command you run was something like 'show me all the data in misc' "
                              "data tables, then you may try to be more specific by explicitly defining your target "
                              "data table. If you think you have already been as sepecific as you could be, then anvi'o "
                              "is as frustrated as you are right now :(" % \
                                    (self.target_table, self.target_data_group, ', '.join(['"%s"' % d for d in self.available_group_names])))


    def get_available_data_keys(self):
        """Will only return the additional data keys so the client can do some controls."""

        if not self.target_data_group:
            raise ConfigError("The target data group is not set. Which should never be the case at this stage. You shall "
                              "not break anvi'o and go back to where you came from, devil :(")

        self.progress.new('Recovering additional keys and data for %s' % self.target_table)
        self.progress.update('...')

        database = db.DB(self.db_path, get_required_version_for_db(self.db_path))
        additional_data_keys_in_db = database.get_single_column_from_table(self.table_name, 'data_key', unique=True, \
                        where_clause="""data_group LIKE '%s'""" % self.target_data_group)
        database.disconnect()

        self.progress.end()

        return additional_data_keys_in_db


    def get(self, additional_data_keys_requested=[]):
        """Get the additional data

        Parameters
        ==========
        additional_data_keys_requested : list, []
            Which keys are requested? If [], all are assumed. Read as "data_keys_requested", they
            "additional" because they are fetched from 'layer_additional_data', 'item_additional_data',
            'amino_acid_additional_data', or 'nucleotide_additional_data' tables.

        Returns
        =======
        (additional_data_keys, additional_data_dict) : tuple
        """

        if not self.target_data_group:
            raise ConfigError("It seems the target data group is not set, which makes zero sense and should never happen "
                              "unless you are doing some hacker stuff. Are you doing hacker stuff? Awesome! Tell us about "
                              "it!")

        if not isinstance(additional_data_keys_requested, list):
            raise ConfigError("The `get` function in AdditionalDataBaseClass is upset with you. You could change that "
                              "by making sure you request additional data keys with a variable of type `list`.")

        additional_data_keys_in_db = self.get_available_data_keys()

        self.progress.new('Recovering additional data for %s' % self.target_table)
        self.progress.update('...')

        database = db.DB(self.db_path, get_required_version_for_db(self.db_path))
        if not len(additional_data_keys_requested):
            additional_data_keys = additional_data_keys_in_db
            additional_data = database.get_some_rows_from_table_as_dict(self.table_name,
                                                                        where_clause = """data_group LIKE '%s'""" % self.target_data_group,
                                                                        error_if_no_data=False)
        else:
            if not len(additional_data_keys_in_db):
                raise ConfigError("The %s database at %s does not contain any additional data for its %s to return. Usually this "
                                  "would not have resulted in an exception, and anvi'o would simply returned an empty data "
                                  "dictionary. But since you are requested to get a specific list of keys ('%s'), we are raising "
                                  "an exception and break the flow just to make sure you are aware of the fact that we don't "
                                  "see the stuff you're requesting :(" \
                                        % (self.db_type, self.db_path, self.target_table, ' ,'.join(additional_data_keys_requested)))

            if set(additional_data_keys_requested) - set(additional_data_keys_in_db):
                raise ConfigError("The keys you requested does not seem to appear in the additional data table of this %s db "
                                  "at '%s' :/ Here is the list of keys you requested: '%s'. And here is the list of keys that anvi'o "
                                  "knows about: '%s'." % (self.db_type, self.db_path,
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

            if entry['data_type'] == 'int':
                # our predictor predicts 'int' if all values are int-ish floats. Ex: [37.0, 36.0, 1.0, ...]
                # but python int() function fails to parse strings like '1.0' while float() can parse them correctly.
                # so here we first parse them with float() and then pass them to int(). Which should not affect anything in theory.
                d[additional_data_item_name][key] = int(float(value or self.nulls_per_type[entry['data_type']]))
            elif entry['data_type'] == 'float':
                d[additional_data_item_name][key] = float(value or self.nulls_per_type[entry['data_type']])
            else:
                d[additional_data_item_name][key] = value

        for additional_data_item_name in d:
            for key in additional_data_keys:
                if key not in d[additional_data_item_name]:
                    d[additional_data_item_name][key] = None

        self.progress.end()

        return additional_data_keys, d


    def update_input_data_dict_with_split_names(self, data_dict):
        """Takes data dict with contig names, turns it into data dict with split names"""

        if not self.contigs_db_path:
            raise ConfigError("Anvi'o was about to attempt to translate item names in the additonal data to "
                              "split names (probably due to the `--contigs-mode`), yet, the environment does not "
                              "seem to know about a contigs database. Did you forget to provide a contigs-db? "
                              "Because it is kind of necessary :/")

        is_profile_db_and_contigs_db_compatible(self.db_path, self.contigs_db_path)

        contig_name_to_split_names_dict = get_contig_name_to_splits_dict(self.contigs_db_path)

        new_data_dict = {}
        for contig_name in data_dict:
            for split_name in contig_name_to_split_names_dict[contig_name]:
                new_data_dict[split_name] = data_dict[contig_name]

        return new_data_dict


    def add(self, data_dict, data_keys_list, skip_check_names=False):
        """Function to add data into the item additional data table.

        Parameters
        ==========
        data_dict : dict
            A dictionary for items or layers additional should follow this format:

                data_dict = {
                        'item_or_layer_name_01': {'data_key_01': value,
                                                  'data_key_02': value,
                                                  'data_key_03': value
                                                  },
                        'item_or_layer_name_02': {'data_key_01': value,
                                                  'data_key_03': value,
                                                  },
                        (...)
                    }

        data_keys_list : list
            A list of keys one or more of which should appear for each item in `data_dict`.
        """

        if self.target_table not in ['items', 'layers', 'nucleotides', 'amino_acids']:
            raise ConfigError("You are using an AdditionalDataBaseClass instance to add %s data into your %s database. But "
                              "you know what? You can't do that :/ Someone made a mistake somewhere. If you are a user, "
                              "check your flags to make sure you are targeting the right data table. If you are a programmer, "
                              "you are fired." % (self.target_table, self.db_type))

        # here we will check for a name conversion request between contig names to split names.
        # this is only necessery if we are adding `items` data and if we are in contigs mode.
        if self.target_table == 'items' and self.contigs_mode:
            data_dict = self.update_input_data_dict_with_split_names(data_dict)

        self.data_dict_sanity_check(data_dict, data_keys_list=data_keys_list)

        self.run.warning(None, "New data for '%s' in data group '%s'" % (self.target_table, self.target_data_group), lc="yellow")
        key_types = {}
        for key in data_keys_list:
            if '!' in key:
                predicted_key_type = "stackedbar"
            else:
                type_class = get_predicted_type_of_items_in_a_dict(data_dict, key)
                predicted_key_type = type_class.__name__ if type_class else None

            key_types[key] = predicted_key_type
            self.run.info('Data key "%s"' % key, 'Predicted type: %s' % (key_types[key]), \
                                            nl_after = 1 if key == data_keys_list[-1] else 0)

        # we be responsible here.
        database = db.DB(self.db_path, get_required_version_for_db(self.db_path))
        all_keys_for_group = database.get_single_column_from_table(self.table_name,
            'data_key', unique=True, where_clause="""data_group='%s'""" % self.target_data_group)
        database.disconnect()

        keys_already_in_db = [c for c in data_keys_list if c in all_keys_for_group]
        if len(keys_already_in_db):
            if self.just_do_it:
                self.run.warning('The following keys in your data dict will replace the ones that are already '
                                 'in your %s database %s table and %s data group: %s.' \
                                                % (self.db_type, self.target_table,
                                                   self.target_data_group, ', '.join(keys_already_in_db)))

                self.remove(keys_already_in_db)
            else:
                self.run.info('Database', self.db_type, nl_before=1)
                self.run.info('Data group', self.target_data_group)
                self.run.info('Data table', self.target_table)
                self.run.info('Data keys already in db', ', '.join(keys_already_in_db), mc='red')

                raise ConfigError("Some of the data keys in your new data are already in the database. If you "
                                  "want to replace them with the ones in your new data input use the `--just-do-it` flag, "
                                  "and watch anvi'o make an exception just for you and complain about nothin' for this "
                                  "once.")

        if skip_check_names or self.skip_check_names:
            self.run.warning("You (or the programmer) asked anvi'o to NOT check the consistency of the names of your %s "
                             "between your additional data and the %s database you are attempting to update. So be it. "
                             "Anvi'o will not check anything, but if things don't look the way you expected them to look, "
                             "you will not blame anvi'o for your poorly prepared data, but choose between yourself or "
                             "Obama." % (self.target_table, self.db_type))
        else:
            if self.target_table not in table_classes:
                raise ConfigError("Congratulations, you managed to hit an uncharted are in anvi'o. It is certainly very "
                                  "curious how you got here unless you are trying to implement a new functionality. Are "
                                  "you? What *IS* it? IS IT FUN?")

            table_classes[self.target_table].check_names(self, data_dict)

        num_entries = len(data_dict)

        self.progress.new("Adding data to DB", progress_total_items=num_entries)

        for i, item_name in enumerate(data_dict):
            if (i % 100000) == 0:
                self.progress.increment(increment_to=i)
                self.progress.update('%d / %d Rows ⚙  | Writing to DB 💾 ...' % (i, num_entries))

                self.store_buffer()
                self.storage_buffer = []

                self.progress.update('%d / %d Rows ⚙  ...' % (i, num_entries))

            for key in data_keys_list:
                self.storage_buffer.append(tuple([item_name,
                                                  key,
                                                  data_dict[item_name][key],
                                                  key_types[key],
                                                  self.target_data_group]))

        self.progress.increment(increment_to=num_entries)
        self.progress.update('%d / %d Rows ⚙  | Writing to DB 💾 ...' % (num_entries, num_entries))

        self.store_buffer()
        self.storage_buffer = []

        self.progress.end()

        self.run.warning('', 'NEW DATA', lc='green')
        self.run.info('Database', self.db_type)
        self.run.info('Data group', self.target_data_group)
        self.run.info('Data table', self.target_table)
        self.run.info('New data keys', '%s.' % (', '.join(data_keys_list)), nl_after=1)


    def store_buffer(self):
        if not len(self.storage_buffer):
            return

        database = db.DB(self.db_path, get_required_version_for_db(self.db_path))
        database._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?)''' % self.table_name, self.storage_buffer)
        database.disconnect()


    def get_all(self):
        keys_dict, data_dict = {}, {}

        for group_name in self.available_group_names:
            self.target_data_group = group_name
            keys_dict[group_name], data_dict[group_name] = self.get()

        return keys_dict, data_dict


    def get_group_names(self):
        database = db.DB(self.db_path, get_required_version_for_db(self.db_path))
        group_names = database.get_single_column_from_table(self.table_name, 'data_group', unique=True)
        database.disconnect()

        return group_names


class TableForItemAdditionalData(AdditionalDataBaseClass):
    """Maintains the item additional data table in anvi'o pan and profile databases.

    Notes
    =====
    - Related issue: https://github.com/merenlab/anvio/issues/662.
    """

    def __init__(self, args, r=terminal.Run(), p=terminal.Progress()):
        self.run = r
        self.progress = p

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.table_name = A('table_name') or t.item_additional_data_table_name

        self.target_table = 'items'

        AdditionalDataBaseClass.__init__(self, args)


    def check_names(self, data_dict):
        """Compares item names found in the data dict to the ones in the db"""

        items_in_db = get_all_item_names_from_the_database(self.db_path)
        items_in_data = set(data_dict.keys())

        # this step is essential for genes database, where item naems are actually
        # gene caller ids that are meant to be integers:
        try:
            items_in_db = set([int(i) for i in items_in_db])
            items_in_data = set([int(i) for i in data_dict.keys()])
        except:
            pass

        # nothing in db?
        if not len(items_in_db):
            raise ConfigError("No item names found in your target database. This is normal if you are working with "
                              "a blank profile database. But in that case, you need to include `skip_check_names=True` "
                              "in your `args` instance. This message is meant for a developer since this part of the "
                              "code only accessible thorugh API calls, if you are a user and seeing this please get in "
                              "touch with anvi'o developers.")

        # there are items in data but not in db?
        items_in_data_but_not_in_db = items_in_data.difference(items_in_db)
        if len(items_in_data_but_not_in_db):
            if len(items_in_data_but_not_in_db) == len(items_in_data):
                items_desc_text = (f"None of the {pp(len(items_in_data))} item names in your additional data are matching "
                                   f"to those that are in the {self.db_type} you are working with")
            else:
                items_desc_text = (f"{pp(len(items_in_data_but_not_in_db))} of {pp(len(items_in_data))} entries in your "
                                   f"additional data seem to be ONLY in your additional data, and not in your database")

            raise ConfigError(f"Well. {items_desc_text} :/ Since there is no reason to add additional data for items "
                              f"that do not exist in your database, anvi'o will stop you right there. Please fix your "
                              f"data and come again. In case you want to see a random item that is only in your data, "
                              f"here is one: '{items_in_data_but_not_in_db.pop()}'. Item names in your db generally "
                              f"looks like this: '{items_in_db.pop()}'. If you are working with contigs names instead "
                              f"of split names, perhaps you need the `--contigs-mode` flag?")

        # there are items in db but not in data?
        items_in_db_but_not_in_data = items_in_db.difference(items_in_data)
        if len(items_in_db_but_not_in_data):
            self.run.warning(f"Your input contains additional data for only {pp(len(items_in_db) - len(items_in_db_but_not_in_data))} "
                             f"of {pp(len(items_in_db))} total number of items in your {self.db_type} database. "
                             f"Just wanted to make sure you know what's up, but we cool.")


class TableForLayerAdditionalData(AdditionalDataBaseClass):
    """Maintains the layer additional data table in anvi'o pan and profile databases.

    Notes
    =====
    - Once upon a time there was something called an anvi'o samples database. This
      is one of two tables that made it irrelevant. Related issue on the topic:
      https://github.com/merenlab/anvio/issues/674.
    """

    def __init__(self, args, r=terminal.Run(), p=terminal.Progress()):
        self.run = r
        self.progress = p

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.table_name = A('table_name') or t.layer_additional_data_table_name

        self.target_table = 'layers'

        AdditionalDataBaseClass.__init__(self, args)


    def check_names(self, data_dict):
        """Compares layer names found in the data dict to the ones in the db"""

        layers_in_db = get_all_sample_names_from_the_database(self.db_path)
        layers_in_data = set(data_dict.keys())

        layers_in_data_but_not_in_db = layers_in_data.difference(layers_in_db)
        if len(layers_in_data_but_not_in_db) and not self.just_do_it:
            raise ConfigError("Grande problemo. %d of %d layers in your additional data are *only* in your data (which "
                              "means they are not in the %s database you are working with). Since there is no reason to "
                              "add additional data for layers that do not exist in your database, anvi'o refuses to "
                              "continue, and hopes that you will try again. In case you want to see a random layer name "
                              "that is only in your data, here is one: %s. In comparison, here is a random layer name "
                              "from your database: %s. If you don't want to deal with this, you could use the flag "
                              "`--just-do-it`, and anvi'o would do something." \
                                    % (len(layers_in_data_but_not_in_db), len(layers_in_data), self.db_type, \
                                       layers_in_data_but_not_in_db.pop(), layers_in_db.pop()))
        elif len(layers_in_data_but_not_in_db) and self.just_do_it:
            self.run.warning("Listen up! %d of %d layers in your additional data were *only* in your data (which "
                             "means they are not in the %s database you are working with). But since you asked anvi'o to "
                             "keep its mouth shut, it removed the ones that were not in your database from your input "
                             "data, hoping that the rest of your probably very dubious operation will go just fine :/" \
                                   % (len(layers_in_data_but_not_in_db), len(layers_in_data), self.db_type))
            for layer_name in layers_in_data_but_not_in_db:
                data_dict.pop(layer_name)

        layers_in_db_but_not_in_data = layers_in_db.difference(layers_in_data)
        if len(layers_in_db_but_not_in_data):
            self.run.warning("Your input contains additional data for only %d of %d total number of layers in your %s "
                             "database. Just wanted to make sure you know what's up, but we cool." \
                                % (len(layers_in_db) - len(layers_in_db_but_not_in_data), len(layers_in_db), self.db_type))


class TableForLayerOrders(OrderDataBaseClass):
    """Maintains the layer order data table in anvi'o pan and profile databases."""

    def __init__(self, args, r=terminal.Run(), p=terminal.Progress()):
        self.run = r
        self.progress = p

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.table_name = A('table_name') or t.layer_orders_table_name

        self.allowde_types = ['newick', 'basic']
        self.target_table = 'layer_orders'

        OrderDataBaseClass.__init__(self, args)


    def check_names(self, data_dict):
        """Compares layer names found in the data dict to the ones in the db"""

        layers_in_db = set(get_all_sample_names_from_the_database(self.db_path))
        layers_in_data = self.get_layer_names(data_dict)

        for data_key in data_dict:
            layers_in_data_for_key = set(layers_in_data[data_key])

            if len(layers_in_data_for_key.difference(layers_in_db)) > 0:
                raise ConfigError("The incoming layer orders data for %s include layer names that do not match the ones in the database :/ "
                                  "Here they are: '%s'" % (data_key, ', '.join(list(layers_in_data_for_key.difference(layers_in_db)))))

            if len(layers_in_data_for_key) != len(layers_in_db):
                self.run.warning("The incoming layer orders data for %s of type %s include layer names your %s database does not know "
                                 "about :/ Anvi'o will let you get away with a warning, but if things go south later on, you know who "
                                 "to blame. For your records, here are the layer names for %s: '%s'. And in contrast, here are the "
                                 "layer names your db recognizes: '%s'. :/" % (data_key,
                                                                               data_dict[data_key]['data_type'],
                                                                               self.db_type,
                                                                               data_key,
                                                                               ', '.join(sorted(layers_in_data_for_key)),
                                                                               ', '.join(layers_in_db)))


class TableForNucleotideAdditionalData(AdditionalDataBaseClass):
    """Maintains 'nucleotide_additional_data' table in anvi'o contigs databases"""

    def __init__(self, args, r=terminal.Run(), p=terminal.Progress()):
        self.run = r
        self.progress = p

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.table_name = A('table_name') or t.nucleotide_additional_data_table_name

        self.target_table = 'nucleotides'

        AdditionalDataBaseClass.__init__(self, args)


    def check_names(self, data_dict):
        """Compares data key values found in the data dict to the ones in the db"""

        database = db.DB(self.db_path, get_required_version_for_db(self.db_path))

        # A dict of {<contig_name>: <contig_length>, ...}
        contig_lengths = dict(database.get_some_columns_from_table(t.contigs_info_table_name, 'contig,length'))

        is_format_good = lambda key: len(key.split(':')) == 2

        not_in_db = set()

        # These tables could be millions of entries. We do not bother compiling a prettified list of
        # their badly formatted table. If it's wrong we complain immediately.
        for i, data_key in enumerate(data_dict):
            if not is_format_good(data_key):
                raise ConfigError("The data key '%s' (entry #%d) is not properly formatted. It should "
                                  "have the format <contig_name>:<position_in_contig>, e.g. "
                                  "contig_000001:125 would be an entry for contig_000001 at the "
                                  "125th position of the sequence (The first letter in the sequence "
                                  "has the position 0, not 1)." % (data_key, i+1))

            contig, pos = data_key.split(':')
            pos = int(pos)

            if contig not in contig_lengths:
                not_in_db.add(data_key)

            elif pos < 0 or pos >= contig_lengths[contig]:
                raise ConfigError("The data key '%s' (entry #%d) is not properly formatted. The contig "
                                  "'%s' exists in the database, but the position %d falls outside of "
                                  "the contig, which has a length of %d." % \
                                  (data_key, i+1, contig, pos, contig_lengths[contig]))

        if len(not_in_db):
            example = next(iter(not_in_db))
            msg = ("Listen up! %d of %d entries in your additional data specified contigs that are *only* "
                   "in your data--that means these contigs don't exist in the %s database you are "
                   "working with. For example, the data key '%s' corresponds to the contig '%s', which "
                   "is non-existent in this database." % (len(not_in_db),
                                                          len(data_dict),
                                                          self.db_type,
                                                          example,
                                                          example.split(':')[0]))

            if self.just_do_it:
                self.run.warning(msg + "... But since you asked anvi'o to keep its mouth shut, it removed the ones that "
                                 "were not in your database from your input data, hoping that the rest of your "
                                 "probably very dubious operation will go just fine :/")

                for data_key in not_in_db:
                    data_dict.pop(data_key)
            else:
                raise ConfigError(msg + " If you want to proceed anyways, rerun with the --just-do-it flag.")

        database.disconnect()


class TableForAminoAcidAdditionalData(AdditionalDataBaseClass):
    """Maintains 'amino_acid_additional_data' table in anvi'o contigs databases"""

    def __init__(self, args, r=terminal.Run(), p=terminal.Progress()):
        self.run = r
        self.progress = p

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.table_name = A('table_name') or t.amino_acid_additional_data_table_name

        self.target_table = 'amino_acids'

        AdditionalDataBaseClass.__init__(self, args)


    def get_multigene_dataframe(self, gene_caller_ids=set([]), keys_of_interest=set([]), group_name=None):
        """Fetches and pivots table data specific to a group of gene_callers_ids

        Parameters
        ==========
        gene_caller_ids : set
            If empty, all are assumed.

        Returns
        =======
        output : pd.DataFrame
            Outputs a dataframe that looks like this:

                gene_callers_id  codon_order_in_gene   MG                SM_               UDP
                1                149                  NaN    0.0349008254091               NaN
                1                150                  NaN    0.0116165975232               NaN
                1                167                  NaN    0.1009224794937   0.1127902398775
                1                168                  NaN    0.0110966620296  0.03150767340385
                1                169                  NaN     0.435271159751   0.4638648995855
                1                171                  NaN    0.0239539337947  0.03150767340385
                1                172                  NaN    0.0347266126837               NaN
                1                173                  NaN     0.244954708876   0.1630916348525
                1                174                  NaN    0.5049807590385    0.604658463037
                1                175                  NaN    0.0237121180997  0.06609853613995
                1                176                  NaN    0.0359138454456  0.10012280860365
                1                177                  NaN    0.0252309121402   0.0691817254722
                1                178                  NaN   0.02856525986075  0.03150767340385
                1                179                  NaN   0.01191453799525               NaN
                ...              ...                  ...               ...                ...

            Only keys that had at least one non-NaN value for the genes are included
        """

        if self.df is None:
            self.init_table_as_dataframe()

        if not len(gene_caller_ids):
            df = self.df
        elif len(gene_caller_ids) == 1:
            df = self.df[self.df['gene_callers_id'] == gene_caller_ids.pop()]
        else:
            df = self.df[self.df['gene_callers_id'].isin(gene_caller_ids)]

        if group_name:
            df = df[df['data_group'] == group_name]

        if len(keys_of_interest):
            df = df[df['data_key'].isin(set(keys_of_interest))]

        # Assumes a data_key can possess only 1 data_type
        dtypes_convert = {
            'str': str,
            'int': int,
            'float': float,
            'stackedbar': str,
            'unknown': str,
        }

        dtypes = {}
        for data_key, subset in df.groupby('data_key'):
            dtypes[data_key] = dtypes_convert[subset['data_type'].iloc[0]]

        if df.empty:
            return pd.DataFrame({}, columns=('codon_order_in_gene',))

        pivot_gene_df = df.pipe(
            self._multi_index_pivot,
            index=['gene_callers_id', 'codon_order_in_gene'],
            columns='data_key',
            values='data_value'
        ).astype(dtypes).reset_index()

        return pivot_gene_df


    def _multi_index_pivot(self, df, index = None, columns = None, values = None):
        # https://github.com/pandas-dev/pandas/issues/23955
        output_df = df.copy(deep = True)
        if index is None:
            names = list(output_df.index.names)
            output_df = output_df.reset_index()
        else:
            names = index
        output_df = output_df.assign(tuples_index = [tuple(i) for i in output_df[names].values])
        if isinstance(columns, list):
            output_df = output_df.assign(tuples_columns = [tuple(i) for i in output_df[columns].values])  # hashable
            output_df = output_df.pivot(index = 'tuples_index', columns = 'tuples_columns', values = values)
            output_df.columns = pd.MultiIndex.from_tuples(output_df.columns, names = columns)  # reduced
        else:
            output_df = output_df.pivot(index = 'tuples_index', columns = columns, values = values)
        output_df.index = pd.MultiIndex.from_tuples(output_df.index, names = names)
        return output_df


    def get_gene_dataframe(self, gene_callers_id, keys_of_interest=set([]), group_name=None):
        """Fetches and pivots table data specific to a gene_callers_id

        Returns
        =======
        output : pd.DataFrame
            Outputs a dataframe that looks like this:

                codon_order_in_gene               MG                SM_               UDP
                149                              NaN    0.0349008254091               NaN
                150                              NaN    0.0116165975232               NaN
                167                              NaN    0.1009224794937   0.1127902398775
                168                              NaN    0.0110966620296  0.03150767340385
                169                              NaN     0.435271159751   0.4638648995855
                171                              NaN    0.0239539337947  0.03150767340385
                172                              NaN    0.0347266126837               NaN
                173                              NaN     0.244954708876   0.1630916348525
                174                              NaN    0.5049807590385    0.604658463037
                175                              NaN    0.0237121180997  0.06609853613995
                176                              NaN    0.0359138454456  0.10012280860365
                177                              NaN    0.0252309121402   0.0691817254722
                178                              NaN   0.02856525986075  0.03150767340385
                179                              NaN   0.01191453799525               NaN
                ...                              ...                ...               ...

            Only keys that had at least one non-NaN value for the gene are included
        """

        pivot_gene_df = self.get_multigene_dataframe(set([gene_callers_id]))

        if not pivot_gene_df.empty:
            pivot_gene_df.drop('gene_callers_id', axis=1, inplace=True)

        return pivot_gene_df


    def init_table_as_dataframe(self):
        """Call AdditionalAndOrderDataBaseClass.init_table_as_dataframe and split `item_name`

        Splits the item_name column (format <gene_callers_id>:<codon_order_in_gene>) into two
        separate columns, `gene_callers_id` and `codon_order_in_gene`, thereby recovering this
        information for parseability, etc.
        """

        super(AdditionalDataBaseClass, self).init_table_as_dataframe()

        if self.df.empty:
            self.df['gene_callers_id'] = None
            self.df['codon_order_in_gene'] = None
            return

        self.df[['gene_callers_id', 'codon_order_in_gene']] = self.df['item_name'].str.split(':', expand=True)
        self.df['gene_callers_id'] = self.df['gene_callers_id'].astype(int)
        self.df['codon_order_in_gene'] = self.df['codon_order_in_gene'].astype(int)


    def check_names(self, data_dict):
        """Compares data key values found in the data dict to the ones in the db

        FIXME This is an unforgiving function, i.e. --just-do-it does nothing
        """

        database = db.DB(self.db_path, get_required_version_for_db(self.db_path))
        gene_lengths = database.get_table_as_dataframe(
            t.genes_in_contigs_table_name,
            columns_of_interest=('gene_callers_id', 'start', 'stop'),
        )
        gene_lengths = dict(zip(gene_lengths['gene_callers_id'], gene_lengths['stop'] - gene_lengths['start']))

        is_format_good = lambda key: len(key.split(':')) == 2

        # These tables could be millions of entries. We do not bother compiling a prettified list of
        # their badly formatted table. If it's wrong we complain immediately.
        for i, data_key in enumerate(data_dict):
            if not is_format_good(data_key):
                raise ConfigError("The data key '%s' (entry #%d) is not properly formatted. It should "
                                  "have the format <gene_callers_id>:<codon_order_in_gene>, "
                                  "e.g. 1248:122 would be an entry for the 122nd amino acid "
                                  "(0-indexed) in gene 1248." % (data_key, i+1))

            gene_callers_id, codon_order_in_gene = data_key.split(':')
            gene_callers_id = int(gene_callers_id)
            codon_order_in_gene = int(codon_order_in_gene)

            if gene_callers_id not in gene_lengths:
                raise ConfigError("There is a problem with data key '%s' (entry #%d). It specifies gene callers "
                                  "ID %d, which does not exist." % \
                                  (data_key, i+1, gene_callers_id))

            if codon_order_in_gene < 0 or codon_order_in_gene >= gene_lengths[gene_callers_id]:
                raise ConfigError("There is a problem with data key '%s' (entry #%d). It specifies codon_order_in_gene "
                                  "%d, which falls outside of gene %d (its length is %d)." % \
                                  (data_key, i+1, codon_order_in_gene, gene_callers_id, gene_lengths[gene_callers_id]))

        database.disconnect()


class MiscDataTableFactory(TableForItemAdditionalData, TableForLayerAdditionalData, TableForLayerOrders):
    """Gives seamless access to additional data or order tables in pan/profile/contigs databases.

    Create an instance with args.target_data_table =
    [items|layers|layer_orders|nucleotides|amino_acids], and you will be golden.
    """

    def __init__(self, args, r=terminal.Run(), p=terminal.Progress()):
        self.run = r
        self.progress = p

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        target_data_table = A('target_data_table')

        if not target_data_table:
            raise ConfigError("When creating an instance from the MiscDataTableFactory class, the `args` object "
                              "must contain the `target_data_table` variable.")

        try:
            table_classes[target_data_table].__init__(self, args, r=self.run, p=self.progress)
        except KeyError:
            raise ConfigError("MiscDataTableFactory does not know about target data tables for '%s' :( "
                              "You can go to the online documentation, or you can try any of %s" % \
                              (target_data_table, ', '.join(["'%s'" % table for table in table_classes])))


table_classes = {
    'items': TableForItemAdditionalData,
    'layers': TableForLayerAdditionalData,
    'layer_orders': TableForLayerOrders,
    'nucleotides': TableForNucleotideAdditionalData,
    'amino_acids': TableForAminoAcidAdditionalData,
}
