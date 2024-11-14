#!/usr/bin/env python
# -*- coding: utf-8

import os
import sys
import shutil
import tempfile
import argparse

import anvio.db as db
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError

current_version, next_version = [x[1:] for x in __name__.split("_to_")]

run = terminal.Run()
progress = terminal.Progress()


item_additional_data_table_name = "item_additional_data"
item_additional_data_table_structure = [
    "entry_id",
    "item_name",
    "data_key",
    "data_value",
    "data_type",
]
item_additional_data_table_types = ["numeric", "text", "text", "text", "text"]

layer_orders_table_name = "layer_orders"
layer_orders_table_structure = ["data_key", "data_type", "data_value"]
layer_orders_table_types = ["text", "text", "text"]

layer_additional_data_table_name = "layer_additional_data"
layer_additional_data_table_structure = [
    "entry_id",
    "item_name",
    "data_key",
    "data_value",
    "data_type",
]
layer_additional_data_table_types = ["numeric", "text", "text", "text", "text"]


def get_temp_file_path():
    f = tempfile.NamedTemporaryFile(delete=False)
    temp_file_name = f.name
    f.close()
    return temp_file_name


class Table(object):
    def __init__(
        self, db_path, version, run=run, progress=progress, quiet=False, simple=False
    ):
        self.quiet = quiet
        self.db_type = None
        self.db_path = db_path
        self.version = version
        self.next_available_id = {}

        self.splits_info = None
        self.contigs_info = None
        self.split_length = None
        self.genes_are_called = None

        self.run = run
        self.progress = progress

        database = db.DB(self.db_path, version, ignore_version=True)
        self.db_type = database.get_meta_value("db_type")

    def next_id(self, table):
        if table not in self.next_available_id:
            raise ConfigError(
                "If you need unique ids, you must call 'set_next_available_id' first"
            )

        self.next_available_id[table] += 1
        return self.next_available_id[table] - 1

    def set_next_available_id(self, table):
        database = db.DB(self.db_path, self.version, ignore_version=True)
        table_content = database.get_table_as_dict(table)
        if table_content:
            self.next_available_id[table] = max(table_content.keys()) + 1
        else:
            self.next_available_id[table] = 0

        database.disconnect()


class SamplesInformationDatabase:
    def __init__(self, db_path):
        self.db_path = db_path

        self.meta = {}
        self.init()

    def init(self):
        self.db = db.DB(self.db_path, "2")
        meta_table = self.db.get_table_as_dict("self")
        self.meta = dict([(k, meta_table[k]["value"]) for k in meta_table])
        self.samples = set([s.strip() for s in self.meta["samples"].split(",")])
        self.sample_names_for_order = (
            set([s.strip() for s in self.meta["sample_names_for_order"].split(",")])
            if self.meta["sample_names_for_order"]
            else self.samples
        )
        self.samples_information_default_layer_order = self.meta[
            "samples_information_default_layer_order"
        ].split(",")

    def recover_samples_information_dict(
        self, samples_information_dict_from_db, aliases_to_attributes_dict
    ):
        samples_information_dict_with_attributes = {}

        for sample_name in samples_information_dict_from_db:
            samples_information_dict_with_attributes[sample_name] = {}

        for alias in aliases_to_attributes_dict:
            attribute = aliases_to_attributes_dict[alias]["attribute"]
            for sample_name in samples_information_dict_with_attributes:
                samples_information_dict_with_attributes[sample_name][attribute] = (
                    samples_information_dict_from_db[sample_name][alias]
                )

        return samples_information_dict_with_attributes

    def get_samples_information_and_order_dicts(self):
        samples_information_dict = self.recover_samples_information_dict(
            self.db.get_table_as_dict("samples_information", error_if_no_data=False),
            self.db.get_table_as_dict(
                "samples_attribute_aliases", error_if_no_data=False
            ),
        )
        samples_order_dict = self.db.get_table_as_dict("samples_order")

        for key in [k for k in list(samples_order_dict.keys()) if k.startswith(">>")]:
            samples_order_dict.pop(key)

        new_order = {}
        for i in samples_order_dict:
            data_type = "newick" if samples_order_dict[i]["newick"] else "basic"
            new_order[i] = {
                "data_type": data_type,
                "data_value": samples_order_dict[i][data_type],
            }

        return samples_information_dict, new_order

    def export_samples_db_files(self):
        """Export whatever information is stored in a ginve anvi'o samples database"""

        order_output_path = get_temp_file_path()
        information_output_path = get_temp_file_path()

        samples_information_dict, samples_order_dict = (
            self.get_samples_information_and_order_dicts()
        )

        utils.store_dict_as_TAB_delimited_file(
            samples_order_dict,
            order_output_path,
            headers=["attributes", "data_type", "data_value"],
        )
        utils.store_dict_as_TAB_delimited_file(
            samples_information_dict,
            information_output_path,
            headers=["samples"]
            + sorted(list(list(samples_information_dict.values())[0].keys())),
        )

        return information_output_path, order_output_path


class AdditionalAndOrderDataBaseClass(Table, object):
    """This is a base class for common operations between order and additional data classes."""

    def __init__(self, args):
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.db_path = A("pan_or_profile_db") or A("profile_db") or A("pan_db")
        self.just_do_it = A("just_do_it")

        if not self.db_path:
            raise ConfigError(
                "The AdditionalAndOrderDataBaseClass is inherited with an args object that did not "
                "contain any database path :/ Even though any of the following would "
                "have worked: `pan_or_profile_db`, `profile_db`, `pan_db` :("
            )

        if not self.table_name:
            raise ConfigError(
                "The AdditionalAndOrderDataBaseClass does not know anything about the table it should "
                "be working with."
            )

        database = db.DB(self.db_path, None, ignore_version=True)
        self.additional_data_keys = database.get_single_column_from_table(
            self.table_name, "data_key"
        )
        database.disconnect()

        Table.__init__(self, self.db_path, None, self.run, self.progress)

    def populate_from_file(self, additional_data_file_path, skip_check_names=None):
        data_keys = utils.get_columns_of_TAB_delim_file(additional_data_file_path)
        data_dict = utils.get_TAB_delimited_file_as_dictionary(
            additional_data_file_path
        )

        if not len(data_keys):
            raise ConfigError(
                "There is something wrong with the additional data file for %s at %s. "
                "It does not seem to have any additional keys for data :/"
                % (self.target, additional_data_file_path)
            )

        if self.target == "layer_orders":
            OrderDataBaseClass.add(self, data_dict, skip_check_names)
        else:
            AdditionalDataBaseClass.add(self, data_dict, data_keys, skip_check_names)


class OrderDataBaseClass(AdditionalAndOrderDataBaseClass, object):
    """Implements a base class to deal with tables that keep order data."""

    def __init__(self, args):
        AdditionalAndOrderDataBaseClass.__init__(self, args)

    def add(self, data_dict, skip_check_names=False):
        data_keys_list = list(data_dict.keys())
        data_key_types = {}
        for key in data_keys_list:
            predicted_key_type = data_dict[key]["data_type"]
            data_key_types[key] = predicted_key_type

        db_entries = []
        for item_name in data_dict:
            db_entries.append(
                tuple(
                    [
                        item_name,
                        data_dict[item_name]["data_type"],
                        data_dict[item_name]["data_value"],
                    ]
                )
            )

        database = db.DB(self.db_path, None, ignore_version=True)
        database._exec_many(
            """INSERT INTO %s VALUES (?,?,?)""" % self.table_name, db_entries
        )
        database.disconnect()


class AdditionalDataBaseClass(AdditionalAndOrderDataBaseClass, object):
    """Implements additional data ops base class.

    See TableForItemAdditionalData or TableForLayerAdditionalData for usage example.

    See AdditionalAndOrderDataBaseClass for inherited functionality.
    """

    def __init__(self, args):
        AdditionalAndOrderDataBaseClass.__init__(self, args)

    def add(self, data_dict, data_keys_list, skip_check_names=False):
        key_types = {}
        for key in data_keys_list:
            if "!" in key:
                predicted_key_type = "stackedbar"
            else:
                type_class = utils.get_predicted_type_of_items_in_a_dict(data_dict, key)
                predicted_key_type = type_class.__name__ if type_class else None

            key_types[key] = predicted_key_type

        db_entries = []
        self.set_next_available_id(self.table_name)
        for item_name in data_dict:
            for key in data_dict[item_name]:
                db_entries.append(
                    tuple(
                        [
                            self.next_id(self.table_name),
                            item_name,
                            key,
                            data_dict[item_name][key],
                            key_types[key],
                        ]
                    )
                )

        database = db.DB(self.db_path, None, ignore_version=True)
        database._exec_many(
            """INSERT INTO %s VALUES (?,?,?,?,?)""" % self.table_name, db_entries
        )
        database.disconnect()


class TableForLayerAdditionalData(AdditionalDataBaseClass):
    def __init__(self, args, r=run, p=progress):
        self.run = r
        self.progress = p

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.table_name = A("table_name") or layer_additional_data_table_name

        self.target = "layers"

        AdditionalDataBaseClass.__init__(self, args)


class TableForLayerOrders(OrderDataBaseClass):
    def __init__(self, args, r=run, p=progress):
        self.run = r
        self.progress = p

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.table_name = A("table_name") or layer_orders_table_name

        self.allowde_types = ["newick", "basic"]
        self.target = "layer_orders"

        OrderDataBaseClass.__init__(self, args)


def check_samples_db_status():
    if "ANVIO_SAMPLES_DB" not in os.environ or os.environ["ANVIO_SAMPLES_DB"] == "SKIP":
        samples_db_path = None
    else:
        samples_db_path = os.environ["ANVIO_SAMPLES_DB"]

        if not os.path.exists(samples_db_path):
            raise ConfigError(
                "Your migration did not finish, and your profile database is still at %s. Although anvi'o found "
                "the environmental variable ANVIO_SAMPLES_DB, the path it pointed, '%s', was nowhere to be found. "
                "If you don't want to incorporate the information in the samples database associated with this "
                "profile datbase you can simply call the migration script this way: 'ANVIO_SAMPLES_DB=SKIP anvi-migrate YOUR_PROFILE_DB_PATH'. "
                "Otherwise, try again with a proper path."
                % (current_version, samples_db_path)
            )

        try:
            database = db.DB(samples_db_path, None, ignore_version=True)
        except:
            raise ConfigError(
                "The file at %s does not look like an anvi'o database" % samples_db_path
            )

        tables = database.get_table_names()
        if "self" not in tables:
            database.disconnect()
            raise ConfigError(
                "'%s' does not seem to be a anvi'o database..." % samples_db_path
            )

        if database.get_meta_value("db_type") != "samples_information":
            raise ConfigError(
                "'%s' does not seem to be a anvi'o samples database..."
                % samples_db_path
            )

        database.disconnect()

    return samples_db_path


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    utils.is_profile_db(db_path)

    profile_db = db.DB(db_path, None, ignore_version=True)
    if str(profile_db.get_version()) != current_version:
        raise ConfigError(
            "Version of this profile database is not %s (hence, this script cannot really do anything)."
            % current_version
        )

    # check samples db
    samples_db_path = check_samples_db_status()

    # start by adding new tables...
    profile_db.create_table(
        layer_orders_table_name, layer_orders_table_structure, layer_orders_table_types
    )
    profile_db.create_table(
        layer_additional_data_table_name,
        layer_additional_data_table_structure,
        layer_additional_data_table_types,
    )

    # update the item_additional_data table
    profile_db.cursor.execute(
        "ALTER TABLE item_additional_data RENAME TO item_additional_data_TEMP;"
    )
    profile_db.cursor.execute(
        "CREATE TABLE item_additional_data (entry_id numeric, item_name text, data_key text, data_value text, data_type text);"
    )
    profile_db.cursor.execute(
        "INSERT INTO item_additional_data(entry_id, item_name, data_key, data_value, data_type) SELECT entry_id, item_name, key, value, type FROM item_additional_data_TEMP;"
    )
    profile_db.cursor.execute("DROP TABLE item_additional_data_TEMP;")

    profile_db.remove_meta_key_value_pair("version")
    profile_db.set_version(next_version)
    profile_db.disconnect()

    if samples_db_path:
        try:
            samples_db = SamplesInformationDatabase(samples_db_path)
            layers_info_path, layers_order_path = samples_db.export_samples_db_files()

            args = argparse.Namespace(profile_db=db_path, target_data_table="layers")
            TableForLayerAdditionalData(args).populate_from_file(layers_info_path)

            args = argparse.Namespace(profile_db=db_path, target_data_table="layers")
            TableForLayerOrders(args).populate_from_file(layers_order_path)

            os.remove(layers_info_path)
            os.remove(layers_order_path)

            fully_upgraded = True
        except Exception as e:
            run.warning(
                "Something went wrong adding the data found in samples database into the profile database. This is what "
                'we know: "%s".' % e
            )
            fully_upgraded = False
    else:
        fully_upgraded = False

    if fully_upgraded:
        shutil.move(samples_db_path, samples_db_path + ".OBSOLETE")
        run.info_single(
            "Your profile db is now version %s. You no longer need your old samples database (which is now "
            "renamed to something ugly so you can see it." % next_version,
            nl_after=1,
            nl_before=1,
            mc="green",
        )
    elif samples_db_path:
        run.info_single(
            "Your profile db is now version %s. BUT THERE WAS THIS: the actual purpose of this script was to "
            "incorporate the data in your samples database into your profile database. But for some reason it "
            "has failed. Probably everything is still alright, but you may have to do that step manually. The "
            "Error messsage should be somewhere above." % next_version,
            nl_after=1,
            nl_before=1,
            mc="green",
        )
    else:
        run.info_single(
            "Your profile db is now version %s. BUT WITHOUT the samples database incorporation as you wished."
            % next_version,
            nl_after=1,
            nl_before=1,
            mc="green",
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="A simple script to upgrade profile database and AUXILIARY-DATA.h5 from version %s to version %s"
        % (current_version, next_version)
    )
    parser.add_argument(
        "profile_db",
        metavar="PROFILE_DB",
        help="An anvi'o profile database of version %s" % current_version,
    )
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.profile_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
