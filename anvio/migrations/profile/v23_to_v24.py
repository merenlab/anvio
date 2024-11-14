#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio.db as db
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError

run = terminal.Run()
progress = terminal.Progress()

current_version, next_version = [x[1:] for x in __name__.split("_to_")]


# Following classess are from anvio.tables.miscdata, and they are trimmed
# to the minimum necessary sizes. they are not suitable to use for anything
# but this upgrade script
class AdditionalAndOrderDataBaseClass(object):
    """This is a base class for common operations between order and additional data classes."""

    def __init__(self, args):
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.db_path = A("profile_db")
        self.just_do_it = True

        database = db.DB(self.db_path, None, ignore_version=True)
        self.additional_data_keys = database.get_single_column_from_table(
            self.table_name, "data_key"
        )

        self.entry_id = database.get_max_value_in_column(
            "layer_additional_data", "entry_id"
        )
        if not self.entry_id:
            self.entry_id = 0
        else:
            self.entry_id = int(self.entry_id) + 1

        database.disconnect()


class AdditionalDataBaseClass(AdditionalAndOrderDataBaseClass, object):
    """Implements additional data ops base class.

    See TableForItemAdditionalData or TableForLayerAdditionalData for usage example.

    See AdditionalAndOrderDataBaseClass for inherited functionality.
    """

    def __init__(self, args):
        AdditionalAndOrderDataBaseClass.__init__(self, args)

    def add(self, data_dict, data_keys_list, skip_check_names=False):
        self.run.warning(None, "New %s additional data..." % self.target, lc="yellow")
        keys_already_in_db = [
            c for c in data_keys_list if c in self.additional_data_keys
        ]
        if len(keys_already_in_db):
            self.run.warning(
                "The following keys in your data dict will replace the ones that are already "
                "in your database: %s." % (", ".join(keys_already_in_db))
            )

            self.remove(keys_already_in_db)

        db_entries = []
        for item_name in data_dict:
            for key in data_keys_list:
                db_entries.append(
                    tuple(
                        [
                            self.entry_id,
                            item_name,
                            key,
                            data_dict[item_name][key],
                            "int",
                        ]
                    )
                )
            self.entry_id += 1

        database = db.DB(self.db_path, None, ignore_version=True)
        database._exec_many(
            """INSERT INTO %s VALUES (?,?,?,?,?)""" % self.table_name, db_entries
        )
        database.disconnect()

        self.run.info(
            "New data added to the db for your %s" % self.target,
            "%s." % (", ".join(data_keys_list)),
            nl_after=1,
        )

    def remove(self, data_keys_list):
        database = db.DB(self.db_path, None, ignore_version=True)
        for key in data_keys_list:
            database._exec(
                '''DELETE from %s WHERE data_key="%s"''' % (self.table_name, key)
            )
        database.disconnect()


class TableForLayerAdditionalData(AdditionalDataBaseClass):
    def __init__(self, args, r=run, p=progress):
        self.run = r
        self.progress = p

        self.table_name = "layer_additional_data"
        self.target = "layers"

        AdditionalDataBaseClass.__init__(self, args)


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    # make sure someone is not being funny
    utils.is_profile_db(db_path)

    # make sure the version is accurate
    profile_db = db.DB(db_path, None, ignore_version=True)
    if str(profile_db.get_version()) != current_version:
        raise ConfigError(
            "Version of this profile database is not %s (hence, this script cannot really do anything)."
            % current_version
        )

    is_merged = profile_db.get_meta_value("merged")
    tables_in_db = profile_db.get_table_names()
    profile_db.disconnect()
    is_full_profile = (
        "mean_coverage_Q2Q3_splits" in tables_in_db
        or "atomic_data_splits" in tables_in_db
    )

    run.info("Profile db type", "Merged" if is_merged else "Single")
    run.info("Full profile", is_full_profile)

    if is_full_profile:
        profile_db = db.DB(db_path, None, ignore_version=True)
        total_reads_mapped = profile_db.get_meta_value("total_reads_mapped")
        samples = profile_db.get_meta_value("samples")
        profile_db.disconnect()

        layer_additional_data_table = TableForLayerAdditionalData(
            argparse.Namespace(profile_db=db_path)
        )

        # we will do this only for full merged or single profiles
        if is_merged:
            full_upgrade = True
            total_reads_mapped = [int(m) for m in total_reads_mapped.split(",")]
            samples = [s.strip() for s in samples.split(",")]
            d = dict(zip(samples, total_reads_mapped))
            data = {}
            for sample in samples:
                data[sample] = {"total_reads_mapped": d[sample]}

            layer_additional_data_table.add(data, ["total_reads_mapped"])
        else:
            total_reads_mapped = int(total_reads_mapped)
            layer_additional_data_table.add(
                {samples: {"total_reads_mapped": total_reads_mapped}},
                ["total_reads_mapped"],
            )

        full_upgrade = True
    else:
        full_upgrade = False

    progress.new("Finalizing profile database upgrade")
    progress.update("...")

    profile_db = db.DB(db_path, None, ignore_version=True)

    if full_upgrade:
        # remove stuff no longer necessary
        profile_db.remove_meta_key_value_pair("total_reads_mapped")

    # set the version
    profile_db.remove_meta_key_value_pair("version")
    profile_db.set_version(next_version)

    # bye
    profile_db.disconnect()
    progress.end()

    if full_upgrade:
        run.info_single(
            "Your profile db is now version %s. You can learn more about what happened here "
            "by taking a look at this issue: https://github.com/merenlab/anvio/issues/800"
            % next_version,
            nl_after=1,
            nl_before=1,
            mc="green",
        )
    else:
        run.info_single(
            "Your profile db is now version %s. But essentially nothing really happened to your "
            "database since it was a blank profile (which is OK, move along)."
            % next_version,
            nl_after=1,
            nl_before=1,
            mc="green",
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="A simple script to upgrade profile database from version %s to version %s"
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
