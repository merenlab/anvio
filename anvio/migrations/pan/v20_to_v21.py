#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio.dbinfo as dbinfo
import anvio.terminal as terminal

from anvio.errors import ConfigError

run = terminal.Run()
progress = terminal.Progress()

current_version, next_version = [x[1:] for x in __name__.split('_to_')]


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    pan_db_info = dbinfo.PanDBInfo(db_path)
    if str(pan_db_info.version) != current_version:
        raise ConfigError(
            f"The version of the provided pan database is {pan_db_info.version}, not the required "
            f"version, {current_version}, so this script cannot upgrade the database."
        )

    pan_db = pan_db_info.load_db()

    progress.new("Migrating")
    progress.update("Updating the self table with two variables if not already there")

    threshold_added = False
    if 'reaction_network_consensus_threshold' not in pan_db_info.get_self_table():
        pan_db.set_meta_value('reaction_network_consensus_threshold', None)
        threshold_added = True

    discard_ties_added = False
    if 'reaction_network_discard_ties' not in pan_db_info.get_self_table():
        pan_db.set_meta_value('reaction_network_discard_ties', None)
        discard_ties_added = True

    progress.update("Updating version")
    pan_db.remove_meta_key_value_pair('version')
    pan_db.set_version(next_version)

    progress.update("Committing changes")
    pan_db.disconnect()

    progress.end()

    if threshold_added and discard_ties_added:
        change_message = (
            "Two placeholder variables were added to the self table. These are filled in when a "
            "reaction network is generated via `anvi-reaction-network`."
        )
    elif threshold_added and not discard_ties_added:
        change_message = (
            "A placeholder variable, 'reaction_network_consensus_threshold', was added to the self "
            "table. Strangely, the variable, 'reaction_network_discard_ties', which goes "
            "hand-in-hand with the other variable, was already present, and since we trust that "
            "you know what's up, we left it alone. These variables are filled in when a reaction "
            "network is generated via `anvi-reaction-network`."
        )
    elif discard_ties_added and not threshold_added:
        change_message = (
            "A placeholder variable, 'reaction_network_discard_ties', was added to the self table. "
            "Strangely, the variable, 'reaction_network_consensus_threshold', which goes "
            "hand-in-hand with the other variable, was already present, and since we trust that "
            "you know what's up, we left it alone. These variables are filled in when a reaction "
            "network is generated via `anvi-reaction-network`."
        )
    else:
        change_message = (
            "The variables, 'reaction_network_consensus_threshold' and "
            "'reaction_network_discard_ties', were already found in the self table, likely because "
            "you have run `anvi-reaction-network`, so we left them alone and didn't change "
            "anything in the database besides the version number."
        )
    message = f"Done! Your pan database is now version {current_version}. This wasn't a biggie. {change_message}"
    run.info_single(message, nl_after=1, nl_before=1, mc='green')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade the pan database from version %s to version %s' % (current_version, next_version))
    parser.add_argument('pan_db', metavar = 'PAN_DB', help = "An anvi'o pan database of version %s" % current_version)
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.pan_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
