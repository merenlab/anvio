#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio.dbinfo as dbinfo
import anvio.terminal as terminal

from anvio.errors import ConfigError

current_version, next_version = [x[1:] for x in __name__.split('_to_')]

protein_abundances_table_name = 'protein_abundances'
protein_abundances_table_structure = [
    'protein_id',
    'reference_source',
    'reference_id',
    'sample_name',
    'abundance_value'
]
protein_abundances_table_types = [
    'numeric',
    'text',
    'text',
    'text',
    'numeric'
]

metabolite_abundances_table_name = 'metabolite_abundances'
metabolite_abundances_table_structure = [
    'reference_source',
    'reference_id',
    'sample_name',
    'abundance_value'
]
metabolite_abundances_table_types = [
    'text',
    'text',
    'text',
    'numeric'
]

run = terminal.Run()
progress = terminal.Progress()

def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    profile_db_info = dbinfo.ProfileDBInfo(db_path)
    if str(profile_db_info.version) != current_version:
        raise ConfigError(
            f"""The version of the provided profile database is {profile_db_info.version}, not the \
            required version, {current_version}, so this script cannot upgrade the database."""
        )

    progress.new("Migrating")
    progress.update("Creating two new tables for protein and metabolite abundances")
    profile_db = profile_db_info.load_db()

    # To be on the safe side, remove any protein and metabolite abundance tables that might exist.
    try:
        profile_db.drop_table(protein_abundances_table_name)
        profile_db.drop_table(metabolite_abundances_table_name)
    except:
        pass

    profile_db.create_table(
        protein_abundances_table_name,
        protein_abundances_table_structure,
        protein_abundances_table_types
    )
    profile_db.create_table(
        metabolite_abundances_table_name,
        metabolite_abundances_table_structure,
        metabolite_abundances_table_types
    )

    progress.update("Updating version")
    profile_db.remove_meta_key_value_pair('version')
    profile_db.set_version(next_version)

    progress.update("Committing changes")
    profile_db.disconnect()

    progress.end()

    run.info_single(
        """Congratulations! Your profile database is now version 40, which means it now contains two \
        new empty tables. Now you can run `anvi-import-protein-profile` and \
        `anvi-import-metabolite-profile` to store protein and metabolite abundances across \
        samples.""",
        nl_after=1, nl_before=1, mc='green'
    )

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=
        f"""A simple script to upgrade an anvi'o profile database from version \
        {current_version} to version {next_version}"""
    )
    parser.add_argument(
        "profile_db",
        metavar="PROFILE_DB",
        help=f"An anvi'o profile database of version {current_version}"
    )
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.profile_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
