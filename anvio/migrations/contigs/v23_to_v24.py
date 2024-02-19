#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio.db as db
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError

current_version, next_version = [x[1:] for x in __name__.split('_to_')]

reaction_network_kegg_table_name      = 'reaction_network_kegg'
reaction_network_kegg_table_structure = ['kegg_id', 'name', 'modules', 'pathways', 'brite_categorization']
reaction_network_kegg_table_types     = [ 'text'  , 'text',  'text'  ,   'text'  ,         'text'        ]

run = terminal.Run()
progress = terminal.Progress()

def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    utils.is_contigs_db(db_path)

    contigs_db = db.DB(db_path, None, ignore_version=True)
    if str(contigs_db.get_version()) != current_version:
        raise ConfigError(
            f"The version of the provided contigs database is {contigs_db.get_version}, not the "
            f"required version, {current_version}, so this script cannot upgrade the database."
        )

    progress.new("Migrating")
    progress.update("Creating a new table for KEGG KO information in a reaction network")

    # To be on the safe side, remove any KEGG table that may already exist.
    try:
        contigs_db.drop_table(reaction_network_kegg_table_name)
    except:
        pass

    contigs_db.create_table(
        reaction_network_kegg_table_name,
        reaction_network_kegg_table_structure,
        reaction_network_kegg_table_types
    )

    progress.update("Renaming other reaction network tables")

    # To be on the safe side, remove any tables with the new names that may already exist.
    try:
        contigs_db.drop_table('reaction_network_reactions')
        contigs_db.drop_table('reaction_network_metabolites')
    except:
        pass

    contigs_db._exec('ALTER TABLE gene_function_reactions RENAME TO reaction_network_reactions')
    contigs_db._exec('ALTER TABLE gene_function_metabolites RENAME TO reaction_network_metabolites')

    progress.update("Updating version")
    contigs_db.remove_meta_key_value_pair('version')
    contigs_db.set_version(next_version)

    progress.update("Committing changes")
    contigs_db.disconnect()

    progress.end()

    message = (
        "Congratulations! Your contigs database is now version 23. An empty table has been added "
        "to improve the functionality of reaction networks, particularly their portability and "
        "reproducibility. The two existing tables for storing reaction networks have also been "
        "renamed for the sake of clarity."
    )
    run.info_single(message, nl_after=1, nl_before=1, mc='green')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade CONTIGS.db from version %s to version %s' % (current_version, next_version))
    parser.add_argument('contigs_db', metavar = 'CONTIGS_DB', help = 'Contigs database at version %s' % current_version)
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.contigs_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
