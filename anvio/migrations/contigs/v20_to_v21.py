#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio.db as db
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError

current_version, next_version = [x[1:] for x in __name__.split('_to_')]

gene_function_reactions_table_name        = 'gene_function_reactions'
gene_function_reactions_table_structure   = ['modelseed_reaction_id', 'modelseed_reaction_name', 'ko_kegg_reaction_source', 'ko_ec_number_source', 'metabolite_modelseed_ids', 'stoichiometry', 'compartments', 'reversibility']
gene_function_reactions_table_types       = [         'text'        ,            'text'        ,           'text'         ,         'text'       ,           'text'          ,      'text'    ,      'str'    ,      'bool'    ]

gene_function_metabolites_table_name      = 'gene_function_metabolites'
gene_function_metabolites_table_structure = ['modelseed_compound_id', 'modelseed_compound_name', 'formula', 'charge']
gene_function_metabolites_table_types     = [         'text'        ,             'text'       ,   'text' , 'numeric']

run = terminal.Run()
progress = terminal.Progress()

def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    utils.is_contigs_db(db_path)

    contigs_db = db.DB(db_path, None, ignore_version = True)
    if str(contigs_db.get_version()) != current_version:
        raise ConfigError("Version of this contigs database is not %s (hence, this script cannot really do anything)." % current_version)

    progress.new("Migrating")
    progress.update("Creating two new tables for reactions and metabolites")

    # to be on the safe side, remove any reaction network tables and metadata that might exist
    try:
        contigs_db.drop_table(gene_function_reactions_table_name)
        contigs_db.drop_table(gene_function_metabolites_table_name)
    except:
        pass

    try:
        contigs_db.remove_meta_key_value_pair('reaction_network_was_run')
        contigs_db.remove_meta_key_value_pair('reaction_network_kegg_database_version')
        contigs_db.remove_meta_key_value_pair('reaction_network_modelseed_database_hash')
    except:
        pass

    contigs_db.set_meta_value('reaction_network_was_run', False)
    contigs_db.set_meta_value('reaction_network_kegg_database_version', None)
    contigs_db.set_meta_value('reaction_network_modelseed_database_hash', None)
    contigs_db.create_table(gene_function_reactions_table_name, gene_function_reactions_table_structure, gene_function_reactions_table_types)
    contigs_db.create_table(gene_function_metabolites_table_name, gene_function_metabolites_table_structure, gene_function_metabolites_table_types)

    progress.update("Updating version")
    contigs_db.remove_meta_key_value_pair('version')
    contigs_db.set_version(next_version)

    progress.update("Committing changes")
    contigs_db.disconnect()

    progress.end()

    message = (
        "Congratulations! Your contigs database is now version 21, which means it now contains two new empty tables. "
        "These tables are not as boring as they first appear, because you can now run `anvi-reaction-network` to store "
        "a network of the metabolic reactions that may be encoded by genes. A metabolic model representing the network "
        "can be exported from the database using `anvi-get-metabolic-model-file`."
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
