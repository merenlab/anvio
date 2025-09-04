#!/usr/bin/env python
# -*- coding: utf-8

import os
import sys
import argparse
import pandas as pd

import anvio.db as db
import anvio.terminal as terminal

from anvio.errors import ConfigError
from anvio.reactionnetwork import ModelSEEDDatabase
from anvio.dbinfo import is_contigs_db

current_version, next_version = [x[1:] for x in __name__.split('_to_')]

reaction_network_kegg_table_name      = 'reaction_network_kegg'
reaction_network_kegg_table_structure = ['kegg_id', 'name', 'modules', 'pathways', 'brite_categorization']
reaction_network_kegg_table_types     = [ 'text'  , 'text',  'text'  ,   'text'  ,         'text'        ]

reaction_network_metabolites_table_structure = ['modelseed_compound_id', 'modelseed_compound_name', 'kegg_aliases', 'formula', 'charge' , 'smiles']
reaction_network_metabolites_table_types     = [         'text'        ,           'text'         ,     'text'    ,   'text' , 'numeric',  'text' ]

run = terminal.Run()
progress = terminal.Progress()

def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    is_contigs_db(db_path)

    contigs_db = db.DB(db_path, None, ignore_version=True)
    if str(contigs_db.get_version()) != current_version:
        raise ConfigError(
            f"The version of the provided contigs database is {contigs_db.get_version}, not the "
            f"required version, {current_version}, so this script cannot upgrade the database."
        )

    progress.new("Creating a new table for KEGG KO information in a reaction network")
    progress.update("...")
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
    progress.end()

    added_smiles_strings = add_smiles_column(contigs_db)

    progress.new("Renaming reaction network tables")
    progress.update("...")
    # To be on the safe side, remove any tables with the new names that may already exist.
    try:
        contigs_db.drop_table('reaction_network_reactions')
        contigs_db.drop_table('reaction_network_metabolites')
    except:
        pass

    contigs_db._exec('ALTER TABLE gene_function_reactions RENAME TO reaction_network_reactions')
    contigs_db._exec('ALTER TABLE gene_function_metabolites RENAME TO reaction_network_metabolites')
    progress.end()

    progress.new("Updating version")
    progress.update("...")
    contigs_db.remove_meta_key_value_pair('version')
    contigs_db.set_version(next_version)

    contigs_db.disconnect()
    progress.end()

    if added_smiles_strings:
        smiles_message = "SMILES string structural data was added to the existing reaction network."
    else:
        smiles_message = (
            "A new column for storage of SMILES string structural data was added to the "
            "metabolites table."
        )
    message = (
        f"Congratulations! Your contigs database is now version {next_version}. An empty table has "
        "been added to improve the functionality of reaction networks, particularly their "
        "portability and reproducibility. The two existing tables for storing reaction networks "
        f"have been renamed for the sake of clarity. {smiles_message} "
    )
    run.info_single(message, nl_after=1, nl_before=1, mc='green')

def add_smiles_column(contigs_db: db.DB) -> bool:
    """
    Add a SMILES string column to the reaction network metabolites table.

    Parameters
    ==========
    contigs_db : db.DB
        Database being migrated.

    Returns
    =======
    bool
        True if SMILES strings could be retrieved from a ModelSEED Biochemistry reference database
        and added. False if only an empty column was added.
    """
    modelseed_db_available = check_modelseed_database(contigs_db)

    progress.new("Adding metabolite SMILES string column")
    progress.update("...")

    contigs_db._exec('ALTER TABLE gene_function_metabolites ADD smiles text')

    if not modelseed_db_available:
        progress.end()
        return False

    metabolites_table = contigs_db.get_table_as_dataframe('gene_function_metabolites')
    modelseed_table = ModelSEEDDatabase().compounds_table
    smiles_strings = []
    for row in metabolites_table.itertuples():
        smiles = modelseed_table.loc[row.modelseed_compound_id]['smiles']
        if pd.isna(smiles):
            smiles_strings.append('')
        else:
            smiles_strings.append(smiles)
    metabolites_table['smiles'] = smiles_strings

    contigs_db.drop_table('gene_function_metabolites')
    contigs_db.create_table(
        'gene_function_metabolites',
        reaction_network_metabolites_table_structure,
        reaction_network_metabolites_table_types
    )
    contigs_db.insert_rows_from_dataframe('gene_function_metabolites', metabolites_table)
    progress.end()

    return True

def check_modelseed_database(contigs_db: db.DB) -> bool:
    """
    Check if the contigs database contains a reaction network that was constructed with a ModelSEED
    Biochemistry database installed at the default anvi'o location.

    Parameters
    ==========
    contigs_db : db.DB
        Database being migrated.

    Returns
    =======
    bool
        True if the contigs database contains a reaction network constructed with a ModelSEED
        database installed at the default anvi'o location.
    """
    network_sha: str = contigs_db.get_meta_value('reaction_network_modelseed_database_sha')
    if not network_sha:
        return False

    sha_txt_path = os.path.join(ModelSEEDDatabase.default_dir, 'sha.txt')
    compounds_db_path = os.path.join(ModelSEEDDatabase.default_dir, 'compounds.tsv')
    if not os.path.isfile(sha_txt_path) or not os.path.isfile(compounds_db_path):
        run.warning(
            "A ModelSEED Biochemistry database was not found to be set up in the default anvi'o "
            f"directory, '{ModelSEEDDatabase.default_dir}', preventing this script from adding "
            "SMILES string structural data to the table of reaction network metabolites in the "
            "contigs database. The reaction network can be reconstructed and overwritten to store "
            "SMILES strings."
        )
        return False

    with open(sha_txt_path) as sha_txt:
        ref_sha = sha_txt.read()
    if network_sha != ref_sha:
        run.warning(
            f"The ID ('{network_sha}') of the ModelSEED Biochemistry database used to build the "
            f"reaction network stored in the contigs database does not match the ID ('{ref_sha}') "
            "of the reference database set up in the default anvi'o directory, indicating that the "
            "network was built with a different version of ModelSEED. This script therefore cannot "
            "add SMILES string structural data to the table of reaction network metabolites in the "
            "contigs database. The reaction network can be reconstructed and overwritten to store "
            "SMILES strings."
        )
        return False

    return True

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade CONTIGS.db from version %s to version %s' % (current_version, next_version))
    parser.add_argument('contigs_db', metavar = 'CONTIGS_DB', help = 'Contigs database at version %s' % current_version)
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.contigs_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
