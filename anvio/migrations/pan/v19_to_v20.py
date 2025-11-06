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

current_version, next_version = [x[1:] for x in __name__.split('_to_')]

pan_reaction_network_metabolites_table_structure = ['modelseed_compound_id', 'modelseed_compound_name', 'kegg_aliases', 'formula', 'charge' , 'smiles']
pan_reaction_network_metabolites_table_types     = [         'text'        ,           'text'         ,     'text'    ,   'text' , 'numeric',  'text' ]

run = terminal.Run()
progress = terminal.Progress()

def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    pan_db = db.DB(db_path, None, ignore_version=True)
    if str(pan_db.get_version()) != current_version:
        raise ConfigError(
            f"The version of the provided pan database is {pan_db.get_version}, not the required "
            f"version, {current_version}, so this script cannot upgrade the database."
        )

    added_smiles_strings = add_smiles_column(pan_db)

    progress.new("Updating version")
    progress.update("...")
    pan_db.remove_meta_key_value_pair('version')
    pan_db.set_version(next_version)
    pan_db.disconnect()
    progress.end()

    if added_smiles_strings:
        smiles_message = (
            "SMILES string compound structural data was added to the existing reaction network."
        )
    else:
        smiles_message = (
            "A new column for storage of SMILES string compound structural data was added to the "
            "metabolites table."
        )
    message = (
        f"Congratulations! Your pan database is now version {next_version}. {smiles_message}"
    )
    run.info_single(message, nl_after=1, nl_before=1, mc='green')

def add_smiles_column(pan_db: db.DB) -> bool:
    """
    Add a SMILES string column to the pan reaction network metabolites table.

    Parameters
    ==========
    pan_db : db.DB
        Database being migrated.

    Returns
    =======
    bool
        True if SMILES strings could be retrieved from a ModelSEED Biochemistry reference database
        and added. False if only an empty column was added.
    """
    modelseed_db_available = check_modelseed_database(pan_db)

    metabolites_table = pan_db.get_table_as_dataframe('pan_reaction_network_metabolites', error_if_no_data=False)
    
    if 'smiles' in metabolites_table.columns:
        pan_db._exec('ALTER TABLE pan_reaction_network_metabolites DROP smiles')
        run.warning("An existing column by the name of 'smiles' was dropped.")

    progress.new("Adding metabolite SMILES string column")
    progress.update("...")
    pan_db._exec('ALTER TABLE pan_reaction_network_metabolites ADD smiles text')

    if not modelseed_db_available:
        progress.end()
        return False

    modelseed_table = ModelSEEDDatabase().compounds_table
    smiles_strings = []
    for row in metabolites_table.itertuples():
        smiles = modelseed_table.loc[row.modelseed_compound_id]['smiles']
        if pd.isna(smiles):
            smiles_strings.append('')
        else:
            smiles_strings.append(smiles)
    metabolites_table['smiles'] = smiles_strings

    pan_db.drop_table('pan_reaction_network_metabolites')
    pan_db.create_table(
        'pan_reaction_network_metabolites',
        pan_reaction_network_metabolites_table_structure,
        pan_reaction_network_metabolites_table_types
    )
    pan_db.insert_rows_from_dataframe('pan_reaction_network_metabolites', metabolites_table)
    progress.end()

    return True

def check_modelseed_database(pan_db: db.DB) -> bool:
    """
    Check if the pan database contains a reaction network that was constructed with a ModelSEED
    Biochemistry database installed at the default anvi'o location.

    Parameters
    ==========
    pan_db : db.DB
        Database being migrated.

    Returns
    =======
    bool
        True if the pan database contains a reaction network constructed with a ModelSEED database
        installed at the default anvi'o location.
    """
    network_sha: str = pan_db.get_meta_value('reaction_network_modelseed_database_sha')
    if not network_sha:
        return False

    sha_txt_path = os.path.join(ModelSEEDDatabase.default_dir, 'sha.txt')
    compounds_db_path = os.path.join(ModelSEEDDatabase.default_dir, 'compounds.tsv')
    if not os.path.isfile(sha_txt_path) or not os.path.isfile(compounds_db_path):
        run.warning(
            "A ModelSEED Biochemistry database was not found to be set up in the default anvi'o "
            f"directory, '{ModelSEEDDatabase.default_dir}', preventing this script from adding "
            "SMILES string structural data to the table of reaction network metabolites in the "
            "pan database. The reaction network can be reconstructed and overwritten to store "
            "SMILES strings."
        )
        return False

    with open(sha_txt_path) as sha_txt:
        ref_sha = sha_txt.read()
    if network_sha != ref_sha:
        run.warning(
            f"The ID ('{network_sha}') of the ModelSEED Biochemistry database used to build the "
            f"reaction network stored in the pan database does not match the ID ('{ref_sha}') of "
            "the reference database set up in the default anvi'o directory, indicating that the "
            "network was built with a different version of ModelSEED. This script therefore cannot "
            "add SMILES string structural data to the table of reaction network metabolites in the "
            "pan database. The reaction network can be reconstructed and overwritten to store "
            "SMILES strings."
        )
        return False

    return True

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade the pan database from version %s to version %s' % (current_version, next_version))
    parser.add_argument('pan_db', metavar = 'PAN_DB', help = "An anvi'o pan database of version %s" % current_version)
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.pan_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
