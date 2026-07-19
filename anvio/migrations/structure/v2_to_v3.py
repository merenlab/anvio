#!/usr/bin/env python

import sys
import argparse

import anvio.db as db
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError

run = terminal.Run()
progress = terminal.Progress()

current_version, next_version = [x[1:] for x in __name__.split('_to_')]

# new columns introduced in v3 to accommodate a second structure-prediction engine (ColabFold).
# ColabFold reports confidence metrics (per-residue pLDDT and model-level pTM) rather than
# MODELLER's DOPE/GA341/molpdf scores, so these columns are simply left NULL for existing
# (MODELLER or external) structures.
new_columns = {
    'models':       [('mean_plddt', 'real'), ('ptm', 'real')],
    'residue_info': [('plddt', 'real')],
}


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    # make sure someone is not being funny
    utils.is_structure_db(db_path)

    structure_db = db.DB(db_path, None, ignore_version = True)
    if str(structure_db.get_version()) != current_version:
        structure_db.disconnect()
        raise ConfigError("Version of this structure database is not %s (hence, this script cannot really do anything)." % current_version)

    progress.new("Migrating structure database to v%s" % next_version)

    progress.update("Adding new columns for ColabFold confidence metrics")
    for table_name, columns in new_columns.items():
        for column_name, column_type in columns:
            structure_db._exec('''ALTER TABLE %s ADD COLUMN %s %s''' % (table_name, column_name, column_type))

    progress.update("Recording the structure-prediction engine")
    # all pre-v3 databases were created with MODELLER (external structures were also stored through
    # the MODELLER code path), so we tag them as such
    structure_db.set_meta_value('engine', 'modeller')

    progress.update("Updating the version")
    structure_db.remove_meta_key_value_pair('version')
    structure_db.set_version(next_version)

    structure_db.disconnect()
    progress.end()

    run.info_single("Your structure db is now version %s. It gained a few columns to store ColabFold "
                    "confidence metrics, and now records which engine was used to predict its "
                    "structures." % next_version, nl_after=1, nl_before=1, mc='green')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade structure database from version %s to version %s' % (current_version, next_version))
    parser.add_argument('structure_db', metavar = 'STRUCTURE_DB', help = "An anvi'o structure database of version %s" % current_version)
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.structure_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
