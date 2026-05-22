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

new_table_name      = 'gene_level_normalized_coverages'
new_table_structure = ['gene_callers_id', 'sample_name', 'log1p', 'rpm', 'zscore_raw', 'zscore_log1p', 'zscore_rpm']
new_table_types     = ['numeric', 'text', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric']


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    utils.is_genes_db(db_path)

    genes_db = db.DB(db_path, None, ignore_version=True)
    if str(genes_db.get_version()) != current_version:
        genes_db.disconnect()
        raise ConfigError("Version of this genes database is not %s (hence, this script cannot really do anything)." % current_version)
    genes_db.disconnect()

    progress.new("Upgrading genes database")
    progress.update("Adding gene_level_normalized_coverages table...")

    genes_db = db.DB(db_path, None, ignore_version=True)

    db_fields = ', '.join(['%s %s' % (col, typ) for col, typ in zip(new_table_structure, new_table_types)])
    genes_db._exec(f'''CREATE TABLE {new_table_name} ({db_fields})''')

    genes_db.set_meta_value('gene_level_normalized_coverages_stored', False)

    genes_db.remove_meta_key_value_pair('version')
    genes_db.set_version(next_version)
    genes_db.disconnect()

    progress.end()

    run.info_single(f"Your genes db is now version {next_version}. A new empty table "
                    f"'gene_level_normalized_coverages' has been added to store normalized "
                    f"gene coverage values (log1p, rpm, zscore_raw, zscore_log1p, zscore_rpm). "
                    f"Run anvi-interactive with --gene-mode --compute-gene-level-normalized-coverages "
                    f"to populate it.", nl_after=1, nl_before=1, mc='green')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade genes database from version %s to version %s' % (current_version, next_version))
    parser.add_argument('genes_db', metavar='GENES_DB', help="An anvi'o genes database of version %s" % current_version)
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.genes_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
