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

# in v4 the per-protein key column `corresponding_gene_call` is renamed to `protein_id` across all four
# structure tables (it is a surrogate key now, equal to the gene caller id for a contigs-db input but an
# opaque id for other input types), and a new `proteins` table records each protein's provenance. These
# are the v4 (post-rename) column definitions for the four tables. SQLite cannot rename a column in place
# on the versions anvi'o supports, so each table is rebuilt with the create-temp / copy / drop / rename
# recipe. The copy selects columns BY NAME (not position) because a v3 database created fresh and one
# migrated up from v2 can have the ColabFold columns (plddt, mean_plddt, ptm) in different physical
# positions -- a positional copy would silently misalign them.
tables = {
    'structures':   ['protein_id', 'pdb_content'],
    'templates':    ['protein_id', 'pdb_id', 'chain_id', 'proper_percent_similarity', 'percent_similarity', 'align_fraction'],
    'models':       ['protein_id', 'molpdf', 'GA341_score', 'DOPE_score', 'mean_plddt', 'ptm', 'picked_as_best'],
    'residue_info': ['protein_id', 'codon_order_in_gene', 'contact_numbers', 'codon', 'amino_acid', 'codon_number', 'sec_struct', 'rel_solvent_acc', 'phi', 'psi', 'plddt', 'NH_O_1_index', 'NH_O_1_energy', 'O_NH_1_index', 'O_NH_1_energy', 'NH_O_2_index', 'NH_O_2_energy', 'O_NH_2_index', 'O_NH_2_energy'],
}
table_types = {
    'structures':   ['integer', 'blob'],
    'templates':    ['integer', 'text', 'text', 'real', 'real', 'real'],
    'models':       ['integer', 'real', 'real', 'real', 'real', 'real', 'integer'],
    'residue_info': ['integer', 'integer', 'text', 'text', 'text', 'integer', 'text', 'real', 'real', 'real', 'real', 'integer', 'real', 'integer', 'real', 'integer', 'real', 'integer', 'real'],
}

# the new provenance table
proteins_structure = ['protein_id', 'input_type', 'genome_name', 'gene_callers_id', 'gene_cluster_id', 'source_key', 'has_structure']
proteins_types     = ['integer', 'text', 'text', 'integer', 'text', 'text', 'integer']


def rename_key_column(db_path, table_name):
    """Rebuild a table with its key column renamed from `corresponding_gene_call` to `protein_id`."""

    progress.new("Renaming the key column in '%s'" % table_name)

    new_structure = tables[table_name]
    types = table_types[table_name]
    db_fields = ', '.join(['%s %s' % (name, column_type) for name, column_type in zip(new_structure, types)])
    temp_table_name = table_name + '_TEMP'

    # the source columns, in the same order as new_structure but with the old name for the key column, so
    # the copy matches columns by name regardless of their physical order in the existing table
    source_columns = ['corresponding_gene_call'] + new_structure[1:]

    _db = db.DB(db_path, None, ignore_version = True)

    progress.update("Creating a temporary table")
    _db._exec('''CREATE TABLE %s (%s)''' % (temp_table_name, db_fields))

    progress.update("Copying data into the temporary table")
    _db._exec('''INSERT INTO %s (%s) SELECT %s FROM %s''' % (temp_table_name,
                                                             ', '.join(new_structure),
                                                             ', '.join(source_columns),
                                                             table_name))

    progress.update("Dropping the original table")
    _db._exec('''DROP TABLE IF EXISTS %s''' % (table_name))

    progress.update("Renaming temporary table to become the original")
    _db._exec('''ALTER TABLE %s RENAME TO %s''' % (temp_table_name, table_name))

    _db.disconnect()
    progress.end()


def parse_gene_id_csv(value):
    """Parse a comma-joined string of gene caller ids (the pre-v4 self-table bookkeeping) into a list of ints."""

    if not value:
        return []

    return [int(x) for x in value.split(',') if not x == '']


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    # make sure someone is not being funny
    utils.is_structure_db(db_path)

    structure_db = db.DB(db_path, None, ignore_version = True)
    if str(structure_db.get_version()) != current_version:
        structure_db.disconnect()
        raise ConfigError("Version of this structure database is not %s (hence, this script cannot really do anything)." % current_version)

    # read the pre-v4 bookkeeping (comma-joined gene caller ids in the self table) before we touch anything,
    # so we can rebuild it as rows in the new proteins table
    genes_with_structure = parse_gene_id_csv(structure_db.get_meta_value('genes_with_structure', try_as_type_int=False, return_none_if_not_in_table=True))
    genes_without_structure = parse_gene_id_csv(structure_db.get_meta_value('genes_without_structure', try_as_type_int=False, return_none_if_not_in_table=True))
    structure_db.disconnect()

    progress.new("Migrating structure database to v%s" % next_version)
    progress.update("...")
    progress.end()

    # rename the key column in each of the four structure tables
    for table_name in tables:
        rename_key_column(db_path, table_name)

    structure_db = db.DB(db_path, None, ignore_version = True)

    progress.new("Migrating structure database to v%s" % next_version)

    # create the new proteins table and enforce protein_id uniqueness with the same index the live code
    # creates for a freshly-made database, so migrated and fresh databases behave identically
    progress.update("Creating the proteins table")
    proteins_db_fields = ', '.join(['%s %s' % (name, column_type) for name, column_type in zip(proteins_structure, proteins_types)])
    structure_db._exec('''CREATE TABLE proteins (%s)''' % proteins_db_fields)
    structure_db._exec('''CREATE UNIQUE INDEX IF NOT EXISTS proteins_protein_id_idx ON proteins(protein_id)''')

    # backfill one row per queried protein. This database was necessarily made from a contigs database
    # (the only input type before v4), so protein_id == gene_callers_id and the other identity columns are
    # left NULL
    progress.update("Backfilling protein provenance")
    entries = []
    for gene_callers_id in genes_with_structure:
        entries.append((gene_callers_id, 'contigs_db', None, gene_callers_id, None, str(gene_callers_id), 1))
    for gene_callers_id in genes_without_structure:
        entries.append((gene_callers_id, 'contigs_db', None, gene_callers_id, None, str(gene_callers_id), 0))

    if entries:
        structure_db._exec_many('''INSERT INTO proteins VALUES (%s)''' % ','.join(['?'] * len(proteins_structure)), entries)

    # record the input type and drop the now-obsolete comma-joined bookkeeping from the self table
    progress.update("Updating the self table")
    structure_db.set_meta_value('input_type', 'contigs_db')
    for key in ['genes_queried', 'genes_with_structure', 'genes_without_structure']:
        structure_db.remove_meta_key_value_pair(key)

    progress.update("Updating the version")
    structure_db.remove_meta_key_value_pair('version')
    structure_db.set_version(next_version)

    structure_db.disconnect()
    progress.end()

    run.info_single("Your structure db is now version %s. Its per-protein key column was renamed to "
                    "'protein_id', and a new 'proteins' table now records where each protein came from "
                    "(all of them from your contigs database, for a database of this vintage)." % next_version,
                    nl_after=1, nl_before=1, mc='green')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade structure database from version %s to version %s' % (current_version, next_version))
    parser.add_argument('structure_db', metavar = 'STRUCTURE_DB', help = "An anvi'o structure database of version %s" % current_version)
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.structure_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
