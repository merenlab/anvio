# -*- coding: utf-8
"""Schema creation and write API for the contigs-database sequence-features tables.

This module is the single source of truth for the SQL schema of the five tables
introduced in contigs-db version 25:

    - contigs_sequence_features
    - feature_types
    - feature_relationships
    - feature_qualifiers
    - CDS_features

`create_sequence_features_tables(db)` is called both from
`anvio.dbops.ContigsDatabase.touch` (when a fresh contigs database is being
generated) and from `anvio.migrations.contigs.v24_to_v25` (when an existing
v24 database is being migrated). Using one function for both paths guarantees
that fresh and migrated databases end up with identical table definitions,
identical indexes, and an identical builtin `feature_types` registry.

The `TablesForSequenceFeatures` class (added in a later commit) is the public
bulk-write entry point used by `anvi-import-genbank-features`.
"""

import anvio
import anvio.db as db
import anvio.utils as utils
import anvio.tables as t
import anvio.terminal as terminal


__copyright__ = "Copyleft 2015-2026, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['semiller10']
__status__ = "Development"


run = terminal.Run()
progress = terminal.Progress()


BUILTIN_FEATURE_TYPES = [
    ('gene',   1, 0, 'A region transcribed as a unit; SO:0000704'),
    ('mRNA',   1, 0, 'Mature messenger RNA; SO:0000234'),
    ('CDS',    1, 1, 'Coding sequence; SO:0000316'),
    ('exon',   1, 0, 'A region of a transcript that remains after splicing; SO:0000147'),
    ('intron', 1, 0, 'A region removed during splicing; SO:0000188'),
]


def create_sequence_features_tables(db):
    """Create the five sequence-features tables, their indexes, and populate `feature_types`.

    Idempotent at the index level (uses `CREATE INDEX IF NOT EXISTS`) and at the
    `feature_types` insert (uses `INSERT OR IGNORE`), but expects the tables
    themselves not to pre-exist on the database — the caller is responsible for
    invoking this only against a database where the five tables have not yet been
    created (either a fresh database during `touch`, or a v24 database being
    upgraded to v25 by the migration script).

    Parameters
    ==========
    db : anvio.db.DB
        An open `anvio.db.DB` connection. The caller manages opening, version
        handling, and disconnection.
    """

    db.create_table(t.contigs_sequence_features_table_name, t.contigs_sequence_features_table_structure, t.contigs_sequence_features_table_types)
    db.create_table(t.feature_types_table_name,             t.feature_types_table_structure,             t.feature_types_table_types)
    db.create_table(t.feature_relationships_table_name,     t.feature_relationships_table_structure,     t.feature_relationships_table_types)
    db.create_table(t.feature_qualifiers_table_name,        t.feature_qualifiers_table_structure,        t.feature_qualifiers_table_types)
    db.create_table(t.CDS_features_table_name,              t.CDS_features_table_structure,              t.CDS_features_table_types)

    # primary-key uniqueness — enforced via UNIQUE INDEX since anvi'o's
    # triplet-driven create_table() does not declare PRIMARY KEY columns.
    db._exec(f'''CREATE UNIQUE INDEX IF NOT EXISTS idx_pk_{t.contigs_sequence_features_table_name} ON {t.contigs_sequence_features_table_name} (feature_id)''')
    db._exec(f'''CREATE UNIQUE INDEX IF NOT EXISTS idx_pk_{t.feature_types_table_name} ON {t.feature_types_table_name} (feature_type)''')
    db._exec(f'''CREATE UNIQUE INDEX IF NOT EXISTS idx_pk_{t.feature_relationships_table_name} ON {t.feature_relationships_table_name} (child_feature_id, parent_feature_id, relationship)''')
    db._exec(f'''CREATE UNIQUE INDEX IF NOT EXISTS idx_pk_{t.feature_qualifiers_table_name} ON {t.feature_qualifiers_table_name} (feature_id, key, position)''')
    db._exec(f'''CREATE UNIQUE INDEX IF NOT EXISTS idx_pk_{t.CDS_features_table_name} ON {t.CDS_features_table_name} (feature_id)''')

    # secondary indexes on contigs_sequence_features for the dominant query
    # patterns: range scans on contigs, type filters, and joins via the
    # group/GCID columns.
    db._exec(f'''CREATE INDEX IF NOT EXISTS idx_{t.contigs_sequence_features_table_name}_contig_start ON {t.contigs_sequence_features_table_name} (contig, start)''')
    db._exec(f'''CREATE INDEX IF NOT EXISTS idx_{t.contigs_sequence_features_table_name}_contig_stop  ON {t.contigs_sequence_features_table_name} (contig, stop)''')
    db._exec(f'''CREATE INDEX IF NOT EXISTS idx_{t.contigs_sequence_features_table_name}_feature_type ON {t.contigs_sequence_features_table_name} (feature_type)''')
    db._exec(f'''CREATE INDEX IF NOT EXISTS idx_{t.contigs_sequence_features_table_name}_group_id     ON {t.contigs_sequence_features_table_name} (feature_group_id)''')
    db._exec(f'''CREATE INDEX IF NOT EXISTS idx_{t.contigs_sequence_features_table_name}_gcid         ON {t.contigs_sequence_features_table_name} (gene_callers_id)''')

    # parent_feature_id lookups are common ("what are the children of feature X");
    # child_feature_id is already covered by the leading column of the PK.
    db._exec(f'''CREATE INDEX IF NOT EXISTS idx_{t.feature_relationships_table_name}_parent ON {t.feature_relationships_table_name} (parent_feature_id)''')

    # builtin feature types — OR IGNORE keeps this safe under future
    # re-introductions (the migration script and touch() each call this once
    # against a fresh table, so the OR IGNORE is defensive, not load-bearing).
    db._exec_many(f'''INSERT OR IGNORE INTO {t.feature_types_table_name} VALUES (?, ?, ?, ?)''', BUILTIN_FEATURE_TYPES)


class TablesForSequenceFeatures:
    """Bulk-write entry point for the sequence-features tables (contigs-db v25).

    Mirrors the `TablesForGeneCalls` pattern: the constructor takes a contigs-db
    path plus optional `args` / `run` / `progress`; the public surface is the
    single `populate_features` method. All inserts (and any pre-import deletion
    triggered by `force=True`) happen within one transaction so a failed import
    leaves the database unchanged.

    Single-row insert methods are deliberately not exposed; the importer is
    expected to assemble its full in-memory representation in pass 1 and write
    it all in one shot. Adding finer-grained writes is a future-PR concern.
    """

    def __init__(self, db_path, args=None, run=terminal.Run(), progress=terminal.Progress()):
        self.db_path = db_path
        self.args = args
        self.run = run
        self.progress = progress

        utils.is_contigs_db(self.db_path)


    def populate_features(self, features, relationships, qualifiers, cds_specific_data, source_name, force=False):
        """Bulk-write the four data tables in one transaction.

        Parameters
        ==========
        features : list of dict
            One dict per row of `contigs_sequence_features`. Every key in
            `contigs_sequence_features_table_structure` must be present on every
            dict; the caller is responsible for value-domain validation
            (direction in {'f','r',None}, half-open coordinates, etc.).
        relationships : list of dict
            Each dict has keys `child_feature_id`, `parent_feature_id`, `relationship`.
        qualifiers : list of dict
            Each dict has keys `feature_id`, `key`, `value`, `position`.
        cds_specific_data : list of dict
            Each dict has keys `feature_id`, `codon_start`, `translation`, `transl_table`.
            For multi-segment CDSs the caller already nulls these on non-canonical rows.
        source_name : str
            Used both to delete pre-existing rows (when `force=True`) and as the
            documented source of these rows in `contigs_sequence_features.source`.
            The caller is responsible for the regex validation required by the CLI.
        force : bool, default False
            If True, delete every row belonging to `source_name` from the four
            target tables before inserting. If False, the caller has already
            verified no such rows exist (and aborts the import otherwise).

        Notes
        =====
        - Any `feature_type` values present in `features` that are not yet
          registered in `feature_types` are inserted with `is_builtin = 0` and
          `has_dedicated_table = 0`. This precedes the main feature insert so
          the database remains consistent without enforced foreign keys.
        - The transaction is opened against the underlying SQLite connection
          via `DB.transaction()`: every helper call within the `with` block
          defers its commit until the context exits, and any exception inside
          the block rolls everything back.
        """

        anvio_db = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))

        try:
            with anvio_db.transaction():
                if force:
                    self._delete_existing_source(anvio_db, source_name)

                self._register_non_builtin_feature_types(anvio_db, features)

                self._insert_dicts(anvio_db, t.contigs_sequence_features_table_name, t.contigs_sequence_features_table_structure, features)
                self._insert_dicts(anvio_db, t.feature_relationships_table_name,     t.feature_relationships_table_structure,     relationships)
                self._insert_dicts(anvio_db, t.feature_qualifiers_table_name,        t.feature_qualifiers_table_structure,        qualifiers)
                self._insert_dicts(anvio_db, t.CDS_features_table_name,              t.CDS_features_table_structure,              cds_specific_data)
        finally:
            anvio_db.disconnect()


    def _insert_dicts(self, anvio_db, table_name, column_order, dicts):
        """Convert dict-rows to tuples in `column_order` and batch-insert."""

        if not dicts:
            return

        tuples = [tuple(d[col] for col in column_order) for d in dicts]
        anvio_db.insert_many(table_name, tuples)


    def _register_non_builtin_feature_types(self, anvio_db, features):
        """Register any non-builtin feature_type values encountered in `features`.

        New rows are inserted with `is_builtin = 0`, `has_dedicated_table = 0`,
        and `description = NULL`. `INSERT OR IGNORE` makes this idempotent:
        running the importer again with the same non-builtin types is a no-op
        on the registry.
        """

        if not features:
            return

        seen_types = {f['feature_type'] for f in features}

        existing_rows = anvio_db._fetchall(anvio_db._exec(f'''SELECT feature_type FROM {t.feature_types_table_name}'''),
                                           t.feature_types_table_name)
        existing = {row[0] for row in existing_rows}

        missing = seen_types - existing
        if not missing:
            return

        rows_to_insert = [(name, 0, 0, None) for name in sorted(missing)]
        anvio_db._exec_many(f'''INSERT OR IGNORE INTO {t.feature_types_table_name} VALUES (?, ?, ?, ?)''',
                            rows_to_insert)


    def _delete_existing_source(self, anvio_db, source_name):
        """Remove every row associated with `source_name` from the four target tables.

        Only `contigs_sequence_features` carries a `source` column, so we first
        collect the relevant `feature_id`s and then cascade through the dependent
        tables. SQLite's default parameter limit is small (defaults to 999 in the
        wheels we ship against), so the IN clauses are batched.
        """

        rows = anvio_db._fetchall(
            anvio_db._exec(f'''SELECT feature_id FROM {t.contigs_sequence_features_table_name} WHERE source = ?''', (source_name,)),
            t.contigs_sequence_features_table_name,
        )
        feature_ids = [row[0] for row in rows]
        if not feature_ids:
            return

        BATCH_SIZE = 500
        for i in range(0, len(feature_ids), BATCH_SIZE):
            chunk = tuple(feature_ids[i:i + BATCH_SIZE])
            placeholders = ','.join(['?'] * len(chunk))

            anvio_db._exec(f'''DELETE FROM {t.feature_relationships_table_name}     WHERE child_feature_id  IN ({placeholders})''', chunk)
            anvio_db._exec(f'''DELETE FROM {t.feature_relationships_table_name}     WHERE parent_feature_id IN ({placeholders})''', chunk)
            anvio_db._exec(f'''DELETE FROM {t.feature_qualifiers_table_name}        WHERE feature_id        IN ({placeholders})''', chunk)
            anvio_db._exec(f'''DELETE FROM {t.CDS_features_table_name}              WHERE feature_id        IN ({placeholders})''', chunk)
            anvio_db._exec(f'''DELETE FROM {t.contigs_sequence_features_table_name} WHERE feature_id        IN ({placeholders})''', chunk)
