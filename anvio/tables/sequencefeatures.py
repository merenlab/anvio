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
