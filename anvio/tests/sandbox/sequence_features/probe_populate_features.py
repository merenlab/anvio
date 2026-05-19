#!/usr/bin/env python
"""End-to-end probe for TablesForSequenceFeatures.populate_features and the new
db.DB.transaction() context manager. Used by run_component_tests_for_sequence_features.sh
in lieu of full-pipeline integration testing for commit 3, since the GenBank importer
that would normally drive populate_features does not yet exist at this commit.
"""

import sys
import sqlite3

import anvio.db as db
import anvio.tables as t
import anvio.terminal as terminal
import anvio.utils as utils
from anvio.tables.sequencefeatures import TablesForSequenceFeatures
from anvio.errors import ConfigError

run = terminal.Run(verbose=False)
progress = terminal.Progress(verbose=False)


def _assert(cond, label):
    if not cond:
        print(f"  FAIL [{label}]")
        sys.exit(1)
    print(f"  OK   [{label}]")


def feature(feature_id, contig, ftype, source, start, stop, direction='f', fg=None, so=None, ext=None, derivation=None, derived_from=None):
    return {
        'feature_id': feature_id, 'contig': contig, 'feature_type': ftype, 'source': source,
        'start': start, 'stop': stop, 'direction': direction,
        'partial_fiveprime': 0, 'partial_threeprime': 0,
        'feature_group_id': fg, 'segment_order': so,
        'external_id': ext, 'gene_callers_id': None,
        'derivation': derivation, 'derived_from_feature_id': derived_from,
    }


def count(db_path, table, where=''):
    con = sqlite3.connect(db_path)
    try:
        return con.execute(f"SELECT COUNT(*) FROM {table} {where}").fetchone()[0]
    finally:
        con.close()


def main():
    if len(sys.argv) != 2:
        print("usage: probe_populate_features.py <contigs-db-path>"); sys.exit(2)
    db_path = sys.argv[1]

    # discover a real contig name we can use for feature.contig
    con = sqlite3.connect(db_path)
    contig_a = con.execute("SELECT contig FROM contigs_basic_info LIMIT 1").fetchone()[0]
    con.close()

    tfsf = TablesForSequenceFeatures(db_path, run=run, progress=progress)

    # -------- first insert: builtins + a brand-new (non-builtin) feature_type --------
    feats_a = [
        feature('A_gene_1',  contig_a, 'gene', 'src_alpha', 100, 500),
        feature('A_cds_1a',  contig_a, 'CDS',  'src_alpha', 100, 200, fg='A_cds_1a', so=0),
        feature('A_cds_1b',  contig_a, 'CDS',  'src_alpha', 300, 400, fg='A_cds_1a', so=1),
        feature('A_reg_1',   contig_a, 'regulatory', 'src_alpha', 600, 700),  # non-builtin
    ]
    rels_a = [
        {'child_feature_id': 'A_cds_1a', 'parent_feature_id': 'A_gene_1', 'relationship': 'derives_from'},
        {'child_feature_id': 'A_cds_1b', 'parent_feature_id': 'A_gene_1', 'relationship': 'derives_from'},
    ]
    quals_a = [
        {'feature_id': 'A_gene_1', 'key': 'locus_tag', 'value': 'LT_A', 'position': 0},
        {'feature_id': 'A_gene_1', 'key': 'db_xref',   'value': 'ref:1', 'position': 0},
        {'feature_id': 'A_gene_1', 'key': 'db_xref',   'value': 'ref:2', 'position': 1},
    ]
    cds_a = [
        {'feature_id': 'A_cds_1a', 'codon_start': 1,    'translation': 'MAB', 'transl_table': 11},
        {'feature_id': 'A_cds_1b', 'codon_start': None, 'translation': None, 'transl_table': None},
    ]

    tfsf.populate_features(feats_a, rels_a, quals_a, cds_a, source_name='src_alpha', force=False)

    _assert(count(db_path, 'contigs_sequence_features') == 4, "4 rows in contigs_sequence_features after first import")
    _assert(count(db_path, 'feature_relationships') == 2, "2 rows in feature_relationships")
    _assert(count(db_path, 'feature_qualifiers') == 3, "3 rows in feature_qualifiers")
    _assert(count(db_path, 'CDS_features') == 2, "2 rows in CDS_features")
    _assert(count(db_path, 'feature_types', "WHERE feature_type='regulatory'") == 1, "non-builtin 'regulatory' was registered in feature_types")
    _assert(count(db_path, 'feature_types', "WHERE feature_type='regulatory' AND is_builtin=0 AND has_dedicated_table=0") == 1, "'regulatory' has is_builtin=0 and has_dedicated_table=0")

    # -------- second insert: different source coexists --------
    feats_b = [
        feature('B_gene_1', contig_a, 'gene', 'src_beta', 1000, 1500),
    ]
    tfsf.populate_features(feats_b, [], [], [], source_name='src_beta', force=False)
    _assert(count(db_path, 'contigs_sequence_features') == 5, "5 rows total after adding second source")
    _assert(count(db_path, 'contigs_sequence_features', "WHERE source='src_alpha'") == 4, "src_alpha rows untouched by src_beta import")

    # -------- force=True on src_alpha: only alpha rows are removed and re-inserted --------
    # rerun the alpha set with force=True; expect alpha row count to match the same 4, and src_beta to remain
    tfsf.populate_features(feats_a, rels_a, quals_a, cds_a, source_name='src_alpha', force=True)
    _assert(count(db_path, 'contigs_sequence_features', "WHERE source='src_alpha'") == 4, "force=True replaced src_alpha rows (count is 4 again)")
    _assert(count(db_path, 'contigs_sequence_features', "WHERE source='src_beta'") == 1, "force=True on src_alpha left src_beta untouched")
    _assert(count(db_path, 'feature_relationships') == 2, "feature_relationships rebuilt for src_alpha (count is 2)")
    _assert(count(db_path, 'feature_qualifiers') == 3, "feature_qualifiers rebuilt for src_alpha (count is 3)")
    _assert(count(db_path, 'CDS_features') == 2, "CDS_features rebuilt for src_alpha (count is 2)")

    # -------- transaction rollback: trigger PK violation mid-batch and confirm nothing landed --------
    pre_count = count(db_path, 'contigs_sequence_features')
    bad_feats = [
        feature('X_new_1', contig_a, 'gene', 'src_gamma', 5000, 5100),
        # the next row collides with the FIRST in the same batch — should trigger UNIQUE INDEX violation
        feature('X_new_1', contig_a, 'gene', 'src_gamma', 5200, 5300),
    ]
    try:
        tfsf.populate_features(bad_feats, [], [], [], source_name='src_gamma', force=False)
        _assert(False, "expected a transaction rollback on duplicate feature_id")
    except sqlite3.IntegrityError:
        pass  # exactly what we want — the transaction context manager re-raised after rollback
    except Exception as e:
        _assert(False, f"unexpected exception type from transaction rollback: {type(e).__name__}: {e}")
    post_count = count(db_path, 'contigs_sequence_features')
    _assert(pre_count == post_count, f"rolled-back insert left contigs_sequence_features unchanged ({pre_count} == {post_count})")
    _assert(count(db_path, 'contigs_sequence_features', "WHERE source='src_gamma'") == 0, "no src_gamma rows landed after rollback")

    # -------- nested transaction is rejected --------
    anvio_db = db.DB(db_path, utils.get_required_version_for_db(db_path))
    try:
        with anvio_db.transaction():
            try:
                with anvio_db.transaction():
                    _assert(False, "expected ConfigError for nested transaction")
            except ConfigError:
                pass
    finally:
        anvio_db.disconnect()
    print("  OK   [nested transaction raises ConfigError]")

    print("All populate_features / transaction probes PASSED.")


if __name__ == '__main__':
    main()
