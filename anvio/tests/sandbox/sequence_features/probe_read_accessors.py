#!/usr/bin/env python
"""Hand-populate the sequence-features tables of a v25 contigs database and probe each
read accessor on ContigsSuperclass.

Used by run_component_tests_for_sequence_features.sh in lieu of full-pipeline integration
testing for commit 2, since the GenBank importer that would normally exercise these
accessors does not yet exist at this commit.
"""

import sys
import sqlite3
from argparse import Namespace

import anvio.tables as t
import anvio.terminal as terminal
from anvio.dbops import ContigsSuperclass

run = terminal.Run(verbose=False)
progress = terminal.Progress(verbose=False)


def _assert(cond, label):
    if not cond:
        print(f"  FAIL [{label}]")
        sys.exit(1)
    print(f"  OK   [{label}]")


def hand_populate(db_path, contig_a, contig_b):
    """Insert deterministic test rows into the five new sequence-features tables."""

    feats = [
        # single-segment gene on contig A
        ('feat_gene_ss',      contig_a, 'gene', 'probe', 100,  500,  'f', 0, 0, None,            None, 'GENE_SS', None,    None, None),
        # multi-segment CDS on contig A, forward (three segments, transcription = genomic order)
        ('feat_cds_3seg_a',   contig_a, 'CDS',  'probe', 100,  200,  'f', 0, 0, 'feat_cds_3seg_a', 0,   'CDS_3SEG', None,  None, None),
        ('feat_cds_3seg_b',   contig_a, 'CDS',  'probe', 300,  400,  'f', 0, 0, 'feat_cds_3seg_a', 1,   'CDS_3SEG', None,  None, None),
        ('feat_cds_3seg_c',   contig_a, 'CDS',  'probe', 500,  600,  'f', 0, 0, 'feat_cds_3seg_a', 2,   'CDS_3SEG', None,  None, None),
        # multi-segment mRNA on contig A, forward (parent of the CDS above)
        ('feat_mrna_2seg_a',  contig_a, 'mRNA', 'probe', 100,  200,  'f', 0, 0, 'feat_mrna_2seg_a', 0,  'MRNA_2SEG', None, None, None),
        ('feat_mrna_2seg_b',  contig_a, 'mRNA', 'probe', 300,  600,  'f', 0, 0, 'feat_mrna_2seg_a', 1,  'MRNA_2SEG', None, None, None),
        # multi-segment CDS on contig A, reverse — segment_order=0 must be the genomic rightmost
        ('feat_cds_rev_a',    contig_a, 'CDS',  'probe', 900,  1000, 'r', 0, 0, 'feat_cds_rev_a',  0,   'CDS_REV',  None, None, None),
        ('feat_cds_rev_b',    contig_a, 'CDS',  'probe', 700,  800,  'r', 0, 0, 'feat_cds_rev_a',  1,   'CDS_REV',  None, None, None),
        # orphan single-segment intron (no parent expected) on contig A
        ('feat_intron_orphan', contig_a, 'intron', 'probe', 2000, 2500, 'f', 0, 0, None,           None, None,     None,  None, None),
        # single-segment feature on contig B (used to confirm contig filtering)
        ('feat_gene_b',       contig_b, 'gene', 'probe', 50,   400,  'f', 0, 0, None,            None, None,      None,   None, None),
    ]

    rels = [
        # CDS three segments all link to canonical mRNA row (part_of)
        ('feat_cds_3seg_a', 'feat_mrna_2seg_a', 'part_of'),
        ('feat_cds_3seg_b', 'feat_mrna_2seg_a', 'part_of'),
        ('feat_cds_3seg_c', 'feat_mrna_2seg_a', 'part_of'),
        # mRNA two segments link to gene (part_of)
        ('feat_mrna_2seg_a', 'feat_gene_ss', 'part_of'),
        ('feat_mrna_2seg_b', 'feat_gene_ss', 'part_of'),
        # reverse CDS two segments link to the gene with derives_from (no mRNA in between)
        ('feat_cds_rev_a', 'feat_gene_ss', 'derives_from'),
        ('feat_cds_rev_b', 'feat_gene_ss', 'derives_from'),
    ]

    quals = [
        ('feat_gene_ss', 'product',   'Test gene',   0),
        ('feat_gene_ss', 'locus_tag', 'TST_001',     0),
        ('feat_gene_ss', 'db_xref',   'ref:first',   0),
        ('feat_gene_ss', 'db_xref',   'ref:second',  1),
        ('feat_cds_3seg_a', 'product', 'Test CDS',   0),
    ]

    cds = [
        ('feat_cds_3seg_a', 1, 'MAKE',  11),
        ('feat_cds_3seg_b', None, None, None),
        ('feat_cds_3seg_c', None, None, None),
        ('feat_cds_rev_a',  1, 'REVS',  11),
        ('feat_cds_rev_b',  None, None, None),
    ]

    con = sqlite3.connect(db_path)
    cur = con.cursor()
    cur.executemany(f"INSERT INTO {t.contigs_sequence_features_table_name} VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", feats)
    cur.executemany(f"INSERT INTO {t.feature_relationships_table_name}      VALUES (?,?,?)", rels)
    cur.executemany(f"INSERT INTO {t.feature_qualifiers_table_name}         VALUES (?,?,?,?)", quals)
    cur.executemany(f"INSERT INTO {t.CDS_features_table_name}               VALUES (?,?,?,?)", cds)
    con.commit()
    con.close()


def main():
    if len(sys.argv) != 2:
        print("usage: probe_read_accessors.py <contigs-db-path>")
        sys.exit(2)
    db_path = sys.argv[1]

    # pick two contig names from the DB so the test data correlates with reality
    con = sqlite3.connect(db_path)
    contigs = [r[0] for r in con.execute("SELECT contig FROM contigs_basic_info ORDER BY contig LIMIT 2")]
    con.close()
    if len(contigs) < 2:
        print("Test db needs at least two contigs in contigs_basic_info."); sys.exit(2)
    contig_a, contig_b = contigs

    hand_populate(db_path, contig_a, contig_b)

    csc = ContigsSuperclass(Namespace(contigs_db=db_path), r=run, p=progress)

    # LazyProperty accessors
    sfd = csc.sequence_features_dict
    _assert(isinstance(sfd, dict) and len(sfd) == 10, f"sequence_features_dict has 10 entries (got {len(sfd)})")
    _assert('feat_gene_ss' in sfd, "sequence_features_dict contains feat_gene_ss")

    ftd = csc.feature_types_dict
    _assert(isinstance(ftd, dict) and len(ftd) >= 5, f"feature_types_dict has >=5 entries (got {len(ftd)})")
    _assert(ftd.get('CDS', {}).get('has_dedicated_table') == 1, "feature_types_dict knows CDS has a dedicated table")

    # get_sequence_features_in_contig — contig A, no type filter
    a_all = csc.get_sequence_features_in_contig(contig_a)
    _assert(len(a_all) == 9, f"get_sequence_features_in_contig(A) returns 9 rows (got {len(a_all)})")
    _assert('feat_gene_b' not in a_all, "contig filter excludes features on contig B")
    _assert('qualifiers' in a_all['feat_gene_ss'], "result rows carry a 'qualifiers' sub-dict")
    _assert(a_all['feat_gene_ss']['qualifiers']['db_xref'] == ['ref:first', 'ref:second'], "qualifier values come back in position order")
    _assert(a_all['feat_gene_ss']['qualifiers']['product'] == ['Test gene'], "single-value qualifier comes back as a one-element list")

    # contig A, type=CDS
    a_cds = csc.get_sequence_features_in_contig(contig_a, feature_type='CDS')
    _assert(len(a_cds) == 5, f"contig A has 5 CDS rows (3 + 2 segments) (got {len(a_cds)})")
    _assert(all(r['feature_type'] == 'CDS' for r in a_cds.values()), "feature_type filter is honored")

    # contig B has only a single gene
    b_all = csc.get_sequence_features_in_contig(contig_b)
    _assert(set(b_all.keys()) == {'feat_gene_b'}, "contig B contains only feat_gene_b")

    # contig with no features returns empty
    empty = csc.get_sequence_features_in_contig('__nonexistent_contig__')
    _assert(empty == {}, "absent contig returns empty dict")

    # overlap query — half-open [start, stop)
    # feat_gene_ss is [100, 500). Query [300, 400) overlaps. Query [500, 600) does NOT (boundary).
    hit = csc.get_sequence_features_overlapping_range(contig_a, 300, 400)
    _assert('feat_gene_ss' in hit, "overlap query catches feat_gene_ss for [300,400)")
    no_hit = csc.get_sequence_features_overlapping_range(contig_a, 500, 600)
    _assert('feat_gene_ss' not in no_hit, "overlap query excludes feat_gene_ss for [500,600) (half-open)")
    # but [499, 500) should still touch feat_gene_ss because the feature stops at 500 exclusive
    hit_edge = csc.get_sequence_features_overlapping_range(contig_a, 499, 500)
    _assert('feat_gene_ss' in hit_edge, "overlap query includes feat_gene_ss for [499,500)")

    # overlap query with feature_type filter
    cds_in_range = csc.get_sequence_features_overlapping_range(contig_a, 0, 700, feature_type='CDS')
    _assert(len(cds_in_range) == 3, f"three forward-CDS segments overlap [0,700) (got {len(cds_in_range)})")

    # get_children_of_sequence_feature — children of the mRNA canonical row
    children = csc.get_children_of_sequence_feature('feat_mrna_2seg_a')
    _assert(set(children.keys()) == {'feat_cds_3seg_a', 'feat_cds_3seg_b', 'feat_cds_3seg_c'}, "mRNA has three CDS segment children")

    # auto-resolution: passing a non-canonical segment should resolve and return the same children
    children_via_seg1 = csc.get_children_of_sequence_feature('feat_mrna_2seg_b')
    _assert(set(children_via_seg1.keys()) == {'feat_cds_3seg_a', 'feat_cds_3seg_b', 'feat_cds_3seg_c'}, "passing a non-canonical segment resolves to same children")

    # filter by relationship
    derives = csc.get_children_of_sequence_feature('feat_gene_ss', relationship='derives_from')
    _assert(set(derives.keys()) == {'feat_cds_rev_a', 'feat_cds_rev_b'}, "derives_from filter returns reverse CDS segments")
    part_of = csc.get_children_of_sequence_feature('feat_gene_ss', relationship='part_of')
    _assert(set(part_of.keys()) == {'feat_mrna_2seg_a', 'feat_mrna_2seg_b'}, "part_of filter returns mRNA segments")

    # nonexistent feature → empty
    nothing = csc.get_children_of_sequence_feature('does_not_exist')
    _assert(nothing == {}, "absent feature_id returns empty children dict")

    # get_parents_of_sequence_feature — works for any segment
    parents_from_seg1 = csc.get_parents_of_sequence_feature('feat_cds_3seg_b')
    _assert(set(parents_from_seg1.keys()) == {'feat_mrna_2seg_a'}, "parents of CDS seg1 is canonical mRNA")
    parents_from_seg2 = csc.get_parents_of_sequence_feature('feat_cds_3seg_c')
    _assert(set(parents_from_seg2.keys()) == {'feat_mrna_2seg_a'}, "parents of CDS seg2 is canonical mRNA")
    # mRNA's parent is the gene
    parents_mrna = csc.get_parents_of_sequence_feature('feat_mrna_2seg_a')
    _assert(set(parents_mrna.keys()) == {'feat_gene_ss'}, "parents of mRNA seg0 is the gene")
    # orphan has no parents
    orphan_parents = csc.get_parents_of_sequence_feature('feat_intron_orphan')
    _assert(orphan_parents == {}, "orphan intron has no parents")

    # get_sequence_feature_group — ordered list, transcription order
    forward_group = csc.get_sequence_feature_group('feat_cds_3seg_b')
    _assert([row['feature_id'] for row in forward_group] == ['feat_cds_3seg_a', 'feat_cds_3seg_b', 'feat_cds_3seg_c'], "forward CDS group returned in transcription order")
    # reverse CDS: segment_order=0 must be the genomic RIGHTMOST
    rev_group = csc.get_sequence_feature_group('feat_cds_rev_b')
    _assert([row['feature_id'] for row in rev_group] == ['feat_cds_rev_a', 'feat_cds_rev_b'], "reverse CDS group: segment_order=0 is genomic rightmost (5' in transcription)")
    _assert(rev_group[0]['start'] == 900 and rev_group[1]['start'] == 700, "reverse CDS group coordinates are right")
    # single-segment feature
    solo = csc.get_sequence_feature_group('feat_gene_ss')
    _assert(len(solo) == 1 and solo[0]['feature_id'] == 'feat_gene_ss', "single-segment group has one row")

    # get_sequence_feature_qualifiers
    qs = csc.get_sequence_feature_qualifiers('feat_gene_ss')
    _assert(qs == {'product': ['Test gene'], 'locus_tag': ['TST_001'], 'db_xref': ['ref:first', 'ref:second']}, "qualifier dict round-trips correctly")
    _assert(csc.get_sequence_feature_qualifiers('feat_intron_orphan') == {}, "feature without qualifiers returns empty dict")

    # get_cds_translation_for_feature — canonical resolution
    _assert(csc.get_cds_translation_for_feature('feat_cds_3seg_a') == 'MAKE', "canonical CDS returns translation directly")
    _assert(csc.get_cds_translation_for_feature('feat_cds_3seg_c') == 'MAKE', "non-canonical CDS segment auto-resolves to seg0 translation")
    _assert(csc.get_cds_translation_for_feature('feat_cds_rev_b') == 'REVS', "reverse-strand non-canonical resolves correctly too")
    _assert(csc.get_cds_translation_for_feature('feat_gene_ss') is None, "non-CDS feature returns None for translation")
    _assert(csc.get_cds_translation_for_feature('does_not_exist') is None, "absent feature_id returns None")

    print("All sequence-features read accessors PASSED.")


if __name__ == '__main__':
    main()
