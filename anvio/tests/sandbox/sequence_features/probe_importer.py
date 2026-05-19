#!/usr/bin/env python
"""Parsing-correctness probes for GenbankFeatureImporter. The shell test invokes this
to assert tricky behaviour (segment ordering on both strands, origin-crossing, canonical
parent relationships) that is hard to express in pure sqlite3 assertions. These probes
complement (rather than replace) the full integration tests that arrive in commit 5.
"""

import os
import sys
import sqlite3
import subprocess
from argparse import Namespace

HERE = os.path.dirname(os.path.abspath(__file__))

import anvio.terminal as terminal
from anvio.genbankimporter import GenbankFeatureImporter

run = terminal.Run(verbose=False)
progress = terminal.Progress(verbose=False)


def _assert(cond, label):
    if not cond:
        print(f"  FAIL [{label}]"); sys.exit(1)
    print(f"  OK   [{label}]")


def build_contigs_db(fasta, db_out):
    """Invoke anvi-gen-contigs-database --skip-gene-calling to produce a tiny v25 db."""

    subprocess.run([
        'anvi-gen-contigs-database', '-f', fasta, '-o', db_out,
        '--project-name', 'ImporterProbe',
        '--skip-gene-calling', '-L', '1000', '--no-progress',
    ], check=True, capture_output=True)


def import_into(db_path, gb_path, source_name):
    args = Namespace(contigs_db=db_path, input_genbank=gb_path, source_name=source_name)
    GenbankFeatureImporter(args, run=run, progress=progress).process()


def query_one(db, sql, params=()):
    con = sqlite3.connect(db)
    try:
        row = con.execute(sql, params).fetchone()
        return row
    finally:
        con.close()


def query_all(db, sql, params=()):
    con = sqlite3.connect(db)
    try:
        return con.execute(sql, params).fetchall()
    finally:
        con.close()


def probe_bacterial(work):
    db = os.path.join(work, 'BACT.db')
    build_contigs_db(os.path.join(HERE, 'bacterial.fa'), db)
    import_into(db, os.path.join(HERE, 'bacterial.gb'), 'bact_probe')

    total = query_one(db, "SELECT COUNT(*) FROM contigs_sequence_features")[0]
    _assert(total == 8, f"bacterial import inserted 8 features (got {total})")

    n_gene = query_one(db, "SELECT COUNT(*) FROM contigs_sequence_features WHERE feature_type='gene'")[0]
    n_cds  = query_one(db, "SELECT COUNT(*) FROM contigs_sequence_features WHERE feature_type='CDS'")[0]
    _assert(n_gene == 4 and n_cds == 4, f"bacterial import has 4 genes + 4 CDS (got {n_gene} / {n_cds})")

    # CDS → gene derives_from
    derives = query_one(db, "SELECT COUNT(*) FROM feature_relationships WHERE relationship='derives_from'")[0]
    _assert(derives == 4, f"4 CDS→gene derives_from relationships (got {derives})")

    # locus_tag round-trips into external_id
    ext_count = query_one(db, "SELECT COUNT(*) FROM contigs_sequence_features WHERE external_id IS NOT NULL")[0]
    _assert(ext_count == 8, f"every bacterial feature has external_id populated from locus_tag (got {ext_count})")

    # translation stored on every CDS row (all single-segment so all are canonical)
    trans_count = query_one(db, "SELECT COUNT(*) FROM CDS_features WHERE translation IS NOT NULL")[0]
    _assert(trans_count == 4, f"4 translations stored (got {trans_count})")

    # qualifier ordering
    db_xrefs = query_all(db, "SELECT value FROM feature_qualifiers WHERE key='db_xref' ORDER BY position")
    _assert([r[0] for r in db_xrefs] == ['REF:first', 'REF:second'], f"db_xref preserves position order (got {db_xrefs})")

    # the note qualifier preserves its whitespace
    note = query_one(db, "SELECT value FROM feature_qualifiers WHERE key='note'")[0]
    _assert(note == 'this is a note with whitespace', f"note qualifier preserves whitespace (got {note!r})")


def probe_eukaryotic(work):
    db = os.path.join(work, 'EUK.db')
    build_contigs_db(os.path.join(HERE, 'eukaryotic.fa'), db)
    import_into(db, os.path.join(HERE, 'eukaryotic.gb'), 'euk_probe')

    # forward CDS group: 3 segments, segment_order 0=leftmost, 1, 2=rightmost
    fwd_cds = query_all(db, """SELECT start, stop, segment_order
                               FROM contigs_sequence_features
                               WHERE feature_type='CDS' AND direction='f'
                               ORDER BY segment_order""")
    expected_fwd = [(99, 500, 0), (900, 1500, 1), (2000, 2500, 2)]
    _assert(fwd_cds == expected_fwd, f"forward multi-CDS segment order matches genomic order (got {fwd_cds})")

    # reverse CDS group: segment_order=0 must be the genomic RIGHTMOST
    rev_cds = query_all(db, """SELECT start, stop, segment_order
                               FROM contigs_sequence_features
                               WHERE feature_type='CDS' AND direction='r'
                               ORDER BY segment_order""")
    expected_rev = [(3999, 4500, 0), (2999, 3500, 1)]
    _assert(rev_cds == expected_rev, f"reverse multi-CDS segment_order=0 is genomic rightmost (got {rev_cds})")

    # group canonical equals segment_order=0 row
    canon = query_one(db, """SELECT feature_id FROM contigs_sequence_features
                             WHERE feature_type='CDS' AND direction='f' AND segment_order=0""")[0]
    siblings = query_all(db, """SELECT feature_id, segment_order FROM contigs_sequence_features
                                WHERE feature_type='CDS' AND direction='f' ORDER BY segment_order""")
    _assert(all(r[0] != canon for r in siblings if r[1] != 0), "forward CDS canonical feature_id differs from non-canonical")
    fgids = query_all(db, """SELECT DISTINCT feature_group_id FROM contigs_sequence_features
                             WHERE feature_type='CDS' AND direction='f'""")
    _assert(len(fgids) == 1 and fgids[0][0] == canon, f"every forward CDS segment shares feature_group_id == canonical (got {fgids})")

    # translation only on canonical (segment_order=0) row
    cds_with_trans = query_one(db, """SELECT COUNT(*) FROM CDS_features cf
                                      JOIN contigs_sequence_features csf USING(feature_id)
                                      WHERE cf.translation IS NOT NULL AND csf.feature_type='CDS' AND csf.direction='f'""")[0]
    _assert(cds_with_trans == 1, f"forward CDS group has translation on only 1 row (got {cds_with_trans})")

    # canonical-parent convention: CDS segments link to canonical mRNA only
    cds_parents = query_all(db, """SELECT DISTINCT fr.parent_feature_id
                                   FROM feature_relationships fr
                                   JOIN contigs_sequence_features csf ON fr.child_feature_id = csf.feature_id
                                   WHERE csf.feature_type='CDS' AND csf.direction='f' AND fr.relationship='part_of'""")
    _assert(len(cds_parents) == 1, f"all forward CDS segments link to the same canonical mRNA (got {len(cds_parents)} distinct parents)")

    # canonical mRNA row check: parent_feature_id must be the row where feature_group_id = feature_id
    mrna_canon = query_one(db, """SELECT feature_id FROM contigs_sequence_features
                                  WHERE feature_type='mRNA' AND direction='f' AND segment_order=0""")[0]
    _assert(cds_parents[0][0] == mrna_canon, f"forward CDS parent is the canonical mRNA (got {cds_parents[0][0]} vs canonical {mrna_canon})")

    # every CDS segment generates its own row in feature_relationships
    cds_rel_count = query_one(db, """SELECT COUNT(*) FROM feature_relationships fr
                                     JOIN contigs_sequence_features csf ON fr.child_feature_id = csf.feature_id
                                     WHERE csf.feature_type='CDS' AND csf.direction='f'""")[0]
    _assert(cds_rel_count == 3, f"3 forward-CDS relationship rows (one per segment) (got {cds_rel_count})")


def probe_origin_crossing(work):
    db = os.path.join(work, 'CIRC.db')
    build_contigs_db(os.path.join(HERE, 'circular.fa'), db)
    import_into(db, os.path.join(HERE, 'circular.gb'), 'circ_probe')

    rows = query_all(db, """SELECT start, stop, segment_order, direction, feature_group_id
                            FROM contigs_sequence_features
                            ORDER BY segment_order""")
    _assert(len(rows) == 2, f"origin-crossing feature stored as exactly 2 segments (got {len(rows)})")
    seg0, seg1 = rows
    _assert(seg0[0] == 899 and seg0[1] == 1000 and seg0[2] == 0, f"segment_order=0 is [899,1000) (got {seg0})")
    _assert(seg1[0] == 0  and seg1[1] == 200  and seg1[2] == 1, f"segment_order=1 is [0,200) (got {seg1})")
    _assert(seg0[4] == seg1[4], f"both segments share feature_group_id (got {seg0[4]} vs {seg1[4]})")

    # external_id replicates across both segments (locus_tag)
    distinct = query_all(db, "SELECT DISTINCT external_id FROM contigs_sequence_features")
    _assert(len(distinct) == 1 and distinct[0][0] == 'ORIGIN_CROSSER_001', f"both segments share external_id (got {distinct})")


def probe_malformed_linear(work):
    db = os.path.join(work, 'LIN.db')
    build_contigs_db(os.path.join(HERE, 'linear_malformed.fa'), db)
    import_into(db, os.path.join(HERE, 'linear_malformed.gb'), 'lin_probe')

    # the malformed feature should be absent; the healthy one should be present
    rows = query_all(db, "SELECT external_id, start, stop FROM contigs_sequence_features")
    ids = {r[0] for r in rows}
    _assert('LINEAR_BAD_001' not in ids, f"malformed origin-crossing on linear contig was skipped (rows: {rows})")
    _assert('LINEAR_OK_001'  in ids,    f"healthy feature on the same contig was still imported (rows: {rows})")


def probe_unknown_type(work):
    db = os.path.join(work, 'UNK.db')
    build_contigs_db(os.path.join(HERE, 'unknown_type.fa'), db)
    import_into(db, os.path.join(HERE, 'unknown_type.gb'), 'unk_probe')

    reg = query_one(db, "SELECT is_builtin, has_dedicated_table FROM feature_types WHERE feature_type='regulatory'")
    _assert(reg == (0, 0), f"non-builtin 'regulatory' registered with is_builtin=0, has_dedicated_table=0 (got {reg})")
    n_reg = query_one(db, "SELECT COUNT(*) FROM contigs_sequence_features WHERE feature_type='regulatory'")[0]
    _assert(n_reg == 1, f"1 regulatory feature imported (got {n_reg})")


def main():
    if len(sys.argv) != 2:
        print("usage: probe_importer.py <workdir>"); sys.exit(2)
    work = sys.argv[1]
    os.makedirs(work, exist_ok=True)

    probe_bacterial(work)
    probe_eukaryotic(work)
    probe_origin_crossing(work)
    probe_malformed_linear(work)
    probe_unknown_type(work)
    print("All GenbankFeatureImporter parsing probes PASSED.")


if __name__ == '__main__':
    main()
