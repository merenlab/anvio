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

    # 8 literal (4 genes + 4 CDS) + 8 synthesized (4 case-2 transcripts + 4 case-2 exons) = 16
    total = query_one(db, "SELECT COUNT(*) FROM contigs_sequence_features")[0]
    _assert(total == 16, f"bacterial import inserted 16 features (8 literal + 8 synthesized) (got {total})")

    # literal counts (derivation IS NULL)
    n_literal = query_one(db, "SELECT COUNT(*) FROM contigs_sequence_features WHERE derivation IS NULL")[0]
    _assert(n_literal == 8, f"8 literal rows (4 gene + 4 CDS) (got {n_literal})")

    n_gene = query_one(db, "SELECT COUNT(*) FROM contigs_sequence_features WHERE feature_type='gene'")[0]
    n_cds  = query_one(db, "SELECT COUNT(*) FROM contigs_sequence_features WHERE feature_type='CDS'")[0]
    _assert(n_gene == 4 and n_cds == 4, f"bacterial import has 4 genes + 4 CDS (got {n_gene} / {n_cds})")

    # synthesis: every bacterial gene gets one transcript + one exon, both with derivation='CDS'
    n_synth_t = query_one(db, "SELECT COUNT(*) FROM contigs_sequence_features WHERE feature_type='transcript' AND derivation='CDS'")[0]
    n_synth_e = query_one(db, "SELECT COUNT(*) FROM contigs_sequence_features WHERE feature_type='exon' AND derivation='CDS'")[0]
    _assert(n_synth_t == 4, f"4 synthesized transcripts with derivation='CDS' (got {n_synth_t})")
    _assert(n_synth_e == 4, f"4 synthesized exons with derivation='CDS' (got {n_synth_e})")

    # CDS → gene derives_from (unchanged by synthesis)
    derives = query_one(db, "SELECT COUNT(*) FROM feature_relationships WHERE relationship='derives_from'")[0]
    _assert(derives == 4, f"4 CDS→gene derives_from relationships (got {derives})")

    # locus_tag round-trips into external_id. Synthesized rows inherit external_id from
    # source (case 2 transcripts inherit from the gene; case 2 exons inherit from the CDS
    # segment). So all 16 rows are non-NULL.
    ext_count = query_one(db, "SELECT COUNT(*) FROM contigs_sequence_features WHERE external_id IS NOT NULL")[0]
    _assert(ext_count == 16, f"every row has external_id populated (literal + inherited by synthesis) (got {ext_count})")
    ext_literal = query_one(db, "SELECT COUNT(*) FROM contigs_sequence_features WHERE external_id IS NOT NULL AND derivation IS NULL")[0]
    _assert(ext_literal == 8, f"every literal feature has external_id populated from locus_tag (got {ext_literal})")

    # translation stored on every CDS row (all single-segment so all are canonical). Synthesis
    # doesn't add CDS_features rows.
    trans_count = query_one(db, "SELECT COUNT(*) FROM CDS_features WHERE translation IS NOT NULL")[0]
    _assert(trans_count == 4, f"4 translations stored (got {trans_count})")

    # qualifier ordering. Synthesis doesn't add qualifier rows.
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

    # canonical-parent convention: CDS segments link to the canonical literal mRNA via part_of.
    # Step 13.5b(i) also re-parents every CDS segment to the synthesized transcript with part_of,
    # so DISTINCT part_of parents is now 2. We restrict to literal-mRNA parents to verify the
    # original canonical convention still holds.
    cds_parents_literal_mrna = query_all(db, """SELECT DISTINCT fr.parent_feature_id
                                                FROM feature_relationships fr
                                                JOIN contigs_sequence_features csf  ON fr.child_feature_id  = csf.feature_id
                                                JOIN contigs_sequence_features pcsf ON fr.parent_feature_id = pcsf.feature_id
                                                WHERE csf.feature_type='CDS' AND csf.direction='f' AND fr.relationship='part_of'
                                                  AND pcsf.feature_type='mRNA' AND pcsf.derivation IS NULL""")
    _assert(len(cds_parents_literal_mrna) == 1, f"all forward CDS segments link to the same canonical literal mRNA (got {len(cds_parents_literal_mrna)} distinct mRNA parents)")
    cds_all_part_of_parents = query_all(db, """SELECT DISTINCT fr.parent_feature_id
                                               FROM feature_relationships fr
                                               JOIN contigs_sequence_features csf ON fr.child_feature_id = csf.feature_id
                                               WHERE csf.feature_type='CDS' AND csf.direction='f' AND fr.relationship='part_of'""")
    _assert(len(cds_all_part_of_parents) == 2, f"after synthesis, forward CDS part_of parents are 2 (literal mRNA + synthesized transcript) (got {len(cds_all_part_of_parents)})")

    # canonical mRNA row check: parent_feature_id must be the row where feature_group_id = feature_id
    mrna_canon = query_one(db, """SELECT feature_id FROM contigs_sequence_features
                                  WHERE feature_type='mRNA' AND direction='f' AND segment_order=0""")[0]
    _assert(cds_parents_literal_mrna[0][0] == mrna_canon, f"forward CDS literal-mRNA parent is the canonical mRNA (got {cds_parents_literal_mrna[0][0]} vs canonical {mrna_canon})")

    # every CDS segment generates its own row in feature_relationships pointing to the
    # canonical literal mRNA. Synthesis additionally re-parents every CDS segment to the
    # synthesized transcript with `part_of` (step 13.5b(i)) — so total relationships are 6.
    cds_rel_to_literal_mrna = query_one(db, """SELECT COUNT(*) FROM feature_relationships fr
                                               JOIN contigs_sequence_features csf  ON fr.child_feature_id  = csf.feature_id
                                               JOIN contigs_sequence_features pcsf ON fr.parent_feature_id = pcsf.feature_id
                                               WHERE csf.feature_type='CDS' AND csf.direction='f'
                                                 AND pcsf.feature_type='mRNA' AND pcsf.derivation IS NULL""")[0]
    _assert(cds_rel_to_literal_mrna == 3, f"3 forward-CDS → literal-mRNA relationships (one per segment) (got {cds_rel_to_literal_mrna})")
    cds_rel_total = query_one(db, """SELECT COUNT(*) FROM feature_relationships fr
                                     JOIN contigs_sequence_features csf ON fr.child_feature_id = csf.feature_id
                                     WHERE csf.feature_type='CDS' AND csf.direction='f'""")[0]
    _assert(cds_rel_total == 6, f"6 total forward-CDS relationship rows (3 to literal mRNA + 3 to synthesized transcript) (got {cds_rel_total})")


def probe_origin_crossing(work):
    db = os.path.join(work, 'CIRC.db')
    build_contigs_db(os.path.join(HERE, 'circular.fa'), db)
    import_into(db, os.path.join(HERE, 'circular.gb'), 'circ_probe')

    # the literal gene is stored as exactly 2 origin-crossing segments
    rows = query_all(db, """SELECT start, stop, segment_order, direction, feature_group_id
                            FROM contigs_sequence_features
                            WHERE feature_type='gene' AND derivation IS NULL
                            ORDER BY segment_order""")
    _assert(len(rows) == 2, f"origin-crossing literal gene stored as exactly 2 segments (got {len(rows)})")
    seg0, seg1 = rows
    _assert(seg0[0] == 899 and seg0[1] == 1000 and seg0[2] == 0, f"literal segment_order=0 is [899,1000) (got {seg0})")
    _assert(seg1[0] == 0  and seg1[1] == 200  and seg1[2] == 1, f"literal segment_order=1 is [0,200) (got {seg1})")
    _assert(seg0[4] == seg1[4], f"both literal segments share feature_group_id (got {seg0[4]} vs {seg1[4]})")

    # synthesis (case 3 — origin-crossing lone gene): the transcript mirrors the source's
    # 2-segment structure (per the origin-crossing rule), and exons are one per segment.
    synth_t = query_all(db, """SELECT start, stop, segment_order FROM contigs_sequence_features
                               WHERE feature_type='transcript' AND derivation='gene'
                               ORDER BY segment_order""")
    _assert(len(synth_t) == 2, f"synthesized transcript mirrors origin-crossing source (2 segments) (got {synth_t})")
    _assert(synth_t[0][0] == 899 and synth_t[0][1] == 1000 and synth_t[0][2] == 0, f"synth transcript seg0 = [899,1000) (got {synth_t[0]})")
    _assert(synth_t[1][0] == 0   and synth_t[1][1] == 200  and synth_t[1][2] == 1, f"synth transcript seg1 = [0,200) (got {synth_t[1]})")

    synth_e = query_all(db, """SELECT start, stop FROM contigs_sequence_features
                               WHERE feature_type='exon' AND derivation='gene'
                               ORDER BY start""")
    _assert(len(synth_e) == 2, f"two synthesized exons (one per gene segment) (got {synth_e})")

    # external_id replicates across literal segments and is inherited by every synthesized row
    distinct = query_all(db, "SELECT DISTINCT external_id FROM contigs_sequence_features")
    _assert(len(distinct) == 1 and distinct[0][0] == 'ORIGIN_CROSSER_001', f"every row shares external_id (got {distinct})")


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
