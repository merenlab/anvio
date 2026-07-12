#!/usr/bin/env python
"""Regression checks for anvi-report-dgrs VR primer construction and read recovery.

All expectations are DERIVED FROM THE INPUTS (reference.fa, mutations.tsv) and the
program's own output (DGRs_found.tsv, DGR_Primers_used_for_VR_diversity.tsv,
PRIMER_MATCHES/*.fa) -- nothing about mutation positions, frequencies, frames or
sequences is hardcoded. Editing mutations.tsv or the read-generation preset does not
require touching this script; the checks re-derive themselves from whatever is there.

For each variable region (VR) discovered by the run, it verifies:

  [1] REFERENCE MATCH  -- the assembled Whole_Primer (a regex, '.' = any base) matches
      the reference contig on at least one strand. A malformed primer (e.g. a duplicated
      base at the initial<->masked junction, a bad flank offset, or the wrong strand
      orientation) will not match the reference and is caught here.

  [2] MASK JUSTIFIED   -- every masked ('.') position in a primer's Masked_Primer is
      explained either by a TR/VR mismatch (from the alignment in DGRs_found.tsv) or by a
      real SNV in the profile.db at that position; and every unmasked position equals the
      reference base there. Catches the frame/rev-comp index off-by-one, which shifts
      masks onto conserved neighbours.

  [3] READS RECOVERED  -- each VR's PRIMER_MATCHES output contains at least one read.
      Catches primers that silently match nothing end to end.

  [4] CREATED SITES MASKED -- every multi-allele site in mutations.tsv that falls inside a
      VR and is detected as an SNV in the profile is masked in that sample's primer. This
      ties the primer back to the variation actually injected into the reads.

Usage:
  check_dgr_primers.py --dgr-output-dir DIR --reference reference.fa --contigs-db CONTIGS.db
                       [--profile-db MERGED_PROFILE/PROFILE.db] [--mutations mutations.tsv]

Exits non-zero (and prints every failure) if any check fails.
"""

import os
import re
import csv
import glob
import sqlite3
import argparse

import anvio.utils as utils


def read_fasta(path):
    """Minimal FASTA reader -> {header: sequence} (no external deps)."""
    seqs, name, chunks = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = ''.join(chunks)
    return seqs


def count_fasta_reads(path):
    n = 0
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                n += 1
    return n


def load_dgrs(dgr_output_dir):
    # the results file is '<output-dir-basename>_DGRs_found.tsv'; exclude the
    # similarly-suffixed 'Parameters_used_in_DGRs_found.tsv'.
    basename = os.path.basename(os.path.normpath(dgr_output_dir))
    preferred = os.path.join(dgr_output_dir, f"{basename}_DGRs_found.tsv")
    if os.path.exists(preferred):
        hits = [preferred]
    else:
        hits = [p for p in glob.glob(os.path.join(dgr_output_dir, '*_DGRs_found.tsv'))
                if not os.path.basename(p).startswith('Parameters')]
    if not hits:
        raise FileNotFoundError(f"No *_DGRs_found.tsv in {dgr_output_dir}")
    vrs = []
    with open(hits[0]) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            if not r.get('VR_sequence'):
                continue
            vrs.append(r)
    return vrs


def load_primers(dgr_output_dir):
    """Return {primer_id: {'*': row, sample: row, ...}} from the compact primers TSV."""
    path = os.path.join(dgr_output_dir, 'DGR_Primers_used_for_VR_diversity.tsv')
    by_primer = {}
    with open(path) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            by_primer.setdefault(r['Primer_ID'], {})[r['Sample_ID']] = r
    return by_primer


def load_mutation_sites(mutations_path):
    """{contig: {position: n_alleles_with_nonzero_freq}} from mutations.tsv (forward, 0-based)."""
    sites = {}
    if not mutations_path or not os.path.exists(mutations_path):
        return sites
    with open(mutations_path) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            freqs = [float(r[c]) for c in ('freq_A', 'freq_T', 'freq_C', 'freq_G')]
            n_alleles = sum(1 for x in freqs if x > 0)
            sites.setdefault(r['contig_name'], {})[int(r['position'])] = n_alleles
    return sites


def snv_positions_by_sample(profile_db, contig, start, end):
    """{sample_id: set(pos_in_contig)} for SNVs in [start, end] inclusive. Empty if no db."""
    out = {}
    if not profile_db or not os.path.exists(profile_db):
        return out
    con = sqlite3.connect(profile_db)
    try:
        q = ("SELECT sample_id, pos_in_contig FROM variable_nucleotides "
             "WHERE split_name LIKE ? AND pos_in_contig >= ? AND pos_in_contig <= ?")
        for sample_id, pos in con.execute(q, (contig + '%', int(start), int(end))):
            out.setdefault(sample_id, set()).add(int(pos))
    finally:
        con.close()
    return out


def mismatch_forward_indices(tr_align, vr_align, base, frame):
    """Forward-orientation indices masked purely by the TR/VR alignment (base mask)."""
    L = len(vr_align)
    idx = set()
    for j in range(L):
        a = j if frame == 1 else (L - 1 - j)
        tr, vr = tr_align[a], vr_align[a]
        if tr == base or tr != vr:
            idx.add(j)
    return idx


def main():
    ap = argparse.ArgumentParser(description="Validate anvi-report-dgrs VR primers against the inputs.")
    ap.add_argument('--dgr-output-dir', required=True)
    ap.add_argument('--reference', required=True)
    ap.add_argument('--contigs-db', required=True)  # kept for interface symmetry / future use
    ap.add_argument('--profile-db', default=None)
    ap.add_argument('--mutations', default=None)
    args = ap.parse_args()

    reference = read_fasta(args.reference)
    vrs = load_dgrs(args.dgr_output_dir)
    primers = load_primers(args.dgr_output_dir)
    mut_sites = load_mutation_sites(args.mutations)
    matches_dir = os.path.join(args.dgr_output_dir, 'PRIMER_MATCHES')

    failures = []
    checked = 0

    for vr in vrs:
        dgr, vr_id = vr['DGR'], vr['VR']
        contig = vr['VR_contig']
        vs, ve = int(vr['VR_start_position']), int(vr['VR_end_position'])
        frame = int(vr['VR_frame_reported'])
        tr_align, vr_align = vr['TR_sequence'], vr['VR_sequence']
        base = vr['Base']
        tag = f"{dgr}_{vr_id}"

        # locate this VR's primer entry (compact TSV: rows keyed by Primer_ID)
        pid = next((p for p in primers if p.startswith(f"{dgr}_{vr_id}_Primer")), None)
        if pid is None:
            failures.append(f"{tag}: no primer row found in primers TSV")
            continue
        rows = primers[pid]

        contig_seq = reference.get(contig)
        if contig_seq is None:
            failures.append(f"{tag}: contig '{contig}' absent from reference {args.reference}")
            continue
        contig_rc = utils.rev_comp(contig_seq)
        fwd_vr = contig_seq[vs:ve + 1]

        mism = mismatch_forward_indices(tr_align, vr_align, base, frame)
        snvs_by_sample = snv_positions_by_sample(args.profile_db, contig, vs, ve)

        # ---- [1] every Whole_Primer (base + overrides) matches the reference on some strand
        for sample_id, row in rows.items():
            whole = row['Whole_Primer']
            if not (re.search(whole, contig_seq) or re.search(whole, contig_rc)):
                failures.append(f"{tag} [{sample_id}]: Whole_Primer does not match the reference "
                                f"contig on either strand (malformed junction/flank/orientation).")

        # ---- [2] masked positions justified; unmasked positions equal the reference
        for sample_id, row in rows.items():
            masked = row['Masked_Primer']
            sample_snvs = set() if sample_id == '*' else snvs_by_sample.get(sample_id, set())
            if len(masked) != len(fwd_vr):
                failures.append(f"{tag} [{sample_id}]: Masked_Primer length {len(masked)} != "
                                f"VR length {len(fwd_vr)}")
                continue
            for j, ch in enumerate(masked):
                p = vs + j
                if ch == '.':
                    if j not in mism and p not in sample_snvs:
                        failures.append(f"{tag} [{sample_id}]: spurious mask at genomic {p} "
                                        f"(forward index {j}) -- not a TR/VR mismatch and not an SNV.")
                elif ch != fwd_vr[j]:
                    failures.append(f"{tag} [{sample_id}]: conserved base at genomic {p} is "
                                    f"'{ch}' but reference is '{fwd_vr[j]}'.")

        # ---- [3] reads recovered for this VR (sum over its per-sample match files)
        vr_match_files = glob.glob(os.path.join(matches_dir, f"*-{pid}-PRIMER-MATCHES.fa"))
        total_reads = sum(count_fasta_reads(f) for f in vr_match_files)
        if not vr_match_files:
            failures.append(f"{tag}: no PRIMER_MATCHES file for primer '{pid}'.")
        elif total_reads == 0:
            failures.append(f"{tag}: primer '{pid}' recovered 0 reads across all samples "
                            f"(check the initial<->masked junction).")

        # ---- [4] created multi-allele sites that became SNVs are masked in that sample
        contig_muts = {p: n for p, n in mut_sites.get(contig, {}).items() if vs <= p <= ve and n >= 2}
        for sample_id, snvset in snvs_by_sample.items():
            row = rows.get(sample_id) or rows.get('*')
            if row is None:
                continue
            masked = row['Masked_Primer']
            if len(masked) != len(fwd_vr):
                continue
            for p in contig_muts:
                if p in snvset and masked[p - vs] != '.':
                    failures.append(f"{tag} [{sample_id}]: created SNV site at genomic {p} "
                                    f"is not masked in the primer (forward index {p - vs}).")

        checked += 1
        n_files = len(vr_match_files)
        print(f"  {tag:<22} frame={frame:>2}  primer={pid.split('_Primer_')[-1] if '_Primer_' in pid else 'NA':<2}  "
              f"reads={total_reads:<4} ({n_files} file(s))  OK")

    print(f"\nChecked {checked} VR(s).")
    if failures:
        print(f"\n{len(failures)} FAILURE(S):")
        for m in failures:
            print(f"  - {m}")
        raise SystemExit(1)
    print("All DGR primer checks passed.")


if __name__ == '__main__':
    main()
