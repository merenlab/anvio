#!/bin/bash
source 00.sh
set -e

# Standalone integration test for read-edge clip profiling. Builds a fully
# synthetic reference + several read sources in-script so the test exercises
# every branch of the clip detector (JUNCTION, UNMAPPED with sequence
# from own SEQ, UNMAPPED with primary-fetched sequence, reverse-strand
# sibling) without depending on external test data.
#
# A slimmed-down version of this test (a single chimera_tail BAM with assertions
# on JUNCTION + UNMAPPED rows) is also part of run_component_tests_for_metagenomics.sh
# under the "Testing read-edge clip profiling" section, so the clip code path
# is exercised whenever the metagenomics suite runs.
#
# Reference design:
#   ref.fa = two unrelated 20 kb contigs (contig_X, contig_Y).
#
# Read sources:
#   chimera_source.fa         contig_X[5000:15000] + contig_Y[5000:15000]
#                             2-way split alignment; the read is FULLY covered
#                             so every clip event is JUNCTION.
#
#   extended_x_source.fa      contig_X + 2 kb random tail
#                             Tail aligns nowhere → primary's right SOFT clip
#                             has no SA → state='UNMAPPED', sequence carries
#                             the tail bases (read directly from primary's SEQ).
#
#   chimera_rc_source.fa      contig_X[5000:15000] + reverse_complement(contig_Y[5000:15000])
#                             Same 2-way split, but the contig_Y portion maps to
#                             the - strand. Exercises strand-flip math in the
#                             alignment-map computation.
#
#   chimera_tail_source.fa    contig_X[5000:15000] + contig_Y[5000:15000] + 2 kb random tail
#                             Reads spanning all three regions produce a HARD-
#                             right clip on the contig_Y supplementary where
#                             the bases are not in this record's SEQ; the
#                             detector fetches the primary and extracts them.
#
# Scenarios:
#   A  chimera                              → all events JUNCTION, sequence=''
#   C  chimera + --skip-clip-profiling      → clippings table is empty
#   D  extended_x                           → SOFT UNMAPPED with sequence
#   E  chimera_rc                           → reverse-strand sibling, JUNCTION
#   F  chimera_tail                         → HARD UNMAPPED via primary fetch
#   M  anvi-merge A + D                     → merged table has both samples

INFO "Creating output directory"
cd sandbox
rm -rf test-output-clippings
mkdir test-output-clippings
cd test-output-clippings

INFO "Anvi'o version"
anvi-profile --version

INFO "Building synthetic reference and read sources"
python - <<'PYEOF'
import random
random.seed(42)


def rand_seq(n):
    return ''.join(random.choices('ACGT', k=n))


def rev_comp(s):
    return s.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


contig_x = rand_seq(20000)
contig_y = rand_seq(20000)
tail_for_extended_x = rand_seq(2000)
tail_for_chimera = rand_seq(2000)

with open('ref.fa', 'w') as f:
    f.write(f'>contig_X\n{contig_x}\n')
    f.write(f'>contig_Y\n{contig_y}\n')

chimera = contig_x[5000:15000] + contig_y[5000:15000]
with open('chimera_source.fa', 'w') as f:
    f.write(f'>chimera\n{chimera}\n')

# Append a 2 kb random tail to contig_X. Reads spanning the boundary will have
# a soft clip with no SA tag — the canonical 'truly unmapped' case.
extended_x = contig_x + tail_for_extended_x
with open('extended_x_source.fa', 'w') as f:
    f.write(f'>extended_x\n{extended_x}\n')

# Plain contig_X used as a short-read source for the bowtie2 scenarios.
with open('pure_x_source.fa', 'w') as f:
    f.write(f'>contig_X\n{contig_x}\n')

# Chimera with contig_Y portion reverse-complemented. minimap2 will place the
# supplementary on contig_Y on the '-' strand, so the SA tag carries strand '-'
# and our detector has to flip the SA's read coordinates.
chimera_rc = contig_x[5000:15000] + rev_comp(contig_y[5000:15000])
with open('chimera_rc_source.fa', 'w') as f:
    f.write(f'>chimera_rc\n{chimera_rc}\n')

# Chimera with a 2 kb tail appended. Reads spanning all three regions create
# a supplementary with right HARD clip whose 'outside' is the unmapped tail.
# The bases for that HARD UNMAPPED row are recovered by fetching the
# primary record from the BAM by query_name.
chimera_tail = contig_x[5000:15000] + contig_y[5000:15000] + tail_for_chimera
with open('chimera_tail_source.fa', 'w') as f:
    f.write(f'>chimera_tail\n{chimera_tail}\n')

# Remember the tails so the assertion script can match the stored bases
# against the truth.
with open('TAILS.txt', 'w') as f:
    f.write(f'extended_x_tail\t{tail_for_extended_x}\n')
    f.write(f'chimera_tail\t{tail_for_chimera}\n')
PYEOF

INFO "Generating HiFi reads from each source"
anvi-script-gen-reads -f chimera_source.fa       --preset pacbio-hifi --coverage 50 --seed 42 -o chimera_reads
anvi-script-gen-reads -f extended_x_source.fa    --preset pacbio-hifi --coverage 50 --seed 17 -o extended_x_reads
anvi-script-gen-reads -f chimera_rc_source.fa    --preset pacbio-hifi --coverage 50 --seed 31 -o chimera_rc_reads
anvi-script-gen-reads -f chimera_tail_source.fa  --preset pacbio-hifi --coverage 50 --seed 53 -o chimera_tail_reads

INFO "Mapping each read set with minimap2 (-Y off so supplementaries hard-clip)"
for sample in chimera extended_x chimera_rc chimera_tail; do
    minimap2 -ax map-hifi --MD ref.fa ${sample}_reads.fastq 2>/dev/null | samtools sort -o ${sample}.bam -
    samtools index ${sample}.bam
done

INFO "Sanity check: chimera.bam carries both S and H clip operations"
soft_count=$(samtools view chimera.bam | awk '{print $6}' | grep -c 'S' || true)
hard_count=$(samtools view chimera.bam | awk '{print $6}' | grep -c 'H' || true)
echo "Reads with S in CIGAR: $soft_count"
echo "Reads with H in CIGAR: $hard_count"
if [ "$soft_count" -eq 0 ] || [ "$hard_count" -eq 0 ]; then
    echo "ERROR: expected both soft- and hard-clipped alignments in chimera.bam"
    exit 1
fi

INFO "Generating anvi'o contigs database"
anvi-gen-contigs-database -f ref.fa -o ref.db -L -1 --skip-gene-calling

INFO "Scenario A: anvi-profile on chimera (default clip profiling)"
anvi-profile -i chimera.bam -c ref.db -o profile_A -M 0 \
             --skip-hierarchical-clustering --sample-name chimera_A

INFO "Scenario C: anvi-profile on chimera with --skip-clip-profiling"
anvi-profile -i chimera.bam -c ref.db -o profile_C -M 0 \
             --skip-hierarchical-clustering --skip-clip-profiling \
             --sample-name chimera_C

INFO "Scenario D: anvi-profile on extended_x (truly unmapped tail)"
anvi-profile -i extended_x.bam -c ref.db -o profile_D -M 0 \
             --skip-hierarchical-clustering --sample-name extended_x_D

INFO "Scenario E: anvi-profile on chimera_rc (reverse-strand sibling)"
anvi-profile -i chimera_rc.bam -c ref.db -o profile_E -M 0 \
             --skip-hierarchical-clustering --sample-name chimera_rc_E

INFO "Scenario F: anvi-profile on chimera_tail (HARD UNMAPPED → primary fetch)"
anvi-profile -i chimera_tail.bam -c ref.db -o profile_F -M 0 \
             --skip-hierarchical-clustering --sample-name chimera_tail_F

INFO "Asserting per-scenario clippings table contents"
python - <<'PYEOF'
import sys
import anvio.tables as t
from anvio.dbinfo import ProfileDBInfo


def load_clips(profile_path):
    db = ProfileDBInfo(profile_path).load_db()
    rows = list(db.get_table_as_dict(t.clippings_table_name).values())
    db.disconnect()
    return rows


def by_state(rows):
    out = {'JUNCTION': 0, 'JUNCTION_WITH_GAP': 0, 'UNMAPPED': 0}
    for r in rows:
        out[r['state']] = out.get(r['state'], 0) + 1
    return out


def show_first(rows, label, n=3):
    print(f"  {label}: {len(rows)} rows")
    for r in rows[:n]:
        seq_preview = r['sequence'][:30] + ('…' if len(r['sequence']) > 30 else '')
        print(f"    - split={r['split_name']} pos={r['pos']} type={r['type']} side={r['side']} "
              f"state={r['state']} len(seq)={len(r['sequence'])} length={r['length']} "
              f"count={r['count']} partner={r['partner_contig']}:{r['partner_junction_pos']}({r['partner_strand']}) "
              f"seq[:30]='{seq_preview}'")


# Pull the canonical tails so we can verify recovered bases.
with open('TAILS.txt') as f:
    tails = dict(line.strip().split('\t') for line in f)

fail = []

# ---------------- Scenario A: chimera (fully covered → all JUNCTION) ----------------
print("\n[Scenario A] chimera reads → all JUNCTION, sequence=''")
clips_A = load_clips('profile_A/PROFILE.db')
print(f"  total clip events: {len(clips_A)}")
print(f"  state distribution: {by_state(clips_A)}")
if len(clips_A) == 0:
    fail.append("A: clippings table is empty")
if any(r['state'] != 'JUNCTION' for r in clips_A):
    bad = [r for r in clips_A if r['state'] != 'JUNCTION'][:3]
    fail.append(f"A: expected every row to be JUNCTION; found {len(clips_A) - sum(1 for r in clips_A if r['state']=='JUNCTION')} non-JUNCTION. e.g. {bad}")
if any(r['sequence'] != '' for r in clips_A):
    fail.append("A: JUNCTION rows should have sequence='' (the bases live in the partner alignment)")
soft_X_right = [r for r in clips_A if r['type'] == 'SOFT' and r['side'] == 'R' and r['split_name'].startswith('contig_X')]
hard_Y_left  = [r for r in clips_A if r['type'] == 'HARD' and r['side'] == 'L' and r['split_name'].startswith('contig_Y')]
show_first(soft_X_right, "SOFT right-clips on contig_X")
show_first(hard_Y_left,  "HARD left-clips on contig_Y")
if not soft_X_right:
    fail.append("A: no SOFT right-clips on contig_X at the chimera junction")
if not hard_Y_left:
    fail.append("A: no HARD left-clips on contig_Y at the chimera junction")
if soft_X_right and {r['partner_contig'] for r in soft_X_right} != {'contig_Y'}:
    fail.append(f"A: SOFT right partners not all 'contig_Y': {set(r['partner_contig'] for r in soft_X_right)}")
if hard_Y_left and {r['partner_contig'] for r in hard_Y_left} != {'contig_X'}:
    fail.append(f"A: HARD left partners not all 'contig_X': {set(r['partner_contig'] for r in hard_Y_left)}")

# ---------------- Scenario C: --skip-clip-profiling → empty table ----------------
print("\n[Scenario C] --skip-clip-profiling → empty table")
clips_C = load_clips('profile_C/PROFILE.db')
print(f"  total clip events: {len(clips_C)}")
if len(clips_C) != 0:
    fail.append(f"C: clippings table is not empty ({len(clips_C)} rows)")

# ---------------- Scenario D: SOFT UNMAPPED (extended_x, no SA) ----------------
print("\n[Scenario D] extended_x reads → SOFT UNMAPPED with sequence, no partner")
clips_D = load_clips('profile_D/PROFILE.db')
print(f"  total clip events: {len(clips_D)}")
print(f"  state distribution: {by_state(clips_D)}")
unexp_D = [r for r in clips_D if r['state'] == 'UNMAPPED']
exp_D = [r for r in clips_D if r['state'] == 'JUNCTION']
show_first(unexp_D, "UNMAPPED rows")

if not unexp_D:
    fail.append("D: expected UNMAPPED rows (right soft clip at contig_X end has no SA)")
else:
    # All UNMAPPED rows here should be SOFT, side R, on contig_X, with no partner.
    for r in unexp_D:
        if r['type'] != 'SOFT':
            fail.append(f"D: UNMAPPED row has type={r['type']}; expected SOFT")
            break
        if r['side'] != 'R':
            fail.append(f"D: UNMAPPED row has side={r['side']}; expected R")
            break
        if not r['split_name'].startswith('contig_X'):
            fail.append(f"D: UNMAPPED row on unexpected split {r['split_name']}")
            break
        if r['partner_contig'] != '' or r['partner_junction_pos'] != -1:
            fail.append(f"D: UNMAPPED row has unexpected partner: {r['partner_contig']}:{r['partner_junction_pos']}")
            break
        if len(r['sequence']) == 0:
            fail.append("D: UNMAPPED row has empty sequence; expected the unmapped tail bases")
            break

    # The bases of at least one UNMAPPED row should appear at the start of the
    # appended random tail (some reads will end mid-tail so the stored sequence
    # is a prefix of the appended tail).
    true_tail = tails['extended_x_tail']
    matches = sum(1 for r in unexp_D if r['sequence'] and true_tail.startswith(r['sequence']))
    print(f"  UNMAPPED rows whose sequence is a prefix of the truth tail: {matches}/{len(unexp_D)}")
    if matches == 0:
        fail.append("D: no UNMAPPED row's sequence matches the truth-tail prefix")

# ---------------- Scenario E: reverse-strand sibling ----------------
print("\n[Scenario E] chimera_rc reads → reverse-strand sibling, JUNCTION")
clips_E = load_clips('profile_E/PROFILE.db')
print(f"  total clip events: {len(clips_E)}")
print(f"  state distribution: {by_state(clips_E)}")
soft_X_right_E = [r for r in clips_E if r['type'] == 'SOFT' and r['side'] == 'R' and r['split_name'].startswith('contig_X')]
show_first(soft_X_right_E, "SOFT right-clips on contig_X (RC chimera)")

if not soft_X_right_E:
    fail.append("E: no SOFT right-clips on contig_X for the RC chimera")
elif any(r['sequence'] != '' for r in soft_X_right_E):
    fail.append("E: JUNCTION rows should have sequence=''")
else:
    # The chimera_rc reads are fully covered (just with the supp on the - strand on contig_Y),
    # so every row must be JUNCTION. UNMAPPED would mean the strand-flip math is wrong —
    # the SA-derived interval would land off-side from the primary's clip and the merged
    # alignment map would have a gap covering it, flipping the row's state.
    unexp_E = sum(1 for r in soft_X_right_E if r['state'] != 'JUNCTION')
    if unexp_E > 0:
        fail.append(f"E: {unexp_E}/{len(soft_X_right_E)} SOFT right-clips on contig_X are not JUNCTION — strand-flip math is likely wrong")

    # We expect *some* partner_strand='r' rows: those reads where minimap2 placed the primary
    # on contig_X (+) with the supplementary on contig_Y (-). The mix of 'r' and 'f' is normal —
    # minimap2 picks the higher-scoring half as primary, so reads can have either contig_X or
    # contig_Y as primary. What matters is that 'r' shows up at all, since only the strand-flip
    # math produces an JUNCTION row on contig_X when the partner is on the opposite strand.
    strand_counts = {}
    for r in soft_X_right_E:
        strand_counts[r['partner_strand']] = strand_counts.get(r['partner_strand'], 0) + 1
    print(f"  partner_strand distribution among SOFT-right on contig_X: {strand_counts}")
    if strand_counts.get('r', 0) == 0:
        fail.append("E: expected some SOFT right-clips on contig_X with partner_strand='r' (reverse-strand sibling)")

# ---------------- Scenario F: HARD UNMAPPED → primary fetch ----------------
print("\n[Scenario F] chimera_tail reads → HARD UNMAPPED with bases from primary fetch")
clips_F = load_clips('profile_F/PROFILE.db')
print(f"  total clip events: {len(clips_F)}")
print(f"  state distribution: {by_state(clips_F)}")
hard_unexp_F = [r for r in clips_F if r['type'] == 'HARD' and r['state'] == 'UNMAPPED']
soft_unexp_F = [r for r in clips_F if r['type'] == 'SOFT' and r['state'] == 'UNMAPPED']
show_first(hard_unexp_F, "HARD UNMAPPED rows")
show_first(soft_unexp_F, "SOFT UNMAPPED rows")

# Either type can carry the tail signal depending on which contig minimap2 picks
# as primary, but at least one must exist.
if not hard_unexp_F and not soft_unexp_F:
    fail.append("F: no UNMAPPED rows; expected the tail to surface as either SOFT or HARD UNMAPPED")

# The unmapped sequence must be a prefix (or full) of the truth tail. This is
# the critical check that the primary-fetch + (potentially) strand-conversion
# code recovers the right bases.
true_tail = tails['chimera_tail']
all_unexp_F = hard_unexp_F + soft_unexp_F
matches = sum(1 for r in all_unexp_F if r['sequence'] and (true_tail.startswith(r['sequence']) or r['sequence'].startswith(true_tail[:200])))
print(f"  UNMAPPED rows whose sequence matches the tail (prefix): {matches}/{len(all_unexp_F)}")
if all_unexp_F and matches == 0:
    print("  (Showing first few sequences vs truth tail prefix for diagnosis)")
    print(f"  truth tail [:60]: {true_tail[:60]}")
    for r in all_unexp_F[:3]:
        print(f"  observed [:60]:   {r['sequence'][:60]} (type={r['type']})")
    fail.append("F: no UNMAPPED row's sequence matches the truth tail; primary fetch may not be recovering bases correctly")

# Regression guard for the per-clip partner attribution fix: UNMAPPED rows must
# NOT carry partner info — those bases don't relate to any SA sibling, so the
# partner_* columns should be empty. Before this fix the code blindly used
# sa_entries[0] as the partner for every clip event regardless of whether that
# entry's M region actually covered the clip's outside.
unmapped_with_partner = [r for r in all_unexp_F if r['partner_contig'] != '' or r['partner_junction_pos'] != -1]
if unmapped_with_partner:
    bad = unmapped_with_partner[:3]
    fail.append(f"F: {len(unmapped_with_partner)} UNMAPPED row(s) have a partner set; expected empty for UNMAPPED. e.g. partner={bad[0]['partner_contig']}:{bad[0]['partner_junction_pos']}")

# Wrap-up
if fail:
    print("\nFAILED:")
    for msg in fail:
        print(f"  - {msg}")
    sys.exit(1)

print("\nAll single-profile assertions passed")
PYEOF

INFO "Merging profiles A and D (mixes JUNCTION + UNMAPPED states)"
anvi-merge profile_A/PROFILE.db profile_D/PROFILE.db -c ref.db -o profile_MERGED

INFO "Asserting merged clippings table"
python - <<'PYEOF'
import sys
import anvio.tables as t
from anvio.dbinfo import ProfileDBInfo

db = ProfileDBInfo('profile_MERGED/PROFILE.db').load_db()
rows = list(db.get_table_as_dict(t.clippings_table_name).values())
db.disconnect()

samples = {r['sample_id'] for r in rows}
states_per_sample = {}
for r in rows:
    states_per_sample.setdefault(r['sample_id'], set()).add(r['state'])

print(f"  merged clip events: {len(rows)}")
print(f"  samples present:    {sorted(samples)}")
print(f"  states per sample:  { {k: sorted(v) for k, v in states_per_sample.items()} }")

fail = []
if 'chimera_A' not in samples:
    fail.append("merge: sample 'chimera_A' missing")
if 'extended_x_D' not in samples:
    fail.append("merge: sample 'extended_x_D' missing")
# After merge we should still see UNMAPPED rows from D and JUNCTION rows from A.
if 'extended_x_D' in states_per_sample and 'UNMAPPED' not in states_per_sample['extended_x_D']:
    fail.append("merge: extended_x_D's UNMAPPED rows did not survive the merge")
if 'chimera_A' in states_per_sample and 'JUNCTION' not in states_per_sample['chimera_A']:
    fail.append("merge: chimera_A's JUNCTION rows did not survive the merge")

if fail:
    print("FAILED:")
    for msg in fail:
        print(f"  - {msg}")
    sys.exit(1)

print("Merge assertions passed")
PYEOF

#
# Scenarios G + H verify the BAM @PG defensive check. Build short reads from
# contig_X and map them with bowtie2 in default end-to-end mode (which does NOT
# emit soft clips). The @PG record will say `bowtie2`, anvi-profile should
# detect this and auto-skip clip profiling with a CLIP PROFILING AUTO-SKIPPED
# warning. Scenario H reruns the same BAM with --force-clip-profiling and
# expects the warning NOT to fire (and clippings table remains empty since
# there genuinely are no clips in the BAM).
#

INFO "Building Illumina short reads from contig_X for the bowtie2 scenarios"
anvi-script-gen-reads -f pure_x_source.fa --preset illumina-paired --coverage 30 --seed 71 -o bt2_reads

INFO "Mapping with bowtie2 default (end-to-end mode — should not emit clips)"
bowtie2-build ref.fa bt2_ref >/dev/null 2>&1
bowtie2 -x bt2_ref -1 bt2_reads-R1.fastq -2 bt2_reads-R2.fastq -p 1 2>/dev/null | samtools sort -o bt2.bam -
samtools index bt2.bam

INFO "Sanity check: bowtie2 default BAM should NOT carry S/H clips"
clip_count=$(samtools view bt2.bam | awk '{print $6}' | grep -c -E '[SH]' || true)
echo "Reads with S or H in CIGAR: $clip_count"
if [ "$clip_count" -gt 0 ]; then
    echo "ERROR: bowtie2 default-mode BAM unexpectedly carries clip CIGAR ops"
    exit 1
fi

INFO "Scenario G: anvi-profile on bowtie2 BAM (expect AUTO-SKIPPED warning)"
anvi-profile -i bt2.bam -c ref.db -o profile_G -M 0 \
             --skip-hierarchical-clustering --sample-name bt2_G 2>&1 | tee profile_G_run.log
if ! grep -q "CLIP PROFILING AUTO-SKIPPED" profile_G_run.log; then
    echo "ERROR: Scenario G did not emit the CLIP PROFILING AUTO-SKIPPED warning for a bowtie2 default BAM"
    exit 1
fi

INFO "Scenario H: anvi-profile on bowtie2 BAM with --force-clip-profiling (override)"
anvi-profile -i bt2.bam -c ref.db -o profile_H -M 0 \
             --skip-hierarchical-clustering --force-clip-profiling \
             --sample-name bt2_H 2>&1 | tee profile_H_run.log
if grep -q "CLIP PROFILING AUTO-SKIPPED" profile_H_run.log; then
    echo "ERROR: Scenario H with --force-clip-profiling unexpectedly emitted AUTO-SKIPPED"
    exit 1
fi

INFO "Asserting bowtie2 scenarios"
python - <<'PYEOF'
import sys
import anvio.tables as t
from anvio.dbinfo import ProfileDBInfo

fail = []

def load_clips_and_flag(profile_path):
    info = ProfileDBInfo(profile_path)
    db = info.load_db()
    rows = list(db.get_table_as_dict(t.clippings_table_name).values())
    meta = db.get_table_as_dict('self')
    db.disconnect()
    clips_profiled = meta.get('clips_profiled', {}).get('value')
    return rows, clips_profiled

# Scenario G: auto-skipped → clippings table empty, clips_profiled='False'
rows_G, flag_G = load_clips_and_flag('profile_G/PROFILE.db')
print(f"\n[Scenario G] bowtie2 default → expect auto-skip")
print(f"  clippings rows: {len(rows_G)}")
print(f"  self.clips_profiled: {flag_G!r}")
if len(rows_G) != 0:
    fail.append(f"G: clippings table should be empty after auto-skip; got {len(rows_G)} rows")
# clips_profiled is stored as the Python repr of the bool, ie. 'False' (string) or 0
if str(flag_G) not in ('False', '0'):
    fail.append(f"G: self.clips_profiled should record False after auto-skip; got {flag_G!r}")

# Scenario H: --force-clip-profiling → no auto-skip; clips_profiled=True; clippings table
# is still empty (the aligner just didn't emit any clip CIGAR ops).
rows_H, flag_H = load_clips_and_flag('profile_H/PROFILE.db')
print(f"\n[Scenario H] bowtie2 default + --force-clip-profiling")
print(f"  clippings rows: {len(rows_H)}")
print(f"  self.clips_profiled: {flag_H!r}")
if str(flag_H) not in ('True', '1'):
    fail.append(f"H: self.clips_profiled should record True when --force-clip-profiling is used; got {flag_H!r}")
# We don't assert rows_H == 0 strictly — bowtie2 might emit a stray clip in edge cases;
# but the BAM-CIGAR sanity check above already guards against that.

if fail:
    print("\nFAILED:")
    for msg in fail:
        print(f"  - {msg}")
    sys.exit(1)

print("\nbowtie2 scenarios passed")
PYEOF

INFO "All automated clipping component tests passed"
echo
echo "  To visually verify the inspect-page clip rendering:"
echo "    cd anvio/tests/sandbox/test-output-clippings"
echo "    anvi-inspect -c ref.db -p profile_F/PROFILE.db --split-name contig_Y_split_00001 --just-do-it"
echo "  Recommended profile DBs:"
echo "    profile_F  → both JUNCTION and UNMAPPED rows; best signal exercise"
echo "    profile_D  → only UNMAPPED rows (truly unmapped tail)"
echo "    profile_A  → only JUNCTION rows (clean 2-way split)"
