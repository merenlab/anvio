#!/bin/bash
source 00.sh

# Component test for the lean profiling flags of anvi-profile (--min-indel-length, --skip-SNV-storage).
#
# Uses the mock-MGE data of run_component_tests_for_report_indels.sh (a transposase inserted at distinct
# sites in three HiFi variants). Every BAM is profiled twice with the same parameters: once with the new
# flags off ("full") and once with `--min-indel-length 100 --skip-SNV-storage` ("lean"). The lean profile
# must hold exactly the full profile's indels of >= 100 bp, no SNV rows, and otherwise identical data;
# anvi-report-indels must report the same thing from both; and the guards for these flags must refuse.

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2 $3
#####################################

# Skip gracefully if long-read mapping tools are missing.
if ! command -v minimap2 >/dev/null 2>&1; then
    INFO "minimap2 not installed -- skipping the lean-profiling component test"
    exit 0
fi
if ! command -v samtools >/dev/null 2>&1; then
    INFO "samtools not installed -- skipping the lean-profiling component test"
    exit 0
fi

INFO "Setting up the lean-profiling mock-MGE analysis directory"
mkdir -p $output_dir
cp -r $files/mock_data_for_indels/* $output_dir/
cd $output_dir/

INFO "Building per-sample variant FASTAs (shared backbone + MGE at distinct sites)"
python3 build_variants.py

INFO "Generating PacBio-HiFi reads for each variant (20 kb mean, ~Q40)"
mkdir -p 01_FASTQ
for sample in S01 S02 S03
do
    anvi-script-gen-reads -f 01_FASTA/${sample}.fa \
                          --preset pacbio-hifi \
                          --read-length 20000 \
                          --error-rate 0.0001 \
                          --coverage 50 \
                          --seed 42 \
                          -o 01_FASTQ/${sample}
done

INFO "Building reference contigs-db from S01"
mkdir -p 02_CONTIGS
anvi-gen-contigs-database -f 01_FASTA/S01.fa \
                          -o 02_CONTIGS/CONTIGS.db \
                          -n MGE_TEST \
                          --quiet \
                          $thread_controller

INFO "Mapping HiFi reads of every variant to the S01 reference with minimap2"
mkdir -p 03_MAPPING
for sample in S01 S02 S03
do
    minimap2 -ax map-hifi 01_FASTA/S01.fa 01_FASTQ/${sample}.fastq 2>/dev/null \
        | samtools sort -o 03_MAPPING/${sample}.bam -
    samtools index 03_MAPPING/${sample}.bam
done

INFO "Profiling each BAM twice: new flags off (full) and --min-indel-length 100 --skip-SNV-storage (lean)"
for sample in S01 S02 S03
do
    anvi-profile -c 02_CONTIGS/CONTIGS.db \
                 -i 03_MAPPING/${sample}.bam \
                 -o 04_PROFILE_FULL/${sample} \
                 --sample-name ${sample} \
                 --min-contig-length 1000 \
                 --report-variability-full \
                 --quiet \
                 $thread_controller

    anvi-profile -c 02_CONTIGS/CONTIGS.db \
                 -i 03_MAPPING/${sample}.bam \
                 -o 04_PROFILE_LEAN/${sample} \
                 --sample-name ${sample} \
                 --min-contig-length 1000 \
                 --report-variability-full \
                 --min-indel-length 100 \
                 --skip-SNV-storage \
                 --quiet \
                 $thread_controller
done

INFO "Comparing lean and full profiles table by table"
python3 - <<'PYEOF'
import collections, gzip, sqlite3, sys

def rows(db, sql):
    con = sqlite3.connect(db)
    try:
        return con.execute(sql).fetchall()
    finally:
        con.close()

def fail(msg):
    print(f"FAIL: {msg}")
    sys.exit(1)

n_short_full = 0
for s in ['S01', 'S02', 'S03']:
    full, lean = f'04_PROFILE_FULL/{s}/PROFILE.db', f'04_PROFILE_LEAN/{s}/PROFILE.db'
    mf, ml = dict(rows(full, 'SELECT key, value FROM self')), dict(rows(lean, 'SELECT key, value FROM self'))
    if (mf.get('SNVs_stored'), mf.get('min_indel_length')) != ('1', '0'):
        fail(f"{s}: full profile meta SNVs_stored/min_indel_length = {mf.get('SNVs_stored')}/{mf.get('min_indel_length')}, expected 1/0")
    if (ml.get('SNVs_stored'), ml.get('min_indel_length')) != ('0', '100'):
        fail(f"{s}: lean profile meta SNVs_stored/min_indel_length = {ml.get('SNVs_stored')}/{ml.get('min_indel_length')}, expected 0/100")

    snv_full = rows(full, 'SELECT COUNT(*) FROM variable_nucleotides')[0][0]
    snv_lean = rows(lean, 'SELECT COUNT(*) FROM variable_nucleotides')[0][0]
    if snv_full == 0:
        fail(f"{s}: the full profile has no SNV rows, so the test would not show anything")
    if snv_lean != 0:
        fail(f"{s}: the lean profile holds {snv_lean} SNV rows (expected none)")

    ind_full = rows(full, 'SELECT * FROM indels')
    ind_lean = rows(lean, 'SELECT * FROM indels')
    length_ix = [r[1] for r in sqlite3.connect(full).execute("PRAGMA table_info(indels)")].index('length')
    n_short_full += sum(1 for r in ind_full if r[length_ix] < 100)
    if collections.Counter(r for r in ind_full if r[length_ix] >= 100) != collections.Counter(ind_lean):
        fail(f"{s}: lean indels are not exactly the full profile's indels >= 100 bp")

    for table in ['variability_splits', 'mean_coverage_splits', 'detection_splits']:
        if collections.Counter(rows(full, f'SELECT * FROM {table}')) != collections.Counter(rows(lean, f'SELECT * FROM {table}')):
            fail(f"{s}: table {table} differs between lean and full profiles")

    un = lambda c: gzip.decompress(c) if isinstance(c, bytes) else c
    cov = lambda p: collections.Counter((a, b, un(c)) for a, b, c in rows(p, 'SELECT split_name, sample_name, coverages FROM split_coverages'))
    if cov(f'04_PROFILE_FULL/{s}/AUXILIARY-DATA.db') != cov(f'04_PROFILE_LEAN/{s}/AUXILIARY-DATA.db'):
        fail(f"{s}: per-nucleotide coverages differ between lean and full profiles")
    print(f"{s}: lean == full (>= 100 bp): {len(ind_lean)} indel rows, {snv_full} SNV rows dropped")

if n_short_full == 0:
    fail("no indel shorter than 100 bp in any full profile, so --min-indel-length was not exercised")
print(f"ok: {n_short_full} indel rows < 100 bp were dropped by --min-indel-length")
PYEOF

INFO "Running anvi-report-indels on the full and on the lean profiles"
for kind in FULL LEAN
do
    anvi-report-indels -c 02_CONTIGS/CONTIGS.db \
                       -p 04_PROFILE_${kind}/S01/PROFILE.db \
                       -p 04_PROFILE_${kind}/S02/PROFILE.db \
                       -p 04_PROFILE_${kind}/S03/PROFILE.db \
                       --min-samples 1 \
                       --min-length 400 \
                       -o INDELS_${kind} \
                       $thread_controller
done

INFO "Comparing the two reports (row order may differ between profiles, so compare as multisets)"
python3 - <<'PYEOF'
import collections, os, sys

def records(path):
    if not path.endswith('.fa'):
        return collections.Counter(open(path).read().splitlines())
    recs, header, seq = [], None, []
    for line in open(path):
        line = line.rstrip('\n')
        if line.startswith('>'):
            if header is not None:
                recs.append((header, ''.join(seq)))
            header, seq = line, []
        elif line:
            seq.append(line)
    if header is not None:
        recs.append((header, ''.join(seq)))
    return collections.Counter(recs)

full_files, lean_files = set(os.listdir('INDELS_FULL')), set(os.listdir('INDELS_LEAN'))
if not full_files or full_files != lean_files:
    print(f"FAIL: report files differ: {sorted(full_files ^ lean_files)}")
    sys.exit(1)
for f in sorted(full_files):
    if records(f'INDELS_FULL/{f}') != records(f'INDELS_LEAN/{f}'):
        print(f"FAIL: {f} differs between the reports made from full and lean profiles")
        sys.exit(1)
    print(f"ok: {f} identical")
PYEOF

# expect_refusal NAME "MESSAGE FRAGMENT" COMMAND...: the command must fail and print the fragment
# (anvi'o wraps long messages across lines, so whitespace is ignored when matching)
expect_refusal() {
    local name=$1 fragment=$2; shift 2
    if "$@" > guard_${name}.log 2>&1
    then
        echo "FAIL: ${name}: the command succeeded but should have been refused"
        exit 1
    fi
    if ! tr -d ' \n\t' < guard_${name}.log | grep -qF -- "$(tr -d ' \n\t' <<< "$fragment")"
    then
        echo "FAIL: ${name}: refused, but not with the expected message ('$fragment'); see guard_${name}.log"
        exit 1
    fi
    echo "ok: ${name} refused as expected"
}

INFO "Guards: anvi-report-indels refuses --min-length below the profiles' --min-indel-length"
expect_refusal report_min_length "min-indel-length 100" \
    anvi-report-indels -c 02_CONTIGS/CONTIGS.db -p 04_PROFILE_LEAN/S01/PROFILE.db -p 04_PROFILE_LEAN/S02/PROFILE.db \
                       --min-samples 1 --min-length 50 -o INDELS_GUARD

INFO "Guards: anvi-profile refuses --profile-SCVs with --skip-SNV-storage"
expect_refusal profile_scvs "can't be combined with --skip-SNV-storage" \
    anvi-profile -c 02_CONTIGS/CONTIGS.db -i 03_MAPPING/S01.bam -o 04_PROFILE_GUARD --report-variability-full \
                 --profile-SCVs --skip-SNV-storage

INFO "Guards: anvi-merge refuses to merge lean and full profiles"
expect_refusal merge_mixed "skip-SNV-storage) are not identical" \
    anvi-merge 04_PROFILE_LEAN/S01/PROFILE.db 04_PROFILE_FULL/S02/PROFILE.db -c 02_CONTIGS/CONTIGS.db -o 04_PROFILE_MIXED

INFO "Guards: anvi-gen-variability-profile refuses a profile without stored SNV rows"
python3 -c "import sqlite3; print('\n'.join(r[0] for r in sqlite3.connect('02_CONTIGS/CONTIGS.db').execute('SELECT split FROM splits_basic_info')))" > splits_of_interest.txt
expect_refusal variability_lean "were not stored" \
    anvi-gen-variability-profile -c 02_CONTIGS/CONTIGS.db -p 04_PROFILE_LEAN/S01/PROFILE.db \
                                 --splits-of-interest splits_of_interest.txt -o variability_lean.txt

INFO "Control: the same anvi-gen-variability-profile call works on the full profile"
anvi-gen-variability-profile -c 02_CONTIGS/CONTIGS.db -p 04_PROFILE_FULL/S01/PROFILE.db \
                             --splits-of-interest splits_of_interest.txt -o variability_full.txt

INFO "Lean-profiling component test passed"
