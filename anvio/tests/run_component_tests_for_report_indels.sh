#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2 $3
#####################################

# Skip gracefully if long-read mapping tools are missing.
if ! command -v minimap2 >/dev/null 2>&1; then
    INFO "minimap2 not installed -- skipping anvi-report-indels component test"
    exit 0
fi
if ! command -v samtools >/dev/null 2>&1; then
    INFO "samtools not installed -- skipping anvi-report-indels component test"
    exit 0
fi

INFO "Setting up the report-indels mock-MGE analysis directory"
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

INFO "Building reference contigs-db from S01 (transposase at its original location)"
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

INFO "Profiling each BAM (indels are profiled by default)"
for sample in S01 S02 S03
do
    anvi-profile -c 02_CONTIGS/CONTIGS.db \
                 -i 03_MAPPING/${sample}.bam \
                 -o 04_PROFILE/${sample} \
                 --sample-name ${sample} \
                 --min-contig-length 1000 \
                 --cluster-contigs \
                 --report-variability-full \
                 --quiet \
                 $thread_controller
done

INFO "Merging the per-sample profile databases"
anvi-merge 04_PROFILE/S01/PROFILE.db \
           04_PROFILE/S02/PROFILE.db \
           04_PROFILE/S03/PROFILE.db \
           -c 02_CONTIGS/CONTIGS.db \
           -o 04_PROFILE/MERGED \
           -W \
           --quiet

INFO "Running anvi-report-indels with --min-length 400 (the mock transposase is 435 bp) and --min-samples 1 (each INS site is unique to one sample in this 3-variant mock)"
anvi-report-indels -c 02_CONTIGS/CONTIGS.db \
                   -p 04_PROFILE/MERGED/PROFILE.db \
                   --min-samples 1 \
                   --min-length 400 \
                   -o INDELS_DEFAULT \
                   $thread_controller

if [ -e INDELS_DEFAULT/.anvi-report-indels-tmp ]
then
    echo "anvi-report-indels left its work directory (the on-disk indel sequences) behind after a successful run"
    exit 1
fi

SHOW_FILE INDELS_DEFAULT/INDELS-VIEW-DATA.txt
SHOW_FILE INDELS_DEFAULT/INDELS-ITEMS-ADDITIONAL.txt
SHOW_FILE INDELS_DEFAULT/INDELS-LOCI-LONG.txt
SHOW_FILE INDELS_DEFAULT/INDELS-CLUSTERS.txt

INFO "The cluster representatives FASTA must be usable as is by anvi-gen-contigs-database"
anvi-gen-contigs-database -f INDELS_DEFAULT/INDELS-CLUSTERS.fa \
                          -o INDELS_DEFAULT/CLUSTERS-CONTIGS.db \
                          --skip-gene-calling \
                          --quiet

INFO "Running anvi-report-indels against multiple single profile-dbs via -p (repeated)"
anvi-report-indels -c 02_CONTIGS/CONTIGS.db \
                   -p 04_PROFILE/S01/PROFILE.db \
                   -p 04_PROFILE/S02/PROFILE.db \
                   -p 04_PROFILE/S03/PROFILE.db \
                   --min-samples 1 \
                   --min-length 400 \
                   -o INDELS_MULTI_P \
                   $thread_controller

INFO "Running anvi-report-indels with --event-type DEL only"
anvi-report-indels -c 02_CONTIGS/CONTIGS.db \
                   -p 04_PROFILE/MERGED/PROFILE.db \
                   --event-type DEL \
                   --min-samples 1 \
                   --min-length 400 \
                   -o INDELS_DEL_ONLY \
                   $thread_controller

INFO "Running anvi-report-indels with --event-type INS only"
anvi-report-indels -c 02_CONTIGS/CONTIGS.db \
                   -p 04_PROFILE/MERGED/PROFILE.db \
                   --event-type INS \
                   --min-samples 1 \
                   --min-length 400 \
                   -o INDELS_INS_ONLY \
                   $thread_controller

INFO "Reproducibility: the report must not depend on the order of rows in the indels tables"
# The profiler stores indel rows in whatever order its worker processes finish contigs, so two
# profiles of the same BAM can hold the same rows in different orders. Shuffle the rows in copies
# of the single-sample profile-dbs and require byte-identical report outputs.
mkdir -p 04_PROFILE_SHUFFLED
for sample in S01 S02 S03
do
    mkdir -p 04_PROFILE_SHUFFLED/${sample}
    cp 04_PROFILE/${sample}/PROFILE.db 04_PROFILE/${sample}/AUXILIARY-DATA.db 04_PROFILE_SHUFFLED/${sample}/
done
python3 - <<'PYEOF'
import random, sqlite3
random.seed(1)
for sample in ['S01', 'S02', 'S03']:
    db = sqlite3.connect(f'04_PROFILE_SHUFFLED/{sample}/PROFILE.db')
    rows = db.execute('SELECT * FROM indels').fetchall()
    original = list(rows)
    random.shuffle(rows)
    if len(rows) > 1 and rows == original:
        rows = rows[1:] + rows[:1]   # make sure the order actually changed
    db.execute('DELETE FROM indels')
    if rows:
        db.executemany(f"INSERT INTO indels VALUES ({','.join('?' * len(rows[0]))})", rows)
    db.commit()
    db.close()
PYEOF
anvi-report-indels -c 02_CONTIGS/CONTIGS.db \
                   -p 04_PROFILE_SHUFFLED/S01/PROFILE.db \
                   -p 04_PROFILE_SHUFFLED/S02/PROFILE.db \
                   -p 04_PROFILE_SHUFFLED/S03/PROFILE.db \
                   --min-samples 1 \
                   --min-length 400 \
                   -o INDELS_MULTI_P_SHUFFLED \
                   $thread_controller

for f in INDELS-VIEW-DATA.txt INDELS-ITEMS-ADDITIONAL.txt INDELS-LOCI-LONG.txt INDELS-CLUSTERS.fa INDELS-CLUSTERS.txt INDELS-ITEMS-ORDER.newick
do
    if ! cmp -s INDELS_MULTI_P/$f INDELS_MULTI_P_SHUFFLED/$f
    then
        echo "anvi-report-indels output $f depends on the order of indel rows in the profile-dbs"
        exit 1
    fi
done
INFO "Reproducibility check passed: identical outputs after shuffling indel rows"

INFO "Scale options: --max-loci-for-tree skips only the tree"
anvi-report-indels -c 02_CONTIGS/CONTIGS.db \
                   -p 04_PROFILE/MERGED/PROFILE.db \
                   --min-samples 1 \
                   --min-length 400 \
                   --max-loci-for-tree 1 \
                   -o INDELS_NO_TREE \
                   $thread_controller

if [ -e INDELS_NO_TREE/INDELS-ITEMS-ORDER.newick ]
then
    echo "anvi-report-indels wrote a tree although there are more loci than --max-loci-for-tree"
    exit 1
fi
for f in INDELS-VIEW-DATA.txt INDELS-ITEMS-ADDITIONAL.txt INDELS-LOCI-LONG.txt INDELS-CLUSTERS.fa INDELS-CLUSTERS.txt
do
    if ! cmp -s INDELS_DEFAULT/$f INDELS_NO_TREE/$f
    then
        echo "--max-loci-for-tree changed $f, but it must only skip the tree"
        exit 1
    fi
done

INFO "Scale options: --exclude-contig-edges sets aside the INS rows near contig ends, and only those"
# the mock genomes have no INS near a contig end, so the edge distance is set just past the INS row
# closest to an end: that row (and any other as close) must be set aside, all others kept
EDGE=$(python3 - <<'PYEOF'
import csv, sqlite3
lengths = dict(sqlite3.connect('02_CONTIGS/CONTIGS.db').execute('SELECT contig, length FROM contigs_basic_info'))
rows = [r for r in csv.DictReader(open('INDELS_DEFAULT/INDELS-LOCI-LONG.txt'), delimiter='\t') if r['type'] == 'INS']
print(min(min(int(r['pos_in_contig']), lengths[r['contig_name']] - int(r['pos_in_contig'])) for r in rows) + 1)
PYEOF
)
INFO "Edge distance used: ${EDGE} bp"
anvi-report-indels -c 02_CONTIGS/CONTIGS.db \
                   -p 04_PROFILE/MERGED/PROFILE.db \
                   --min-samples 1 \
                   --min-length 400 \
                   --exclude-contig-edges ${EDGE} \
                   -o INDELS_NO_EDGES \
                   $thread_controller

EDGE=${EDGE} python3 - <<'PYEOF'
import csv, os, sqlite3, sys

EDGE = int(os.environ['EDGE'])
db = sqlite3.connect('02_CONTIGS/CONTIGS.db')
lengths = dict(db.execute('SELECT contig, length FROM contigs_basic_info'))
near_end = lambda r: int(r['pos_in_contig']) < EDGE or int(r['pos_in_contig']) >= lengths[r['contig_name']] - EDGE

read = lambda path: list(csv.DictReader(open(path), delimiter='\t')) if os.path.exists(path) else []
default_rows = read('INDELS_DEFAULT/INDELS-LOCI-LONG.txt')
kept_rows = read('INDELS_NO_EDGES/INDELS-LOCI-LONG.txt')
edge_rows = read('INDELS_NO_EDGES/INDELS-CONTIG-EDGE-INS.txt')

errors = []
if any(r['type'] != 'INS' or not near_end(r) for r in edge_rows):
    errors.append(f'a row set aside as contig-edge INS is not an INS within {EDGE} bp of a contig end')
if any(r['type'] == 'INS' and near_end(r) for r in kept_rows):
    errors.append(f'an INS row within {EDGE} bp of a contig end is still in the report')
if not edge_rows:
    errors.append('no row was set aside, but the edge distance was chosen so that at least one must be')
if len(kept_rows) + len(edge_rows) != len(default_rows):
    errors.append(f'rows kept ({len(kept_rows)}) + rows set aside ({len(edge_rows)}) != rows of the default run ({len(default_rows)})')
expected_edge = sum(1 for r in default_rows if r['type'] == 'INS' and near_end(r))
if len(edge_rows) != expected_edge:
    errors.append(f'{len(edge_rows)} rows set aside, but the default run has {expected_edge} INS rows near contig ends')

if errors:
    print('\n'.join(errors))
    sys.exit(1)
print(f'{len(edge_rows)} contig-edge INS rows set aside, {len(kept_rows)} rows kept')
PYEOF
INFO "Scale option checks passed"

INFO "Scale options: --sparse-matrix holds exactly the non-zero cells of the dense matrix"
anvi-report-indels -c 02_CONTIGS/CONTIGS.db \
                   -p 04_PROFILE/MERGED/PROFILE.db \
                   --min-samples 1 \
                   --min-length 400 \
                   --sparse-matrix \
                   -o INDELS_SPARSE \
                   $thread_controller

python3 - <<'PYEOF'
import csv, sys
dense = list(csv.reader(open('INDELS_DEFAULT/INDELS-VIEW-DATA.txt'), delimiter='\t'))
samples = dense[0][1:]
expected = sorted((r[0], s, v) for r in dense[1:] for s, v in zip(samples, r[1:]) if float(v) > 0)
sparse = list(csv.reader(open('INDELS_SPARSE/INDELS-VIEW-DATA-SPARSE.txt'), delimiter='\t'))
if sparse[0] != ['item', 'sample', 'presence']:
    sys.exit(f'unexpected header of the sparse matrix: {sparse[0]}')
observed = sorted(tuple(r) for r in sparse[1:])
if observed != expected:
    sys.exit(f'the sparse matrix ({len(observed)} cells) differs from the non-zero cells of the dense matrix ({len(expected)})')
print(f'{len(observed)} non-zero cells, identical in the sparse and the dense matrix')
PYEOF
for f in INDELS-ITEMS-ADDITIONAL.txt INDELS-LOCI-LONG.txt INDELS-CLUSTERS.fa INDELS-CLUSTERS.txt INDELS-ITEMS-ORDER.newick
do
    if ! cmp -s INDELS_DEFAULT/$f INDELS_SPARSE/$f
    then
        echo "--sparse-matrix changed $f, but it must only change how the matrix is written"
        exit 1
    fi
done
INFO "Sparse matrix checks passed"

INFO "Scale options: --mmseqs-split-memory-limit and --mmseqs-extra-params reach mmseqs2"
anvi-report-indels -c 02_CONTIGS/CONTIGS.db \
                   -p 04_PROFILE/MERGED/PROFILE.db \
                   --min-samples 1 \
                   --min-length 400 \
                   --mmseqs-split-memory-limit 30G \
                   --mmseqs-extra-params "--cluster-reassign 0" \
                   -o INDELS_MMSEQS_PARAMS \
                   $thread_controller

for f in INDELS-VIEW-DATA.txt INDELS-ITEMS-ADDITIONAL.txt INDELS-LOCI-LONG.txt INDELS-CLUSTERS.fa INDELS-CLUSTERS.txt INDELS-ITEMS-ORDER.newick
do
    if ! cmp -s INDELS_DEFAULT/$f INDELS_MMSEQS_PARAMS/$f
    then
        echo "a memory limit and a parameter set to its default value changed $f"
        exit 1
    fi
done

# a parameter mmseqs2 does not know must make the run fail, which shows that the parameters are passed on
if anvi-report-indels -c 02_CONTIGS/CONTIGS.db \
                      -p 04_PROFILE/MERGED/PROFILE.db \
                      --min-samples 1 \
                      --min-length 400 \
                      --mmseqs-extra-params "--not-an-mmseqs-parameter 1" \
                      -o INDELS_MMSEQS_BAD_PARAM \
                      $thread_controller > /dev/null 2>&1
then
    echo "anvi-report-indels ran although --mmseqs-extra-params held a parameter mmseqs2 does not know"
    exit 1
fi
INFO "MMseqs2 parameter checks passed"
