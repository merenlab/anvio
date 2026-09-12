#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2 $3
#####################################

# This tests the metagenomics workflow's ability to download its own reads from the SRA and
# delete them as it goes. Almost everything here runs offline: a metadata file describing the SRA
# runs is shipped as a test fixture, so anvi'o never needs to ask NCBI anything and the scenarios
# below are all dry runs that finish in seconds. Only the last section actually downloads reads.

INFO "Setting up the metagenomics SRA workflow test directory"
mkdir $output_dir/workflow_test
cp -r $files/workflows/metagenomics_sra/* $output_dir/workflow_test/
cp -r $files/workflows/metagenomics/three_samples_example $output_dir/workflow_test/
cd $output_dir/workflow_test

# `anvi-run-workflow --dry-run` does not show which job depends on which, and that is precisely
# what has to be checked here, so some of the tests below ask snakemake directly. The interpreter
# comes from the anvi-run-workflow entry point rather than from whatever `python` happens to be
# on the PATH, which is not always the one anvi'o was installed with.
ANVIO_PYTHON=$(head -1 "$(command -v anvi-run-workflow)" | sed 's/^#!//')
ANVIO_SNAKEMAKE="$(dirname $ANVIO_PYTHON)/snakemake"
METAGENOMICS_SNAKEFILE=$($ANVIO_PYTHON -c "import anvio.workflows as w; print(w.get_workflow_snake_file_path('metagenomics'))")

DRY_RUN() {
    anvi-run-workflow -w metagenomics -c $1 --dry-run
}

SNAKEMAKE_DAG() {
    $ANVIO_SNAKEMAKE --snakefile $METAGENOMICS_SNAKEFILE --configfile $1 --dryrun --cores 4
}

INFO "A samples-txt with SRA accessions builds a workflow, even though no reads exist yet"
DRY_RUN config-references.json

INFO "The reads anvi'o will download are described in a samples-txt of its own"
SHOW_FILE 01_SRA-references/samples-txt-with-downloaded-reads.txt

# The disk budget is the point of the whole feature. With a budget big enough for both samples
# there is nothing to wait for, and with a budget that only fits one, the second sample's
# download must wait for the first sample's reads to have been released.
INFO "With room for every sample, no download waits for another"
SNAKEMAKE_DAG config-references.json > dag-roomy.txt
if grep -q "input: 01_SRA-references/released" dag-roomy.txt; then
    echo "FAIL: a download is waiting for another even though the budget has room for every sample."
    exit 1
fi

INFO "With room for one sample at a time, the second download waits for the first"
$ANVIO_PYTHON -c "
import json
config = json.load(open('config-references.json'))
config['download_reads']['max_disk_gb'] = 1.0
json.dump(config, open('config-one-at-a-time.json', 'w'), indent=4)
"
SNAKEMAKE_DAG config-one-at-a-time.json > dag-tight.txt
ASSERT_FILE_CONTAINS dag-tight.txt "01_SRA-references/released/S01.released"

INFO "Long reads from a nanopore run are recognized, and so is their technology"
DRY_RUN config-nanopore.json
ASSERT_FILE_CONTAINS 01_SRA-nanopore/samples-txt-with-downloaded-reads.txt "ont"
ASSERT_FILE_CONTAINS 01_SRA-nanopore/samples-txt-with-downloaded-reads.txt "S04_LR.fastq.gz"

INFO "A PacBio run whose chemistry cannot be told from its metadata gets a warning, not an error"
DRY_RUN config-pacbio.json > pacbio-output.txt 2>&1
ASSERT_FILE_CONTAINS pacbio-output.txt "SRR11951439"

INFO "One sample can be made of a short-read run and a long-read run at once"
DRY_RUN config-hybrid.json
ASSERT_FILE_CONTAINS 01_SRA-hybrid/samples-txt-with-downloaded-reads.txt "S06_R1.fastq.gz"
ASSERT_FILE_CONTAINS 01_SRA-hybrid/samples-txt-with-downloaded-reads.txt "S06_LR.fastq.gz"

INFO "Downloaded reads and reads that were already on disk can be mixed in one samples-txt"
DRY_RUN config-mixed-with-local.json
ASSERT_FILE_CONTAINS 01_SRA-mixed-with-local/samples-txt-with-downloaded-reads.txt "sample-01-R1.fastq.gz"

# Samples that are co-assembled have to be on disk at the same time, so they are released
# together rather than one by one.
INFO "Co-assembled samples share a single release unit"
SNAKEMAKE_DAG config-co-assembly.json > dag-co-assembly.txt
ASSERT_FILE_CONTAINS dag-co-assembly.txt "01_SRA-co-assembly/released/CO.released"
if grep -c "sra_reads_released *1$" dag-co-assembly.txt > /dev/null; then
    echo "Both samples of the co-assembly group are released together, as they should be."
fi

# Mapping every sample against every assembly leaves anvi'o no room to download a few samples at
# a time: nothing can be mapped until every assembly exists, and no assembly exists until its own
# reads have been downloaded. That is allowed — the reads are still deleted afterwards — but all
# the samples end up in a single release unit, and the budget can then only say yes or no.
INFO "Mapping everything against everything while assembling puts every sample in one release unit"
DRY_RUN config-all-against-all.json > all-against-all-output.txt 2>&1
ASSERT_FILE_CONTAINS all-against-all-output.txt "EVERY SAMPLE WILL BE ON DISK AT ONCE"
SNAKEMAKE_DAG config-all-against-all.json > dag-all-against-all.txt
ASSERT_FILE_CONTAINS dag-all-against-all.txt "01_SRA-all-against-all/released/every-sample-at-once.released"

EXPECT_FAIL "single-end short reads, which this workflow cannot process" \
    anvi-run-workflow -w metagenomics -c config-single-end.json --dry-run

EXPECT_FAIL "a disk budget too small to fit even one sample" \
    anvi-run-workflow -w metagenomics -c config-tiny-budget.json --dry-run

EXPECT_FAIL "a budget too small for a run where every sample has to be downloaded at once" \
    anvi-run-workflow -w metagenomics -c config-all-against-all-tight.json --dry-run

EXPECT_FAIL "removing reads based on references while quality filtering is off" \
    anvi-run-workflow -w metagenomics -c config-ref-removal-without-qc.json --dry-run

# ---------------------------------------------------------------------------------------------
# Everything below actually downloads reads from NCBI. Two small metagenomes are mapped against
# a small reference and profiled, with a disk budget that only has room for one of them at a
# time — so the second sample cannot be downloaded until the first one's reads are gone.
# ---------------------------------------------------------------------------------------------

ASSERT_FILE_EXISTS() {
    if [ ! -f "$1" ]; then
        echo "FAIL: expected the file '$1' to exist, but it does not."
        exit 1
    fi
}

# The claim this whole feature makes is that only so many metagenomes are ever on disk at once.
# Nothing about the workflow's output can show that after the fact, so we watch while it runs.
INFO "Watching how many samples have reads on disk while the workflow runs"
# Only files that hold reads are counted. Both directories are also full of small per-sample
# files that stay around for the whole run (.ini files, statistics, read counts), and counting
# those would mean always seeing every sample no matter what the workflow was doing.
cat > watch_residency.sh << 'EOF'
#!/bin/bash
while true; do
    ls 01_SRA-real/reads/ 01_QC-real/ 2>/dev/null \
        | grep '\.fastq' | grep -oE "^S[0-9]+" | sort -u | wc -l >> residency.txt
    sleep 1
done
EOF
chmod +x watch_residency.sh
touch residency.txt
./watch_residency.sh &
watcher_pid=$!

INFO "Downloading, mapping and profiling two metagenomes from the SRA"
anvi-run-workflow -w metagenomics -c config-real.json

kill $watcher_pid 2>/dev/null || true

INFO "Making sure only one sample's reads were ever on disk at a time"
most_at_once=$(sort -n residency.txt | tail -1)
echo "The most samples that had reads on disk at any one moment: $most_at_once"
if [ "$most_at_once" -gt 1 ]; then
    echo "FAIL: the disk budget only had room for one sample at a time, and yet $most_at_once"
    echo "samples had reads on disk at once."
    exit 1
fi

INFO "Making sure no downloaded reads were left behind"
if [ -n "$(find 01_SRA-real -name '*.fastq*' 2>/dev/null)" ]; then
    echo "FAIL: downloaded reads are still on disk:"
    find 01_SRA-real -name '*.fastq*'
    exit 1
fi
if [ -n "$(find 01_SRA-real -name '*.sra' -o -name '*.sralite' 2>/dev/null)" ]; then
    echo "FAIL: downloaded SRA archives are still on disk."
    exit 1
fi

INFO "Making sure the workflow produced what it was supposed to produce"
ASSERT_FILE_EXISTS 05_ANVIO_PROFILE-real/SRAREF/S01/PROFILE.db
ASSERT_FILE_EXISTS 05_ANVIO_PROFILE-real/SRAREF/S02/PROFILE.db
ASSERT_FILE_EXISTS 06_MERGED-real/SRAREF/PROFILE.db

INFO "Making sure a second run does not download anything all over again"
anvi-run-workflow -w metagenomics -c config-real.json
if [ -n "$(find 01_SRA-real -name '*.fastq*' 2>/dev/null)" ]; then
    echo "FAIL: the second run of the workflow downloaded the reads again."
    exit 1
fi

INFO "The SRA metadata anvi'o worked with"
SHOW_FILE SRA-METADATA.txt
