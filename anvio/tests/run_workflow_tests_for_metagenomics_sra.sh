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
config['download_reads']['max_disk_gb'] = 0.5
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

# Anvi'o either knows the technology of every long-read sample or of none of them. Knowing
# some would mean processing the rest with another sample's presets, so it is refused -- and
# when it knows none, it must not leave behind an lr_technology column that merely looks filled
# in, because the checks that make sure config presets are set would then skip themselves.
INFO "A long-read technology anvi'o cannot determine leaves no lr_technology column behind"
DRY_RUN config-pacbio.json
if grep -q "lr_technology" 01_SRA-pacbio/samples-txt-with-downloaded-reads.txt; then
    echo "FAIL: an lr_technology column was declared when anvi'o knows no technologies at all."
    cat 01_SRA-pacbio/samples-txt-with-downloaded-reads.txt
    exit 1
fi

INFO "Knowing some long-read technologies but not all of them is refused"
printf 'sample\tlr\tlr_technology\tsra_accession\n' > samples-half-known.txt
printf 'S_local\tthree_samples_example/sample-01-LR.fastq.gz\tont\t\n' >> samples-half-known.txt
printf 'S05\t\t\tSRR11951439\n' >> samples-half-known.txt
$ANVIO_PYTHON -c "
import json
config = json.load(open('config-pacbio.json'))
config['samples_txt'] = 'samples-half-known.txt'
json.dump(config, open('config-half-known.json', 'w'), indent=4)
"
EXPECT_FAIL "a set of long-read samples whose technologies are only partly known" \
    anvi-run-workflow -w metagenomics -c config-half-known.json --dry-run

# ... and once the gap is filled, what the user declared by hand is what gets used.
INFO "A technology declared by hand survives into the derived samples-txt"
cp SRA-METADATA.txt SRA-METADATA-filled.txt
$ANVIO_PYTHON -c "
import io
path = 'SRA-METADATA-filled.txt'
text = io.open(path).read().replace('SRR11951439\tLR\t\t', 'SRR11951439\tLR\tpb-hifi\t')
io.open(path, 'w').write(text)

import json
config = json.load(open('config-half-known.json'))
config['download_reads']['metadata_cache'] = 'SRA-METADATA-filled.txt'
config['output_dirs'] = {k: v.replace('-pacbio', '-halfknown') for k, v in config['output_dirs'].items()}
json.dump(config, open('config-half-known-filled.json', 'w'), indent=4)
"
DRY_RUN config-half-known-filled.json
ASSERT_FILE_CONTAINS 01_SRA-halfknown/samples-txt-with-downloaded-reads.txt "ont"
ASSERT_FILE_CONTAINS 01_SRA-halfknown/samples-txt-with-downloaded-reads.txt "pb-hifi"

INFO "One sample can be made of a short-read run and a long-read run at once"
DRY_RUN config-hybrid.json
ASSERT_FILE_CONTAINS 01_SRA-hybrid/samples-txt-with-downloaded-reads.txt "S06_R1.fastq.gz"
ASSERT_FILE_CONTAINS 01_SRA-hybrid/samples-txt-with-downloaded-reads.txt "S06_LR.fastq.gz"

INFO "Downloaded reads and reads that were already on disk can be mixed in one samples-txt"
DRY_RUN config-mixed-with-local.json > mixed-output.txt 2>&1
ASSERT_FILE_CONTAINS 01_SRA-mixed-with-local/samples-txt-with-downloaded-reads.txt "sample-01-R1.fastq.gz"

# Quality filtering is one rule for the whole run, and snakemake decides whether a rule's output
# is temporary before it knows which sample it is working on. So the filtered reads of a sample
# that was never downloaded go along with everyone else's, and anvi'o says so rather than letting
# it be discovered afterwards.
ASSERT_FILE_CONTAINS mixed-output.txt "FILTERED READS GO FOR YOUR OWN SAMPLES TOO"

# A sample with nothing to download has nothing to release, so it gets no release unit of its own.
INFO "A sample that was already on disk does not become a release unit"
SNAKEMAKE_DAG config-mixed-with-local.json > dag-mixed.txt
ASSERT_FILE_CONTAINS dag-mixed.txt "01_SRA-mixed-with-local/released/S01.released"
if grep -q "released/S_local.released" dag-mixed.txt; then
    echo "FAIL: a sample whose reads were already on disk was given a release unit."
    exit 1
fi

# Reads that are kept are not deleted as the workflow goes, so they pile up for the whole run
# and a disk budget has nothing to say about them. Anvi'o says so, with a number.
INFO "Asking anvi'o to keep the reads gets a warning about what they will cost"
$ANVIO_PYTHON -c "
import json
config = json.load(open('config-references.json'))
config['download_reads']['keep_reads'] = 'both'
json.dump(config, open('config-keep-reads.json', 'w'), indent=4)
"
DRY_RUN config-keep-reads.json > keep-reads-output.txt 2>&1
ASSERT_FILE_CONTAINS keep-reads-output.txt "THE READS YOU ARE KEEPING NEED ROOM OF THEIR OWN"

# Downloading reads is what makes quality-filtered reads disposable. Reference-based read
# removal makes that true of SHORT reads on its own, but it is a short-read step and has no say
# over what filtlong makes, so a run with no SRA in it anywhere must keep its filtered long
# reads -- which is what it did before this workflow learned to download anything.
INFO "A run with no SRA accessions keeps its filtered long reads, even with reference removal"
printf 'sample\tr1\tr2\tlr\n' > samples-no-sra.txt
printf 'S_sr\tthree_samples_example/sample-01-R1.fastq.gz\tthree_samples_example/sample-01-R2.fastq.gz\t\n' >> samples-no-sra.txt
printf 'S_lr\t\t\tthree_samples_example/sample-01-LR.fastq.gz\n' >> samples-no-sra.txt
$ANVIO_PYTHON -c "
import json
config = json.load(open('config-references.json'))
config['samples_txt'] = 'samples-no-sra.txt'
config['references_mode'] = True
config['filtlong'] = {'run': True, 'use_anvio_conda_yaml': True, '--min-length': 500}
config['minimap2'] = dict(config.get('minimap2', {}), preset='map-ont')
config['remove_short_reads_based_on_references'] = {
    'threads': 1, 'dont_remove_just_map': None,
    'references_for_removal_txt': 'references-for-removal.txt',
    'delimiter-for-iu-remove-ids-from-fastq': ' '}
config['output_dirs'] = {k: v.replace('-references', '-nosra') for k, v in config['output_dirs'].items()}
json.dump(config, open('config-no-sra-lr.json', 'w'), indent=4)
"
SNAKEMAKE_DAG config-no-sra-lr.json > dag-no-sra.txt
if grep -q "Would remove temporary output.*FILTERED_LR" dag-no-sra.txt; then
    echo "FAIL: filtered long reads are deleted in a run that downloads nothing."
    grep "Would remove temporary output.*FILTERED_LR" dag-no-sra.txt
    exit 1
fi
# ... while the short reads of that same run are superseded by the removal step, as always.
ASSERT_FILE_CONTAINS dag-no-sra.txt "Would remove temporary output 01_QC-nosra"

# Samples that are co-assembled have to be on disk at the same time, so they are released
# together rather than one by one.
INFO "Co-assembled samples share a single release unit"
SNAKEMAKE_DAG config-co-assembly.json > dag-co-assembly.txt
ASSERT_FILE_CONTAINS dag-co-assembly.txt "01_SRA-co-assembly/released/CO.released"
# ... and released together means released once: neither of them gets a unit of its own.
if grep -qE "released/S0[12]\.released" dag-co-assembly.txt; then
    echo "FAIL: a co-assembled sample was given a release unit of its own."
    grep -E "released/S0[12]\.released" dag-co-assembly.txt
    exit 1
fi

# An accession can be named by more than one sample. The download happens once, so it has to
# wait for the earliest unit that wants it -- waiting on a later one deadlocks, because that
# unit's gate can be the earlier unit, whose reads this very download is needed to produce.
INFO "An accession shared by samples in different groups does not deadlock the workflow"
printf 'sample\tgroup\tsra_accession\n' > samples-shared-accession.txt
printf 'A\tG1\tERR6450080\n' >> samples-shared-accession.txt
printf 'B\tG2\tERR6450081\n' >> samples-shared-accession.txt
printf 'C\tG1\tERR6450081\n' >> samples-shared-accession.txt
$ANVIO_PYTHON -c "
import json
config = json.load(open('config-co-assembly.json'))
config['samples_txt'] = 'samples-shared-accession.txt'
config['download_reads']['max_disk_gb'] = 0.9
config['output_dirs'] = {k: v.replace('-co-assembly', '-shared') for k, v in config['output_dirs'].items()}
json.dump(config, open('config-shared-accession.json', 'w'), indent=4)
"
# A cyclic graph makes snakemake exit non-zero, which `set -e` would turn into an abort with no
# explanation, so let it fail and let the check below say what actually went wrong.
SNAKEMAKE_DAG config-shared-accession.json > dag-shared.txt 2>&1 || true
if grep -qi "cyclic" dag-shared.txt; then
    echo "FAIL: an accession named by two samples produced a cyclic workflow."
    grep -i -A3 "cyclic" dag-shared.txt
    exit 1
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

# A platform anvi'o has not learned about does not stop the lookup: the run is written down
# with everything NCBI did say and an empty read_type, so that the row this error asks you to
# edit is one that exists. The workflow then refuses until someone fills it in.
INFO "An unfamiliar sequencing platform is written down with an empty read_type"
cp SRA-METADATA.txt SRA-METADATA-unknown.txt
printf 'SRR9999999\t\t\tHELICOS\tHeliscope\tSINGLE\t1000\t300000\t1\tncbi\n' >> SRA-METADATA-unknown.txt
printf 'sample\tsra_accession\nS07\tSRR9999999\n' > samples-unknown-platform.txt
$ANVIO_PYTHON -c "
import json
config = json.load(open('config-references.json'))
config['samples_txt'] = 'samples-unknown-platform.txt'
config['download_reads']['metadata_cache'] = 'SRA-METADATA-unknown.txt'
json.dump(config, open('config-unknown-platform.json', 'w'), indent=4)
"
EXPECT_FAIL "an accession whose sequencing platform anvi'o does not recognize" \
    anvi-run-workflow -w metagenomics -c config-unknown-platform.json --dry-run

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
