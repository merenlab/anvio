#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2 $3
#####################################

# A deletion that crosses a split boundary must be recorded once, whole, in the split where it starts.

for program in minimap2 samtools
do
    if ! command -v $program >/dev/null 2>&1; then
        INFO "$program is not installed -- skipping this component test"
        exit 0
    fi
done

mkdir -p $output_dir
cp $files/mock_data_for_indels_across_splits/mock.py $output_dir/
cd $output_dir/

INFO "A 60 kb contig and its contigs database"
python3 mock.py reference
anvi-gen-contigs-database -f ref.fa -o CONTIGS.db -n SPLITS --quiet $thread_controller

INFO "Reads of a genome that lacks 1,500 bp across the first split boundary"
python3 mock.py reads CONTIGS.db
minimap2 -ax map-ont ref.fa reads.fa 2>/dev/null | samtools sort -o reads.bam - 2>/dev/null
samtools index reads.bam

INFO "Profiling"
anvi-profile -c CONTIGS.db -i reads.bam -o PROFILE -S mock --min-contig-length 1000 --report-variability-full --quiet $thread_controller
anvi-export-table PROFILE/PROFILE.db --table indels -o indels.txt --quiet

INFO "Checking the indels table"
python3 mock.py check indels.txt CONTIGS.db

INFO "Splitting the profile into bins that cut the contig at the boundary, and into a bin that keeps both sides"
python3 mock.py collections CONTIGS.db
anvi-import-collection cut.txt -p PROFILE/PROFILE.db -c CONTIGS.db -C CUT --quiet
anvi-import-collection kept.txt -p PROFILE/PROFILE.db -c CONTIGS.db -C KEPT --quiet
anvi-split -p PROFILE/PROFILE.db -c CONTIGS.db -C CUT -o SPLIT_CUT --quiet
anvi-split -p PROFILE/PROFILE.db -c CONTIGS.db -C KEPT -o SPLIT_KEPT --quiet

for bin in SPLIT_CUT/BEFORE SPLIT_CUT/AFTER SPLIT_KEPT/BOTH
do
    anvi-export-table $bin/PROFILE.db --table indels -o $bin/indels.txt --quiet
done

INFO "Checking the indels tables of the split profiles"
python3 mock.py check-split indels.txt CONTIGS.db SPLIT_CUT/BEFORE/indels.txt SPLIT_CUT/AFTER/indels.txt SPLIT_KEPT/BOTH/indels.txt
