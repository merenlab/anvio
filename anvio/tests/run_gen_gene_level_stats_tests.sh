#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1
#####################################

# Set INSEQ files path
inseq_files=$files/example_files_for_inseq_tnseq

INFO "Initializing raw BAM files"
# init raw bam files.
for f in 01 02 03
do
    anvi-init-bam $inseq_files/INSEQ-$f-RAW.bam --output-file $output_dir/INSEQ-$f.bam
    echo
done

INFO "Generating a contigs database from a reference genome"
anvi-gen-contigs-database -f $inseq_files/contigs_inseq.fa \
                          -o $output_dir/CONTIGS.db \
                          --project-name "Contigs DB for generate gene level stats."

for f in 01 02 03
do
    INFO "Profiling sample SAMPLE-$f"
    anvi-profile -i $output_dir/INSEQ-$f.bam -o $output_dir/INSEQ-$f -c $output_dir/CONTIGS.db

    echo
done

INFO "Merging profiles"
anvi-merge $output_dir/INSEQ-*/PROFILE.db \
           -o $output_dir/SAMPLES-MERGED \
           -c $output_dir/CONTIGS.db

# Add a default collection
INFO "Adding default collection to merged profile"
anvi-script-add-default-collection -p $output_dir/SAMPLES-MERGED/PROFILE.db

# Generate gene-level-stats database
INFO "Computing gene level stats database"
anvi-gen-gene-level-stats-databases -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                                    -c $output_dir/CONTIGS.db \
                                    -C DEFAULT -b EVERYTHING
# Delete the GENE database
rm -r $output_dir/SAMPLES-MERGED/GENES/

# Compute INSeq stats database
INFO "Computing INSeq stats database"
anvi-gen-gene-level-stats-databases -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                                    -c $output_dir/CONTIGS.db  \
                                    -C DEFAULT -b EVERYTHING \
                                    --inseq-stats --just-do-it
