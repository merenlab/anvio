#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2
#####################################

INFO "Writing tRNA feature parameterization file"
anvi-trnaseq --default-feature-param-file ${output_dir}/DEFAULT-PARAMS.txt

INFO "Analyzing each sample FASTA of merged tRNA reads"
for sample in S01 S02
do
    for split in untreated demethylase
    do
        anvi-trnaseq -f ${files}/TRNASEQ-SAMPLE-${sample}_${split}.fa \
                     -S ${sample}_${split} \
                     -o ${output_dir}/${sample}_${split} \
                     --treatment ${split} \
                     --write-checkpoints \
                     --skip-fasta-check
        echo
    done
done

INFO "Converting tRNA-seq databases to contigs and profile databases"
anvi-convert-trnaseq-database ${output_dir}/S01_untreated/S01_untreated-TRNASEQ.db \
                              ${output_dir}/S01_demethylase/S01_demethylase-TRNASEQ.db \
                              ${output_dir}/S02_untreated/S02_untreated-TRNASEQ.db \
                              ${output_dir}/S02_demethylase/S02_demethylase-TRNASEQ.db \
                              -o ${output_dir}/CONVERTED \
                              -n TEST_EXPERIMENT \
                              --max-reported-trna-seeds 1000 \
                              --preferred-treatment demethylase \
                              --nonspecific-output nonspecific_db,combined_db,summed_db

INFO "Assigning taxonomy to tRNA seeds"
anvi-run-trna-taxonomy -c ${output_dir}/CONVERTED/CONTIGS.db \
                       --all-hits-output-file ${output_dir}/CONVERTED/TAXONOMY-HITS.txt

INFO "Firing up the interactive interface with the \"combined\" profile database"
anvi-interactive -p ${output_dir}/CONVERTED/COMBINED_COVERAGE/PROFILE.db \
                 -c ${output_dir}/CONVERTED/CONTIGS.db \
                 ${dry_run_controller}
