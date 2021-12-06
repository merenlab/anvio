#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2
#####################################

INFO "Setting up the ecophylo workflow test directory"
mkdir $output_dir/workflow_test
cp -r $files/workflows/ecophylo/*                                       $output_dir/workflow_test/
cp $files/data/metagenomes/human_gut/IGD_SUBSET/CONTIGS.db              $output_dir/workflow_test
cp $files/data/genomes/bacteria/*.db                                    $output_dir/workflow_test
cp $files/data/genomes/archaea/*.db                                     $output_dir/workflow_test
cp $files/data/input_files/metagenomes.txt                              $output_dir/workflow_test
cp $files/data/input_files/external-genomes.txt                         $output_dir/workflow_test
cp $files/data/input_files/hmm_list.txt                                 $output_dir/workflow_test
cd $output_dir/workflow_test

INFO "Creating a default config for ecophylo workflow"
anvi-run-workflow -w ecophylo --get-default-config default-config.json

INFO "Generating r1 and r2 short reads for samples"
mkdir -p output
samples="sampleB_thetaiotamicron_Ribosomal_L16 sample_IGD_SUBSET_Ribosomal_L16 sample_P_marinus_Ribosomal_L16"

for sample in $samples
do
    ~/github/reads-for-assembly/gen-paired-end-reads samples/$sample.ini
done

echo -e "sample\tr1\tr2" > samples.txt
for reads in `ls output/*-R1.fastq`;
do 
    name=$(basename $reads -R1.fastq | sed 's|-|_|g')
    echo -e "${name}\t${reads}\t${reads/-R1.fastq/-R2.fastq}" >> samples.txt
done

INFO "Listing dependencies for ecophylo workflow"
anvi-run-workflow -w ecophylo -c default-config.json --list-dependencies

INFO "Saving a workflow graph"
anvi-run-workflow -w ecophylo -c default-config.json --save-workflow-graph

INFO "Running ecophylo workflow with ecophylo dry-run"
anvi-run-workflow -w ecophylo -c default-config.json --dry-run

INFO "Running ecophylo workflow"
anvi-run-workflow -w ecophylo -c default-config.json

INFO "Running ecophylo workflow interactive"
HMM="Ribosomal_L16"
anvi-interactive -c ECOPHYLO_WORKFLOW/METAGENOMICS_WORKFLOW/03_CONTIGS/"${HMM}"-contigs.db \
                 -p ECOPHYLO_WORKFLOW/METAGENOMICS_WORKFLOW/06_MERGED/"${HMM}"/PROFILE.db \
                 --state-autoload "${HMM}"