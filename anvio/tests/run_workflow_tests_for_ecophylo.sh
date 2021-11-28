#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2
#####################################

INFO "Setting up the ecophylo workflow test directory"
mkdir $output_dir/workflow_test
cp -r $files/workflows/metagenomics/*                                       $output_dir/workflow_test/
cp -r $files/workflows/ecophylo/*                                       $output_dir/workflow_test/
cp $files/data/input_files/metagenomes.txt                          $output_dir/workflow_test
cp $files/data/input_files/external-genomes.txt                     $output_dir/workflow_test
cp $files/data/input_files/hmm_list.txt                     $output_dir/workflow_test
cp $files/data/metagenomes/human_gut/IGD_SUBSET/CONTIGS.db                     $output_dir/workflow_test
cp $files/data/genomes/bacteria/*.db                    $output_dir/workflow_test
cp $files/data/genomes/archaea/*.db                     $output_dir/workflow_test
cd $output_dir/workflow_test

INFO "Creating a default config for ecophylo workflow"
anvi-run-workflow -w ecophylo --get-default-config default-config.json

# EcoPhylo does not use samples with groups
# sed 's|"samples_txt": ""|"samples_txt": "samples-no-groups.txt"|' default-config.json > ecophylo-config.json

INFO "Listing dependencies for ecophylo workflow"
anvi-run-workflow -w ecophylo -c default-config.json --list-dependencies

INFO "Saving a workflow graph"
anvi-run-workflow -w ecophylo -c default-config.json --save-workflow-graph

INFO "Running ecophylo workflow with ecophylo dry-run"
anvi-run-workflow -w ecophylo -c default-config.json --dry-run

INFO "Generating r1 and r2 short reads for samples"
for contig in sampleB_thetaiotamicron_Ribosomal_L16 sample_IGD_SUBSET_Ribosomal_L16 sample_P_marinus_Ribosomal_L16;
do 
    anvi-script-gen-short-reads "${contig}".ini --output-file-path $output_dir/workflow_test/"${contig}".fa
    anvi-script-gen-pseudo-paired-reads-from-fastq -f $output_dir/workflow_test/"${contig}".fa -O $output_dir/workflow_test/"${contig}";
done

rm samples.txt
echo -e "sample\tr1\tr2" > samples.txt
for contig in sampleB_thetaiotamicron_Ribosomal_L16 sample_IGD_SUBSET_Ribosomal_L16 sample_P_marinus_Ribosomal_L16;
do 
    echo -e "${contig}\t${contig}_1.fastq\t${contig}_2.fastq" >> samples.txt
done

INFO "Running ecophylo workflow"
anvi-run-workflow -w ecophylo -c default-config.json