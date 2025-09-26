#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2 $3
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
cp $files/data/input_files/hmm_list_group.txt                           $output_dir/workflow_test
cp $files/data/input_files/hmm_list_external.txt                        $output_dir/workflow_test
cd $output_dir/workflow_test

INFO "Migrating all databases"
anvi-migrate *db --migrate-quickly

INFO "Creating a default config for ecophylo workflow"
anvi-run-workflow -w ecophylo --get-default-config default-config.json
sed 's|external-genomes.txt||' default-config.json > only-metagenomes-txt-config.json
sed 's|metagenomes.txt||' default-config.json > only-external-genomes-txt-config.json
sed 's|samples\.txt||' default-config.json > no-samples-txt-config.json
sed 's|\"AA_mode\"\: false|\"AA_mode\"\: true|' no-samples-txt-config.json > AA-mode-config.json
sed 's|external-genomes.txt||' default-config.json | sed 's|samples\.txt||' > no-samples-only-metagenomes-txt-config.json
sed 's|metagenomes.txt||' default-config.json | sed 's|samples\.txt||' > no-samples-only-external-genomes-txt-config.json
sed 's|"run_genomes_sanity_check": true|"run_genomes_sanity_check": false|' default-config.json > no-genomes-sanity-check-config.json
awk '{
    if ($0 ~ /hmm_list.txt/) {
        sub(/hmm_list.txt/, "hmm_list_group.txt", $0)
    }
    if ($1 ~ /anvi_run_scg_taxonomy/) {
        print
        getline
        sub(/true/, "false", $0)
    }
    print
}' default-config.json > merge-by-group-config.json


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

INFO "Listing dependencies for ecophylo workflow (profile-mode)"
anvi-run-workflow -w ecophylo -c default-config.json --list-dependencies

INFO "Saving a workflow graph (profile-mode)"
anvi-run-workflow -w ecophylo -c default-config.json --save-workflow-graph

INFO "Running ecophylo workflow with ecophylo dry-run (profile-mode)"
anvi-run-workflow -w ecophylo -c default-config.json -A --dry-run

INFO "Running ecophylo workflow with ecophylo dry-run: only metagenomes.txt (profile-mode)"
anvi-run-workflow -w ecophylo -c only-metagenomes-txt-config.json -A --dry-run

INFO "Running ecophylo workflow with ecophylo dry-run: only external-genomes.txt (profile-mode)"
anvi-run-workflow -w ecophylo -c only-external-genomes-txt-config.json -A --dry-run

INFO "Running ecophylo workflow with ecophylo dry-run: no GenomeDescriptions and MetagenomeDescriptions"
anvi-run-workflow -w ecophylo -c no-genomes-sanity-check-config.json -A --dry-run

INFO "Running ecophylo workflow (profile-mode)"
anvi-run-workflow -w ecophylo -c default-config.json

INFO "Running ecophylo workflow interactive (profile-mode)"
HMM=`awk 'NR==2{print $2 "_" $1}' hmm_list.txt`
echo $HMM
anvi-interactive -c ECOPHYLO_WORKFLOW/METAGENOMICS_WORKFLOW/03_CONTIGS/${HMM}-contigs.db \
                 -p ECOPHYLO_WORKFLOW/METAGENOMICS_WORKFLOW/06_MERGED/${HMM}/PROFILE.db \
                 $dry_run_controller

rm -rf $output_dir/workflow_test/ECOPHYLO_WORKFLOW/

INFO "Saving a workflow graph - no samples.txt (tree-mode)"
anvi-run-workflow -w ecophylo -c no-samples-txt-config.json --save-workflow-graph

INFO "Running ecophylo workflow with ecophylo dry-run - no samples.txt (tree-mode)"
anvi-run-workflow -w ecophylo -c no-samples-txt-config.json -A --dry-run

INFO "Running ecophylo workflow with ecophylo dry-run - no samples.txt, only metagenomes.txt (tree-mode)"
anvi-run-workflow -w ecophylo -c no-samples-only-metagenomes-txt-config.json -A --dry-run

INFO "Running ecophylo workflow with ecophylo dry-run - no samples.txt, only external-genomes.txt (tree-mode)"
anvi-run-workflow -w ecophylo -c no-samples-only-external-genomes-txt-config.json -A --dry-run

INFO "Running ecophylo workflow - no samples.txt and AA-mode (tree-mode)"
anvi-run-workflow -w ecophylo -c AA-mode-config.json -A --dry-run

INFO "Running ecophylo workflow - no samples.txt (tree-mode)"
rm -rf $output_dir/workflow_test/ECOPHYLO_WORKFLOW/
anvi-run-workflow -w ecophylo -c no-samples-txt-config.json

INFO "Running ecophylo workflow interactive (tree-mode)"
anvi-interactive -t ECOPHYLO_WORKFLOW/05_TREES/${HMM}/${HMM}_renamed.nwk \
                 -p ECOPHYLO_WORKFLOW/05_TREES/${HMM}/${HMM}-PROFILE.db \
                 --manual \
                 $dry_run_controller

rm -rf $output_dir/workflow_test/ECOPHYLO_WORKFLOW/

INFO "Running ecophylo workflow - external HMM (tree-mode)"
anvi-run-workflow -w ecophylo -c no-samples-only-external-genomes-txt-config.json

INFO "Running ecophylo workflow interactive from external HMM (tree-mode)"
anvi-interactive -t ECOPHYLO_WORKFLOW/05_TREES/${HMM}/${HMM}_renamed.nwk \
                 -p ECOPHYLO_WORKFLOW/05_TREES/${HMM}/${HMM}-PROFILE.db \
                 --manual \
                 $dry_run_controller

rm -rf $output_dir/workflow_test/ECOPHYLO_WORKFLOW/

INFO "Running ecophylo workflow - merge by group (profile mode)"
anvi-run-workflow -w ecophylo -c merge-by-group-config.json

INFO "Running ecophylo workflow interactive (merge by group - profile mode)"
GROUP=`awk 'NR==2{print $4}' hmm_list_group.txt`
anvi-interactive -c ECOPHYLO_WORKFLOW/METAGENOMICS_WORKFLOW/03_CONTIGS/${GROUP}-contigs.db \
                 -p ECOPHYLO_WORKFLOW/METAGENOMICS_WORKFLOW/06_MERGED/${GROUP}/PROFILE.db \
                 $dry_run_controller

