
#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2 $3
#####################################

INFO "Setting up the anvi-export-locus test directory"
mkdir $output_dir/export_locus_test
cp $files/data/genomes/bacteria/*.db                    $output_dir/export_locus_test
cp $files/data/genomes/archaea/*.db                     $output_dir/export_locus_test
cp $files/data/metagenomes/human_gut/IGD_SUBSET/*.db    $output_dir/export_locus_test
cp $files/data/input_files/*.txt                        $output_dir/export_locus_test

cd $output_dir/export_locus_test

mkdir test

INFO "Migrating all databases"
anvi-migrate *db --migrate-quickly

INFO "Running anvi-export-locus in default-mode"
anvi-export-locus -c CONTIGS.db --gene-caller-id 68 -n 7,7 -o test/ -O metagenome_68

INFO "Running anvi-export-locus in default-mode find gene with --search-term"
anvi-export-locus -c CONTIGS.db --search-term "L29" -n 7,7 -o test/ -O metagenome_RP_L29

INFO "Running anvi-export-locus in default-mode find gene with --search-term and --case-sensitive"
# should return nothing because all Ribosomal_L29 genes have capital L
anvi-export-locus -c CONTIGS.db --search-term "l29" --case-sensitive -n 7,7 -o test/ -O metagenome_RP_l29

INFO "Running anvi-export-locus in default-mode find gene with --search-term and --exact-match"
# KOFAM accession exact match
anvi-export-locus -c CONTIGS.db --search-term 'ko02044!!!ko02035' --exact-match -n 7,7 -o test/ -O metagenome_Type_II_secretion

INFO "Running anvi-export-locus in default-mode find gene with --search-term and --annotation-sources"
# Find Ribosomal_L29 gene with COG14_FUNCTION annotation
anvi-export-locus -c CONTIGS.db --search-term "L29" --annotation-sources 'COG14_FUNCTION' -n 7,7 -o test/ -O metagenome_COG14_FUNCTION_RP_L29 

INFO "Running anvi-export-locus in flank-mode with premature contig end"
anvi-export-locus -c CONTIGS.db --gene-caller-id 8616 -o test/ -n 7,30 -O metagenome_premature 

INFO "Running anvi-export-locus in flank-mode"
anvi-export-locus -c CONTIGS.db --gene-caller-id 68,78 -o test/ -O metagenome_68_78 --flank-mode

INFO "List hmm-sources in the contigs-db"
anvi-export-locus -c P_marinus_CCMP1375.db --use-hmm --hmm-sources Bacteria_71 --search-term "Exonuc_VII_L" -n 500,500 -o test/ -O P_marinus_CCMP1375 --list-hmm-sources

INFO "Running anvi-export-locus in default-mode using HMM to find a BIG locus in a genome"
anvi-export-locus -c P_marinus_CCMP1375.db --use-hmm  --search-term "Exonuc_VII_L" -n 500,500 -o test/ -O P_marinus_CCMP1375

INFO "Running anvi-export-locus in default-mode using HMM to find a BIG locus in a genome specificaly searching the hmm-source Bacteria_71"
anvi-export-locus -c P_marinus_CCMP1375.db --use-hmm --hmm-sources Bacteria_71 --search-term "Exonuc_VII_L" -n 500,500 -o test/ -O P_marinus_CCMP1375_Bacteria_71