
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

INFO "Running anvi-export-locus in flank-mode with premature contig end"
anvi-export-locus -c CONTIGS.db --gene-caller-id 8616 -o test/ -n 7,30 -O metagenome_premature 

INFO "Running anvi-export-locus in flank-mode"
anvi-export-locus -c CONTIGS.db --gene-caller-id 68,78 -o test/ -O metagenome_68_78 --flank-mode