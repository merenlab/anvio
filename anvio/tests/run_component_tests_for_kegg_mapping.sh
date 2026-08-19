#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2 $3
#####################################

rn_python_script=$(readlink -f run_component_tests_for_reaction_network)

INFO "Checking for the required KEGG database set up by anvi'o in a default location"
$rn_python_script --check-default-kegg-database

INFO "Setting up the KEGG mapping analysis directory"
cp $files/mock_data_for_kegg_mapping/*.txt                $output_dir/
# These databases should already contain KO annotations.
cp $files/mock_data_for_pangenomics/*.db                  $output_dir/
cp $files/mock_data_for_pangenomics/external-genomes.txt  $output_dir/
cp $files/mock_data_for_pangenomics/group-information.txt $output_dir/pan-group-information.txt
cp $files/example_description.md                          $output_dir/
cd $output_dir/

# anvi-draw-kegg-pathways writes maps rather than tables, so there is nothing to show from its
# output. These are the inputs the tests below draw from.
SHOW_FILE draw_kos_samples.reaction.txt
SHOW_FILE draw_kos_samples_coverage.reaction.txt
SHOW_FILE draw-sample-group-information.txt

INFO "Migrating all databases"
anvi-migrate *db --migrate-quickly

INFO "Generating an anvi'o genomes storage"
anvi-gen-genomes-storage -e external-genomes.txt -o TEST-GENOMES.db --no-progress

INFO "Running the pangenome analysis with default parameters"
anvi-pan-genome -g TEST-GENOMES.db \
                -n TEST \
                --use-ncbi-blast \
                --description example_description.md \
                --no-progress \
                $thread_controller

# A permanent ModelSEED database can be installed with 'anvi-setup-modelseed-database'. Without
# one, a temporary database serves for this run, and the reaction network below is told where it is.
modelseed_controller=""
use_default_modelseed_db=$($rn_python_script --check-default-modelseed-database)
if [ "$use_default_modelseed_db" == "True" ]
then
    INFO "Using the ModelSEED Biochemistry database already set up in a default location"
else
    INFO "Setting up the ModelSEED Biochemistry database in a temporary directory"
    data_dir=$(mktemp -d)
    anvi-setup-modelseed-database --dir $data_dir
    modelseed_controller="--modelseed-dir $data_dir/MODELSEED"
fi

INFO "Generating a pangenomic reaction network"
anvi-reaction-network --pan-db TEST-PAN.db \
                      --genomes-storage TEST-GENOMES.db \
                      --no-progress \
                      $modelseed_controller

# 00010 is a standard map and 01200 an overview map, and every test below draws both. 01100, the
# global metabolism map, renders about ten times slower than those two together, so only one test
# per input engine adds it, to cover the gray a global map reserves for unhighlighted elements.
pathway_numbers="00010 01200"
global_pathway_numbers="$pathway_numbers 01100"

## TEXT FILE INPUT
INFO "Drawing KOs from a reaction-txt file by presence in a single color"
anvi-draw-kegg-pathways --reaction-txt draw_kos.reaction.txt \
                        --output-dir draw_txt_kos \
                        --pathway-numbers $global_pathway_numbers \
                        --no-progress

INFO "Drawing KOs across the samples of a 'sample' column, as grids and a map per sample"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
                        --output-dir draw_txt_samples_kos \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

INFO "Drawing KOs across groups of samples, as grids and a map per group"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
                        --groups-txt draw-sample-group-information.txt \
                        --group-threshold 0 \
                        --output-dir draw_txt_samples_kos_groups \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

INFO "Drawing KEGG reaction IDs from a reaction-txt file, colored by a value column"
anvi-draw-kegg-pathways --reaction-txt draw_reactions.reaction.txt \
                        --output-dir draw_txt_reactions \
                        --pathway-numbers $pathway_numbers \
                        --no-progress

INFO "Drawing compounds from a compound-txt file, colored by a value column"
anvi-draw-kegg-pathways --compound-txt draw_compounds.compound.txt \
                        --output-dir draw_txt_compounds \
                        --pathway-numbers $pathway_numbers \
                        --no-progress

INFO "Drawing the reaction and compound layers together, each on its own value scale"
anvi-draw-kegg-pathways --reaction-txt draw_kos_and_compounds.reaction.txt \
                        --compound-txt draw_compounds.compound.txt \
                        --output-dir draw_txt_reactions_and_compounds \
                        --pathway-numbers $pathway_numbers \
                        --no-progress

## VALUE COLORING AND COLOR SCALES
INFO "Coloring reactions by a value column"
anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
                        --output-dir draw_txt_kos_quantitative \
                        --pathway-numbers $pathway_numbers \
                        --no-progress

# A layer drawn across samples has two color scales, set separately: one for the 'unified' map
# that summarizes the samples, and one shared by the maps of the individual samples.
INFO "Coloring reactions by a value column across samples, as grids and a map per sample"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
                        --output-dir draw_txt_samples_kos_quantitative \
                        --reaction-gene-aggregation mean \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

INFO "Summarizing samples by the mean of their values on the 'unified' map"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
                        --output-dir draw_txt_samples_kos_quantitative_mean \
                        --reaction-gene-aggregation mean \
                        --reaction-sample-summary mean \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

# Truncating both ends of both scales keeps a few extreme values from stretching the colors over
# a range in which the rest cannot be told apart. Each colorbar marks its truncated ends.
INFO "Limiting both color scales of a value column, truncating each at both ends"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
                        --output-dir draw_txt_samples_kos_value_limits \
                        --reaction-gene-aggregation mean \
                        --reaction-sample-summary mean \
                        --reaction-value-limits 5 40 \
                        --reaction-category-value-limits 5 40 \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

# These values lean well to the positive side, so centering widens each scale downwards until it
# runs equally far either side of zero, and the colormap's middle color comes to mean zero.
INFO "Centering the two color scales of a signed value column on zero"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_log_ratio.reaction.txt \
                        --output-dir draw_txt_samples_kos_value_center \
                        --reaction-gene-aggregation mean \
                        --reaction-accession-aggregation mean \
                        --reaction-sample-summary mean \
                        --reaction-colormap RdBu_r \
                        --reaction-value-center \
                        --reaction-category-value-center \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

# Log-transformed abundances are negative throughout, so their midpoint is nowhere near zero and
# each layer's scale has to be centered on a number of its own.
INFO "Centering the reaction and compound scales on numbers of their own"
anvi-draw-kegg-pathways --reaction-txt draw_kos_negative.reaction.txt \
                        --compound-txt draw_compounds.compound.txt \
                        --output-dir draw_txt_both_layers_value_center \
                        --reaction-gene-aggregation mean \
                        --reaction-accession-aggregation mean \
                        --reaction-colormap coolwarm \
                        --compound-colormap PuOr \
                        --reaction-value-center -4 \
                        --compound-value-center 9 \
                        --pathway-numbers $pathway_numbers \
                        --no-progress

# Summarizing samples by their standard deviation makes the 'unified' map show how much they
# disagree rather than how large they are. A colormap of its own keeps that spread from being
# read as the same quantity, and each scale writes its own colorbar.
INFO "Giving the per-sample color scale a colormap of its own, with a colorbar for each"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
                        --output-dir draw_txt_samples_kos_category_colormap \
                        --reaction-gene-aggregation mean \
                        --reaction-sample-summary std \
                        --reaction-colormap viridis \
                        --reaction-category-colormap RdYlGn 0.1 0.9 \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

if [ ! -s draw_txt_samples_kos_category_colormap/colorbar_reactions.pdf ] \
    || [ ! -s draw_txt_samples_kos_category_colormap/colorbar_reactions_samples.pdf ]
then
    echo "ERROR: two colormaps for one layer should have written a colorbar for each scale."
    exit 1
fi

INFO "Drawing both layers across the same samples, each with its own aggregation"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
                        --compound-txt draw_compounds_samples.compound.txt \
                        --reaction-gene-aggregation mean \
                        --compound-accession-aggregation sum \
                        --reaction-sample-summary mean \
                        --compound-sample-summary max \
                        --output-dir draw_txt_both_layers_samples \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

## COUNTS AND CATEGORY COLORS
INFO "Summarizing samples by count rather than by membership"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
                        --output-dir draw_txt_samples_kos_count \
                        --reaction-sample-summary count \
                        --pathway-numbers $pathway_numbers \
                        --no-progress

INFO "Summarizing samples by count on a continuous color scale"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples4.reaction.txt \
                        --output-dir draw_txt_samples_kos_count_continuous \
                        --reaction-sample-summary count_continuous \
                        --pathway-numbers $pathway_numbers \
                        --no-progress

INFO "Coloring sample membership from a color given per sample"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
                        --reaction-category-colors draw-sample-colors.txt \
                        --output-dir draw_txt_samples_category_colors \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

INFO "Coloring group membership from a color given per group"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
                        --groups-txt draw-sample-group-information.txt \
                        --group-threshold 0.5 \
                        --reaction-category-colors draw-group-colors.txt \
                        --output-dir draw_txt_groups_category_colors \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

# Each group's ramp runs from white to that group's own color, so a group's maps are read on a
# scale of their own while staying recognizable as that group's.
INFO "Building each group's color scale as a ramp with '--group-colormap category'"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
                        --groups-txt draw-sample-group-information.txt \
                        --group-threshold 0.5 \
                        --reaction-category-colors draw-group-colors.txt \
                        --group-colormap category \
                        --output-dir draw_txt_groups_category_ramps \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

INFO "Giving the per-group color scale a colormap of its own"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
                        --output-dir draw_txt_groups_kos_category_colormap \
                        --groups-txt draw-sample-group-information.txt \
                        --reaction-sample-summary mean \
                        --reaction-group-summary std \
                        --reaction-colormap cividis \
                        --reaction-category-colormap RdBu_r \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

if [ ! -s draw_txt_groups_kos_category_colormap/colorbar_reactions_groups.pdf ]
then
    echo "ERROR: a grouped run with two colormaps should have written a per-group colorbar."
    exit 1
fi

## PRESENCE COLORING
INFO "Mixing presence coloring of reactions with value coloring of compounds"
anvi-draw-kegg-pathways --reaction-txt draw_kos.reaction.txt \
                        --compound-txt draw_compounds.compound.txt \
                        --output-dir draw_txt_mixed_presence_quant \
                        --pathway-numbers $pathway_numbers \
                        --no-progress

INFO "Mixing membership coloring across samples with value coloring of compounds"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
                        --compound-txt draw_compounds.compound.txt \
                        --output-dir draw_txt_mixed_samples \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

INFO "Coloring the reaction and compound presence layers with colors of their own"
anvi-draw-kegg-pathways --reaction-txt draw_kos.reaction.txt \
                        --compound-txt draw_compounds_presence.compound.txt \
                        --reaction-color "#FF8800" \
                        --compound-color "#0000FF" \
                        --output-dir draw_txt_custom_colors \
                        --pathway-numbers $pathway_numbers \
                        --no-progress

INFO "Highlighting a reaction presence layer in the reference map's original colors"
anvi-draw-kegg-pathways --reaction-txt draw_kos.reaction.txt \
                        --original-color \
                        --output-dir draw_txt_original_color \
                        --pathway-numbers $pathway_numbers \
                        --no-progress

## OUTPUT LAYOUT
INFO "Gathering the maps of individual samples into a directory per map"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
                        --output-dir draw_txt_collate_by_map \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --collate-files-by-map \
                        --no-progress

for sample in SAMPLE_1 SAMPLE_2 SAMPLE_3
do
    if [ ! -s draw_txt_collate_by_map/by_map/ko00010/$sample.pdf ]
    then
        echo "ERROR: the map gathered by map is missing the file of $sample."
        exit 1
    fi
done

# '-ef' asks whether two paths are the same file on disk. A map gathered by map should be a hard
# link to the map already drawn for the sample, not a second copy of it.
if ! [ draw_txt_collate_by_map/by_map/ko00010/SAMPLE_1.pdf \
       -ef draw_txt_collate_by_map/individual/SAMPLE_1/ko00010.pdf ]
then
    echo "ERROR: a map gathered by map should be a link to the map drawn for the sample."
    exit 1
fi

INFO "Gathering by map while keeping the BRITE subdirectories and the named files"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
                        --output-dir draw_txt_collate_by_map_categorized \
                        --pathway-numbers $pathway_numbers \
                        --name-files \
                        --categorize-files \
                        --draw-individual-files \
                        --draw-grid \
                        --collate-files-by-map \
                        --no-progress

collated_dir=draw_txt_collate_by_map_categorized/by_map
brite_dir=$collated_dir/Metabolism/Carbohydrate_metabolism
if [ ! -s $brite_dir/ko00010_Glycolysis_Gluconeogenesis/SAMPLE_1.pdf ]
then
    echo "ERROR: gathering by map did not keep the BRITE subdirectories of the map."
    exit 1
fi
if [ -e $collated_dir/all_maps ]
then
    echo "ERROR: the flat directory of links was gathered as if it held maps of its own."
    exit 1
fi

# The flat directory holds links to the categorized maps rather than copies of them.
categorized_dir=draw_txt_collate_by_map_categorized/unified
brite_dir=$categorized_dir/Metabolism/Carbohydrate_metabolism
if ! [ $categorized_dir/all_maps/ko00010_Glycolysis_Gluconeogenesis.pdf \
       -ef $brite_dir/ko00010_Glycolysis_Gluconeogenesis.pdf ]
then
    echo "ERROR: the flat directory should hold a link to each categorized map, not a copy."
    exit 1
fi

INFO "Gathering by map only the maps of the groups asked for individually"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
                        --groups-txt draw-sample-group-information.txt \
                        --group-threshold 0 \
                        --output-dir draw_txt_collate_by_map_groups \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files G1 \
                        --draw-grid \
                        --collate-files-by-map \
                        --no-progress

collated_dir=draw_txt_collate_by_map_groups/by_map
if [ ! -s $collated_dir/ko00010/G1.pdf ]
then
    echo "ERROR: the group asked for individually was not gathered by map."
    exit 1
fi
if [ -e $collated_dir/ko00010/G2.pdf ]
then
    echo "ERROR: a group drawn as a grid panel alone was gathered by map."
    exit 1
fi
if [ -e $collated_dir/colorbar ]
then
    echo "ERROR: a group's colorbar was gathered by map as if it were one of its maps."
    exit 1
fi

## CONTIGS DATABASE INPUT
INFO "Mapping KOs from a genomic contigs database"
anvi-draw-kegg-pathways --contigs-dbs E_faecalis_6240.db \
                        --output-dir contigs_db_kos \
                        --name-files \
                        --categorize-files \
                        --pathway-numbers $pathway_numbers 05310 \
                        --draw-bare-maps \
                        --no-progress

INFO "Mapping KOs from multiple contigs databases, displaying database counts"
anvi-draw-kegg-pathways --external-genomes external-genomes.txt \
                        --output-dir contigs_dbs_kos_count \
                        --categorize-files \
                        --pathway-numbers $pathway_numbers \
                        --draw-grid \
                        --no-progress

INFO "Mapping KOs from multiple contigs databases, displaying database membership"
anvi-draw-kegg-pathways --contigs-dbs E_faecalis_6240.db E_faecalis_6255.db E_faecalis_6512.db \
                        --output-dir contigs_dbs_kos_membership \
                        --pathway-numbers $global_pathway_numbers \
                        --no-progress

# Each database's own map takes that database's color and matches its band on the 'unified' map.
INFO "Coloring contigs database membership from a color given per database"
anvi-draw-kegg-pathways --contigs-dbs E_faecalis_6240.db E_faecalis_6255.db E_faecalis_6512.db \
                        --reaction-category-colors draw-contigs-db-colors.txt \
                        --output-dir contigs_dbs_kos_category_colors \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

## PANGENOMIC DATABASE INPUT
INFO "Mapping KOs from a pangenomic database, displaying genome counts"
anvi-draw-kegg-pathways --pan-db TEST-PAN.db \
                        --genomes-storage TEST-GENOMES.db \
                        --output-dir pan_db_kos_genome_count_emphasize_shared \
                        --pathway-numbers $global_pathway_numbers \
                        --no-progress

INFO "Mapping KOs from a pangenomic database, displaying presence and absence"
anvi-draw-kegg-pathways --pan-db TEST-PAN.db \
                        --genomes-storage TEST-GENOMES.db \
                        --output-dir pan_db_kos_presence_absence \
                        --pathway-numbers $pathway_numbers \
                        --reaction-color \
                        --no-progress

INFO "Mapping KOs from a pangenomic database, displaying counts for a subset of groups"
anvi-draw-kegg-pathways --pan-db TEST-PAN.db \
                        --genomes-storage TEST-GENOMES.db \
                        --groups-txt pan-group-information.txt \
                        --group-threshold 0 \
                        --output-dir pan_db_kos_group_count_emphasize_shared \
                        --categorize-files \
                        --draw-individual-files G1 \
                        --draw-grid G1 \
                        --pathway-numbers $pathway_numbers \
                        --no-progress

INFO "Coloring pangenome reactions by genome count on a continuous scale"
anvi-draw-kegg-pathways --pan-db TEST-PAN.db \
                        --genomes-storage TEST-GENOMES.db \
                        --output-dir pan_db_kos_count_continuous \
                        --presence-colormap-scheme by_count_continuous \
                        --pathway-numbers $pathway_numbers \
                        --no-progress
