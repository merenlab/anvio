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
                        --reaction-value-limits 5 35 \
                        --reaction-individual-value-limits 5 35 \
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
                        --reaction-individual-value-center \
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
                        --output-dir draw_txt_samples_kos_individual_colormap \
                        --reaction-gene-aggregation mean \
                        --reaction-sample-summary std \
                        --reaction-colormap viridis \
                        --reaction-individual-colormap RdYlGn 0.1 0.9 \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

if [ ! -s draw_txt_samples_kos_individual_colormap/colorbar_reactions.pdf ] \
    || [ ! -s draw_txt_samples_kos_individual_colormap/colorbar_reactions_samples.pdf ]
then
    echo "ERROR: two colormaps for one layer should have written a colorbar for each scale."
    exit 1
fi

# Clock times repeat every 24 hours. A period of 24 makes each scale of values span one whole day,
# from 0 to 24, so both ends of the scale are midnight. The scales take the cyclic 'clocktime'
# colormap, whose two ends have the same color. The KOs of a map element are averaged on a circle,
# so that times on either side of midnight average to a time near midnight.
INFO "Coloring the clock times of each sample with the cyclic 'clocktime' colormap"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_clocktime.reaction.txt \
                        --output-dir draw_txt_samples_kos_clocktime \
                        --reaction-value-period 24 \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

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

## COUNTS AND MEMBERSHIP COLORS
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
                        --reaction-membership-colors draw-sample-colors.txt \
                        --output-dir draw_txt_samples_membership_colors \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

INFO "Coloring group membership from a color given per group"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
                        --groups-txt draw-sample-group-information.txt \
                        --group-threshold 0.5 \
                        --reaction-membership-colors draw-group-colors.txt \
                        --output-dir draw_txt_groups_membership_colors \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

# Each group's ramp runs from white to that group's own color, so a group's maps are read on a
# scale of their own while staying recognizable as that group's.
INFO "Building each group's color scale as a ramp with '--group-colormap tint'"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
                        --groups-txt draw-sample-group-information.txt \
                        --group-threshold 0.5 \
                        --reaction-membership-colors draw-group-colors.txt \
                        --group-colormap tint \
                        --output-dir draw_txt_groups_tint_ramps \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

INFO "Giving the per-group color scale a colormap of its own"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
                        --output-dir draw_txt_groups_kos_individual_colormap \
                        --groups-txt draw-sample-group-information.txt \
                        --reaction-sample-summary mean \
                        --reaction-group-summary std \
                        --reaction-colormap cividis \
                        --reaction-individual-colormap RdBu_r \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

if [ ! -s draw_txt_groups_kos_individual_colormap/colorbar_reactions_groups.pdf ]
then
    echo "ERROR: a grouped run with two colormaps should have written a per-group colorbar."
    exit 1
fi

## NORMALIZING ELEMENT VALUES ACROSS SAMPLES
# Each sample's map shows how much more or less than usual an element is in that sample, rather than
# how much of it there is. The element values being sums of the values of the KOs an element stands
# for, this cannot be worked out per accession in the input file. The 'unified' map is unaffected:
# it summarizes the very samples the normalization compares, so it goes on showing coverage itself.
INFO "Rescaling each sample's element values against the element's mean across samples"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
                        --output-dir draw_txt_samples_kos_element_normalization \
                        --reaction-element-normalization relative_to_mean \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

if [ ! -s draw_txt_samples_kos_element_normalization/colorbar_reactions_samples.pdf ]
then
    echo "ERROR: a normalized run should have written a colorbar for the per-sample scale."
    exit 1
fi

# A label of the normalization's own goes on the colorbar in place of the one anvi'o derives, and
# both layers can be rescaled at once, each on its own scale.
INFO "Rescaling both layers, with a colorbar label given for one of them"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
                        --compound-txt draw_compounds_samples.compound.txt \
                        --output-dir draw_txt_both_layers_element_normalization \
                        --reaction-element-normalization log2_ratio_to_mean "log2 fold change" \
                        --compound-element-normalization z_score \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --no-progress

if [ ! -s draw_txt_both_layers_element_normalization/colorbar_compounds_samples.pdf ]
then
    echo "ERROR: a normalized compound layer should have written a per-sample colorbar."
    exit 1
fi

# Under a groups-txt the rescaling happens across the GROUPS rather than across the samples, each
# group's value being what the sample summary pooled from its samples, so both summaries have to
# pool values: by presence the group maps would be colored by sample counts, with nothing to
# rescale, and the 'unified' map would need a presence threshold.
INFO "Rescaling each group's element values against the element's mean across groups"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
                        --groups-txt draw-sample-group-information.txt \
                        --output-dir draw_txt_groups_kos_element_normalization \
                        --reaction-sample-summary mean \
                        --reaction-group-summary mean \
                        --reaction-element-normalization z_score "coverage z-score across groups" \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

if [ ! -s draw_txt_groups_kos_element_normalization/colorbar_reactions_groups.pdf ]
then
    echo "ERROR: a normalized grouped run should have written a colorbar for the per-group scale."
    exit 1
fi

## CATEGORICAL VALUES
# A value column can hold names, such as 'up', 'down' and 'same', instead of numbers. Each name
# takes a color of its own. The map of each sample shows that sample's values. The 'unified' map
# shows presence across the samples, since a value need not mean the same thing in every sample.
INFO "Coloring each sample's reactions by categorical values, as grids and a map per sample"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_status.reaction.txt \
                        --reaction-value-type categorical \
                        --output-dir draw_txt_samples_kos_categorical \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

if [ ! -s draw_txt_samples_kos_categorical/colorbar_reactions_samples.pdf ]
then
    echo "ERROR: a run of categorical values should have written a legend for the sample maps."
    exit 1
fi

# A colors file gives each value a color of its own choosing. The order of its rows also ranks the
# values, so 'order' gives a map element whose KOs have several values the first of them.
INFO "Coloring categorical values with colors of your choosing, which also rank the values"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_status.reaction.txt \
                        --reaction-value-type categorical \
                        --reaction-value-colors draw-status-colors.txt \
                        --reaction-accession-aggregation order \
                        --output-dir draw_txt_samples_kos_categorical_colors \
                        --pathway-numbers $pathway_numbers \
                        --draw-grid \
                        --no-progress

# Where a value means the same thing in every sample, a rule chooses one value across the samples
# for the 'unified' map. Here a tie between samples is broken at random, from a given seed.
INFO "Choosing one categorical value across the samples, breaking ties at random"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_status.reaction.txt \
                        --reaction-value-type categorical \
                        --reaction-sample-summary plurality_random \
                        --seed 1 \
                        --output-dir draw_txt_samples_kos_categorical_unified \
                        --pathway-numbers $pathway_numbers \
                        --no-progress

# With groups, a rule among the samples of each group colors that group's map. A rule across the
# groups colors the 'unified' map. Each group gives one vote, whatever its size.
INFO "Choosing one categorical value in each group of samples, and one across the groups"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_status.reaction.txt \
                        --reaction-value-type categorical \
                        --groups-txt draw-sample-group-information.txt \
                        --reaction-sample-summary plurality_order \
                        --reaction-group-summary plurality_order \
                        --output-dir draw_txt_groups_kos_categorical \
                        --pathway-numbers $pathway_numbers \
                        --draw-grid \
                        --no-progress

## SLICING REACTION ELEMENTS
# A map element often stands for several KOs. With '--reaction-slice-by accession', each KO with a
# value on a map takes a slice of the element in the color of its own value. A box is cut into
# columns, and a line into repeating dashes. Only maps that show values are sliced. Here these are
# the maps of the samples. The 'unified' map shows presence, so its elements are drawn whole.
INFO "Slicing each element by KO on the maps of individual samples, as grids and a map per sample"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
                        --reaction-slice-by accession \
                        --output-dir draw_txt_samples_kos_slices \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

# Slices have no labels, so a table beside the colorbars lists the slices of every map.
if [ ! -s draw_txt_samples_kos_slices/slice_key_reactions.txt ]
then
    echo "ERROR: a run that slices reaction elements should have written a key table of slices."
    exit 1
fi

# A column of weights sizes the slices. Here the color of a slice shows the expression of its KO,
# and its width shows the coverage of the KO's genes. The map of each sample takes its shares from
# that sample. Under 'share_mean', the default, the map of a group averages the shares of its
# samples, and the 'unified' map averages the shares of the groups.
INFO "Sizing slices by a column of weights, with shares averaged within groups and across them"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_weights.reaction.txt \
                        --groups-txt draw-sample-group-information.txt \
                        --reaction-sample-summary mean \
                        --reaction-group-summary mean \
                        --reaction-slice-by accession \
                        --reaction-slice-weights coverage \
                        --reaction-slice-share-summary share_mean \
                        --output-dir draw_txt_groups_kos_slices_share_mean \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

# 'pooled_sum' adds up the weights of each KO across the samples and divides each sum by the total
# of the element. A sample with more weight therefore counts for more.
INFO "Pooling the weights of slices across samples by their sum on the 'unified' map"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_weights.reaction.txt \
                        --reaction-sample-summary mean \
                        --reaction-slice-by accession \
                        --reaction-slice-weights coverage \
                        --reaction-slice-share-summary pooled_sum \
                        --output-dir draw_txt_samples_kos_slices_pooled_sum \
                        --pathway-numbers $pathway_numbers \
                        --draw-grid \
                        --no-progress

# Categorical values are sliced too. Each slice takes the color of its KO's value. On the 'unified'
# map, a rule across the samples chooses the value of each slice.
INFO "Slicing elements by the categorical values of their KOs"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_status.reaction.txt \
                        --reaction-value-type categorical \
                        --reaction-value-colors draw-status-colors.txt \
                        --reaction-sample-summary plurality_order \
                        --reaction-slice-by accession \
                        --output-dir draw_txt_samples_kos_categorical_slices \
                        --pathway-numbers $pathway_numbers \
                        --draw-grid \
                        --no-progress

# A file of KEGG reaction IDs is sliced by reaction rather than by KO. R01061 and R01063 share a map
# element on both maps.
INFO "Slicing elements by KEGG reaction"
anvi-draw-kegg-pathways --reaction-txt draw_reactions_shared_element.reaction.txt \
                        --reaction-slice-by accession \
                        --output-dir draw_txt_reactions_slices \
                        --pathway-numbers $pathway_numbers \
                        --no-progress

# With '--reaction-slice-by gene', each gene with a value on a map takes a slice of each element
# whose KOs it has rows for. A gene with rows for two KOs of one element is one slice. One gene ID
# names the same gene in every sample, so the 'unified' map pools the values of each gene across the
# samples. Here the color of a slice shows the expression of its gene, and its width shows the
# coverage of the gene.
INFO "Slicing each element by gene, with the values of each gene pooled across the samples"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_genes.reaction.txt \
                        --reaction-sample-summary mean \
                        --reaction-slice-by gene \
                        --reaction-slice-weights coverage \
                        --output-dir draw_txt_samples_genes_slices \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

# With groups, the map of each group pools the genes of its samples. Under 'pooled_sum', the map of
# a group adds up the weights of each gene over the samples of the group, and the 'unified' map adds
# them up over every sample.
INFO "Slicing each element by gene on the maps of groups, with weights pooled by their sum"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_genes.reaction.txt \
                        --groups-txt draw-sample-group-information.txt \
                        --reaction-sample-summary mean \
                        --reaction-group-summary mean \
                        --reaction-slice-by gene \
                        --reaction-slice-weights coverage \
                        --reaction-slice-share-summary pooled_sum \
                        --output-dir draw_txt_groups_genes_slices_pooled_sum \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

## BORDER VALUES
# A second column can color the borders of map elements, while the value column still colors their
# fills. '--reaction-border-column' names it. A border is a ring inside the outline of each box on
# standard maps, and of each compound circle on standard and overview maps. Here the fills show the
# coverage of each sample, and the borders show its status. The 'unified' map chooses one status
# across the samples, so it shows a border too. The borders take keys of their own.
INFO "Coloring borders by a column of categorical values beside a value column"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage_status.reaction.txt \
                        --reaction-border-column status \
                        --reaction-border-value-type categorical \
                        --reaction-border-value-colors draw-status-colors.txt \
                        --reaction-border-sample-summary plurality_order \
                        --output-dir draw_txt_samples_kos_border_status \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

if [ ! -s draw_txt_samples_kos_border_status/colorbar_reactions_border.pdf ] \
    || [ ! -s draw_txt_samples_kos_border_status/colorbar_reactions_border_samples.pdf ]
then
    echo "ERROR: a run with borders should have written a key for each scale of the borders."
    exit 1
fi

# A file whose only column of values is its border column colors its fills by presence. With
# groups, the border of each group's map takes the status chosen among the group's samples, and the
# 'unified' map takes the status chosen among the groups.
INFO "Coloring borders across groups of samples, with the fills showing presence"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_status.reaction.txt \
                        --reaction-border-column status \
                        --reaction-border-value-type categorical \
                        --groups-txt draw-sample-group-information.txt \
                        --group-threshold 0.5 \
                        --reaction-border-sample-summary plurality_order \
                        --reaction-border-group-summary plurality_order \
                        --output-dir draw_txt_groups_kos_border_status \
                        --pathway-numbers $pathway_numbers \
                        --draw-grid \
                        --no-progress

if [ ! -s draw_txt_groups_kos_border_status/colorbar_reactions_border.pdf ] \
    || [ ! -s draw_txt_groups_kos_border_status/colorbar_reactions_border_groups.pdf ]
then
    echo "ERROR: a grouped run with borders should have written a key for each border scale."
    exit 1
fi

# Numeric border values take 'viridis' by default, so their colors differ from those of the fills.
# Each layer has a border of its own. The reaction file here has no other column of values, so its
# fills show presence and get no colorbar.
INFO "Coloring the borders of both layers by numbers"
anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
                        --reaction-border-column coverage \
                        --compound-txt draw_compounds_samples.compound.txt \
                        --compound-border-column concentration \
                        --compound-border-sample-summary mean \
                        --output-dir draw_txt_both_layers_border \
                        --pathway-numbers $pathway_numbers \
                        --no-progress

if [ ! -s draw_txt_both_layers_border/colorbar_reactions_border.pdf ] \
    || [ ! -s draw_txt_both_layers_border/colorbar_compounds_border.pdf ]
then
    echo "ERROR: a run with borders on both layers should have written a key for each border."
    exit 1
fi
if [ -e draw_txt_both_layers_border/colorbar_reactions.pdf ]
then
    echo "ERROR: fills that show presence without samples should have no colorbar."
    exit 1
fi

# '--compound-color none' leaves the fills of a compound file without a value column uncolored, so
# that its circles show only their rings. On the overview map 01200, the circles keep the fills that
# they take from the reactions that touch them. A circle there shows the coverage of its reactions
# in its fill and the mean concentration of its compound in its ring. The compound layer writes a
# key of its rings and none of its fills.
INFO "Drawing compound rings alone, over the fills that compounds take from reactions"
anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
                        --compound-txt draw_compounds_samples.compound.txt \
                        --compound-border-column concentration \
                        --compound-border-sample-summary mean \
                        --compound-color none \
                        --output-dir draw_txt_both_layers_compound_rings_alone \
                        --pathway-numbers $pathway_numbers \
                        --no-progress

if [ ! -s draw_txt_both_layers_compound_rings_alone/colorbar_compounds_border.pdf ]
then
    echo "ERROR: a compound layer of rings alone should have written a key of its rings."
    exit 1
fi
if [ -e draw_txt_both_layers_compound_rings_alone/colorbar_compounds.pdf ]
then
    echo "ERROR: a compound layer of rings alone should have written no key of its fills."
    exit 1
fi

# '--reaction-border-slice-by' divides the ring of each box among its KOs or genes. Here the fill
# is sliced by KO as well, so the rings are aligned with it. Each column of the fill takes the
# status of its own KO, as bands above and below the column. The coverage of each KO sizes its
# column, and its expression colors the column. A second key table lists the slices of the rings.
INFO "Slicing borders by KO, aligned with fills that are sliced by KO and sized by weights"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_weights_status.reaction.txt \
                        --reaction-sample-summary mean \
                        --reaction-slice-by accession \
                        --reaction-slice-weights coverage \
                        --reaction-border-column status \
                        --reaction-border-value-type categorical \
                        --reaction-border-value-colors draw-status-colors.txt \
                        --reaction-border-sample-summary plurality_order \
                        --reaction-border-slice-by accession \
                        --output-dir draw_txt_samples_kos_border_slices_aligned \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

if [ ! -s draw_txt_samples_kos_border_slices_aligned/slice_key_reactions_border.txt ]
then
    echo "ERROR: a run that slices borders should have written a key table of the slices of rings."
    exit 1
fi

# Without '--reaction-slice-by', each ring is divided on its own, clockwise from the top left corner
# of its box, while the fill stays whole. The weights size the slices of the rings.
INFO "Slicing borders by KO on their own, with weights, around fills drawn whole"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_weights_status.reaction.txt \
                        --reaction-sample-summary mean \
                        --reaction-slice-weights coverage \
                        --reaction-border-column status \
                        --reaction-border-value-type categorical \
                        --reaction-border-value-colors draw-status-colors.txt \
                        --reaction-border-sample-summary plurality_order \
                        --reaction-border-slice-by accession \
                        --output-dir draw_txt_samples_kos_border_slices_independent \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

if [ ! -s draw_txt_samples_kos_border_slices_independent/slice_key_reactions_border.txt ]
then
    echo "ERROR: a run that slices borders on their own should have written their key table."
    exit 1
fi
if [ -e draw_txt_samples_kos_border_slices_independent/slice_key_reactions.txt ]
then
    echo "ERROR: fills drawn whole should have no key table of slices."
    exit 1
fi

## PRESENCE COLORING
# A bar column sets the length of a bar in each reaction box on standard maps, while the value
# column still sets the color of its fill. A bar is drawn in the color of the fill over a gray track,
# and a thin dark line marks its end. Here the fills show the expression of each sample, and the
# bars show its coverage. The 'unified' map shows bars too, since a bar sample summary pools the
# coverage of the samples. The bars take keys of their own.
INFO "Drawing bars of a second column in boxes whose fills are colored by values"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_genes.reaction.txt \
                        --reaction-bar-column coverage \
                        --reaction-bar-sample-summary mean \
                        --output-dir draw_txt_samples_kos_bars \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

if [ ! -s draw_txt_samples_kos_bars/bar_key_reactions.pdf ] \
    || [ ! -s draw_txt_samples_kos_bars/bar_key_reactions_samples.pdf ]
then
    echo "ERROR: a run with bars should have written a key for each bar scale."
    exit 1
fi

# A file whose only other column of values is its bar column colors its fills by presence, and its
# bars take the colors of presence. A ring of a border column lies around the bar.
INFO "Drawing bars inside rings, in boxes whose fills show presence"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage_status.reaction.txt \
                        --reaction-bar-column coverage \
                        --reaction-bar-sample-summary mean \
                        --reaction-border-column status \
                        --reaction-border-value-type categorical \
                        --reaction-border-value-colors draw-status-colors.txt \
                        --reaction-border-sample-summary plurality_order \
                        --output-dir draw_txt_samples_kos_bars_rings \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --no-progress

# The strip style draws the fill of each box whole and the bar as a dark bar on a white strip along
# the bottom of the box. Both values then stay readable: expression by color, and coverage by
# length. A ring of a third column lies around both.
INFO "Drawing bars as strips under whole fills, inside rings"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_weights_status.reaction.txt \
                        --reaction-sample-summary mean \
                        --reaction-bar-column coverage \
                        --reaction-bar-style strip \
                        --reaction-bar-sample-summary mean \
                        --reaction-border-column status \
                        --reaction-border-value-type categorical \
                        --reaction-border-value-colors draw-status-colors.txt \
                        --reaction-border-sample-summary plurality_order \
                        --output-dir draw_txt_samples_kos_strips_rings \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --no-progress

if [ ! -s draw_txt_samples_kos_strips_rings/bar_key_reactions.pdf ] \
    || [ ! -s draw_txt_samples_kos_strips_rings/bar_key_reactions_samples.pdf ]
then
    echo "ERROR: a run with strips should have written a key for each bar scale."
    exit 1
fi

# Bars can run up the height of their boxes in place of across them. The top of a vertical bar would
# run behind the label of its box, so these boxes drop their labels, and the key table of unlabeled
# boxes names them on each drawn map. The bar keys have their rulers beside their boxes.
INFO "Drawing vertical bars"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_genes.reaction.txt \
                        --reaction-bar-column coverage \
                        --reaction-bar-sample-summary mean \
                        --reaction-bar-orientation vertical \
                        --output-dir draw_txt_samples_kos_vertical_bars \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --no-progress

if [ "$(sed 1d draw_txt_samples_kos_vertical_bars/label_key_reactions.txt | cut -f2 | uniq \
    | tr '\n' ' ')" != "unified SAMPLE_1 SAMPLE_2 SAMPLE_3 " ] \
    || [ ! -s draw_txt_samples_kos_vertical_bars/bar_key_reactions_samples.pdf ]
then
    echo "ERROR: a run with vertical bars should have listed its boxes and keyed its bars."
    exit 1
fi

# A bar column with values on both sides of 0 draws bars that grow both ways from a gray zero line.
# The boxes on such a scale drop their labels, which would hide the zero line, and a key table names
# them by their accessions.
INFO "Drawing bars that grow both ways from 0"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_log_ratio.reaction.txt \
                        --reaction-bar-column log_ratio \
                        --reaction-bar-sample-summary mean \
                        --output-dir draw_txt_samples_kos_signed_bars \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --no-progress

# Each drawn map lists its own boxes, the 'unified' map first.
if [ "$(sed 1d draw_txt_samples_kos_signed_bars/label_key_reactions.txt | cut -f2 | uniq \
    | tr '\n' ' ')" != "unified SAMPLE_1 SAMPLE_2 SAMPLE_3 " ]
then
    echo "ERROR: a run with bars on both sides of 0 should have listed the boxes without labels."
    exit 1
fi

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
                        --animate-files-by-map \
                        --no-progress

collated_dir=draw_txt_collate_by_map_categorized/by_map
brite_dir=$collated_dir/Metabolism/Carbohydrate_metabolism
if [ ! -s $brite_dir/ko00010_Glycolysis_Gluconeogenesis/SAMPLE_1.pdf ]
then
    echo "ERROR: gathering by map did not keep the BRITE subdirectories of the map."
    exit 1
fi
if [ ! -s $brite_dir/ko00010_Glycolysis_Gluconeogenesis.html ]
then
    echo "ERROR: the player of a map did not keep the BRITE subdirectories and the map's name."
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

INFO "Writing a player of each map gathered by map"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
                        --output-dir draw_txt_animate_by_map \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --collate-files-by-map \
                        --animate-files-by-map \
                        --no-progress

player=draw_txt_animate_by_map/by_map/ko00010.html
if [ ! -s $player ]
then
    echo "ERROR: the map gathered by map has no player."
    exit 1
fi
# The page describes each frame with a label.
if [ "$(grep -o '"label":' $player | wc -l)" -ne 3 ]
then
    echo "ERROR: the player of a map should hold a frame for each of the three samples."
    exit 1
fi
# The player holds everything it shows, so it can be opened and passed around as one file.
if grep -q '<script src=\|<link \|<img [^>]*src="[^d]' $player
then
    echo "ERROR: the player of a map loads a file of its own."
    exit 1
fi

INFO "Skipping the 'unified' map so that only the individual samples are drawn"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
                        --output-dir draw_txt_skip_unified \
                        --pathway-numbers $pathway_numbers \
                        --skip-unified-maps \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

if [ -e draw_txt_skip_unified/unified ]
then
    echo "ERROR: the 'unified' map was drawn for a run that skipped it."
    exit 1
fi
# The colorbar keying the summary sits at the top of the output directory rather than beside the
# map it keys, so not drawing that map has to take the bar with it.
if [ -e draw_txt_skip_unified/colorbar_reactions.pdf ]
then
    echo "ERROR: the colorbar keying the skipped 'unified' map was drawn."
    exit 1
fi
if [ ! -s draw_txt_skip_unified/individual/SAMPLE_1/ko00010.pdf ]
then
    echo "ERROR: skipping the 'unified' map should leave the maps of the individual samples."
    exit 1
fi
# The grid leads with the 'unified' map, so a run skipping it has to go on drawing grids from the
# individual maps alone rather than quietly drawing none.
if [ ! -s draw_txt_skip_unified/grid/ko00010.pdf ]
then
    echo "ERROR: skipping the 'unified' map should leave a grid of the individual samples."
    exit 1
fi

## LABELS
# By default, only the reaction elements of standard maps are labeled. Here the reaction elements
# and compound circles of standard and overview maps are all labeled.
INFO "Labeling reaction elements and compound circles on standard and overview maps"
anvi-draw-kegg-pathways --reaction-txt draw_kos_and_compounds.reaction.txt \
                        --compound-txt draw_compounds.compound.txt \
                        --reaction-labels standard overview \
                        --compound-labels standard overview \
                        --output-dir draw_txt_labels \
                        --pathway-numbers $pathway_numbers \
                        --no-progress

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
                        --reaction-membership-colors draw-contigs-db-colors.txt \
                        --output-dir contigs_dbs_kos_membership_colors \
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

## REACTION NETWORK JSON INPUT
INFO "Exporting the pangenomic reaction network to a metabolic model JSON file"
anvi-get-metabolic-model-file --pan-db TEST-PAN.db \
                              --genomes-storage TEST-GENOMES.db \
                              --output-file pan-network.json

INFO "Mapping KOs from a pangenomic reaction network JSON file"
anvi-draw-kegg-pathways --reaction-network-json pan-network.json \
                        --output-dir json_pan_kos \
                        --pathway-numbers $pathway_numbers \
                        --no-progress

# The genes of a pangenomic JSON are gene clusters. Each has one KO, under 'id'. Read the way the
# genes of a genomic JSON are read, its KOs would match no map, and nothing would be drawn.
if [ ! -s json_pan_kos/unified/ko00010.pdf ]
then
    echo "ERROR: the KOs of a pangenomic reaction network JSON should have drawn map 00010."
    exit 1
fi
