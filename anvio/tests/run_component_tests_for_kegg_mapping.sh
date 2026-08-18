#!/bin/bash
source 00.sh

# Several draw commands below are piped into 'tee' so that a warning can be matched against the log
# afterwards. Without this, such a pipeline exits with the status of 'tee', which always succeeds,
# so the 'set -e' of 00.sh would not see a failed draw and the test would pass on a crash.
set -o pipefail

SETUP_WITH_OUTPUT_DIR $1 $2 $3

rn_python_script=$(readlink -f run_component_tests_for_reaction_network)

INFO "Checking for the required KEGG database set up by anvi'o in a default location"
${rn_python_script} --check-default-kegg-database

INFO "Setting up the KEGG mapping analysis directory"
mkdir -p ${output_dir}/
# These databases should already contain KO annotations.
cp ${files}/mock_data_for_pangenomics/*.db ${output_dir}/
cp ${files}/mock_data_for_pangenomics/external-genomes.txt ${output_dir}/
cp ${files}/example_description.md ${output_dir}/
cp ${files}/mock_data_for_pangenomics/group-information.txt ${output_dir}/pan-group-information.txt
cp ${files}/data/input_files/minimal_enzymes_input.txt ${output_dir}/
awk -F'\t' 'NR==1 {print $0} NR>1 {$1=$1".db"; print $0}' OFS='\t' \
${output_dir}/pan-group-information.txt > ${output_dir}/contigs-db-group-information.txt

cd ${output_dir}/
# Build per-layer draw-kegg-pathways text files from the enzymes input's KOfam rows (gene ID and KO
# accession), plus a groups file assigning samples to groups. Variants add a 'sample' column and/or
# a numeric 'coverage' value column to test comparison across samples/groups and quantitative
# coloring. Reaction-layer files ('.reaction.txt') hold KO or KEGG reaction accessions;
# compound-layer files ('.compound.txt') hold KEGG compound accessions.
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "accession","gene_id"} NR>1{print $2,$1}' \
minimal_enzymes_input.txt > draw_kos.reaction.txt
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "accession","gene_id","sample"} NR>1{s=((NR-2)%3)+1; print $2,$1,"SAMPLE_"s}' \
minimal_enzymes_input.txt > draw_kos_samples.reaction.txt
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "accession","gene_id","sample"} NR>1{s=((NR-2)%4)+1; print $2,$1,"SAMPLE_"s}' \
minimal_enzymes_input.txt > draw_kos_samples4.reaction.txt
printf 'sample\tgroup\nSAMPLE_1\tG1\nSAMPLE_2\tG2\nSAMPLE_3\tG2\n' > draw-sample-group-information.txt
# More samples than a colormap has distinguishable colors (a colormap holds 256), so that coloring
# elements by sample count cannot give each count its own color. Presence there is drawn on a
# continuous scale rather than in discrete bands.
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "accession","gene_id","sample"} NR>1{s=((NR-2)%300)+1; print $2,$1,sprintf("SAMPLE_%03d",s)}' \
minimal_enzymes_input.txt > draw_kos_samples300.reaction.txt
# Enough samples that a discrete colorbar cannot label one band per count in type large enough to
# read, but not so many that the colormap runs out of distinguishable colors: the count is drawn on
# a continuous scale for the second reason rather than the first.
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "accession","gene_id","sample"} NR>1{s=((NR-2)%60)+1; print $2,$1,sprintf("SAMPLE_%02d",s)}' \
minimal_enzymes_input.txt > draw_kos_samples60.reaction.txt
# Groups over those samples, nearly all of them in one group, so that the group's own scale runs
# over more counts than its colorbar can label.
{
    printf 'sample\tgroup\n'
    for sample in $(seq -f 'SAMPLE_%02g' 1 59)
    do
        printf '%s\tWIDE\n' "${sample}"
    done
    printf 'SAMPLE_60\tTINY\n'
} > draw-sample-group-information-60.txt
# The same 300 samples, but with each accession in only a handful of them: a scale running to 300
# would squeeze every element into the bottom of the colormap, which is what '--count-scale-max' is
# for.
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "accession","gene_id","sample"} NR>1{for (i = 1; i <= ((NR-2)%5)+1; i++) print $2,$1,sprintf("SAMPLE_%03d",(((NR-2)*7+i)%300)+1)}' \
minimal_enzymes_input.txt > draw_kos_sparse300.reaction.txt
# A compound layer over the same samples, so that both layers fall back to a continuous count scale
# in one run and each says which layer it is talking about.
{
    printf 'accession\tsample\n'
    for accession in C00031 C00668 C00074 C00267
    do
        for sample in $(seq -f 'SAMPLE_%03g' 1 300)
        do
            printf '%s\t%s\n' "${accession}" "${sample}"
        done
    done
} > draw_compounds_samples300.compound.txt
# Groups over those 300 samples, nearly all of them in one group, so that the group's own count scale
# runs over more counts than its colormap has distinguishable colors. That is what makes a group's
# maps fall back to a continuous count scale, as the 'unified' map's can.
{
    printf 'sample\tgroup\n'
    for sample in $(seq -f 'SAMPLE_%03g' 1 299)
    do
        printf '%s\tBIG\n' "${sample}"
    done
    printf 'SAMPLE_300\tTINY\n'
} > draw-sample-group-information-300.txt
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "accession","gene_id","coverage"} NR>1{print $2,$1,((NR-2)%37)+0.5}' \
minimal_enzymes_input.txt > draw_kos_coverage.reaction.txt
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "accession","gene_id","sample","coverage"} NR>1{s=((NR-2)%3)+1; print $2,$1,"SAMPLE_"s,((NR-2)%37)+0.5}' \
minimal_enzymes_input.txt > draw_kos_samples_coverage.reaction.txt

# A constant 'coverage' column. With 'max' at BOTH reduction levels, every per-reaction value reduces
# to the same constant, exercising the degenerate (vmin == vmax) branch of quantitative coloring and
# its single-value colorbar. (Both levels are named explicitly in the test below so that it stays
# degenerate whatever the defaults are; with 'sum', an accession with more genes, or an element with
# more accessions, would differ from the rest.)
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "accession","gene_id","coverage"} NR>1{print $2,$1,1}' \
minimal_enzymes_input.txt > draw_kos_constant_coverage.reaction.txt

# A value column of negative numbers, cycling from -2 down to -6, as log-transformed abundances are.
# Used to limit a color scale with negative limits, which is what such a column asks for, and which
# the command line has to carry through argparse without mistaking them for option names.
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "accession","gene_id","log_abundance"} NR>1{print $2,$1,-((NR-2)%9)/2-2}' \
minimal_enzymes_input.txt > draw_kos_negative.reaction.txt

# A signed value column with a 'sample' column, cycling from -2 up to +6, as a log ratio between two
# conditions is: it can come out either way, but here it leans well to the positive side. Both of
# its color scales therefore sit off zero, which is what centering them puts right, and both are
# lopsided enough that centering visibly moves them.
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "accession","gene_id","sample","log_ratio"} NR>1{s=((NR-2)%3)+1; print $2,$1,"SAMPLE_"s,(((NR-2)%17)-4)/2}' \
minimal_enzymes_input.txt > draw_kos_samples_log_ratio.reaction.txt

# Reaction-level (KEGG reaction IDs) and compound-level fixtures (KEGG IDs present on the tested maps).
printf 'accession\tvalue\nR01068\t5\nR00710\t10\nR00229\t15\nR00755\t25\nR00431\t40\nR00014\t30\n' > draw_reactions.reaction.txt
# A reaction file and a compound file drawn together as two quantitative layers.
{ printf 'accession\tgene_id\tvalue\n'; awk -F'\t' 'NR>1{print $2"\t"$1"\t"(((NR-2)%37)+0.5)}' minimal_enzymes_input.txt; } > draw_kos_and_compounds.reaction.txt
printf 'accession\tvalue\nC00031\t2\nC00668\t5\nC00074\t9\nC00267\t14\n' > draw_compounds.compound.txt
# A compound file with a 'sample' column and a value column, for comparing metabolomes across samples.
{
    printf 'accession\tsample\tconcentration\n'
    printf 'C00031\tSAMPLE_1\t55\nC00031\tSAMPLE_2\t31\nC00668\tSAMPLE_1\t12\n'
    printf 'C00074\tSAMPLE_2\t30.5\nC00074\tSAMPLE_3\t8\nC00267\tSAMPLE_3\t14\n'
} > draw_compounds_samples.compound.txt
# A compound file with no value column (presence), for testing a custom compound presence color.
printf 'accession\nC00031\nC00668\nC00074\nC00267\n' > draw_compounds_presence.compound.txt

# Colors chosen per category rather than sampled from a colormap: one file per level, since a run's
# categories are its samples when they are ungrouped and its groups when they are not. The
# combination file additionally overrides the color that blending two samples' colors would derive.
printf 'category\tcolor\nSAMPLE_1\t#1f77b4\nSAMPLE_2\t#ff7f0e\nSAMPLE_3\t#2ca02c\n' > draw-sample-colors.txt
printf 'category\tcolor\nSAMPLE_1\t#1f77b4\nSAMPLE_2\t#ff7f0e\nSAMPLE_3\t#2ca02c\nSAMPLE_1, SAMPLE_2\t#7d5ba6\n' > draw-sample-colors-combo.txt
printf 'category\tcolor\nG1\t#1f77b4\nG2\t#d62728\n' > draw-group-colors.txt
printf 'category\tcolor\nSAMPLE_1\t#1f77b4\nSAMPLE_2\t#ff7f0e\n' > draw-sample-colors-incomplete.txt
# Two samples given the same color, so that two of the seven combinations of three samples come out
# the same and the membership scale cannot label its bands.
printf 'category\tcolor\nSAMPLE_1\t#1f77b4\nSAMPLE_2\t#ff7f0e\nSAMPLE_3\t#1f77b4\n' > draw-sample-colors-clashing.txt
# Contigs databases are named by their project names, which is what the maps compare them by. The
# groups of them are named G1 and G2, as the sample groups above are, so one file of group colors
# serves both inputs.
printf 'genome\tcolor\nE_faecalis_6240\t#1f77b4\nE_faecalis_6255\t#ff7f0e\nE_faecalis_6512\t#2ca02c\n' > draw-contigs-db-colors.txt

# Malformed files exercising the input-validation guards (each drawn command below is expected to fail).
printf 'accession\nK00844\nR00200\n' > draw_bad_mixed_kr.reaction.txt
printf 'accession\tgene_id\nC00031\tgene_1\n' > draw_bad_compound_gene_id.compound.txt
printf 'accession\tcol_a\tcol_b\nK00844\t1\t2\n' > draw_bad_two_values.reaction.txt
printf 'accession\tgene_id\tcoverage\nK00844\tgene_1\t10\nK00844\tgene_1\t30\n' > draw_bad_repeated_rows.reaction.txt
printf 'accession\tconcentration\nC00031\t10\nC00031\t30\n' > draw_bad_repeated_compounds.compound.txt

# Samples named after the output directory's own subdirectories and after a BRITE category: these
# are drawn under 'individual', so they cannot collide with anything anvi'o creates for itself.
printf 'accession\tsample\nK00844\tMetabolism\nK00845\tsymlink\nK02358\tgrid\nK00844\tunified\n' > draw_kos_awkward_sample_names.reaction.txt

# A sample whose name cannot be a directory name, alongside one that can. Only a sample drawn on its
# own maps becomes a subdirectory, so this file is fine until such maps are asked for it.
printf 'accession\tsample\nK00844\tSAMPLE_1\nK00845\tBAD/NAME\nK02358\tSAMPLE_1\n' > draw_kos_unusable_sample_names.reaction.txt

mkdir contigs_db_kos
mkdir contigs_dbs_kos_count
mkdir pan_db_kos_genome_count_emphasize_shared
mkdir pan_db_kos_genome_count_emphasize_unshared
mkdir pan_db_kos_presence_absence

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
                ${thread_controller}

use_default_modelseed_db=$(${rn_python_script} --check-default-modelseed-database)
if [ "${use_default_modelseed_db}" == "True" ]
then
    INFO "Using the ModelSEED Biochemistry database already set up by anvi'o in a default location"
else
    INFO "Setting up the ModelSEED Biochemistry database in a temporary directory (a permanent \
ModelSEED database can be installed in the default location with \
'anvi-setup-modelseed-database')"
    data_dir=$(mktemp -d)
    anvi-setup-modelseed-database --dir ${data_dir}
    modelseed_data_dir=${data_dir}/MODELSEED
fi

INFO "Generating a pangenomic reaction network"
args=()
args+=( "--pan-db" "TEST-PAN.db" )
args+=( "--genomes-storage" "TEST-GENOMES.db" )
if [ ${use_default_modelseed_db} == "False" ]
then
    args+=( "--modelseed-dir" ${modelseed_data_dir} )
fi
args+=( "--no-progress" )
anvi-reaction-network "${args[@]}"

# The map classes anvi'o treats differently: 00010 is a standard map and 01200 an overview map, both
# of which every test below draws. 01100, the global metabolism map, holds the whole of metabolism
# and takes about ten times as long to render as the other two together, so only the handful of tests
# that exercise something specific to a global map add it. What is specific to them is the color they
# reserve for their own unhighlighted elements -- gray, where overview maps use black and white and
# standard maps use white -- so the tests that add it are those staging colors through a distinct
# path: each coloring mode, each drawer, and each input engine, once.
pathway_numbers=( "00010" "01200" )
global_pathway_numbers=( "${pathway_numbers[@]}" "01100" )

INFO "Testing drawing KOs from a kegg-reaction-txt file by presence (single color)"
args=()
args+=( "--reaction-txt" "draw_kos.reaction.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_kos )
# Draws the global map too, for presence in a single color.
args+=( "--pathway-numbers" "${global_pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing drawing KOs across samples defined by a 'sample' column, drawing map grids and map \
files for each sample"
args=()
args+=( "--reaction-txt" "draw_kos_samples.reaction.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_samples_kos )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-individual-files" )
args+=( "--draw-grid" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing drawing KOs across groups of samples, where KOs are in any samples in the group, \
drawing map grids and map files for each group"
args=()
args+=( "--reaction-txt" "draw_kos_samples.reaction.txt" )
args+=( "--groups-txt" "draw-sample-group-information.txt" )
args+=( "--group-threshold" "0" )
args+=( "--output-dir" ${output_dir}/draw_txt_samples_kos_groups )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-individual-files" )
args+=( "--draw-grid" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing drawing reactions from a kegg-reaction-txt file of KEGG reaction IDs, colored by a \
value column"
args=()
args+=( "--reaction-txt" "draw_reactions.reaction.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_reactions )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing drawing compounds from a kegg-compound-txt file, colored by a value column"
args=()
args+=( "--compound-txt" "draw_compounds.compound.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_compounds )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing drawing the reaction and compound layers together from two files, each colored by \
its own value column on its own scale"
args=()
args+=( "--reaction-txt" "draw_kos_and_compounds.reaction.txt" )
args+=( "--compound-txt" "draw_compounds.compound.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_reactions_and_compounds )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing coloring reactions by a value column from a kegg-reaction-txt file"
args=()
args+=( "--reaction-txt" "draw_kos_coverage.reaction.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_kos_quantitative )
# Draws the global map too, for a value column, whose reactions also derive compound colors.
args+=( "--pathway-numbers" "${global_pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing coloring reactions by a value column across samples, where the per-sample maps show \
each sample's own values continuously and the default sample summary makes the 'unified' map show \
which samples contain each reaction, drawing map grids and map files for each sample"
args=()
args+=( "--reaction-txt" "draw_kos_samples_coverage.reaction.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_samples_kos_quantitative )
args+=( "--reaction-gene-aggregation" "mean" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-individual-files" )
args+=( "--draw-grid" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing summarizing samples by pooling their values, so that the 'unified' map is colored by \
the mean of the samples' values"
args=()
args+=( "--reaction-txt" "draw_kos_samples_coverage.reaction.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_samples_kos_quantitative_mean )
args+=( "--reaction-gene-aggregation" "mean" )
args+=( "--reaction-sample-summary" "mean" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-individual-files" )
args+=( "--draw-grid" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing limits on the two color scales of a value column, which are set separately: the \
'unified' map's scale, of means across samples, and the scale its per-sample maps share, of the \
samples' own values. Both are truncated at both ends here, so that a few extreme elements stop \
stretching the colors over a range in which the rest cannot be told apart, and both colorbars mark \
the truncated ends '<=' and '>='"
args=()
args+=( "--reaction-txt" "draw_kos_samples_coverage.reaction.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_samples_kos_value_limits )
args+=( "--reaction-gene-aggregation" "mean" )
args+=( "--reaction-sample-summary" "mean" )
args+=( "--reaction-value-limits" "5" "40" )
args+=( "--reaction-category-value-limits" "5" "40" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-individual-files" )
args+=( "--draw-grid" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing limits that no value crosses, which therefore truncate nothing: both scales still \
span their own values and neither colorbar is marked, which is what lets one set of limits be \
carried from dataset to dataset without rescaling the ones it was never needed for"
args=()
args+=( "--reaction-txt" "draw_kos_samples_coverage.reaction.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_samples_kos_value_limits_inert )
args+=( "--reaction-gene-aggregation" "mean" )
args+=( "--reaction-sample-summary" "mean" )
args+=( "--reaction-value-limits" "0" "200" )
args+=( "--reaction-category-value-limits" "0" "200" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-individual-files" )
args+=( "--draw-grid" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing negative limits, which is what a log-transformed value column asks for, with one end \
of the scale left open by 'none'. A layer with no 'sample' column has the one scale, so \
'--reaction-value-limits' bounds it and no second pair of limits is wanted"
args=()
args+=( "--reaction-txt" "draw_kos_negative.reaction.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_kos_negative_value_limits )
args+=( "--reaction-gene-aggregation" "mean" )
args+=( "--reaction-accession-aggregation" "mean" )
args+=( "--reaction-value-limits" "-5" "none" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing that limiting one of two scales that both color by value warns about the other, \
which goes on spanning whatever its own values reach"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --output-dir ${output_dir}/draw_txt_samples_kos_value_limits_one \
    --reaction-gene-aggregation mean --reaction-sample-summary mean \
    --reaction-value-limits 5 40 --pathway-numbers "${pathway_numbers[@]}" --draw-grid \
    --no-progress 2>&1 | tee draw_txt_value_limits_one.log
if ! tr '\n' ' ' < draw_txt_value_limits_one.log \
    | grep -q "still spans whatever its own values happen to reach"
then
    echo "ERROR: limiting one of two scales colored by value did not warn about the other."
    exit 1
fi

INFO "Testing centering the two color scales of a signed value column on zero, given as a bare \
flag, with the diverging colormap that centering is for. The values lean well to the positive \
side, so each scale is widened downwards until it runs the same distance either side of zero, and \
the neutral middle color of the colormap comes to stand for zero on both the 'unified' map and the \
maps of the individual samples"
args=()
args+=( "--reaction-txt" "draw_kos_samples_log_ratio.reaction.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_samples_kos_value_center )
args+=( "--reaction-gene-aggregation" "mean" )
args+=( "--reaction-accession-aggregation" "mean" )
args+=( "--reaction-sample-summary" "mean" )
args+=( "--reaction-colormap" "RdBu_r" )
args+=( "--reaction-value-center" )
args+=( "--reaction-category-value-center" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-individual-files" )
args+=( "--draw-grid" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing centering a color scale on a number of its own rather than on zero, which is what a \
quantity read against some other reference asks for, and centering the compound layer of a second \
file at the same time: the two layers are colored on independent scales, so each takes its own \
center"
args=()
args+=( "--reaction-txt" "draw_kos_negative.reaction.txt" )
args+=( "--compound-txt" "draw_compounds.compound.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_both_layers_value_center )
args+=( "--reaction-gene-aggregation" "mean" )
args+=( "--reaction-accession-aggregation" "mean" )
args+=( "--reaction-colormap" "coolwarm" )
args+=( "--compound-colormap" "PuOr" )
args+=( "--reaction-value-center" "-4" )
args+=( "--compound-value-center" "9" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing a centered scale whose longer side is trimmed back by a limit, which is how a \
centered scale is kept from wasting most of its colormap on a lopsided tail: the limits truncate \
what the values reach and the center then widens whichever side falls short, so the scale comes \
out centered on zero, over a good deal less than the values span, with its top marked '>='"
args=()
args+=( "--reaction-txt" "draw_kos_samples_log_ratio.reaction.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_samples_kos_value_center_limits )
args+=( "--reaction-gene-aggregation" "mean" )
args+=( "--reaction-accession-aggregation" "mean" )
args+=( "--reaction-sample-summary" "mean" )
args+=( "--reaction-colormap" "RdBu_r" )
args+=( "--reaction-value-center" )
args+=( "--reaction-category-value-center" )
args+=( "--reaction-value-limits" "none" "3" )
args+=( "--reaction-category-value-limits" "none" "3" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-grid" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing that centering one of two scales that both color by value warns about the other, \
which goes on sitting wherever its own values leave it while one colormap colors both"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_log_ratio.reaction.txt \
    --output-dir ${output_dir}/draw_txt_samples_kos_value_center_one \
    --reaction-gene-aggregation mean --reaction-sample-summary mean \
    --reaction-colormap RdBu_r --reaction-value-center \
    --pathway-numbers "${pathway_numbers[@]}" --draw-grid \
    --no-progress 2>&1 | tee draw_txt_value_center_one.log
if ! tr '\n' ' ' < draw_txt_value_center_one.log \
    | grep -q "still sits wherever its own values leave it"
then
    echo "ERROR: centering one of two scales colored by value did not warn about the other."
    exit 1
fi

INFO "Testing that centering a scale drawn from a sequential colormap warns, that colormap having \
no middle for a centered value to take, and that values lying entirely on one side of the center \
warn too, half of the colormap then going unused"
anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --output-dir ${output_dir}/draw_txt_kos_value_center_warnings \
    --reaction-gene-aggregation mean --reaction-accession-aggregation mean \
    --reaction-value-center --pathway-numbers "${pathway_numbers[@]}" \
    --no-progress 2>&1 | tee draw_txt_value_center_warnings.log
for phrase in "is not a diverging one" "all lie on one side of it"
do
    if ! tr '\n' ' ' < draw_txt_value_center_warnings.log | grep -q "${phrase}"
    then
        echo "ERROR: centering a sequential colormap over one-sided values did not warn (${phrase})."
        exit 1
    fi
done

INFO "Testing that centering a scale drawn from a diverging colormap trimmed off-center warns, the \
neutral color no longer being in the middle of what is drawn"
anvi-draw-kegg-pathways --reaction-txt draw_kos_negative.reaction.txt \
    --output-dir ${output_dir}/draw_txt_kos_value_center_trim \
    --reaction-gene-aggregation mean --reaction-accession-aggregation mean \
    --reaction-colormap RdBu_r 0.2 1.0 --reaction-value-center -4 \
    --pathway-numbers "${pathway_numbers[@]}" \
    --no-progress 2>&1 | tee draw_txt_value_center_trim.log
if ! tr '\n' ' ' < draw_txt_value_center_trim.log | grep -q "That trimming is not symmetric"
then
    echo "ERROR: centering a scale on an off-center trim of a diverging colormap did not warn."
    exit 1
fi

INFO "Testing a centered scale with nothing left to span, every value being the center itself: the \
colorbar is a single band, and it takes the middle color of the colormap rather than its top color, \
that middle being what the centering asked for"
args=()
args+=( "--reaction-txt" "draw_kos_constant_coverage.reaction.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_kos_value_center_degenerate )
args+=( "--reaction-gene-aggregation" "max" )
args+=( "--reaction-accession-aggregation" "max" )
args+=( "--reaction-colormap" "RdBu_r" )
args+=( "--reaction-value-center" "1" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing a colormap for the scale the per-sample maps share that differs from the one \
coloring the 'unified' map. The samples are summarized by their standard deviation, so the \
'unified' map shows how much they disagree rather than how large they are, and giving that \
spread its own colormap keeps the two from being read as the same quantity. The second \
colormap is trimmed, as the first can be, and each scale writes its own colorbar"
args=()
args+=( "--reaction-txt" "draw_kos_samples_coverage.reaction.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_samples_kos_category_colormap )
args+=( "--reaction-gene-aggregation" "mean" )
args+=( "--reaction-sample-summary" "std" )
args+=( "--reaction-colormap" "viridis" )
args+=( "--reaction-category-colormap" "RdYlGn" "0.1" "0.9" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-individual-files" )
args+=( "--draw-grid" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

if [ ! -s ${output_dir}/draw_txt_samples_kos_category_colormap/colorbar_reactions.pdf ] \
    || [ ! -s ${output_dir}/draw_txt_samples_kos_category_colormap/colorbar_reactions_samples.pdf ]
then
    echo "ERROR: two colormaps for one layer should have written a colorbar for each scale."
    exit 1
fi

INFO "Testing summarizing samples by count rather than by membership, overriding the default for \
three or fewer samples"
args=()
args+=( "--reaction-txt" "draw_kos_samples.reaction.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_samples_kos_count )
args+=( "--reaction-sample-summary" "count" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing summarizing samples by count on a continuous scale, which is asked for explicitly \
here even though four samples would fit in discrete bands: with the default sequential colormap \
the colors are the same as 'count' would assign, but the colorbar is a gradient from the lowest \
count to the highest"
args=()
args+=( "--reaction-txt" "draw_kos_samples4.reaction.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_samples_kos_count_continuous )
args+=( "--reaction-sample-summary" "count_continuous" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing that with a scale running to more samples than the colormap has distinguishable \
colors, the default presence summary falls back from discrete bands to a continuous count scale \
rather than failing"
args=()
args+=( "--reaction-txt" "draw_kos_samples300.reaction.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_samples300_kos_count_continuous )
args+=( "--count-scale-max" "total" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing that a scale with colors enough for every count but more bands than a colorbar can \
label falls back to the continuous scale too, and that asking for the discrete bands outright \
draws them with a warning that the labels will be too small to read"
args=()
args+=( "--reaction-txt" "draw_kos_samples60.reaction.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_samples60_bands_unreadable )
args+=( "--count-scale-max" "total" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}" 2>&1 | tee draw_txt_samples60_fallback.log
# Warnings are wrapped to the terminal width, so a phrase is matched against the log with its line
# breaks flattened to spaces rather than line by line.
if ! tr '\n' ' ' < draw_txt_samples60_fallback.log \
    | grep -q "more labels than it can set in type large enough to read"
then
    echo "ERROR: a count scale of unlabelable bands did not fall back to the continuous scale."
    exit 1
fi
# For a text file the scheme comes from the layer's summary rather than from
# '--presence-colormap-scheme', which such a run rejects.
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples60.reaction.txt \
    --output-dir ${output_dir}/draw_txt_samples60_bands_asked_for --count-scale-max total \
    --reaction-sample-summary count --pathway-numbers "${pathway_numbers[@]}" --no-progress 2>&1 \
    | tee draw_txt_samples60_bands.log
if ! tr '\n' ' ' < draw_txt_samples60_bands.log | grep -q "The bands are drawn as asked"
then
    echo "ERROR: insisting on more bands than a colorbar can label did not warn about the labels."
    exit 1
fi

INFO "Testing the same fallback with a reaction layer and a compound layer together: each layer \
decides for itself, so the run reports the fallback once per layer, saying which layer each time"
args=()
args+=( "--reaction-txt" "draw_kos_samples300.reaction.txt" )
args+=( "--compound-txt" "draw_compounds_samples300.compound.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_samples300_two_layers_count_continuous )
args+=( "--count-scale-max" "total" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}" 2>&1 | tee draw_txt_samples300_two_layers.log
for element_type in reaction compound
do
    if ! grep -q "The ${element_type} layer is colored by counts" \
        draw_txt_samples300_two_layers.log
    then
        echo "ERROR: the continuous count scale fallback did not name the ${element_type} layer."
        exit 1
    fi
done

INFO "Testing that a count scale stops at the highest count in the data by default, so that a \
sparse layer over 300 samples spreads its colors over the counts that occur instead of squeezing \
them into the bottom of the colormap"
args=()
args+=( "--reaction-txt" "draw_kos_sparse300.reaction.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_sparse300_count_scale_observed )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing running a count scale to every sample there is rather than to the highest count \
observed"
args=()
args+=( "--reaction-txt" "draw_kos_sparse300.reaction.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_sparse300_count_scale_total )
args+=( "--count-scale-max" "total" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing pinning a count scale to a given number, so that separate figures share a scale, \
with counts above it taking the top color"
args=()
args+=( "--reaction-txt" "draw_kos_sparse300.reaction.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_sparse300_count_scale_number )
args+=( "--count-scale-max" "3" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing a static presence color on a layer with a 'sample' column: it overrides coloring by \
sample on the 'unified' map, which shows presence in any sample in that one color, while the \
individual sample maps are unaffected"
args=()
args+=( "--reaction-txt" "draw_kos_samples.reaction.txt" )
args+=( "--reaction-color" "#FF8800" )
args+=( "--output-dir" ${output_dir}/draw_txt_samples_static_color )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-individual-files" )
args+=( "--draw-grid" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing coloring by membership from a color given per sample rather than sampled from a \
colormap: each sample alone takes its own color, each combination of them the blend of their \
colors, and each sample's own map that sample's color"
args=()
args+=( "--reaction-txt" "draw_kos_samples.reaction.txt" )
args+=( "--reaction-category-colors" "draw-sample-colors.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_samples_category_colors )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-individual-files" )
args+=( "--draw-grid" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing a row of the colors file that names a combination of samples, which sets that \
combination's color directly in place of the blend of its members' colors"
args=()
args+=( "--reaction-txt" "draw_kos_samples.reaction.txt" )
args+=( "--reaction-category-colors" "draw-sample-colors-combo.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_samples_category_colors_combo )
# Draws the global map too, for colors given per category, blends included.
args+=( "--pathway-numbers" "${global_pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}" 2>&1 | tee draw_txt_samples_combo.log
# A combination whose members do not parse back to the run's samples is reported and ignored rather
# than refused, so the run would exit 0 with the row silently unused and the combination drawn in
# the blend it was written to override. Its absence from the log is what says the row was taken.
if tr '\n' ' ' < draw_txt_samples_combo.log \
    | grep -q "name at least one thing that is not a sample of this run"
then
    echo "ERROR: the combination row was ignored rather than recoloring the combination."
    exit 1
fi

INFO "Testing that colors given per sample also color a compound layer, independently of the \
reaction layer, so that the two layers drawn on one map keep their own scales"
args=()
args+=( "--reaction-txt" "draw_kos_samples.reaction.txt" )
args+=( "--compound-txt" "draw_compounds_samples.compound.txt" )
args+=( "--reaction-category-colors" "draw-sample-colors.txt" )
args+=( "--compound-category-colors" "draw-sample-colors-combo.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_two_layers_category_colors )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing colors given per GROUP, which color the 'unified' map by group membership while each \
group's own maps keep the default group colormap: the identity and magnitude channels stay \
separate unless the group colormap is asked to follow the colors"
args=()
args+=( "--reaction-txt" "draw_kos_samples.reaction.txt" )
args+=( "--groups-txt" "draw-sample-group-information.txt" )
args+=( "--group-threshold" "0.5" )
args+=( "--reaction-category-colors" "draw-group-colors.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_groups_category_colors )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-individual-files" )
args+=( "--draw-grid" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing '--group-colormap category', which builds each group's own scale as a ramp running \
from a pale tint to that group's own color instead of sampling one colormap shared by every group"
args=()
args+=( "--reaction-txt" "draw_kos_samples.reaction.txt" )
args+=( "--groups-txt" "draw-sample-group-information.txt" )
args+=( "--group-threshold" "0.5" )
args+=( "--reaction-category-colors" "draw-group-colors.txt" )
args+=( "--group-colormap" "category" )
args+=( "--output-dir" ${output_dir}/draw_txt_groups_category_ramps )
# Draws the global map too, for a scale per group, built from that group's own color.
args+=( "--pathway-numbers" "${global_pathway_numbers[@]}" )
args+=( "--draw-individual-files" )
args+=( "--draw-grid" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing the same group ramps with an explicit span, which says how far from white each ramp \
starts and stops rather than which fraction of a colormap to use"
args=()
args+=( "--reaction-txt" "draw_kos_samples.reaction.txt" )
args+=( "--groups-txt" "draw-sample-group-information.txt" )
args+=( "--group-threshold" "0.5" )
args+=( "--reaction-category-colors" "draw-group-colors.txt" )
args+=( "--group-colormap" "category" "0.4" "0.95" )
args+=( "--group-reverse-overlay" )
args+=( "--output-dir" ${output_dir}/draw_txt_groups_category_ramps_span )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-individual-files" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing '--group-colormap-scheme by_count_continuous', which draws every group's count scale \
as a gradient rather than in discrete bands of one color per source, here for a ramp built from \
each group's own color, so that a group's colorbar spans that group's own hue"
args=()
args+=( "--reaction-txt" "draw_kos_samples.reaction.txt" )
args+=( "--groups-txt" "draw-sample-group-information.txt" )
args+=( "--group-threshold" "0.5" )
args+=( "--reaction-category-colors" "draw-group-colors.txt" )
args+=( "--group-colormap" "category" )
args+=( "--group-colormap-scheme" "by_count_continuous" )
args+=( "--output-dir" ${output_dir}/draw_txt_groups_scheme_continuous )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-individual-files" )
args+=( "--draw-grid" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing that a group whose count scale runs over more counts than its colorbar can label \
falls back to the continuous scale too, even though its colormap has a distinguishable color for \
every one of those counts"
args=()
args+=( "--reaction-txt" "draw_kos_samples60.reaction.txt" )
args+=( "--groups-txt" "draw-sample-group-information-60.txt" )
args+=( "--group-threshold" "0.5" )
args+=( "--count-scale-max" "total" )
args+=( "--output-dir" ${output_dir}/draw_txt_groups60_bands_unreadable )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-individual-files" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}" 2>&1 | tee draw_txt_groups60_fallback.log
flattened_groups60_log=$(tr '\n' ' ' < draw_txt_groups60_fallback.log)
if ! grep -q "more labels than it can set in type large enough to read" \
    <<< "${flattened_groups60_log}"
then
    echo "ERROR: a group scale of unlabelable bands did not fall back to the continuous scale."
    exit 1
fi
if grep -q "could supply only" <<< "${flattened_groups60_log}"
then
    echo "ERROR: the fallback blamed the colors when it was the labels that did not fit."
    exit 1
fi

INFO "Testing that asking for the discrete bands outright keeps them for a group whose scale runs \
over more counts than its colorbar can label, warning that the labels will be small — the colormap \
has a distinguishable color for every one of those counts, so nothing stops them being drawn"
args=()
args+=( "--reaction-txt" "draw_kos_samples60.reaction.txt" )
args+=( "--groups-txt" "draw-sample-group-information-60.txt" )
args+=( "--group-threshold" "0.5" )
args+=( "--count-scale-max" "total" )
args+=( "--group-colormap-scheme" "by_count" )
args+=( "--output-dir" ${output_dir}/draw_txt_groups60_scheme_bands )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-individual-files" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}" 2>&1 | tee draw_txt_groups60_scheme_bands.log
flattened_groups60_bands_log=$(tr '\n' ' ' < draw_txt_groups60_scheme_bands.log)
if ! grep -q "The bands are drawn as asked" <<< "${flattened_groups60_bands_log}"
then
    echo "ERROR: insisting on bands too many to label did not warn that they were drawn anyway."
    exit 1
fi
if grep -q "drawn on a CONTINUOUS color scale" <<< "${flattened_groups60_bands_log}"
then
    echo "ERROR: the discrete bands were asked for but the continuous scale was drawn."
    exit 1
fi
# The colormap has a color for every count here, so the OTHER limit must not be what is reported.
if grep -q "could supply only" <<< "${flattened_groups60_bands_log}"
then
    echo "ERROR: the warning blamed the colors when it was the labels that did not fit."
    exit 1
fi

INFO "Testing that a group whose count scale runs over more counts than its colormap can \
distinguish falls back to the continuous scale on its own, with a warning, and that asking for the \
discrete bands outright keeps them and warns that neighboring counts share a color"
args=()
args+=( "--reaction-txt" "draw_kos_samples300.reaction.txt" )
args+=( "--groups-txt" "draw-sample-group-information-300.txt" )
args+=( "--group-threshold" "0.5" )
args+=( "--count-scale-max" "total" )
args+=( "--output-dir" ${output_dir}/draw_txt_groups300_scheme_fallback )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-individual-files" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}" 2>&1 | tee draw_txt_groups300_scheme_fallback.log
# Warnings are wrapped to the terminal width, so a phrase is matched against the log with its line
# breaks flattened to spaces rather than line by line.
if ! tr '\n' ' ' < draw_txt_groups300_scheme_fallback.log \
    | grep -q "drawn on a CONTINUOUS color scale"
then
    echo "ERROR: a group scale short of colors did not fall back to the continuous scale."
    exit 1
fi
args=()
args+=( "--reaction-txt" "draw_kos_samples300.reaction.txt" )
args+=( "--groups-txt" "draw-sample-group-information-300.txt" )
args+=( "--group-threshold" "0.5" )
args+=( "--count-scale-max" "total" )
args+=( "--group-colormap-scheme" "by_count" )
args+=( "--output-dir" ${output_dir}/draw_txt_groups300_scheme_bands )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-individual-files" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}" 2>&1 | tee draw_txt_groups300_scheme_bands.log
flattened_bands_log=$(tr '\n' ' ' < draw_txt_groups300_scheme_bands.log)
if ! grep -q "Neighboring counts therefore share a color" <<< "${flattened_bands_log}"
then
    echo "ERROR: insisting on the discrete bands did not warn that counts share a color."
    exit 1
fi
if grep -q "drawn on a CONTINUOUS color scale" <<< "${flattened_bands_log}"
then
    echo "ERROR: the discrete bands were asked for but the continuous scale was drawn."
    exit 1
fi

INFO "Testing a qualitative colormap with the drawing order reversed, which samples the colormap at \
whole positions rather than at fractions of its range"
args=()
args+=( "--reaction-txt" "draw_kos_samples.reaction.txt" )
args+=( "--reaction-colormap" "tab10" )
args+=( "--reaction-reverse-overlay" )
args+=( "--output-dir" ${output_dir}/draw_txt_samples_qualitative_reverse )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing both layers across the same samples, each with its own aggregation and its own \
sample summary: reactions pooled by their mean and compounds by their maximum concentration"
args=()
args+=( "--reaction-txt" "draw_kos_samples_coverage.reaction.txt" )
args+=( "--compound-txt" "draw_compounds_samples.compound.txt" )
args+=( "--reaction-gene-aggregation" "mean" )
args+=( "--compound-accession-aggregation" "sum" )
args+=( "--reaction-sample-summary" "mean" )
args+=( "--compound-sample-summary" "max" )
args+=( "--output-dir" ${output_dir}/draw_txt_both_layers_samples )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-individual-files" )
args+=( "--draw-grid" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing that each layer's scales are limited on their own, four of them in one run. The \
reaction layer's two scales and the compound layer's per-sample scale are each truncated at both \
ends. The compound layer's unified scale is instead given a minimum equal to the highest value it \
has, which leaves it a single value: its colorbar is one band marked <=, standing for that value \
and everything below it."
args=()
args+=( "--reaction-txt" "draw_kos_samples_coverage.reaction.txt" )
args+=( "--compound-txt" "draw_compounds_samples.compound.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_both_layers_value_limits )
args+=( "--reaction-gene-aggregation" "mean" )
args+=( "--reaction-sample-summary" "mean" )
args+=( "--compound-sample-summary" "max" )
args+=( "--reaction-value-limits" "5" "40" )
args+=( "--reaction-category-value-limits" "5" "40" )
args+=( "--compound-value-limits" "55" "none" )
args+=( "--compound-category-value-limits" "20" "40" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-individual-files" )
args+=( "--draw-grid" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing that each layer's TWO colormaps are chosen on their own as well, giving the run four \
color scales in all: a 'unified' map and a set of per-sample maps for each of the two layers, each \
scale with its own colormap and its own colorbar"
args=()
args+=( "--reaction-txt" "draw_kos_samples_coverage.reaction.txt" )
args+=( "--compound-txt" "draw_compounds_samples.compound.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_both_layers_category_colormap )
args+=( "--reaction-gene-aggregation" "mean" )
args+=( "--reaction-sample-summary" "std" )
args+=( "--compound-sample-summary" "std" )
args+=( "--reaction-colormap" "viridis" )
args+=( "--reaction-category-colormap" "RdYlGn" )
args+=( "--compound-colormap" "cividis" )
args+=( "--compound-category-colormap" "PuOr" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-individual-files" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

for colorbar in reactions reactions_samples compounds compounds_samples
do
    if [ ! -s ${output_dir}/draw_txt_both_layers_category_colormap/colorbar_${colorbar}.pdf ]
    then
        echo "ERROR: two layers with two colormaps each should have written four colorbars."
        exit 1
    fi
done

INFO "Testing an aggregation outside the validated names, resolved through pandas ('var' within \
a sample), together with one that maps disagreement among samples ('std' across them), including \
accessions for which the standard deviation is undefined and are therefore left uncolored"
args=()
args+=( "--reaction-txt" "draw_kos_samples_coverage.reaction.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_samples_kos_std )
args+=( "--reaction-gene-aggregation" "var" )
args+=( "--reaction-sample-summary" "std" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-individual-files" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing the two summary levels together for replicates: each group's map shows the mean of \
its samples' values, while the 'unified' map shows the number of groups containing each reaction, \
drawing map grids and map files for each group"
args=()
args+=( "--reaction-txt" "draw_kos_samples_coverage.reaction.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_groups_kos_quantitative )
args+=( "--groups-txt" "draw-sample-group-information.txt" )
args+=( "--group-threshold" "0.5" )
args+=( "--reaction-sample-summary" "mean" )
args+=( "--reaction-group-summary" "count" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-individual-files" )
args+=( "--draw-grid" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing pooling values at both summary levels, so that the 'unified' map is colored by the \
mean across groups of each group's mean"
args=()
args+=( "--reaction-txt" "draw_kos_samples_coverage.reaction.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_groups_kos_quantitative_mean )
args+=( "--groups-txt" "draw-sample-group-information.txt" )
args+=( "--reaction-sample-summary" "mean" )
args+=( "--reaction-group-summary" "mean" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-individual-files" )
args+=( "--draw-grid" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing a colormap of their own for the maps of the individual groups, whose scale is the \
groups' pooled sample values, alongside a different colormap for the 'unified' map, whose scale \
pools those group values in turn"
args=()
args+=( "--reaction-txt" "draw_kos_samples_coverage.reaction.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_groups_kos_category_colormap )
args+=( "--groups-txt" "draw-sample-group-information.txt" )
args+=( "--reaction-sample-summary" "mean" )
args+=( "--reaction-group-summary" "std" )
args+=( "--reaction-colormap" "cividis" )
args+=( "--reaction-category-colormap" "RdBu_r" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-individual-files" )
args+=( "--draw-grid" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

if [ ! -s ${output_dir}/draw_txt_groups_kos_category_colormap/colorbar_reactions_groups.pdf ]
then
    echo "ERROR: a grouped run with two colormaps should have written a per-group colorbar."
    exit 1
fi

INFO "Testing summarizing the samples of a value-column layer by presence at both levels, where \
the group maps show within-group sample counts and the value column goes unused"
args=()
args+=( "--reaction-txt" "draw_kos_samples_coverage.reaction.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_groups_kos_presence )
args+=( "--groups-txt" "draw-sample-group-information.txt" )
args+=( "--group-threshold" "0" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-individual-files" )
args+=( "--draw-grid" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing value-column coloring with a constant column and 'max' aggregation, exercising the \
degenerate value range and its single-value colorbar"
args=()
args+=( "--reaction-txt" "draw_kos_constant_coverage.reaction.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_kos_quantitative_constant )
args+=( "--reaction-gene-aggregation" "max" )
args+=( "--reaction-accession-aggregation" "max" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing MIXED coloring: presence of KOs (reaction layer) with quantitative compounds \
(compound layer) on the same maps"
args=()
args+=( "--reaction-txt" "draw_kos.reaction.txt" )
args+=( "--compound-txt" "draw_compounds.compound.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_mixed_presence_quant )
# Draws the global map too, for two layers at once, which suppresses derived compound colors.
args+=( "--pathway-numbers" "${global_pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing MIXED coloring across samples: KO membership across samples (reaction layer) with a \
quantitative compound layer held constant across the per-sample maps"
args=()
args+=( "--reaction-txt" "draw_kos_samples.reaction.txt" )
args+=( "--compound-txt" "draw_compounds.compound.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_mixed_samples )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-individual-files" )
args+=( "--draw-grid" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing custom presence colors for the reaction and compound layers of draw-kegg-pathways \
text files"
args=()
args+=( "--reaction-txt" "draw_kos.reaction.txt" )
args+=( "--compound-txt" "draw_compounds_presence.compound.txt" )
args+=( "--reaction-color" "#FF8800" )
args+=( "--compound-color" "#0000FF" )
args+=( "--output-dir" ${output_dir}/draw_txt_custom_colors )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing highlighting a reaction presence layer in the reference map's original colors"
args=()
args+=( "--reaction-txt" "draw_kos.reaction.txt" )
args+=( "--original-color" )
args+=( "--output-dir" ${output_dir}/draw_txt_original_color )
# Draws the global map too, for the reference-color drawer, a separate path from the element engine.
args+=( "--pathway-numbers" "${global_pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing the reference map's original colors across samples: the 'unified' map is the union \
of the samples and each sample also gets its own map"
args=()
args+=( "--reaction-txt" "draw_kos_samples.reaction.txt" )
args+=( "--original-color" )
args+=( "--output-dir" ${output_dir}/draw_txt_original_color_samples )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-individual-files" )
args+=( "--draw-grid" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing the reference map's original colors with sample groups: the 'unified' map is still \
the union of the samples, while each group's map falls back to within-group sample counts, since \
one color cannot distinguish a group's samples"
args=()
args+=( "--reaction-txt" "draw_kos_samples.reaction.txt" )
args+=( "--original-color" )
args+=( "--groups-txt" "draw-sample-group-information.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_original_color_groups )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-individual-files" )
args+=( "--draw-grid" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing samples whose names match the output directory's own subdirectories ('unified', \
'grid', 'symlink') and a BRITE category ('Metabolism'), which are drawn under 'individual' and so \
cannot collide with them"
args=()
args+=( "--reaction-txt" "draw_kos_awkward_sample_names.reaction.txt" )
args+=( "--categorize-files" )
args+=( "--output-dir" ${output_dir}/draw_txt_awkward_sample_names )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-individual-files" )
args+=( "--draw-grid" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing a sample whose name cannot be a directory name, which is summarized on the 'unified' \
map by color alone and so never becomes a path"
args=()
args+=( "--reaction-txt" "draw_kos_unusable_sample_names.reaction.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_unusable_sample_names )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing that individual maps can still be drawn for a subset of samples that leaves out the \
sample whose name cannot be a directory name"
args=()
args+=( "--reaction-txt" "draw_kos_unusable_sample_names.reaction.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_unusable_sample_names_subset )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-individual-files" "SAMPLE_1" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing gathering the maps of individual samples into a directory per map, so that one map \
can be stepped through sample by sample"
args=()
args+=( "--reaction-txt" "draw_kos_samples.reaction.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_collate_by_map )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-individual-files" )
args+=( "--collate-files-by-map" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

for sample in SAMPLE_1 SAMPLE_2 SAMPLE_3
do
    if [ ! -s ${output_dir}/draw_txt_collate_by_map/by_map/ko00010/${sample}.pdf ]
    then
        echo "ERROR: the map gathered by map is missing the file of ${sample}."
        exit 1
    fi
done

# The gathered files are links to the maps already drawn rather than copies of them.
if ! [ ${output_dir}/draw_txt_collate_by_map/by_map/ko00010/SAMPLE_1.pdf \
    -ef ${output_dir}/draw_txt_collate_by_map/individual/SAMPLE_1/ko00010.pdf ]
then
    echo "ERROR: a map gathered by map should be a link to the map drawn for the sample."
    exit 1
fi

INFO "Testing that gathering by map keeps the BRITE subdirectories and the named files used \
everywhere else in the output, and leaves the flat directory of links out of it"
args=()
args+=( "--reaction-txt" "draw_kos_samples.reaction.txt" )
args+=( "--output-dir" ${output_dir}/draw_txt_collate_by_map_categorized )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--name-files" )
args+=( "--categorize-files" )
args+=( "--draw-individual-files" )
args+=( "--draw-grid" )
args+=( "--collate-files-by-map" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

collated_dir=${output_dir}/draw_txt_collate_by_map_categorized/by_map
if [ ! -s ${collated_dir}/Metabolism/Carbohydrate_metabolism/ko00010_Glycolysis_Gluconeogenesis/SAMPLE_1.pdf ]
then
    echo "ERROR: gathering by map did not keep the BRITE subdirectories of the map."
    exit 1
fi
if [ -e ${collated_dir}/symlink ]
then
    echo "ERROR: the flat directory of links was gathered as if it held maps of its own."
    exit 1
fi

INFO "Testing that gathering by map takes only the maps of the groups asked for individually, and \
neither the colorbar beside them nor the maps drawn as grid panels alone"
args=()
args+=( "--reaction-txt" "draw_kos_samples.reaction.txt" )
args+=( "--groups-txt" "draw-sample-group-information.txt" )
args+=( "--group-threshold" "0" )
args+=( "--output-dir" ${output_dir}/draw_txt_collate_by_map_groups )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-individual-files" "G1" )
args+=( "--draw-grid" )
args+=( "--collate-files-by-map" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

collated_dir=${output_dir}/draw_txt_collate_by_map_groups/by_map
if [ ! -s ${collated_dir}/ko00010/G1.pdf ]
then
    echo "ERROR: the group asked for individually was not gathered by map."
    exit 1
fi
if [ -e ${collated_dir}/ko00010/G2.pdf ]
then
    echo "ERROR: a group drawn as a grid panel alone was gathered by map."
    exit 1
fi
if [ -e ${collated_dir}/colorbar ]
then
    echo "ERROR: a group's colorbar was gathered by map as if it were one of its maps."
    exit 1
fi

INFO "Testing that the draw-kegg-pathways text file input rejects invalid files and argument \
combinations (each command below is expected to fail)"

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --collate-files-by-map \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: gathering by map with no individual maps to gather should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos.reaction.txt \
    --draw-individual-files --collate-files-by-map \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: gathering by map for a file with no 'sample' column should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_unusable_sample_names.reaction.txt \
    --draw-individual-files \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: individual maps for a sample whose name cannot be a directory name should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_unusable_sample_names.reaction.txt \
    --draw-grid "BAD/NAME" "SAMPLE_1" \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a grid panel for a sample whose name cannot be a directory name should have failed."
    exit 1
fi
if anvi-draw-kegg-pathways --reaction-txt draw_bad_mixed_kr.reaction.txt \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a reaction file mixing 'K' and 'R' accessions should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --compound-txt draw_bad_compound_gene_id.compound.txt \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a compound file with a 'gene_id' column should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_bad_two_values.reaction.txt \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a file with two candidate value columns should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt --reaction-color "#FF8800" \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--reaction-color' on a value-colored reaction layer should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos.reaction.txt \
    --compound-txt draw_compounds_presence.compound.txt --original-color \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--original-color' with a compound file (reaction-only) should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt --original-color \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--original-color' on a reaction layer with a value column should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt --original-color \
    --reaction-sample-summary count \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a sample summary with '--original-color' should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt --original-color \
    --groups-txt draw-sample-group-information.txt --group-threshold 0 --draw-individual-files \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--group-threshold' with '--original-color' should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos.reaction.txt --original-color --reaction-color "#FF8800" \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--original-color' with '--reaction-color' should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos.reaction.txt --reaction-colormap plasma \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--reaction-colormap' with a single-color presence reaction layer should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --reaction-category-colors draw-sample-colors-incomplete.txt \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a colors file leaving a sample of the run uncolored should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --reaction-category-colors draw-sample-colors-clashing.txt \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a colors file whose combinations cannot be told apart should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --reaction-category-colors draw-sample-colors.txt --reaction-colormap tab10 \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: colors given per category together with a colormap should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --reaction-category-colors draw-sample-colors.txt --reaction-sample-summary count \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: colors given per category with a count summary should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos.reaction.txt \
    --reaction-category-colors draw-sample-colors.txt \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: colors given per category for a layer with no 'sample' column should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --groups-txt draw-sample-group-information.txt --group-threshold 0.5 --draw-individual-files \
    --group-colormap category \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: group ramps built from each group's own color, with no colors given, should have failed."
    exit 1
fi

# The pale end of a group's ramp is white at a span starting from 0, and standard and overview maps
# keep white for their own unhighlighted elements. Caught before any file is written, rather than by
# the drawing code partway through.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --groups-txt draw-sample-group-information.txt --group-threshold 0.5 --draw-individual-files \
    --reaction-category-colors draw-group-colors.txt --group-colormap category 0.0 1.0 \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a group ramp running all the way to white should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt --compound-colormap cividis \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--compound-colormap' with no compound file should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt --count-scale-max nonsense \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--count-scale-max' with a value that is neither a rule nor a number should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt --count-scale-max 0 \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--count-scale-max' below 1 should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples300.reaction.txt \
    --reaction-sample-summary count --count-scale-max total \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: 'count' in discrete bands over a scale of more samples than the colormap has distinct colors should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt --presence-colormap-scheme by_count \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--presence-colormap-scheme' on a text run should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --groups-txt draw-sample-group-information.txt --group-threshold 0.5 \
    --reaction-sample-summary mean --reaction-group-summary mean \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--group-threshold' with no group summary of presence should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --groups-txt draw-sample-group-information.txt \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a group summary of presence without '--group-threshold' should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --reaction-sample-summary mean \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a value sample summary on a layer with no value column should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-sample-summary mean \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a sample summary on a layer with no 'sample' column should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --compound-sample-summary count \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a sample summary naming a layer with no file should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-group-summary count \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a group summary without '--groups-txt' should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --groups-txt draw-sample-group-information.txt --group-threshold 0.5 \
    --reaction-sample-summary count --reaction-group-summary mean \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a value group summary with a presence sample summary should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --groups-txt draw-sample-group-information.txt --group-threshold 0.5 \
    --reaction-sample-summary mean --group-colormap plasma --draw-individual-files \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--group-colormap' with value-colored group maps should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --group-colormap-scheme by_count_continuous \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--group-colormap-scheme' without individual group maps should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_bad_repeated_rows.reaction.txt \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a reaction file repeating an accession and gene should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --compound-txt draw_bad_repeated_compounds.compound.txt \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a compound file repeating a compound should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_reactions.reaction.txt \
    --reaction-gene-aggregation mean \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a gene aggregation on a file with no 'gene_id' column should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-gene-aggregation cumsum \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: an aggregation that transforms rather than reduces should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-gene-aggregation idxmax \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: an aggregation meaning different things at different levels should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples4.reaction.txt \
    --reaction-sample-summary membership \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: 'membership' needing more colors than the colormap has should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos.reaction.txt --draw-individual-files \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--draw-individual-files' with no 'sample' column should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt --reaction-color "#FF8800" \
    --reaction-sample-summary count \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a sample summary with a static presence color should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-colormap plasma 0.9 0.2 \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: reversed colormap limits should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-colormap Greys \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a colormap reaching a color reserved for unhighlighted elements should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos.reaction.txt \
    --reaction-color "#ffffffff" \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: another spelling of a reserved color should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-colormap not_a_colormap \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: an unrecognized colormap name should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos.reaction.txt \
    --genomes-storage TEST-GENOMES.db \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a pangenome option on a text run should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-gene-aggregation not_an_aggregation \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: an unrecognized aggregation name should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-value-limits 40 5 \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: value limits given the wrong way around should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-value-limits none none \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: value limits leaving both ends open should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-value-limits not_a_number 40 \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a value limit that is neither a number nor 'none' should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --compound-txt draw_compounds.compound.txt \
    --reaction-value-limits 5 40 \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: limits for a layer with no file should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos.reaction.txt \
    --reaction-value-limits 5 40 \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: limits for a layer with no value column should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-category-value-limits 5 40 \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: limits on the individual maps' scale where the file has no 'sample' column should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-category-value-limits 5 40 --pathway-numbers "${pathway_numbers[@]}" \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: limits on the individual maps' scale where no individual map is drawn should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-sample-summary count --reaction-value-limits 5 40 \
    --pathway-numbers "${pathway_numbers[@]}" \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: limits on a 'unified' map colored by presence rather than value should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-gene-aggregation mean --reaction-sample-summary mean \
    --reaction-value-limits 1000 2000 --pathway-numbers "${pathway_numbers[@]}" \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: limits lying entirely outside the values should have failed but did not."
    exit 1
fi

# A center is refused in all the same situations a pair of value limits is: no file, no value
# column, no 'sample' column, no individual maps to draw, or a scale that is not colored by value.
# It is refused again where it clashes with the limits on the same scale.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-value-center not_a_number \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a value center that is not a number should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --compound-txt draw_compounds.compound.txt \
    --reaction-value-center \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a center for a layer with no file should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos.reaction.txt \
    --reaction-value-center \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a center for a layer with no value column should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-category-value-center \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a center on the individual maps' scale where the file has no 'sample' column should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-category-value-center --pathway-numbers "${pathway_numbers[@]}" \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a center on the individual maps' scale where no individual map is drawn should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-sample-summary count --reaction-value-center \
    --pathway-numbers "${pathway_numbers[@]}" \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a center on a 'unified' map colored by presence rather than value should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-value-center 100 --reaction-value-limits 5 40 \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a center lying outside the limits on the same scale should have failed but did not."
    exit 1
fi

# The limits truncate first and the center widens afterwards, so a center can push an end of a scale back past
# a limit that was holding a tail in check. That undoes what the limit was for, and is refused rather
# than carried out quietly.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_log_ratio.reaction.txt \
    --reaction-gene-aggregation mean --reaction-accession-aggregation mean \
    --reaction-sample-summary mean --reaction-value-center --reaction-value-limits -1 none \
    --pathway-numbers "${pathway_numbers[@]}" \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: centering a scale past a limit that was truncating should have failed but did not."
    exit 1
fi

# The colormap for the maps of the individual samples or groups is refused wherever there is no such
# scale to color: no file, no value column, no 'sample' column, or a scale that is not colored by
# value.
if anvi-draw-kegg-pathways --compound-txt draw_compounds.compound.txt \
    --reaction-category-colormap viridis \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a category colormap for a layer with no file should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --reaction-category-colormap viridis --pathway-numbers "${pathway_numbers[@]}" \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a category colormap for a layer with no value column should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-category-colormap viridis \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a category colormap where the file has no 'sample' column should have failed."
    exit 1
fi

# Grouped, the per-group maps are colored by value only where the sample summary pools values;
# summarized by presence they show within-group sample counts, which '--group-colormap' styles.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --groups-txt draw-sample-group-information.txt --group-threshold 0.5 \
    --reaction-category-colormap viridis --pathway-numbers "${pathway_numbers[@]}" \
    --draw-individual-files \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a category colormap for group maps colored by presence should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-category-colormap viridis 0.9 0.1 --pathway-numbers "${pathway_numbers[@]}" \
    --draw-individual-files \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: category colormap fractions the wrong way around should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-category-colormap viridis 0.1 0.5 0.9 --pathway-numbers "${pathway_numbers[@]}" \
    --draw-individual-files \
    --output-dir ${output_dir}/draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a category colormap given three numbers should have failed but did not."
    exit 1
fi

INFO "Testing that a database run rejects the group-map coloring options where no individual group \
maps are drawn for them to style, as a text run does (each command below is expected to fail)"

# Without a grouping there are no group maps at all, whatever else is drawn.
if anvi-draw-kegg-pathways --external-genomes external-genomes.txt --draw-grid \
    --group-colormap-scheme by_count_continuous \
    --output-dir ${output_dir}/draw_db_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--group-colormap-scheme' without a grouping should have failed."
    exit 1
fi

# Grouped, but nothing asks for the individual group maps the options would style.
if anvi-draw-kegg-pathways --external-genomes external-genomes.txt \
    --groups-txt contigs-db-group-information.txt --group-threshold 0.5 \
    --group-colormap viridis \
    --output-dir ${output_dir}/draw_db_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--group-colormap' without individual group maps should have failed."
    exit 1
fi

# The whole family is covered, not just the option this check was added for.
if anvi-draw-kegg-pathways --external-genomes external-genomes.txt --draw-grid \
    --group-reverse-overlay \
    --output-dir ${output_dir}/draw_db_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--group-reverse-overlay' without a grouping should have failed."
    exit 1
fi

# A database run colors elements by how many databases or genomes contain them rather than by a
# value, so there is no scale of values for limits to bound.
if anvi-draw-kegg-pathways --external-genomes external-genomes.txt \
    --reaction-value-limits 5 40 \
    --output-dir ${output_dir}/draw_db_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: value limits on a database run should have failed but did not."
    exit 1
fi

# Nor is there one for a second colormap to color: a database run's individual maps show one
# database or genome each, in a single color rather than along a scale.
if anvi-draw-kegg-pathways --external-genomes external-genomes.txt --draw-individual-files \
    --reaction-category-colormap viridis \
    --output-dir ${output_dir}/draw_db_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a category colormap on a database run should have failed but did not."
    exit 1
fi

# A pangenome run takes the same check.
if anvi-draw-kegg-pathways --pan-db TEST-PAN.db --genomes-storage TEST-GENOMES.db \
    --group-colormap-scheme by_count \
    --output-dir ${output_dir}/draw_db_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--group-colormap-scheme' on an ungrouped pangenome should have failed."
    exit 1
fi

# A single database is drawn on one map per pathway, leaving no individual maps to gather by map.
if anvi-draw-kegg-pathways --contigs-dbs E_faecalis_6240.db \
    --draw-individual-files --collate-files-by-map \
    --output-dir ${output_dir}/draw_db_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: gathering by map for a single contigs database should have failed."
    exit 1
fi

INFO "Testing mapping KOs from a genomic contigs database"
args=()
args+=( "--contigs-dbs" "E_faecalis_6240.db" )
args+=( "--output-dir" ${output_dir}/contigs_db_kos )
args+=( "--name-files" )
args+=( "--categorize-files" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" "05310" )
args+=( "--draw-bare-maps" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing mapping KOs from multiple contigs databases, displaying database counts emphasizing \
shared reactions, drawing map grids"
args=()
args+=( "--external-genomes" "external-genomes.txt" )
args+=( "--output-dir" ${output_dir}/contigs_dbs_kos_count )
args+=( "--categorize-files" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-grid" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing mapping KOs from multiple contigs databases, displaying database membership"
args=()
args+=( "--contigs-dbs" \
"E_faecalis_6240.db" \
"E_faecalis_6255.db" \
"E_faecalis_6512.db"
)
args+=( "--output-dir" ${output_dir}/contigs_dbs_kos_membership )
# Draws the global map too, for the contigs database engine.
args+=( "--pathway-numbers" "${global_pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing mapping KOs from multiple contigs databases, coloring database membership from a \
color given per database rather than sampled from a colormap, so that each database's own map \
takes that database's color and matches its band on the 'unified' map"
args=()
args+=( "--contigs-dbs" \
"E_faecalis_6240.db" \
"E_faecalis_6255.db" \
"E_faecalis_6512.db"
)
args+=( "--reaction-category-colors" "draw-contigs-db-colors.txt" )
args+=( "--output-dir" ${output_dir}/contigs_dbs_kos_category_colors )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-individual-files" )
args+=( "--draw-grid" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing colors given per group of contigs databases, with each group's own maps colored by a \
ramp running to that group's color rather than by one colormap shared across the groups"
args=()
args+=( "--contigs-dbs" \
"E_faecalis_6240.db" \
"E_faecalis_6255.db" \
"E_faecalis_6512.db" \
"E_faecalis_6557.db" \
"E_faecalis_6563.db" \
)
args+=( "--groups-txt" ${output_dir}/contigs-db-group-information.txt)
args+=( "--group-threshold" "0" )
args+=( "--reaction-category-colors" "draw-group-colors.txt" )
args+=( "--group-colormap" "category" )
args+=( "--output-dir" ${output_dir}/contigs_dbs_kos_group_category_ramps )
args+=( "--draw-individual-files" )
args+=( "--draw-grid" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing mapping KOs from multiple contigs databases, displaying group membership where KOs \
are in any databases in the group, emphasizing shared reactions, drawing map grids and map files \
for each group"
args=()
args+=( "--contigs-dbs" \
"E_faecalis_6240.db" \
"E_faecalis_6255.db" \
"E_faecalis_6512.db" \
"E_faecalis_6557.db" \
"E_faecalis_6563.db" \
)
args+=( "--groups-txt" ${output_dir}/contigs-db-group-information.txt)
args+=( "--group-threshold" "0" )
args+=( "--output-dir" ${output_dir}/contigs_dbs_kos_group_membership_min_threshold )
args+=( "--categorize-files" )
args+=( "--draw-individual-files" )
args+=( "--draw-grid" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing mapping KOs from multiple contigs databases, displaying group membership where KOs \
are in all databases in the group, emphasizing shared reactions, drawing map grids and map files \
for a subset of the groups"
args=()
args+=( "--contigs-dbs" \
"E_faecalis_6240.db" \
"E_faecalis_6255.db" \
"E_faecalis_6512.db" \
"E_faecalis_6557.db" \
"E_faecalis_6563.db" \
)
args+=( "--groups-txt" ${output_dir}/contigs-db-group-information.txt)
args+=( "--group-threshold" "1" )
args+=( "--output-dir" ${output_dir}/contigs_dbs_kos_group_membership_max_threshold )
args+=( "--draw-individual-files" "G2" )
args+=( "--draw-grid" "G2" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing mapping KOs from a pangenomic database, displaying genome counts emphasizing shared \
reactions"
args=()
args+=( "--pan-db" "TEST-PAN.db" )
args+=( "--genomes-storage" "TEST-GENOMES.db" )
args+=( "--output-dir" ${output_dir}/pan_db_kos_genome_count_emphasize_shared )
# Draws the global map too, for the pangenome engine.
args+=( "--pathway-numbers" "${global_pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing mapping KOs from a pangenomic database, displaying genome counts emphasizing \
unshared reactions, drawing map grids and map files for each genome"
args=()
args+=( "--pan-db" "TEST-PAN.db" )
args+=( "--genomes-storage" "TEST-GENOMES.db" )
args+=( "--output-dir" ${output_dir}/pan_db_kos_genome_count_emphasize_unshared )
args+=( "--categorize-files" )
args+=( "--draw-individual-files" )
args+=( "--draw-grid" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--reaction-colormap" "plasma" "0.1" "0.9")
args+=( "--reaction-reverse-overlay" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing mapping KOs from a pangenomic database, displaying presence/absence"
args=()
args+=( "--pan-db" "TEST-PAN.db" )
args+=( "--genomes-storage" "TEST-GENOMES.db" )
args+=( "--output-dir" ${output_dir}/pan_db_kos_presence_absence )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--reaction-color" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing mapping KOs from a pangenomic database, displaying group count where KOs are in \
any genomes in the group, emphasizing shared reactions, drawing map grids and map files for a \
subset of groups"
args=()
args+=( "--pan-db" "TEST-PAN.db" )
args+=( "--genomes-storage" "TEST-GENOMES.db" )
args+=( "--groups-txt" "pan-group-information.txt" )
args+=( "--group-threshold" "0" )
args+=( "--output-dir" ${output_dir}/pan_db_kos_group_count_emphasize_shared )
args+=( "--categorize-files" )
args+=( "--draw-individual-files" "G1" )
args+=( "--draw-grid" "G1" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing mapping KOs from a pangenomic database, displaying group membership where KOs are \
in any genomes in the group, emphasizing unshared reactions, drawing map grids for each group"
args=()
args+=( "--pan-db" "TEST-PAN.db" )
args+=( "--genomes-storage" "TEST-GENOMES.db" )
args+=( "--groups-txt" "pan-group-information.txt" )
args+=( "--group-threshold" "0" )
args+=( "--output-dir" ${output_dir}/pan_db_kos_group_membership_emphasize_unshared )
args+=( "--draw-grid" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--reaction-colormap" "plasma" "0.1" "0.9")
args+=( "--presence-colormap-scheme" "by_count" )
args+=( "--reaction-reverse-overlay" )
args+=( "--group-colormap" "plasma" "0.1" "0.9" )
args+=( "--group-reverse-overlay" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing mapping KOs from a pangenomic database, coloring reactions by genome count on a \
continuous scale rather than in discrete bands"
args=()
args+=( "--pan-db" "TEST-PAN.db" )
args+=( "--genomes-storage" "TEST-GENOMES.db" )
args+=( "--output-dir" ${output_dir}/pan_db_kos_count_continuous )
args+=( "--presence-colormap-scheme" "by_count_continuous" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"
