#!/bin/bash
# Development-only companion to run_component_tests_for_kegg_mapping.sh. It runs that suite first
# and then adds the exhaustive tests the packaged one deliberately leaves out: the checks that
# anvi-draw-kegg-pathways refuses invalid files and argument combinations, the checks on the wording
# of the warnings it gives, and the colormap and colorbar capacity fallbacks.
#
#   bash run_component_tests_for_kegg_mapping_exhaustive.sh OUTPUT_DIR no_interactive [NUM_THREADS]
#
# This script is deliberately absent from anvio/tests/__init__.py, so neither anvi-self-test nor CI
# runs it. It needs two idioms that the packaged suite does not contain: commands expected to fail,
# and warnings matched against a log. Neither has a shared helper in 00.sh yet.

# 00.sh sets 'files' relative to the tests directory, and the suite sourced below ends by changing
# the working directory to the output directory, so remember where the input files are first.
tests_dir=$(pwd)

source run_component_tests_for_kegg_mapping.sh $1 $2 $3

# Several draws below are piped into 'tee' so that a warning can be matched against the log
# afterwards. Without this, such a pipeline exits with the status of 'tee', which always succeeds,
# so the 'set -e' of 00.sh would not see a failed draw and the test would pass on a crash.
set -o pipefail

INFO "Setting up the input files that only the development tests use"
cp $tests_dir/$files/data/input_files/minimal_enzymes_input.txt .
awk -F'\t' 'NR==1 {print $0} NR>1 {$1=$1".db"; print $0}' OFS='\t' \
pan-group-information.txt > contigs-db-group-information.txt
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
        printf '%s\tWIDE\n' "$sample"
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
            printf '%s\t%s\n' "$accession" "$sample"
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
        printf '%s\tBIG\n' "$sample"
    done
    printf 'SAMPLE_300\tTINY\n'
} > draw-sample-group-information-300.txt
# A constant 'coverage' column. With 'max' at BOTH reduction levels, every per-reaction value reduces
# to the same constant, exercising the degenerate (vmin == vmax) branch of quantitative coloring and
# its single-value colorbar. (Both levels are named explicitly in the test below so that it stays
# degenerate whatever the defaults are; with 'sum', an accession with more genes, or an element with
# more accessions, would differ from the rest.)
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "accession","gene_id","coverage"} NR>1{print $2,$1,1}' \
minimal_enzymes_input.txt > draw_kos_constant_coverage.reaction.txt
printf 'category\tcolor\nSAMPLE_1\t#1f77b4\nSAMPLE_2\t#ff7f0e\nSAMPLE_3\t#2ca02c\nSAMPLE_1, SAMPLE_2\t#7d5ba6\n' > draw-sample-colors-combo.txt
printf 'category\tcolor\nSAMPLE_1\t#1f77b4\nSAMPLE_2\t#ff7f0e\n' > draw-sample-colors-incomplete.txt
# Two samples given the same color, so that two of the seven combinations of three samples come out
# the same and the membership scale cannot label its bands.
printf 'category\tcolor\nSAMPLE_1\t#1f77b4\nSAMPLE_2\t#ff7f0e\nSAMPLE_3\t#1f77b4\n' > draw-sample-colors-clashing.txt
# Samples named after the output directory's own subdirectories and after a BRITE category: these
# are drawn under 'individual', so they cannot collide with anything anvi'o creates for itself.
printf 'accession\tsample\nK00844\tMetabolism\nK00845\tall_maps\nK02358\tgrid\nK00844\tunified\nK02358\tby_map\n' > draw_kos_awkward_sample_names.reaction.txt
# A sample whose name cannot be a directory name, alongside one that can. Only a sample drawn on its
# own maps becomes a subdirectory, so this file is fine until such maps are asked for it.
printf 'accession\tsample\nK00844\tSAMPLE_1\nK00845\tBAD/NAME\nK02358\tSAMPLE_1\n' > draw_kos_unusable_sample_names.reaction.txt
# Malformed files exercising the input-validation guards (each drawn command below is expected to fail).
printf 'accession\nK00844\nR00200\n' > draw_bad_mixed_kr.reaction.txt
printf 'accession\tgene_id\tcoverage\nK00844\tgene_1\t10\nK00844\tgene_1\t30\n' > draw_bad_repeated_rows.reaction.txt
printf 'accession\tcol_a\tcol_b\nK00844\t1\t2\n' > draw_bad_two_values.reaction.txt
printf 'accession\tgene_id\nC00031\tgene_1\n' > draw_bad_compound_gene_id.compound.txt
printf 'accession\tconcentration\nC00031\t10\nC00031\t30\n' > draw_bad_repeated_compounds.compound.txt

# A test that expects a refusal discards the output of the draw, so a fixture that was never written
# would make that test pass for the wrong reason. Fail loudly here instead.
for fixture in contigs-db-group-information.txt draw_kos_samples300.reaction.txt \
              draw_kos_samples60.reaction.txt draw-sample-group-information-60.txt \
              draw_kos_sparse300.reaction.txt draw_compounds_samples300.compound.txt \
              draw-sample-group-information-300.txt draw_kos_constant_coverage.reaction.txt \
              draw-sample-colors-combo.txt draw-sample-colors-incomplete.txt draw-sample-colors-clashing.txt \
              draw_kos_awkward_sample_names.reaction.txt draw_kos_unusable_sample_names.reaction.txt \
              draw_bad_mixed_kr.reaction.txt draw_bad_repeated_rows.reaction.txt \
              draw_bad_two_values.reaction.txt draw_bad_compound_gene_id.compound.txt \
              draw_bad_repeated_compounds.compound.txt
do
    if [ ! -s $fixture ]
    then
        echo "ERROR: the development fixture $fixture was not generated."
        exit 1
    fi
done

# Limits that no value crosses, which therefore truncate nothing: both scales still span their own
# values and neither colorbar is marked, which is what lets one set of limits be carried from
# dataset to dataset without rescaling the ones it was never needed for.
INFO "Testing limits that no value crosses, which therefore truncate nothing"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
                        --output-dir draw_txt_samples_kos_value_limits_inert \
                        --reaction-gene-aggregation mean \
                        --reaction-sample-summary mean \
                        --reaction-value-limits 0 200 \
                        --reaction-category-value-limits 0 200 \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

# Negative limits, which is what a log-transformed value column asks for, with one end of the scale
# left open by 'none'. A layer with no 'sample' column has the one scale, so
# '--reaction-value-limits' bounds it and no second pair of limits is wanted.
INFO "Testing negative limits, which is what a log-transformed value column asks for"
anvi-draw-kegg-pathways --reaction-txt draw_kos_negative.reaction.txt \
                        --output-dir draw_txt_kos_negative_value_limits \
                        --reaction-gene-aggregation mean \
                        --reaction-accession-aggregation mean \
                        --reaction-value-limits -5 none \
                        --pathway-numbers $pathway_numbers \
                        --no-progress

# Limiting one of two scales that both color by value warns about the other, which goes on spanning
# whatever its own values reach.
INFO "Testing that limiting one of two scales that both color by value warns about the other"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --output-dir draw_txt_samples_kos_value_limits_one \
    --reaction-gene-aggregation mean --reaction-sample-summary mean \
    --reaction-value-limits 5 40 --pathway-numbers $pathway_numbers --draw-grid \
    --no-progress 2>&1 | tee draw_txt_value_limits_one.log
if ! tr '\n' ' ' < draw_txt_value_limits_one.log \
    | grep -q "still spans whatever its own values happen to reach"
then
    echo "ERROR: limiting one of two scales colored by value did not warn about the other."
    exit 1
fi

# A centered scale whose longer side is trimmed back by a limit, which is how a centered scale is
# kept from wasting most of its colormap on a lopsided tail: the limits truncate what the values
# reach and the center then widens whichever side falls short, so the scale comes out centered on
# zero, over a good deal less than the values span, with its top marked '>='.
INFO "Testing a centered scale whose longer side is trimmed back by a limit"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_log_ratio.reaction.txt \
                        --output-dir draw_txt_samples_kos_value_center_limits \
                        --reaction-gene-aggregation mean \
                        --reaction-accession-aggregation mean \
                        --reaction-sample-summary mean \
                        --reaction-colormap RdBu_r \
                        --reaction-value-center \
                        --reaction-category-value-center \
                        --reaction-value-limits none 3 \
                        --reaction-category-value-limits none 3 \
                        --pathway-numbers $pathway_numbers \
                        --draw-grid \
                        --no-progress

# Centering one of two scales that both color by value warns about the other, which goes on sitting
# wherever its own values leave it while one colormap colors both.
INFO "Testing that centering one of two scales that both color by value warns about the other"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_log_ratio.reaction.txt \
    --output-dir draw_txt_samples_kos_value_center_one \
    --reaction-gene-aggregation mean --reaction-sample-summary mean \
    --reaction-colormap RdBu_r --reaction-value-center \
    --pathway-numbers $pathway_numbers --draw-grid \
    --no-progress 2>&1 | tee draw_txt_value_center_one.log
if ! tr '\n' ' ' < draw_txt_value_center_one.log \
    | grep -q "still sits wherever its own values leave it"
then
    echo "ERROR: centering one of two scales colored by value did not warn about the other."
    exit 1
fi

# Centering a scale drawn from a sequential colormap warns, that colormap having no middle for a
# centered value to take, and that values lying entirely on one side of the center warn too, half of
# the colormap then going unused.
INFO "Testing that centering a scale drawn from a sequential colormap warns"
anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --output-dir draw_txt_kos_value_center_warnings \
    --reaction-gene-aggregation mean --reaction-accession-aggregation mean \
    --reaction-value-center --pathway-numbers $pathway_numbers \
    --no-progress 2>&1 | tee draw_txt_value_center_warnings.log
for phrase in "is not a diverging one" "all lie on one side of it"
do
    if ! tr '\n' ' ' < draw_txt_value_center_warnings.log | grep -q "$phrase"
    then
        echo "ERROR: centering a sequential colormap over one-sided values did not warn ($phrase)."
        exit 1
    fi
done

# Centering a scale drawn from a diverging colormap trimmed off-center warns, the neutral color no
# longer being in the middle of what is drawn.
INFO "Testing that centering a scale drawn from a diverging colormap trimmed off-center warns"
anvi-draw-kegg-pathways --reaction-txt draw_kos_negative.reaction.txt \
    --output-dir draw_txt_kos_value_center_trim \
    --reaction-gene-aggregation mean --reaction-accession-aggregation mean \
    --reaction-colormap RdBu_r 0.2 1.0 --reaction-value-center -4 \
    --pathway-numbers $pathway_numbers \
    --no-progress 2>&1 | tee draw_txt_value_center_trim.log
if ! tr '\n' ' ' < draw_txt_value_center_trim.log | grep -q "That trimming is not symmetric"
then
    echo "ERROR: centering a scale on an off-center trim of a diverging colormap did not warn."
    exit 1
fi

# A centered scale with nothing left to span, every value being the center itself: the colorbar is a
# single band, and it takes the middle color of the colormap rather than its top color, that middle
# being what the centering asked for.
INFO "Testing a centered scale with nothing left to span, every value being the center itself"
anvi-draw-kegg-pathways --reaction-txt draw_kos_constant_coverage.reaction.txt \
                        --output-dir draw_txt_kos_value_center_degenerate \
                        --reaction-gene-aggregation max \
                        --reaction-accession-aggregation max \
                        --reaction-colormap RdBu_r \
                        --reaction-value-center 1 \
                        --pathway-numbers $pathway_numbers \
                        --no-progress

# With a scale running to more samples than the colormap has distinguishable colors, the default
# presence summary falls back from discrete bands to a continuous count scale rather than failing.
INFO "Testing that with a scale running to more samples than the colormap has distinguishable colors"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples300.reaction.txt \
                        --output-dir draw_txt_samples300_kos_count_continuous \
                        --count-scale-max total \
                        --pathway-numbers $pathway_numbers \
                        --no-progress

# A scale with colors enough for every count but more bands than a colorbar can label falls back to
# the continuous scale too, and that asking for the discrete bands outright draws them with a
# warning that the labels will be too small to read.
INFO "Testing that a scale with colors enough for every count but more bands than a colorbar can"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples60.reaction.txt \
                        --output-dir draw_txt_samples60_bands_unreadable \
                        --count-scale-max total \
                        --pathway-numbers $pathway_numbers \
                        --no-progress 2>&1 | tee draw_txt_samples60_fallback.log
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
    --output-dir draw_txt_samples60_bands_asked_for --count-scale-max total \
    --reaction-sample-summary count --pathway-numbers $pathway_numbers --no-progress 2>&1 \
    | tee draw_txt_samples60_bands.log
if ! tr '\n' ' ' < draw_txt_samples60_bands.log | grep -q "The bands are drawn as asked"
then
    echo "ERROR: insisting on more bands than a colorbar can label did not warn about the labels."
    exit 1
fi

# The same fallback with a reaction layer and a compound layer together: each layer decides for
# itself, so the run reports the fallback once per layer, saying which layer each time.
INFO "Testing the same fallback with a reaction layer and a compound layer together"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples300.reaction.txt \
                        --compound-txt draw_compounds_samples300.compound.txt \
                        --output-dir draw_txt_samples300_two_layers_count_continuous \
                        --count-scale-max total \
                        --pathway-numbers $pathway_numbers \
                        --no-progress 2>&1 | tee draw_txt_samples300_two_layers.log
for element_type in reaction compound
do
    if ! grep -q "The $element_type layer is colored by counts" \
        draw_txt_samples300_two_layers.log
    then
        echo "ERROR: the continuous count scale fallback did not name the $element_type layer."
        exit 1
    fi
done

# A count scale stops at the highest count in the data by default, so that a sparse layer over 300
# samples spreads its colors over the counts that occur instead of squeezing them into the bottom of
# the colormap.
INFO "Testing that a count scale stops at the highest count in the data by default"
anvi-draw-kegg-pathways --reaction-txt draw_kos_sparse300.reaction.txt \
                        --output-dir draw_txt_sparse300_count_scale_observed \
                        --pathway-numbers $pathway_numbers \
                        --no-progress

# Running a count scale to every sample there is rather than to the highest count observed.
INFO "Testing running a count scale to every sample there is rather than to the highest count"
anvi-draw-kegg-pathways --reaction-txt draw_kos_sparse300.reaction.txt \
                        --output-dir draw_txt_sparse300_count_scale_total \
                        --count-scale-max total \
                        --pathway-numbers $pathway_numbers \
                        --no-progress

# Pinning a count scale to a given number, so that separate figures share a scale, with counts above
# it taking the top color.
INFO "Testing pinning a count scale to a given number, so that separate figures share a scale"
anvi-draw-kegg-pathways --reaction-txt draw_kos_sparse300.reaction.txt \
                        --output-dir draw_txt_sparse300_count_scale_number \
                        --count-scale-max 3 \
                        --pathway-numbers $pathway_numbers \
                        --no-progress

# A static presence color on a layer with a 'sample' column: it overrides coloring by sample on the
# 'unified' map, which shows presence in any sample in that one color, while the individual sample
# maps are unaffected.
INFO "Testing a static presence color on a layer with a 'sample' column"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
                        --reaction-color "#FF8800" \
                        --output-dir draw_txt_samples_static_color \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

# A row of the colors file that names a combination of samples, which sets that combination's color
# directly in place of the blend of its members' colors.
INFO "Testing a row of the colors file that names a combination of samples"
# Draws the global map too, for colors given per category, blends included.
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
                        --reaction-category-colors draw-sample-colors-combo.txt \
                        --output-dir draw_txt_samples_category_colors_combo \
                        --pathway-numbers $global_pathway_numbers \
                        --no-progress 2>&1 | tee draw_txt_samples_combo.log
# A combination whose members do not parse back to the run's samples is reported and ignored rather
# than refused, so the run would exit 0 with the row silently unused and the combination drawn in
# the blend it was written to override. Its absence from the log is what says the row was taken.
if tr '\n' ' ' < draw_txt_samples_combo.log \
    | grep -q "name at least one thing that is not a sample of this run"
then
    echo "ERROR: the combination row was ignored rather than recoloring the combination."
    exit 1
fi

# Colors given per sample also color a compound layer, independently of the reaction layer, so that
# the two layers drawn on one map keep their own scales.
INFO "Testing that colors given per sample also color a compound layer"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
                        --compound-txt draw_compounds_samples.compound.txt \
                        --reaction-category-colors draw-sample-colors.txt \
                        --compound-category-colors draw-sample-colors-combo.txt \
                        --output-dir draw_txt_two_layers_category_colors \
                        --pathway-numbers $pathway_numbers \
                        --no-progress

# The same group ramps with an explicit span, which says how far from white each ramp starts and
# stops rather than which fraction of a colormap to use.
INFO "Testing the same group ramps with an explicit span"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
                        --groups-txt draw-sample-group-information.txt \
                        --group-threshold 0.5 \
                        --reaction-category-colors draw-group-colors.txt \
                        --group-colormap category 0.4 0.95 \
                        --group-reverse-overlay \
                        --output-dir draw_txt_groups_category_ramps_span \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --no-progress

# '--group-colormap-scheme by_count_continuous', which draws every group's count scale as a gradient
# rather than in discrete bands of one color per source, here for a ramp built from each group's own
# color, so that a group's colorbar spans that group's own hue.
INFO "Testing '--group-colormap-scheme by_count_continuous'"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
                        --groups-txt draw-sample-group-information.txt \
                        --group-threshold 0.5 \
                        --reaction-category-colors draw-group-colors.txt \
                        --group-colormap category \
                        --group-colormap-scheme by_count_continuous \
                        --output-dir draw_txt_groups_scheme_continuous \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

# A group whose count scale runs over more counts than its colorbar can label falls back to the
# continuous scale too, even though its colormap has a distinguishable color for every one of those
# counts.
INFO "Testing that a group whose count scale runs over more counts than its colorbar can label"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples60.reaction.txt \
                        --groups-txt draw-sample-group-information-60.txt \
                        --group-threshold 0.5 \
                        --count-scale-max total \
                        --output-dir draw_txt_groups60_bands_unreadable \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --no-progress 2>&1 | tee draw_txt_groups60_fallback.log
flattened_groups60_log=$(tr '\n' ' ' < draw_txt_groups60_fallback.log)
if ! grep -q "more labels than it can set in type large enough to read" \
    <<< "$flattened_groups60_log"
then
    echo "ERROR: a group scale of unlabelable bands did not fall back to the continuous scale."
    exit 1
fi
if grep -q "could supply only" <<< "$flattened_groups60_log"
then
    echo "ERROR: the fallback blamed the colors when it was the labels that did not fit."
    exit 1
fi

# Asking for the discrete bands outright keeps them for a group whose scale runs over more counts
# than its colorbar can label, warning that the labels will be small -- the colormap has a
# distinguishable color for every one of those counts, so nothing stops them being drawn.
INFO "Testing that asking for the discrete bands outright keeps them for a group whose scale runs"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples60.reaction.txt \
                        --groups-txt draw-sample-group-information-60.txt \
                        --group-threshold 0.5 \
                        --count-scale-max total \
                        --group-colormap-scheme by_count \
                        --output-dir draw_txt_groups60_scheme_bands \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --no-progress 2>&1 | tee draw_txt_groups60_scheme_bands.log
flattened_groups60_bands_log=$(tr '\n' ' ' < draw_txt_groups60_scheme_bands.log)
if ! grep -q "The bands are drawn as asked" <<< "$flattened_groups60_bands_log"
then
    echo "ERROR: insisting on bands too many to label did not warn that they were drawn anyway."
    exit 1
fi
if grep -q "drawn on a CONTINUOUS color scale" <<< "$flattened_groups60_bands_log"
then
    echo "ERROR: the discrete bands were asked for but the continuous scale was drawn."
    exit 1
fi
# The colormap has a color for every count here, so the OTHER limit must not be what is reported.
if grep -q "could supply only" <<< "$flattened_groups60_bands_log"
then
    echo "ERROR: the warning blamed the colors when it was the labels that did not fit."
    exit 1
fi

# A group whose count scale runs over more counts than its colormap can distinguish falls back to
# the continuous scale on its own, with a warning, and that asking for the discrete bands outright
# keeps them and warns that neighboring counts share a color.
INFO "Testing that a group whose count scale runs over more counts than its colormap can"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples300.reaction.txt \
                        --groups-txt draw-sample-group-information-300.txt \
                        --group-threshold 0.5 \
                        --count-scale-max total \
                        --output-dir draw_txt_groups300_scheme_fallback \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --no-progress 2>&1 | tee draw_txt_groups300_scheme_fallback.log
# Warnings are wrapped to the terminal width, so a phrase is matched against the log with its line
# breaks flattened to spaces rather than line by line.
if ! tr '\n' ' ' < draw_txt_groups300_scheme_fallback.log \
    | grep -q "drawn on a CONTINUOUS color scale"
then
    echo "ERROR: a group scale short of colors did not fall back to the continuous scale."
    exit 1
fi
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples300.reaction.txt \
                        --groups-txt draw-sample-group-information-300.txt \
                        --group-threshold 0.5 \
                        --count-scale-max total \
                        --group-colormap-scheme by_count \
                        --output-dir draw_txt_groups300_scheme_bands \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --no-progress 2>&1 | tee draw_txt_groups300_scheme_bands.log
flattened_bands_log=$(tr '\n' ' ' < draw_txt_groups300_scheme_bands.log)
if ! grep -q "Neighboring counts therefore share a color" <<< "$flattened_bands_log"
then
    echo "ERROR: insisting on the discrete bands did not warn that counts share a color."
    exit 1
fi
if grep -q "drawn on a CONTINUOUS color scale" <<< "$flattened_bands_log"
then
    echo "ERROR: the discrete bands were asked for but the continuous scale was drawn."
    exit 1
fi

# A qualitative colormap with the drawing order reversed, which samples the colormap at whole
# positions rather than at fractions of its range.
INFO "Testing a qualitative colormap with the drawing order reversed"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
                        --reaction-colormap tab10 \
                        --reaction-reverse-overlay \
                        --output-dir draw_txt_samples_qualitative_reverse \
                        --pathway-numbers $pathway_numbers \
                        --no-progress

# Each layer's scales are limited on their own, four of them in one run. The reaction layer's two
# scales and the compound layer's per-sample scale are each truncated at both ends. The compound
# layer's unified scale is instead given a minimum equal to the highest value it has, which leaves
# it a single value: its colorbar is one band marked <=, standing for that value and everything
# below it.
INFO "Testing that each layer's scales are limited on their own"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
                        --compound-txt draw_compounds_samples.compound.txt \
                        --output-dir draw_txt_both_layers_value_limits \
                        --reaction-gene-aggregation mean \
                        --reaction-sample-summary mean \
                        --compound-sample-summary max \
                        --reaction-value-limits 5 40 \
                        --reaction-category-value-limits 5 40 \
                        --compound-value-limits 55 none \
                        --compound-category-value-limits 20 40 \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

# Each layer's TWO colormaps are chosen on their own as well, giving the run four color scales in
# all: a 'unified' map and a set of per-sample maps for each of the two layers, each scale with its
# own colormap and its own colorbar.
INFO "Testing that each layer's TWO colormaps are chosen on their own as well"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
                        --compound-txt draw_compounds_samples.compound.txt \
                        --output-dir draw_txt_both_layers_category_colormap \
                        --reaction-gene-aggregation mean \
                        --reaction-sample-summary std \
                        --compound-sample-summary std \
                        --reaction-colormap viridis \
                        --reaction-category-colormap RdYlGn \
                        --compound-colormap cividis \
                        --compound-category-colormap PuOr \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --no-progress

for colorbar in reactions reactions_samples compounds compounds_samples
do
    if [ ! -s draw_txt_both_layers_category_colormap/colorbar_$colorbar.pdf ]
    then
        echo "ERROR: two layers with two colormaps each should have written four colorbars."
        exit 1
    fi
done

# An aggregation outside the validated names, resolved through pandas ('var' within a sample),
# together with one that maps disagreement among samples ('std' across them), including accessions
# for which the standard deviation is undefined and are therefore left uncolored.
INFO "Testing an aggregation outside the validated names"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
                        --output-dir draw_txt_samples_kos_std \
                        --reaction-gene-aggregation var \
                        --reaction-sample-summary std \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --no-progress

# The two summary levels together for replicates: each group's map shows the mean of its samples'
# values, while the 'unified' map shows the number of groups containing each reaction, drawing map
# grids and map files for each group.
INFO "Testing the two summary levels together for replicates"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
                        --output-dir draw_txt_groups_kos_quantitative \
                        --groups-txt draw-sample-group-information.txt \
                        --group-threshold 0.5 \
                        --reaction-sample-summary mean \
                        --reaction-group-summary count \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

# Pooling values at both summary levels, so that the 'unified' map is colored by the mean across
# groups of each group's mean.
INFO "Testing pooling values at both summary levels"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
                        --output-dir draw_txt_groups_kos_quantitative_mean \
                        --groups-txt draw-sample-group-information.txt \
                        --reaction-sample-summary mean \
                        --reaction-group-summary mean \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

# Summarizing the samples of a value-column layer by presence at both levels, where the group maps
# show within-group sample counts and the value column goes unused.
INFO "Testing summarizing the samples of a value-column layer by presence at both levels"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
                        --output-dir draw_txt_groups_kos_presence \
                        --groups-txt draw-sample-group-information.txt \
                        --group-threshold 0 \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

# Value-column coloring with a constant column and 'max' aggregation, exercising the degenerate
# value range and its single-value colorbar.
INFO "Testing value-column coloring with a constant column and 'max' aggregation"
anvi-draw-kegg-pathways --reaction-txt draw_kos_constant_coverage.reaction.txt \
                        --output-dir draw_txt_kos_quantitative_constant \
                        --reaction-gene-aggregation max \
                        --reaction-accession-aggregation max \
                        --pathway-numbers $pathway_numbers \
                        --no-progress

# The reference map's original colors across samples: the 'unified' map is the union of the samples
# and each sample also gets its own map.
INFO "Testing the reference map's original colors across samples"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
                        --original-color \
                        --output-dir draw_txt_original_color_samples \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

# The reference map's original colors with sample groups: the 'unified' map is still the union of
# the samples, while each group's map falls back to within-group sample counts, since one color
# cannot distinguish a group's samples.
INFO "Testing the reference map's original colors with sample groups"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
                        --original-color \
                        --groups-txt draw-sample-group-information.txt \
                        --output-dir draw_txt_original_color_groups \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

# Samples whose names match the output directory's own subdirectories ('unified', 'grid', 'by_map',
# 'all_maps') and a BRITE category ('Metabolism'), which are drawn under 'individual' and so cannot
# collide with them.
INFO "Testing samples whose names match the output directory's own subdirectories ('unified', 'grid'"
anvi-draw-kegg-pathways --reaction-txt draw_kos_awkward_sample_names.reaction.txt \
                        --categorize-files \
                        --output-dir draw_txt_awkward_sample_names \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --no-progress

# A sample whose name cannot be a directory name, which is summarized on the 'unified' map by color
# alone and so never becomes a path.
INFO "Testing a sample whose name cannot be a directory name"
anvi-draw-kegg-pathways --reaction-txt draw_kos_unusable_sample_names.reaction.txt \
                        --output-dir draw_txt_unusable_sample_names \
                        --pathway-numbers $pathway_numbers \
                        --no-progress

# Individual maps can still be drawn for a subset of samples that leaves out the sample whose name
# cannot be a directory name.
INFO "Testing that individual maps can still be drawn for a subset of samples that leaves out the"
anvi-draw-kegg-pathways --reaction-txt draw_kos_unusable_sample_names.reaction.txt \
                        --output-dir draw_txt_unusable_sample_names_subset \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files SAMPLE_1 \
                        --no-progress

# The draw-kegg-pathways text file input rejects invalid files and argument combinations (each
# command below is expected to fail).
INFO "Testing that the draw-kegg-pathways text file input rejects invalid files and argument"
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --collate-files-by-map \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: gathering by map with no individual maps to gather should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos.reaction.txt \
    --draw-individual-files --collate-files-by-map \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: gathering by map for a file with no 'sample' column should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_unusable_sample_names.reaction.txt \
    --draw-individual-files \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: individual maps for a sample whose name cannot be a directory name should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_unusable_sample_names.reaction.txt \
    --draw-grid "BAD/NAME" "SAMPLE_1" \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a grid panel for a sample whose name cannot be a directory name should have failed."
    exit 1
fi
if anvi-draw-kegg-pathways --reaction-txt draw_bad_mixed_kr.reaction.txt \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a reaction file mixing 'K' and 'R' accessions should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --compound-txt draw_bad_compound_gene_id.compound.txt \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a compound file with a 'gene_id' column should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_bad_two_values.reaction.txt \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a file with two candidate value columns should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt --reaction-color "#FF8800" \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--reaction-color' on a value-colored reaction layer should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos.reaction.txt \
    --compound-txt draw_compounds_presence.compound.txt --original-color \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--original-color' with a compound file (reaction-only) should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt --original-color \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--original-color' on a reaction layer with a value column should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt --original-color \
    --reaction-sample-summary count \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a sample summary with '--original-color' should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt --original-color \
    --groups-txt draw-sample-group-information.txt --group-threshold 0 --draw-individual-files \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--group-threshold' with '--original-color' should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos.reaction.txt --original-color --reaction-color "#FF8800" \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--original-color' with '--reaction-color' should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos.reaction.txt --reaction-colormap plasma \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--reaction-colormap' with a single-color presence reaction layer should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --reaction-category-colors draw-sample-colors-incomplete.txt \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a colors file leaving a sample of the run uncolored should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --reaction-category-colors draw-sample-colors-clashing.txt \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a colors file whose combinations cannot be told apart should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --reaction-category-colors draw-sample-colors.txt --reaction-colormap tab10 \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: colors given per category together with a colormap should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --reaction-category-colors draw-sample-colors.txt --reaction-sample-summary count \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: colors given per category with a count summary should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos.reaction.txt \
    --reaction-category-colors draw-sample-colors.txt \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: colors given per category for a layer with no 'sample' column should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --groups-txt draw-sample-group-information.txt --group-threshold 0.5 --draw-individual-files \
    --group-colormap category \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
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
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a group ramp running all the way to white should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt --compound-colormap cividis \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--compound-colormap' with no compound file should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt --count-scale-max nonsense \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--count-scale-max' with a value that is neither a rule nor a number should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt --count-scale-max 0 \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--count-scale-max' below 1 should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples300.reaction.txt \
    --reaction-sample-summary count --count-scale-max total \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: 'count' in discrete bands over a scale of more samples than the colormap has distinct colors should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt --presence-colormap-scheme by_count \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--presence-colormap-scheme' on a text run should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --groups-txt draw-sample-group-information.txt --group-threshold 0.5 \
    --reaction-sample-summary mean --reaction-group-summary mean \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--group-threshold' with no group summary of presence should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --groups-txt draw-sample-group-information.txt \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a group summary of presence without '--group-threshold' should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --reaction-sample-summary mean \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a value sample summary on a layer with no value column should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-sample-summary mean \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a sample summary on a layer with no 'sample' column should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --compound-sample-summary count \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a sample summary naming a layer with no file should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-group-summary count \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a group summary without '--groups-txt' should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --groups-txt draw-sample-group-information.txt --group-threshold 0.5 \
    --reaction-sample-summary count --reaction-group-summary mean \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a value group summary with a presence sample summary should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --groups-txt draw-sample-group-information.txt --group-threshold 0.5 \
    --reaction-sample-summary mean --group-colormap plasma --draw-individual-files \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--group-colormap' with value-colored group maps should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --group-colormap-scheme by_count_continuous \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--group-colormap-scheme' without individual group maps should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_bad_repeated_rows.reaction.txt \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a reaction file repeating an accession and gene should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --compound-txt draw_bad_repeated_compounds.compound.txt \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a compound file repeating a compound should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_reactions.reaction.txt \
    --reaction-gene-aggregation mean \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a gene aggregation on a file with no 'gene_id' column should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-gene-aggregation cumsum \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: an aggregation that transforms rather than reduces should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-gene-aggregation idxmax \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: an aggregation meaning different things at different levels should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples4.reaction.txt \
    --reaction-sample-summary membership \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: 'membership' needing more colors than the colormap has should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos.reaction.txt --draw-individual-files \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--draw-individual-files' with no 'sample' column should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt --reaction-color "#FF8800" \
    --reaction-sample-summary count \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a sample summary with a static presence color should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-colormap plasma 0.9 0.2 \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: reversed colormap limits should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-colormap Greys \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a colormap reaching a color reserved for unhighlighted elements should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos.reaction.txt \
    --reaction-color "#ffffffff" \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: another spelling of a reserved color should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-colormap not_a_colormap \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: an unrecognized colormap name should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos.reaction.txt \
    --genomes-storage TEST-GENOMES.db \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a pangenome option on a text run should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-gene-aggregation not_an_aggregation \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: an unrecognized aggregation name should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-value-limits 40 5 \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: value limits given the wrong way around should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-value-limits none none \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: value limits leaving both ends open should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-value-limits not_a_number 40 \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a value limit that is neither a number nor 'none' should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --compound-txt draw_compounds.compound.txt \
    --reaction-value-limits 5 40 \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: limits for a layer with no file should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos.reaction.txt \
    --reaction-value-limits 5 40 \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: limits for a layer with no value column should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-category-value-limits 5 40 \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: limits on the individual maps' scale where the file has no 'sample' column should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-category-value-limits 5 40 --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: limits on the individual maps' scale where no individual map is drawn should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-sample-summary count --reaction-value-limits 5 40 \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: limits on a 'unified' map colored by presence rather than value should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-gene-aggregation mean --reaction-sample-summary mean \
    --reaction-value-limits 1000 2000 --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: limits lying entirely outside the values should have failed but did not."
    exit 1
fi

# A center is refused in all the same situations a pair of value limits is: no file, no value
# column, no 'sample' column, no individual maps to draw, or a scale that is not colored by value.
# It is refused again where it clashes with the limits on the same scale.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-value-center not_a_number \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a value center that is not a number should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --compound-txt draw_compounds.compound.txt \
    --reaction-value-center \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a center for a layer with no file should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos.reaction.txt \
    --reaction-value-center \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a center for a layer with no value column should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-category-value-center \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a center on the individual maps' scale where the file has no 'sample' column should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-category-value-center --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a center on the individual maps' scale where no individual map is drawn should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-sample-summary count --reaction-value-center \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a center on a 'unified' map colored by presence rather than value should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-value-center 100 --reaction-value-limits 5 40 \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
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
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: centering a scale past a limit that was truncating should have failed but did not."
    exit 1
fi

# The colormap for the maps of the individual samples or groups is refused wherever there is no such
# scale to color: no file, no value column, no 'sample' column, or a scale that is not colored by
# value.
if anvi-draw-kegg-pathways --compound-txt draw_compounds.compound.txt \
    --reaction-category-colormap viridis \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a category colormap for a layer with no file should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --reaction-category-colormap viridis --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a category colormap for a layer with no value column should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-category-colormap viridis \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a category colormap where the file has no 'sample' column should have failed."
    exit 1
fi

# Grouped, the per-group maps are colored by value only where the sample summary pools values;
# summarized by presence they show within-group sample counts, which '--group-colormap' styles.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --groups-txt draw-sample-group-information.txt --group-threshold 0.5 \
    --reaction-category-colormap viridis --pathway-numbers $pathway_numbers \
    --draw-individual-files \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a category colormap for group maps colored by presence should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-category-colormap viridis 0.9 0.1 --pathway-numbers $pathway_numbers \
    --draw-individual-files \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: category colormap fractions the wrong way around should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-category-colormap viridis 0.1 0.5 0.9 --pathway-numbers $pathway_numbers \
    --draw-individual-files \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a category colormap given three numbers should have failed but did not."
    exit 1
fi

# A database run rejects the group-map coloring options where no individual group maps are drawn for
# them to style, as a text run does (each command below is expected to fail).
INFO "Testing that a database run rejects the group-map coloring options where no individual"
# Without a grouping there are no group maps at all, whatever else is drawn.
if anvi-draw-kegg-pathways --external-genomes external-genomes.txt --draw-grid \
    --group-colormap-scheme by_count_continuous \
    --output-dir draw_db_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--group-colormap-scheme' without a grouping should have failed."
    exit 1
fi

# Grouped, but nothing asks for the individual group maps the options would style.
if anvi-draw-kegg-pathways --external-genomes external-genomes.txt \
    --groups-txt contigs-db-group-information.txt --group-threshold 0.5 \
    --group-colormap viridis \
    --output-dir draw_db_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--group-colormap' without individual group maps should have failed."
    exit 1
fi

# The whole family is covered, not just the option this check was added for.
if anvi-draw-kegg-pathways --external-genomes external-genomes.txt --draw-grid \
    --group-reverse-overlay \
    --output-dir draw_db_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--group-reverse-overlay' without a grouping should have failed."
    exit 1
fi

# A database run colors elements by how many databases or genomes contain them rather than by a
# value, so there is no scale of values for limits to bound.
if anvi-draw-kegg-pathways --external-genomes external-genomes.txt \
    --reaction-value-limits 5 40 \
    --output-dir draw_db_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: value limits on a database run should have failed but did not."
    exit 1
fi

# Nor is there one for a second colormap to color: a database run's individual maps show one
# database or genome each, in a single color rather than along a scale.
if anvi-draw-kegg-pathways --external-genomes external-genomes.txt --draw-individual-files \
    --reaction-category-colormap viridis \
    --output-dir draw_db_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a category colormap on a database run should have failed but did not."
    exit 1
fi

# A pangenome run takes the same check.
if anvi-draw-kegg-pathways --pan-db TEST-PAN.db --genomes-storage TEST-GENOMES.db \
    --group-colormap-scheme by_count \
    --output-dir draw_db_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--group-colormap-scheme' on an ungrouped pangenome should have failed."
    exit 1
fi

# A single database is drawn on one map per pathway, leaving no individual maps to gather by map.
if anvi-draw-kegg-pathways --contigs-dbs E_faecalis_6240.db \
    --draw-individual-files --collate-files-by-map \
    --output-dir draw_db_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: gathering by map for a single contigs database should have failed."
    exit 1
fi

# Colors given per group of contigs databases, with each group's own maps colored by a ramp running
# to that group's color rather than by one colormap shared across the groups.
INFO "Testing colors given per group of contigs databases"
anvi-draw-kegg-pathways --contigs-dbs E_faecalis_6240.db E_faecalis_6255.db E_faecalis_6512.db E_faecalis_6557.db E_faecalis_6563.db \
                        --groups-txt contigs-db-group-information.txt \
                        --group-threshold 0 \
                        --reaction-category-colors draw-group-colors.txt \
                        --group-colormap category \
                        --output-dir contigs_dbs_kos_group_category_ramps \
                        --draw-individual-files \
                        --draw-grid \
                        --pathway-numbers $pathway_numbers \
                        --no-progress

# Mapping KOs from multiple contigs databases, displaying group membership where KOs are in any
# databases in the group, emphasizing shared reactions, drawing map grids and map files for each
# group.
INFO "Testing mapping KOs from multiple contigs databases"
anvi-draw-kegg-pathways --contigs-dbs E_faecalis_6240.db E_faecalis_6255.db E_faecalis_6512.db E_faecalis_6557.db E_faecalis_6563.db \
                        --groups-txt contigs-db-group-information.txt \
                        --group-threshold 0 \
                        --output-dir contigs_dbs_kos_group_membership_min_threshold \
                        --categorize-files \
                        --draw-individual-files \
                        --draw-grid \
                        --pathway-numbers $pathway_numbers \
                        --no-progress

# Mapping KOs from multiple contigs databases, displaying group membership where KOs are in all
# databases in the group, emphasizing shared reactions, drawing map grids and map files for a subset
# of the groups.
INFO "Testing mapping KOs from multiple contigs databases"
anvi-draw-kegg-pathways --contigs-dbs E_faecalis_6240.db E_faecalis_6255.db E_faecalis_6512.db E_faecalis_6557.db E_faecalis_6563.db \
                        --groups-txt contigs-db-group-information.txt \
                        --group-threshold 1 \
                        --output-dir contigs_dbs_kos_group_membership_max_threshold \
                        --draw-individual-files G2 \
                        --draw-grid G2 \
                        --pathway-numbers $pathway_numbers \
                        --no-progress

# Mapping KOs from a pangenomic database, displaying genome counts emphasizing unshared reactions,
# drawing map grids and map files for each genome.
INFO "Testing mapping KOs from a pangenomic database"
anvi-draw-kegg-pathways --pan-db TEST-PAN.db \
                        --genomes-storage TEST-GENOMES.db \
                        --output-dir pan_db_kos_genome_count_emphasize_unshared \
                        --categorize-files \
                        --draw-individual-files \
                        --draw-grid \
                        --pathway-numbers $pathway_numbers \
                        --reaction-colormap plasma 0.1 0.9 \
                        --reaction-reverse-overlay \
                        --no-progress

# Mapping KOs from a pangenomic database, displaying group membership where KOs are in any genomes
# in the group, emphasizing unshared reactions, drawing map grids for each group.
INFO "Testing mapping KOs from a pangenomic database"
anvi-draw-kegg-pathways --pan-db TEST-PAN.db \
                        --genomes-storage TEST-GENOMES.db \
                        --groups-txt pan-group-information.txt \
                        --group-threshold 0 \
                        --output-dir pan_db_kos_group_membership_emphasize_unshared \
                        --draw-grid \
                        --pathway-numbers $pathway_numbers \
                        --reaction-colormap plasma 0.1 0.9 \
                        --presence-colormap-scheme by_count \
                        --reaction-reverse-overlay \
                        --group-colormap plasma 0.1 0.9 \
                        --group-reverse-overlay \
                        --no-progress
