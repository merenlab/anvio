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
# Log abundances spread over three samples, negative throughout. Every element's mean is therefore
# negative, and a ratio measured against it is undefined for all of them, so a normalization leaves
# every map of an individual sample blank.
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "accession","gene_id","sample","log_abundance"} NR>1{s=((NR-2)%3)+1; print $1,$2,sprintf("SAMPLE_%d",s),$3}' \
draw_kos_negative.reaction.txt > draw_kos_samples_negative.reaction.txt
printf 'sample\tcolor\nSAMPLE_1\t#1f77b4\nSAMPLE_2\t#ff7f0e\nSAMPLE_3\t#2ca02c\nSAMPLE_1, SAMPLE_2\t#7d5ba6\n' > draw-sample-colors-combo.txt
printf 'sample\tcolor\nSAMPLE_1\t#1f77b4\nSAMPLE_2\t#ff7f0e\n' > draw-sample-colors-incomplete.txt
# Two samples given the same color, so that two of the seven combinations of three samples come out
# the same and the membership scale cannot label its bands.
printf 'sample\tcolor\nSAMPLE_1\t#1f77b4\nSAMPLE_2\t#ff7f0e\nSAMPLE_3\t#1f77b4\n' > draw-sample-colors-clashing.txt
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
# Two columns of one name, which pandas would read as 'coverage' and 'coverage.1', or as 'sample'
# and 'sample.1'. A column name that starts with '__', like the columns that anvi'o adds as it reads
# a file.
printf 'accession\tcoverage\tcoverage\nK00844\t1\t10\nK01810\t2\t20\n' > draw_bad_repeated_weights.reaction.txt
printf 'accession\tsample\tsample\nK00844\tS1\tup\nK01810\tS2\tdown\n' > draw_bad_repeated_sample.reaction.txt
printf 'accession\t__accession\nK00844\t5\nK01810\t7\n' > draw_bad_reserved_header.reaction.txt
# A number with a stray opening quote.
printf 'accession\tcoverage\nK00844\t"10\nK01810\t20\n' > draw_bad_stray_quote_number.reaction.txt
# Categorical values without a 'sample' column, so that one map shows them. Each row cycles through
# three classes, so the genes of an accession and the KOs of an element often differ.
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "accession","gene_id","class"} NR>1{split("alpha beta gamma", c, " "); print $1,$2,c[((NR-2)%3)+1]}' \
draw_kos.reaction.txt > draw_kos_class.reaction.txt
# Categorical values of compounds across samples.
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "accession","sample","status"} NR>1{split("up down", s, " "); print $1,$2,s[((NR-2)%2)+1]}' \
draw_compounds_samples.compound.txt > draw_compounds_samples_status.compound.txt
# Eighteen different values, as many as anvi'o has default colors for.
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "accession","gene_id","class"} NR>1{print $1,$2,sprintf("class %d",((NR-2)%18)+1)}' \
draw_kos.reaction.txt > draw_kos_class18.reaction.txt
# Nineteen different values, one more than anvi'o has default colors for.
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "accession","gene_id","class"} NR>1{print $1,$2,sprintf("class %d",((NR-2)%19)+1)}' \
draw_kos.reaction.txt > draw_kos_class19.reaction.txt
# The same nineteen values across samples.
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "accession","gene_id","sample","class"} NR>1{print $1,$2,$3,sprintf("class %d",((NR-2)%19)+1)}' \
draw_kos_samples.reaction.txt > draw_kos_samples_class19.reaction.txt
# A colors file missing a value, and one giving two values the same color.
printf 'status\tcolor\nup\t#d62728\ndown\t#1f77b4\n' > draw-status-colors-incomplete.txt
printf 'status\tcolor\nup\t#d62728\ndown\t#d62728\nsame\t#bcbd22\n' > draw-status-colors-shared.txt
# A row without a categorical value, and two rows giving one accession and gene different values.
printf 'accession\tstatus\nK00844\tup\nK00845\t\n' > draw_bad_blank_status.reaction.txt
printf 'accession\tgene_id\tstatus\nK00844\tgene_1\tup\nK00844\tgene_1\tdown\n' > draw_bad_repeated_status.reaction.txt
# A file of categorical values with a header and no rows.
printf 'accession\tstatus\n' > draw_kos_header_only_status.reaction.txt
# One KO with a different value in each of two samples, so that 'plurality' ties across them.
printf 'accession\tsample\tstatus\nK00845\tS2\tup\nK00845\tS3\tdown\n' > draw_kos_tie.reaction.txt
# Two KOs that share a map element in 00010. One KO has two genes with the value 'up'. The other has
# one gene with the value 'down'.
printf 'accession\tgene_id\tstatus\nK01623\tgene_1\tup\nK01623\tgene_2\tup\nK01624\tgene_3\tdown\n' > draw_kos_gene_votes.reaction.txt
# A value with a stray opening quote and a value with a stray closing quote, and a colors file that
# names the values with their quotes.
printf 'accession\tstatus\nK00844\t"high\nK01810\tlow"\nK00850\tmid\n' > draw_kos_stray_quote.reaction.txt
printf 'status\tcolor\n"high\t#d62728\nlow"\t#1f77b4\nmid\t#bcbd22\n' > draw-stray-quote-colors.txt
# Every field in double quotes, as R's write.table writes a file by default.
printf '"accession"\t"status"\n"K00844"\t"up"\n"K01810"\t"down"\n' > draw_kos_quoted_fields.reaction.txt
# A stray opening quote in the name of a column of categorical values and of a column of numbers.
printf 'accession\t"status\nK00844\tup\nK01810\tdown\n' > draw_kos_header_quote_status.reaction.txt
printf 'accession\t"coverage\nK00844\t10\nK01810\t20\n' > draw_kos_header_quote_coverage.reaction.txt
# Weights without a 'sample' column. Element 18 of map 00010 is sliced between K01624 and K01622,
# since K16306 has a weight of 0. Both KOs of element 154 have a weight of 0, so that element is
# left uncolored.
printf 'accession\tcoverage\tweight\nK00174\t1.5\t0\nK00175\t3.0\t0\nK01624\t2.0\t3\nK01622\t4.0\t1\nK16306\t6.0\t0\n' > draw_kos_weights.reaction.txt
# Four of the five KOs of entry 6838 of the global map, a line too short for four slices, the three
# KOs of entry 1327, a line whose visible part is too short for three slices, and KOs that share
# longer lines of that map.
printf 'accession\tcoverage\nK08907\t1\nK08908\t2\nK08909\t3\nK08910\t4\nK16043\t5\nK16044\t6\nK22230\t7\nK00174\t8\nK00175\t9\nK01624\t10\nK01622\t11\nK16306\t12\n' > draw_kos_global_slices.reaction.txt
# A negative weight.
printf 'accession\tcoverage\tweight\nK00174\t1.5\t-1\nK00175\t3.0\t2\n' > draw_bad_negative_weight.reaction.txt
# A row without a gene ID.
printf 'accession\tgene_id\tcoverage\nK00844\tgene_1\t10\nK00845\t\t30\n' > draw_bad_blank_gene_id.reaction.txt
# Gene caller IDs of two assemblies, used as they are. Gene 0 has K01623 in S1 and K01624 in S2, so
# the ID names two genes. Gene 1 has K01622 in both samples.
printf 'accession\tgene_id\tsample\tcoverage\nK01623\t0\tS1\t1\nK01624\t0\tS2\t9\nK01622\t1\tS1\t2\nK01622\t1\tS2\t4\n' > draw_kos_gene_id_collision.reaction.txt
# Genes with rows for one or two KOs of element 18 of map 00010. The values are numbers, categorical
# values across samples, and clock times.
printf 'accession\tgene_id\tcoverage\nK01623\tgene_1\t2\nK01624\tgene_1\t5\nK01622\tgene_2\t4\n' > draw_kos_genes_coverage.reaction.txt
printf 'accession\tgene_id\tsample\tstatus\nK01623\tgene_1\tS1\tup\nK01624\tgene_1\tS1\tdown\nK01622\tgene_2\tS1\tup\nK01622\tgene_2\tS2\tup\nK01624\tgene_2\tS3\tdown\nK01622\tgene_2\tS3\tdown\nK01622\tgene_3\tS1\tdown\nK01624\tgene_3\tS1\tup\n' > draw_kos_samples_genes_status.reaction.txt
printf 'accession\tgene_id\tclocktime\nK01623\tgene_1\t23\nK01624\tgene_1\t3\nK01622\tgene_2\t12\n' > draw_kos_genes_clocktime.reaction.txt
# A value column and a border column, with an empty border cell, and a file with one more column.
printf 'accession\tcoverage\tstatus\nK00844\t10\tup\nK01810\t20\tdown\nK00850\t30\t\n' > draw_kos_coverage_status.reaction.txt
printf 'accession\tcoverage\tstatus\textra\nK00844\t10\tup\t1\n' > draw_bad_three_columns.reaction.txt
# A border column named twice, a border column named with '__', and a border column with no values.
printf 'accession\tstatus\tstatus\nK00844\tup\tdown\n' > draw_bad_repeated_border.reaction.txt
printf 'accession\t__status\nK00844\tup\n' > draw_bad_reserved_border.reaction.txt
printf 'accession\tcoverage\tstatus\nK00844\t10\t\nK01810\t20\t\n' > draw_bad_empty_border.reaction.txt
# 'NA' in a column of numeric border values, and two rows of one KO with only a border column.
printf 'accession\tcoverage\tlevel\nK00844\t10\t2.5\nK01810\t20\tNA\n' > draw_bad_na_border.reaction.txt
printf 'accession\tstatus\nK00844\tup\nK00844\tdown\n' > draw_bad_repeated_border_rows.reaction.txt
# Eight categorical border values, one more than the default border colors.
printf 'accession\tstatus\nK00844\tv1\nK00845\tv2\nK01810\tv3\nK00850\tv4\nK16370\tv5\nK01623\tv6\nK01624\tv7\nK01803\tv8\n' > draw_kos_border_class8.reaction.txt
# A stray quote in a categorical border value of a file whose fill values are numbers.
printf 'accession\tcoverage\tstatus\nK00844\t10\t"high\nK01810\t20\tlow\n' > draw_kos_border_stray_quote.reaction.txt
# K00844 has one gene, so the standard deviation of its genes gives it no fill value. Its border
# value is the mean of its genes. It is the only KO of the file on map 00010.
printf 'accession\tgene_id\tcoverage\tlevel\nK00844\tgene_1\t10\t4\n' > draw_kos_border_no_fill.reaction.txt
# K03850 is drawn as a line on the standard map 00510, and K00844 as a box on 00010.
printf 'accession\tstatus\nK03850\tup\nK00844\tdown\n' > draw_kos_border_lines.reaction.txt
# Twelve and eleven genes of K00844, which boxes 69 and 75 of map 00010 stand for, with border
# values. A ring of 12 slices is over the cap of 11, and a fill of 11 slices inside a ring is over
# the cap of 10.
awk 'BEGIN {OFS = "\t"; print "accession", "gene_id", "coverage", "status"; split("up down same", s, " "); for (i = 1; i <= 12; i++) print "K00844", "gene_" i, i, s[i % 3 + 1]}' > draw_kos_genes12_status.reaction.txt
awk 'NR <= 12' draw_kos_genes12_status.reaction.txt > draw_kos_genes11_status.reaction.txt
# The mock of weights and border values without its value column, and without its weights too.
cut -f 1-3,5,6 draw_kos_samples_weights_status.reaction.txt > draw_kos_samples_weights_border.reaction.txt
cut -f 1-3,6 draw_kos_samples_weights_status.reaction.txt > draw_kos_samples_genes_border.reaction.txt
# Weights of 0 with border values. K16306 of box 18 of map 00010 has no slice of the fill, but it
# has a border value. Both KOs of box 154 have a weight of 0.
printf 'accession\tcoverage\tweight\tstatus\nK00174\t1.5\t0\tup\nK00175\t3.0\t0\tdown\nK01624\t2.0\t3\tup\nK01622\t4.0\t1\tdown\nK16306\t6.0\t0\tsame\n' > draw_kos_weights_status.reaction.txt
# A border column alone, with an empty border cell. K00844, K01810, and K00850 are boxes of map
# 00010.
printf 'accession\tstatus\nK00844\tup\nK01810\tdown\nK00850\t\n' > draw_kos_status_empty_cell.reaction.txt
# A bar column beside a value column, with an empty bar cell. K01810 then has a fill and no bar value.
printf 'accession\tcoverage\tlevel\nK00844\t10\t2.5\nK01810\t20\t\n' > draw_kos_bar_empty_cell.reaction.txt
# Bar values of 0 only, 'NA' among bar values, a negative bar value, and a bar column with no values.
printf 'accession\tlevel\nK00844\t0\nK01810\t0\n' > draw_kos_bar_zero.reaction.txt
printf 'accession\tcoverage\tlevel\nK00844\t10\tNA\nK01810\t20\t3\n' > draw_bad_na_bar.reaction.txt
printf 'accession\tcoverage\tlevel\nK00844\t10\t-2\nK01810\t20\t3\n' > draw_bad_negative_bar.reaction.txt
printf 'accession\tcoverage\tlevel\nK00844\t10\t\nK01810\t20\t\n' > draw_bad_empty_bar.reaction.txt
# K03850 is drawn as a line on the standard map 00510, and K00844 as a box on 00010.
printf 'accession\tlevel\nK03850\t2\nK00844\t3\n' > draw_kos_bar_lines.reaction.txt

# A test that expects a refusal discards the output of the draw, so a fixture that was never written
# would make that test pass for the wrong reason. Fail loudly here instead.
for fixture in contigs-db-group-information.txt draw_kos_samples300.reaction.txt \
              draw_kos_samples60.reaction.txt draw-sample-group-information-60.txt \
              draw_kos_sparse300.reaction.txt draw_compounds_samples300.compound.txt \
              draw-sample-group-information-300.txt draw_kos_constant_coverage.reaction.txt \
              draw_kos_samples_negative.reaction.txt \
              draw-sample-colors-combo.txt draw-sample-colors-incomplete.txt draw-sample-colors-clashing.txt \
              draw_kos_awkward_sample_names.reaction.txt draw_kos_unusable_sample_names.reaction.txt \
              draw_bad_mixed_kr.reaction.txt draw_bad_repeated_rows.reaction.txt \
              draw_bad_two_values.reaction.txt draw_bad_compound_gene_id.compound.txt \
              draw_bad_repeated_compounds.compound.txt draw_kos_class.reaction.txt \
              draw_bad_repeated_weights.reaction.txt draw_bad_repeated_sample.reaction.txt \
              draw_bad_reserved_header.reaction.txt draw_bad_stray_quote_number.reaction.txt \
              draw_compounds_samples_status.compound.txt draw_kos_class19.reaction.txt \
              draw-status-colors-incomplete.txt draw-status-colors-shared.txt \
              draw_bad_blank_status.reaction.txt draw_kos_tie.reaction.txt \
              draw_kos_class18.reaction.txt draw_kos_samples_class19.reaction.txt \
              draw_bad_repeated_status.reaction.txt draw_kos_header_only_status.reaction.txt \
              draw_kos_gene_votes.reaction.txt draw_kos_stray_quote.reaction.txt \
              draw-stray-quote-colors.txt draw_kos_quoted_fields.reaction.txt \
              draw_kos_header_quote_status.reaction.txt \
              draw_kos_header_quote_coverage.reaction.txt \
              draw_kos_weights.reaction.txt draw_kos_global_slices.reaction.txt \
              draw_bad_negative_weight.reaction.txt draw_bad_blank_gene_id.reaction.txt \
              draw_kos_gene_id_collision.reaction.txt draw_kos_genes_coverage.reaction.txt \
              draw_kos_samples_genes_status.reaction.txt draw_kos_genes_clocktime.reaction.txt \
              draw_kos_coverage_status.reaction.txt draw_bad_three_columns.reaction.txt \
              draw_bad_repeated_border.reaction.txt draw_bad_reserved_border.reaction.txt \
              draw_bad_empty_border.reaction.txt draw_bad_na_border.reaction.txt \
              draw_bad_repeated_border_rows.reaction.txt draw_kos_border_class8.reaction.txt \
              draw_kos_border_stray_quote.reaction.txt draw_kos_border_no_fill.reaction.txt \
              draw_kos_border_lines.reaction.txt draw_kos_genes12_status.reaction.txt \
              draw_kos_genes11_status.reaction.txt draw_kos_samples_weights_border.reaction.txt \
              draw_kos_samples_genes_border.reaction.txt draw_kos_weights_status.reaction.txt \
              draw_kos_status_empty_cell.reaction.txt draw_kos_bar_empty_cell.reaction.txt \
              draw_kos_bar_zero.reaction.txt draw_bad_na_bar.reaction.txt \
              draw_bad_negative_bar.reaction.txt draw_bad_empty_bar.reaction.txt \
              draw_kos_bar_lines.reaction.txt
do
    if [ ! -s $fixture ]
    then
        echo "ERROR: the development fixture $fixture was not generated."
        exit 1
    fi
done

# Limits that no value reaches still set the ends of both scales. Each colorbar labels both limits
# with no '≤' or '≥' mark.
INFO "Testing limits that no value reaches, which still set the ends of both scales"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
                        --output-dir draw_txt_samples_kos_value_limits_unreached \
                        --reaction-gene-aggregation mean \
                        --reaction-sample-summary mean \
                        --reaction-value-limits 0 200 \
                        --reaction-individual-value-limits 0 200 \
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
# kept from wasting most of its colormap on a lopsided tail. The limit sets the top of the scale,
# and the center then widens the bottom to match. The scale comes out centered on zero, over a good
# deal less than the values span, with its top marked '>='.
INFO "Testing a centered scale whose longer side is trimmed back by a limit"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_log_ratio.reaction.txt \
                        --output-dir draw_txt_samples_kos_value_center_limits \
                        --reaction-gene-aggregation mean \
                        --reaction-accession-aggregation mean \
                        --reaction-sample-summary mean \
                        --reaction-colormap RdBu_r \
                        --reaction-value-center \
                        --reaction-individual-value-center \
                        --reaction-value-limits none 3 \
                        --reaction-individual-value-limits none 3 \
                        --pathway-numbers $pathway_numbers \
                        --draw-grid \
                        --no-progress

# Two limits set both ends of a scale, so a center is accepted only at their midpoint. These values
# run from -2 to 6, so the bottom limit is labeled with no mark and the top is marked '>='.
INFO "Testing a centered scale set at both ends by limits the same distance from the center"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_log_ratio.reaction.txt \
                        --output-dir draw_txt_samples_kos_value_center_both_limits \
                        --reaction-gene-aggregation mean \
                        --reaction-accession-aggregation mean \
                        --reaction-sample-summary mean \
                        --reaction-colormap RdBu_r \
                        --reaction-value-center \
                        --reaction-individual-value-center \
                        --reaction-value-limits -3 3 \
                        --reaction-individual-value-limits -3 3 \
                        --pathway-numbers $pathway_numbers \
                        --draw-grid \
                        --no-progress

# Centering only one of two scales that both color by value warns that the other is not centered.
# Both scales use one colormap here.
INFO "Testing that centering one of two scales that both color by value warns about the other"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_log_ratio.reaction.txt \
    --output-dir draw_txt_samples_kos_value_center_one \
    --reaction-gene-aggregation mean --reaction-sample-summary mean \
    --reaction-colormap RdBu_r --reaction-value-center \
    --pathway-numbers $pathway_numbers --draw-grid \
    --no-progress 2>&1 | tee draw_txt_value_center_one.log
if ! tr '\n' ' ' < draw_txt_value_center_one.log \
    | grep -q "is not centered"
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
# Also draws the global map, in the colors of the samples and their blends.
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
                        --reaction-membership-colors draw-sample-colors-combo.txt \
                        --output-dir draw_txt_samples_membership_colors_combo \
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
                        --reaction-membership-colors draw-sample-colors.txt \
                        --compound-membership-colors draw-sample-colors-combo.txt \
                        --output-dir draw_txt_two_layers_membership_colors \
                        --pathway-numbers $pathway_numbers \
                        --no-progress

# The same group ramps with an explicit span, which says how far from white each ramp starts and
# stops rather than which fraction of a colormap to use.
INFO "Testing the same group ramps with an explicit span"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
                        --groups-txt draw-sample-group-information.txt \
                        --group-threshold 0.5 \
                        --reaction-membership-colors draw-group-colors.txt \
                        --group-colormap tint 0.4 0.95 \
                        --group-reverse-overlay \
                        --output-dir draw_txt_groups_tint_ramps_span \
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
                        --reaction-membership-colors draw-group-colors.txt \
                        --group-colormap tint \
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
                        --reaction-value-limits 5 35 \
                        --reaction-individual-value-limits 5 35 \
                        --compound-value-limits 55 none \
                        --compound-individual-value-limits 20 40 \
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
                        --output-dir draw_txt_both_layers_individual_colormap \
                        --reaction-gene-aggregation mean \
                        --reaction-sample-summary std \
                        --compound-sample-summary std \
                        --reaction-colormap viridis \
                        --reaction-individual-colormap RdYlGn \
                        --compound-colormap cividis \
                        --compound-individual-colormap PuOr \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --no-progress

for colorbar in reactions reactions_samples compounds compounds_samples
do
    if [ ! -s draw_txt_both_layers_individual_colormap/colorbar_$colorbar.pdf ]
    then
        echo "ERROR: two layers with two colormaps each should have written four colorbars."
        exit 1
    fi
done

# An '_r' suffix reverses a colormap that anvi'o defines. It does the same for a Matplotlib
# colormap.
INFO "Testing the reversed anvi'o colormap 'clocktime_r'"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_clocktime.reaction.txt \
                        --output-dir draw_txt_samples_kos_clocktime_r \
                        --reaction-individual-colormap clocktime_r \
                        --reaction-individual-value-limits 0 24 \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --no-progress

# A circular mean of clock times that cancel out is undefined. The two samples of K01251 are exactly
# 12 hours apart. K01251 stands alone on the global map, so its element is left uncolored on the
# 'unified' map, and a warning names it. Limits of exactly 0 and the period are accepted.
INFO "Testing that a circular mean of clock times that cancel out is undefined"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_clocktime.reaction.txt \
    --output-dir draw_txt_samples_kos_circular_mean \
    --reaction-value-period 24 --reaction-sample-summary circular_mean \
    --reaction-value-limits 0 24 --pathway-numbers $global_pathway_numbers \
    --no-progress 2>&1 | tee draw_txt_circular_mean.log
for phrase in "A circular mean is undefined" "K01251"
do
    if ! tr '\n' ' ' < draw_txt_circular_mean.log | grep -q "$phrase"
    then
        echo "ERROR: a circular mean of clock times that cancel out did not warn ($phrase)."
        exit 1
    fi
done

# The angular deviation colors the 'unified' map by how far the samples' clock times are spread
# apart. It needs times in at least 3 samples. Elements in fewer samples are left uncolored, and a
# warning says so.
INFO "Testing the angular deviation of clock times on the 'unified' map"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_clocktime.reaction.txt \
    --output-dir draw_txt_samples_kos_angular_deviation \
    --reaction-value-period 24 --reaction-sample-summary angular_deviation \
    --pathway-numbers $pathway_numbers --draw-individual-files \
    --no-progress 2>&1 | tee draw_txt_angular_deviation.log
if ! tr '\n' ' ' < draw_txt_angular_deviation.log | grep -q "needs at least 3 values"
then
    echo "ERROR: an angular deviation of fewer than 3 clock times did not warn."
    exit 1
fi

# Each sample's map shows each element's clock time as its offset from the element's circular mean
# across the samples. The two samples of K01251 are exactly 12 hours apart, so they have no circular
# mean. K01251 stands alone on the global map, where it is left uncolored, and a warning names it.
INFO "Testing offsets of clock times from their circular mean"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_clocktime.reaction.txt \
    --output-dir draw_txt_samples_kos_circular_offsets \
    --reaction-value-period 24 --reaction-element-normalization difference_from_circular_mean \
    --pathway-numbers $global_pathway_numbers --draw-individual-files \
    --no-progress 2>&1 | tee draw_txt_circular_offsets.log
for phrase in "no circular mean to be offset from" "K01251"
do
    if ! tr '\n' ' ' < draw_txt_circular_offsets.log | grep -q "$phrase"
    then
        echo "ERROR: an offset from a circular mean that is undefined did not warn ($phrase)."
        exit 1
    fi
done

# The compound file has a value column and no 'sample' column, so the compound layer has no scale
# of its own over the samples: every individual map takes its compound colors from the scale of the
# 'unified' map. Skipping that map therefore has to leave its colorbar behind, that bar being the
# only key to the colors on the maps that are drawn.
INFO "Testing that skipping the 'unified' map keeps the colorbar of the scale its maps borrow"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
                        --compound-txt draw_compounds.compound.txt \
                        --output-dir draw_txt_skip_unified_borrowed_scale \
                        --pathway-numbers $pathway_numbers \
                        --compound-value-limits 1 6 \
                        --skip-unified-maps \
                        --draw-individual-files \
                        --no-progress

if [ ! -s draw_txt_skip_unified_borrowed_scale/colorbar_compounds.pdf ]
then
    echo "ERROR: skipping the 'unified' map should keep the colorbar its individual maps borrow."
    exit 1
fi
if [ -e draw_txt_skip_unified_borrowed_scale/colorbar_reactions.pdf ]
then
    echo "ERROR: the colorbar of the skipped 'unified' map's own scale should not be drawn."
    exit 1
fi

# An aggregation outside the validated names, resolved through pandas ('var' within a sample),
# together with one that maps disagreement among samples ('std' across them). The standard deviation
# is undefined for some map elements. These are left uncolored.
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

# The players of the same samples. Samples 'grid' and 'by_map' have no data on the map. The grid
# draws blank maps for them and deletes them, and the player draws its own, so it still has a frame
# for each of the five samples.
INFO "Writing players for samples whose names match the output directory's own subdirectories"
anvi-draw-kegg-pathways --reaction-txt draw_kos_awkward_sample_names.reaction.txt \
                        --categorize-files \
                        --output-dir draw_txt_awkward_sample_names_animated \
                        --pathway-numbers $pathway_numbers \
                        --draw-individual-files \
                        --draw-grid \
                        --collate-files-by-map \
                        --animate-files-by-map \
                        --no-progress

player=draw_txt_awkward_sample_names_animated/by_map/Metabolism/Carbohydrate_metabolism/ko00010.html
if [ "$(grep -o '"label":' $player | wc -l)" -ne 5 ]
then
    echo "ERROR: the player should hold a frame for each of the five samples."
    exit 1
fi
if [ "$(grep -o 'No data on this map' $player | wc -l)" -ne 2 ]
then
    echo "ERROR: the player should show the bare map for the two samples without data."
    exit 1
fi

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

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --draw-individual-files --animate-files-by-map \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_animate_no_collate.log 2>&1
then
    echo "ERROR: a player without gathering by map should have failed."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_animate_no_collate.log | grep -q "but that flag was not given"
then
    echo "ERROR: a player without gathering by map was refused for another reason."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --draw-individual-files SAMPLE_1 --collate-files-by-map --animate-files-by-map \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_animate_one_sample.log 2>&1
then
    echo "ERROR: a player of one sample's maps should have failed."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_animate_one_sample.log | grep -q "so each player would hold"
then
    echo "ERROR: a player of one sample's maps was refused for another reason."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos.reaction.txt \
    --reaction-labels none standard \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_labels_none_and_kind.log 2>&1
then
    echo "ERROR: labels given as 'none' and a kind of map should have failed."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_labels_none_and_kind.log | grep -q "'none' labels no map"
then
    echo "ERROR: labels given as 'none' and a kind of map were refused for another reason."
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

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --skip-unified-maps \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: skipping the 'unified' map with no other maps asked for should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos.reaction.txt \
    --skip-unified-maps --draw-individual-files \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: skipping the 'unified' map of a file with no 'sample' column should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --skip-unified-maps --draw-grid --reaction-sample-summary count \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a sample summary for a skipped 'unified' map should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --groups-txt draw-sample-group-information.txt --reaction-group-summary count \
    --skip-unified-maps --draw-grid \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a group summary for a skipped 'unified' map should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --groups-txt draw-sample-group-information.txt --group-threshold 0 \
    --skip-unified-maps --draw-grid \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a group threshold for a skipped 'unified' map should have failed."
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

# A column name that starts with '__' could clash with a column that anvi'o adds as it reads the
# file, such as '__accession'.
if anvi-draw-kegg-pathways --reaction-txt draw_bad_reserved_header.reaction.txt \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_reserved_header.log 2>&1
then
    echo "ERROR: a column named '__accession' should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_reserved_header.log \
    | grep -q "Names that start with '__' are kept for"
then
    echo "ERROR: a column named '__accession' was refused for another reason."
    exit 1
fi

# A stray quote stays in its field, so a number that holds one is refused. The quote does not join
# the rows up to the next quote into one field.
if anvi-draw-kegg-pathways --reaction-txt draw_bad_stray_quote_number.reaction.txt \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_stray_quote_number.log 2>&1
then
    echo "ERROR: a number with a stray quote should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_stray_quote_number.log | grep -q "must contain a finite number"
then
    echo "ERROR: a number with a stray quote was refused for another reason."
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
    --reaction-membership-colors draw-sample-colors-incomplete.txt \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a colors file leaving a sample of the run uncolored should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --reaction-membership-colors draw-sample-colors-clashing.txt \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a colors file whose combinations cannot be told apart should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --reaction-membership-colors draw-sample-colors.txt --reaction-colormap tab10 \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a membership colors file together with a colormap should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --reaction-membership-colors draw-sample-colors.txt --reaction-sample-summary count \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a membership colors file with a count summary should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos.reaction.txt \
    --reaction-membership-colors draw-sample-colors.txt \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a membership colors file for a layer with no 'sample' column should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --groups-txt draw-sample-group-information.txt --group-threshold 0.5 --draw-individual-files \
    --group-colormap tint \
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
    --reaction-membership-colors draw-group-colors.txt --group-colormap tint 0.0 1.0 \
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

# A circular mean needs to know the period after which the values repeat. A period in turn rules out
# every aggregation, summary of values, and normalization that treats the values as lying on a line.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_clocktime.reaction.txt \
    --reaction-sample-summary circular_mean \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a circular mean without '--reaction-value-period' should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_clocktime.reaction.txt \
    --reaction-value-period 24 --reaction-accession-aggregation sum \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: an aggregation other than a circular mean with a period should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_clocktime.reaction.txt \
    --reaction-value-period 24 --reaction-sample-summary mean \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a sample summary other than a circular mean with a period should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_clocktime.reaction.txt \
    --reaction-value-period 24 --reaction-element-normalization difference_from_mean \
    --draw-individual-files \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a normalization on a line, of values with a period, should have failed."
    exit 1
fi

# The offset from a circular mean needs a period. Its scale is fixed from -12 to 12, with 0 in the
# middle, so it takes no other limits and no center.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_clocktime.reaction.txt \
    --reaction-element-normalization difference_from_circular_mean --draw-individual-files \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: an offset from a circular mean without a period should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_clocktime.reaction.txt \
    --reaction-value-period 24 --reaction-element-normalization difference_from_circular_mean \
    --reaction-individual-value-limits 0 24 --draw-individual-files \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: limits other than -12 and 12 for offsets from a circular mean should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_clocktime.reaction.txt \
    --reaction-value-period 24 --reaction-element-normalization difference_from_circular_mean \
    --reaction-individual-value-center 0 --draw-individual-files \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a center for offsets from a circular mean should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_clocktime.reaction.txt \
    --reaction-value-period 0 \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a period that is not a positive number should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --reaction-value-period 24 \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a period for a file with no value column should have failed."
    exit 1
fi

# A period fixes each scale of values to run from 0 to the period, so it takes no other limits and
# no center.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_clocktime.reaction.txt \
    --reaction-value-period 24 --reaction-individual-value-limits 6 18 --draw-individual-files \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: limits other than 0 and the period, with a period, should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_clocktime.reaction.txt \
    --reaction-value-period 24 --reaction-sample-summary circular_mean \
    --reaction-value-center 12 \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a center with a period should have failed."
    exit 1
fi

# An angular deviation needs a period as well. Since angular deviation is a spread rather than a
# point in the period, only the summary that colors the 'unified' map takes it, and that map's scale
# is fixed, with a colormap that is not cyclic.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_clocktime.reaction.txt \
    --reaction-sample-summary angular_deviation \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: an angular deviation without '--reaction-value-period' should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_clocktime.reaction.txt \
    --reaction-value-period 24 --reaction-accession-aggregation angular_deviation \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: an angular deviation as an accession aggregation should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_clocktime.reaction.txt \
    --groups-txt draw-sample-group-information.txt --reaction-value-period 24 \
    --reaction-sample-summary angular_deviation --reaction-group-summary circular_mean \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: an angular deviation as the sample summary of a grouped run should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_clocktime.reaction.txt \
    --reaction-value-period 24 --reaction-sample-summary angular_deviation \
    --reaction-colormap clocktime \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a cyclic colormap for an angular deviation should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_clocktime.reaction.txt \
    --reaction-value-period 24 --reaction-sample-summary angular_deviation \
    --reaction-value-limits 0 24 \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: limits on the 'unified' scale of an angular deviation should have failed."
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
    --reaction-color "#FFFFFF" \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: another spelling of a reserved color should have failed but did not."
    exit 1
fi

# The refusal of an unrecognized name lists the colormaps that anvi'o defines.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-colormap not_a_colormap \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_bad_colormap.log 2>&1
then
    echo "ERROR: an unrecognized colormap name should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_bad_colormap.log | grep -q "'clocktime'"
then
    echo "ERROR: the refusal of an unrecognized colormap name did not list 'clocktime'."
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
    --reaction-individual-value-limits 5 40 \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: limits on the individual maps' scale where the file has no 'sample' column should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-individual-value-limits 5 40 --pathway-numbers $pathway_numbers \
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
    --reaction-individual-value-center \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a center on the individual maps' scale where the file has no 'sample' column should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-individual-value-center --pathway-numbers $pathway_numbers \
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

# Neither end of a centered scale can lie on its center, so a center on a limit is refused.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-value-center 5 --reaction-value-limits 5 none \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_center_on_limit.log 2>&1
then
    echo "ERROR: a center lying on a limit of the same scale should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_center_on_limit.log | grep -q "Neither of its ends"
then
    echo "ERROR: a center lying on a limit of the same scale was refused for another reason."
    exit 1
fi

# Two limits set both ends of a scale, which leaves centering nothing to widen. A center between
# them is refused unless it is their midpoint, and this is caught before any value is read.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-value-center 10 --reaction-value-limits 5 40 \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_center_off_midpoint.log 2>&1
then
    echo "ERROR: a center other than the midpoint of two limits should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_center_off_midpoint.log | grep -q "centered on their midpoint"
then
    echo "ERROR: a center other than the midpoint of two limits was refused for another reason."
    exit 1
fi

# A limit sets its end of a scale, and the center widens the shorter side afterwards. Centering can
# therefore ask to move a limit. That undoes what the limit was for, and is refused rather than
# carried out quietly. These values run from -2 to 6. The first limit cuts off values below -1. The
# second limit is reached by no value, and is refused all the same.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_log_ratio.reaction.txt \
    --reaction-gene-aggregation mean --reaction-accession-aggregation mean \
    --reaction-sample-summary mean --reaction-value-center --reaction-value-limits -1 none \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: centering a scale past a limit that values lie past should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_log_ratio.reaction.txt \
    --reaction-gene-aggregation mean --reaction-accession-aggregation mean \
    --reaction-sample-summary mean --reaction-value-center --reaction-value-limits -2.5 none \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_center_moves_limit.log 2>&1
then
    echo "ERROR: centering a scale past a limit that no value reaches should have failed."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_center_moves_limit.log | grep -q "would then have its minimum"
then
    echo "ERROR: centering past a limit that no value reaches was refused for another reason."
    exit 1
fi

# The colormap for the maps of the individual samples or groups is refused wherever there is no such
# scale to color: no file, no value column, no 'sample' column, or a scale that is not colored by
# value.
if anvi-draw-kegg-pathways --compound-txt draw_compounds.compound.txt \
    --reaction-individual-colormap viridis \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--reaction-individual-colormap' for a layer with no file should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --reaction-individual-colormap viridis --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--reaction-individual-colormap' for a layer with no value column should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-individual-colormap viridis \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--reaction-individual-colormap' where the file has no 'sample' column should have failed."
    exit 1
fi

# Grouped, the per-group maps are colored by value only where the sample summary pools values;
# summarized by presence they show within-group sample counts, which '--group-colormap' styles.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --groups-txt draw-sample-group-information.txt --group-threshold 0.5 \
    --reaction-individual-colormap viridis --pathway-numbers $pathway_numbers \
    --draw-individual-files \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--reaction-individual-colormap' for group maps colored by presence should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-individual-colormap viridis 0.9 0.1 --pathway-numbers $pathway_numbers \
    --draw-individual-files \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--reaction-individual-colormap' fractions the wrong way around should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-individual-colormap viridis 0.1 0.5 0.9 --pathway-numbers $pathway_numbers \
    --draw-individual-files \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--reaction-individual-colormap' given three numbers should have failed but did not."
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

# Nor are there values that could repeat after a period.
if anvi-draw-kegg-pathways --external-genomes external-genomes.txt \
    --reaction-value-period 24 \
    --output-dir draw_db_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a period on a database run should have failed but did not."
    exit 1
fi

# The same goes for a normalization, which rescales the values of a value column that only the input
# text files carry.
if anvi-draw-kegg-pathways --external-genomes external-genomes.txt --draw-grid \
    --reaction-element-normalization relative_to_mean \
    --output-dir draw_db_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a normalization on a database run should have failed but did not."
    exit 1
fi

# Nor is there one for a second colormap to color: a database run's individual maps show one
# database or genome each, in a single color rather than along a scale.
if anvi-draw-kegg-pathways --external-genomes external-genomes.txt --draw-individual-files \
    --reaction-individual-colormap viridis \
    --output-dir draw_db_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--reaction-individual-colormap' on a database run should have failed but did not."
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

# A reaction network JSON file is a single source drawn in one color, so it has no individual maps
# for a player to show.
if anvi-draw-kegg-pathways --reaction-network-json pan-network.json \
    --animate-files-by-map --pathway-numbers $pathway_numbers \
    --output-dir draw_db_bad --overwrite-output-destinations \
    --no-progress > draw_json_animate.log 2>&1
then
    echo "ERROR: a player of a reaction network JSON should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_json_animate.log \
    | grep -q "single source drawn in one color): '--animate-files-by-map'"
then
    echo "ERROR: a player of a reaction network JSON was refused for another reason."
    exit 1
fi

INFO "Testing that a normalization of element values is refused where it has nothing to act on"
# A normalization rescales only the maps of the individual samples or groups, so asking for one
# without asking for those maps would go nowhere.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-element-normalization relative_to_mean \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a normalization with no individual maps should have failed."
    exit 1
fi

# Without a 'sample' column there is a single set of values with nothing to rescale them against.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-element-normalization relative_to_mean \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a normalization without a 'sample' column should have failed."
    exit 1
fi

# Without a value column the layer is colored by presence, so there is nothing to rescale.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --reaction-element-normalization relative_to_mean --draw-individual-files \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a normalization without a value column should have failed."
    exit 1
fi

# A normalization can be undefined for every element it is asked about, which is not a refusal: the
# maps of the individual samples simply come out blank. The 'unified' map is drawn from the
# unnormalized values and is unaffected, so the warning has to name the normalization rather than
# the accessions or the aggregation, which are both fine here.
INFO "Testing the warning that a normalization was undefined for every element"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_negative.reaction.txt \
    --output-dir draw_txt_samples_kos_normalization_undefined \
    --reaction-element-normalization relative_to_mean \
    --pathway-numbers $pathway_numbers --draw-individual-files \
    --no-progress 2>&1 | tee draw_txt_normalization_undefined.log
if ! tr '\n' ' ' < draw_txt_normalization_undefined.log \
    | grep -q "the normalization is undefined for every map element"
then
    echo "ERROR: a normalization undefined everywhere did not warn about the normalization."
    exit 1
fi
if [ -s draw_txt_samples_kos_normalization_undefined/colorbar_reactions_samples.pdf ]
then
    echo "ERROR: a normalization undefined everywhere should have drawn no per-sample colorbar."
    exit 1
fi

# A pandas Series method that gives each value a new value of its own works as a normalization too,
# and takes the label anvi'o composes from the method's own name. This is the incidental way to ask
# for one, the named normalizations being the way it is meant to be used, so it is checked here
# rather than in the packaged suite.
INFO "Testing a pandas transform used in place of a named normalization"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --output-dir draw_txt_samples_kos_element_normalization_transform \
    --reaction-element-normalization abs \
    --pathway-numbers $pathway_numbers --draw-individual-files \
    --no-progress
if [ ! -s draw_txt_samples_kos_element_normalization_transform/colorbar_reactions_samples.pdf ]
then
    echo "ERROR: a pandas transform used as a normalization should have written a colorbar."
    exit 1
fi

INFO "Testing the names a normalization cannot be given"
# Each of these is a name that does the wrong thing with a set of values, and each is caught by a
# different probe. A reduction such as 'mean' summarizes samples rather than rescaling each of them,
# which is what '--reaction-sample-summary' is for. A name that reorders the samples would give a
# sample another sample's value. A name that compares one sample with its neighbors depends on their
# order. That order comes only from how the samples happen to be named. 'duplicated' answers for one
# of two equal samples according to which came first, which values that all differ never reveal.
# 'squeeze' returns a bare number for the single value that an element found in one sample alone
# gives it, which is commonplace and would otherwise fail in the middle of drawing. 'convert_dtypes'
# answers in a nullable dtype, whose missing value no array of plain floats holds.
for normalization in not_a_pandas_method mean idxmax sort_values cumsum pct_change duplicated \
                     squeeze convert_dtypes
do
    if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
        --reaction-element-normalization $normalization --draw-individual-files \
        --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
    then
        echo "ERROR: '$normalization' should have failed as a normalization but did not."
        exit 1
    fi
done

# A colorbar labeled with nothing says nothing about what it shows.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-element-normalization relative_to_mean "   " --draw-individual-files \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a blank colorbar label should have failed but did not."
    exit 1
fi

# A normalization whose neutral value is zero centers the rescaled scale there, which means the run
# holds a center that nobody asked for. Wherever that center is reported, the message has to name
# the normalization that supplied it rather than the option for giving one by hand: naming
# '--reaction-individual-value-center' would send the reader to an option they never used. It is
# reported in four places. The first is a per-group scale that is not colored by value at all. The
# second is a scale whose limits leave no room for the center. The third is a scale whose two limits
# do not have the center as their midpoint. The fourth is a limit that centering would move.
INFO "Testing that a center supplied by a normalization is reported against the normalization"
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --groups-txt draw-sample-group-information.txt --group-threshold 0.5 --draw-grid \
    --reaction-element-normalization relative_to_mean \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_normalization_grouped.log 2>&1
then
    echo "ERROR: a normalization over groups summarized by presence should have failed."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_normalization_grouped.log \
    | grep -q "reaction-element-normalization"
then
    echo "ERROR: the refusal over presence-summarized groups did not name the normalization."
    exit 1
fi
if tr '\n' ' ' < draw_txt_normalization_grouped.log \
    | grep -q "reaction-individual-value-center"
then
    echo "ERROR: the refusal over presence-summarized groups blamed a center nobody gave."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-element-normalization relative_to_mean \
    --reaction-individual-value-limits 1 5 --draw-individual-files \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_normalization_center.log 2>&1
then
    echo "ERROR: a centered normalization under limits excluding zero should have failed."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_normalization_center.log \
    | grep -q "reaction-element-normalization"
then
    echo "ERROR: a center from a normalization was not reported against the normalization."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-element-normalization relative_to_mean \
    --reaction-individual-value-limits -0.5 1 --draw-individual-files \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_normalization_midpoint.log 2>&1
then
    echo "ERROR: a centered normalization under limits not centered on zero should have failed."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_normalization_midpoint.log | grep -q "normalization' centers"
then
    echo "ERROR: a center from a normalization off the midpoint of two limits was misreported."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-element-normalization relative_to_mean \
    --reaction-individual-value-limits -0.5 none --draw-individual-files \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_normalization_moved_limit.log 2>&1
then
    echo "ERROR: centering a normalization's scale past one of its limits should have failed."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_normalization_moved_limit.log | grep -q "normalization' centers"
then
    echo "ERROR: a center from a normalization that would move a limit was misreported."
    exit 1
fi
if tr '\n' ' ' < draw_txt_normalization_moved_limit.log | grep -q "value-center' centers"
then
    echo "ERROR: the refusal of a center that would move a limit blamed a center nobody gave."
    exit 1
fi

# Categorical values on a global map. A compound there takes the value that the accession rule
# chooses among the reactions touching it. The run reports how many choices the rule made.
INFO "Testing categorical values on a global map, where compounds take a value from their reactions"
anvi-draw-kegg-pathways --reaction-txt draw_kos_class.reaction.txt \
    --reaction-value-type categorical \
    --output-dir draw_txt_kos_categorical_global \
    --pathway-numbers $global_pathway_numbers \
    --no-progress 2>&1 | tee draw_txt_categorical_global.log
if ! tr '\n' ' ' < draw_txt_categorical_global.log \
    | grep -q "Choices of a reaction value among reactions touching a compound"
then
    echo "ERROR: a run of categorical values on a global map did not report choices for compounds."
    exit 1
fi

INFO "Testing categorical values of compounds across samples"
anvi-draw-kegg-pathways --compound-txt draw_compounds_samples_status.compound.txt \
    --compound-value-type categorical \
    --output-dir draw_txt_compounds_categorical \
    --pathway-numbers $pathway_numbers --draw-grid \
    --no-progress

# The compound file has the values 'up' and 'down'. The colors file also names 'same', which the run
# reports and ignores.
INFO "Testing a colors file that names a categorical value the compound file lacks"
anvi-draw-kegg-pathways --compound-txt draw_compounds_samples_status.compound.txt \
    --compound-value-type categorical \
    --compound-value-colors draw-status-colors.txt \
    --output-dir draw_txt_compounds_categorical_colors \
    --pathway-numbers $pathway_numbers --draw-grid \
    --no-progress 2>&1 | tee draw_txt_compounds_categorical_colors.log
if ! tr '\n' ' ' < draw_txt_compounds_categorical_colors.log \
    | grep -q "The value is 'same'. It is ignored."
then
    echo "ERROR: a value in a colors file that the compound file lacks was not reported."
    exit 1
fi

# Each value takes one of the default colors.
INFO "Testing as many categorical values as anvi'o has default colors for"
anvi-draw-kegg-pathways --reaction-txt draw_kos_class18.reaction.txt \
    --reaction-value-type categorical \
    --output-dir draw_txt_kos_categorical_default_colors \
    --pathway-numbers $pathway_numbers \
    --no-progress

# Where map elements overlap, the first value is drawn on top. '--reaction-reverse-overlay' draws
# the last value on top instead.
INFO "Testing categorical values with the last value drawn on top"
anvi-draw-kegg-pathways --reaction-txt draw_kos_class.reaction.txt \
    --reaction-value-type categorical \
    --reaction-reverse-overlay \
    --output-dir draw_txt_kos_categorical_reverse_overlay \
    --pathway-numbers $pathway_numbers \
    --no-progress

# A stray quote is part of a categorical value, as it is in a colors file, so it does not join the
# rows up to the next quote into one value. Joined rows would make a value that the colors file
# lacks, and the run would fail.
INFO "Testing that a stray quote is part of a categorical value"
anvi-draw-kegg-pathways --reaction-txt draw_kos_stray_quote.reaction.txt \
    --reaction-value-type categorical \
    --reaction-value-colors draw-stray-quote-colors.txt \
    --output-dir draw_txt_kos_categorical_quotes \
    --pathway-numbers $pathway_numbers \
    --no-progress

# A pair of quotes enclosing a whole field is dropped, so a file written with every field quoted
# still has an 'accession' column.
INFO "Testing a file of categorical values with every field in quotes"
anvi-draw-kegg-pathways --reaction-txt draw_kos_quoted_fields.reaction.txt \
    --reaction-value-type categorical \
    --output-dir draw_txt_kos_categorical_quoted_fields \
    --pathway-numbers $pathway_numbers \
    --no-progress

# A stray quote in a column name is part of the name, so it does not join the rows up to the next
# quote into the header line. This holds for a column of categorical values and of numbers.
INFO "Testing a stray quote in the name of a column of categorical values"
anvi-draw-kegg-pathways --reaction-txt draw_kos_header_quote_status.reaction.txt \
    --reaction-value-type categorical \
    --output-dir draw_txt_kos_categorical_header_quote \
    --pathway-numbers $pathway_numbers \
    --no-progress

INFO "Testing a stray quote in the name of a column of numbers"
anvi-draw-kegg-pathways --reaction-txt draw_kos_header_quote_coverage.reaction.txt \
    --output-dir draw_txt_kos_header_quote \
    --pathway-numbers $pathway_numbers \
    --no-progress

# A file of categorical values with a header and no rows has nothing to color, and the run says so.
INFO "Testing a file of categorical values with a header and no rows"
anvi-draw-kegg-pathways --reaction-txt draw_kos_header_only_status.reaction.txt \
    --reaction-value-type categorical \
    --output-dir draw_txt_categorical_header_only \
    --pathway-numbers $pathway_numbers \
    --no-progress 2>&1 | tee draw_txt_categorical_header_only.log
if ! tr '\n' ' ' < draw_txt_categorical_header_only.log | grep -q "has a header but no rows"
then
    echo "ERROR: a file of categorical values with a header and no rows was not reported."
    exit 1
fi

# 'plurality' gives no value where values tie, and 'majority' none where no value has more than
# half of the votes. Elements given no value are left uncolored, and a warning says how many.
INFO "Testing categorical rules that can give no value"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_status.reaction.txt \
    --reaction-value-type categorical \
    --reaction-accession-aggregation plurality \
    --reaction-sample-summary majority \
    --output-dir draw_txt_categorical_no_value \
    --pathway-numbers $pathway_numbers --draw-grid \
    --no-progress 2>&1 | tee draw_txt_categorical_no_value.log
if ! tr '\n' ' ' < draw_txt_categorical_no_value.log | grep -q "plurality' gave no value to"
then
    echo "ERROR: elements given no value by 'plurality' were not reported."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_categorical_no_value.log | grep -q "majority' gave no value to"
then
    echo "ERROR: elements given no value by 'majority' were not reported."
    exit 1
fi

# A map element gets one vote from each KO, not from each gene. The KO with two genes valued 'up'
# and the KO with one gene valued 'down' each give one vote, so 'plurality' gives the element no
# value.
INFO "Testing that a map element gets one vote from each KO, not from each gene"
anvi-draw-kegg-pathways --reaction-txt draw_kos_gene_votes.reaction.txt \
    --reaction-value-type categorical \
    --reaction-accession-aggregation plurality \
    --output-dir draw_txt_categorical_ko_votes \
    --pathway-numbers 00010 \
    --no-progress 2>&1 | tee draw_txt_categorical_ko_votes.log
if ! tr '\n' ' ' < draw_txt_categorical_ko_votes.log \
    | grep -q "plurality' gave no value to 1 map element"
then
    echo "ERROR: a map element did not get one vote from each KO."
    exit 1
fi

# A map where a rule gives no value to any element has nothing to color, so it is not drawn, as a
# map lacking data is not.
INFO "Testing that a map with no categorical value to color is not drawn"
anvi-draw-kegg-pathways --reaction-txt draw_kos_tie.reaction.txt \
    --reaction-value-type categorical \
    --reaction-sample-summary plurality \
    --output-dir draw_txt_categorical_tie \
    --pathway-numbers 00010 \
    --no-progress
if [ -e draw_txt_categorical_tie/unified/ko00010.pdf ]
then
    echo "ERROR: a map with no categorical value to color should not have been drawn."
    exit 1
fi

# A map grid is drawn even where the 'unified' map is not. The grid has a blank panel in its place.
# The blank map is drawn only for the grid, so it is not kept.
INFO "Testing that a map grid is drawn where the 'unified' map has no categorical value to color"
anvi-draw-kegg-pathways --reaction-txt draw_kos_tie.reaction.txt \
    --reaction-value-type categorical \
    --reaction-sample-summary plurality \
    --output-dir draw_txt_categorical_tie_grid \
    --pathway-numbers 00010 --draw-grid \
    --no-progress
if [ ! -s draw_txt_categorical_tie_grid/grid/ko00010.pdf ]
then
    echo "ERROR: a map grid should have been drawn with a blank panel for the 'unified' map."
    exit 1
fi
if [ -e draw_txt_categorical_tie_grid/unified/ko00010.pdf ]
then
    echo "ERROR: the blank 'unified' map drawn for a map grid should have been removed."
    exit 1
fi

# Grouped, each group's map counts the group's samples unless a rule chooses among them. The
# values then color nothing, which is worth a warning. No legend of the values is written.
INFO "Testing grouped categorical values without a rule among the samples of each group"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_status.reaction.txt \
    --reaction-value-type categorical \
    --groups-txt draw-sample-group-information.txt --group-threshold 0 \
    --output-dir draw_txt_groups_categorical_counts \
    --pathway-numbers $pathway_numbers --draw-grid \
    --no-progress 2>&1 | tee draw_txt_groups_categorical_counts.log
if ! tr '\n' ' ' < draw_txt_groups_categorical_counts.log \
    | grep -q "With groups of samples, each group's map counts"
then
    echo "ERROR: grouped categorical values coloring nothing did not warn."
    exit 1
fi
if [ -e draw_txt_groups_categorical_counts/colorbar_reactions_groups.pdf ]
then
    echo "ERROR: grouped categorical values without a rule should have written no legend."
    exit 1
fi

# Without a rule among each group's samples, no map shows the values. They need no colors, so there
# can be more of them than anvi'o has default colors for.
INFO "Testing grouped categorical values without a rule, more of them than there are default colors"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_class19.reaction.txt \
    --reaction-value-type categorical \
    --groups-txt draw-sample-group-information.txt --group-threshold 0 \
    --output-dir draw_txt_groups_categorical_counts_many \
    --pathway-numbers $pathway_numbers --draw-grid \
    --no-progress

# The choices across the groups are reported with the others.
INFO "Testing that choices of a categorical value across groups are reported"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_status.reaction.txt \
    --reaction-value-type categorical \
    --groups-txt draw-sample-group-information.txt \
    --reaction-sample-summary plurality_order --reaction-group-summary plurality_order \
    --output-dir draw_txt_groups_categorical_report \
    --pathway-numbers $pathway_numbers \
    --no-progress 2>&1 | tee draw_txt_groups_categorical_report.log
if ! tr '\n' ' ' < draw_txt_groups_categorical_report.log | grep -q "value among groups"
then
    echo "ERROR: choices of a categorical value across groups were not reported."
    exit 1
fi

# The 'unified' map of a layer with samples shows presence by default. Drawing it alone colors
# nothing by the values, which is worth a warning, whether the values are names or numbers.
INFO "Testing the warning that a value column colors nothing where only the 'unified' map is drawn"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_status.reaction.txt \
    --reaction-value-type categorical \
    --output-dir draw_txt_categorical_unified_only \
    --pathway-numbers $pathway_numbers \
    --no-progress 2>&1 | tee draw_txt_categorical_unified_only.log
if ! tr '\n' ' ' < draw_txt_categorical_unified_only.log \
    | grep -q "nothing on the maps is colored by it"
then
    echo "ERROR: categorical values coloring nothing did not warn."
    exit 1
fi
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --output-dir draw_txt_values_unified_only \
    --pathway-numbers $pathway_numbers \
    --no-progress 2>&1 | tee draw_txt_values_unified_only.log
if ! tr '\n' ' ' < draw_txt_values_unified_only.log | grep -q "nothing on the maps is colored by it"
then
    echo "ERROR: numeric values coloring nothing did not warn."
    exit 1
fi

# A column of names read as numbers is refused, and the refusal says how to read it as names.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_status.reaction.txt \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_names_as_numbers.log 2>&1
then
    echo "ERROR: a column of names read as numbers should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_names_as_numbers.log | grep -q "If the column holds categorical values"
then
    echo "ERROR: refusing a column of names read as numbers did not say how to read it as names."
    exit 1
fi

# Categorical values are refused wherever an option expects numbers, or expects something they do
# not have (each command below is expected to fail).
INFO "Testing that categorical values are refused where they cannot apply"
# Categorical values take colors of their own, not a scale of numbers.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_status.reaction.txt \
    --reaction-value-type categorical --reaction-value-limits 0 1 \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: value limits on categorical values should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_status.reaction.txt \
    --reaction-value-type categorical --reaction-value-period 24 \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a period for categorical values should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_status.reaction.txt \
    --reaction-value-type categorical --reaction-element-normalization z_score --draw-grid \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a normalization of categorical values should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_status.reaction.txt \
    --reaction-value-type categorical --reaction-individual-colormap viridis --draw-grid \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: '--reaction-individual-colormap' for categorical values should have failed but did not."
    exit 1
fi

# Categorical values cannot be averaged, across samples or within one.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_status.reaction.txt \
    --reaction-value-type categorical --reaction-sample-summary mean \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: averaging categorical values across samples should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_status.reaction.txt \
    --reaction-value-type categorical --reaction-accession-aggregation mean \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: averaging categorical values within an element should have failed but did not."
    exit 1
fi

# With groups, a rule across the groups needs a rule among each group's samples. Without one,
# the groups have no values to choose among.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_status.reaction.txt \
    --reaction-value-type categorical \
    --groups-txt draw-sample-group-information.txt --group-threshold 0 \
    --reaction-group-summary plurality_order \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a rule across groups without a rule among their samples should have failed."
    exit 1
fi

# A rule across the groups leaves no group presence for a threshold to decide.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_status.reaction.txt \
    --reaction-value-type categorical \
    --groups-txt draw-sample-group-information.txt --group-threshold 0 \
    --reaction-sample-summary plurality_order --reaction-group-summary plurality_order \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a group threshold with a rule across the groups should have failed."
    exit 1
fi

# Grouped, the values reach a map only through a rule among each group's samples, so a colors
# file of values without one would color nothing.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_status.reaction.txt \
    --reaction-value-type categorical --reaction-value-colors draw-status-colors.txt \
    --groups-txt draw-sample-group-information.txt --group-threshold 0 \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: value colors for grouped samples without a rule should have failed."
    exit 1
fi

# A seed with no rule that chooses at random would change nothing. An error indents its wrapped
# lines, so the refusal is matched by a phrase from its first line. The test of a random rule for
# numbers below checks that its refusal lacks this phrase.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_status.reaction.txt \
    --reaction-value-type categorical --seed 1 \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_seed_without_random_rule.log 2>&1
then
    echo "ERROR: a seed without a random rule should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_seed_without_random_rule.log | grep -q "'--seed' sets the seed"
then
    echo "ERROR: the refusal of a seed without a random rule did not name the seed."
    exit 1
fi

# Categorical values take colors of their own, not one presence color.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_status.reaction.txt \
    --reaction-value-type categorical --reaction-color \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a static color for categorical values should have failed but did not."
    exit 1
fi

# Every value needs a color, and no two values can share one.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_status.reaction.txt \
    --reaction-value-type categorical --reaction-value-colors draw-status-colors-incomplete.txt \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a colors file missing a value should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_status.reaction.txt \
    --reaction-value-type categorical --reaction-value-colors draw-status-colors-shared.txt \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a colors file giving two values one color should have failed but did not."
    exit 1
fi

# A layer without samples shows its values on its one map, so no colormap colors presence.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_class.reaction.txt \
    --reaction-value-type categorical --reaction-colormap viridis \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a colormap for categorical values without samples should have failed."
    exit 1
fi

# More values than anvi'o has default colors for need a colors file.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_class19.reaction.txt \
    --reaction-value-type categorical \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: nineteen categorical values without a colors file should have failed."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_bad_blank_status.reaction.txt \
    --reaction-value-type categorical \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a row without a categorical value should have failed but did not."
    exit 1
fi

# An accession and gene can have only one categorical value.
if anvi-draw-kegg-pathways --reaction-txt draw_bad_repeated_status.reaction.txt \
    --reaction-value-type categorical \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: two categorical values for one accession and gene should have failed but did not."
    exit 1
fi

# Pandas would read a second 'sample' column as 'sample.1' and take it for the value column.
if anvi-draw-kegg-pathways --reaction-txt draw_bad_repeated_sample.reaction.txt \
    --reaction-value-type categorical \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_repeated_sample.log 2>&1
then
    echo "ERROR: two columns named 'sample' should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_repeated_sample.log | grep -q "A name may appear only once"
then
    echo "ERROR: two columns named 'sample' were refused for another reason."
    exit 1
fi

# A categorical value type needs a value column to read.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --reaction-value-type categorical \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a categorical value type for a file without values should have failed."
    exit 1
fi

# Numbers are reduced with an aggregation. A rule chooses only among categorical values.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-sample-summary plurality_order \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a categorical rule for numeric values across samples should have failed but did not."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-accession-aggregation order \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a categorical rule for numeric values within an element should have failed but did not."
    exit 1
fi

# Colors per value need categorical values.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-value-colors draw-status-colors.txt \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: colors per value for numeric values should have failed but did not."
    exit 1
fi

# A random rule given to numbers is refused for the rule. The seed given with it is not blamed.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-accession-aggregation random --seed 5 \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_random_rule_numbers.log 2>&1
then
    echo "ERROR: a random rule for numeric values should have failed but did not."
    exit 1
fi
if tr '\n' ' ' < draw_txt_random_rule_numbers.log | grep -q "'--seed' sets the seed"
then
    echo "ERROR: the refusal of a random rule for numeric values blamed the seed."
    exit 1
fi

# Only a text file has a value column, so a database run has no values to read as names.
if anvi-draw-kegg-pathways --external-genomes external-genomes.txt \
    --reaction-value-type categorical \
    --output-dir draw_db_bad --overwrite-output-destinations --no-progress > /dev/null 2>&1
then
    echo "ERROR: a value type on a database run should have failed but did not."
    exit 1
fi

# A column of weights gives the slices their shares, and the key table lists them. Element 154 of
# map 00010 stands for K00174 and K00175. Their weights are 90 and 10 in SAMPLE_1, where K00174 has
# two genes of 60 and 30. They are 1 and 9 in SAMPLE_2, and SAMPLE_3 has K00174 alone. The shares of
# the two KOs are therefore 0.9 and 0.1 in SAMPLE_1, 0.1 and 0.9 in SAMPLE_2, and 1 and 0 in
# SAMPLE_3. Group G1 holds SAMPLE_1 alone and takes its shares. Group G2 averages the shares of
# SAMPLE_2 and SAMPLE_3 to 0.55 and 0.45. The 'unified' map averages the groups to 0.725 and 0.275.
INFO "Testing the shares of weighted slices averaged within groups and across them"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_weights.reaction.txt \
                        --groups-txt draw-sample-group-information.txt \
                        --reaction-sample-summary mean \
                        --reaction-group-summary mean \
                        --reaction-slice-by accession \
                        --reaction-slice-weights coverage \
                        --output-dir draw_txt_slices_share_mean \
                        --pathway-numbers 00010 \
                        --draw-individual-files \
                        --no-progress

slice_key=draw_txt_slices_share_mean/slice_key_reactions.txt
if [ "$(head -n 1 $slice_key | tr '\t' ' ')" \
     != "map drawn_map entry_id x y label slice accession gene_id value share color drawn_whole" ]
then
    echo "ERROR: the key table of slices does not have the expected header."
    exit 1
fi
# A slice of a KO names no gene, so its 'gene_id' is empty, though the file has a 'gene_id' column.
if awk -F'\t' 'NR > 1 && $9 != "" {found = 1} END {exit !found}' $slice_key
then
    echo "ERROR: a slice of a KO named a gene in the key table."
    exit 1
fi
shares=$(awk -F'\t' '$3 == "154" {printf "%s %s %.3f\n", $2, $8, $11}' $slice_key)
expected_shares=$(printf '%s\n' "unified K00174 0.725" "unified K00175 0.275" \
                  "G1 K00174 0.900" "G1 K00175 0.100" "G2 K00174 0.550" "G2 K00175 0.450")
if [ "$shares" != "$expected_shares" ]
then
    echo "ERROR: the shares of weighted slices were not averaged within groups and across them."
    exit 1
fi
# Slices follow the order of the KOs in the KGML file, where K00172 comes before K00171.
if [ "$(awk -F'\t' '$2 == "G2" && $3 == "133" {print $8}' $slice_key | tr '\n' ' ')" \
     != "K00169 K00170 K00172 K00171 " ]
then
    echo "ERROR: the slices of an element did not follow the order of its KOs in the KGML file."
    exit 1
fi

# 'pooled_sum' adds up the weights of each KO across the samples and divides each sum by the total
# of the element. On the 'unified' map, K00174 has 96 of the 115 units of weight of element 154, and
# K00175 has 19.
INFO "Testing the shares of weighted slices pooled by the sum of their weights"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_weights.reaction.txt \
                        --reaction-sample-summary mean \
                        --reaction-slice-by accession \
                        --reaction-slice-weights coverage \
                        --reaction-slice-share-summary pooled_sum \
                        --output-dir draw_txt_slices_pooled_sum \
                        --pathway-numbers 00010 \
                        --no-progress

shares=$(awk -F'\t' '$3 == "154" {printf "%s %s %.3f\n", $2, $8, $11}' \
         draw_txt_slices_pooled_sum/slice_key_reactions.txt)
if [ "$shares" != "$(printf '%s\n' "unified K00174 0.835" "unified K00175 0.165")" ]
then
    echo "ERROR: the shares of weighted slices were not pooled by the sum of their weights."
    exit 1
fi

# The weights of the genes of a KO are added up by default. With 'mean', K00174 has a weight of 45
# in SAMPLE_1, the mean of its two genes, and K00175 has a weight of 10.
INFO "Testing the weights of the genes of a KO reduced by their mean"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_weights.reaction.txt \
                        --reaction-slice-by accession \
                        --reaction-slice-weights coverage \
                        --reaction-slice-weight-aggregation mean \
                        --output-dir draw_txt_slices_weight_mean \
                        --pathway-numbers 00010 \
                        --draw-individual-files \
                        --no-progress

shares=$(awk -F'\t' '$2 == "SAMPLE_1" && $3 == "154" {printf "%s %.3f\n", $8, $11}' \
         draw_txt_slices_weight_mean/slice_key_reactions.txt)
if [ "$shares" != "$(printf '%s\n' "K00174 0.818" "K00175 0.182")" ]
then
    echo "ERROR: the weights of the genes of a KO were not reduced by their mean."
    exit 1
fi

# A KO with a weight of 0 gets no slice. Element 18 of map 00010 is sliced between K01624 and K01622
# in the shares 0.75 and 0.25. Both KOs of element 154 have a weight of 0, so the element is left
# uncolored. It has no rows in the key table, and a warning counts it.
INFO "Testing that a KO with a weight of 0 gets no slice"
anvi-draw-kegg-pathways --reaction-txt draw_kos_weights.reaction.txt \
    --reaction-slice-by accession --reaction-slice-weights weight \
    --output-dir draw_txt_slices_zero_weights --pathway-numbers 00010 \
    --no-progress 2>&1 | tee draw_txt_slices_zero_weights.log
if ! tr '\n' ' ' < draw_txt_slices_zero_weights.log | grep -q "were left uncolored, since every"
then
    echo "ERROR: an element whose KOs all have a weight of 0 was not reported as uncolored."
    exit 1
fi
shares=$(awk -F'\t' '$3 == "18" || $3 == "154" {printf "%s %s %.3f\n", $3, $8, $11}' \
         draw_txt_slices_zero_weights/slice_key_reactions.txt)
if [ "$shares" != "$(printf '%s\n' "18 K01624 0.750" "18 K01622 0.250")" ]
then
    echo "ERROR: a KO with a weight of 0 was given a slice."
    exit 1
fi

# Lines on the global map keep the widths of the KGML file. Entry 6838 is a line 20 points long and
# 6 points wide. A quarter of it is shorter than it is wide, so its four KOs with values do not fit
# as slices. The line is drawn whole, and a warning says so. The key table lists the slices the line
# would have had, with the reason. Other KOs share longer lines of the map, which are sliced.
INFO "Testing slices on the global map, where a line too short for its slices is drawn whole"
anvi-draw-kegg-pathways --reaction-txt draw_kos_global_slices.reaction.txt \
    --reaction-slice-by accession \
    --output-dir draw_txt_slices_global --pathway-numbers $global_pathway_numbers \
    --no-progress 2>&1 | tee draw_txt_slices_global.log
if ! tr '\n' ' ' < draw_txt_slices_global.log | grep -q "had more slices than fit"
then
    echo "ERROR: a line too short for its slices was not reported as drawn whole."
    exit 1
fi
slice_key=draw_txt_slices_global/slice_key_reactions.txt
if [ "$(awk -F'\t' '$3 == "6838" && $7 == "" && $13 == "line too short" {n++} END {print n + 0}' \
       $slice_key)" != 4 ]
then
    echo "ERROR: the key table did not list the four slices of a line drawn whole."
    exit 1
fi
# Entry 1327 is a line 32 points long. Compound circles cover 19 points of its ends, and 13 points
# stay visible. A third of its whole length would fit as a slice, but a third of its visible part is
# shorter than the line is wide, so the line is drawn whole.
if [ "$(awk -F'\t' '$1 == "01100" && $3 == "1327" && $13 == "line too short" {n++}
                    END {print n + 0}' $slice_key)" != 3 ]
then
    echo "ERROR: a line whose visible part is too short for its slices was not drawn whole."
    exit 1
fi
if ! awk -F'\t' '$1 == "01100" && $7 == "2" {found = 1} END {exit !found}' $slice_key
then
    echo "ERROR: no line of the global map was sliced."
    exit 1
fi

# A summary can be undefined for a slice, as the standard deviation is for a KO found in one sample.
# Such a slice is not drawn, and the warning about undefined summaries counts it.
INFO "Testing the warning that a summary was undefined for slices"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-slice-by accession --reaction-sample-summary std \
    --output-dir draw_txt_slices_std --pathway-numbers $pathway_numbers \
    --no-progress 2>&1 | tee draw_txt_slices_std.log
if ! tr '\n' ' ' < draw_txt_slices_std.log | grep -q "These slices are not drawn on the"
then
    echo "ERROR: slices with an undefined summary were not reported."
    exit 1
fi

# Element 18 of map 00010 stands for K01623, K01624, K01622, K11645, K16305, and K16306. Its slices
# of genes are ordered by the first of these KOs that each gene has a row for, and then by gene ID,
# with digits compared as numbers. gene_2 therefore comes before gene_10, and gene_10 before gene_1.
# gene_5 has rows for K01624 and K01622, so it is one slice that lists both KOs. In SAMPLE_1, these
# rows have values of 2 and 6 and weights of 30 and 10. The slice takes their means, a value of 4
# and a weight of 20. The slices of the element there weigh 10, 20, 20, and 50.
INFO "Testing slices of genes, ordered by KO and then by gene ID"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_genes.reaction.txt \
    --reaction-sample-summary mean --reaction-slice-by gene --reaction-slice-weights coverage \
    --output-dir draw_txt_slices_genes --pathway-numbers 00010 --draw-individual-files \
    --no-progress 2>&1 | tee draw_txt_slices_genes.log
slice_key=draw_txt_slices_genes/slice_key_reactions.txt
slices=$(awk -F'\t' '$2 == "unified" && $3 == "18" {print $9, $8}' $slice_key)
expected_slices=$(printf '%s\n' "gene_2 K01623" "gene_10 K01623" "gene_1 K01624" \
                  "gene_5 K01624,K01622" "gene_3 K01622" "gene_4 K16306")
if [ "$slices" != "$expected_slices" ]
then
    echo "ERROR: the slices of genes were not ordered by KO and then by gene ID."
    exit 1
fi
slices=$(awk -F'\t' '$2 == "SAMPLE_1" && $3 == "18" {printf "%s %.3f %.3f\n", $9, $10, $11}' \
         $slice_key)
expected_slices=$(printf '%s\n' "gene_2 1.000 0.100" "gene_10 5.000 0.200" "gene_5 4.000 0.200" \
                  "gene_3 8.000 0.500")
if [ "$slices" != "$expected_slices" ]
then
    echo "ERROR: a gene with rows for two KOs of an element did not take the means of its rows."
    exit 1
fi
# Element 64 lists K01622 but not K01624, so gene_5 there stands for its rows of K01622 alone. The
# 'unified' map takes the mean of the values of a gene in the samples where it has rows for the
# element. gene_10 has rows in SAMPLE_1 and SAMPLE_3 alone, with values of 5 and 7, so its mean is
# 6. In element 18, gene_5 has values of 4, 6, and 8 in the three samples, so its mean is 6. It has
# no row for K01622 in SAMPLE_3, so SAMPLE_3 is skipped in element 64. There its values of 6 and 9
# give a mean of 7.5.
values=$(awk -F'\t' '$2 == "unified" && ($9 == "gene_5" || $9 == "gene_10") {
                     printf "%s %s %s %.3f\n", $3, $9, $8, $10}' $slice_key)
expected_values=$(printf '%s\n' "18 gene_10 K01623 6.000" "18 gene_5 K01624,K01622 6.000" \
                  "64 gene_5 K01622 7.500")
if [ "$values" != "$expected_values" ]
then
    echo "ERROR: the values of a gene were not pooled from its own rows in each element."
    exit 1
fi
# Several genes of this file occur in two or more samples, so the warning that none do is not given.
# gene_5 has K01624 and K01622 in SAMPLE_1 and SAMPLE_2, but K01624 alone in SAMPLE_3. Its rows in
# any two samples share K01624, so its ID is not named as one that may name two genes.
for phrase in "Each slice of a gene there shows the value of a single sample" \
              "suggest that one ID names two genes"
do
    if tr '\n' ' ' < draw_txt_slices_genes.log | grep -q "$phrase"
    then
        echo "ERROR: the gene IDs of a file whose genes are shared were reported ($phrase)."
        exit 1
    fi
done

# '--reaction-accession-aggregation' reduces the rows of a gene for the KOs of an element when it
# is given. gene_1 has rows for K01623 and K01624 of element 18, with values of 2 and 5. 'sum' gives
# its slice a value of 7, where the default, 'mean', would give 3.5. gene_2 has one row, for K01622,
# and takes its value of 4.
INFO "Testing the rows of a gene reduced by their sum"
anvi-draw-kegg-pathways --reaction-txt draw_kos_genes_coverage.reaction.txt \
    --reaction-slice-by gene --reaction-accession-aggregation sum \
    --output-dir draw_txt_slices_genes_sum --pathway-numbers 00010 --no-progress
values=$(awk -F'\t' '$3 == "18" {printf "%s %s %.3f\n", $9, $8, $10}' \
         draw_txt_slices_genes_sum/slice_key_reactions.txt)
if [ "$values" != "$(printf '%s\n' "gene_1 K01623,K01624 7.000" "gene_2 K01622 4.000")" ]
then
    echo "ERROR: the rows of a gene for the KOs of an element were not reduced by their sum."
    exit 1
fi

# With '--reaction-value-period', the rows of a gene for the KOs of an element are reduced by their
# circular mean by default. gene_1 has rows for K01623 and K01624 of element 18 at 23 and 3 hours,
# an hour before and three hours after midnight. Their circular mean is 1, the hour halfway between
# them across midnight. Their arithmetic mean would be 13. The circular mean can differ from 1 by a
# tiny amount, so it is read to three decimal places. gene_2 has one row, for K01622, and takes its
# value of 12.
INFO "Testing the rows of a gene reduced by their circular mean by default"
anvi-draw-kegg-pathways --reaction-txt draw_kos_genes_clocktime.reaction.txt \
    --reaction-value-period 24 --reaction-slice-by gene \
    --output-dir draw_txt_slices_genes_period --pathway-numbers 00010 --no-progress
values=$(awk -F'\t' '$3 == "18" {printf "%s %s %.3f\n", $9, $8, $10}' \
         draw_txt_slices_genes_period/slice_key_reactions.txt)
if [ "$values" != "$(printf '%s\n' "gene_1 K01623,K01624 1.000" "gene_2 K01622 12.000")" ]
then
    echo "ERROR: the rows of a gene with a period were not reduced by their circular mean."
    exit 1
fi

# A categorical slice of a gene takes the value that the accession rule chooses among the rows of
# the gene for the element in a sample. In S1, gene_1 has rows for K01623 and K01624 of element 18,
# valued 'up' and 'down'. Each row gives one vote, and 'plurality_order', the default rule, breaks
# the tie by the order of the values. Without a colors file, the values are sorted, so 'down' comes
# first and gene_1 takes it. gene_1 has no rows in S2 and S3, so it has no slice there. gene_2 takes
# the value of its one row in S1 and S2, 'up'. Its two rows in S3 are both 'down'. gene_3 ties in S1
# as gene_1 does, but its 'down' row comes first, so the check fails if a gene takes its last row.
INFO "Testing the categorical values of slices of genes in samples and across them"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_genes_status.reaction.txt \
    --reaction-value-type categorical --reaction-sample-summary plurality_order \
    --reaction-slice-by gene --output-dir draw_txt_slices_genes_categorical \
    --pathway-numbers 00010 --draw-individual-files --no-progress
slice_key=draw_txt_slices_genes_categorical/slice_key_reactions.txt
values=$(awk -F'\t' '$2 != "unified" && $3 == "18" {print $2, $9, $10}' $slice_key)
expected_values=$(printf '%s\n' "S1 gene_1 down" "S1 gene_2 up" "S1 gene_3 down" "S2 gene_2 up" \
                  "S3 gene_2 down")
if [ "$values" != "$expected_values" ]
then
    echo "ERROR: the accession rule did not choose the value of a gene among its rows."
    exit 1
fi
# On the 'unified' map, '--reaction-sample-summary' chooses among the values of a gene in the
# samples where it has rows for the element, and each sample gives one vote. gene_2 has 'up' in S1
# and S2 and 'down' in S3, so 'plurality_order' gives it 'up', though 'down' comes first in the
# order of the values. If each row voted, 'up' and 'down' would tie at two votes, and the tie would
# go to 'down'. gene_1 and gene_3 have rows in S1 alone and keep 'down'.
values=$(awk -F'\t' '$2 == "unified" && $3 == "18" {print $9, $10}' $slice_key)
if [ "$values" != "$(printf '%s\n' "gene_1 down" "gene_2 up" "gene_3 down")" ]
then
    echo "ERROR: the sample rule did not choose the value of a gene among the samples."
    exit 1
fi

# No gene ID of this file occurs in two or more samples. On a map that pools samples, each slice of
# a gene then shows the value of a single sample, and a warning names these maps. Group G1 holds
# SAMPLE_1 alone, so its map pools nothing and is not named.
INFO "Testing the warning that no gene ID occurs in two or more samples of a map that pools them"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --groups-txt draw-sample-group-information.txt --reaction-sample-summary mean \
    --reaction-group-summary mean --reaction-slice-by gene \
    --output-dir draw_txt_slices_genes_unshared --pathway-numbers 00010 --draw-individual-files \
    --no-progress 2>&1 | tee draw_txt_slices_genes_unshared.log
for phrase in "Each slice of a gene there shows the value of a single sample" \
              "'unified' map and the map of group 'G2'"
do
    if ! tr '\n' ' ' < draw_txt_slices_genes_unshared.log | grep -q "$phrase"
    then
        echo "ERROR: maps that pool no gene across samples were misreported ($phrase)."
        exit 1
    fi
done
if tr '\n' ' ' < draw_txt_slices_genes_unshared.log | grep -q "'G1'"
then
    echo "ERROR: the map of a group of one sample was named in the warning about gene IDs."
    exit 1
fi

# Gene caller IDs start at 0 in every contigs database. Here gene 0 of one assembly has K01623 in
# S1, and gene 0 of another has K01624 in S2, so its rows in the two samples share no KO. The
# 'unified' map would pool them as one gene, so a warning names the ID. Gene 1 has the same KO in
# both samples and is not named.
INFO "Testing the warning that a gene ID has no KO in common between two samples"
anvi-draw-kegg-pathways --reaction-txt draw_kos_gene_id_collision.reaction.txt \
    --reaction-sample-summary mean --reaction-slice-by gene \
    --output-dir draw_txt_slices_genes_collision --pathway-numbers 00010 \
    --no-progress 2>&1 | tee draw_txt_slices_genes_collision.log
for phrase in "1 gene ID of the reaction file" \
              "no accession in common between its rows in two samples" \
              "'0'. A gene of one contigs database is annotated once"
do
    if ! tr '\n' ' ' < draw_txt_slices_genes_collision.log | grep -q "$phrase"
    then
        echo "ERROR: a gene ID with no KO in common between two samples was misreported ($phrase)."
        exit 1
    fi
done

# Slicing is refused where it has no values to slice or where another option rules it out (each
# command below is expected to fail).
INFO "Testing that slicing is refused where it cannot apply"
# Only a text file has values to slice, not a database run or a reaction network JSON.
if anvi-draw-kegg-pathways --external-genomes external-genomes.txt --reaction-slice-by accession \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_db_bad --overwrite-output-destinations \
    --no-progress > draw_db_slice_by.log 2>&1
then
    echo "ERROR: '--reaction-slice-by' on a database run should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_db_slice_by.log | grep -q "These options were given: '--reaction-slice-by'"
then
    echo "ERROR: '--reaction-slice-by' on a database run was refused for another reason."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-network-json pan-network.json \
    --reaction-slice-weights coverage --pathway-numbers $pathway_numbers \
    --output-dir draw_db_bad --overwrite-output-destinations \
    --no-progress > draw_json_slice_weights.log 2>&1
then
    echo "ERROR: slice weights on a reaction network JSON should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_json_slice_weights.log \
    | grep -q "These options were given: '--reaction-slice-weights'"
then
    echo "ERROR: slice weights on a reaction network JSON were refused for another reason."
    exit 1
fi

# Slicing applies to the reaction layer alone.
if anvi-draw-kegg-pathways --compound-txt draw_compounds.compound.txt \
    --reaction-slice-by accession --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_slice_no_reaction_txt.log 2>&1
then
    echo "ERROR: '--reaction-slice-by' with no reaction file should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_slice_no_reaction_txt.log \
    | grep -q "no reaction file ('--reaction-txt') was provided"
then
    echo "ERROR: '--reaction-slice-by' with no reaction file was refused for another reason."
    exit 1
fi

# The options that size and pool slices need '--reaction-slice-by'.
for slice_option in "--reaction-slice-weights coverage" "--reaction-slice-weight-aggregation mean" \
                    "--reaction-slice-share-summary pooled_sum"
do
    if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
        $slice_option --pathway-numbers $pathway_numbers \
        --output-dir draw_txt_bad --overwrite-output-destinations \
        --no-progress > draw_txt_slice_option_alone.log 2>&1
    then
        echo "ERROR: '$slice_option' without '--reaction-slice-by' should have failed but did not."
        exit 1
    fi
    if ! tr '\n' ' ' < draw_txt_slice_option_alone.log | grep -q "was not given, so no element"
    then
        echo "ERROR: '$slice_option' without '--reaction-slice-by' was refused for another reason."
        exit 1
    fi
done

# A column of weights cannot be one of the columns that identify a row. It has to exist, and it
# cannot be the only column that could hold values, since the slices then have no values to color.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-slice-by accession --reaction-slice-weights sample \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_slice_weights_key_column.log 2>&1
then
    echo "ERROR: weights in the 'sample' column should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_slice_weights_key_column.log \
    | grep -q "one of the columns that identify a row"
then
    echo "ERROR: weights in the 'sample' column were refused for another reason."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-slice-by accession --reaction-slice-weights weight \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_slice_weights_missing.log 2>&1
then
    echo "ERROR: weights in a column the file lacks should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_slice_weights_missing.log | grep -q "has no column with that header"
then
    echo "ERROR: weights in a column the file lacks were refused for another reason."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-slice-by accession --reaction-slice-weights coverage --draw-individual-files \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_slice_weights_only_column.log 2>&1
then
    echo "ERROR: weights in the only column that could hold values should have failed."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_slice_weights_only_column.log \
    | grep -q "so it would otherwise be the value column"
then
    echo "ERROR: weights in the only column that could hold values were refused for another reason."
    exit 1
fi

# Every weight is a number of 0 or more.
if anvi-draw-kegg-pathways --reaction-txt draw_bad_negative_weight.reaction.txt \
    --reaction-slice-by accession --reaction-slice-weights weight \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_slice_negative_weight.log 2>&1
then
    echo "ERROR: a negative weight should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_slice_negative_weight.log \
    | grep -q "A weight is a finite number of 0 or more"
then
    echo "ERROR: a negative weight was refused for another reason."
    exit 1
fi

# Pandas would read a second 'coverage' column as 'coverage.1'. The weights would then come from the
# first column and the values from the second.
if anvi-draw-kegg-pathways --reaction-txt draw_bad_repeated_weights.reaction.txt \
    --reaction-slice-by accession --reaction-slice-weights coverage \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_slice_repeated_weights.log 2>&1
then
    echo "ERROR: two columns named 'coverage' should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_slice_repeated_weights.log | grep -q "A name may appear only once"
then
    echo "ERROR: two columns named 'coverage' were refused for another reason."
    exit 1
fi

# Reducing the weights of the genes of a KO needs weights and genes.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-slice-by accession --reaction-slice-weight-aggregation mean \
    --draw-individual-files --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_slice_weight_aggregation_no_weights.log 2>&1
then
    echo "ERROR: a weight aggregation without weights should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_slice_weight_aggregation_no_weights.log \
    | grep -q "no column of weights was named"
then
    echo "ERROR: a weight aggregation without weights was refused for another reason."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_weights.reaction.txt \
    --reaction-slice-by accession --reaction-slice-weights weight \
    --reaction-slice-weight-aggregation mean --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_slice_weight_aggregation_no_genes.log 2>&1
then
    echo "ERROR: a weight aggregation without a 'gene_id' column should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_slice_weight_aggregation_no_genes.log \
    | grep -q "has no 'gene_id' column naming"
then
    echo "ERROR: a weight aggregation without a 'gene_id' column was refused for another reason."
    exit 1
fi

# Pooling the shares of slices needs weights, since without them the slices are equal on every map.
# It also needs a drawn map that pools samples and shows values. A file without a 'sample' column
# has no such map. Without groups, only the 'unified' map pools samples, so it has to be drawn and
# show values.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-slice-by accession --reaction-slice-share-summary pooled_sum \
    --reaction-sample-summary mean --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_slice_share_summary_no_weights.log 2>&1
then
    echo "ERROR: a share summary without weights should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_slice_share_summary_no_weights.log \
    | grep -q "no column of weights was named"
then
    echo "ERROR: a share summary without weights was refused for another reason."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_weights.reaction.txt \
    --reaction-slice-by accession --reaction-slice-weights weight \
    --reaction-slice-share-summary pooled_sum --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_slice_share_summary_no_samples.log 2>&1
then
    echo "ERROR: a share summary for a file with no 'sample' column should have failed."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_slice_share_summary_no_samples.log \
    | grep -q "has no 'sample' column, so no map pools samples"
then
    echo "ERROR: a share summary for a file with no 'sample' column was refused for another reason."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_weights.reaction.txt \
    --reaction-slice-by accession --reaction-slice-weights coverage \
    --reaction-slice-share-summary pooled_sum --draw-individual-files \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_slice_share_summary_presence.log 2>&1
then
    echo "ERROR: a share summary for a 'unified' map of presence should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_slice_share_summary_presence.log \
    | grep -q "are therefore drawn whole"
then
    echo "ERROR: a share summary for a 'unified' map of presence was refused for another reason."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_weights.reaction.txt \
    --reaction-slice-by accession --reaction-slice-weights coverage \
    --reaction-slice-share-summary pooled_sum --draw-individual-files --skip-unified-maps \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_slice_share_summary_skip_unified.log 2>&1
then
    echo "ERROR: a share summary for a skipped 'unified' map should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_slice_share_summary_skip_unified.log \
    | grep -q "draws no 'unified' map"
then
    echo "ERROR: a share summary for a skipped 'unified' map was refused for another reason."
    exit 1
fi

# Each slice takes the color of its own value, so the layer cannot take one color.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --original-color --reaction-slice-by accession --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_slice_original_color.log 2>&1
then
    echo "ERROR: '--reaction-slice-by' with '--original-color' should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_slice_original_color.log \
    | grep -q "These options were given along with '--original-color': '--reaction-slice-by'"
then
    echo "ERROR: '--reaction-slice-by' with '--original-color' was refused for another reason."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-color --reaction-slice-by accession --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_slice_reaction_color.log 2>&1
then
    echo "ERROR: '--reaction-slice-by' with '--reaction-color' should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_slice_reaction_color.log \
    | grep -q "sets a single presence color"
then
    echo "ERROR: '--reaction-slice-by' with '--reaction-color' was refused for another reason."
    exit 1
fi

# A file without a value column is colored by presence, which is not sliced yet. Rescaling each
# slice on its own is not supported yet either.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples.reaction.txt \
    --reaction-slice-by accession --draw-individual-files --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_slice_presence.log 2>&1
then
    echo "ERROR: slicing a file with no value column should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_slice_presence.log | grep -q "of presence are not supported yet"
then
    echo "ERROR: slicing a file with no value column was refused for another reason."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-slice-by accession --reaction-element-normalization relative_to_mean \
    --draw-individual-files --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_slice_normalization.log 2>&1
then
    echo "ERROR: slicing with a normalization should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_slice_normalization.log \
    | grep -q "Rescaling each slice on its own is not supported"
then
    echo "ERROR: slicing with a normalization was refused for another reason."
    exit 1
fi

# A slice of a KO shows the value of that KO. An element drawn whole takes the weighted mean of its
# slice values in place of an aggregation of numbers across its KOs.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-slice-by accession --reaction-accession-aggregation mean \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_slice_accession_aggregation.log 2>&1
then
    echo "ERROR: slicing with an accession aggregation of numbers should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_slice_accession_aggregation.log \
    | grep -q "in place of this reduction"
then
    echo "ERROR: slicing with an accession aggregation of numbers was refused for another reason."
    exit 1
fi

# Slicing needs a drawn map that shows values. The refusal names an option that would draw one, in
# words that fit the run: without groups, with group maps of values, and with group maps of
# presence, drawn or not.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --reaction-slice-by accession --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_slice_no_value_maps.log 2>&1
then
    echo "ERROR: slicing with no drawn map of values should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_slice_no_value_maps.log \
    | grep -q "The map of each sample shows values, but none were"
then
    echo "ERROR: slicing with no drawn map of values was refused for another reason."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --groups-txt draw-sample-group-information.txt --group-threshold 0.5 \
    --reaction-sample-summary mean --reaction-slice-by accession \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_slice_no_value_maps_groups.log 2>&1
then
    echo "ERROR: slicing with no drawn group map of values should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_slice_no_value_maps_groups.log \
    | grep -q "The map of each group shows values"
then
    echo "ERROR: slicing with no drawn group map of values was refused for another reason."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --groups-txt draw-sample-group-information.txt --group-threshold 0.5 \
    --draw-individual-files --reaction-slice-by accession \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_slice_group_presence.log 2>&1
then
    echo "ERROR: slicing group maps of presence should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_slice_group_presence.log \
    | grep -q "The map of each group counts its samples"
then
    echo "ERROR: slicing group maps of presence was refused for another reason."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage.reaction.txt \
    --groups-txt draw-sample-group-information.txt --group-threshold 0.5 \
    --reaction-slice-by accession --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_slice_group_presence_no_maps.log 2>&1
then
    echo "ERROR: slicing undrawn group maps of presence should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_slice_group_presence_no_maps.log \
    | grep -q "Then either draw those maps"
then
    echo "ERROR: slicing undrawn group maps of presence was refused for another reason."
    exit 1
fi

# The refusal names an option that fits the values. Categorical values need a rule, and values with
# a period need a circular mean.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_status.reaction.txt \
    --reaction-value-type categorical --reaction-slice-by accession \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_slice_no_value_maps_categorical.log 2>&1
then
    echo "ERROR: slicing categorical values with no drawn map of values should have failed."
    exit 1
fi
for phrase in "no drawn map of the reaction layer" "a rule such as 'plurality_order'"
do
    if ! tr '\n' ' ' < draw_txt_slice_no_value_maps_categorical.log | grep -q "$phrase"
    then
        echo "ERROR: slicing categorical values with no map of values was misreported ($phrase)."
        exit 1
    fi
done

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_clocktime.reaction.txt \
    --reaction-value-period 24 --reaction-slice-by accession \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_slice_no_value_maps_period.log 2>&1
then
    echo "ERROR: slicing values with a period and no drawn map of values should have failed."
    exit 1
fi
for phrase in "no drawn map of the reaction layer" "an aggregation such as 'circular_mean'"
do
    if ! tr '\n' ' ' < draw_txt_slice_no_value_maps_period.log | grep -q "$phrase"
    then
        echo "ERROR: slicing values with a period and no map of values was misreported ($phrase)."
        exit 1
    fi
done

# Slices of genes need a 'gene_id' column, and a gene ID in every row.
if anvi-draw-kegg-pathways --reaction-txt draw_reactions.reaction.txt \
    --reaction-slice-by gene --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_slice_genes_no_gene_id.log 2>&1
then
    echo "ERROR: slicing by gene without a 'gene_id' column should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_slice_genes_no_gene_id.log \
    | grep -q "draws each gene of a reaction element as a slice"
then
    echo "ERROR: slicing by gene without a 'gene_id' column was refused for another reason."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_bad_blank_gene_id.reaction.txt \
    --reaction-slice-by gene --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_slice_genes_blank_gene_id.log 2>&1
then
    echo "ERROR: slicing by gene with a row without a gene ID should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_slice_genes_blank_gene_id.log | grep -q "in its 'gene_id' column"
then
    echo "ERROR: slicing by gene with a row without a gene ID was refused for another reason."
    exit 1
fi

# Each gene is a slice of its own, so the values and the weights of the genes of a KO are never
# reduced to those of the KO.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-slice-by gene --reaction-gene-aggregation mean \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_slice_genes_gene_aggregation.log 2>&1
then
    echo "ERROR: slicing by gene with a gene aggregation should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_slice_genes_gene_aggregation.log \
    | grep -q "so no values of genes are reduced to an accession"
then
    echo "ERROR: slicing by gene with a gene aggregation was refused for another reason."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_weights.reaction.txt \
    --reaction-slice-by gene --reaction-slice-weights coverage \
    --reaction-slice-weight-aggregation mean --draw-individual-files \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_slice_genes_weight_aggregation.log 2>&1
then
    echo "ERROR: slicing by gene with a weight aggregation should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_slice_genes_weight_aggregation.log \
    | grep -q "so no weights are reduced to an accession"
then
    echo "ERROR: slicing by gene with a weight aggregation was refused for another reason."
    exit 1
fi

# A gene with one row for an element shows the value of that row. 'std' does not return a single
# value unchanged, so it cannot reduce the rows of a gene.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-slice-by gene --reaction-accession-aggregation std \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_slice_genes_std.log 2>&1
then
    echo "ERROR: slicing by gene with 'std' over a gene's rows should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_slice_genes_std.log \
    | grep -q "return a single value unchanged, so genes with several rows"
then
    echo "ERROR: slicing by gene with 'std' over a gene's rows was refused for another reason."
    exit 1
fi

# A border column colors a ring inside the outline of each box on standard maps, and of each
# compound circle on standard and overview maps. A row with an empty border cell has no border
# value, and the other rows still color their rings.
INFO "Testing a border column with an empty cell beside a value column"
anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage_status.reaction.txt \
    --reaction-border-column status --reaction-border-value-type categorical \
    --output-dir draw_txt_border_empty_cell --pathway-numbers $pathway_numbers \
    --no-progress 2>&1 | tee draw_txt_border_empty_cell.log
if ! tr '\n' ' ' < draw_txt_border_empty_cell.log | grep -q "Reaction rows without a border value"
then
    echo "ERROR: the rows without a border value were not counted."
    exit 1
fi
if [ ! -s draw_txt_border_empty_cell/colorbar_reactions_border.pdf ]
then
    echo "ERROR: a categorical border should have written a legend of its own."
    exit 1
fi

# A stray quote in a categorical border value is part of the value. The file is read as text, even
# though its fill values are numbers.
INFO "Testing a stray quote in a categorical border value of a file of numbers"
anvi-draw-kegg-pathways --reaction-txt draw_kos_border_stray_quote.reaction.txt \
    --reaction-border-column status --reaction-border-value-type categorical \
    --output-dir draw_txt_border_stray_quote --pathway-numbers $pathway_numbers \
    --no-progress

# K00844 has no fill value, since the standard deviation of one gene is undefined. Its border value
# still draws a ring, so map 00010 is drawn though nothing on it has a fill color.
INFO "Testing an element with a border and no fill color"
anvi-draw-kegg-pathways --reaction-txt draw_kos_border_no_fill.reaction.txt \
    --reaction-gene-aggregation std --reaction-border-column level \
    --output-dir draw_txt_border_no_fill --pathway-numbers 00010 \
    --no-progress
if [ ! -s draw_txt_border_no_fill/unified/ko00010.pdf ]
then
    echo "ERROR: a map with a ring and no fill color should have been drawn."
    exit 1
fi

# Rings are drawn on circles of standard and overview maps, not on the global map. The run counts
# the maps that show no compound border.
INFO "Testing a compound border on maps of every kind"
anvi-draw-kegg-pathways --compound-txt draw_compounds_samples.compound.txt \
    --compound-border-column concentration --compound-border-sample-summary mean \
    --output-dir draw_txt_compound_border_global --pathway-numbers $global_pathway_numbers \
    --no-progress 2>&1 | tee draw_txt_compound_border_global.log
if ! tr '\n' ' ' < draw_txt_compound_border_global.log \
    | grep -q "Maps that show no compound borders"
then
    echo "ERROR: the maps without compound borders were not counted."
    exit 1
fi

# A player shows the key of a border only beside a map that can show its rings. The player of the
# standard map shows the keys of the coverage, the reaction border, and the compound border, in that
# order. The overview map shows compound rings but no reaction rings, and the global map shows no
# rings, so their players leave out the keys of the borders they lack. A key is an image that holds
# its data, so the same key is the same element in every player.
INFO "Testing the keys of borders in the players of maps of every kind"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage_status.reaction.txt \
    --reaction-border-column status --reaction-border-value-type categorical \
    --compound-txt draw_compounds_samples.compound.txt --compound-border-column concentration \
    --output-dir draw_txt_border_players --pathway-numbers $global_pathway_numbers \
    --draw-individual-files --collate-files-by-map --animate-files-by-map \
    --no-progress
players=draw_txt_border_players/by_map
standard_keys=$(grep -o '<img class="key"[^>]*' $players/ko00010.html)
overview_keys=$(grep -o '<img class="key"[^>]*' $players/ko01200.html)
global_keys=$(grep -o '<img class="key"[^>]*' $players/ko01100.html)
if [ "$(echo "$standard_keys" | wc -l)" -ne 3 ] \
    || [ "$overview_keys" != "$(echo "$standard_keys" | sed -n '1p;3p')" ] \
    || [ "$global_keys" != "$(echo "$standard_keys" | sed -n 1p)" ]
then
    echo "ERROR: a player should show the key of a border only beside maps that show its rings."
    exit 1
fi

# Lines show no ring. A reaction element that a standard map draws as a line is reported.
INFO "Testing the warning that reaction lines on a standard map show no border"
anvi-draw-kegg-pathways --reaction-txt draw_kos_border_lines.reaction.txt \
    --reaction-border-column status --reaction-border-value-type categorical \
    --output-dir draw_txt_border_lines --pathway-numbers 00010 00510 \
    --no-progress 2>&1 | tee draw_txt_border_lines.log
if ! tr '\n' ' ' < draw_txt_border_lines.log | grep -q "Lines show no border"
then
    echo "ERROR: reaction lines with border values on a standard map were not reported."
    exit 1
fi

# A border summary that is undefined for an element leaves it without a ring.
INFO "Testing the warning that a border summary was undefined"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage_status.reaction.txt \
    --reaction-border-column coverage --reaction-value-type categorical \
    --reaction-border-sample-summary std \
    --output-dir draw_txt_border_std --pathway-numbers $pathway_numbers \
    --no-progress 2>&1 | tee draw_txt_border_std.log
if ! tr '\n' ' ' < draw_txt_border_std.log | grep -q "show no border on the 'unified' map"
then
    echo "ERROR: elements with an undefined border summary were not reported."
    exit 1
fi

# Choices among categorical border values are reported apart from those of the fill.
INFO "Testing that choices of a categorical border value are reported"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_status.reaction.txt \
    --reaction-border-column status --reaction-border-value-type categorical \
    --reaction-border-sample-summary plurality_random --seed 2 \
    --output-dir draw_txt_border_choices --pathway-numbers $pathway_numbers \
    --no-progress 2>&1 | tee draw_txt_border_choices.log
if ! tr '\n' ' ' < draw_txt_border_choices.log | grep -q "Choices of a reaction border value among"
then
    echo "ERROR: the choices of categorical border values were not reported."
    exit 1
fi

# Skipping the 'unified' map drops the key of its border and keeps the key of the sample maps.
INFO "Testing that skipping the 'unified' map drops its border key"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage_status.reaction.txt \
    --reaction-border-column status --reaction-border-value-type categorical \
    --output-dir draw_txt_border_skip_unified --pathway-numbers $pathway_numbers \
    --skip-unified-maps --draw-grid \
    --no-progress
if [ -e draw_txt_border_skip_unified/colorbar_reactions_border.pdf ]
then
    echo "ERROR: the border key of the skipped 'unified' map was drawn."
    exit 1
fi
if [ ! -s draw_txt_border_skip_unified/colorbar_reactions_border_samples.pdf ]
then
    echo "ERROR: skipping the 'unified' map should keep the border key of the sample maps."
    exit 1
fi

# Border keys stay at the top of the output directory, so gathering maps by map leaves them out.
INFO "Testing that gathering maps by map gathers no border key"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage_status.reaction.txt \
    --reaction-border-column status --reaction-border-value-type categorical \
    --output-dir draw_txt_border_collate --pathway-numbers $pathway_numbers \
    --draw-individual-files --collate-files-by-map \
    --no-progress
if ls -R draw_txt_border_collate/by_map | grep -q "colorbar"
then
    echo "ERROR: a border key was gathered by map as if it were a map."
    exit 1
fi

INFO "Testing that border columns are refused where they cannot apply"
# A border column is read like the column of weights. It must be a column of the file, and it
# cannot be a key column, a name given twice, or a name that starts with '__'.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage_status.reaction.txt \
    --reaction-border-column nothing --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_border_missing.log 2>&1
then
    echo "ERROR: a border column that the file lacks should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_border_missing.log | grep -q "as the border column, but"
then
    echo "ERROR: a border column that the file lacks was refused for another reason."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage_status.reaction.txt \
    --reaction-border-column sample --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_border_key_column.log 2>&1
then
    echo "ERROR: a key column as the border column should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_border_key_column.log | grep -q "key column, which cannot hold border values"
then
    echo "ERROR: a key column as the border column was refused for another reason."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_bad_repeated_border.reaction.txt \
    --reaction-border-column status --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_border_repeated.log 2>&1
then
    echo "ERROR: a border column named twice should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_border_repeated.log | grep -q "more than once"
then
    echo "ERROR: a border column named twice was refused for another reason."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_bad_reserved_border.reaction.txt \
    --reaction-border-column __status --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_border_reserved.log 2>&1
then
    echo "ERROR: a border column named with '__' should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_border_reserved.log | grep -q "that start with '__'"
then
    echo "ERROR: a border column named with '__' was refused for another reason."
    exit 1
fi

# One column cannot hold both the weights of slices and the border values.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_weights.reaction.txt \
    --reaction-slice-by accession --reaction-slice-weights coverage \
    --reaction-border-column coverage --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_border_weights.log 2>&1
then
    echo "ERROR: one column as both weights and border should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_border_weights.log | grep -q "can either color borders or hold"
then
    echo "ERROR: one column as both weights and border was refused for another reason."
    exit 1
fi

# A border column leaves two other candidate value columns in this file.
if anvi-draw-kegg-pathways --reaction-txt draw_bad_three_columns.reaction.txt \
    --reaction-border-column status --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_border_three_columns.log 2>&1
then
    echo "ERROR: a file with two value columns besides its border column should have failed."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_border_three_columns.log | grep -q "is ambiguous which to color elements by"
then
    echo "ERROR: a file with two value columns besides its border column failed for another reason."
    exit 1
fi

# The values of a border column: none at all, 'NA' among numbers, more categorical values than
# default colors, and repeated rows in a file whose only values are border values.
if anvi-draw-kegg-pathways --reaction-txt draw_bad_empty_border.reaction.txt \
    --reaction-border-column status --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_border_empty.log 2>&1
then
    echo "ERROR: a border column without values should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_border_empty.log | grep -q "has no values"
then
    echo "ERROR: a border column without values was refused for another reason."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_bad_na_border.reaction.txt \
    --reaction-border-column level --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_border_na.log 2>&1
then
    echo "ERROR: 'NA' in a column of numeric border values should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_border_na.log | grep -q "such as text or 'NA'"
then
    echo "ERROR: 'NA' in a column of numeric border values was refused for another reason."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_border_class8.reaction.txt \
    --reaction-border-column status --reaction-border-value-type categorical \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_border_class8.log 2>&1
then
    echo "ERROR: more categorical border values than default colors should have failed."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_border_class8.log | grep -q "default colors for"
then
    echo "ERROR: more categorical border values than default colors failed for another reason."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_bad_repeated_border_rows.reaction.txt \
    --reaction-border-column status --reaction-border-value-type categorical \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_border_repeated_rows.log 2>&1
then
    echo "ERROR: repeated rows of a file with only border values should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_border_repeated_rows.log | grep -q "only one border value"
then
    echo "ERROR: repeated rows of a file with only border values failed for another reason."
    exit 1
fi

# Every other border option needs the border column of its layer.
for border_option in "--reaction-border-value-type categorical" \
                     "--reaction-border-value-colors draw-status-colors.txt" \
                     "--reaction-border-value-period 24" \
                     "--reaction-border-gene-aggregation mean" \
                     "--reaction-border-accession-aggregation mean" \
                     "--reaction-border-sample-summary mean" \
                     "--reaction-border-group-summary mean" \
                     "--reaction-border-colormap viridis" \
                     "--reaction-border-individual-colormap viridis" \
                     "--reaction-border-value-limits 0 1" \
                     "--reaction-border-individual-value-limits 0 1" \
                     "--reaction-border-value-center" \
                     "--reaction-border-individual-value-center" \
                     "--reaction-border-element-normalization z_score" \
                     "--reaction-border-slice-by accession" \
                     "--compound-border-sample-summary mean"
do
    if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage_status.reaction.txt \
        $border_option --pathway-numbers $pathway_numbers \
        --output-dir draw_txt_bad --overwrite-output-destinations \
        --no-progress > draw_txt_border_option_alone.log 2>&1
    then
        echo "ERROR: '$border_option' without a border column should have failed but did not."
        exit 1
    fi
    if ! tr '\n' ' ' < draw_txt_border_option_alone.log | grep -q "no border column was named"
    then
        echo "ERROR: '$border_option' without a border column was refused for another reason."
        exit 1
    fi
done

# A border shows values only, so its summaries take no presence names.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage_status.reaction.txt \
    --reaction-border-column status --reaction-border-value-type categorical \
    --reaction-border-sample-summary count --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_border_presence_summary.log 2>&1
then
    echo "ERROR: a presence summary of a border should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_border_presence_summary.log | grep -q "A border shows values only"
then
    echo "ERROR: a presence summary of a border was refused for another reason."
    exit 1
fi

# The colormap of a border colors numbers only.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage_status.reaction.txt \
    --reaction-border-column status --reaction-border-value-type categorical \
    --reaction-border-colormap viridis --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_border_categorical_colormap.log 2>&1
then
    echo "ERROR: a colormap of categorical border values should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_border_categorical_colormap.log \
    | grep -q "colors numeric border values"
then
    echo "ERROR: a colormap of categorical border values was refused for another reason."
    exit 1
fi

# A colors file needs categorical border values.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage_status.reaction.txt \
    --reaction-value-type categorical --reaction-border-column coverage \
    --reaction-border-value-colors draw-status-colors.txt --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_border_numeric_colors.log 2>&1
then
    echo "ERROR: colors for numeric border values should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_border_numeric_colors.log | grep -q "Give it to read the border column as categorical"
then
    echo "ERROR: colors for numeric border values were refused for another reason."
    exit 1
fi

# The border value of a group comes from the border sample summary, so a group summary needs one.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage_status.reaction.txt \
    --reaction-value-type categorical --reaction-border-column coverage \
    --groups-txt draw-sample-group-information.txt --group-threshold 0.5 \
    --reaction-border-group-summary mean --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_border_group_summary.log 2>&1
then
    echo "ERROR: a border group summary without a border sample summary should have failed."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_border_group_summary.log | grep -q "which was not given"
then
    echo "ERROR: a border group summary without a border sample summary failed for another reason."
    exit 1
fi

# The reference map's colors leave no colors for values, and they are drawn without rings.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt --original-color \
    --reaction-border-column coverage --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_border_original.log 2>&1
then
    echo "ERROR: '--original-color' with a border column should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_border_original.log | grep -q "along with '--original-color'"
then
    echo "ERROR: '--original-color' with a border column was refused for another reason."
    exit 1
fi

# A border that no drawn map can show. Reaction borders need a standard map, and compound borders a
# standard or an overview map.
for border_map in 01200 01100
do
    if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage_status.reaction.txt \
        --reaction-border-column status --reaction-border-value-type categorical \
        --reaction-border-sample-summary plurality_order --pathway-numbers $border_map \
        --output-dir draw_txt_bad --overwrite-output-destinations \
        --no-progress > draw_txt_border_map_kind.log 2>&1
    then
        echo "ERROR: a reaction border drawn only on map $border_map should have failed."
        exit 1
    fi
    if ! tr '\n' ' ' < draw_txt_border_map_kind.log | grep -q "This run draws no such map"
    then
        echo "ERROR: a reaction border drawn only on map $border_map failed for another reason."
        exit 1
    fi
done

if anvi-draw-kegg-pathways --compound-txt draw_compounds_samples.compound.txt \
    --compound-border-column concentration --compound-border-sample-summary mean \
    --pathway-numbers 01100 \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_compound_border_map_kind.log 2>&1
then
    echo "ERROR: a compound border drawn only on map 01100 should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_compound_border_map_kind.log | grep -q "This run draws no such map"
then
    echo "ERROR: a compound border drawn only on map 01100 was refused for another reason."
    exit 1
fi

# The options alone can leave the border off every map. Here the 'unified' map shows no border
# without a border summary, and no map of a sample is drawn.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage_status.reaction.txt \
    --reaction-border-column status --reaction-border-value-type categorical \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_border_no_map.log 2>&1
then
    echo "ERROR: a border that no map shows should have failed but did not."
    exit 1
fi

# A database run has no text file and no border column.
if anvi-draw-kegg-pathways --contigs-dbs E_faecalis_6240.db --reaction-border-column status \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_db_bad --overwrite-output-destinations \
    --no-progress > draw_db_border.log 2>&1
then
    echo "ERROR: a border column on a database run should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_db_border.log \
    | grep -q "These options were given: '--reaction-border-column'"
then
    echo "ERROR: a border column on a database run was refused for another reason."
    exit 1
fi

# Messages from four places name the options of the border, not those of the fill: the builder of
# a border model, the check of the maps that show a border, a guard of the command line looped over
# the borders, and the check of a normalization. ConfigError breaks its lines after hyphens, so the
# spaces and line breaks are deleted before the flags are matched.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage_status.reaction.txt \
    --reaction-value-type categorical --reaction-border-column coverage \
    --reaction-border-value-limits 0 1 --draw-individual-files \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_border_message_builder.log 2>&1
then
    echo "ERROR: limits on a 'unified' map that shows no border should have failed but did not."
    exit 1
fi
if ! tr -d ' \n' < draw_txt_border_message_builder.log \
    | grep -q "since'--reaction-border-sample-summary'wasnotgiven"
then
    echo "ERROR: limits on a 'unified' map without a border named the wrong option."
    exit 1
fi
if ! tr -d ' \n' < draw_txt_border_no_map.log \
    | grep -q "'--reaction-border-sample-summary'wasnotgiven,sothe'unified'mapshowsnoborder"
then
    echo "ERROR: a border that no map shows named the wrong option."
    exit 1
fi
if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage_status.reaction.txt \
    --reaction-border-column status --reaction-border-value-type categorical \
    --reaction-border-individual-value-limits 0 1 --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_border_message_guard.log 2>&1
then
    echo "ERROR: limits of categorical border values should have failed but did not."
    exit 1
fi
if ! tr -d ' \n' < draw_txt_border_message_guard.log \
    | grep -q "('--reaction-border-value-typecategorical')"
then
    echo "ERROR: limits of categorical border values named the wrong option."
    exit 1
fi
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage_status.reaction.txt \
    --reaction-value-type categorical --reaction-border-column coverage \
    --reaction-border-element-normalization mean --draw-individual-files \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_border_message_normalization.log 2>&1
then
    echo "ERROR: a border normalization that reduces values should have failed but did not."
    exit 1
fi
if ! tr -d ' \n' < draw_txt_border_message_normalization.log \
    | grep -q "belongsto'--reaction-border-sample-summary'instead"
then
    echo "ERROR: a border normalization that reduces values named the wrong option."
    exit 1
fi

# '--reaction-border-slice-by gene' divides each ring among the genes of its box. With
# '--reaction-slice-by gene', the rings are aligned with the fill, so each slice of a ring is a
# column of the fill. Its position, accession, gene, and share are those of the column.
INFO "Testing rings sliced by gene and aligned with fills sliced by gene"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_weights_status.reaction.txt \
    --reaction-sample-summary mean --reaction-slice-by gene --reaction-slice-weights coverage \
    --reaction-border-column status --reaction-border-value-type categorical \
    --reaction-border-sample-summary plurality_order --reaction-border-slice-by gene \
    --output-dir draw_txt_border_slices_genes_aligned --pathway-numbers 00010 \
    --draw-individual-files \
    --no-progress
if [ "$(awk -F'\t' 'FNR > 1 && $9 != ""' \
    draw_txt_border_slices_genes_aligned/slice_key_reactions_border.txt | wc -l)" -eq 0 ]
then
    echo "ERROR: rings sliced by gene should have listed the genes of their slices."
    exit 1
fi
if awk -F'\t' 'NR == FNR {fill[$1 FS $2 FS $3 FS $7 FS $8 FS $9 FS $11]; next}
    FNR > 1 && !(($1 FS $2 FS $3 FS $7 FS $8 FS $9 FS $11) in fill) {bad++} END {exit !bad}' \
    draw_txt_border_slices_genes_aligned/slice_key_reactions.txt \
    draw_txt_border_slices_genes_aligned/slice_key_reactions_border.txt
then
    echo "ERROR: an aligned ring has a slice that no column of the fill has."
    exit 1
fi

# A file without a value column colors its fills by presence. Its rings can still be sliced by
# gene, each on its own.
INFO "Testing rings sliced by gene in a file without a value column"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_genes_border.reaction.txt \
    --reaction-border-column status --reaction-border-value-type categorical \
    --reaction-border-slice-by gene \
    --output-dir draw_txt_border_slices_genes_no_fill --pathway-numbers 00010 \
    --draw-individual-files \
    --no-progress 2>&1 | tee draw_txt_border_slices_genes_no_fill.log
if ! tr '\n' ' ' < draw_txt_border_slices_genes_no_fill.log \
    | grep -q "by gene, independent of the fill"
then
    echo "ERROR: rings sliced by gene without a value column were not sliced on their own."
    exit 1
fi
if [ "$(awk -F'\t' 'FNR > 1 && $9 != ""' \
    draw_txt_border_slices_genes_no_fill/slice_key_reactions_border.txt | wc -l)" -eq 0 ]
then
    echo "ERROR: rings sliced by gene without a value column should have listed their genes."
    exit 1
fi

# The weights size the slices of rings that are sliced on their own, so a file needs no value
# column for them. 'pooled_sum' adds up the weights of the rows with a border value. On the
# 'unified' map, K00174 weighs 96 and K00175 weighs 19 in box 154, so K00174 takes 96 / 115 of the
# ring.
INFO "Testing weights that size the slices of rings in a file without a value column"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_weights_border.reaction.txt \
    --reaction-slice-weights coverage --reaction-slice-share-summary pooled_sum \
    --reaction-border-column status --reaction-border-value-type categorical \
    --reaction-border-sample-summary plurality_order --reaction-border-slice-by accession \
    --output-dir draw_txt_border_slices_weights --pathway-numbers 00010 \
    --no-progress
if [ "$(awk -F'\t' '$2 == "unified" && $3 == "154" && $8 == "K00174" {print $11}' \
    draw_txt_border_slices_weights/slice_key_reactions_border.txt)" != "0.8347826086956521" ]
then
    echo "ERROR: the weights of the rows with border values did not size the slice of a ring."
    exit 1
fi

# The weights of the genes of an accession are reduced where the rings are sliced by accession,
# even where the fill is sliced by gene.
INFO "Testing a weight aggregation for rings sliced by accession around fills sliced by gene"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_weights_status.reaction.txt \
    --reaction-sample-summary mean --reaction-slice-by gene --reaction-slice-weights coverage \
    --reaction-slice-weight-aggregation mean --reaction-border-column status \
    --reaction-border-value-type categorical --reaction-border-sample-summary plurality_order \
    --reaction-border-slice-by accession \
    --output-dir draw_txt_border_slices_weight_aggregation --pathway-numbers 00010 \
    --draw-individual-files \
    --no-progress
if [ ! -s draw_txt_border_slices_weight_aggregation/slice_key_reactions_border.txt ]
then
    echo "ERROR: rings sliced by accession around fills sliced by gene should have been drawn."
    exit 1
fi

# The fill of a box inside an aligned ring holds up to 10 slices, so a fill of 11 genes is drawn
# whole, and so is its ring. A ring sliced on its own holds up to 11 slices.
INFO "Testing the caps of sliced rings"
anvi-draw-kegg-pathways --reaction-txt draw_kos_genes11_status.reaction.txt \
    --reaction-slice-by gene --reaction-border-column status \
    --reaction-border-value-type categorical --reaction-border-slice-by gene \
    --output-dir draw_txt_border_slices_cap_aligned --pathway-numbers 00010 \
    --no-progress
if [ "$(awk -F'\t' '$13 == "over the cap"' \
    draw_txt_border_slices_cap_aligned/slice_key_reactions_border.txt | wc -l)" -ne 22 ]
then
    echo "ERROR: aligned rings around fills over the cap were not drawn whole."
    exit 1
fi
anvi-draw-kegg-pathways --reaction-txt draw_kos_genes11_status.reaction.txt \
    --reaction-border-column status --reaction-border-value-type categorical \
    --reaction-border-slice-by gene \
    --output-dir draw_txt_border_slices_cap_fits --pathway-numbers 00010 \
    --no-progress
if [ "$(awk -F'\t' '$7 == "11"' \
    draw_txt_border_slices_cap_fits/slice_key_reactions_border.txt | wc -l)" -ne 2 ]
then
    echo "ERROR: rings of 11 slices sliced on their own should have been drawn in slices."
    exit 1
fi
anvi-draw-kegg-pathways --reaction-txt draw_kos_genes12_status.reaction.txt \
    --reaction-border-column status --reaction-border-value-type categorical \
    --reaction-border-slice-by gene \
    --output-dir draw_txt_border_slices_cap_over --pathway-numbers 00010 \
    --no-progress 2>&1 | tee draw_txt_border_slices_cap_over.log
if ! tr '\n' ' ' < draw_txt_border_slices_cap_over.log | grep -q "have had more than 11 slices"
then
    echo "ERROR: rings of 12 slices were not reported as drawn whole."
    exit 1
fi

# The 'unified' map shows presence in its fills here, so its fills are not sliced. Its aligned
# rings are drawn whole, while the rings of the maps of the samples are sliced with the fills.
INFO "Testing whole rings around fills that are not sliced"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_weights_status.reaction.txt \
    --reaction-slice-by accession --reaction-slice-weights coverage \
    --reaction-border-column status --reaction-border-value-type categorical \
    --reaction-border-sample-summary plurality_order --reaction-border-slice-by accession \
    --output-dir draw_txt_border_slices_not_sliced --pathway-numbers 00010 \
    --draw-individual-files \
    --no-progress
if [ "$(awk -F'\t' '$2 == "unified" && $13 != "fill not sliced"' \
    draw_txt_border_slices_not_sliced/slice_key_reactions_border.txt | wc -l)" -ne 0 ]
then
    echo "ERROR: the rings of a 'unified' map whose fills are not sliced were not drawn whole."
    exit 1
fi
if [ "$(awk -F'\t' '$2 == "SAMPLE_1" && $7 == "2"' \
    draw_txt_border_slices_not_sliced/slice_key_reactions_border.txt | wc -l)" -eq 0 ]
then
    echo "ERROR: the rings of the map of a sample should have been sliced with the fills."
    exit 1
fi

# Weights of 0 leave out slices. An aligned ring cannot show a border value whose accession has
# no slice of the fill, and a ring all of whose slices weigh 0 is not drawn.
INFO "Testing border values that no column shows, and rings left out by weights of 0"
anvi-draw-kegg-pathways --reaction-txt draw_kos_weights_status.reaction.txt \
    --reaction-slice-by accession --reaction-slice-weights weight \
    --reaction-border-column status --reaction-border-value-type categorical \
    --reaction-border-slice-by accession \
    --output-dir draw_txt_border_slices_zero_aligned --pathway-numbers 00010 \
    --no-progress 2>&1 | tee draw_txt_border_slices_zero_aligned.log
for phrase in "so no band shows their border values" "every slice of the fill had a weight"
do
    if ! tr '\n' ' ' < draw_txt_border_slices_zero_aligned.log | grep -q "$phrase"
    then
        echo "ERROR: the warning that contains '$phrase' was not given."
        exit 1
    fi
done
anvi-draw-kegg-pathways --reaction-txt draw_kos_weights_status.reaction.txt \
    --reaction-slice-weights weight --reaction-border-column status \
    --reaction-border-value-type categorical --reaction-border-slice-by accession \
    --output-dir draw_txt_border_slices_zero_independent --pathway-numbers 00010 \
    --no-progress 2>&1 | tee draw_txt_border_slices_zero_independent.log
if ! tr '\n' ' ' < draw_txt_border_slices_zero_independent.log \
    | grep -q "every accession with a border value"
then
    echo "ERROR: a ring all of whose slices weigh 0 was not reported."
    exit 1
fi

INFO "Testing that sliced rings are refused where they cannot apply"
# Rings sliced by gene need a 'gene_id' column.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage_status.reaction.txt \
    --reaction-border-column status --reaction-border-value-type categorical \
    --reaction-border-slice-by gene --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_border_slices_no_genes.log 2>&1
then
    echo "ERROR: rings sliced by gene without a 'gene_id' column should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_border_slices_no_genes.log \
    | grep -q "has no 'gene_id' column naming genes"
then
    echo "ERROR: rings sliced by gene without a 'gene_id' column were refused for another reason."
    exit 1
fi

# Each gene is a slice of its own, so no border values of genes are reduced to an accession.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_genes_border.reaction.txt \
    --reaction-border-column status --reaction-border-value-type categorical \
    --reaction-border-gene-aggregation plurality_order --reaction-border-slice-by gene \
    --draw-individual-files --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_border_slices_gene_aggregation.log 2>&1
then
    echo "ERROR: a border gene aggregation for rings sliced by gene should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_border_slices_gene_aggregation.log \
    | grep -q "so no border values of genes are reduced"
then
    echo "ERROR: a border gene aggregation for rings sliced by gene was refused for another reason."
    exit 1
fi

# Each slice of a ring sliced by accession shows the border value of one accession.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage_status.reaction.txt \
    --reaction-value-type categorical --reaction-border-column coverage \
    --reaction-border-sample-summary mean --reaction-border-accession-aggregation max \
    --reaction-border-slice-by accession --draw-individual-files \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_border_slices_accession_aggregation.log 2>&1
then
    echo "ERROR: a border accession aggregation for rings sliced by accession did not fail."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_border_slices_accession_aggregation.log \
    | grep -q "each slice of a ring shows the border value of one accession"
then
    echo "ERROR: a border accession aggregation for rings sliced by accession failed otherwise."
    exit 1
fi

# The slices of rings cannot be normalized yet.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage_status.reaction.txt \
    --reaction-value-type categorical --reaction-border-column coverage \
    --reaction-border-element-normalization z_score --reaction-border-slice-by accession \
    --draw-individual-files --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_border_slices_normalization.log 2>&1
then
    echo "ERROR: a border normalization with sliced rings should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_border_slices_normalization.log \
    | grep -q "Rescaling each slice of a ring on its own"
then
    echo "ERROR: a border normalization with sliced rings was refused for another reason."
    exit 1
fi

# A gene with one row shows the border value of that row, so an aggregation of the rows of a gene
# must return a single value unchanged. The mapper refuses 'std'.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage_status.reaction.txt \
    --reaction-value-type categorical --reaction-border-column coverage \
    --reaction-border-sample-summary mean --reaction-border-accession-aggregation std \
    --reaction-border-slice-by gene --draw-individual-files \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_border_slices_genes_std.log 2>&1
then
    echo "ERROR: rings sliced by gene with 'std' over a gene's rows should have failed but did not."
    exit 1
fi

# A reaction border slices the rings of boxes, which only standard maps draw.
for pathway_number in 01200 01100
do
    if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_genes_border.reaction.txt \
        --reaction-border-column status --reaction-border-value-type categorical \
        --reaction-border-slice-by accession --draw-individual-files \
        --pathway-numbers $pathway_number \
        --output-dir draw_txt_bad --overwrite-output-destinations \
        --no-progress > draw_txt_border_slices_map_kind.log 2>&1
    then
        echo "ERROR: rings sliced on map $pathway_number alone should have failed but did not."
        exit 1
    fi
    if ! tr '\n' ' ' < draw_txt_border_slices_map_kind.log | grep -q "This run draws no such map"
    then
        echo "ERROR: rings sliced on map $pathway_number alone were refused for another reason."
        exit 1
    fi
done

# The weights size the slices of elements or of rings, so they need one kind of slicing.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_weights_status.reaction.txt \
    --reaction-slice-weights coverage --reaction-border-column status \
    --reaction-border-value-type categorical --draw-individual-files \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_border_slices_weights_neither.log 2>&1
then
    echo "ERROR: weights with neither kind of slicing should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_border_slices_weights_neither.log \
    | grep -q "so no element or ring is sliced"
then
    echo "ERROR: weights with neither kind of slicing were refused for another reason."
    exit 1
fi

# The weights of the genes of an accession are reduced only where something is sliced by accession.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_weights_status.reaction.txt \
    --reaction-slice-by gene --reaction-slice-weights coverage \
    --reaction-slice-weight-aggregation mean --reaction-border-column status \
    --reaction-border-value-type categorical --reaction-border-slice-by gene \
    --draw-individual-files --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_border_slices_weight_aggregation_genes.log 2>&1
then
    echo "ERROR: a weight aggregation with all slices by gene should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_border_slices_weight_aggregation_genes.log \
    | grep -q "so no weights are reduced to an"
then
    echo "ERROR: a weight aggregation with all slices by gene was refused for another reason."
    exit 1
fi

# A summary of shares needs a drawn map that pools samples and shows slices. Here the rings show no
# border on the 'unified' map, and no map of a group is drawn.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_weights_status.reaction.txt \
    --reaction-slice-weights coverage --reaction-slice-share-summary pooled_sum \
    --reaction-border-column status --reaction-border-value-type categorical \
    --reaction-border-slice-by accession --draw-individual-files \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_border_slices_share_summary.log 2>&1
then
    echo "ERROR: a summary of shares with no map that pools slices should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_border_slices_share_summary.log \
    | grep -q "slices of elements or of rings"
then
    echo "ERROR: a summary of shares with no map that pools slices was refused for another reason."
    exit 1
fi

# Messages of the command line and of the mapper name the options of the border where rings are
# sliced. ConfigError breaks its lines after hyphens, so the spaces and line breaks are deleted
# before the flags are matched.
if ! tr -d ' \n' < draw_txt_border_slices_accession_aggregation.log \
    | grep -q "With'--reaction-border-slice-byaccession',eachsliceofaring"
then
    echo "ERROR: a border accession aggregation for sliced rings named the wrong option."
    exit 1
fi
if ! tr -d ' \n' < draw_txt_border_slices_genes_std.log \
    | grep -q "With'--reaction-border-slice-bygene',itreducestherows"
then
    echo "ERROR: 'std' over the rows of a gene in a ring named the wrong option."
    exit 1
fi
if ! tr -d ' \n' < draw_txt_border_slices_weight_aggregation_genes.log \
    | grep -q "With'--reaction-slice-bygene'and'--reaction-border-slice-bygene',eachgene"
then
    echo "ERROR: a weight aggregation with everything sliced by gene named the wrong options."
    exit 1
fi

# '--reaction-color none' leaves the fills of a file without a value column uncolored, so that its
# boxes show only their rings. Its entries do not count as data, so a map is drawn only where a ring
# shows. The overview map 01200 shows no rings of reactions, so it is not drawn. A KO whose border
# cell is empty shows nothing, as if it were not in the file, and the run counts such rows.
INFO "Testing a reaction layer of rings alone, with an empty border cell"
anvi-draw-kegg-pathways --reaction-txt draw_kos_status_empty_cell.reaction.txt \
    --reaction-border-column status --reaction-border-value-type categorical \
    --reaction-color none \
    --output-dir draw_txt_rings_alone_empty_cell --pathway-numbers $pathway_numbers \
    --no-progress 2>&1 | tee draw_txt_rings_alone_empty_cell.log
if ! tr '\n' ' ' < draw_txt_rings_alone_empty_cell.log \
    | grep -q "Reaction rows without a border value \.*: 1 "
then
    echo "ERROR: the row without a border value of a layer of rings alone was not counted."
    exit 1
fi
if [ ! -s draw_txt_rings_alone_empty_cell/unified/ko00010.pdf ]
then
    echo "ERROR: a map with rings and no colored fills should have been drawn."
    exit 1
fi
if [ -e draw_txt_rings_alone_empty_cell/unified/ko01200.pdf ]
then
    echo "ERROR: a map that shows no ring of a layer of rings alone should not have been drawn."
    exit 1
fi

# Compounds take fills from reactions only where a reaction layer colors fills. A compound layer of
# rings alone without a reaction layer draws white circles with rings, on the overview map 01200 as
# on the standard map 00010.
INFO "Testing a compound layer of rings alone, without a reaction layer"
anvi-draw-kegg-pathways --compound-txt draw_compounds.compound.txt \
    --compound-border-column value --compound-color none \
    --output-dir draw_txt_compound_rings_alone --pathway-numbers $pathway_numbers \
    --no-progress
if [ ! -s draw_txt_compound_rings_alone/unified/ko01200.pdf ]
then
    echo "ERROR: an overview map with compound rings and no reaction layer should have been drawn."
    exit 1
fi
if [ -e draw_txt_compound_rings_alone/colorbar_compounds.pdf ] \
    || [ ! -s draw_txt_compound_rings_alone/colorbar_compounds_border.pdf ]
then
    echo "ERROR: a compound layer of rings alone should have keyed its rings and not its fills."
    exit 1
fi

# With a static color, the map of each group counts the samples of the group and takes a key of the
# counts. A layer of rings alone counts nothing there, so neither the maps of the groups nor their
# panels in the grid take a key of counts. The run needs no '--group-threshold'.
INFO "Testing the maps of groups of a layer of rings alone"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_status.reaction.txt \
    --reaction-border-column status --reaction-border-value-type categorical \
    --reaction-color none --groups-txt draw-sample-group-information.txt \
    --reaction-border-sample-summary plurality_order \
    --reaction-border-group-summary plurality_order \
    --output-dir draw_txt_groups_rings_alone --pathway-numbers $pathway_numbers \
    --draw-individual-files --draw-grid \
    --no-progress
for group in G1 G2
do
    if [ -e draw_txt_groups_rings_alone/individual/$group/colorbar.pdf ] \
        || [ -e draw_txt_groups_rings_alone/grid/colorbar_$group.pdf ]
    then
        echo "ERROR: the maps of group $group took a key of counts for a layer of rings alone."
        exit 1
    fi
    if [ ! -s draw_txt_groups_rings_alone/individual/$group/ko00010.pdf ]
    then
        echo "ERROR: the map of group $group should have been drawn for its rings."
        exit 1
    fi
done

# Without a border sample summary, the 'unified' map shows no ring, and a layer of rings alone shows
# nothing there. The maps of the samples show their rings.
INFO "Testing the warning that the 'unified' map shows nothing of a layer of rings alone"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_status.reaction.txt \
    --reaction-border-column status --reaction-border-value-type categorical \
    --reaction-color none \
    --output-dir draw_txt_rings_alone_no_summary --pathway-numbers $pathway_numbers \
    --draw-individual-files \
    --no-progress 2>&1 | tee draw_txt_rings_alone_no_summary.log
if ! tr '\n' ' ' < draw_txt_rings_alone_no_summary.log \
    | grep -q "The 'unified' map shows nothing of the reaction layer"
then
    echo "ERROR: a 'unified' map that shows nothing of a layer of rings alone was not reported."
    exit 1
fi

INFO "Testing that the colors of layers are refused unless they are hex codes or 'none'"
# '--reaction-color' and '--compound-color' take a six-digit color hex code, as the colors files do.
# A color name, such as 'orange', is refused for every input.
for color_option in "--reaction-txt draw_kos.reaction.txt --reaction-color" \
                    "--compound-txt draw_compounds_presence.compound.txt --compound-color"
do
    if anvi-draw-kegg-pathways $color_option orange --pathway-numbers $pathway_numbers \
        --output-dir draw_txt_bad --overwrite-output-destinations \
        --no-progress > draw_txt_color_name.log 2>&1
    then
        echo "ERROR: '$color_option orange' should have failed but did not."
        exit 1
    fi
    if ! tr '\n' ' ' < draw_txt_color_name.log | grep -q "It takes a six-digit color hex code"
    then
        echo "ERROR: '$color_option orange' was refused for another reason."
        exit 1
    fi
done
for color_input in "--reaction-network-json pan-network.json" \
                   "--contigs-dbs E_faecalis_6240.db" \
                   "--external-genomes external-genomes.txt" \
                   "--pan-db TEST-PAN.db --genomes-storage TEST-GENOMES.db"
do
    if anvi-draw-kegg-pathways $color_input --reaction-color orange \
        --pathway-numbers $pathway_numbers \
        --output-dir draw_db_bad --overwrite-output-destinations \
        --no-progress > draw_db_color_name.log 2>&1
    then
        echo "ERROR: '--reaction-color orange' with '$color_input' should have failed but did not."
        exit 1
    fi
    if ! tr '\n' ' ' < draw_db_color_name.log | grep -q "It takes a six-digit color hex code"
    then
        echo "ERROR: '--reaction-color orange' with '$color_input' was refused for another reason."
        exit 1
    fi
done

# The mapper takes 'original' as its way of asking for the colors of the reference map. Every input
# refuses it as a color, in favor of '--original-color'.
if anvi-draw-kegg-pathways --contigs-dbs E_faecalis_6240.db --reaction-color original \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_db_bad --overwrite-output-destinations \
    --no-progress > draw_db_color_original.log 2>&1
then
    echo "ERROR: '--reaction-color original' on a contigs database run should have failed."
    exit 1
fi
if ! tr '\n' ' ' < draw_db_color_original.log | grep -q "use the '--original-color' flag instead"
then
    echo "ERROR: '--reaction-color original' on a contigs database run failed for another reason."
    exit 1
fi

INFO "Testing that rings alone are refused where they cannot apply"
# Only a text file has a border column, so only a text file can draw rings alone.
if anvi-draw-kegg-pathways --contigs-dbs E_faecalis_6240.db --reaction-color none \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_db_bad --overwrite-output-destinations \
    --no-progress > draw_db_color_none.log 2>&1
then
    echo "ERROR: '--reaction-color none' on a contigs database run should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_db_color_none.log | grep -q "Only a text file given with"
then
    echo "ERROR: '--reaction-color none' on a contigs database run was refused for another reason."
    exit 1
fi

# Without a border column, a layer of rings alone would draw nothing.
if anvi-draw-kegg-pathways --reaction-txt draw_kos.reaction.txt --reaction-color none \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_color_none_no_border.log 2>&1
then
    echo "ERROR: '--reaction-color none' without a border column should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_color_none_no_border.log | grep -q "so the layer would draw nothing"
then
    echo "ERROR: '--reaction-color none' without a border column was refused for another reason."
    exit 1
fi

# The options that color fills have no fills to color. A summary of presence colors the fills of the
# 'unified' map, and the rings take summaries of their own. Slices of the fill need values, and the
# rings can still be sliced with '--reaction-border-slice-by'. '--reaction-reverse-overlay' changes
# only the drawing order of the fills.
for fill_option in "--reaction-sample-summary count" "--reaction-slice-by accession" \
                   "--reaction-reverse-overlay"
do
    if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_status.reaction.txt \
        --reaction-border-column status --reaction-border-value-type categorical \
        --reaction-color none $fill_option --pathway-numbers $pathway_numbers \
        --output-dir draw_txt_bad --overwrite-output-destinations \
        --no-progress > draw_txt_color_none_fill_option.log 2>&1
    then
        echo "ERROR: '$fill_option' for rings alone should have failed but did not."
        exit 1
    fi
    if ! tr '\n' ' ' < draw_txt_color_none_fill_option.log \
        | grep -q "These options color the fills of the reaction layer"
    then
        echo "ERROR: '$fill_option' for rings alone was refused for another reason."
        exit 1
    fi
done

# Global maps draw no rings, so a run of rings alone that draws only global maps would draw nothing.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_status_empty_cell.reaction.txt \
    --reaction-border-column status --reaction-border-value-type categorical \
    --reaction-color none --pathway-numbers 01100 \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_color_none_global.log 2>&1
then
    echo "ERROR: rings alone on the global map alone should have failed but did not."
    exit 1
fi
if ! tr '\n' ' ' < draw_txt_color_none_global.log | grep -q "This run draws no such map"
then
    echo "ERROR: rings alone on the global map alone were refused for another reason."
    exit 1
fi

# The group options style the counts on the maps of groups. A layer of rings alone counts nothing
# there, and no other layer is drawn.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_status.reaction.txt \
    --reaction-border-column status --reaction-border-value-type categorical \
    --reaction-color none --groups-txt draw-sample-group-information.txt \
    --reaction-border-sample-summary plurality_order --group-colormap viridis \
    --draw-individual-files --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_color_none_group_colormap.log 2>&1
then
    echo "ERROR: '--group-colormap' for a grouped layer of rings alone should have failed."
    exit 1
fi
if ! tr -d ' \n' < draw_txt_color_none_group_colormap.log \
    | grep -q "countnothingforthereactionlayer,since'--reaction-colornone'"
then
    echo "ERROR: '--group-colormap' for a grouped layer of rings alone failed for another reason."
    exit 1
fi

# A bar column sets the length of a bar in each reaction box of a standard map, in the color of the
# fill. A box with a fill and no bar value is left uncolored on a map that shows bars, and the run
# counts such boxes.
INFO "Testing a bar column with an empty cell beside a value column"
anvi-draw-kegg-pathways --reaction-txt draw_kos_bar_empty_cell.reaction.txt \
    --reaction-bar-column level --output-dir draw_txt_bar_empty_cell \
    --pathway-numbers $pathway_numbers --no-progress 2>&1 | tee draw_txt_bar_empty_cell.log
for phrase in "Reactionrowswithoutabarvalue" "havenobarvaluefromthe'level'barcolumn"
do
    if ! tr -d ' \n' < draw_txt_bar_empty_cell.log | grep -q "$phrase"
    then
        echo "ERROR: the warning or report that contains '$phrase' was not given."
        exit 1
    fi
done
if [ ! -s draw_txt_bar_empty_cell/bar_key_reactions.pdf ]
then
    echo "ERROR: a run with bars should have written a bar key."
    exit 1
fi

# Bar values of 0 alone draw bars of no length, and a key labeled 0 alone.
INFO "Testing a bar scale whose only value is 0"
anvi-draw-kegg-pathways --reaction-txt draw_kos_bar_zero.reaction.txt --reaction-bar-column level \
    --output-dir draw_txt_bar_zero --pathway-numbers $pathway_numbers --no-progress
if [ ! -s draw_txt_bar_zero/bar_key_reactions.pdf ] \
    || [ ! -s draw_txt_bar_zero/unified/ko00010.pdf ]
then
    echo "ERROR: bars of 0 alone should have been drawn with a key."
    exit 1
fi

# Lines take no bars. A reaction element that a standard map draws as a line is reported.
INFO "Testing the warning that reaction lines on a standard map show no bars"
anvi-draw-kegg-pathways --reaction-txt draw_kos_bar_lines.reaction.txt --reaction-bar-column level \
    --output-dir draw_txt_bar_lines --pathway-numbers 00010 00510 \
    --no-progress 2>&1 | tee draw_txt_bar_lines.log
if ! tr -d ' \n' < draw_txt_bar_lines.log | grep -q "Linesshownobars"
then
    echo "ERROR: reaction lines with bar values on a standard map were not reported."
    exit 1
fi

# A player shows the key of bars only beside the maps that show bars, which are standard maps.
# Skipping the 'unified' map drops its bar key.
INFO "Testing the keys of bars in players and without the 'unified' map"
anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_genes.reaction.txt \
    --reaction-bar-column coverage --skip-unified-maps \
    --output-dir draw_txt_bar_players --pathway-numbers $global_pathway_numbers \
    --draw-individual-files --collate-files-by-map --animate-files-by-map \
    --no-progress
players=draw_txt_bar_players/by_map
standard_keys=$(grep -o '<img class="key"[^>]*' $players/ko00010.html)
overview_keys=$(grep -o '<img class="key"[^>]*' $players/ko01200.html)
if [ "$(echo "$standard_keys" | wc -l)" -ne 2 ] \
    || [ "$overview_keys" != "$(echo "$standard_keys" | sed -n 1p)" ]
then
    echo "ERROR: a player should show the key of bars only beside maps that show bars."
    exit 1
fi
if [ -e draw_txt_bar_players/bar_key_reactions.pdf ] \
    || [ ! -s draw_txt_bar_players/bar_key_reactions_samples.pdf ]
then
    echo "ERROR: skipping the 'unified' map should drop its bar key and keep that of the samples."
    exit 1
fi

INFO "Testing that bar columns are refused where they cannot apply"

# Every bar option needs the bar column.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_genes.reaction.txt \
    --reaction-bar-sample-summary mean --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_bar_option_without_column.log 2>&1
then
    echo "ERROR: a bar option without the bar column should have failed but did not."
    exit 1
fi
if ! tr -d ' \n' < draw_txt_bar_option_without_column.log \
    | grep -q "appliestothebarsofthereactionlayer,butnobarcolumnwasnamed"
then
    echo "ERROR: a bar option without the bar column was refused for another reason."
    exit 1
fi

# A key column cannot hold bar values.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-bar-column sample --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_bar_key_column.log 2>&1
then
    echo "ERROR: a key column as the bar column should have failed but did not."
    exit 1
fi
if ! tr -d ' \n' < draw_txt_bar_key_column.log \
    | grep -q "soitcannotholdbarvalues"
then
    echo "ERROR: a key column as the bar column was refused for another reason."
    exit 1
fi

# The bar column must be a column of the file.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_coverage.reaction.txt \
    --reaction-bar-column nothing --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_bar_missing_column.log 2>&1
then
    echo "ERROR: a missing bar column should have failed but did not."
    exit 1
fi
if ! tr -d ' \n' < draw_txt_bar_missing_column.log \
    | grep -q "asthebarcolumn,butthereactionfilehasnocolumnwiththatheader"
then
    echo "ERROR: a missing bar column was refused for another reason."
    exit 1
fi

# One column cannot set the lengths of bars and color borders.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage_status.reaction.txt \
    --reaction-bar-column status --reaction-border-column status \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_bar_is_border.log 2>&1
then
    echo "ERROR: one column as the bar and border column should have failed but did not."
    exit 1
fi
if ! tr -d ' \n' < draw_txt_bar_is_border.log \
    | grep -q "eithersetthelengthsofbarsorcolorborders"
then
    echo "ERROR: one column as the bar and border column was refused for another reason."
    exit 1
fi

# One column cannot set the lengths of bars and hold the weights of slices.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_weights.reaction.txt \
    --reaction-bar-column coverage --reaction-slice-weights coverage \
    --reaction-border-column expression --reaction-border-slice-by accession \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_bar_is_weights.log 2>&1
then
    echo "ERROR: one column as the bar and weight column should have failed but did not."
    exit 1
fi
if ! tr -d ' \n' < draw_txt_bar_is_weights.log \
    | grep -q "eithersetthelengthsofbarsorholdtheweightsofslices"
then
    echo "ERROR: one column as the bar and weight column was refused for another reason."
    exit 1
fi

# A bar shows values only, so its summaries take no presence names.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_genes.reaction.txt \
    --reaction-bar-column coverage --reaction-bar-sample-summary count \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_bar_presence_summary.log 2>&1
then
    echo "ERROR: a presence summary of bars should have failed but did not."
    exit 1
fi
if ! tr -d ' \n' < draw_txt_bar_presence_summary.log \
    | grep -q "Abarshowsvaluesonly"
then
    echo "ERROR: a presence summary of bars was refused for another reason."
    exit 1
fi

# Bars take numbers only, so their summaries take no rules for categorical values.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_genes.reaction.txt \
    --reaction-bar-column coverage --reaction-bar-sample-summary plurality_order \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_bar_rule_summary.log 2>&1
then
    echo "ERROR: a rule as a bar summary should have failed but did not."
    exit 1
fi
if ! tr -d ' \n' < draw_txt_bar_rule_summary.log \
    | grep -q "Poolthemwithanaggregation"
then
    echo "ERROR: a rule as a bar summary was refused for another reason."
    exit 1
fi

# Bars take no period.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_genes.reaction.txt \
    --reaction-bar-column coverage --reaction-bar-sample-summary circular_mean \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_bar_circular_summary.log 2>&1
then
    echo "ERROR: a circular summary of bars should have failed but did not."
    exit 1
fi
if ! tr -d ' \n' < draw_txt_bar_circular_summary.log \
    | grep -q "takenoperiod,sincethelengthofabarstartsat0"
then
    echo "ERROR: a circular summary of bars was refused for another reason."
    exit 1
fi

# A bar group summary needs a bar sample summary.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_genes.reaction.txt \
    --reaction-bar-column coverage --groups-txt draw-sample-group-information.txt \
    --group-threshold 0.5 --reaction-bar-group-summary mean --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_bar_group_without_sample.log 2>&1
then
    echo "ERROR: a bar group summary without a bar sample summary should have failed but did not."
    exit 1
fi
if ! tr -d ' \n' < draw_txt_bar_group_without_sample.log \
    | grep -q "whichpoolsthebarvaluesofthegroups"
then
    echo "ERROR: a bar group summary without a bar sample summary was refused for another reason."
    exit 1
fi

# A bar scale always includes 0.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_genes.reaction.txt \
    --reaction-bar-column coverage --reaction-bar-sample-summary mean \
    --reaction-bar-value-limits 10 100 --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_bar_limits_without_zero.log 2>&1
then
    echo "ERROR: bar limits that leave out 0 should have failed but did not."
    exit 1
fi
if ! tr -d ' \n' < draw_txt_bar_limits_without_zero.log \
    | grep -q "limitsthatleaveout0"
then
    echo "ERROR: bar limits that leave out 0 was refused for another reason."
    exit 1
fi

# Bars that grow both ways from 0 are not supported yet.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_genes.reaction.txt \
    --reaction-bar-column coverage --reaction-bar-sample-summary mean \
    --reaction-bar-value-limits -5 100 --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_bar_limits_negative.log 2>&1
then
    echo "ERROR: bar limits below 0 should have failed but did not."
    exit 1
fi
if ! tr -d ' \n' < draw_txt_bar_limits_negative.log \
    | grep -q "limitsthatreachbelow0"
then
    echo "ERROR: bar limits below 0 was refused for another reason."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_bad_negative_bar.reaction.txt \
    --reaction-bar-column level --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_bar_negative_values.log 2>&1
then
    echo "ERROR: a negative bar value should have failed but did not."
    exit 1
fi
if ! tr -d ' \n' < draw_txt_bar_negative_values.log \
    | grep -q "mustbe0ormore,sincebarsthatgrowbothwaysfrom0arenotsupportedyet"
then
    echo "ERROR: a negative bar value was refused for another reason."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_genes.reaction.txt \
    --reaction-bar-column coverage --reaction-bar-sample-summary skew \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_bar_aggregation_negative.log 2>&1
then
    echo "ERROR: bar values reduced below 0 should have failed but did not."
    exit 1
fi
if ! tr -d ' \n' < draw_txt_bar_aggregation_negative.log \
    | grep -q "reduceitsvaluestonumbersbelow0"
then
    echo "ERROR: bar values reduced below 0 was refused for another reason."
    exit 1
fi

# Each bar value must be a number, and a bar column needs a value.
if anvi-draw-kegg-pathways --reaction-txt draw_bad_na_bar.reaction.txt --reaction-bar-column level \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_bar_na_value.log 2>&1
then
    echo "ERROR: 'NA' as a bar value should have failed but did not."
    exit 1
fi
if ! tr -d ' \n' < draw_txt_bar_na_value.log \
    | grep -q "somethingotherthanafinitenumber"
then
    echo "ERROR: 'NA' as a bar value was refused for another reason."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_bad_empty_bar.reaction.txt \
    --reaction-bar-column level --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_bar_empty_column.log 2>&1
then
    echo "ERROR: a bar column with no values should have failed but did not."
    exit 1
fi
if ! tr -d ' \n' < draw_txt_bar_empty_column.log \
    | grep -q "itwouldsetnobar"
then
    echo "ERROR: a bar column with no values was refused for another reason."
    exit 1
fi

# A bar takes the color of its fill, which the reference map's colors and 'none' do not give.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_bar_zero.reaction.txt \
    --reaction-bar-column level --original-color --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_bar_original_color.log 2>&1
then
    echo "ERROR: bars with '--original-color' should have failed but did not."
    exit 1
fi
if ! tr -d ' \n' < draw_txt_bar_original_color.log \
    | grep -q "givenalongwith'--original-color':'--reaction-bar-column'"
then
    echo "ERROR: bars with '--original-color' was refused for another reason."
    exit 1
fi

if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_coverage_status.reaction.txt \
    --reaction-bar-column coverage --reaction-border-column status \
    --reaction-border-value-type categorical --reaction-color none \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_bar_color_none.log 2>&1
then
    echo "ERROR: bars with '--reaction-color none' should have failed but did not."
    exit 1
fi
if ! tr -d ' \n' < draw_txt_bar_color_none.log \
    | grep -q "sothebarswouldhavenocolor"
then
    echo "ERROR: bars with '--reaction-color none' was refused for another reason."
    exit 1
fi

# Bars are not drawn in slices yet.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_genes.reaction.txt \
    --reaction-bar-column coverage --reaction-slice-by accession \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_bar_slice_by.log 2>&1
then
    echo "ERROR: bars with '--reaction-slice-by' should have failed but did not."
    exit 1
fi
if ! tr -d ' \n' < draw_txt_bar_slice_by.log \
    | grep -q "cannotbecombinedwith'--reaction-bar-column'yet"
then
    echo "ERROR: bars with '--reaction-slice-by' was refused for another reason."
    exit 1
fi

# Only a text file has a bar column.
if anvi-draw-kegg-pathways --contigs-dbs E_faecalis_6240.db --reaction-bar-column coverage \
    --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_bar_contigs_db.log 2>&1
then
    echo "ERROR: a bar column for a contigs database should have failed but did not."
    exit 1
fi
if ! tr -d ' \n' < draw_txt_bar_contigs_db.log \
    | grep -q "Theseoptionsweregiven:'--reaction-bar-column'"
then
    echo "ERROR: a bar column for a contigs database was refused for another reason."
    exit 1
fi

# Bars are drawn only in the reaction boxes of standard maps.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_genes.reaction.txt \
    --reaction-bar-column coverage --pathway-numbers 01200 \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_bar_only_overview.log 2>&1
then
    echo "ERROR: bars on an overview map alone should have failed but did not."
    exit 1
fi
if ! tr -d ' \n' < draw_txt_bar_only_overview.log \
    | grep -q "drawsbarsinreactionboxes,whichonlystandardmapsdraw"
then
    echo "ERROR: bars on an overview map alone was refused for another reason."
    exit 1
fi

# Without a bar sample summary, the 'unified' map of a file with samples shows no bars.
if anvi-draw-kegg-pathways --reaction-txt draw_kos_samples_genes.reaction.txt \
    --reaction-bar-column coverage --pathway-numbers $pathway_numbers \
    --output-dir draw_txt_bad --overwrite-output-destinations \
    --no-progress > draw_txt_bar_unified_without_summary.log 2>&1
then
    echo "ERROR: bars that no map shows should have failed but did not."
    exit 1
fi
if ! tr -d ' \n' < draw_txt_bar_unified_without_summary.log \
    | grep -q "butnomapofthisrunwouldshowitsbars"
then
    echo "ERROR: bars that no map shows was refused for another reason."
    exit 1
fi

# Colors given per group of contigs databases, with each group's own maps colored by a ramp running
# to that group's color rather than by one colormap shared across the groups.
INFO "Testing colors given per group of contigs databases"
anvi-draw-kegg-pathways --contigs-dbs E_faecalis_6240.db E_faecalis_6255.db E_faecalis_6512.db E_faecalis_6557.db E_faecalis_6563.db \
                        --groups-txt contigs-db-group-information.txt \
                        --group-threshold 0 \
                        --reaction-membership-colors draw-group-colors.txt \
                        --group-colormap tint \
                        --output-dir contigs_dbs_kos_group_tint_ramps \
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
