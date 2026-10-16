#!/usr/bin/env python
DESCRIPTION = "Write KEGG pathway map files incorporating data sourced from external inputs"

import os
import re
import sys
import functools
import importlib
import pandas as pd

from argparse import Namespace

import anvio.filesnpaths as filesnpaths

from anvio.argparse import ArgumentParser
from anvio import A, K, __version__ as VERSION
from anvio.metabolism.context import KeggContext
from anvio.errors import ConfigError, FilesNPathsError


def check_package_dependencies():
    missing_packages = []
    for package_name in ('Bio', 'reportlab', 'fitz'):
        try:
            importlib.import_module(package_name)
        except ImportError:
            missing_packages.append(package_name)

    if missing_packages:
        message = ', '.join(f"'{package_name}'" for package_name in missing_packages)
        # Printed rather than raised: this runs at import time, where there is no handler to turn
        # a ConfigError into the message alone rather than into a traceback.
        print(ConfigError(
            f"The following Python packages required to run `anvi-draw-kegg-pathways` could not be "
            f"imported: {message}. All Python dependencies of anvi'o are described in its "
            f"`pyproject.toml`, and can be installed by running the command `pip install -e .` in "
            f"the top directory of the anvi'o codebase, which may have a location like "
            f"'~/github/anvio' in your file system if you followed the anvio.org installation "
            f"instructions."
        ))
        sys.exit(-1)

# Checked before the import below rather than from within main: 'anvio.keggmapping' imports 'fitz'
# and, through 'anvio.kgml', 'reportlab', so by the time main runs a missing one of those has
# already reached the user as an ImportError traceback.
check_package_dependencies()

from anvio.keggmapping import (
    AGGREGATION_FUNCTIONS, CATEGORICAL_RULES, CATEGORICAL_RULES_PHRASE, DEFAULT_CATEGORICAL_SEED,
    DEFAULT_CENTERED_COLORMAP, DEFAULT_GROUP_TINT_SPAN, DEFAULT_VALUE_COLORS,
    ELEMENT_NORMALIZATION_PHRASE, GROUP_COLORMAP_TINT, GROUP_SCHEME_OPTIONS,
    MAX_DISCRETE_COUNT_BANDS, MIN_ANGULAR_DEVIATION_VALUES, MIN_BOX_SLICE_WIDTH,
    SLICE_SHARE_SUMMARIES, SLICE_UNITS, SLICE_WEIGHT_AGGREGATIONS, SUMMARY_PRESENCE_PHRASE,
    SUMMARY_PRESENCE_SCHEMES, VALUE_TYPES, Mapper
)


__authors__ = ["semiller10"]
__copyright__ = "Copyleft 2015-2026, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = VERSION
__requires__ = ['kegg-data']
__can_use__ = ['contigs-db', 'external-genomes', 'pan-db', 'genomes-storage-db',
               'reaction-network-json', 'kegg-reaction-txt', 'kegg-compound-txt',
               'kegg-membership-colors-txt', 'kegg-value-colors-txt']
__provides__ = ['kegg-pathway-map', 'kegg-slice-key-txt']
__description__ = DESCRIPTION

# Default presence/absence colors, matching the Mapper method defaults. Used as the
# '--reaction-color' and '--compound-color' values when those flags are given without an explicit
# hex code.
DEFAULT_REACTION_COLOR = '#2ca02c'
DEFAULT_COMPOUND_COLOR = '#e239af'

# The '--*-sample-summary'/'--*-group-summary' values that summarize presence across samples or
# groups, taken from the mapper's schemes for coloring them so that the two cannot drift. Any other
# value of those options, and any value of the '--*-aggregation' options, names an aggregation that
# pools a layer's value column, so it requires the layer to have one. These names take precedence
# over the aggregations, so pandas 'count' (which counts rows) is not reachable through a summary
# option.
PRESENCE_SUMMARIES = tuple(SUMMARY_PRESENCE_SCHEMES)

# The '--count-scale-max' values naming a rule rather than giving a count.
COUNT_SCALE_MAX_RULES = ('observed', 'total')

# The word standing in for a number at either end of the '--*-value-limits' options, leaving that
# end of a color scale wherever the data puts it. Both options take two values, so an end can be
# left open without being left out.
VALUE_LIMIT_NONE = 'none'

# The value the '--*-value-center' options center a color scale on when they are given as a bare
# flag. Zero is what a signed quantity is read against — a log ratio, a difference of means — so it
# is what a bare flag is nearly always meant to ask for. It stays a string, as a value typed after
# the option would be, so that 'Mapper._resolve_value_center' converts and checks every center in
# one place.
VALUE_CENTER_DEFAULT = '0'

# Aggregation names with fast paths in the mapper, listed in help text as the recommended ones. Any
# other pandas aggregation reducing values to one number is accepted too ('_resolve_aggregation').
RECOMMENDED_AGGREGATIONS = ', '.join(f"'{name}'" for name in AGGREGATION_FUNCTIONS)

# The modes of a text layer with a value column ('_detect_txt_layer'). Its values are numbers or
# categorical values.
VALUE_MODES = ('quantitative', 'categorical')


def get_args() -> Namespace:
    parser = ArgumentParser(description=DESCRIPTION)
    # Disable prefix abbreviation: with the many symmetric '--reaction-*'/'--compound-*'/'--group-*'
    # options, abbreviations are mostly ambiguous anyway, and turning them off keeps a removed or
    # mistyped flag a clear "unrecognized arguments" error rather than a silent partial match.
    parser.allow_abbrev = False

    groupTXT = parser.add_argument_group(
        "INPUT TEXT FILES",
        "Draw pathway maps directly from tab-delimited text files, coloring reaction and/or "
        "compound map elements. Provide a reaction-layer file, a compound-layer file, or both of "
        "these draw-kegg-pathway text files; each layer is colored independently by its own "
        "auto-detected mode."
    )
    groupTXT.add_argument(
        '--reaction-txt', type=str, metavar='FILE', help=
        "Path to a tab-delimited kegg-reaction-txt file coloring the reaction elements of maps. It "
        "must have an 'accession' column whose values are ALL KO IDs ('K' followed by digits) or "
        "ALL KEGG reaction IDs ('R' followed by digits): the two cannot be mixed in one file. An "
        "optional 'gene_id' column names the gene that a row comes from. Gene IDs affect coloring "
        "only with '--reaction-slice-by gene', which draws each gene as a slice. Otherwise, the "
        "values of the genes of an accession are reduced to one value "
        "('--reaction-gene-aggregation'). An optional 'sample' column assigns each row to a sample "
        "of origin, enabling comparison across samples (and, with '--groups-txt', groups of "
        "samples) as when using multiple contigs databases. The single column that is not "
        "'accession', 'gene_id', 'sample', or a column of weights ('--reaction-slice-weights') is "
        "auto-detected as the value column. With it, reactions are colored by value, and without "
        "it, by presence. Its values are numbers, or categorical values with "
        "'--reaction-value-type categorical'."
    )
    groupTXT.add_argument(
        '--compound-txt', type=str, metavar='FILE', help=
        "Path to a tab-delimited kegg-compound-txt file coloring the compound elements of maps. It "
        "must have an 'accession' column whose values are all KEGG compound IDs ('C' followed by "
        "digits). It must NOT have a 'gene_id' column, since compounds are not gene products. An "
        "optional 'sample' column works as for '--reaction-txt'. The single column that is not "
        "'accession' or 'sample' is auto-detected as the value column: with it, compounds are "
        "colored by value; without it, by presence. Its values are numbers, or categorical values "
        "with '--compound-value-type categorical'."
    )
    groupTXT.add_argument(
        '--reaction-gene-aggregation', metavar='NAME', help=
        f"How to reduce the values of the GENES annotated with one accession to that accession's "
        f"own value, when the reaction layer is colored by a value column and its file names genes "
        f"in a 'gene_id' column. Verified names are {RECOMMENDED_AGGREGATIONS}, and the default is "
        f"'sum': the genes carrying a KO contribute to that KO's total. Any other pandas "
        f"aggregation that reduces a series of numbers to a single number works too, such as 'var' "
        f"or 'sem'; a name that transforms rather than reduces, like 'cumsum', or that only a "
        f"grouping offers, like 'first', is rejected. Note that two rows describing the same gene "
        f"and accession are refused rather than combined, so this reduces genuinely different "
        f"genes. Where an aggregation is undefined for the values available, as 'std' is for a "
        f"single value, the element is left uncolored and a warning reports how many accessions "
        f"were affected. With '--reaction-value-period', the default is 'circular_mean', and no "
        f"other name is accepted. With '--reaction-value-type categorical', this takes a rule that "
        f"chooses one value among the genes instead: {CATEGORICAL_RULES_PHRASE}. The default is "
        f"'plurality_order' (see '--reaction-value-type'). With '--reaction-slice-by gene', each "
        f"gene is a slice of its own, so this option is refused."
    )
    groupTXT.add_argument(
        '--reaction-accession-aggregation', metavar='NAME', help=
        "How to reduce the values of the several accessions that ONE MAP ELEMENT stands for to "
        "that element's color, when the reaction layer is colored by a value column. Nearly a "
        "quarter of the reaction elements in KEGG's maps stand for more than one KO, so this is a "
        "common case. The default is 'sum': the KOs of a reaction contribute to that reaction's "
        "total. It takes the same names as '--reaction-gene-aggregation', and the two can differ: "
        "averaging across the genes of each KO while totaling across the KOs of a reaction, say. "
        "With '--reaction-value-period', the default is 'circular_mean', and no other name is "
        "accepted. With '--reaction-value-type categorical', this takes a rule that chooses one "
        "value among the accessions of a map element instead, with the same names and default as "
        "'--reaction-gene-aggregation'. With '--reaction-slice-by gene', this reduces the rows of "
        "one gene for the accessions of an element, and its default for numbers is 'mean' rather "
        "than 'sum'. A gene with one such row shows the value of that row, so a name that does not "
        "return a single value unchanged, such as 'std' or 'count', is refused. With "
        "'--reaction-slice-by accession', it is refused for numbers, since each slice shows the "
        "value of one accession."
    )
    groupTXT.add_argument(
        '--compound-accession-aggregation', metavar='NAME', help=
        "Like '--reaction-accession-aggregation', but for the compound layer, and with the same "
        "default of 'sum': a single circle on a map occasionally stands for more than one KEGG "
        "compound, such as the circle covering four forms of cobalamin, and this reduces the "
        "values of whichever of them the file lists to that circle's color. It applies to about "
        "one compound circle in fifty, so it rarely has anything to do. There is no reduction "
        "below it: a compound file names no genes, and rows repeating a compound are refused "
        "rather than combined, so each compound has one value per sample. With "
        "'--compound-value-period', the default is 'circular_mean', and no other name is accepted. "
        "With '--compound-value-type categorical', this takes a rule that chooses one value among "
        "the compounds of a circle instead, with the same names and default as "
        "'--reaction-gene-aggregation'."
    )
    groupTXT.add_argument(
        '--reaction-value-period', metavar='PERIOD', help=
        f"The period after which the values of the reaction layer repeat, in the units of its "
        f"value column. For clock times in hours, give 24. Values just below the period and just "
        f"above 0 are then close together — with a period of 24, the circular mean of 23.5 and 0.5 "
        f"is 0, not 12. Every reduction of the values averages them on a circle with "
        f"'circular_mean'. It is the default of '--reaction-gene-aggregation' and "
        f"'--reaction-accession-aggregation'. It is also the only name that these options accept, "
        f"and the only average that the sample and group summaries accept. Values that cancel out "
        f"have no mean, such as 6 and 18; elements with undefined means are left uncolored. A "
        f"period does not affect presence summaries. With this argument, "
        f"'--reaction-element-normalization' takes only 'difference_from_circular_mean', since the "
        f"other normalizations compare values on a line. Each color scale of the values themselves "
        f"runs from 0 to the period, so a color always means the same point in the period. Limits "
        f"other than 0 and the period are refused, and so is a center. Without a colormap, the "
        f"scales use the cyclic 'clocktime'. On global and overview maps drawn from reactions "
        f"alone, compound colors are also averaged on a circle. The summary that colors the "
        f"'unified' map ('--reaction-sample-summary' without groups and '--reaction-group-summary' "
        f"with groups) can instead be 'angular_deviation', with each element colored by how far "
        f"its values in the samples or groups are spread apart. With a period of 24, the scale "
        f"ranges from 0 to 5.4. This summary requires values from at least "
        f"{MIN_ANGULAR_DEVIATION_VALUES} samples or groups, and its scale takes a sequential "
        f"colormap ('plasma_r' unless '--reaction-colormap' names another)."
    )
    groupTXT.add_argument(
        '--compound-value-period', metavar='PERIOD', help=
        "Like '--reaction-value-period', but for the compound layer."
    )
    groupTXT.add_argument(
        '--reaction-value-type', choices=list(VALUE_TYPES), help=
        f"What the value column of the reaction file holds. 'numeric', the default, means numbers. "
        f"They are colored along a scale. 'categorical' means named values, such as 'up', 'down' "
        f"and 'same'. Each value takes a color of its own. Values written as numbers, such as 1 "
        f"and 2, are names too. A value is read as text, so 'NA' is a value like any other. Every "
        f"row needs a value. Each accession can have one value per sample, or one per gene and "
        f"sample. The values take {len(DEFAULT_VALUE_COLORS)} default colors in sorted order. Up "
        f"to 9 values take the colors of Matplotlib's 'tab10' without its gray. The 10th to 18th "
        f"values take the light colors of 'tab20' without its light gray. "
        f"'--reaction-value-colors' gives colors of your choosing instead. Where several values "
        f"meet, a rule chooses one of them. '--reaction-gene-aggregation' chooses among the genes "
        f"of an accession, and '--reaction-accession-aggregation' among the KOs or reactions of a "
        f"map element. The rules are {CATEGORICAL_RULES_PHRASE}. 'plurality_order', the default, "
        f"takes the most common value, even with less than half of the votes. It breaks a tie by "
        f"the order of the values. 'plurality_random' breaks a tie at random. 'plurality' gives no "
        f"value where there is a tie. 'majority' takes the value with more than half of the votes. "
        f"Without one, it gives no value. 'order' takes the value that comes first in the order of "
        f"the values. 'random' takes one of the values at random. An accession given no value does "
        f"not vote in its map elements. An element given no value is left uncolored. It does not "
        f"vote across the samples or groups. A warning says how many elements have no value. The "
        f"order is the row order of '--reaction-value-colors', or else sorted order. Where every "
        f"value is a number, the values are sorted as numbers. Otherwise the digits in the values "
        f"are compared as numbers, so 'cluster 2' comes before 'cluster 10'. A random choice "
        f"depends on '--seed'. The map of each sample shows that sample's values. A value in one "
        f"sample need not mean the same thing in another sample, so by default the 'unified' map "
        f"counts the samples that contain each element. A rule given to "
        f"'--reaction-sample-summary' instead chooses one value across the samples. On global and "
        f"overview maps drawn from reactions alone, each compound takes the value that the "
        f"accession rule chooses among the reactions that touch it. A compound given no value "
        f"keeps the color of a compound that no colored reaction touches. These options for a "
        f"scale of numbers are refused: '--reaction-value-limits', "
        f"'--reaction-individual-value-limits', '--reaction-value-center', "
        f"'--reaction-individual-value-center', '--reaction-value-period', "
        f"'--reaction-element-normalization', and '--reaction-individual-colormap'. With "
        f"'--groups-txt', each group's map counts its samples, and the 'unified' map counts the "
        f"groups. A rule in '--reaction-sample-summary' colors each group's map instead. "
        f"'--reaction-group-summary' can then give a rule for the 'unified' map. A group given no "
        f"value does not vote for the 'unified' map."
    )
    groupTXT.add_argument(
        '--compound-value-type', choices=list(VALUE_TYPES), help=
        "Like '--reaction-value-type', but for the value column of the compound file. A compound "
        "file has no genes, so '--compound-accession-aggregation' is the only rule below the "
        "samples."
    )
    groupTXT.add_argument(
        '--seed', type=int, metavar='INT', help=
        f"The seed of random choices among categorical values. The rules 'plurality_random' and "
        f"'random' make these choices (see '--reaction-value-type'). The default is "
        f"{DEFAULT_CATEGORICAL_SEED}. Each choice comes from the seed and from what the value is "
        f"chosen for, such as the KOs of a map element in one sample. A map element therefore gets "
        f"the same value on every map, and in every run with the same seed. This option is refused "
        f"where no rule chooses at random."
    )
    groupTXT.add_argument(
        '--reaction-slice-by', choices=list(SLICE_UNITS), help=
        f"Draw each reaction element in slices, one for each of its accessions or genes with a "
        f"value on the map. Each slice takes the color of its own value. This needs a "
        f"'--reaction-txt' file with a value column. Maps that show values are sliced, but maps "
        f"that show presence are drawn whole. With 'accession', a KO file is sliced by KO, and a "
        f"reaction file by KEGG reaction. Slices follow the order of the accessions in the map's "
        f"KGML file. With 'gene', an element is sliced by the genes of a 'gene_id' column, and "
        f"every row needs a gene ID. A gene with rows for two or more accessions of an element is "
        f"one slice. In each sample, the values of these rows are reduced to their mean, or to "
        f"their circular mean with '--reaction-value-period'. Without a period, "
        f"'--reaction-accession-aggregation' can name another aggregation that returns a single "
        f"value unchanged, such as 'sum' or 'max'. For categorical values, its rule chooses among "
        f"the rows. Slices of genes are ordered by the first of their accessions in the KGML file, "
        f"and then by gene ID, with digits compared as numbers. A gene ID must name the same gene "
        f"in every sample, since the maps of groups and the 'unified' map pool the values of each "
        f"gene across samples. Gene caller IDs start at 0 in every contigs database, so IDs from "
        f"separate assemblies must be made unique, for example with a sample prefix. A warning "
        f"names the gene IDs that have no accession in common between their rows in two samples of "
        f"a map that pools samples. Such an ID likely names two genes. Another warning names the "
        f"maps that pool samples where no gene ID occurs in two or more samples. Each gene is a "
        f"slice of its own, so '--reaction-gene-aggregation' and "
        f"'--reaction-slice-weight-aggregation' are refused. A box is cut into columns from left "
        f"to right, and a line into repeating dashes. The dashes fill the part of the line that "
        f"compound circles leave visible. The slices of an element are equal unless "
        f"'--reaction-slice-weights' names a column of weights. A box is drawn whole where an "
        f"equal slice would be narrower than {MIN_BOX_SLICE_WIDTH} points, and a line where an "
        f"equal slice of its visible part would be shorter than the line is wide. An element drawn "
        f"whole takes the mean of its slice values, weighted by their shares. With "
        f"'--reaction-value-period', this is a circular mean. With '--reaction-value-type "
        f"categorical', the element takes the value that '--reaction-accession-aggregation' "
        f"chooses among the accessions or genes of its slices. Each slice shows the value of one "
        f"accession or gene, so the color scales span the values of single accessions or genes. "
        f"When numbers are sliced by accession, the weighted mean takes the place of "
        f"'--reaction-accession-aggregation', which is refused. An element keeps its usual label, "
        f"and its slices are not labeled. The file 'slice_key_reactions.txt' in the output "
        f"directory lists each slice on each map, and each element drawn whole."
    )
    groupTXT.add_argument(
        '--reaction-slice-weights', metavar='COLUMN', help=
        "The column of the reaction file that holds the weights of slices ('--reaction-slice-by'). "
        "Without this option, the slices of an element are equal. This column is never taken for "
        "the value column. Every row needs a weight, which is a number of 0 or more. Where a map "
        "shows one sample, or the file has no 'sample' column, the share of a slice is its weight "
        "divided by the total weight of the slices drawn. An accession or gene with a weight of 0 "
        "gets no slice. When slicing by accession with a 'gene_id' column, the weights of the "
        "genes of an accession are added up ('--reaction-slice-weight-aggregation'). When slicing "
        "by gene, a gene takes the mean weight of its rows for the accessions of an element, so a "
        "gene with two KOs does not take twice the space. Group maps and the 'unified' map pool "
        "the shares of their samples ('--reaction-slice-share-summary'). The width of a slice in a "
        "box, or its length on a line, shows its share of the element, not its weight. The weights "
        "do not change the colors of slices."
    )
    groupTXT.add_argument(
        '--reaction-slice-weight-aggregation', choices=list(SLICE_WEIGHT_AGGREGATIONS), help=
        "How to reduce the weights of the genes annotated with one accession to the weight of that "
        "accession. This needs '--reaction-slice-by accession', '--reaction-slice-weights', and a "
        "'gene_id' column. 'sum', the default, suits amounts that add up over genes, such as "
        "counts and coverage. 'mean' suits rates and scores of single genes. With 'mean', a KO "
        "with ten genes takes the same share as a KO with one gene, if all of their genes have the "
        "same weight."
    )
    groupTXT.add_argument(
        '--reaction-slice-share-summary', choices=list(SLICE_SHARE_SUMMARIES), help=
        "How group maps and the 'unified' map pool the shares of weighted slices across samples "
        "('--reaction-slice-weights'). The map of a single sample takes its shares from the "
        "weights of that sample. 'share_mean', the default, finds the share of each accession in "
        "each sample and averages these shares over the samples in which the element has weight. "
        "An accession absent from a sample has a share of 0 there. Each sample counts once, so the "
        "weights need to be comparable only within a sample. 'pooled_sum' adds up the weights of "
        "each accession over the samples and divides each sum by the total over the accessions of "
        "the element. A sample with more weight counts for more, so use it where the weights are "
        "in the same units in every sample, such as counts at comparable sequencing depth. Slices "
        "of genes are pooled by the same rules, with genes in place of accessions. With "
        "'--groups-txt', the map of each group pools its samples, and the 'unified' map pools the "
        "groups by the same rule. Under 'share_mean', each group counts once. This option needs "
        "'--reaction-slice-weights', since slices without weights are equal on every map."
    )

    groupJSON = parser.add_argument_group(
        "REACTION NETWORK JSON",
        "Display KO data from a metabolic model JSON file produced by 'anvi-reaction-network "
        "--enzymes-txt' or 'anvi-get-metabolic-model-file'."
    )
    groupJSON.add_argument(
        '--reaction-network-json', type=str, metavar='FILE', help=
        "Path to a metabolic model JSON file (produced by 'anvi-reaction-network --enzymes-txt' or "
        "'anvi-get-metabolic-model-file'). KO IDs are extracted from the gene annotations in the "
        "JSON and used to highlight reactions in pathway maps."
    )

    groupCONTIGS = parser.add_argument_group(
        "CONTIGS DATABASE",
        "Display KO data from one or more contigs databases, e.g., for genomes and metagenomes."
    )
    groupCONTIGS.add_argument(
        '--contigs-dbs', type=str, nargs='+', help=
        "One or more anvi'o contigs databases generated by 'anvi-gen-contigs-database'. Contigs "
        "databases can alternatively be provided by '--external-genomes'."
    )
    groupCONTIGS.add_argument(
        '--external-genomes', type=str, help=
        "A two-column tab-delimited text file that lists anvi'o contigs databases generated by "
        "'anvi-gen-contigs-database'. Contigs databases can alternatively be provided by "
        "'--contigs-dbs'. The header row must be 'name' and 'contigs_db_path' in that order. Each "
        "line in the file should describe a single entry, where the first column contains a "
        "database name, and the second contains the path to the database."
    )

    groupPAN = parser.add_argument_group("PANGENOMIC DATABASE", "Display KO data from a pangenome.")
    groupPAN.add_argument(*A('pan-db'), **K('pan-db', {'required': False}))
    groupPAN.add_argument(*A('genomes-storage'), **K('genomes-storage', {'required': False}))
    groupPAN.add_argument(
        '--consensus-threshold', default=None, type=float, metavar='FLOAT', help=
        "If this argument is provided, then a KO annotation must be assigned to this minimum "
        "proportion of genes in a cluster to be imputed to the cluster as a whole. By default, "
        "without this argument, the annotation assigned to the most genes becomes the annotation "
        "of the cluster (see '--discard-ties'). The consensus threshold must be on [0, 1]."
    )
    groupPAN.add_argument(
        '--discard-ties', action='store_true', default=False, help=
        "By default, a gene cluster is assigned a KO annotation by finding the KO that occurs in "
        "the greatest number of genes in the cluster (see '--consensus-threshold') and arbitrarily "
        "choosing one KO in case of a tie. With this flag, a tie instead results in a KO "
        "annotation not being assigned to the cluster."
    )

    groupGROUP = parser.add_argument_group(
        "GROUPS",
        "Display data from groups of contigs databases, genomes in a pangenome, or samples in a "
        "draw-kegg-pathways text file."
    )
    groupGROUP.add_argument(
        *A('groups-txt'), **K('groups-txt', {'help':
        "Path to a tab-delimited text file specifying which group each item belongs to. Items may "
        "be: contigs database files if provided with '--contigs-dbs', the names of genomes in a "
        "pangenomic database if provided with '--pan-db', or the names of samples in the 'sample' "
        "column of a draw-kegg-pathways text file if provided with "
        "'--reaction-txt'/'--compound-txt'. The first column can have any header, e.g., "
        "'contigs_db', 'pangenome', 'sample', or 'source'. The second column, headed 'group', must "
        "contain the group name for each item. Each item can only be associated with a single "
        "group. It is recommended that group names be single words without fancy characters, such "
        "as 'HIGH_TEMPERATURE' or 'LOW_FITNESS' rather than 'my group #1' or 'IS-THIS-OK?'. With "
        "groups, maps are drawn for individual groups rather than individual items "
        "('--draw-individual-files', '--draw-grid'), and the 'unified' map summarizes the groups "
        "rather than the items. For a draw-kegg-pathways text file, how each group's samples and "
        "how the groups themselves are summarized are set by "
        "'--reaction-sample-summary'/'--compound-sample-summary' and "
        "'--reaction-group-summary'/'--compound-group-summary'. Whenever the groups are summarized "
        "by presence, which is the default, '--group-threshold' must also be given to define when "
        "a group contains a map element, unless '--skip-unified-maps' is also used, in which case "
        "'--group-threshold' is invalid, as it only applies to the 'unified' map."
    }))
    groupGROUP.add_argument(
        '--group-threshold', type=float, metavar='FLOAT', help=
        "The proportion of items in a group containing data of interest for it to be represented "
        "in terms of presence/absence in a map feature. Here is a concrete example. Say the "
        "'groups-txt' file defines three groups of genome items in a pangenome representing "
        "different species, 'A', 'B', and 'C'. You wish to understand the distribution of "
        "metabolic capabilities across the 3 species using KO annotations of genes. Reaction "
        "colors are assigned based on the groups rather than individual genomes containing the "
        "reaction. Thresholds on [0, 1] can be set to define group membership: a threshold of '0' "
        "would mean that ANY genome in the group can contain the reaction via KOs for the reaction "
        "to be considered present in the group; a threshold of '0.75' means at least 75%% of the "
        "genomes in the group must contain the reaction for it to be present; a threshold of '1' "
        "means that ALL genomes in the group must contain the reaction for it to be present. In "
        "our example, set the threshold to '0.5'. Reaction J on a map corresponds to KO X, and "
        "Reaction K on a map corresponds to KOs Y and Z. 90%% of species A genomes, 50%% of "
        "species B genomes, and 10%% of species C genomes contain KO X, so Reaction J would be "
        "colored to indicate that it is represented in species A and B. 0%% of species A genomes, "
        "15%% of species B genomes, and 40%% of species C genomes contain KO Y or KO Z, so "
        "Reaction K would not be colored. This threshold defines group presence, so it applies "
        "whenever groups are summarized by presence: always for contigs databases and pangenomes, "
        "and for a draw-kegg-pathways text layer whose "
        "'--reaction-group-summary'/'--compound-group-summary' summarizes presence ('count', "
        "'count_continuous', or 'membership'), which is the default. It does not apply when every "
        "text layer's groups are summarized by pooling values, or by choosing a categorical value "
        "with a rule, instead. It also does not apply with '--skip-unified-maps', since only the "
        "'unified' map uses the threshold, while a group's own map counts the group's own sources "
        "whatever the threshold is."
    )

    groupSUMMARY = parser.add_argument_group(
        "SAMPLE AND GROUP SUMMARIES",
        "How a draw-kegg-pathways text layer with a 'sample' column is summarized across its "
        "samples, and across groups of samples defined by '--groups-txt'. These are the two "
        "reductions above the within-sample gene and accession aggregations. Each can "
        "independently summarize presence (in how many, or exactly which, samples or groups an "
        "element occurs) or pool the values of a value column. For categorical values, each can "
        "choose one value with a rule instead."
    )
    groupSUMMARY.add_argument(
        '--reaction-sample-summary', metavar='NAME', help=
        f"How to summarize a SET OF SAMPLES for the reaction layer. {SUMMARY_PRESENCE_PHRASE} "
        f"summarize presence, coloring an element by how many ('count' and 'count_continuous') or "
        f"exactly which ('membership') samples contain it. 'count' gives each count its own band "
        f"of a discrete colorbar, which needs a distinguishable color per sample and a legible "
        f"label per band, and so is refused when there are more samples than the colormap can "
        f"supply and warns when there are more bands than a bar can label; 'count_continuous' "
        f"colors counts from the same colormap — assigning the very same colors, with the "
        f"sequential colormap a count scale calls for — but draws the scale as a gradient from the "
        f"lowest count to the highest, which any number of samples can share, at the cost of a "
        f"color reading as a position along the range rather than as an exact count. Any other "
        f"name is an aggregation that pools the samples' values from the layer's value column, "
        f"which the file must have: the validated ones are {RECOMMENDED_AGGREGATIONS}, and any "
        f"other pandas aggregation reducing numbers to a single number works too (see "
        f"'--reaction-gene-aggregation'). An aggregation pools the values of each map element. "
        f"These are the values that the maps of the individual samples draw. 'std' is an example "
        f"here, mapping how much replicate samples disagree. This summary colors the 'unified' map "
        f"when the samples are not grouped, and each group's map when they are grouped with "
        f"'--groups-txt' — though only an aggregation, or a rule for categorical values, affects a "
        f"group's map, since presence there is always the count of the group's own samples, so the "
        f"presence names are rejected in a grouped run as having nothing to choose between. Maps "
        f"for individual samples are never summarized, always showing that one sample alone: its "
        f"own values with a value column, its presence without one. Without this option, samples "
        f"are summarized by presence, by membership with 3 or fewer samples, by count above that, "
        f"and by a continuous count scale where a discrete one would run out of distinguishable "
        f"colors or of room to label its bands (above {MAX_DISCRETE_COUNT_BANDS} of them). "
        f"Presence is the default because it is meaningful for any set of samples, whereas pooling "
        f"values is only meaningful when the samples are commensurable, such as replicates of one "
        f"condition. If no summary of the samples is meaningful or desired, use "
        f"'--skip-unified-maps' to leave the map out. Without groups, '--reaction-sample-summary' "
        f"colors the 'unified' map alone; an ungrouped run with '--skip-unified-maps' therefore "
        f"refuses this option. With '--reaction-value-period', its aggregations are "
        f"'circular_mean' and, without groups, 'angular_deviation'. With '--reaction-value-type "
        f"categorical', aggregations are refused. Presence stays the default. It counts the "
        f"samples, however few there are. Coloring by which samples contain an element would use "
        f"the 'tab10' colors of the values. A color would then mean a sample on one map and a "
        f"value on another. Without groups, a rule ({CATEGORICAL_RULES_PHRASE}) instead chooses "
        f"one value across the samples for the 'unified' map (see '--reaction-value-type'). With "
        f"'--groups-txt', the rule chooses one value among the samples of each group for that "
        f"group's map. The 'unified' map then counts the groups, unless '--reaction-group-summary' "
        f"gives a rule. Use a rule only where a value means the same thing in every sample."
    )
    groupSUMMARY.add_argument(
        '--compound-sample-summary', metavar='NAME', help=
        "Like '--reaction-sample-summary', but for the compound layer. The two layers can be "
        "summarized differently."
    )
    groupSUMMARY.add_argument(
        '--reaction-group-summary', metavar='NAME', help=
        "How to summarize THE SAMPLE GROUPS of '--groups-txt' for the reaction layer, which colors "
        "the 'unified' map of a grouped run. 'count', 'count_continuous' and 'membership' "
        "summarize presence, coloring an element by how many or exactly which groups contain it, a "
        "group containing it subject to the '--group-threshold'; 'count' and 'count_continuous' "
        "differ only in whether the scale is drawn in discrete bands or as a gradient, as "
        "described for '--reaction-sample-summary'. Any other name is an aggregation pooling the "
        "groups' values, each group's value being its own samples summarized by "
        "'--reaction-sample-summary' (which must therefore also pool values). Without this option, "
        "groups are summarized by presence, by membership with 3 or fewer groups, by count above "
        "that, and by a continuous count scale where a discrete one would run out of "
        "distinguishable colors or of room to label its bands. A useful combination for replicates "
        "is '--reaction-sample-summary mean --reaction-group-summary count': each group's map "
        "shows the mean of its replicate samples, while the 'unified' map shows the number of "
        "groups in which each element occurs. If no summary of the groups is meaningful or "
        "desired, use '--skip-unified-maps' to leave the map out, and without a map to color, "
        "'--reaction-group-summary' is refused. With '--reaction-value-period', its aggregations "
        "are 'circular_mean' and 'angular_deviation'. With '--reaction-value-type categorical', "
        "presence stays the default. It counts the groups, however few there are. A rule chooses "
        "one value across the groups instead (see '--reaction-value-type'). It needs "
        "'--reaction-sample-summary' to give a rule as well, to choose each group's value among "
        "its samples. Each group with a value gives one vote, whatever its size. With few groups, "
        "ties are common. With two groups, every disagreement is a tie, so 'plurality_order' "
        "settles it by the order of the values."
    )
    groupSUMMARY.add_argument(
        '--compound-group-summary', metavar='NAME', help=
        "Like '--reaction-group-summary', but for the compound layer. The two layers can be "
        "summarized differently."
    )

    groupOUT = parser.add_argument_group("OUTPUT", "Output files.")
    groupOUT.add_argument(*A('output-dir'), **K('output-dir'))
    groupOUT.add_argument(*A('overwrite-output-destinations'), **K('overwrite-output-destinations'))
    groupOUT.add_argument(
        '--name-files', action='store_true', default=False, help=
        "Include the pathway name along with the number in output map file names. For example, in "
        "drawing KO presence/absence data, by default, the 'Glycolysis / Gluconeogenesis' map "
        "would be saved to a file named 'ko00010.pdf', but, with this flag, the map would be saved "
        "to a file named 'ko00010_Glycolysis_Gluconeogenesis.pdf'. To further illustrate how "
        "special characters in pathway names are translated to characters in file paths, the file "
        "name for 'Glycosylphosphatidylinositol (GPI)-anchor biosynthesis' would be "
        "'ko00563_Glycosylphosphatidylinositol_(GPI)_anchor_biosynthesis.pdf', and the file name "
        "for 'Biosynthesis of 12-, 14- and 16-membered macrolides' would be "
        "'ko00522_Biosynthesis_of_12_14_and_16_membered_macrolides.pdf' with this flag."
    )
    groupOUT.add_argument(
        '--categorize-files', action='store_true', default=False, help=
        "Categorize output files by pathway map within subdirectories corresponding to the BRITE "
        "hierarchy of maps (see https://www.genome.jp/brite/br08901). Links to these files are "
        "also gathered in one directory for easier browsing. For example, if drawing a map file "
        "for 'Glycolysis / Gluconeogenesis', '00010', then the file will be written to "
        "subdirectory 'Metabolism/Carbohydrate_metabolism' of the output directory, and a link "
        "with the same basename as the file will be created in subdirectory 'all_maps' of the "
        "output directory. If drawing a map grid for 'Quorum sensing', '02024', then the file will "
        "be written to a subdirectory named "
        "'grid/Cellular_Processes/Cellular_community_prokaryotes', and a hard link will be created "
        "in a subdirectory named 'grid/all_maps'."
    )
    groupOUT.add_argument(
        '--draw-individual-files', nargs='*', help=
        "Draw pathway maps for individual contigs databases if multiple ungrouped databases are "
        "provided, for individual genomes of an ungrouped pangenome, or for the individual samples "
        "of a draw-kegg-pathways text file with a 'sample' column. If used as a flag (without "
        "values), save files for all of the individual databases, genomes, or samples. "
        "Alternatively, the project names of a subset of contigs databases, or the names of a "
        "subset of genomes or samples, can be provided. If groups are defined by '--groups-txt', "
        "then maps are drawn for individual database, pan genome, or sample groups instead. A "
        "subset of group names can be provided to draw maps for select groups."
    )
    groupOUT.add_argument(
        '--draw-grid', nargs='*', help=
        "Draw a grid for each pathway map. If using multiple ungrouped contigs databases, the grid "
        "shows the 'unified' map of data from all databases and maps for individual databases. If "
        "using an ungrouped pangenomic database, the grid shows the pangenomic map and maps for "
        "individual genomes. If using a draw-kegg-pathways text file with a 'sample' column, the "
        "grid shows the 'unified' map summarizing all of the samples and maps for individual "
        "samples. The grid view facilitates identification of the contigs databases, genomes, or "
        "samples containing reactions highlighted in the integrative map. If used as a flag "
        "(without values), all of the contigs databases, genomes, or samples are included in the "
        "grid. Alternatively, the project names of a subset of contigs databases, or the names of "
        "a subset of genomes or samples, can be provided. If groups are defined by '--groups-txt', "
        "then the grid instead shows the 'unified' map alongside maps for individual contigs "
        "database, pan genome, or sample groups. A subset of group names can be provided to select "
        "maps in the grid. With '--skip-unified-maps' each grid only shows the individual maps, "
        "leaving out the first 'unified' map panel."
    )
    groupOUT.add_argument(
        '--skip-unified-maps', action='store_true', default=False, help=
        "Skip the 'unified' map. That map summarizes every sample, group, contigs database, or "
        "genome at once. Skipping it drops the files written under 'unified' and the first panel "
        "of each map grid. The summary is not always wanted. Pooling values is only meaningful "
        "across commensurable samples. The default 'unified' map display of the presence of "
        "sample or group sources says little when the point is to compare the samples or groups. A "
        "run normally draws the map even when no summary was asked for. This option does not "
        "affect the maps of individual samples or groups and the color scale they share. The "
        "colorbar for the summary ('colorbar_reactions.pdf' or 'colorbar_compounds.pdf') is kept "
        "only where the same scale also colors the individual maps. '--draw-individual-files' "
        "and/or '--draw-grid' are therefore required with '--skip-unified-maps'. Options that only "
        "shape the summary are refused: the sample and group summaries, "
        "'--presence-colormap-scheme', and the limits and center of the summary's own value scale. "
        "'--group-threshold' is also refused, though a grouped run otherwise requires it."
    )
    groupOUT.add_argument(
        '--collate-files-by-map', action='store_true', default=False, help=
        "Gather the map files drawn by '--draw-individual-files' into a directory per map, each "
        "holding one file per sample, group, contigs database, or pan genome named after it. This "
        "is a second arrangement of the same files, standing beside the directory per source that "
        "'--draw-individual-files' writes rather than replacing it, and it lets a single map be "
        "compared across sources by stepping through one directory, as a file browser's preview "
        "moves from one file to the next. For example, with samples 'A' through 'E', the "
        "'Glycolysis / Gluconeogenesis' map would be gathered into a directory named "
        "'by_map/ko00010' holding files named 'A.pdf' through 'E.pdf', or "
        "'by_map/ko00010_Glycolysis_Gluconeogenesis' with '--name-files'. With "
        "'--categorize-files', the BRITE subdirectories precede each map's directory, just as "
        "they precede each map file elsewhere in the output. The gathered files are links to the "
        "maps already drawn, so they take up no disk space of their own."
    )
    groupOUT.add_argument(
        '--draw-bare-maps', action='store_true', default=False, help=
        "By default, without this flag, only draw maps containing select input data. Even if "
        "pathway maps are given explicitly with '--pathway-numbers' (e.g., \"00010 01100\"), if "
        "they do not contain input data, they are not drawn unless this flag is used."
    )

    groupMAP = parser.add_argument_group("MAP", "Pathway maps to draw.")
    groupMAP.add_argument(
        '--pathway-numbers', type=str, nargs='+', help=
        "Five-digit numbers identify pathway maps to draw. By default, all maps are drawn. Numbers "
        "are five-digits long. This argument accepts regular expression patterns. For example, the "
        "values, 01100 03... , will draw the global 'Metabolic pathways' map '01100' and all of "
        "the 'Genetic Information Processing' maps with numbers starting '03'. See the following "
        "website for a classification of the maps: https://www.genome.jp/kegg/pathway.html"
    )
    groupMAP.add_argument(
        '--kegg-dir', type=str, metavar='PATH', help=
        "Path to KEGG database directory containing map files. If this option is not used, the "
        "program expects a database set up in the default location used by 'anvi-setup-kegg-data'."
    )

    groupCOLOR = parser.add_argument_group("COLOR", "Everything to do with coloring.")
    groupCOLOR.add_argument(
        '--reaction-color', nargs='?', const=DEFAULT_REACTION_COLOR, metavar='HEX', help=
        "Color reaction elements a single static color for presence/absence, rather than "
        "dynamically (across databases, genomes, samples, or groups) or by a value column. Used as "
        "a flag (without a value), reactions are colored the default color (green). Alternatively, "
        "give a color hex code, such as '#FFA500' for orange. A COLOR HEX CODE MUST BE ENCLOSED IN "
        "QUOTES, SINCE AN UNQUOTED '#' MAKES THE SHELL IGNORE THE REST OF THE COMMAND. This "
        "applies to the reaction layer of any input; for multiple databases, genomes, samples, or "
        "groups it overrides dynamic coloring to show simple presence/absence. It cannot be "
        "combined with '--reaction-colormap' (dynamic reaction coloring) or '--original-color'."
    )
    groupCOLOR.add_argument(
        '--compound-color', nargs='?', const=DEFAULT_COMPOUND_COLOR, metavar='HEX', help=
        "Like '--reaction-color', but for the compound layer of a '--compound-txt' file colored by "
        "presence/absence. Used as a flag, compounds are colored the default color (pink); or give "
        "a quoted color hex code. Database, pangenome, and reaction-network-JSON inputs have no "
        "compound layer, so this applies only to a compound text file. It cannot be combined with "
        "'--compound-colormap' or '--original-color'."
    )
    groupCOLOR.add_argument(
        '--original-color', action='store_true', default=False, help=
        "Color reaction elements using the original colors of the reference map, rather than a "
        "single color, dynamic coloring, or a value column. On global and overview maps, the "
        "compounds of those reactions are automatically colored to match, just as they appear in "
        "the reference maps. This draws only from the reaction layer, through a dedicated "
        "renderer, so it cannot also take an explicit compound file ('--compound-txt'), and it "
        "cannot be combined with '--reaction-color', '--compound-color', or the colormaps. It "
        "works across multiple sources: with several contigs databases, the genomes of a "
        "pangenome, or a '--reaction-txt' file with a 'sample' column, the 'unified' map "
        "highlights the union of them all and '--draw-individual-files'/'--draw-grid' add a "
        "reference-colored map per database, genome, or sample. Because the reference map fixes "
        "both the colors and their drawing order, there is no color scale left to give away: a "
        "value column, the sample and group summary options, '--reaction-reverse-overlay' and "
        "'--group-threshold' therefore cannot be combined with it. With '--groups-txt', the "
        "'unified' map is still the union, while each group's own map counts the group's sources "
        "instead, since one color cannot distinguish them."
    )
    groupCOLOR.add_argument(
        '--reaction-colormap', nargs='+', help=
        "This option takes the name of a colormap which is sampled in coloring data. This can be a "
        "Matplotlib colormap or 'clocktime', a colormap defined in anvi'o to show clock times. In "
        "addition to the colormap name, two decimal values between 0.0 and 1.0, with the first "
        "value smaller than the second, can be provided to limit the fraction of the colormap "
        "used. For example, the values, plasma 0.2 0.9 , would extract 70%% of the 'plasma' "
        "colormap, ignoring the darkest 20%% and lightest 10%%. Here is how a colormap is applied "
        "to KO occurrence data. KO reactions can be dynamically colored by occurrence in multiple "
        "contigs databases, the genomes of a pangenome, the samples of a draw-kegg-pathways text "
        "file, or groups of any of these. Pangenomes by default use the sequential colormap, "
        "'plasma_r' ('_r' can be added to colormap names to reverse the order of colors), trimming "
        "the top and bottom 10%%. 'plasma_r' spans yellow (fewer genomes) to blue-violet (more "
        "genomes), which accentuates in darker colors reactions that are shared rather than "
        "unshared across genomes. In contrast, a colormap spanning dark to light, such as "
        "'plasma', can be better for drawing attention to unshared reactions. Multiple contigs "
        "databases and groups can use three 'schemes' for dynamic coloring: 'by_count', "
        "'by_count_continuous', or 'by_membership' (see '--presence-colormap-scheme'). As with "
        "pangenomes, the two count schemes use by default the 'plasma_r' colormap, trimming the "
        "top and bottom 10%%. 'by_membership' uses by default the qualitative colormap, 'tab10', "
        "without trimming. This colormap contains distinct colors suitable for clearly "
        "differentiating the databases or groups containing reactions. All of the above concerns "
        "occurrence coloring. When coloring a reaction layer by a value column instead (a "
        "'--reaction-txt' file with a value column of numbers), this same option selects the "
        "sequential colormap that is sampled continuously to map reaction values (default "
        "'plasma_r'), and '--presence-colormap-scheme' does not apply. That one colormap then "
        "colors the 'unified' map and the maps of the individual samples or groups alike, unless "
        "'--reaction-individual-colormap' gives the latter a colormap of their own. With "
        "'--reaction-value-type categorical', each value takes a color of its own. This option "
        "then colors only presence across the samples or groups on the 'unified' map. It is "
        "refused where the 'unified' map shows the values, or where no 'unified' map is drawn. See "
        "the following webpage for named Matplotlib colormaps: "
        "https://matplotlib.org/stable/users/explain/colors/colormaps.html#classes-of-colormaps"
    )
    groupCOLOR.add_argument(
        '--reaction-individual-colormap', nargs='+', help=
        "Like '--reaction-colormap', but for the single color scale shared by the individual maps "
        "('--draw-individual-files'/'--draw-grid'). These are the maps of the samples, or of the "
        "groups when the samples are grouped with '--groups-txt'. It takes the same values: a "
        "colormap name on its own, or a name followed by two decimals limiting the fraction of the "
        "colormap to sample. Without this option, one colormap colors both of those maps and the "
        "'unified' map, which is what you want when the two show the same quantity — the values of "
        "one sample beside the mean of them all, for example. Give it when they do not. With "
        "'--reaction-sample-summary std', for instance, the 'unified' map shows how much replicate "
        "samples disagree, a spread that is never negative and will not conform to the scale of "
        "the underlying values, so drawing it in the colors that show magnitude on the per-sample "
        "maps invites mistaking one scale for the other. A diverging colormap could make sense for "
        "per-sample values running either side of zero, and a sequential colormap given here for "
        "'unified' map values like standard deviation that only grow. 'Unified' and "
        "per-sample/group maps have their own colorbars regardless ('colorbar_reactions.pdf' for "
        "the 'unified' map and 'colorbar_reactions_samples.pdf'/'colorbar_reactions_groups.pdf'), "
        "so nothing about which colors mean what is left implicit. This option only applies where "
        "those maps are colored by value: the reaction file needs a value column and a 'sample' "
        "column, and, with '--groups-txt', '--reaction-sample-summary' must pool the values of "
        "each group's samples rather than summarize their presence, which '--group-colormap' "
        "colors instead. It is refused for categorical values ('--reaction-value-type "
        "categorical'), which take colors of their own."
    )
    groupCOLOR.add_argument(
        '--compound-individual-colormap', nargs='+', help=
        "Like '--reaction-individual-colormap', but for the compound layer of a '--compound-txt' "
        "file colored by a value column. The two layers are colored on independent scales, each "
        "with its own colorbar, so each takes its own colormaps."
    )
    groupCOLOR.add_argument(
        '--reaction-membership-colors', type=str, metavar='PATH', help=
        "Path to a tab-delimited kegg-membership-colors-txt file. It gives a color hex code to "
        "each sample of a '--reaction-txt' file, each contigs database, each genome of a "
        "pangenome, or each group of any of these when '--groups-txt' is used. The file colors the "
        "reaction layer by membership, in colors you choose. It replaces '--reaction-colormap', so "
        "the two cannot be combined. The first column of the file holds the name of each sample, "
        "database, genome, or group. It can be headed anything. Its 'color' column holds the hex "
        "codes. Coloring by membership colors an element by exactly WHICH samples, databases, "
        "genomes, or groups contain it, so a color is needed for every combination of them. Each "
        "of them on its own takes its own color. A combination takes the blend of the colors of "
        "its members. A row naming the combination overrides the blend. The row gives the names of "
        "the members separated by commas, exactly as the colorbar labels the combination. Without "
        "groups, the same colors also color the individual map of each sample, database, or genome "
        "('--draw-individual-files'/'--draw-grid'). The color of a panel of a grid then says which "
        "sample, database, or genome it is. That color matches the band of the same sample, "
        "database, or genome on the 'unified' map. For a layer with a value column, the map of "
        "each sample shows that sample's values instead. With groups, the individual maps of each "
        "group count the group's samples, databases, or genomes. Only '--group-colormap tint' puts "
        "the group's color on those maps. It colors them with a ramp of tints of the group's "
        "color. Where '--reaction-sample-summary' gives an aggregation or a rule, a group's maps "
        "show values, not counts. Every sample, database, genome, or group of the run needs a "
        "color. Names in the file that are not in the run are reported and ignored, so one file "
        "can serve runs that draw different subsets. Note that giving this file colors the layer "
        "by membership however many samples, databases, genomes, or groups there are. Without the "
        "file, a layer with more than three of them is colored by count. A layer of categorical "
        "values ('--reaction-value-type categorical') is colored by count even with two or three "
        "samples or groups. With the file, six of them take a 63-band scale of every combination, "
        "rather than a six-band scale of counts. No color scale can distinguish the combinations "
        "of more than a handful of them, so anvi'o refuses rather than draw a scale whose bands "
        "cannot be told apart."
    )
    groupCOLOR.add_argument(
        '--compound-membership-colors', type=str, metavar='PATH', help=
        "Like '--reaction-membership-colors', but for the compound layer of a '--compound-txt' "
        "file colored by sample or group presence. The two layers are colored independently, so "
        "each takes its own file. Pointing both options at a single file gives the two layers the "
        "same color for each sample or group, and so one legend. It cannot be combined with "
        "'--compound-colormap'."
    )
    groupCOLOR.add_argument(
        '--reaction-value-colors', type=str, metavar='PATH', help=
        "Path to a tab-delimited kegg-value-colors-txt file giving a color hex code to each "
        "categorical value of the reaction layer ('--reaction-value-type categorical'). Its first "
        "column holds the values and can be headed anything. Its 'color' column holds the hex "
        "codes. The order of the rows is the order of the values. That order sets the legend, the "
        "drawing order of overlapping elements, and the choices of the rules 'order' and "
        "'plurality_order'. Every value of the reaction file needs a color. Values that the "
        "reaction file does not have are reported and ignored, so one colors file can serve "
        "several runs. Two values cannot share a color. With '--groups-txt' and a 'sample' column "
        "in the reaction file, the colors file needs a rule in '--reaction-sample-summary'. "
        "Without one, each group's map counts its samples, and no map shows the values."
    )
    groupCOLOR.add_argument(
        '--compound-value-colors', type=str, metavar='PATH', help=
        "Like '--reaction-value-colors', but for the categorical values of the compound layer."
    )
    groupCOLOR.add_argument(
        '--presence-colormap-scheme',
        choices=['by_count', 'by_count_continuous', 'by_membership'], help=
        f"There are three ways of dynamically coloring elements by occurrence in multiple contigs "
        f"databases ('--contigs-dbs'), the genomes of a pangenome, or groups ('--groups-txt') of "
        f"either: by count in discrete bands, by count on a continuous scale, or by membership. "
        f"For a draw-kegg-pathways text file this choice is part of the sample and group summaries "
        f"instead ('--reaction-sample-summary' and the related options), so this option is "
        f"rejected for a text run. It is also rejected with '--skip-unified-maps', which leaves no "
        f"map for it to color. A group's own map is drawn by '--group-colormap-scheme'. A "
        f"database's or genome's own map draws its reactions in a single color. By default, with 4 "
        f"or more databases or groups, reactions are colored by count of database or group; with 2 "
        f"or 3, reactions are colored explicitly by database or group membership. In coloring by "
        f"count, the colormap should be sequential, such that the color of a reaction changes "
        f"'smoothly' with the count. Because its colorbar labels one band per count, 'by_count' "
        f"additionally needs a distinguishable color for every count, and is refused when there "
        f"are more databases or groups than the colormap can supply — with a few hundred of them "
        f"no color scale could tell them apart anyway — and it needs room to set every one of "
        f"those labels, so above about {MAX_DISCRETE_COUNT_BANDS} counts it warns that they will "
        f"be too small to read. 'by_count_continuous' colors counts from the same colormap in the "
        f"same way but draws its colorbar as a gradient from the lowest count to the highest, "
        f"which needs neither a distinguishable color nor a label per count, so a color there "
        f"reads as a position along the range rather than as an exact count; it is chosen "
        f"automatically, with a warning, wherever 'by_count' would be refused or its labels would "
        f"be unreadable. In contrast, coloring by membership means that the color of a reaction is "
        f"set by its membership in a database or group, or in a combination of databases or "
        f"groups, so a qualitative colormap can be used instead of a sequential colormap, to give "
        f"a distinct color to each combination. This is the default with 2 or 3 databases or "
        f"groups."
    )
    groupCOLOR.add_argument(
        '--count-scale-max', metavar='NAME_OR_NUMBER', default='observed', help=
        "Where a color scale of counts stops, which decides how the colors are spread over the "
        "counts. By default, 'observed' stops at the highest count anything on the drawn maps "
        "actually has, so that the colors span the counts that occur. 'total' runs the scale to "
        "the number of samples, contigs databases, genomes, or groups being counted. That number "
        "is the highest the count could be. A number stops the scale there instead, and counts "
        "above it take the top color, which is how separate figures are given the same scale. The "
        "default matters most for sparse data: with 100 samples and no element in more than 10 of "
        "them, a scale running to 100 draws every element in indistinguishable shades at the "
        "bottom of the colormap, while one running to 10 spreads them across it. Use 'total' or a "
        "number if a fixed scale matters more, since an 'observed' scale depends on which maps "
        "were drawn and so is not comparable between runs that draw different maps. This applies "
        "to every count scale: the 'unified' map's, per layer, and the within-group source counts "
        "of individual group maps. Coloring by membership is unaffected, needing a color per "
        "combination; so is coloring by a value column, which always spans the input values."
    )
    groupCOLOR.add_argument(
        '--reaction-reverse-overlay', action='store_true', default=False, help=
        "By default, without this flag, reactions with a greater numerical value (e.g., in more "
        "contigs databases, pan genomes, samples, or groups, or with a higher reaction value "
        "column) are drawn on top of those with a lesser value. With this flag, the opposite "
        "applies. This drawing order matters most in global and overview maps, where reaction "
        "arrows frequently cross and overlap, so which reaction is drawn on top is actually "
        "visible; especially with a non-default colormap spanning dark to light, reversing the "
        "order there accentuates unshared rather than shared parts of a pathway. (In standard "
        "maps, reactions rarely overlap, so the order seldom matters.) This applies to the "
        "reaction layer across every source: contigs databases, pan genomes, and the reaction "
        "layer of draw-kegg-pathways text files ('--reaction-txt'). The compound layer has its own "
        "'--compound-reverse-overlay'. For categorical values ('--reaction-value-type "
        "categorical'), the value that comes first in their order is drawn on top, and this flag "
        "puts the last value on top instead. It does not change which value a rule chooses."
    )
    groupCOLOR.add_argument(
        '--compound-colormap', nargs='+', help=
        "Like '--reaction-colormap', but for the compound layer when drawing compounds from a "
        "'--compound-txt' file colored by a value column. The reaction and compound layers are "
        "colored on independent scales, each with its own colorbar (both default 'plasma_r'). As "
        "with '--reaction-colormap', this takes a colormap name, optionally followed by two "
        "decimal limits. With '--compound-value-type categorical', it colors only presence across "
        "the samples or groups on the 'unified' map."
    )
    groupCOLOR.add_argument(
        '--compound-reverse-overlay', action='store_true', default=False, help=
        "Like '--reaction-reverse-overlay', but for the compound layer of a '--compound-txt' file. "
        "By default, compounds with a greater value are drawn on top of those with a lesser value; "
        "with this flag, the opposite applies. Note that compound circles rarely overlap on a map, "
        "unlike the crossing reaction arrows of global and overview maps, so the compound drawing "
        "order is usually of little visual consequence. For categorical values "
        "('--compound-value-type categorical'), it puts the last value on top, as "
        "'--reaction-reverse-overlay' does."
    )
    groupCOLOR.add_argument(
        '--reaction-value-limits', nargs=2, metavar='LIMIT', help=
        "Limits on the values that the reaction layer's color scale spans on the 'unified' map, "
        "given as a minimum and then a maximum, either of which can be the word 'none' to leave "
        "that end where the data puts it. For instance, '--reaction-value-limits -6 none' starts "
        "the scale at -6 and lets it run up to whatever the values reach. A limit sets its end of "
        "the scale whether or not any value reaches it. Every element past a limit takes the color "
        "of that end of the scale. The colorbar labels each limit. Where values lie past a limit, "
        "the label at the limit is marked '<=' or '>='. This tells a reader that the color stands "
        "for that value or anything beyond it, rather than for the value alone. The limits are "
        "read in the units of the colorbar, which are the values of map elements once both "
        "reductions have been applied ('--reaction-gene-aggregation' and then "
        "'--reaction-accession-aggregation'), so with the default of 'sum' at the accession level, "
        "a single element standing for a dozen KOs may sit well past any one value in the file. "
        "With a 'sample' column, the 'unified' map shows a summary of these element values across "
        "the samples or groups. With '--reaction-slice-by', the limits are read in the values of "
        "slices instead. Each slice shows the value of one accession or gene."
    )
    groupCOLOR.add_argument(
        '--reaction-individual-value-limits', nargs=2, metavar='LIMIT', help=
        "Like '--reaction-value-limits', but for the single color scale shared by the individual "
        "maps ('--draw-individual-files'/'--draw-grid'). These are the maps of the samples, or of "
        "the groups when the samples are grouped with '--groups-txt'. The two scales are limited "
        "separately because they can span quite different things: where '--reaction-sample-summary "
        "mean' colors the 'unified' map, that map shows means across samples, which cluster far "
        "more tightly than the per-sample values underlying them, so limits that suit one can be "
        "badly wrong for the other. Giving one of the two options without the other is perfectly "
        "legitimate, and anvi'o warns when you do, since the scale left alone goes on spanning "
        "whatever its own values happen to reach."
    )
    groupCOLOR.add_argument(
        '--compound-value-limits', nargs=2, metavar='LIMIT', help=
        "Like '--reaction-value-limits', but for the compound layer of a '--compound-txt' file "
        "colored by a value column. The two layers are colored on independent scales, each with "
        "its own colorbar, so each takes its own limits."
    )
    groupCOLOR.add_argument(
        '--compound-individual-value-limits', nargs=2, metavar='LIMIT', help=
        "Like '--reaction-individual-value-limits', but for the compound layer: the scale shared "
        "by the compound colors of the individual sample or group maps."
    )
    groupCOLOR.add_argument(
        '--reaction-value-center', nargs='?', const=VALUE_CENTER_DEFAULT, metavar='VALUE', help=
        f"Put a value at the middle of the reaction layer's color scale on the 'unified' map. "
        f"Given as a bare flag without an argument, it centers the scale on "
        f"{VALUE_CENTER_DEFAULT}, which is what a signed quantity is read against. Give a number "
        f"instead to center the scale on that. Centering is important in making use of a diverging "
        f"colormap, e.g., '--reaction-colormap RdBu_r'. Such a colormap has a neutral color in the "
        f"middle and two opposed ramps either side — only a centered scale puts the correct value "
        f"at the neutral color. Without it, a scale spans exactly the values it is given, so "
        f"values running from -4 to +6 put the neutral color at +1. The scale is widened, never "
        f"narrowed: it is adjusted to run the same distance on either side of the center, reaching "
        f"as far as the farther of its two ends did — nothing is clipped that was not clipped "
        f"before with '--*-value-limits' arguments. The price is that the shorter side of the "
        f"colormap goes partly unused. Centering never moves an end set by "
        f"'--reaction-value-limits'. A center that would move one is refused. With both limits "
        f"given, the center has to be their midpoint."
    )
    groupCOLOR.add_argument(
        '--reaction-individual-value-center', nargs='?', const=VALUE_CENTER_DEFAULT,
        metavar='VALUE', help=
        "Like '--reaction-value-center', but for the single color scale shared by the individual "
        "maps ('--draw-individual-files'/'--draw-grid'). These are the maps of the samples, or of "
        "the groups when the samples are grouped with '--groups-txt'. The two scales are centered "
        "separately, for the same reason they are limited separately: a summary can put them on "
        "quite different footings. Giving one of the two without the other is legitimate, and "
        "anvi'o warns when you do, since one colormap colors both scales unless "
        "'--reaction-individual-colormap' gives them one each, and the middle color would then "
        "mean the centered value on the one map and whatever the values leave in the middle on the "
        "other."
    )
    groupCOLOR.add_argument(
        '--compound-value-center', nargs='?', const=VALUE_CENTER_DEFAULT, metavar='VALUE', help=
        "Like '--reaction-value-center', but for the compound layer of a '--compound-txt' file "
        "colored by a value column. The two layers are colored on independent scales, each with "
        "its own colorbar, so each takes its own center."
    )
    groupCOLOR.add_argument(
        '--compound-individual-value-center', nargs='?', const=VALUE_CENTER_DEFAULT,
        metavar='VALUE', help=
        "Like '--reaction-individual-value-center', but for the compound layer: the scale shared "
        "by the compound colors of the individual sample or group maps."
    )
    groupCOLOR.add_argument(
        '--reaction-element-normalization', nargs='+', metavar='NAME', help=
        f"Color the reaction elements of sample or group maps after a transformation comparing "
        f"their levels to the other maps: each element's quantitative value is rescaled against "
        f"that same element's values across all of the samples or groups. Given a value of "
        f"'relative_to_mean', a reaction element value of -0.26 means that the is 26%% less "
        f"abundant than the average value across samples. Certain normalizations including "
        f"'relative_to_mean' are defined: {ELEMENT_NORMALIZATION_PHRASE}. A pandas Series method "
        f"that transforms each value to a new value, such as 'abs', can be used instead. A second "
        f"value after the normalization can be given to label the colorbar of the rescaled scale: "
        f"this replaces the default label composed from the normalization and the header of the "
        f"'--reaction-txt' file value column. Remember to quote the argument value if it has "
        f"spaces, e.g., 'relative_to_mean \"change vs. mean\"'. '--draw-individual-files' and/or "
        f"'--draw-grid' are required with this argument; only the maps of individual samples or "
        f"groups are rescaled, while the 'unified' map summarizes unnormalized values. With "
        f"'--groups-txt', rescaling happens across groups, so '--reaction-sample-summary' has to "
        f"pool each group's samples into a value: left at its default, it counts by occurrence "
        f"among the group's samples, which isn't something that should be rescaled. A "
        f"normalization whose neutral value is zero, such as 'relative_to_mean', centers the "
        f"rescaled scale on zero and colors it with the diverging colormap, "
        f"'{DEFAULT_CENTERED_COLORMAP}'; override with '--reaction-individual-value-center' and "
        f"'--reaction-individual-colormap'. An element the normalization cannot value is left "
        f"uncolored — a ratio needs a positive value to measure against, and a z-score needs an "
        f"element found in more than one sample. With '--reaction-value-period', the only "
        f"normalization is 'difference_from_circular_mean'. Each value becomes its offset from the "
        f"element's circular mean across the samples or groups, the shorter way around the circle. "
        f"With a period of 24, the offsets run from -12 to 12, and -12 and 12 are the same point. "
        f"Their scale is fixed, with 0 in the middle. It takes the cyclic 'clocktime' unless "
        f"'--reaction-individual-colormap' names another colormap. An element found in one sample "
        f"has an offset of 0. An element whose values cancel out has no circular mean, so it is "
        f"left uncolored, and a warning says so."
    )
    groupCOLOR.add_argument(
        '--compound-element-normalization', nargs='+', metavar='NAME', help=
        "Like '--reaction-element-normalization', but for the compound layer of a '--compound-txt' "
        "file colored by a value column. The two layers are rescaled independently, each on its "
        "own scale and with its own colorbar, so each takes its own normalization and its own "
        "label."
    )
    groupGROUP.add_argument(
        '--group-colormap', nargs='+', help=
        f"This option is like '--reaction-colormap', but only applies to drawing files for "
        f"individual groups ('--draw-individual-files') and panels for individual groups in map "
        f"grids ('--draw-grid'). These maps for individual groups show data from group sources — "
        f"samples, contigs databases, or pan genomes. Presence may only be colored by count — "
        f"e.g., the number of samples in the group containing the data — not membership. Like "
        f"'--reaction-colormap', this parameter takes the name of a colormap, and optionally, two "
        f"decimal values between 0.0 and 1.0 to limit the fraction of the colormap used. The "
        f"default configuration is the same, with the colormap being 'plasma_r' and the limits "
        f"being 0.1 and 0.9. Note that per-group maps show counts in discrete bands by default, "
        f"one per count in the group's scale, so a group whose scale spans more counts than this "
        f"colormap has distinguishable colors or over more than {MAX_DISCRETE_COUNT_BANDS} of them "
        f"has its count drawn on a continuous scale instead, along with a warning saying so "
        f"('--group-colormap-scheme'). This option colors group maps that show counts. Where a "
        f"layer's samples are pooled into a value instead — "
        f"'--reaction-sample-summary'/'--compound-sample-summary' given an aggregation — its group "
        f"maps show that value, and are colored by "
        f"'--reaction-individual-colormap'/'--compound-individual-colormap', or by the layer's own "
        f"'--reaction-colormap'/'--compound-colormap' if that is not given. In place of a colormap "
        f"name, give the value '{GROUP_COLORMAP_TINT}'. It colors the maps of each group with a "
        f"ramp of tints of the group's color. The ramp runs from a pale tint to the group's color. "
        f"'--reaction-membership-colors'/'--compound-membership-colors' gives the color of each "
        f"group. Every group's ramp is built the same way, so the panels of a grid stay "
        f"comparable, while each panel shows which group it is. The two decimal values then say "
        f"how far from white the ramp starts and stops rather than which fraction of a colormap to "
        f"use, defaulting to {DEFAULT_GROUP_TINT_SPAN[0]} and {DEFAULT_GROUP_TINT_SPAN[1]}: the "
        f"pale end stops short of white because standard and overview maps keep white for their "
        f"own unhighlighted elements. A named colormap remains the default, since one built for "
        f"showing magnitude usually reads better as a count than a color chosen to identify a "
        f"group does."
    )
    groupGROUP.add_argument(
        '--group-colormap-scheme', choices=list(GROUP_SCHEME_OPTIONS), help=
        f"How the count scale of individual group maps is drawn: 'by_count' gives each count a "
        f"band of a discrete colorbar, one per source in the group, while 'by_count_continuous' "
        f"draws a gradient from the lowest count to the highest. Which color a count takes is the "
        f"same either way — it is always sampled at an even fraction of '--group-colormap' — so "
        f"this option changes the colorbar and nothing on the maps themselves. Discrete bands need "
        f"one distinguishable color per count and room to label every band, both of which a group "
        f"of many sources can exhaust (the labels above about {MAX_DISCRETE_COUNT_BANDS} counts); "
        f"a gradient needs neither, since a color on it reads as a position along the range of "
        f"counts rather than as an exact count, and so works for a group of any size. By default "
        f"the bands are drawn while they can be told apart and read, and the gradient takes over "
        f"with a warning where they cannot. Unlike '--presence-colormap-scheme', which each layer "
        f"sets for itself, one scheme covers every layer's group maps, because the reaction and "
        f"compound layers of a group's map share one set of colors and one colorbar; for the same "
        f"reason the choice is made once for the whole run rather than per group, so that the "
        f"panels of a grid all carry the same kind of colorbar. There is no 'by_membership' here: "
        f"a group's own map counts that group's sources and never shows which of them contain an "
        f"element. Like '--group-colormap', this applies only to individual group maps that show "
        f"source counts."
    )
    groupGROUP.add_argument(
        '--group-reverse-overlay', action='store_true', default=False, help=
        "This flag is like '--reaction-reverse-overlay', but only applies to drawing files for "
        "individual groups ('--draw-individual-files') and panels for individual groups in map "
        "grids ('--draw-grid'). With this flag, these maps for individual groups draw reactions "
        "with a lesser numerical value (e.g., in fewer group sources like contigs databases, pan "
        "genomes, or samples) on top of reactions with a greater numerical value, the opposite of "
        "the default drawing order. Like '--group-colormap', this applies only to individual group "
        "maps that show source counts; when a text layer's samples are summarized by pooling "
        "values, its individual group maps reuse "
        "'--reaction-reverse-overlay'/'--compound-reverse-overlay' instead."
    )

    args = parser.get_args(parser)
    return args

def check_kegg_data(args: Namespace) -> None:
    kegg_args = Namespace()
    kegg_args.kegg_data_dir = args.kegg_dir
    kegg_context = KeggContext(kegg_args)

    if not os.path.exists(kegg_context.kegg_map_image_kgml_file):
        raise ConfigError(
            "One of the key files required by KEGG pathway maps is missing in your active anvi'o "
            "installation. If your KEGG data are not stored at the default KEGG data location, "
            "include that path using the `--kegg-dir` parameter. Otherwise, please consider using "
            "the program `anvi-setup-kegg-data` to set up the latest KEGG data that includes the "
            "necessary files for KEGG pathway maps."
        )

def consolidate_contigs_dbs(args: Namespace) -> None:
    """Transfer contigs database paths from an external_genomes file to the contigs_dbs argument."""
    if args.external_genomes is None:
        return

    if args.contigs_dbs is None:
        args.contigs_dbs = []

    filesnpaths.is_file_tab_delimited(args.external_genomes, expected_number_of_fields=2)
    external_genomes_table = pd.read_csv(args.external_genomes, sep='\t', header=0)
    assert external_genomes_table.columns.tolist() == ['name', 'contigs_db_path']
    args.contigs_dbs += external_genomes_table['contigs_db_path'].tolist()

def map_json_network_ko_data(args: Namespace, mapper: Mapper) -> None:
    """Draw KO data from a reaction network JSON file."""
    # A reaction network JSON file is a single, sample-less source drawn in one color, like a
    # draw-kegg-pathways text file without a 'sample' column. Flags for comparing or dynamically
    # coloring multiple samples/groups do not apply. Rather than silently ignoring them, tell the
    # user. '--draw-individual-files'/'--draw-grid' are 'nargs' args defaulting to None, and
    # '--group-threshold' can legitimately be 0, so presence is tested with 'is not None'.
    unsupported_args: list = []
    if args.groups_txt is not None:
        unsupported_args.append('--groups-txt')
    if args.group_threshold is not None:
        unsupported_args.append('--group-threshold')
    if args.draw_individual_files is not None:
        unsupported_args.append('--draw-individual-files')
    if args.draw_grid is not None:
        unsupported_args.append('--draw-grid')
    if args.collate_files_by_map:
        unsupported_args.append('--collate-files-by-map')
    if args.reaction_colormap is not None:
        unsupported_args.append('--reaction-colormap')
    if args.compound_colormap is not None:
        unsupported_args.append('--compound-colormap')
    if args.reaction_individual_colormap is not None:
        unsupported_args.append('--reaction-individual-colormap')
    if args.compound_individual_colormap is not None:
        unsupported_args.append('--compound-individual-colormap')
    if args.reaction_membership_colors is not None:
        unsupported_args.append('--reaction-membership-colors')
    if args.compound_membership_colors is not None:
        unsupported_args.append('--compound-membership-colors')
    if args.compound_reverse_overlay:
        unsupported_args.append('--compound-reverse-overlay')
    for summary_arg, summary_flag in (
        (args.reaction_sample_summary, '--reaction-sample-summary'),
        (args.compound_sample_summary, '--compound-sample-summary'),
        (args.reaction_group_summary, '--reaction-group-summary'),
        (args.compound_group_summary, '--compound-group-summary')
    ):
        if summary_arg is not None:
            unsupported_args.append(summary_flag)
    for aggregation_arg, aggregation_flag in (
        (args.reaction_gene_aggregation, '--reaction-gene-aggregation'),
        (args.reaction_accession_aggregation, '--reaction-accession-aggregation'),
        (args.compound_accession_aggregation, '--compound-accession-aggregation')
    ):
        if aggregation_arg is not None:
            unsupported_args.append(aggregation_flag)
    if args.presence_colormap_scheme is not None:
        unsupported_args.append('--presence-colormap-scheme')
    if args.count_scale_max != 'observed':
        unsupported_args.append('--count-scale-max')
    if args.reaction_reverse_overlay:
        unsupported_args.append('--reaction-reverse-overlay')
    if args.group_colormap is not None:
        unsupported_args.append('--group-colormap')
    if args.group_colormap_scheme is not None:
        unsupported_args.append('--group-colormap-scheme')
    if args.group_reverse_overlay:
        unsupported_args.append('--group-reverse-overlay')

    if unsupported_args:
        message = ', '.join(f"'{flag}'" for flag in unsupported_args)
        raise ConfigError(
            f"The following arguments apply to comparing or dynamically coloring KOs across "
            f"multiple samples or groups, which is not supported for a reaction network JSON file "
            f"(a single source drawn in one color): {message}. Please remove these arguments; use "
            f"'--reaction-color' or '--original-color' to choose the color."
        )

    # Refused on its own rather than as another entry in the list above. That list says these
    # options are about coloring across sources, and tells the user to pick a color instead.
    # Neither fits here. The map being skipped is the only map this run draws.
    if args.skip_unified_maps:
        raise ConfigError(
            "'--skip-unified-maps' skips the map summarizing every source at once. A reaction "
            "network JSON file is a single source drawn in one color. That map is the only one "
            "this run draws. Skipping it would leave an empty output directory. Please drop this "
            "flag."
        )

    map_reaction_network_json_kos = mapper.map_reaction_network_json_kos

    # A reaction-only source: '--original-color' uses the reference colors, '--reaction-color' a
    # single color (its bare-flag default matches the method default, so it need not be forwarded).
    reaction_color = 'original' if args.original_color else args.reaction_color
    if reaction_color is not None:
        map_reaction_network_json_kos = functools.partial(
            map_reaction_network_json_kos, reaction_color=reaction_color
        )

    map_reaction_network_json_kos(
        args.reaction_network_json,
        args.output_dir,
        pathway_numbers=args.pathway_numbers,
        draw_maps_lacking_data=args.draw_bare_maps
    )

def _detect_txt_layer(path: str, value_type: str = None, weight_column: str = None) -> dict:
    """Light structural peek at a per-layer text file for CLI flag validation.

    Returns {'mode': 'quantitative'|'categorical'|'membership'|'single'|None, 'has_sample': bool,
    'has_gene_id': bool, 'has_weights': bool}. 'mode' is None when the file is malformed in a way
    the reader ('Mapper._read_element_txt') should report (no 'accession' column, an ambiguous >1
    value column, a column name given twice, or a name starting with '__'); the caller then skips
    flag guards and lets the reader raise the definitive error. The names are read as the reader
    reads them, with quotes left as text ('Mapper._read_txt_header'). A value column is
    'categorical' where 'value_type' says so, and 'quantitative' otherwise. 'weight_column' names
    the column of slice weights ('--reaction-slice-weights'). It is left out of value detection,
    and 'has_weights' says whether the file has it.
    """
    filesnpaths.is_file_exists(path)
    columns = Mapper._read_txt_header(path)
    if not columns:
        # An empty file has no columns to detect, so defer to the reader for the definitive error.
        return {'mode': None, 'has_sample': False, 'has_gene_id': False, 'has_weights': False}
    has_sample = 'sample' in columns
    has_gene_id = 'gene_id' in columns
    has_weights = weight_column is not None and weight_column in columns
    structure = {'has_sample': has_sample, 'has_gene_id': has_gene_id, 'has_weights': has_weights}
    # Pandas names each blank column itself, such as 'Unnamed: 2', so blank names are not compared.
    named = [column for column in columns if column]
    repeated = len(set(named)) < len(named)
    reserved = any(column.startswith('__') for column in columns)
    if 'accession' not in columns or repeated or reserved:
        return {'mode': None, **structure}
    value_columns = [
        c for c in columns if c not in ('accession', 'gene_id', 'sample', weight_column)
    ]
    if len(value_columns) > 1:
        return {'mode': None, **structure}
    if len(value_columns) == 1:
        mode = 'categorical' if value_type == 'categorical' else 'quantitative'
    elif has_sample:
        mode = 'membership'
    else:
        mode = 'single'
    return {'mode': mode, **structure}


def map_txt_data(args: Namespace, mapper: Mapper) -> None:
    """Draw reaction and/or compound data from per-layer draw-kegg-pathways text files."""
    def _resolve_colormap(colormap_arg, flag):
        # '--reaction-colormap'/'--compound-colormap'/'--group-colormap' are 'nargs' args: a name,
        # or a name plus two fraction limits. None means use the layer default.
        if colormap_arg is None:
            return None, None
        if len(colormap_arg) == 1:
            return colormap_arg[0], None
        # Given '--group-colormap tint' the two numbers are not a fraction of a colormap at all, but
        # how far from white a group's own ramp starts and stops, so the messages say what the
        # numbers actually mean. The constraint on them is the same either way.
        tinted = (
            flag == '--group-colormap' and colormap_arg[0] == GROUP_COLORMAP_TINT
        )
        limits_mean = (
            "how far from white towards a group's own color its ramp starts and stops"
            if tinted else "the fraction of the colormap to use"
        )
        if len(colormap_arg) != 3:
            example = (
                f"'{GROUP_COLORMAP_TINT} {DEFAULT_GROUP_TINT_SPAN[0]} {DEFAULT_GROUP_TINT_SPAN[1]}'"
            ) if tinted else "'plasma 0.1 0.9'"
            raise ConfigError(
                f"'{flag}' takes either a colormap name on its own, or a name followed by TWO "
                f"decimal values giving {limits_mean}, such as {example}. "
                f"{len(colormap_arg)} values were given: {' '.join(colormap_arg)}."
            )
        try:
            limits = (float(colormap_arg[1]), float(colormap_arg[2]))
        except ValueError:
            raise ConfigError(
                f"The two limits given to '{flag}' must be decimal numbers, but these are not: "
                f"{colormap_arg[1]}, {colormap_arg[2]}."
            )
        if not 0.0 <= limits[0] <= limits[1] <= 1.0:
            raise ConfigError(
                f"The two limits given to '{flag}' say {limits_mean}, so they must lie between 0.0 "
                f"and 1.0 with the smaller one first. These do not: {limits[0]}, {limits[1]}."
            )
        return colormap_arg[0], limits

    def _resolve_element_normalization(normalization_arg, flag):
        # '--*-element-normalization' is an 'nargs' arg: a normalization on its own, or one followed
        # by the label its colorbar takes. The normalization itself goes to the mapper as it was
        # typed, so that it is checked by the one place that knows the names,
        # 'Mapper._resolve_element_normalization', and the two cannot drift apart.
        if normalization_arg is None:
            return None, None
        if len(normalization_arg) == 1:
            return normalization_arg[0], None
        if len(normalization_arg) == 2:
            return normalization_arg[0], normalization_arg[1]
        raise ConfigError(
            f"'{flag}' takes either a normalization on its own, or a normalization followed by the "
            f"label its colorbar takes, such as 'relative_to_mean \"% change vs. "
            f"mean\"'. {len(normalization_arg)} values were given: {' '.join(normalization_arg)}. "
            f"If the label contains spaces, it has to be quoted so that it arrives as one value."
        )

    def _resolve_value_limits(limits_arg):
        # '--*-value-limits' are two-value args: a minimum and then a maximum, either of which can
        # be 'VALUE_LIMIT_NONE' to leave that end of the scale where the data puts it. Only that
        # word is resolved here; anything else goes to the mapper as it was typed, so that the
        # numbers are converted and checked in the one place ('Mapper._resolve_value_limits') that
        # programmatic callers reach too, and the two cannot drift apart.
        if limits_arg is None:
            return None
        return tuple(
            None if limit.strip().lower() == VALUE_LIMIT_NONE else limit for limit in limits_arg
        )

    def _nargs_selection(value):
        # Translate a 'nargs' individual-file/grid argument (None, [], or a list of names).
        if value is None:
            return False
        if len(value) == 0:
            return True
        return value

    # Detect the structure of each provided file for flag validation (full validation is deferred to
    # the reader). A layer is 'quantitative' (a value column of numbers), 'categorical' (a value
    # column of categorical values), 'membership' (presence + a 'sample' column), or 'single'
    # (presence, no samples). This says what the file provides. The mapper uses it to choose how the
    # layer colors each map context. A column of slice weights is not a value column.
    reaction = _detect_txt_layer(
        args.reaction_txt, args.reaction_value_type, args.reaction_slice_weights
    ) if args.reaction_txt is not None else None
    compound = _detect_txt_layer(
        args.compound_txt, args.compound_value_type
    ) if args.compound_txt is not None else None
    present = [layer for layer in (reaction, compound) if layer is not None]
    any_sample = any(layer['has_sample'] for layer in present)

    draw_individual_files = _nargs_selection(args.draw_individual_files)
    draw_grid = _nargs_selection(args.draw_grid)
    draw_individual_maps = draw_individual_files is not False or draw_grid is not False

    # Each layer's file flag, detected structure, and the two summaries reducing its samples and its
    # sample groups.
    summaries = (
        (
            'reaction', '--reaction-txt', reaction, args.reaction_sample_summary,
            args.reaction_group_summary
        ),
        (
            'compound', '--compound-txt', compound, args.compound_sample_summary,
            args.compound_group_summary
        )
    )
    # A sampled layer's groups are summarized by presence, needing '--group-threshold' to define
    # when a group contains an element, unless the group summary pools values instead. Its per-group
    # maps show within-group sample counts, which '--group-colormap' styles, unless the sample
    # summary pools values instead. A layer whose colors come from the reference map
    # ('--original-color') or from one static color ('--reaction-color'/'--compound-color')
    # summarizes nothing: its 'unified' map is presence in any sample, so no group-level presence
    # threshold applies to it.
    def _statically_colored(element_type: str) -> bool:
        return args.original_color or getattr(args, f'{element_type}_color') is not None
    group_presence = args.groups_txt is not None and any(
        layer is not None and layer['has_sample'] and not _statically_colored(element_type)
        and (group_summary is None or group_summary in PRESENCE_SUMMARIES)
        for element_type, _, layer, _, group_summary in summaries
    )
    group_presence_maps = args.groups_txt is not None and any(
        layer is not None and layer['has_sample']
        and (sample_summary is None or sample_summary in PRESENCE_SUMMARIES)
        for _, _, layer, sample_summary, _ in summaries
    )

    # Skip the flag guards if a provided file is malformed/ambiguous, so the reader raises the
    # definitive file error rather than a guard reporting a misleading secondary problem.
    if not any(layer['mode'] is None for layer in present):
        # Slicing draws reaction elements in slices, so every slicing option needs the reaction
        # file. The options that size slices also need '--reaction-slice-by'. A column of weights is
        # left out when the value column is detected. Its checks come first, so that a mistake in
        # naming it is not reported as a missing value column.
        for setting, flag, does in (
            (args.reaction_slice_by, '--reaction-slice-by',
             "draws each accession or gene of a reaction element as a slice"),
            (args.reaction_slice_weights, '--reaction-slice-weights',
             "sizes the slices of reaction elements by a column of weights"),
            (args.reaction_slice_weight_aggregation, '--reaction-slice-weight-aggregation',
             "reduces the weights of the genes of an accession to one weight"),
            (args.reaction_slice_share_summary, '--reaction-slice-share-summary',
             "pools the shares of slices across samples or groups")
        ):
            if setting is None:
                continue
            if reaction is None:
                raise ConfigError(
                    f"'{flag}' {does}, but no reaction file ('--reaction-txt') was provided."
                )
            if args.reaction_slice_by is None:
                # Weights of genes are reduced to an accession only when slicing by accession.
                units = (
                    "'--reaction-slice-by accession'"
                    if args.reaction_slice_weight_aggregation is not None
                    else "'--reaction-slice-by accession' or '--reaction-slice-by gene'"
                )
                raise ConfigError(
                    f"'{flag}' {does}. '--reaction-slice-by' was not given, so no element is "
                    f"sliced. To slice elements, give {units}. Otherwise, drop '{flag}'."
                )
        # Slices by gene need the genes of a 'gene_id' column. A file without a value column skips
        # this check, so that the refusal of slices of presence below comes first.
        if (
            args.reaction_slice_by == 'gene' and reaction['mode'] in VALUE_MODES
            and not reaction['has_gene_id']
        ):
            raise ConfigError(
                "'--reaction-slice-by gene' draws each gene of a reaction element as a slice, but "
                "the reaction file has no 'gene_id' column naming genes. Add a 'gene_id' column, "
                "or give '--reaction-slice-by accession' to slice elements by their accessions."
            )
        weight_column = args.reaction_slice_weights
        if weight_column is not None:
            if weight_column in ('accession', 'gene_id', 'sample'):
                raise ConfigError(
                    f"'--reaction-slice-weights' names the column '{weight_column}' as the weights "
                    f"of slices. That column is one of the columns that identify a row "
                    f"('accession', 'gene_id', and 'sample'), so it cannot hold weights. Name the "
                    f"column of weights instead."
                )
            if not reaction['has_weights']:
                raise ConfigError(
                    f"'--reaction-slice-weights' names the column '{weight_column}' as the weights "
                    f"of slices, but the reaction file has no column with that header."
                )
            if reaction['mode'] not in VALUE_MODES:
                raise ConfigError(
                    f"'--reaction-slice-weights' names the column '{weight_column}' as the weights "
                    f"of slices. It is the only column of the reaction file besides 'accession', "
                    f"'gene_id', and 'sample', so it would otherwise be the value column. Slices "
                    f"need a value column for their colors and a separate column for their "
                    f"weights. Add a value column, or drop '--reaction-slice-weights' to color "
                    f"equal slices by '{weight_column}'."
                )
        if args.reaction_slice_weight_aggregation is not None:
            if args.reaction_slice_by == 'gene':
                raise ConfigError(
                    "'--reaction-slice-weight-aggregation' reduces the weights of the genes of an "
                    "accession to one weight. With '--reaction-slice-by gene', each gene is a "
                    "slice of its own, so no weights are reduced to an accession. A gene with rows "
                    "for two or more accessions of an element takes the mean weight of these rows. "
                    "Please drop '--reaction-slice-weight-aggregation'."
                )
            if weight_column is None:
                raise ConfigError(
                    "'--reaction-slice-weight-aggregation' reduces the weights of the genes of an "
                    "accession to one weight, but no column of weights was named with "
                    "'--reaction-slice-weights'. Without weights, the slices of an element are "
                    "equal. Name the column, or drop '--reaction-slice-weight-aggregation'."
                )
            if not reaction['has_gene_id']:
                raise ConfigError(
                    "'--reaction-slice-weight-aggregation' reduces the weights of the genes of an "
                    "accession to one weight, but the reaction file has no 'gene_id' column naming "
                    "genes. Each accession there already has a single weight per sample, since "
                    "rows repeating an accession are refused. Please drop "
                    "'--reaction-slice-weight-aggregation'."
                )
        if args.reaction_slice_share_summary is not None and weight_column is None:
            raise ConfigError(
                "'--reaction-slice-share-summary' pools the shares of weighted slices across "
                "samples or groups, but no column of weights was named with "
                "'--reaction-slice-weights'. Without weights, the slices of an element are equal "
                "on every map, so there are no shares to pool. Name the column, or drop "
                "'--reaction-slice-share-summary'."
            )

        # A value type says what a layer's value column holds, so it needs that layer's file and its
        # value column. Colors given per categorical value need categorical values, and a map that
        # shows them. With groups, only a rule among each group's samples puts the values on a map.
        for element_type, file_flag, layer, sample_summary, _ in summaries:
            value_type = getattr(args, f'{element_type}_value_type')
            type_flag = f'--{element_type}-value-type'
            if value_type is not None and layer is None:
                raise ConfigError(
                    f"'{type_flag}' says what the value column of the {element_type} file holds, "
                    f"but no {element_type} file ('{file_flag}') was provided."
                )
            if value_type is not None and layer['mode'] not in VALUE_MODES:
                raise ConfigError(
                    f"'{type_flag}' says what the value column of the {element_type} file holds, "
                    f"but the file has no value column. Its layer is colored by presence."
                )
            colors_flag = f'--{element_type}-value-colors'
            if getattr(args, f'{element_type}_value_colors') is not None:
                if layer is None:
                    raise ConfigError(
                        f"'{colors_flag}' gives a color to each categorical value of the "
                        f"{element_type} layer, but no {element_type} file ('{file_flag}') was "
                        f"provided."
                    )
                if layer['mode'] not in VALUE_MODES:
                    raise ConfigError(
                        f"'{colors_flag}' gives a color to each categorical value of the "
                        f"{element_type} layer, but the {element_type} file has no value column. "
                        f"Its layer is colored by presence."
                    )
                if value_type != 'categorical':
                    raise ConfigError(
                        f"'{colors_flag}' gives a color to each categorical value of the "
                        f"{element_type} layer, but '{type_flag} categorical' was not given. Give "
                        f"it to read the value column as categorical values."
                    )
                if (
                    layer['has_sample'] and args.groups_txt is not None
                    and sample_summary not in CATEGORICAL_RULES
                ):
                    raise ConfigError(
                        f"'{colors_flag}' gives a color to each categorical value. Here the "
                        f"samples are grouped, and '--{element_type}-sample-summary' gives no "
                        f"rule. Each group's map therefore counts the group's samples, and no map "
                        f"shows the values. Give '--{element_type}-sample-summary' a rule, such as "
                        f"'plurality_order', to color each group's map by a value chosen among its "
                        f"samples. Or drop the colors file."
                    )

        # '--original-color' highlights the reaction layer in the reference map's colors (its
        # associated compounds follow, on global/overview maps) through a dedicated renderer that
        # draws only the reaction layer. Both the colors and their render order come from the
        # reference map, so the flag carries no data: a value column has nowhere to go, and the
        # options that shape dynamic coloring do nothing. A 'sample' column is fine, giving a
        # reference-colored map per sample and a union 'unified' map, as multiple contigs databases
        # or the genomes of a pangenome already do.
        if args.original_color:
            if args.reaction_txt is None:
                raise ConfigError(
                    "'--original-color' colors reaction elements in the reference map's colors, "
                    "but no reaction file ('--reaction-txt') was provided."
                )
            if args.compound_txt is not None:
                raise ConfigError(
                    "'--original-color' draws only the reaction layer, through a dedicated "
                    "renderer that also colors those reactions' compounds in the reference colors "
                    "(on global and overview maps), so it cannot also take an explicit compound "
                    "file ('--compound-txt')."
                )
            if reaction['mode'] in VALUE_MODES:
                raise ConfigError(
                    "'--original-color' colors reactions in the reference map's own colors, which "
                    "leaves no colors for a value column, but the reaction file has one. Drop "
                    "'--original-color' to color the reactions by their values, or drop the value "
                    "column to highlight the reactions in the reference colors."
                )
            # '--group-colormap'/'--group-reverse-overlay' are NOT listed here: with groups, the
            # per-group maps show within-group sample counts (as they do for databases and
            # pangenomes), and those two options style exactly that coloring. '--group-threshold' is
            # listed because nothing consults it here — the 'unified' map is the union of every
            # sample, and the per-group maps count a group's own samples without a threshold.
            original_flags = [flag for present_flag, flag in (
                (args.reaction_reverse_overlay, '--reaction-reverse-overlay'),
                (args.reaction_sample_summary is not None, '--reaction-sample-summary'),
                (args.reaction_group_summary is not None, '--reaction-group-summary'),
                (args.group_threshold is not None, '--group-threshold'),
                (args.reaction_slice_by is not None, '--reaction-slice-by'),
                (args.reaction_slice_weights is not None, '--reaction-slice-weights'),
                (args.reaction_slice_weight_aggregation is not None,
                 '--reaction-slice-weight-aggregation'),
                (args.reaction_slice_share_summary is not None, '--reaction-slice-share-summary')
            ) if present_flag]
            if original_flags:
                message = ', '.join(f"'{flag}'" for flag in original_flags)
                raise ConfigError(
                    f"These options were given along with '--original-color': {message}. They "
                    f"decide how elements are colored by value or across samples and sample "
                    f"groups, but '--original-color' takes both its colors and its drawing order "
                    f"from the reference map, so none of them has anything to act on. With a "
                    f"'sample' column, '--original-color' draws the union of the samples in the "
                    f"'unified' map and a map per sample, which "
                    f"'--draw-individual-files'/'--draw-grid' request."
                )

        # The mapper takes the literal string 'original' as its way of asking for the reference
        # map's colors, which is what '--original-color' sets. Passing it through a color option
        # instead would bypass every guard above, so it is rejected in favor of the flag.
        for color_value, color_flag in (
            (args.reaction_color, '--reaction-color'), (args.compound_color, '--compound-color')
        ):
            if color_value == 'original':
                raise ConfigError(
                    f"'{color_flag}' takes a color hex code, such as '#FFA500'. To use the "
                    f"reference map's own colors, use the '--original-color' flag instead of "
                    f"passing 'original' here."
                )

        # '--reaction-color'/'--compound-color' set a single presence color for a layer; they do not
        # apply to a value-colored layer, whose color comes from its colormap.
        if args.reaction_color is not None:
            if args.reaction_txt is None:
                raise ConfigError(
                    "'--reaction-color' colors the reaction layer, but no reaction file "
                    "('--reaction-txt') was provided."
                )
            if reaction['mode'] == 'quantitative':
                raise ConfigError(
                    "'--reaction-color' sets a single presence color, but the reaction layer is "
                    "colored by a value column; set its colormap with '--reaction-colormap' "
                    "instead."
                )
            if reaction['mode'] == 'categorical':
                raise ConfigError(
                    "'--reaction-color' sets a single presence color, but the reaction layer is "
                    "colored by categorical values. Set the colors of the values with "
                    "'--reaction-value-colors' instead."
                )
        if args.compound_color is not None:
            if args.compound_txt is None:
                raise ConfigError(
                    "'--compound-color' colors the compound layer, but no compound file "
                    "('--compound-txt') was provided."
                )
            if compound['mode'] == 'quantitative':
                raise ConfigError(
                    "'--compound-color' sets a single presence color, but the compound layer is "
                    "colored by a value column; set its colormap with '--compound-colormap' "
                    "instead."
                )
            if compound['mode'] == 'categorical':
                raise ConfigError(
                    "'--compound-color' sets a single presence color, but the compound layer is "
                    "colored by categorical values. Set the colors of the values with "
                    "'--compound-value-colors' instead."
                )

        # A membership colors file colors a layer by which samples or groups contain an element. It
        # also colors their individual maps, so the layer's file needs a 'sample' column, or there
        # are no samples or groups to color. Which map contexts actually take those colors depends
        # on the summaries, which the mapper resolves; it refuses there if none of them does.
        for element_type, file_flag, layer, _, _ in summaries:
            if getattr(args, f'{element_type}_membership_colors') is None:
                continue
            colors_flag = f'--{element_type}-membership-colors'
            if layer is None:
                raise ConfigError(
                    f"'{colors_flag}' colors the {element_type} layer, but no {element_type} file "
                    f"('{file_flag}') was provided."
                )
            if not layer['has_sample']:
                raise ConfigError(
                    f"'{colors_flag}' gives a color to each sample, or to each group of samples. "
                    f"The {element_type} file has no 'sample' column, so it has no samples or "
                    f"groups to color. Set the single color of this layer with "
                    f"'--{element_type}-color' instead."
                )

        # Reaction-layer coloring options apply only when a reaction layer is drawn dynamically (by
        # value or by sample/group membership), not for a single-color presence reaction layer.
        reaction_flags = [flag for present_flag, flag in (
            (args.reaction_colormap is not None, '--reaction-colormap'),
            (args.reaction_reverse_overlay, '--reaction-reverse-overlay')
        ) if present_flag]
        if reaction_flags:
            message = ', '.join(f"'{flag}'" for flag in reaction_flags)
            if args.reaction_txt is None:
                raise ConfigError(
                    f"Reaction-layer coloring options were given ({message}), but no reaction file "
                    f"('--reaction-txt') was provided."
                )
            if reaction['mode'] == 'single':
                raise ConfigError(
                    f"Reaction-layer coloring options were given ({message}), but the reaction "
                    f"layer is colored by single-color presence (the file has no value column and "
                    f"no 'sample' column), which does not use them."
                )

        compound_flags = [flag for present_flag, flag in (
            (args.compound_colormap is not None, '--compound-colormap'),
            (args.compound_reverse_overlay, '--compound-reverse-overlay')
        ) if present_flag]
        if compound_flags:
            message = ', '.join(f"'{flag}'" for flag in compound_flags)
            if args.compound_txt is None:
                raise ConfigError(
                    f"Compound-layer coloring options were given ({message}), but no compound file "
                    f"('--compound-txt') was provided."
                )
            if compound['mode'] == 'single':
                raise ConfigError(
                    f"Compound-layer coloring options were given ({message}), but the compound "
                    f"layer is colored by single-color presence (the file has no value column and "
                    f"no 'sample' column), which does not use them."
                )

        # Categorical values take colors of their own. A colormap of such a layer colors only its
        # presence on the 'unified' map. It has nothing to color where that map shows the values, or
        # where that map is not drawn. The 'unified' map is colored by the group summary with
        # groups, and by the sample summary without them. In a grouped run without a sample rule,
        # each group's map counts the group's samples, and '--group-colormap' colors those counts. A
        # rule across the groups in such a run is refused further below.
        unified_level = 'group' if args.groups_txt is not None else 'sample'
        for element_type, _, layer, sample_summary, group_summary in summaries:
            colormap_flag = f'--{element_type}-colormap'
            unified_summary = group_summary if args.groups_txt is not None else sample_summary
            if layer is None or layer['mode'] != 'categorical':
                continue
            if getattr(args, f'{element_type}_colormap') is None:
                continue
            group_counts = args.groups_txt is not None and sample_summary not in CATEGORICAL_RULES
            if not layer['has_sample']:
                reason = "The file has no 'sample' column, so its one map shows the values."
            elif group_counts and args.skip_unified_maps:
                raise ConfigError(
                    f"'{colormap_flag}' colors presence on the 'unified' map, but "
                    f"'--skip-unified-maps' draws no 'unified' map. Each group's map counts the "
                    f"group's samples, because '--{element_type}-sample-summary' gives no rule. "
                    f"Use '--group-colormap' to color those counts."
                )
            elif group_counts:
                continue
            elif unified_summary in CATEGORICAL_RULES:
                reason = (
                    f"'--{element_type}-{unified_level}-summary {unified_summary}' colors the "
                    f"'unified' map by the values too."
                )
            elif args.skip_unified_maps:
                reason = "'--skip-unified-maps' draws no 'unified' map."
            else:
                continue
            raise ConfigError(
                f"'{colormap_flag}' colors presence on the 'unified' map. Here no map of the "
                f"{element_type} layer is colored by presence. {reason} Categorical values take "
                f"their colors from '--{element_type}-value-colors', or from anvi'o's default "
                f"colors."
            )

        # Skipping the 'unified' map makes the maps of the individual samples or groups the whole
        # of the output. The run therefore needs samples. It also needs those maps to be asked for.
        if args.skip_unified_maps:
            if not any_sample:
                raise ConfigError(
                    "'--skip-unified-maps' skips the map summarizing all of the samples. No input "
                    "file has a 'sample' column, so there is only one map to draw. Skipping it "
                    "would leave an empty output directory. Add a 'sample' column to compare "
                    "origins, or drop this flag."
                )
            if draw_individual_files is False and draw_grid is False:
                raise ConfigError(
                    "'--skip-unified-maps' skips the map summarizing all of the samples. That "
                    "leaves the maps of the individual samples or groups to be drawn. None were "
                    "asked for. Please add '--draw-individual-files' and/or '--draw-grid', or drop "
                    "this flag."
                )

        # Individual maps and grids are drawn per sample, or per group of samples, so without a
        # 'sample' column there is nothing to draw.
        if (draw_individual_files is not False or draw_grid is not False) and not any_sample:
            flags = ', '.join(f"'{flag}'" for present_flag, flag in (
                (draw_individual_files is not False, '--draw-individual-files'),
                (draw_grid is not False, '--draw-grid')
            ) if present_flag)
            raise ConfigError(
                f"{flags} draws a map for each sample, or for each group of samples, but no input "
                f"file has a 'sample' column, so there is only one map to draw. Add a 'sample' "
                f"column to compare origins, or drop these options."
            )

        # Options belonging to a pangenome have nothing to act on in a text run.
        pan_flags = [flag for present_flag, flag in (
            (args.genomes_storage is not None, '--genomes-storage'),
            (args.consensus_threshold is not None, '--consensus-threshold'),
            (args.discard_ties, '--discard-ties')
        ) if present_flag]
        if pan_flags:
            message = ', '.join(f"'{flag}'" for flag in pan_flags)
            raise ConfigError(
                f"These options describe a pangenome: {message}. This run draws draw-kegg-pathways "
                f"text files instead, which have no gene clusters or genomes to interpret, so the "
                f"options cannot apply. Please remove them, or run against a '--pan-db' with its "
                f"'--genomes-storage' in a separate output directory."
            )

        # Each aggregation reduces a layer's values, so it needs that layer's file and its value
        # column; the gene level additionally needs the genes it reduces.
        for aggregation, flag, element_type, file_flag, layer in (
            (args.reaction_gene_aggregation, '--reaction-gene-aggregation', 'reaction',
             '--reaction-txt', reaction),
            (args.reaction_accession_aggregation, '--reaction-accession-aggregation', 'reaction',
             '--reaction-txt', reaction),
            (args.compound_accession_aggregation, '--compound-accession-aggregation', 'compound',
             '--compound-txt', compound)
        ):
            if aggregation is None:
                continue
            if layer is None:
                raise ConfigError(
                    f"'{flag}' reduces the values of the {element_type} layer, but no "
                    f"{element_type} file ('{file_flag}') was provided."
                )
            if layer['mode'] not in VALUE_MODES:
                raise ConfigError(
                    f"'{flag}' reduces the values of the {element_type} layer, but its file has no "
                    f"value column, so the layer is colored by presence and there are no values to "
                    f"reduce."
                )
            if aggregation in CATEGORICAL_RULES and layer['mode'] == 'quantitative':
                raise ConfigError(
                    f"'{flag}' was given as '{aggregation}', which chooses one categorical value "
                    f"among several. The values of the {element_type} layer are numbers, so they "
                    f"are reduced with an aggregation, such as 'sum' or 'mean'. To read them as "
                    f"categorical values, give '--{element_type}-value-type categorical'."
                )
            if flag == '--reaction-gene-aggregation' and not layer['has_gene_id']:
                raise ConfigError(
                    f"'{flag}' reduces the values of the genes annotated with one accession, but "
                    f"the reaction file has no 'gene_id' column naming genes. Each accession there "
                    f"already has a single value per sample, since rows repeating an accession are "
                    f"refused. Use '--reaction-accession-aggregation' to set how a map element "
                    f"combines the several accessions it stands for."
                )

        # A period describes a layer's values, so it needs that layer's file and its value column.
        # The mapper checks the number itself, as it does for programmatic calls.
        for period, flag, element_type, file_flag, layer in (
            (args.reaction_value_period, '--reaction-value-period', 'reaction', '--reaction-txt',
             reaction),
            (args.compound_value_period, '--compound-value-period', 'compound', '--compound-txt',
             compound)
        ):
            if period is None:
                continue
            if layer is None:
                raise ConfigError(
                    f"'{flag}' says how far the values of the {element_type} layer go before they "
                    f"repeat, but no {element_type} file ('{file_flag}') was provided."
                )
            if layer['mode'] == 'categorical':
                raise ConfigError(
                    f"'{flag}' says how far the values of the {element_type} layer go before they "
                    f"repeat. The values of that layer are categorical "
                    f"('--{element_type}-value-type categorical'). They are names, which do not "
                    f"repeat after a period."
                )
            if layer['mode'] != 'quantitative':
                raise ConfigError(
                    f"'{flag}' says how far the values of the {element_type} layer go before they "
                    f"repeat, but its file has no value column, so the layer is colored by "
                    f"presence and has no values."
                )

        # Each pair of value limits bounds, and each value center centers, a color scale of values,
        # so each needs its layer's file and that file's value column; the ones acting on the scale
        # the individual maps share need the samples those maps are drawn for as well. The mapper
        # checks everything past this point like a programmatic call, including the numbers
        # themselves and whether the context in question ends up colored by value once the summaries
        # have had their say.
        for setting, flag, element_type, file_flag, layer, per_individual, verb in (
            (args.reaction_value_limits, '--reaction-value-limits', 'reaction', '--reaction-txt',
             reaction, False, 'bound'),
            (args.reaction_individual_value_limits, '--reaction-individual-value-limits',
             'reaction', '--reaction-txt', reaction, True, 'bound'),
            (args.compound_value_limits, '--compound-value-limits', 'compound', '--compound-txt',
             compound, False, 'bound'),
            (args.compound_individual_value_limits, '--compound-individual-value-limits',
             'compound', '--compound-txt', compound, True, 'bound'),
            (args.reaction_value_center, '--reaction-value-center', 'reaction', '--reaction-txt',
             reaction, False, 'center'),
            (args.reaction_individual_value_center, '--reaction-individual-value-center',
             'reaction', '--reaction-txt', reaction, True, 'center'),
            (args.compound_value_center, '--compound-value-center', 'compound', '--compound-txt',
             compound, False, 'center'),
            (args.compound_individual_value_center, '--compound-individual-value-center',
             'compound', '--compound-txt', compound, True, 'center')
        ):
            if setting is None:
                continue
            if layer is None:
                raise ConfigError(
                    f"'{flag}' {verb}s the color scale of the {element_type} layer, but no "
                    f"{element_type} file ('{file_flag}') was provided."
                )
            if layer['mode'] == 'categorical':
                raise ConfigError(
                    f"'{flag}' {verb}s a color scale of numbers. The values of the {element_type} "
                    f"layer are categorical ('--{element_type}-value-type categorical'). Each "
                    f"takes a color of its own, so there is no scale to {verb}."
                )
            if layer['mode'] != 'quantitative':
                raise ConfigError(
                    f"'{flag}' {verb}s the color scale spanning the values of the {element_type} "
                    f"layer, but its file has no value column, so the layer is colored by presence "
                    f"and there is no scale of values to {verb}."
                )
            if per_individual and not layer['has_sample']:
                raise ConfigError(
                    f"'{flag}' {verb}s the color scale shared by the maps of the individual "
                    f"samples or groups, but the {element_type} file has no 'sample' column, so it "
                    f"draws no such maps and its values take a single scale. {verb.capitalize()} "
                    f"that one with "
                    f"'--{element_type}-value-{'limits' if verb == 'bound' else 'center'}'."
                )

        # A normalization needs the same things: the layer's file, a value column to rescale, and
        # the samples to rescale each value against. It also needs maps for those samples to be
        # drawn at all, which the mapper checks once it knows what this run draws.
        for normalization, element_type, file_flag, layer in (
            (args.reaction_element_normalization, 'reaction', '--reaction-txt', reaction),
            (args.compound_element_normalization, 'compound', '--compound-txt', compound)
        ):
            if normalization is None:
                continue
            flag = f'--{element_type}-element-normalization'
            if layer is None:
                raise ConfigError(
                    f"'{flag}' rescales the values of the {element_type} layer, but no "
                    f"{element_type} file ('{file_flag}') was provided."
                )
            if layer['mode'] == 'categorical':
                raise ConfigError(
                    f"'{flag}' rescales numbers. The values of the {element_type} layer are "
                    f"categorical ('--{element_type}-value-type categorical'). They are names, "
                    f"which cannot be rescaled."
                )
            if layer['mode'] != 'quantitative':
                raise ConfigError(
                    f"'{flag}' rescales the values of the {element_type} layer, but its file has "
                    f"no value column, so the layer is colored by presence and has no values to "
                    f"rescale."
                )
            if not layer['has_sample']:
                raise ConfigError(
                    f"'{flag}' rescales the value of each sample or group against the values of "
                    f"all of them, but the {element_type} file has no 'sample' column, so there is "
                    f"a single set of values with nothing to compare them against. Add a 'sample' "
                    f"column to compare samples."
                )

        # A colormap for the scale the individual maps share needs the same things that pair of
        # limits needs: the layer's file, a value column to put those maps on a scale, and the
        # samples they are drawn for. Whether that scale ends up colored by value once the summaries
        # have had their say is checked by the mapper, as for programmatic callers.
        for individual_colormap, element_type, file_flag, layer in (
            (args.reaction_individual_colormap, 'reaction', '--reaction-txt', reaction),
            (args.compound_individual_colormap, 'compound', '--compound-txt', compound)
        ):
            if individual_colormap is None:
                continue
            flag = f'--{element_type}-individual-colormap'
            if layer is None:
                raise ConfigError(
                    f"'{flag}' colors the maps of the individual samples or groups of the "
                    f"{element_type} layer, but no {element_type} file ('{file_flag}') was "
                    f"provided."
                )
            if layer['mode'] == 'categorical':
                raise ConfigError(
                    f"'{flag}' colors the maps of the individual samples or groups along a scale "
                    f"of numbers. The values of the {element_type} layer are categorical "
                    f"('--{element_type}-value-type categorical'). They take their colors from "
                    f"'--{element_type}-value-colors', or from anvi'o's default colors."
                )
            if layer['mode'] != 'quantitative':
                raise ConfigError(
                    f"'{flag}' colors the maps of the individual samples or groups of the "
                    f"{element_type} layer along a scale of values, but its file has no value "
                    f"column, so the layer is colored by presence instead. "
                    f"'--{element_type}-colormap' is what colors a presence scale."
                )
            if not layer['has_sample']:
                raise ConfigError(
                    f"'{flag}' colors the maps of the individual samples or groups, but the "
                    f"{element_type} file has no 'sample' column, so it draws no such maps and its "
                    f"values take a single scale. Color that one with '--{element_type}-colormap'."
                )

        # '--presence-colormap-scheme' is superseded for text layers by the summaries, which choose
        # the presence scheme per layer and per level.
        if args.presence_colormap_scheme is not None:
            raise ConfigError(
                f"'--presence-colormap-scheme' selects how presence is colored for contigs "
                f"databases, pan genomes, and groups of either. For a draw-kegg-pathways text "
                f"file, that choice is made per layer and per level by the summary options: use "
                f"'--reaction-sample-summary'/'--compound-sample-summary' with one of "
                f"{SUMMARY_PRESENCE_PHRASE} for coloring across samples, and "
                f"'--reaction-group-summary'/'--compound-group-summary' for coloring across "
                f"groups."
            )

        # Sample and group summaries. A summary reduces a set of samples, or the sample groups, to
        # one value or one presence state per map element. It needs its layer's file and a 'sample'
        # column. A summary that pools values also needs a value column.
        for element_type, file_flag, layer, sample_summary, group_summary in summaries:
            for level, summary in (('sample', sample_summary), ('group', group_summary)):
                if summary is None:
                    continue
                flag = f'--{element_type}-{level}-summary'
                if layer is None:
                    raise ConfigError(
                        f"'{flag}' summarizes the {element_type} layer, but no {element_type} file "
                        f"('{file_flag}') was provided."
                    )
                if not layer['has_sample']:
                    raise ConfigError(
                        f"'{flag}' summarizes the samples of the {element_type} layer, but its "
                        f"file has no 'sample' column, so there are no samples to summarize."
                    )
                if level == 'group' and args.groups_txt is None:
                    raise ConfigError(
                        f"'{flag}' summarizes groups of samples, but no groups were defined with "
                        f"'--groups-txt'."
                    )
                # A categorical rule chooses one categorical value across the samples or groups.
                # Every other name that is not presence is an aggregation, which pools numbers. A
                # rule given where it cannot apply is refused with what the user can do instead.
                # That depends on the file.
                if summary in CATEGORICAL_RULES and layer['mode'] not in VALUE_MODES:
                    raise ConfigError(
                        f"'{flag}' was given as '{summary}', which chooses one categorical value "
                        f"across the {level}s. The {element_type} file has no value column, so "
                        f"there are no values to choose among. Summarize presence with one of "
                        f"{SUMMARY_PRESENCE_PHRASE} instead."
                    )
                if summary in CATEGORICAL_RULES and layer['mode'] == 'quantitative':
                    raise ConfigError(
                        f"'{flag}' was given as '{summary}', which chooses one categorical value "
                        f"across the {level}s. The values of the {element_type} file are numbers. "
                        f"To read them as categorical values, give '--{element_type}-value-type "
                        f"categorical'."
                    )
                if (
                    layer['mode'] == 'categorical'
                    and summary not in PRESENCE_SUMMARIES + CATEGORICAL_RULES
                ):
                    # With groups, the samples of a group are counted unless a rule chooses among
                    # them, so presence names are no alternative at the sample level there.
                    if level == 'sample' and args.groups_txt is not None:
                        alternatives = (
                            f"Choose one value among the samples of each group with one of "
                            f"{CATEGORICAL_RULES_PHRASE}."
                        )
                    else:
                        alternatives = (
                            f"Summarize presence across the {level}s with one of "
                            f"{SUMMARY_PRESENCE_PHRASE}. Or choose one value across the {level}s "
                            f"with one of {CATEGORICAL_RULES_PHRASE}."
                        )
                    raise ConfigError(
                        f"'{flag}' was given as '{summary}'. The values of the {element_type} "
                        f"layer are categorical ('--{element_type}-value-type categorical'). They "
                        f"cannot be added up or averaged across the {level}s. {alternatives} "
                        f"Choosing assumes that a value means the same thing in every one of them."
                    )
                if (
                    summary not in PRESENCE_SUMMARIES + CATEGORICAL_RULES
                    and layer['mode'] != 'quantitative'
                ):
                    raise ConfigError(
                        f"'{flag}' was given as '{summary}', which is an aggregation that pools "
                        f"values, but the {element_type} file has no value column, so there are no "
                        f"values to pool. Summarize presence with one of {SUMMARY_PRESENCE_PHRASE} "
                        f"instead, or add a value column to the file."
                    )
                # With groups, presence at the SAMPLE level is always the count of a group's own
                # samples in discrete bands from '--group-colormap', so every presence name
                # describes the same per-group map and none of them changes anything: the choice
                # only shapes the ungrouped 'unified' map.
                if (
                    level == 'sample' and args.groups_txt is not None
                    and summary in PRESENCE_SUMMARIES
                ):
                    alternative = (
                        "give a rule such as 'plurality_order' to color each group's map by a "
                        "value chosen among its samples"
                    ) if layer['mode'] == 'categorical' else (
                        "use an aggregation such as 'mean' to color each group's map by pooled "
                        "values"
                    )
                    raise ConfigError(
                        f"'{flag}' was given as '{summary}', but the samples here are grouped with "
                        f"'--groups-txt', and each group's map shows the count of that group's own "
                        f"samples regardless of how presence is asked for — in discrete bands or "
                        f"on a continuous scale, as '--group-colormap-scheme' says — so "
                        f"{SUMMARY_PRESENCE_PHRASE} would all draw the same thing. Summarizing the "
                        f"samples by presence is already the default, so drop the option; "
                        f"{alternative} instead, and note that '--{element_type}-group-summary' is "
                        f"what chooses how presence is colored across the groups themselves."
                    )
                # A static color overrides summarizing across samples or groups, so asking for both
                # is asking for two different things on the 'unified' map.
                static_flag = f'--{element_type}-color'
                if getattr(args, f'{element_type}_color') is not None:
                    raise ConfigError(
                        f"'{flag}' summarizes the {element_type} layer's samples or groups, but "
                        f"'{static_flag}' colors that layer one static color for presence in any "
                        f"sample, which overrides the summary. Please use only one of them."
                    )
                # A summary reduces a set of samples or groups for the 'unified' map. With groups,
                # the sample summary also colors each group's own map, so only the group summary
                # has nothing left to do. Without groups, a per-sample map shows one sample, so
                # neither summary colors it.
                if args.skip_unified_maps and (level == 'group' or args.groups_txt is None):
                    elsewhere = (
                        f"Each group's own map is colored by '--{element_type}-sample-summary'. "
                        f"That option summarizes the group's own samples."
                    ) if level == 'group' else (
                        "The map of an individual sample shows that one sample alone. It shows "
                        "the sample's own values where the file has a value column, and its "
                        "presence where it does not."
                    )
                    raise ConfigError(
                        f"'{flag}' decides how the 'unified' map summarizes the {level}s of the "
                        f"{element_type} layer. Nothing else uses it. {elsewhere} "
                        f"'--skip-unified-maps' draws no 'unified' map to color. Please use "
                        f"'{flag}' or '--skip-unified-maps', not both."
                    )
            # A group's value comes from the sample summary, so a rule or an aggregation across the
            # groups needs one of the same kind among each group's samples.
            if (
                layer is not None and layer['mode'] == 'categorical'
                and group_summary in CATEGORICAL_RULES and sample_summary not in CATEGORICAL_RULES
            ):
                raise ConfigError(
                    f"'--{element_type}-group-summary' was given as '{group_summary}', which "
                    f"chooses one value among the values of the groups. Each group's value comes "
                    f"from '--{element_type}-sample-summary'. Without a rule, that option counts "
                    f"each group's samples, so the groups have no values to choose among. Please "
                    f"also give '--{element_type}-sample-summary' a rule, such as "
                    f"'plurality_order'."
                )
            if group_summary is not None and group_summary not in PRESENCE_SUMMARIES and (
                sample_summary is None or sample_summary in PRESENCE_SUMMARIES
            ) and group_summary not in CATEGORICAL_RULES:
                raise ConfigError(
                    f"'--{element_type}-group-summary' was given as '{group_summary}', which pools "
                    f"the values of the groups, but each group's value comes from "
                    f"'--{element_type}-sample-summary', which summarizes the group's samples by "
                    f"presence rather than by value. Please also set "
                    f"'--{element_type}-sample-summary' to an aggregation such as 'mean'."
                )

        # Group options.
        if args.groups_txt is not None and not any_sample:
            raise ConfigError(
                "'--groups-txt' groups the samples of the text files, but neither file has a "
                "'sample' column, so there are no samples to group."
            )
        if args.group_threshold is not None and not group_presence:
            raise ConfigError(
                "'--group-threshold' sets the proportion of a group's samples in which an element "
                "must occur for the group to be considered to contain it, so it applies when a "
                "layer's groups are summarized by presence. No layer's groups are summarized that "
                "way here: either no samples are grouped with '--groups-txt', or every sampled "
                "layer's '--*-group-summary' pools values or chooses a categorical value instead."
            )
        if args.group_threshold is not None and args.skip_unified_maps:
            raise ConfigError(
                "'--group-threshold' sets the proportion of a group's samples in which an element "
                "must occur for the group to count as containing it. Only the 'unified' map uses "
                "that threshold. A group's own map counts that group's own samples whatever the "
                "threshold is. '--skip-unified-maps' draws no 'unified' map. Please use "
                "'--group-threshold' or '--skip-unified-maps', not both."
            )
        if group_presence and args.group_threshold is None and not args.skip_unified_maps:
            raise ConfigError(
                "When grouping samples ('--groups-txt') whose groups are summarized by presence, "
                "which is the default, '--group-threshold' must also be given: the proportion of a "
                "group's samples in which an element must be present for the group to be "
                "considered to contain it. Alternatively, summarize the groups by pooling their "
                "values, or by choosing a categorical value with a rule, with "
                "'--reaction-group-summary'/'--compound-group-summary'."
            )
        group_flags = [flag for present_flag, flag in (
            (args.group_colormap is not None, '--group-colormap'),
            (args.group_colormap_scheme is not None, '--group-colormap-scheme'),
            (args.group_reverse_overlay, '--group-reverse-overlay')
        ) if present_flag]
        if group_flags and not (group_presence_maps and draw_individual_maps):
            message = ', '.join(f"'{flag}'" for flag in group_flags)
            raise ConfigError(
                f"Group-map coloring options were given ({message}), but they apply only to "
                f"individual group maps showing within-group sample counts, which are drawn with "
                f"'--groups-txt' plus '--draw-individual-files' and/or '--draw-grid' for a layer "
                f"whose samples are summarized by presence. A layer whose samples are summarized "
                f"by pooling values colors its group maps continuously from "
                f"'--reaction-individual-colormap'/'--compound-individual-colormap' (or the "
                f"layer's own '--reaction-colormap'/'--compound-colormap') and "
                f"'--reaction-reverse-overlay'/'--compound-reverse-overlay' instead. A layer of "
                f"categorical values whose samples a rule summarizes colors its group maps by "
                f"those values."
            )

        # Each slice takes the color of its own value, so slicing needs a reaction layer of values
        # and a drawn map that shows them. Without a 'sample' column, every map shows the values.
        # With one, the map of each sample shows them. The map of a group and the 'unified' map show
        # them only where their summary pools values or chooses a value with a rule. Maps that show
        # presence are drawn whole. The mapper repeats these checks for programmatic callers.
        if args.reaction_slice_by is not None:
            if reaction['mode'] not in VALUE_MODES:
                raise ConfigError(
                    "'--reaction-slice-by' gives each slice the color of its own value, but the "
                    "reaction file has no value column. Its layer is colored by presence, and "
                    "slices of presence are not supported yet. Add a value column, or drop "
                    "'--reaction-slice-by'."
                )
            if args.reaction_element_normalization is not None:
                raise ConfigError(
                    "'--reaction-slice-by' cannot be combined with "
                    "'--reaction-element-normalization' yet. Rescaling each slice on its own is "
                    "not supported. Please use only one."
                )
            # When slicing by gene, the accession aggregation reduces the rows of one gene for the
            # accessions of an element. It is refused for numbers only when slicing by accession.
            if (
                args.reaction_accession_aggregation is not None
                and args.reaction_slice_by == 'accession'
                and reaction['mode'] == 'quantitative'
            ):
                raise ConfigError(
                    "'--reaction-accession-aggregation' reduces the values of the accessions of a "
                    "map element to one value. With '--reaction-slice-by accession', each slice "
                    "shows the value of one accession. An element drawn whole takes the mean of "
                    "its slice values, weighted by their shares, in place of this reduction. "
                    "Please drop '--reaction-accession-aggregation'."
                )
            if args.reaction_gene_aggregation is not None and args.reaction_slice_by == 'gene':
                raise ConfigError(
                    "'--reaction-gene-aggregation' reduces the values of the genes annotated with "
                    "one accession to one value. With '--reaction-slice-by gene', each gene is a "
                    "slice of its own, so no values of genes are reduced to an accession. "
                    "'--reaction-accession-aggregation' reduces the rows of a gene for two or more "
                    "accessions of an element. Please drop '--reaction-gene-aggregation', or give "
                    "'--reaction-slice-by accession' to slice elements by their accessions."
                )

            has_sample = reaction['has_sample']
            grouped = has_sample and args.groups_txt is not None
            sample_summary = args.reaction_sample_summary
            unified_summary = args.reaction_group_summary if grouped else sample_summary
            group_values = sample_summary is not None and sample_summary not in PRESENCE_SUMMARIES
            individual_values = draw_individual_maps and (not grouped or group_values)
            unified_values = not args.skip_unified_maps and (
                not has_sample
                or (unified_summary is not None and unified_summary not in PRESENCE_SUMMARIES)
            )
            if reaction['mode'] == 'categorical':
                value_summary = "a rule such as 'plurality_order'"
            elif args.reaction_value_period is not None:
                # With a period, 'circular_mean' is the only average that the summaries accept.
                value_summary = "an aggregation such as 'circular_mean'"
            else:
                value_summary = "an aggregation such as 'mean'"
            draw_flags = "'--draw-individual-files' or '--draw-grid'"
            if not individual_values and not unified_values:
                if not grouped:
                    reason = (
                        "The 'unified' map shows presence, since '--reaction-sample-summary' "
                        "summarizes presence. The map of each sample shows values, but none were "
                        "asked for."
                    )
                    alternative = (
                        f"Add {draw_flags}, or give '--reaction-sample-summary' {value_summary}."
                    )
                elif group_values:
                    reason = (
                        "The 'unified' map shows presence, since '--reaction-group-summary' "
                        "summarizes presence. The map of each group shows values, but none were "
                        "asked for."
                    )
                    alternative = (
                        f"Add {draw_flags}, or give '--reaction-group-summary' {value_summary}."
                    )
                else:
                    reason = (
                        "The map of each group counts its samples, since "
                        "'--reaction-sample-summary' summarizes presence. The 'unified' map then "
                        "shows presence across the groups."
                    )
                    alternative = (
                        f"To show values on the map of each group, give "
                        f"'--reaction-sample-summary' {value_summary}."
                    )
                    if not draw_individual_maps:
                        alternative += (
                            f" Then either draw those maps with {draw_flags}, or also give "
                            f"'--reaction-group-summary' {value_summary} to show values on the "
                            f"'unified' map."
                        )
                raise ConfigError(
                    f"'--reaction-slice-by' slices the reaction elements of maps that show values. "
                    f"Maps that show presence are drawn whole. Here no drawn map of the reaction "
                    f"layer shows values. {reason} {alternative}"
                )

            # Group maps and the 'unified' map of a layer with samples pool the shares of slices.
            # The map of one sample, and the one map of a file without samples, take their shares
            # from the weights directly. In a grouped run every map pools samples, so the check
            # above already found one that shows values.
            pooled_values = (grouped and individual_values) or (has_sample and unified_values)
            if args.reaction_slice_share_summary is not None and not pooled_values:
                if not has_sample:
                    reason = (
                        "The reaction file has no 'sample' column, so no map pools samples. Each "
                        "map takes its shares from the weights in the file."
                    )
                    alternative = "Please drop '--reaction-slice-share-summary'."
                elif args.skip_unified_maps:
                    reason = (
                        "Without groups, only the 'unified' map pools the samples, and "
                        "'--skip-unified-maps' draws no 'unified' map. The map of each sample "
                        "takes its shares from the weights of that sample."
                    )
                    alternative = "Please drop '--reaction-slice-share-summary'."
                else:
                    reason = (
                        "Without groups, only the 'unified' map pools the samples. It shows "
                        "presence, since '--reaction-sample-summary' summarizes presence. Its "
                        "elements are therefore drawn whole. The map of each sample takes its "
                        "shares from the weights of that sample."
                    )
                    alternative = (
                        f"Give '--reaction-sample-summary' {value_summary} to show values on the "
                        f"'unified' map, or drop '--reaction-slice-share-summary'."
                    )
                raise ConfigError(
                    f"'--reaction-slice-share-summary' says how the maps of groups and the "
                    f"'unified' map pool the shares of slices across samples. Here no such map is "
                    f"drawn with values. {reason} {alternative}"
                )

    reaction_colormap, reaction_colormap_limits = _resolve_colormap(
        args.reaction_colormap, '--reaction-colormap'
    )
    compound_colormap, compound_colormap_limits = _resolve_colormap(
        args.compound_colormap, '--compound-colormap'
    )
    reaction_individual_colormap, reaction_individual_colormap_limits = _resolve_colormap(
        args.reaction_individual_colormap, '--reaction-individual-colormap'
    )
    compound_individual_colormap, compound_individual_colormap_limits = _resolve_colormap(
        args.compound_individual_colormap, '--compound-individual-colormap'
    )
    group_colormap, group_colormap_limits = _resolve_colormap(
        args.group_colormap, '--group-colormap'
    )
    reaction_element_normalization, reaction_element_normalization_label = (
        _resolve_element_normalization(
            args.reaction_element_normalization, '--reaction-element-normalization'
        )
    )
    compound_element_normalization, compound_element_normalization_label = (
        _resolve_element_normalization(
            args.compound_element_normalization, '--compound-element-normalization'
        )
    )

    # Presence colors: '--original-color' routes the reaction layer through the reference-color
    # drawer; otherwise '--reaction-color'/'--compound-color' set each layer's single presence color
    # (None leaves the layer's default). Their bare-flag values are the defaults, so forwarding them
    # is harmless.
    reaction_color = 'original' if args.original_color else args.reaction_color
    compound_color = args.compound_color

    kwargs = {
        'reaction_txt': args.reaction_txt,
        'compound_txt': args.compound_txt,

        'reaction_sample_summary': args.reaction_sample_summary,
        'compound_sample_summary': args.compound_sample_summary,
        'reaction_group_summary': args.reaction_group_summary,
        'compound_group_summary': args.compound_group_summary,
        'groups_txt': args.groups_txt,
        'group_threshold': args.group_threshold,
        'pathway_numbers': args.pathway_numbers,
        'draw_unified_maps': not args.skip_unified_maps,
        'draw_individual_files': draw_individual_files,
        'draw_grid': draw_grid,
        'reaction_reverse_overlay': args.reaction_reverse_overlay,
        'compound_reverse_overlay': args.compound_reverse_overlay,
        'reaction_value_limits': _resolve_value_limits(args.reaction_value_limits),
        'reaction_individual_value_limits': _resolve_value_limits(
            args.reaction_individual_value_limits
        ),
        'compound_value_limits': _resolve_value_limits(args.compound_value_limits),
        'compound_individual_value_limits': _resolve_value_limits(
            args.compound_individual_value_limits
        ),
        # The centers go to the mapper as they were typed, bare flag included
        # ('VALUE_CENTER_DEFAULT' standing in for the value it did not carry), so that the number is
        # converted and checked in the one place programmatic callers reach too.
        'reaction_value_center': args.reaction_value_center,
        'reaction_individual_value_center': args.reaction_individual_value_center,
        'compound_value_center': args.compound_value_center,
        'compound_individual_value_center': args.compound_individual_value_center,
        # The periods go to the mapper as they were typed, as the centers do.
        'reaction_value_period': args.reaction_value_period,
        'compound_value_period': args.compound_value_period,
        'reaction_element_normalization': reaction_element_normalization,
        'reaction_element_normalization_label': reaction_element_normalization_label,
        'compound_element_normalization': compound_element_normalization,
        'compound_element_normalization_label': compound_element_normalization_label,
        'group_reverse_overlay': args.group_reverse_overlay,
        'group_colormap_scheme': args.group_colormap_scheme,
        'draw_maps_lacking_data': args.draw_bare_maps
    }
    for aggregation, key in (
        (args.reaction_gene_aggregation, 'reaction_gene_aggregation'),
        (args.reaction_accession_aggregation, 'reaction_accession_aggregation'),
        (args.compound_accession_aggregation, 'compound_accession_aggregation')
    ):
        if aggregation is not None:
            kwargs[key] = aggregation
    if reaction_colormap is not None:
        kwargs['reaction_colormap'] = reaction_colormap
        kwargs['reaction_colormap_limits'] = reaction_colormap_limits
    if compound_colormap is not None:
        kwargs['compound_colormap'] = compound_colormap
        kwargs['compound_colormap_limits'] = compound_colormap_limits
    # Left unset, each layer's per-sample/per-group scale is colored from that layer's own colormap,
    # which the mapper arranges, so only a colormap actually given is forwarded.
    if reaction_individual_colormap is not None:
        kwargs['reaction_individual_colormap'] = reaction_individual_colormap
        kwargs['reaction_individual_colormap_limits'] = reaction_individual_colormap_limits
    if compound_individual_colormap is not None:
        kwargs['compound_individual_colormap'] = compound_individual_colormap
        kwargs['compound_individual_colormap_limits'] = compound_individual_colormap_limits
    if group_colormap is not None:
        kwargs['group_colormap'] = group_colormap
        # Limits left unset are resolved by the mapper, which defaults them differently for a named
        # colormap (the fraction of it to use) and for a ramp built from a group's own color (how
        # far from white it runs).
        kwargs['group_colormap_limits'] = group_colormap_limits
    for membership_colors, key in (
        (args.reaction_membership_colors, 'reaction_membership_colors'),
        (args.compound_membership_colors, 'compound_membership_colors')
    ):
        if membership_colors is not None:
            kwargs[key] = membership_colors
    # A colormap of False tells the mapper that a static color was chosen for the layer, so presence
    # in any sample is drawn in that one color on the 'unified' map rather than being colored by
    # sample or group, matching how the same choice works for contigs databases and pangenomes.
    if reaction_color is not None:
        kwargs['reaction_color'] = reaction_color
        if not args.original_color:
            kwargs['reaction_colormap'] = False
    if compound_color is not None:
        kwargs['compound_color'] = compound_color
        kwargs['compound_colormap'] = False

    # A value type, the colors of categorical values, and the seed go to the mapper only where
    # given. The mapper refuses a seed where no rule chooses at random.
    for element_type in ('reaction', 'compound'):
        value_type = getattr(args, f'{element_type}_value_type')
        if value_type is not None:
            kwargs[f'{element_type}_value_type'] = value_type
        value_colors = getattr(args, f'{element_type}_value_colors')
        if value_colors is not None:
            kwargs[f'{element_type}_value_colors'] = value_colors
    if args.seed is not None:
        kwargs['seed'] = args.seed

    # The slicing options go to the mapper only where given. The mapper resolves their defaults.
    for setting, key in (
        (args.reaction_slice_by, 'reaction_slice_by'),
        (args.reaction_slice_weights, 'reaction_slice_weights'),
        (args.reaction_slice_weight_aggregation, 'reaction_slice_weight_aggregation'),
        (args.reaction_slice_share_summary, 'reaction_slice_share_summary')
    ):
        if setting is not None:
            kwargs[key] = setting

    kwargs['count_scale_max'] = args.count_scale_max

    mapper.map_kegg_pathways_txt(args.output_dir, **kwargs)

def map_single_contigs_db_ko_data(args: Namespace, mapper: Mapper) -> None:
    """Draw KO data from a single contigs database source in the absence of a colormap."""
    map_contigs_database_kos = mapper.map_contigs_database_kos

    # A single database is always drawn in one static color: '--original-color' uses the reference
    # colors, '--reaction-color' a chosen color (its bare-flag value is the default, so it need not
    # be forwarded), otherwise the default color.
    reaction_color = 'original' if args.original_color else args.reaction_color
    if reaction_color is not None:
        map_contigs_database_kos = functools.partial(
            map_contigs_database_kos, reaction_color=reaction_color
        )

    map_contigs_database_kos(
        args.contigs_dbs[0],
        args.output_dir,
        pathway_numbers=args.pathway_numbers,
        draw_maps_lacking_data=args.draw_bare_maps
    )

def map_multiple_contigs_dbs_ko_data(args: Namespace, mapper: Mapper) -> None:
    """Draw KO data from contigs database sources."""
    map_contigs_databases_kos = mapper.map_contigs_databases_kos

    if args.draw_individual_files is None:
        pass
    elif len(args.draw_individual_files) == 0:
        # Draw maps for all contigs databases or groups.
        map_contigs_databases_kos = functools.partial(
            map_contigs_databases_kos, draw_individual_files=True
        )
    else:
        # Draw maps for select contigs databases or groups.
        map_contigs_databases_kos = functools.partial(
            map_contigs_databases_kos, draw_individual_files=args.draw_individual_files
        )

    if args.draw_grid is None:
        pass
    elif len(args.draw_grid) == 0:
        # Draw a grid of maps including all contigs databases or all groups.
        map_contigs_databases_kos = functools.partial(map_contigs_databases_kos, draw_grid=True)
    else:
        # Include select contigs databases or select groups.
        map_contigs_databases_kos = functools.partial(
            map_contigs_databases_kos, draw_grid=args.draw_grid
        )

    assert not (
        (args.original_color or args.reaction_color is not None)
        and (args.reaction_colormap is not None)
    )

    if args.reaction_colormap is None:
        # Dynamically color reactions by database or group in unified maps using the default
        # colormap.
        pass
    elif len(args.reaction_colormap) == 1:
        # Use the provided colormap name.
        map_contigs_databases_kos = functools.partial(
            map_contigs_databases_kos, reaction_colormap=args.reaction_colormap[0]
        )
    else:
        # Use the provided colormap name and limits.
        assert len(args.reaction_colormap) == 3
        min_limit = float(args.reaction_colormap[1])
        max_limit = float(args.reaction_colormap[2])
        map_contigs_databases_kos = functools.partial(
            map_contigs_databases_kos,
            reaction_colormap=args.reaction_colormap[0],
            reaction_colormap_limits=(min_limit, max_limit)
        )

    if args.presence_colormap_scheme is None:
        # The scheme is determined automatically by the number of contigs databases or groups.
        pass
    else:
        map_contigs_databases_kos = functools.partial(
            map_contigs_databases_kos, colormap_scheme=args.presence_colormap_scheme
        )

    if args.reaction_membership_colors is not None:
        # Color reactions by membership from a color given per database, genome, or group, rather
        # than from a colormap.
        map_contigs_databases_kos = functools.partial(
            map_contigs_databases_kos, reaction_membership_colors=args.reaction_membership_colors
        )

    if args.original_color:
        # Preserve the reference map's colors, statically (overriding dynamic coloring).
        map_contigs_databases_kos = functools.partial(
            map_contigs_databases_kos, reaction_colormap=False, reaction_color='original'
        )
    elif args.reaction_color is not None:
        # Color reactions one static color by presence/absence, overriding dynamic coloring (the
        # bare '--reaction-color' value is the default color).
        map_contigs_databases_kos = functools.partial(
            map_contigs_databases_kos, reaction_colormap=False, reaction_color=args.reaction_color
        )
    # else: dynamically color reactions by database or group count or membership (the default).

    if args.groups_txt:
        if args.group_colormap is None:
            # Dynamically color reactions by group in individual group maps using the default
            # colormap.
            pass
        elif len(args.group_colormap) == 1:
            # Use the provided group colormap name.
            map_contigs_databases_kos = functools.partial(
                map_contigs_databases_kos, group_colormap=args.group_colormap[0]
            )
        else:
            # Use the provided group colormap name and limits.
            assert len(args.group_colormap) == 3
            min_limit = float(args.group_colormap[1])
            max_limit = float(args.group_colormap[2])
            map_contigs_databases_kos = functools.partial(
                map_contigs_databases_kos,
                group_colormap=args.group_colormap[0],
                group_colormap_limits=(min_limit, max_limit)
            )

    map_contigs_databases_kos(
        args.contigs_dbs,
        args.output_dir,
        groups_txt=args.groups_txt,
        group_threshold=args.group_threshold,
        pathway_numbers=args.pathway_numbers,
        draw_unified_maps=not args.skip_unified_maps,
        reaction_reverse_overlay=args.reaction_reverse_overlay,
        group_reverse_overlay=args.group_reverse_overlay,
        group_colormap_scheme=args.group_colormap_scheme,
        count_scale_max=args.count_scale_max,
        draw_maps_lacking_data=args.draw_bare_maps
    )

def map_pan_db_ko_data(args: Namespace, mapper: Mapper) -> None:
    """Draw KO data from a pangenomic database source."""
    map_pan_database_kos = mapper.map_pan_database_kos

    if args.draw_individual_files is None:
        pass
    elif len(args.draw_individual_files) == 0:
        # Draw maps for all genomes or groups.
        map_pan_database_kos = functools.partial(map_pan_database_kos, draw_individual_files=True)
    else:
        # Draw maps for select genomes or groups.
        map_pan_database_kos = functools.partial(
            map_pan_database_kos, draw_individual_files=args.draw_individual_files
        )

    if args.draw_grid is None:
        pass
    elif len(args.draw_grid) == 0:
        # Draw a grid of maps including all genomes or all groups.
        map_pan_database_kos = functools.partial(map_pan_database_kos, draw_grid=True)
    else:
        # Include select genomes or select groups.
        map_pan_database_kos = functools.partial(
            map_pan_database_kos, draw_grid=args.draw_grid
        )

    assert not (
        (args.original_color or args.reaction_color is not None)
        and (args.reaction_colormap is not None)
    )

    if args.reaction_colormap is None:
        # Dynamically color reactions by genome or group in 'unified' maps using the default
        # colormap.
        pass
    elif len(args.reaction_colormap) == 1:
        # Use the provided colormap name.
        map_pan_database_kos = functools.partial(
            map_pan_database_kos, reaction_colormap=args.reaction_colormap[0]
        )
    else:
        # Use the provided colormap name and limits.
        assert len(args.reaction_colormap) == 3
        min_limit = float(args.reaction_colormap[1])
        max_limit = float(args.reaction_colormap[2])
        map_pan_database_kos = functools.partial(
            map_pan_database_kos,
            reaction_colormap=args.reaction_colormap[0],
            reaction_colormap_limits=(min_limit, max_limit)
        )

    if args.presence_colormap_scheme is None:
        # The scheme is determined automatically by the number of genomes or groups.
        pass
    else:
        map_pan_database_kos = functools.partial(
            map_pan_database_kos, colormap_scheme=args.presence_colormap_scheme
        )

    if args.reaction_membership_colors is not None:
        # Color reactions by membership from a color given per database, genome, or group, rather
        # than from a colormap.
        map_pan_database_kos = functools.partial(
            map_pan_database_kos, reaction_membership_colors=args.reaction_membership_colors
        )

    if args.original_color:
        # Preserve the reference map's colors, statically (overriding dynamic coloring).
        map_pan_database_kos = functools.partial(
            map_pan_database_kos, reaction_colormap=False, reaction_color='original'
        )
    elif args.reaction_color is not None:
        # Color reactions one static color by presence/absence, overriding dynamic coloring (the
        # bare '--reaction-color' value is the default color).
        map_pan_database_kos = functools.partial(
            map_pan_database_kos, reaction_colormap=False, reaction_color=args.reaction_color
        )
    # else: dynamically color reactions by genome or group count or membership (the default).

    if args.groups_txt:
        if args.group_colormap is None:
            # Dynamically color reactions by group in individual group maps using the default
            # colormap.
            pass
        elif len(args.group_colormap) == 1:
            # Use the provided group colormap name.
            map_pan_database_kos = functools.partial(
                map_pan_database_kos, group_colormap=args.group_colormap[0]
            )
        else:
            # Use the provided group colormap name and limits.
            assert len(args.group_colormap) == 3
            min_limit = float(args.group_colormap[1])
            max_limit = float(args.group_colormap[2])
            map_pan_database_kos = functools.partial(
                map_pan_database_kos,
                group_colormap=args.group_colormap[0],
                group_colormap_limits=(min_limit, max_limit)
            )

    map_pan_database_kos(
        args.pan_db,
        args.genomes_storage,
        args.output_dir,
        groups_txt=args.groups_txt,
        group_threshold=args.group_threshold,
        consensus_threshold=args.consensus_threshold,
        discard_ties=args.discard_ties,
        pathway_numbers=args.pathway_numbers,
        draw_unified_maps=not args.skip_unified_maps,
        reaction_reverse_overlay=args.reaction_reverse_overlay,
        group_reverse_overlay=args.group_reverse_overlay,
        group_colormap_scheme=args.group_colormap_scheme,
        count_scale_max=args.count_scale_max,
        draw_maps_lacking_data=args.draw_bare_maps
    )

def main() -> None:
    args = get_args()

    try:
        check_kegg_data(args)
        consolidate_contigs_dbs(args)

        # '--count-scale-max' either names a rule or gives the count to stop at. Normalized here,
        # once, so that every dispatcher below hands the mapper either a rule name or an int.
        if args.count_scale_max not in COUNT_SCALE_MAX_RULES:
            try:
                args.count_scale_max = int(args.count_scale_max)
            except ValueError:
                raise ConfigError(
                    f"'--count-scale-max' was given as '{args.count_scale_max}', which is neither "
                    f"{' nor '.join(repr(rule) for rule in COUNT_SCALE_MAX_RULES)} nor a number. "
                    f"Use 'observed' to stop a count scale at the highest count in the data. Use "
                    f"'total' to run it to the number of samples, contigs databases, genomes, or "
                    f"groups being counted. Or give a number to stop it there."
                )
            if args.count_scale_max < 1:
                raise ConfigError(
                    f"'--count-scale-max' was given as {args.count_scale_max}, but a scale of "
                    f"counts has to reach at least 1. An element colored by count is in at least "
                    f"one sample, contigs database, genome, or group."
                )

        # Only one primary data source can be drawn per run. Each source writes its map files (and
        # shared colorbars) into the same output directory with the same names, so combining sources
        # would collide or silently overwrite. The reaction- and compound-layer text files are one
        # source pair (they are drawn together on each map), so they count as a single source.
        is_txt_run = args.reaction_txt is not None or args.compound_txt is not None
        provided_sources: list = []
        if is_txt_run:
            provided_sources.append('--reaction-txt/--compound-txt')
        if args.reaction_network_json is not None:
            provided_sources.append('--reaction-network-json')
        if args.contigs_dbs is not None:
            provided_sources.append('--contigs-dbs/--external-genomes')
        # A genomes storage on its own is not enough to draw a pangenome, but naming one alongside
        # another source is still asking for two sources, so it counts here.
        if args.pan_db is not None or args.genomes_storage is not None:
            provided_sources.append('--pan-db/--genomes-storage')
        if len(provided_sources) > 1:
            message = ', '.join(f"'{source}'" for source in provided_sources)
            raise ConfigError(
                f"Only one data source can be drawn per run, but more than one was provided: "
                f"{message}. Each source writes its maps and colorbar into the same output "
                f"directory, so please run the program once per source, into separate output "
                f"directories."
            )

        # Static presence colors vs. dynamic colormaps: a layer is colored one way or the other.
        # '--original-color' preempts custom colors (it uses the reference map's colors); the color
        # flags preempt their layer's colormap (several dispatchers below also assert this).
        if args.original_color and (
            args.reaction_color is not None or args.compound_color is not None
        ):
            raise ConfigError(
                "'--original-color' uses the reference map's original colors, so it cannot be "
                "combined with '--reaction-color' or '--compound-color'."
            )
        if (args.original_color or args.reaction_color is not None) and (
            args.reaction_colormap is not None
        ):
            raise ConfigError(
                "'--reaction-color'/'--original-color' set a static reaction color, while "
                "'--reaction-colormap' colors reactions dynamically across a range. Please use "
                "only one."
            )
        if args.compound_color is not None and args.compound_colormap is not None:
            raise ConfigError(
                "'--compound-color' sets a static compound color, while '--compound-colormap' "
                "colors compounds dynamically across a range. Please use only one."
            )

        # A membership colors file and a colormap answer the same question. The question is where a
        # layer's colors come from, so a layer takes one or the other. A static color, or the
        # reference map's own colors, answers it a third way.
        for element_type in ('reaction', 'compound'):
            if getattr(args, f'{element_type}_membership_colors') is None:
                continue
            colors_flag = f'--{element_type}-membership-colors'
            if getattr(args, f'{element_type}_colormap') is not None:
                raise ConfigError(
                    f"'{colors_flag}' colors the {element_type} layer with the colors given in the "
                    f"file. '--{element_type}-colormap' has anvi'o choose the colors of the layer "
                    f"from a colormap instead. Please use only one."
                )
            if getattr(args, f'{element_type}_color') is not None:
                raise ConfigError(
                    f"'{colors_flag}' colors the {element_type} layer with the colors given in the "
                    f"file. '--{element_type}-color' colors the layer one static color for "
                    f"presence anywhere. Please use only one."
                )
            if args.original_color:
                raise ConfigError(
                    f"'{colors_flag}' colors the {element_type} layer with the colors given in the "
                    f"file. '--original-color' takes both the colors and the drawing order of the "
                    f"layer from the reference map, so there is no color left to choose. Please "
                    f"use only one."
                )

        # With one contigs database, membership has a single combination, and so a single color.
        # That is a static color rather than a scale. The run takes the single-database path. That
        # path has no set of databases for a membership colors file to give colors to.
        if args.reaction_membership_colors is not None and args.contigs_dbs is not None and (
            len(args.contigs_dbs) == 1
        ):
            raise ConfigError(
                "'--reaction-membership-colors' gives a color to each contigs database or group of "
                "them, which colors reactions by exactly which of them contain a reaction. Only "
                "one contigs database was provided, so every reaction it contains is in the same "
                "single database and there is nothing to tell apart by color. Please provide more "
                "databases, or set the color of this one with '--reaction-color'."
            )

        # A ramp per group runs to that group's own color, so the colors have to come from
        # somewhere. The layer-level checks in 'map_txt_data' and below decide whether the file
        # given is one this run can use; this only catches asking for the ramp without any file.
        if args.group_colormap is not None and args.group_colormap[0] == (
            GROUP_COLORMAP_TINT
        ) and args.reaction_membership_colors is None and args.compound_membership_colors is None:
            # A database, pangenome or JSON run has a reaction layer alone, so pointing it at
            # '--compound-membership-colors' would only earn it a second error.
            colors_flags = (
                "'--reaction-membership-colors'/'--compound-membership-colors'" if is_txt_run
                else "'--reaction-membership-colors'"
            )
            raise ConfigError(
                f"'--group-colormap {GROUP_COLORMAP_TINT}' colors each group's own maps by a ramp "
                f"running from a pale tint to that group's own color, but no color was given for "
                f"any group. Give one per group with {colors_flags}, or name a Matplotlib colormap "
                f"for the group maps instead."
            )
        # Only a compound text file provides a compound layer, so every option that shapes one has
        # nothing to act on without it. The same goes for the sample and group summaries and the
        # value limits, which only the text engine reads: a value column is what puts elements on a
        # continuous scale, and only a text file has one. The options for categorical values and the
        # seed of their random choices need a value column too, as do the options that slice
        # reaction elements.
        if not is_txt_run:
            compound_only_flags = [flag for present_flag, flag in (
                (args.compound_color is not None, '--compound-color'),
                (args.compound_colormap is not None, '--compound-colormap'),
                (args.compound_membership_colors is not None, '--compound-membership-colors'),
                (args.compound_reverse_overlay, '--compound-reverse-overlay'),
                (args.compound_accession_aggregation is not None,
                 '--compound-accession-aggregation'),
                (args.reaction_gene_aggregation is not None, '--reaction-gene-aggregation'),
                (args.reaction_accession_aggregation is not None,
                 '--reaction-accession-aggregation'),
                (args.reaction_value_period is not None, '--reaction-value-period'),
                (args.compound_value_period is not None, '--compound-value-period'),
                (args.reaction_sample_summary is not None, '--reaction-sample-summary'),
                (args.compound_sample_summary is not None, '--compound-sample-summary'),
                (args.reaction_group_summary is not None, '--reaction-group-summary'),
                (args.compound_group_summary is not None, '--compound-group-summary'),
                (args.reaction_value_limits is not None, '--reaction-value-limits'),
                (args.reaction_individual_value_limits is not None,
                 '--reaction-individual-value-limits'),
                (args.compound_value_limits is not None, '--compound-value-limits'),
                (args.compound_individual_value_limits is not None,
                 '--compound-individual-value-limits'),
                (args.reaction_value_center is not None, '--reaction-value-center'),
                (args.reaction_individual_value_center is not None,
                 '--reaction-individual-value-center'),
                (args.compound_value_center is not None, '--compound-value-center'),
                (args.compound_individual_value_center is not None,
                 '--compound-individual-value-center'),
                (args.reaction_individual_colormap is not None, '--reaction-individual-colormap'),
                (args.compound_individual_colormap is not None, '--compound-individual-colormap'),
                (args.reaction_element_normalization is not None,
                 '--reaction-element-normalization'),
                (args.compound_element_normalization is not None,
                 '--compound-element-normalization'),
                (args.reaction_value_type is not None, '--reaction-value-type'),
                (args.compound_value_type is not None, '--compound-value-type'),
                (args.reaction_value_colors is not None, '--reaction-value-colors'),
                (args.compound_value_colors is not None, '--compound-value-colors'),
                (args.seed is not None, '--seed'),
                (args.reaction_slice_by is not None, '--reaction-slice-by'),
                (args.reaction_slice_weights is not None, '--reaction-slice-weights'),
                (args.reaction_slice_weight_aggregation is not None,
                 '--reaction-slice-weight-aggregation'),
                (args.reaction_slice_share_summary is not None, '--reaction-slice-share-summary')
            ) if present_flag]
            if compound_only_flags:
                message = ', '.join(f"'{flag}'" for flag in compound_only_flags)
                raise ConfigError(
                    f"These options were given: {message}. They apply only to a text file given "
                    f"with '--reaction-txt' or '--compound-txt'. Database, pangenome, and "
                    f"reaction-network-JSON inputs have only a reaction layer. That layer is "
                    f"colored by how many databases or genomes contain an element. It has no "
                    f"samples of a text file and no column of values, so these options have "
                    f"nothing to act on."
                )
        # '--original-color' draws only the reaction layer (its compounds follow in the reference
        # colors on global/overview maps); it cannot also stage an explicit compound file.
        if args.original_color and args.compound_txt is not None:
            raise ConfigError(
                "'--original-color' draws only the reaction layer, through a dedicated renderer "
                "that also colors those reactions' compounds in the reference colors (on global "
                "and overview maps), so it cannot also take an explicit compound file "
                "('--compound-txt')."
            )

        # The group-map coloring options style the individual maps of each group, which are drawn
        # only for a grouped run that asks for them. Rather than accept them and do nothing, refuse
        # them here. A text run refuses them in 'map_txt_data'. A reaction-network-JSON run refuses
        # the whole family in its own dispatcher, since it has no samples, databases, genomes, or
        # groups at all. Note that a static reaction color does not make them idle: a group's own
        # map still counts the group's databases or genomes, and still takes a colorbar of its own.
        if args.contigs_dbs is not None or args.pan_db is not None:
            group_flags = [flag for present_flag, flag in (
                (args.group_colormap is not None, '--group-colormap'),
                (args.group_colormap_scheme is not None, '--group-colormap-scheme'),
                (args.group_reverse_overlay, '--group-reverse-overlay')
            ) if present_flag]
            draws_individual_maps = (
                args.draw_individual_files is not None or args.draw_grid is not None
            )
            if group_flags and not (args.groups_txt is not None and draws_individual_maps):
                message = ', '.join(f"'{flag}'" for flag in group_flags)
                missing = (
                    "the databases or genomes are not grouped ('--groups-txt')"
                    if args.groups_txt is None else
                    "no maps for individual groups were asked for "
                    "('--draw-individual-files'/'--draw-grid')"
                )
                raise ConfigError(
                    f"Group-map coloring options were given ({message}), but they style the "
                    f"individual maps of each group, and here {missing}. Add what is missing, or "
                    f"drop these options: the 'unified' map summarizing every group is colored by "
                    f"'--reaction-colormap'/'--reaction-membership-colors' and "
                    f"'--presence-colormap-scheme' instead."
                )

        # To group contigs databases or pan genomes, give both '--groups-txt' and the
        # presence/absence '--group-threshold'. Give neither to leave them ungrouped.
        # '--skip-unified-maps' is the exception, as the threshold only applies to 'unified' maps;
        # '--groups-txt' may be used alone in this case. A text run's grouping depends on each
        # layer's auto-detected mode, so its group-option checks live in 'map_txt_data'.
        if not is_txt_run and not args.skip_unified_maps and (
            (args.groups_txt is None and args.group_threshold is not None) or
            (args.groups_txt is not None and args.group_threshold is None)
        ):
            raise ConfigError(
                "For groups to be used, arguments to both '--groups-txt' and '--group-threshold' "
                "must be provided."
            )

        # One database and no '--reaction-colormap' takes the route below that draws only the
        # 'unified' map: no individual maps, no grids. With a colormap it takes the route for
        # several databases instead, which draws those too.
        is_single_db_run = (
            args.contigs_dbs is not None and
            len(args.contigs_dbs) == 1 and
            args.reaction_colormap is None
        )

        # Gathering the individual maps by map is a second arrangement of files that have to exist,
        # so the flag is refused wherever the run draws none of them and it would otherwise go
        # quietly ignored. Wherever a run does draw them, it is honored however few they are: one
        # database or one sample gives a directory per map holding a single file, which is no less
        # than '--draw-individual-files' gives there on its own. A reaction-network-JSON run is left
        # to its own dispatcher. It refuses this flag along with the rest of the family, since the
        # run has no samples, databases, genomes, or groups to draw individual maps for.
        if args.collate_files_by_map and args.reaction_network_json is None:
            if args.draw_individual_files is None:
                raise ConfigError(
                    "'--collate-files-by-map' gathers the maps drawn for individual contigs "
                    "databases, genomes, samples, or groups into a directory per map, but none "
                    "were asked for. Please add '--draw-individual-files', or drop this flag."
                )
            if is_single_db_run:
                raise ConfigError(
                    "'--collate-files-by-map' gathers the maps drawn for individual contigs "
                    "databases into a directory per map, but a lone database drawn in a single "
                    "color has no individual maps to gather: the one map it draws per pathway is "
                    "already the whole of the data. Please provide the other databases to compare "
                    "it to, or drop this flag."
                )

        # Skipping the 'unified' map makes the maps of the individual databases, genomes, or groups
        # the whole of the output. The run therefore needs those maps, and has nothing to skip where
        # the 'unified' map is all that should be drawn, as in the following cases. A text run is
        # checked in 'map_txt_data', which alone knows whether any file has a 'sample' column, and a
        # reaction-network-JSON run is checked in its own dispatcher.
        if args.skip_unified_maps and not is_txt_run and args.reaction_network_json is None:
            if is_single_db_run:
                raise ConfigError(
                    "'--skip-unified-maps' should be used to skip the map summarizing several "
                    "sources at once, which does not apply to a single contigs database drawn in a "
                    "single color. Please provide other databases to compare it to, or drop this "
                    "flag."
                )
            if args.draw_individual_files is None and args.draw_grid is None:
                raise ConfigError(
                    "'--skip-unified-maps' skips the map summarizing every contigs database, "
                    "genome, or group at once. That leaves the maps of the individual databases, "
                    "genomes, or groups to be drawn. None were asked for. Please add "
                    "'--draw-individual-files' and/or '--draw-grid', or drop this flag."
                )
            if args.group_threshold is not None:
                raise ConfigError(
                    "'--group-threshold' sets the proportion of a group's databases or genomes in "
                    "which an element must occur for the group to count as containing it. Only "
                    "the 'unified' map uses that threshold. A group's own map counts that group's "
                    "own sources whatever the threshold is. '--skip-unified-maps' draws no "
                    "'unified' map. Please use '--group-threshold' or '--skip-unified-maps', not "
                    "both."
                )
            if args.presence_colormap_scheme is not None:
                raise ConfigError(
                    "'--presence-colormap-scheme' selects how the 'unified' map colors presence "
                    "across contigs databases, genomes, or groups of either. Nothing else uses it. "
                    "A group's own map is drawn by '--group-colormap'/'--group-colormap-scheme'. A "
                    "database's or genome's own map draws its reactions in a single color. "
                    "'--skip-unified-maps' draws no 'unified' map. Please use "
                    "'--presence-colormap-scheme' or '--skip-unified-maps', not both."
                )

        mapper = Mapper(kegg_dir=args.kegg_dir,
                        overwrite_output=args.overwrite_output_destinations,
                        name_files=args.name_files,
                        categorize_files=args.categorize_files,
                        collate_files_by_map=args.collate_files_by_map)
        performed = False

        if is_txt_run:
            map_txt_data(args, mapper)
            performed = True

        if args.reaction_network_json is not None:
            map_json_network_ko_data(args, mapper)
            performed = True

        if is_single_db_run:
            map_single_contigs_db_ko_data(args, mapper)
            performed = True
        elif args.contigs_dbs is not None:
            map_multiple_contigs_dbs_ko_data(args, mapper)
            performed = True

        if args.pan_db is not None and args.genomes_storage is not None:
            map_pan_db_ko_data(args, mapper)
            performed = True

        if not performed:
            raise ConfigError(
                "No task was performed! The minimum requirement is a data source: a "
                "draw-kegg-pathways reaction- and/or compound-layer text file "
                "(`--reaction-txt`/`--compound-txt`), a reaction network JSON file "
                "(`--reaction-network-json`), one or more contigs databases "
                "(`--contigs-dbs`/`--external-genomes`), or a pangenomic database (`--pan-db` with "
                "`--genomes-storage`)."
            )

    except ConfigError as e:
        e_str = re.sub(r'\s+', ' ', str(e))
        if (
            "Unprioritized entry graphics cannot be assigned the same combination of foreground "
            "and background colors as prioritized entries of the same entry and graphics types."
        ) in e_str:
            # The colors of a run are checked before anything is drawn. The check covers the colors
            # of the layers ('_check_reserved_colors') and the compound colors derived from
            # reactions ('_check_derived_compound_colors'). Only a clash those checks did not
            # foresee arrives here. The kgml error names the entry type, which says which layer's
            # colors produced it.
            if "'compound' entries" in e_str:
                layer_clause = (
                    "Specifically, a compound would be drawn in a reserved color. On a map drawn "
                    "from reactions alone, each compound takes its color from the reactions it "
                    "touches, so the reaction color or colormap is what to change "
                    "('--reaction-color' or '--reaction-colormap'). "
                )
            elif "'ortholog' entries" in e_str:
                layer_clause = (
                    "Specifically, a color given to the reaction layer collided with a reserved "
                    "color, so adjust the reaction colormap or color ('--reaction-colormap' or "
                    "'--reaction-color'). "
                )
            else:
                layer_clause = ""
            # Printed rather than raised: raising from inside this handler would reach the user as a
            # traceback instead of as the message just assembled.
            print(ConfigError(
                f"{layer_clause}The colors of highlighted reactions and compounds cannot be set to "
                f"reserved colors of other un-highlighted reactions and compounds, respectively. "
                f"In global maps, other reactions and compounds are colored gray ('#E0E0E0'), so "
                f"this should not be used as a static color or dynamic color in a colormap. In "
                f"overview maps, other reactions are colored black ('#000000') and other compounds "
                f"are colored white ('#FFFFFF'), so these should not be used as colors. In "
                f"standard maps, other reactions and compounds are colored white, so this should "
                f"not be used as a color."
            ))
            sys.exit(-1)
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)

if __name__ == '__main__':
    main()
