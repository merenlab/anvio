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
from anvio.keggmapping import (
    AGGREGATION_FUNCTIONS, DEFAULT_GROUP_TINT_SPAN, GROUP_COLORMAP_FROM_CATEGORY,
    GROUP_SCHEME_OPTIONS, MAX_DISCRETE_COUNT_BANDS, SUMMARY_PRESENCE_PHRASE,
    SUMMARY_PRESENCE_SCHEMES, Mapper
)


__authors__ = ["semiller10"]
__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = VERSION
__requires__ = ['kegg-data']
__can_use__ = ['contigs-db', 'external-genomes', 'pan-db', 'genomes-storage-db',
               'reaction-network-json', 'kegg-reaction-txt', 'kegg-compound-txt',
               'kegg-category-colors-txt']
__provides__ = ['kegg-pathway-map']
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
        "optional 'gene_id' column can label the gene a value comes from (it does not affect "
        "coloring; a value is aggregated per accession across rows). An optional 'sample' column "
        "assigns each row to a sample of origin, enabling comparison across samples (and, with "
        "'--groups-txt', groups of samples) as when using multiple contigs databases. The single "
        "column that is not 'accession', 'gene_id', or 'sample' is auto-detected as a numeric "
        "value column: with it, reactions are colored quantitatively by value; without it, by "
        "presence."
    )
    groupTXT.add_argument(
        '--compound-txt', type=str, metavar='FILE', help=
        "Path to a tab-delimited kegg-compound-txt file coloring the compound elements of maps. It "
        "must have an 'accession' column whose values are all KEGG compound IDs ('C' followed by "
        "digits). It must NOT have a 'gene_id' column, since compounds are not gene products. An "
        "optional 'sample' column works as for '--reaction-txt'. The single column that is not "
        "'accession' or 'sample' is auto-detected as a numeric value column: with it, compounds "
        "are colored quantitatively by value; without it, by presence."
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
        f"were affected."
    )
    groupTXT.add_argument(
        '--reaction-accession-aggregation', metavar='NAME', help=
        "How to reduce the values of the several accessions that ONE MAP ELEMENT stands for to "
        "that element's color, when the reaction layer is colored by a value column. Nearly a "
        "quarter of the reaction elements in KEGG's maps stand for more than one KO, so this is a "
        "common case. The default is 'sum': the KOs of a reaction contribute to that reaction's "
        "total. It takes the same names as '--reaction-gene-aggregation', and the two can differ: "
        "averaging across the genes of each KO while totaling across the KOs of a reaction, say."
    )
    groupTXT.add_argument(
        '--compound-accession-aggregation', metavar='NAME', help=
        "Like '--reaction-accession-aggregation', but for the compound layer, and with the same "
        "default of 'sum': a single circle on a map occasionally stands for more than one KEGG "
        "compound, such as the circle covering four forms of cobalamin, and this reduces the "
        "values of whichever of them the file lists to that circle's color. It applies to about "
        "one compound circle in fifty, so it rarely has anything to do. There is no reduction "
        "below it: a compound file names no genes, and rows repeating a compound are refused "
        "rather than combined, so each compound has one value per sample."
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
        "'--reaction-group-summary'/'--compound-group-summary'; whenever the groups are summarized "
        "by presence, which is the default, '--group-threshold' must also be given to define when "
        "a group contains an element."
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
        "text layer's groups are summarized by pooling values instead."
    )

    groupSUMMARY = parser.add_argument_group(
        "SAMPLE AND GROUP SUMMARIES",
        "How a draw-kegg-pathways text layer with a 'sample' column is summarized across its "
        "samples, and across groups of samples defined by '--groups-txt'. These are the two "
        "reductions above the within-sample gene and accession aggregations. Each can "
        "independently summarize presence (in how many, or exactly which, samples or groups an "
        "element occurs) or pool the values of a value column."
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
        f"'--reaction-gene-aggregation'). 'std' is an example here, mapping how much replicate "
        f"samples disagree. This summary colors the 'unified' map when the samples are not "
        f"grouped, and each group's map when they are grouped with '--groups-txt' — though only an "
        f"aggregation affects a group's map, since presence there is always the count of the "
        f"group's own samples, so the presence names are rejected in a grouped run as having "
        f"nothing to choose between. Maps for individual samples are never summarized, always "
        f"showing that one sample alone: its own values with a value column, its presence without "
        f"one. Without this option, samples are summarized by presence, by membership with 3 or "
        f"fewer samples, by count above that, and by a continuous count scale where a discrete one "
        f"would run out of distinguishable colors or of room to label its bands (above "
        f"{MAX_DISCRETE_COUNT_BANDS} of them). Presence is the default because it is meaningful "
        f"for any set of samples, whereas pooling values is only meaningful when the samples are "
        f"commensurable, such as replicates of one condition."
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
        "shows the mean of its replicate samples, while the 'unified' map shows in how many groups "
        "each element occurs."
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
        "shows the unified map of data from all databases and maps for individual databases. If "
        "using an ungrouped pangenomic database, the grid shows the pangenomic map and maps for "
        "individual genomes. If using a draw-kegg-pathways text file with a 'sample' column, the "
        "grid shows the unified map summarizing all of the samples and maps for individual "
        "samples. The grid view facilitates identification of the contigs databases, genomes, or "
        "samples containing reactions highlighted in the integrative map. If used as a flag "
        "(without values), all of the contigs databases, genomes, or samples are included in the "
        "grid. Alternatively, the project names of a subset of contigs databases, or the names of "
        "a subset of genomes or samples, can be provided. If groups are defined by '--groups-txt', "
        "then the grid instead shows the unified map alongside maps for individual contigs "
        "database, pan genome, or sample groups. A subset of group names can be provided to select "
        "maps in the grid."
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
        "This option takes the name of a Matplotlib Colormap which is sampled in coloring data. In "
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
        "'--reaction-txt' file with a value column), this same option selects the sequential "
        "colormap that is sampled continuously to map reaction values (default 'plasma_r'), and "
        "'--presence-colormap-scheme' does not apply. That one colormap then colors the 'unified' "
        "map and the maps of the individual samples or groups alike, unless "
        "'--reaction-category-colormap' gives the latter a colormap of their own. See the "
        "following webpage for named colormaps: "
        "https://matplotlib.org/stable/users/explain/colors/colormaps.html#classes-of-colormaps"
    )
    groupCOLOR.add_argument(
        '--reaction-category-colormap', nargs='+', help=
        "Like '--reaction-colormap', but for the single color scale shared by the maps of the "
        "individual samples ('--draw-individual-files'/'--draw-grid'), or of the individual groups "
        "when the samples are grouped with '--groups-txt'. It takes the same values: a Matplotlib "
        "colormap name on its own, or a name followed by two decimals limiting the fraction of the "
        "colormap to sample. Without it, one colormap colors both of those maps and the 'unified' "
        "map, which is what you want when the two show the same quantity — the values of one "
        "sample beside the mean of them all, for example. Give it when they do not. With "
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
        "colors instead."
    )
    groupCOLOR.add_argument(
        '--compound-category-colormap', nargs='+', help=
        "Like '--reaction-category-colormap', but for the compound layer of a '--compound-txt' "
        "file colored by a value column. The two layers are colored on independent scales, each "
        "with its own colorbar, so each takes its own colormaps."
    )
    groupCOLOR.add_argument(
        '--reaction-category-colors', type=str, metavar='PATH', help=
        "Path to a tab-delimited kegg-category-colors-txt file giving a color hex code to each "
        "category of the run — each sample of a '--reaction-txt' file, each contigs database, each "
        "genome of a pangenome, or each group of any of these when '--groups-txt' is used. This "
        "colors the reaction layer by membership from colors you choose rather than from a "
        "colormap, so it replaces '--reaction-colormap' and cannot be combined with it. Its first "
        "column holds category names and can be headed anything; its 'color' column holds the hex "
        "codes. Because coloring by membership colors an element by exactly WHICH categories "
        "contain it, a color is needed for every combination of them: a category on its own takes "
        "its own color, and a combination takes the blend of its members' colors, which a row "
        "naming the combination — its category names separated by commas, exactly as the colorbar "
        "labels it — overrides. Where the categories are samples, databases or genomes, the same "
        "colors also color each category's own map ('--draw-individual-files'/'--draw-grid'), so a "
        "panel of a grid says which category it is and matches the band it takes on the 'unified' "
        "map; where they are groups, each group's own maps count that group's sources instead, and "
        "'--group-colormap category' is what carries the group's color onto them. Every category "
        "of the run needs a color; names in the file that are not categories of the run are "
        "reported and ignored, so one file can serve runs drawing different subsets. Note that "
        "giving this file colors the layer by membership however many categories there are, rather "
        "than by count above the three that would otherwise switch to counting: six categories "
        "then take a 63-band scale of every combination of them rather than a six-band scale of "
        "counts. Since no color scale can distinguish combinations of more than a handful of "
        "categories, anvi'o refuses rather than draw a scale whose bands cannot be told apart."
    )
    groupCOLOR.add_argument(
        '--compound-category-colors', type=str, metavar='PATH', help=
        "Like '--reaction-category-colors', but for the compound layer of a '--compound-txt' file "
        "colored by sample or group presence. The two layers are colored independently, so each "
        "takes its own file; pointing both options at a single file gives the two layers one set "
        "of category colors and so one legend. It cannot be combined with '--compound-colormap'."
    )
    groupCOLOR.add_argument(
        '--presence-colormap-scheme',
        choices=['by_count', 'by_count_continuous', 'by_membership'], help=
        f"There are three ways of dynamically coloring elements by occurrence in multiple contigs "
        f"databases ('--contigs-dbs'), the genomes of a pangenome, or groups ('--groups-txt') of "
        f"either: by count in discrete bands, by count on a continuous scale, or by membership. "
        f"For a draw-kegg-pathways text file this choice is part of the sample and group summaries "
        f"instead ('--reaction-sample-summary' and the related options), so this option is "
        f"rejected for a text run. By default, with 4 or more databases or groups, reactions are "
        f"colored by count of database or group; with 2 or 3, reactions are colored explicitly by "
        f"database or group membership. In coloring by count, the colormap should be sequential, "
        f"such that the color of a reaction changes 'smoothly' with the count. Because its "
        f"colorbar labels one band per count, 'by_count' additionally needs a distinguishable "
        f"color for every count, and is refused when there are more databases or groups than the "
        f"colormap can supply — with a few hundred of them no color scale could tell them apart "
        f"anyway — and it needs room to set every one of those labels, so above about "
        f"{MAX_DISCRETE_COUNT_BANDS} counts it warns that they will be too small to read. "
        f"'by_count_continuous' colors counts from the same colormap in the same way but draws its "
        f"colorbar as a gradient from the lowest count to the highest, which needs neither a "
        f"distinguishable color nor a label per count, so a color there reads as a position along "
        f"the range rather than as an exact count; it is chosen automatically, with a warning, "
        f"wherever 'by_count' would be refused or its labels would be unreadable. In contrast, "
        f"coloring by membership means reaction color is determined by membership in a "
        f"database/group or combination of databases/groups, so a qualitative colormap can be used "
        f"instead of a sequential colormap, as by default with 2 or 3 categories, to give a "
        f"distinct color to each membership category."
    )
    groupCOLOR.add_argument(
        '--count-scale-max', metavar='NAME_OR_NUMBER', default='observed', help=
        "Where a color scale of counts stops, which decides how the colors are spread over the "
        "counts. By default, 'observed' stops at the highest count anything on the drawn maps "
        "actually has, so that the colors span the counts that occur. 'total' runs the scale to "
        "every category there is — every sample, contigs database, genome, or group — which is "
        "what the count could be at most. A number stops the scale there instead, and counts above "
        "it take the top color, which is how separate figures are given the same scale. The "
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
        "'--compound-reverse-overlay'."
    )
    groupCOLOR.add_argument(
        '--compound-colormap', nargs='+', help=
        "Like '--reaction-colormap', but for the compound layer when drawing compounds from a "
        "'--compound-txt' file colored by a value column. The reaction and compound layers are "
        "colored on independent scales, each with its own colorbar (both default 'plasma_r'). As "
        "with '--reaction-colormap', this takes a Matplotlib colormap name, optionally followed by "
        "two decimal limits."
    )
    groupCOLOR.add_argument(
        '--compound-reverse-overlay', action='store_true', default=False, help=
        "Like '--reaction-reverse-overlay', but for the compound layer of a '--compound-txt' file. "
        "By default, compounds with a greater value are drawn on top of those with a lesser value; "
        "with this flag, the opposite applies. Note that compound circles rarely overlap on a map, "
        "unlike the crossing reaction arrows of global and overview maps, so the compound drawing "
        "order is usually of little visual consequence."
    )
    groupCOLOR.add_argument(
        '--reaction-value-limits', nargs=2, metavar='LIMIT', help=
        "Limits on the values that the reaction layer's color scale spans on the 'unified' map, "
        "given as a minimum and then a maximum, either of which can be the word 'none' to leave "
        "that end where the data puts it. For instance, '--reaction-value-limits -6 none' stops "
        "the scale at -6 from below and lets it run up to whatever the values reach. A limit only "
        "takes effect where the values actually cross it: it does nothing to a scale whose values "
        "all sit inside it. Where a limit does truncate, every element past it takes the color of "
        "that end of the scale, and the colorbar labels on that end are labeled with '<=' or '>=' "
        "so a reader can tell that its color stands for that value or anything beyond it rather "
        "than for the value alone. The limits are read in the units of the colorbar, which are the "
        "values of map elements once both reductions have been applied "
        "('--reaction-gene-aggregation' and then '--reaction-accession-aggregation'), so with the "
        "default of 'sum' at the accession level a single element standing for a dozen KOs may sit "
        "well past any one value in the file. Note that these limits are an entirely different "
        "thing from the two decimals '--reaction-colormap' accepts, which choose what fraction of "
        "the colormap to sample and say nothing at all about the values."
    )
    groupCOLOR.add_argument(
        '--reaction-category-value-limits', nargs=2, metavar='LIMIT', help=
        "Like '--reaction-value-limits', but for the single color scale shared by the maps of the "
        "individual samples ('--draw-individual-files'/'--draw-grid'), or of the individual groups "
        "when the samples are grouped with '--groups-txt'. The two scales are limited separately "
        "because they can span quite different things: where '--reaction-sample-summary mean' "
        "colors the 'unified' map, that map shows means across samples, which cluster far more "
        "tightly than the per-sample values underlying them, so limits that suit one can be badly "
        "wrong for the other. Giving one of the two options without the other is perfectly "
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
        '--compound-category-value-limits', nargs=2, metavar='LIMIT', help=
        "Like '--reaction-category-value-limits', but for the compound layer: the scale shared by "
        "the compound colors of the individual sample or group maps."
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
        f"colormap goes partly unused."
    )
    groupCOLOR.add_argument(
        '--reaction-category-value-center', nargs='?', const=VALUE_CENTER_DEFAULT, metavar='VALUE',
        help=
        "Like '--reaction-value-center', but for the single color scale shared by the maps of the "
        "individual samples ('--draw-individual-files'/'--draw-grid'), or of the individual groups "
        "when the samples are grouped with '--groups-txt'. The two scales are centered separately, "
        "for the same reason they are limited separately: a summary can put them on quite "
        "different footings. Giving one of the two without the other is legitimate, and anvi'o "
        "warns when you do, since one colormap colors both scales unless "
        "'--reaction-category-colormap' gives them one each, and the middle color would then mean "
        "the centered value on the one map and whatever the values leave in the middle on the "
        "other."
    )
    groupCOLOR.add_argument(
        '--compound-value-center', nargs='?', const=VALUE_CENTER_DEFAULT, metavar='VALUE', help=
        "Like '--reaction-value-center', but for the compound layer of a '--compound-txt' file "
        "colored by a value column. The two layers are colored on independent scales, each with "
        "its own colorbar, so each takes its own center."
    )
    groupCOLOR.add_argument(
        '--compound-category-value-center', nargs='?', const=VALUE_CENTER_DEFAULT, metavar='VALUE',
        help=
        "Like '--reaction-category-value-center', but for the compound layer: the scale shared by "
        "the compound colors of the individual sample or group maps."
    )
    groupGROUP.add_argument(
        '--group-colormap', nargs='+', help=
        f"This option is like '--reaction-colormap', but only applies to drawing files for "
        f"individual groups ('--draw-individual-files') and panels for individual groups in map "
        f"grids ('--draw-grid'). These maps for individual groups show data from group sources — "
        f"samples, contigs databases, or pan genomes. Presence may only be colored by count — "
        f"e.g., the number of samples in the group containing the data — not membership. Like "
        f"'--reaction-colormap', this parameter takes the name of a Matplotlib Colormap, and "
        f"optionally, two decimal values between 0.0 and 1.0 to limit the fraction of the colormap "
        f"used. The default configuration is the same, with the colormap being 'plasma_r' and the "
        f"limits being 0.1 and 0.9. Note that per-group maps show counts in discrete bands by "
        f"default, one per count in the group's scale, so a group whose scale spans more counts "
        f"than this colormap has distinguishable colors or over more than "
        f"{MAX_DISCRETE_COUNT_BANDS} of them has its count drawn on a continuous scale instead, "
        f"along with a warning saying so ('--group-colormap-scheme'). This option colors group "
        f"maps that show counts. Where a layer's samples are pooled into a value instead — "
        f"'--reaction-sample-summary'/'--compound-sample-summary' given an aggregation — its group "
        f"maps show that value, and are colored by "
        f"'--reaction-category-colormap'/'--compound-category-colormap', or by the layer's own "
        f"'--reaction-colormap'/'--compound-colormap' if that is not given. In place of a colormap "
        f"name, the value '{GROUP_COLORMAP_FROM_CATEGORY}' colors each group's own maps by a ramp "
        f"running from a pale tint to that group's own color, which "
        f"'--reaction-category-colors'/'--compound-category-colors' gives; every group's ramp is "
        f"built the same way, so the panels of a grid stay comparable while each says which group "
        f"it is. The two decimal values then say how far from white the ramp starts and stops "
        f"rather than which fraction of a colormap to use, defaulting to "
        f"{DEFAULT_GROUP_TINT_SPAN[0]} and {DEFAULT_GROUP_TINT_SPAN[1]}: the pale end stops short "
        f"of white because standard and overview maps keep white for their own unhighlighted "
        f"elements. A named colormap remains the default, since one built for showing magnitude "
        f"usually reads better as a count than a color chosen to identify a group does."
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


def check_package_dependencies():
    missing_packages = []
    for package_name in ('Bio', 'reportlab', 'fitz'):
        try:
            importlib.import_module(package_name)
        except ImportError:
            missing_packages.append(package_name)

    if missing_packages:
        message = ', '.join(f"'{package_name}'" for package_name in missing_packages)
        raise ConfigError(
            f"The following Python packages required to run `anvi-draw-kegg-pathways` could not be "
            f"imported: {message}. All Python dependencies of anvi'o can be installed by running "
            f"the command `pip install -r requirements.txt` in the top directory of the anvi'o "
            f"codebase, which may have a location like '~/github/anvio' in your file system if you "
            f"followed the anvio.org installation instructions."
        )

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
    if args.reaction_category_colormap is not None:
        unsupported_args.append('--reaction-category-colormap')
    if args.compound_category_colormap is not None:
        unsupported_args.append('--compound-category-colormap')
    if args.reaction_category_colors is not None:
        unsupported_args.append('--reaction-category-colors')
    if args.compound_category_colors is not None:
        unsupported_args.append('--compound-category-colors')
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

def _detect_txt_layer(path: str) -> dict:
    """Light structural peek at a per-layer text file for CLI flag validation.

    Returns {'mode': 'quantitative'|'membership'|'single'|None, 'has_sample': bool}. 'mode' is None
    when the file is malformed in a way the reader ('Mapper._read_element_txt') should report (no
    'accession' column, or an ambiguous >1 value column); the caller then skips flag guards and lets
    the reader raise the definitive error.
    """
    filesnpaths.is_file_exists(path)
    try:
        columns = list(pd.read_csv(path, sep='\t', nrows=0).columns)
    except pd.errors.EmptyDataError:
        # An empty file has no columns to detect, so defer to the reader for the definitive error.
        return {'mode': None, 'has_sample': False, 'has_gene_id': False}
    has_sample = 'sample' in columns
    has_gene_id = 'gene_id' in columns
    if 'accession' not in columns:
        return {'mode': None, 'has_sample': has_sample, 'has_gene_id': has_gene_id}
    value_columns = [c for c in columns if c not in ('accession', 'gene_id', 'sample')]
    if len(value_columns) > 1:
        return {'mode': None, 'has_sample': has_sample, 'has_gene_id': has_gene_id}
    if len(value_columns) == 1:
        mode = 'quantitative'
    elif has_sample:
        mode = 'membership'
    else:
        mode = 'single'
    return {'mode': mode, 'has_sample': has_sample, 'has_gene_id': has_gene_id}


def map_txt_data(args: Namespace, mapper: Mapper) -> None:
    """Draw reaction and/or compound data from per-layer draw-kegg-pathways text files."""
    def _resolve_colormap(colormap_arg, flag):
        # '--reaction-colormap'/'--compound-colormap'/'--group-colormap' are 'nargs' args: a name,
        # or a name plus two fraction limits. None means use the layer default.
        if colormap_arg is None:
            return None, None
        if len(colormap_arg) == 1:
            return colormap_arg[0], None
        # Given '--group-colormap category' the two numbers are not a fraction of a colormap at all,
        # but how far from white a group's own ramp starts and stops, so the messages say what the
        # numbers actually mean. The constraint on them is the same either way.
        from_category = (
            flag == '--group-colormap' and colormap_arg[0] == GROUP_COLORMAP_FROM_CATEGORY
        )
        limits_mean = (
            "how far from white towards a group's own color its ramp starts and stops"
            if from_category else "the fraction of the colormap to use"
        )
        if len(colormap_arg) != 3:
            example = (
                f"'{GROUP_COLORMAP_FROM_CATEGORY} {DEFAULT_GROUP_TINT_SPAN[0]} "
                f"{DEFAULT_GROUP_TINT_SPAN[1]}'"
            ) if from_category else "'plasma 0.1 0.9'"
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
    # the reader). A layer is 'quantitative' (a value column), 'membership' (presence + a 'sample'
    # column), or 'single' (presence, no samples). This says what the file provides, from which the
    # mapper derives how the layer colors each map context.
    reaction = _detect_txt_layer(args.reaction_txt) if args.reaction_txt is not None else None
    compound = _detect_txt_layer(args.compound_txt) if args.compound_txt is not None else None
    present = [layer for layer in (reaction, compound) if layer is not None]
    any_sample = any(layer['has_sample'] for layer in present)

    draw_individual_files = _nargs_selection(args.draw_individual_files)
    draw_grid = _nargs_selection(args.draw_grid)
    draw_category_maps = draw_individual_files is not False or draw_grid is not False

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
            if reaction['mode'] == 'quantitative':
                raise ConfigError(
                    "'--original-color' colors reactions in the reference map's own colors, which "
                    "leaves no color scale for a value column, but the reaction file has one. Draw "
                    "the values with '--reaction-colormap' instead, or drop the value column to "
                    "highlight the reactions in the reference colors."
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
                (args.group_threshold is not None, '--group-threshold')
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

        # A color per category colors a layer by which samples or groups contain an element, and
        # colors each of their own maps, so the layer's file needs a 'sample' column for there to be
        # any categories at all. Which map contexts actually take those colors depends on the
        # summaries, which the mapper resolves; it refuses there if none of them does.
        for element_type, file_flag, layer, _, _ in summaries:
            if getattr(args, f'{element_type}_category_colors') is None:
                continue
            colors_flag = f'--{element_type}-category-colors'
            if layer is None:
                raise ConfigError(
                    f"'{colors_flag}' colors the {element_type} layer, but no {element_type} file "
                    f"('{file_flag}') was provided."
                )
            if not layer['has_sample']:
                raise ConfigError(
                    f"'{colors_flag}' gives a color to each sample, or to each group of samples, "
                    f"but the {element_type} file has no 'sample' column, so there are no "
                    f"categories to color. Set the single color of this layer with "
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
            if layer['mode'] != 'quantitative':
                raise ConfigError(
                    f"'{flag}' reduces the values of the {element_type} layer, but its file has no "
                    f"value column, so the layer is colored by presence and there are no values to "
                    f"reduce."
                )
            if flag == '--reaction-gene-aggregation' and not layer['has_gene_id']:
                raise ConfigError(
                    f"'{flag}' reduces the values of the genes annotated with one accession, but "
                    f"the reaction file has no 'gene_id' column naming genes. Each accession there "
                    f"already has a single value per sample, since rows repeating an accession are "
                    f"refused. Use '--reaction-accession-aggregation' to set how a map element "
                    f"combines the several accessions it stands for."
                )

        # Each pair of value limits bounds, and each value center centers, a color scale of values,
        # so each needs its layer's file and that file's value column; the ones acting on the scale
        # the individual maps share need the samples those maps are drawn for as well. The mapper
        # checks everything past this point like a programmatic call, including the numbers
        # themselves and whether the context in question ends up colored by value once the summaries
        # have had their say.
        for setting, flag, element_type, file_flag, layer, per_category, verb in (
            (args.reaction_value_limits, '--reaction-value-limits', 'reaction', '--reaction-txt',
             reaction, False, 'bound'),
            (args.reaction_category_value_limits, '--reaction-category-value-limits', 'reaction',
             '--reaction-txt', reaction, True, 'bound'),
            (args.compound_value_limits, '--compound-value-limits', 'compound', '--compound-txt',
             compound, False, 'bound'),
            (args.compound_category_value_limits, '--compound-category-value-limits', 'compound',
             '--compound-txt', compound, True, 'bound'),
            (args.reaction_value_center, '--reaction-value-center', 'reaction', '--reaction-txt',
             reaction, False, 'center'),
            (args.reaction_category_value_center, '--reaction-category-value-center', 'reaction',
             '--reaction-txt', reaction, True, 'center'),
            (args.compound_value_center, '--compound-value-center', 'compound', '--compound-txt',
             compound, False, 'center'),
            (args.compound_category_value_center, '--compound-category-value-center', 'compound',
             '--compound-txt', compound, True, 'center')
        ):
            if setting is None:
                continue
            if layer is None:
                raise ConfigError(
                    f"'{flag}' {verb}s the color scale of the {element_type} layer, but no "
                    f"{element_type} file ('{file_flag}') was provided."
                )
            if layer['mode'] != 'quantitative':
                raise ConfigError(
                    f"'{flag}' {verb}s the color scale spanning the values of the {element_type} "
                    f"layer, but its file has no value column, so the layer is colored by presence "
                    f"and there is no scale of values to {verb}."
                )
            if per_category and not layer['has_sample']:
                raise ConfigError(
                    f"'{flag}' {verb}s the color scale shared by the maps of the individual "
                    f"samples or groups, but the {element_type} file has no 'sample' column, so it "
                    f"draws no such maps and its values take a single scale. {verb.capitalize()} "
                    f"that one with "
                    f"'--{element_type}-value-{'limits' if verb == 'bound' else 'center'}'."
                )

        # A colormap for the scale the individual maps share needs the same things that pair of
        # limits needs: the layer's file, a value column to put those maps on a scale, and the
        # samples they are drawn for. Whether that scale ends up colored by value once the summaries
        # have had their say is checked by the mapper, as for programmatic callers.
        for category_colormap, element_type, file_flag, layer in (
            (args.reaction_category_colormap, 'reaction', '--reaction-txt', reaction),
            (args.compound_category_colormap, 'compound', '--compound-txt', compound)
        ):
            if category_colormap is None:
                continue
            flag = f'--{element_type}-category-colormap'
            if layer is None:
                raise ConfigError(
                    f"'{flag}' colors the maps of the individual samples or groups of the "
                    f"{element_type} layer, but no {element_type} file ('{file_flag}') was "
                    f"provided."
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
        # one statement per accession, so it needs its layer's file, a 'sample' column to summarize,
        # and, when it pools values, a value column to pool.
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
                if summary not in PRESENCE_SUMMARIES and layer['mode'] != 'quantitative':
                    raise ConfigError(
                        f"'{flag}' was given as '{summary}', which is an aggregation that pools "
                        f"values, but the {element_type} file has no value column, so there are no "
                        f"values to pool. Summarize presence with one of {SUMMARY_PRESENCE_PHRASE} "
                        f"instead, or add a value column to the file."
                    )
                if level == 'group' and args.groups_txt is None:
                    raise ConfigError(
                        f"'{flag}' summarizes groups of samples, but no groups were defined with "
                        f"'--groups-txt'."
                    )
                # With groups, presence at the SAMPLE level is always the count of a group's own
                # samples in discrete bands from '--group-colormap', so every presence name
                # describes the same per-group map and none of them changes anything: the choice
                # only shapes the ungrouped 'unified' map.
                if (
                    level == 'sample' and args.groups_txt is not None
                    and summary in PRESENCE_SUMMARIES
                ):
                    raise ConfigError(
                        f"'{flag}' was given as '{summary}', but the samples here are grouped with "
                        f"'--groups-txt', and each group's map shows the count of that group's own "
                        f"samples regardless of how presence is asked for — in discrete bands or on "
                        f"a continuous scale, as '--group-colormap-scheme' says — so "
                        f"{SUMMARY_PRESENCE_PHRASE} would all draw the same thing. Summarizing the "
                        f"samples by presence is already the default, so drop the option; use an "
                        f"aggregation such as 'mean' to color each group's map by pooled values "
                        f"instead, and note that '--{element_type}-group-summary' is what chooses "
                        f"how presence is colored across the groups themselves."
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
            if group_summary is not None and group_summary not in PRESENCE_SUMMARIES and (
                sample_summary is None or sample_summary in PRESENCE_SUMMARIES
            ):
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
                "layer's '--*-group-summary' pools values instead."
            )
        if group_presence and args.group_threshold is None:
            raise ConfigError(
                "When grouping samples ('--groups-txt') whose groups are summarized by presence, "
                "which is the default, '--group-threshold' must also be given: the proportion of a "
                "group's samples in which an element must be present for the group to be "
                "considered to contain it. Alternatively, summarize the groups by pooling their "
                "values with '--reaction-group-summary'/'--compound-group-summary'."
            )
        group_flags = [flag for present_flag, flag in (
            (args.group_colormap is not None, '--group-colormap'),
            (args.group_colormap_scheme is not None, '--group-colormap-scheme'),
            (args.group_reverse_overlay, '--group-reverse-overlay')
        ) if present_flag]
        if group_flags and not (group_presence_maps and draw_category_maps):
            message = ', '.join(f"'{flag}'" for flag in group_flags)
            raise ConfigError(
                f"Group-map coloring options were given ({message}), but they apply only to "
                f"individual group maps showing within-group sample counts, which are drawn with "
                f"'--groups-txt' plus '--draw-individual-files' and/or '--draw-grid' for a layer "
                f"whose samples are summarized by presence. A layer whose samples are summarized "
                f"by pooling values colors its group maps continuously from "
                f"'--reaction-category-colormap'/'--compound-category-colormap' (or the layer's "
                f"own '--reaction-colormap'/'--compound-colormap') and "
                f"'--reaction-reverse-overlay'/'--compound-reverse-overlay' instead."
            )

    reaction_colormap, reaction_colormap_limits = _resolve_colormap(
        args.reaction_colormap, '--reaction-colormap'
    )
    compound_colormap, compound_colormap_limits = _resolve_colormap(
        args.compound_colormap, '--compound-colormap'
    )
    reaction_category_colormap, reaction_category_colormap_limits = _resolve_colormap(
        args.reaction_category_colormap, '--reaction-category-colormap'
    )
    compound_category_colormap, compound_category_colormap_limits = _resolve_colormap(
        args.compound_category_colormap, '--compound-category-colormap'
    )
    group_colormap, group_colormap_limits = _resolve_colormap(
        args.group_colormap, '--group-colormap'
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
        'draw_individual_files': draw_individual_files,
        'draw_grid': draw_grid,
        'reaction_reverse_overlay': args.reaction_reverse_overlay,
        'compound_reverse_overlay': args.compound_reverse_overlay,
        'reaction_value_limits': _resolve_value_limits(args.reaction_value_limits),
        'reaction_category_value_limits': _resolve_value_limits(
            args.reaction_category_value_limits
        ),
        'compound_value_limits': _resolve_value_limits(args.compound_value_limits),
        'compound_category_value_limits': _resolve_value_limits(
            args.compound_category_value_limits
        ),
        # The centers go to the mapper as they were typed, bare flag included
        # ('VALUE_CENTER_DEFAULT' standing in for the value it did not carry), so that the number is
        # converted and checked in the one place programmatic callers reach too.
        'reaction_value_center': args.reaction_value_center,
        'reaction_category_value_center': args.reaction_category_value_center,
        'compound_value_center': args.compound_value_center,
        'compound_category_value_center': args.compound_category_value_center,
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
    if reaction_category_colormap is not None:
        kwargs['reaction_category_colormap'] = reaction_category_colormap
        kwargs['reaction_category_colormap_limits'] = reaction_category_colormap_limits
    if compound_category_colormap is not None:
        kwargs['compound_category_colormap'] = compound_category_colormap
        kwargs['compound_category_colormap_limits'] = compound_category_colormap_limits
    if group_colormap is not None:
        kwargs['group_colormap'] = group_colormap
        # Limits left unset are resolved by the mapper, which defaults them differently for a named
        # colormap (the fraction of it to use) and for a ramp built from a group's own color (how
        # far from white it runs).
        kwargs['group_colormap_limits'] = group_colormap_limits
    for category_colors, key in (
        (args.reaction_category_colors, 'reaction_category_colors'),
        (args.compound_category_colors, 'compound_category_colors')
    ):
        if category_colors is not None:
            kwargs[key] = category_colors
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

    if args.reaction_category_colors is not None:
        # Color reactions by membership from a color given per database, genome, or group, rather
        # than from a colormap.
        map_contigs_databases_kos = functools.partial(
            map_contigs_databases_kos, reaction_category_colors=args.reaction_category_colors
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
        # Dynamically color reactions by genome or group in unified maps using the default colormap.
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

    if args.reaction_category_colors is not None:
        # Color reactions by membership from a color given per database, genome, or group, rather
        # than from a colormap.
        map_pan_database_kos = functools.partial(
            map_pan_database_kos, reaction_category_colors=args.reaction_category_colors
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
        reaction_reverse_overlay=args.reaction_reverse_overlay,
        group_reverse_overlay=args.group_reverse_overlay,
        group_colormap_scheme=args.group_colormap_scheme,
        count_scale_max=args.count_scale_max,
        draw_maps_lacking_data=args.draw_bare_maps
    )

def main() -> None:
    args = get_args()

    try:
        check_package_dependencies()
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
                    f"Use 'observed' to stop a count scale at the highest count in the data, "
                    f"'total' to run it to every category there is, or a number to stop it there."
                )
            if args.count_scale_max < 1:
                raise ConfigError(
                    f"'--count-scale-max' was given as {args.count_scale_max}, but a scale of "
                    f"counts has to reach at least 1, since an element colored by count is in at "
                    f"least one category."
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

        # A color per category and a colormap are two answers to the same question — where a layer's
        # colors come from — so a layer takes one or the other. A static color, or the reference
        # map's own colors, answers it a third way.
        for element_type in ('reaction', 'compound'):
            if getattr(args, f'{element_type}_category_colors') is None:
                continue
            colors_flag = f'--{element_type}-category-colors'
            if getattr(args, f'{element_type}_colormap') is not None:
                raise ConfigError(
                    f"'{colors_flag}' gives the {element_type} layer a color per category, and "
                    f"'--{element_type}-colormap' has anvi'o choose its colors from a colormap "
                    f"instead. Please use only one."
                )
            if getattr(args, f'{element_type}_color') is not None:
                raise ConfigError(
                    f"'{colors_flag}' gives the {element_type} layer a color per category, while "
                    f"'--{element_type}-color' colors it one static color for presence in any "
                    f"category. Please use only one."
                )
            if args.original_color:
                raise ConfigError(
                    f"'{colors_flag}' gives the {element_type} layer a color per category, but "
                    f"'--original-color' takes both its colors and its drawing order from the "
                    f"reference map, so there is no color left to choose. Please use only one."
                )

        # One contigs database is one category, so its membership has a single combination and a
        # single color: that is a static color rather than a scale, and the run takes the
        # single-database path, which has no category dimension for a colors file to describe.
        if args.reaction_category_colors is not None and args.contigs_dbs is not None and (
            len(args.contigs_dbs) == 1
        ):
            raise ConfigError(
                "'--reaction-category-colors' gives a color to each contigs database or group of "
                "them, which colors reactions by exactly which of them contain a reaction. Only "
                "one contigs database was provided, so every reaction it contains is in the same "
                "single database and there is nothing to tell apart by color. Please provide more "
                "databases, or set the color of this one with '--reaction-color'."
            )

        # A ramp per group runs to that group's own color, so the colors have to come from
        # somewhere. The layer-level checks in 'map_txt_data' and below decide whether the file
        # given is one this run can use; this only catches asking for the ramp without any file.
        if args.group_colormap is not None and args.group_colormap[0] == (
            GROUP_COLORMAP_FROM_CATEGORY
        ) and args.reaction_category_colors is None and args.compound_category_colors is None:
            # A database, pangenome or JSON run has a reaction layer alone, so pointing it at
            # '--compound-category-colors' would only earn it a second error.
            colors_flags = (
                "'--reaction-category-colors'/'--compound-category-colors'" if is_txt_run
                else "'--reaction-category-colors'"
            )
            raise ConfigError(
                f"'--group-colormap {GROUP_COLORMAP_FROM_CATEGORY}' colors each group's own maps "
                f"by a ramp running from a pale tint to that group's own color, but no color was "
                f"given for any group. Give one per group with {colors_flags}, or name a "
                f"Matplotlib colormap for the group maps instead."
            )
        # Only a compound text file provides a compound layer, so every option that shapes one has
        # nothing to act on without it. The same goes for the sample and group summaries and the
        # value limits, which only the text engine reads: a value column is what puts elements on a
        # continuous scale, and only a text file has one.
        if not is_txt_run:
            compound_only_flags = [flag for present_flag, flag in (
                (args.compound_color is not None, '--compound-color'),
                (args.compound_colormap is not None, '--compound-colormap'),
                (args.compound_category_colors is not None, '--compound-category-colors'),
                (args.compound_reverse_overlay, '--compound-reverse-overlay'),
                (args.compound_accession_aggregation is not None,
                 '--compound-accession-aggregation'),
                (args.reaction_gene_aggregation is not None, '--reaction-gene-aggregation'),
                (args.reaction_accession_aggregation is not None,
                 '--reaction-accession-aggregation'),
                (args.reaction_sample_summary is not None, '--reaction-sample-summary'),
                (args.compound_sample_summary is not None, '--compound-sample-summary'),
                (args.reaction_group_summary is not None, '--reaction-group-summary'),
                (args.compound_group_summary is not None, '--compound-group-summary'),
                (args.reaction_value_limits is not None, '--reaction-value-limits'),
                (args.reaction_category_value_limits is not None,
                 '--reaction-category-value-limits'),
                (args.compound_value_limits is not None, '--compound-value-limits'),
                (args.compound_category_value_limits is not None,
                 '--compound-category-value-limits'),
                (args.reaction_value_center is not None, '--reaction-value-center'),
                (args.reaction_category_value_center is not None,
                 '--reaction-category-value-center'),
                (args.compound_value_center is not None, '--compound-value-center'),
                (args.compound_category_value_center is not None,
                 '--compound-category-value-center'),
                (args.reaction_category_colormap is not None, '--reaction-category-colormap'),
                (args.compound_category_colormap is not None, '--compound-category-colormap')
            ) if present_flag]
            if compound_only_flags:
                message = ', '.join(f"'{flag}'" for flag in compound_only_flags)
                raise ConfigError(
                    f"These options were given: {message}. They color a compound layer, or "
                    f"summarize the samples of a draw-kegg-pathways text file, or bound, center or "
                    f"color a scale of values, all of which only '--reaction-txt'/'--compound-txt' "
                    f"provide. Database, pangenome and reaction-network-JSON inputs have a "
                    f"reaction layer alone, drawn from one source per database or genome and "
                    f"colored by how many of them contain an element rather than by a value, so "
                    f"none of these can apply."
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
        # only for a grouped run that asks for them. Rather than accepting them and doing nothing,
        # say so — as a text run does from 'map_txt_data', and a reaction-network-JSON run from its
        # own dispatcher, which refuses the whole family outright for want of any categories at all.
        # Note that a static reaction color does not make them idle: a group's own map still counts
        # the group's databases or genomes, and still takes a colorbar of its own.
        if args.contigs_dbs is not None or args.pan_db is not None:
            group_flags = [flag for present_flag, flag in (
                (args.group_colormap is not None, '--group-colormap'),
                (args.group_colormap_scheme is not None, '--group-colormap-scheme'),
                (args.group_reverse_overlay, '--group-reverse-overlay')
            ) if present_flag]
            draws_category_maps = (
                args.draw_individual_files is not None or args.draw_grid is not None
            )
            if group_flags and not (args.groups_txt is not None and draws_category_maps):
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
                    f"'--reaction-colormap'/'--reaction-category-colors' and "
                    f"'--presence-colormap-scheme' instead."
                )

        # Categorical grouping of contigs databases or pan genomes pairs '--groups-txt' with the
        # presence/absence '--group-threshold' (both or neither). A text run's grouping depends on
        # each layer's auto-detected mode, so its group-option checks live in 'map_txt_data'.
        if not is_txt_run and (
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
        # to its own dispatcher, which refuses this flag along with the rest of the family, for want
        # of any categories at all to draw individually.
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
            # A layer's own colors are checked before anything is drawn ('_check_reserved_colors'),
            # so reaching here means the clash came from a color anvi'o derived rather than one that
            # was asked for: on a reaction-only global or overview map, compounds take a color
            # averaged from the reactions around them. The kgml error names the entry type, which
            # says which layer's colormap produced it.
            if "'compound' entries" in e_str:
                layer_clause = (
                    "Specifically, a color given to the compound layer collided with a reserved "
                    "color. On a map drawn from reactions alone, each compound takes a color "
                    "derived from the reactions it touches, so the reaction colormap "
                    "('--reaction-colormap') is what to adjust; an explicit compound layer takes "
                    "its colors from '--compound-colormap' or '--compound-color'. "
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
