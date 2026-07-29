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
from anvio.keggmapping import AGGREGATION_FUNCTIONS, SUMMARY_DISCRETE_SCHEMES, Mapper


__authors__ = ["semiller10"]
__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = VERSION
__requires__ = ['kegg-data']
__can_use__ = ['contigs-db', 'external-genomes', 'pan-db', 'genomes-storage-db',
               'reaction-network-json', 'kegg-reaction-txt', 'kegg-compound-txt']
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
PRESENCE_SUMMARIES = tuple(SUMMARY_DISCRETE_SCHEMES)

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
        "'--reaction-group-summary'/'--compound-group-summary' is 'count' or 'membership' (the "
        "default). It does not apply when every text layer's groups are summarized by pooling "
        "values instead."
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
        f"How to summarize a SET OF SAMPLES for the reaction layer. 'count' and 'membership' "
        f"summarize presence, coloring an element by how many or exactly which samples contain it. "
        f"Any other name is an aggregation that pools the samples' values from the layer's value "
        f"column, which the file must have: the recommended ones are {RECOMMENDED_AGGREGATIONS}, "
        f"and any other pandas aggregation reducing numbers to a single number works too (see "
        f"'--reaction-gene-aggregation'). 'std' is a natural fit here, mapping how much replicate "
        f"samples disagree. This summary colors the 'unified' map when the samples are not "
        f"grouped, and each group's map when they are grouped with '--groups-txt' — though only an "
        f"aggregation changes a group's map, since presence there is always the count of the "
        f"group's own samples, so 'count' and 'membership' are rejected in a grouped run as having "
        f"nothing to choose between. Maps for individual samples are never summarized, always "
        f"showing that one sample alone: its own values with a value column, its presence without "
        f"one. Without this option, samples are summarized by presence, by membership with 3 or "
        f"fewer samples and by count above that. Presence is the default because it is meaningful "
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
        "the 'unified' map of a grouped run. 'count' and 'membership' summarize presence, coloring "
        "an element by how many or exactly which groups contain it, a group containing it per "
        "'--group-threshold'. Any other name is an aggregation pooling the groups' values, each "
        "group's value being its own samples summarized by '--reaction-sample-summary' (which must "
        "therefore also pool values). Without this option, groups are summarized by presence, by "
        "membership with 3 or fewer groups and by count above that. A useful combination for "
        "replicates is '--reaction-sample-summary mean --reaction-group-summary count': each "
        "group's map shows the mean of its replicate samples, while the 'unified' map shows in how "
        "many groups each element occurs."
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
        "would be saved to a file named 'kos_00010.pdf', but, with this flag, the map would be "
        "saved to a file named 'kos_00010_Glycolysis_Gluconeogenesis.pdf'. To further illustrate "
        "how special characters in pathway names are translated to characters in file paths, the "
        "file name for 'Glycosylphosphatidylinositol (GPI)-anchor biosynthesis' would be "
        "'kos_00563_Glycosylphosphatidylinositol_(GPI)_anchor_biosynthesis.pdf', and the file name "
        "for 'Biosynthesis of 12-, 14- and 16-membered macrolides' would be "
        "'kos_00522_Biosynthesis_of_12_14_and_16_membered_macrolides.pdf' with this flag."
    )
    groupOUT.add_argument(
        '--categorize-files', action='store_true', default=False, help=
        "Categorize output files by pathway map within subdirectories corresponding to the BRITE "
        "hierarchy of maps (see https://www.genome.jp/brite/br08901). Symlinks to these files are "
        "also created for easier browsing. For example, if drawing a map file for 'Glycolysis / "
        "Gluconeogenesis', '00010', then the file will be written to subdirectory "
        "'Metabolism/Carbohydrate_metabolism' of the output directory, and a symlink with the same "
        "basename as the file will be created in subdirectory 'symlink' of the output directory. "
        "If drawing a map grid for 'Quorum sensing', '02024', then the file will be written to a "
        "subdirectory named 'grid/Cellular_Processes/Cellular_community_prokaryotes', and a "
        "symlink will be created in a subdirectory named 'grid/symlink'."
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
        "databases and groups can use two 'schemes' for dynamic coloring, 'by_count' and "
        "'by_membership' (see the argument, '--discrete-colormap-scheme'). As with pangenomes, "
        "'by_count' uses by default the 'plasma_r' colormap, trimming the top and bottom 10%%. "
        "'by_membership' uses by default the qualitative colormap, 'tab10', without trimming. This "
        "colormap contains distinct colors suitable for clearly differentiating the databases or "
        "groups containing reactions. All of the above concerns occurrence coloring. When coloring "
        "a reaction layer by a value column instead (a '--reaction-txt' file with a value column), "
        "this same option selects the sequential colormap that is sampled continuously to map "
        "reaction values (default 'plasma_r'), and '--discrete-colormap-scheme' does not apply. "
        "See the following webpage for named colormaps: "
        "https://matplotlib.org/stable/users/explain/colors/colormaps.html#classes-of-colormaps "
    )
    groupCOLOR.add_argument(
        '--discrete-colormap-scheme', choices=['by_count', 'by_membership'], help=
        "There are two ways of dynamically coloring elements by occurrence in multiple contigs "
        "databases ('--contigs-dbs'), the genomes of a pangenome, or groups ('--groups-txt') of "
        "either: by count or by membership. For a draw-kegg-pathways text file this choice is part "
        "of the sample and group summaries instead ('--reaction-sample-summary' and the related "
        "options), so this option is rejected for a text run. By default, with 4 or more databases "
        "or groups, reactions are colored by count of database or group; with 2 or 3, reactions "
        "are colored explicitly by database or group membership. In coloring by count, the "
        "colormap should be sequential, such that the color of a reaction changes 'smoothly' with "
        "the count. In contrast, coloring by membership means reaction color is determined by "
        "membership in a database/group or combination of databases/groups, so a qualitative "
        "colormap can be used instead of a sequential colormap, as by default with 2 or 3 "
        "categories, to give a distinct color to each membership category."
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
    groupGROUP.add_argument(
        '--group-colormap', nargs='+', help=
        "This option is like '--reaction-colormap', but only applies to drawing files for "
        "individual groups ('--draw-individual-files') and panels for individual groups in map "
        "grids ('--draw-grid'). These maps for individual groups show data membership in group "
        "sources, e.g., contigs databases, pan genomes, or the samples of a draw-kegg-pathways "
        "text file. They are always colored dynamically by count, e.g., the number of databases or "
        "genomes in the group containing the data. Like '--reaction-colormap', this parameter "
        "takes the name of a Matplotlib Colormap, and optionally, two decimal values between 0.0 "
        "and 1.0 to limit the fraction of the colormap used. The default configuration is the "
        "same, with the colormap being 'plasma_r' and the limits being 0.1 and 0.9. This applies "
        "only to individual group maps that show source counts, so for a draw-kegg-pathways text "
        "layer it applies when the layer's samples are summarized by presence "
        "('--reaction-sample-summary'/'--compound-sample-summary' of 'count' or 'membership', the "
        "default). When a layer's samples are summarized by pooling values, its individual group "
        "maps are colored continuously from '--reaction-colormap'/'--compound-colormap' instead, "
        "and no per-group source-count colorbar is written."
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
    if args.reaction_colormap is not None:
        unsupported_args.append('--reaction-colormap')
    if args.compound_colormap is not None:
        unsupported_args.append('--compound-colormap')
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
    if args.discrete_colormap_scheme is not None:
        unsupported_args.append('--discrete-colormap-scheme')
    if args.reaction_reverse_overlay:
        unsupported_args.append('--reaction-reverse-overlay')
    if args.group_colormap is not None:
        unsupported_args.append('--group-colormap')
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
        if len(colormap_arg) != 3:
            raise ConfigError(
                f"'{flag}' takes either a colormap name on its own, or a name followed by TWO "
                f"decimal limits on the fraction of the colormap to use, such as 'plasma 0.1 0.9'. "
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
                f"The two limits given to '{flag}' cut the colormap down to the fraction of it "
                f"between them, so they must lie between 0.0 and 1.0 with the smaller one first. "
                f"These do not: {limits[0]}, {limits[1]}."
            )
        return colormap_arg[0], limits

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

        # '--discrete-colormap-scheme' is superseded for text layers by the summaries, which choose
        # count vs. membership per layer and per level.
        if args.discrete_colormap_scheme is not None:
            raise ConfigError(
                "'--discrete-colormap-scheme' selects count vs. membership coloring for contigs "
                "databases, pangenome genomes, and groups of either. For a draw-kegg-pathways text "
                "file, that choice is made per layer and per level by the summary options: use "
                "'--reaction-sample-summary'/'--compound-sample-summary' with 'count' or "
                "'membership' for coloring across samples, and "
                "'--reaction-group-summary'/'--compound-group-summary' for coloring across groups."
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
                        f"values to pool. Summarize presence with 'count' or 'membership' instead, "
                        f"or add a value column to the file."
                    )
                if level == 'group' and args.groups_txt is None:
                    raise ConfigError(
                        f"'{flag}' summarizes groups of samples, but no groups were defined with "
                        f"'--groups-txt'."
                    )
                # With groups, presence at the SAMPLE level is always the count of a group's own
                # samples, so 'count' and 'membership' describe the same per-group map and neither
                # changes anything: the choice only shapes the ungrouped 'unified' map.
                if (
                    level == 'sample' and args.groups_txt is not None
                    and summary in PRESENCE_SUMMARIES
                ):
                    raise ConfigError(
                        f"'{flag}' was given as '{summary}', but the samples here are grouped with "
                        f"'--groups-txt', and each group's map shows the count of that group's own "
                        f"samples however presence is asked for, so 'count' and 'membership' would "
                        f"draw the same thing. Summarizing the samples by presence is already the "
                        f"default, so drop the option; use an aggregation such as 'mean' to color "
                        f"each group's map by pooled values instead, and note that "
                        f"'--{element_type}-group-summary' is what chooses count versus membership "
                        f"across the groups themselves."
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
                f"'--reaction-colormap'/'--compound-colormap' and "
                f"'--reaction-reverse-overlay'/'--compound-reverse-overlay' instead."
            )

    reaction_colormap, reaction_colormap_limits = _resolve_colormap(
        args.reaction_colormap, '--reaction-colormap'
    )
    compound_colormap, compound_colormap_limits = _resolve_colormap(
        args.compound_colormap, '--compound-colormap'
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
        'group_reverse_overlay': args.group_reverse_overlay,
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
    if group_colormap is not None:
        kwargs['group_colormap'] = group_colormap
        kwargs['group_colormap_limits'] = (
            group_colormap_limits if group_colormap_limits is not None else (0.1, 0.9)
        )
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

    if args.discrete_colormap_scheme is None:
        # The scheme is determined automatically by the number of contigs databases or groups.
        pass
    else:
        map_contigs_databases_kos = functools.partial(
            map_contigs_databases_kos, colormap_scheme=args.discrete_colormap_scheme
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

    if args.discrete_colormap_scheme is None:
        # The scheme is determined automatically by the number of genomes or groups.
        pass
    else:
        map_pan_database_kos = functools.partial(
            map_pan_database_kos, colormap_scheme=args.discrete_colormap_scheme
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
        draw_maps_lacking_data=args.draw_bare_maps
    )

def main() -> None:
    args = get_args()

    try:
        check_package_dependencies()
        check_kegg_data(args)
        consolidate_contigs_dbs(args)

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
        # Only a compound text file provides a compound layer, so every option that shapes one has
        # nothing to act on without it. The same goes for the sample and group summaries, which only
        # the text engine reads.
        if not is_txt_run:
            compound_only_flags = [flag for present_flag, flag in (
                (args.compound_color is not None, '--compound-color'),
                (args.compound_colormap is not None, '--compound-colormap'),
                (args.compound_reverse_overlay, '--compound-reverse-overlay'),
                (args.compound_accession_aggregation is not None,
                 '--compound-accession-aggregation'),
                (args.reaction_gene_aggregation is not None, '--reaction-gene-aggregation'),
                (args.reaction_accession_aggregation is not None,
                 '--reaction-accession-aggregation'),
                (args.reaction_sample_summary is not None, '--reaction-sample-summary'),
                (args.compound_sample_summary is not None, '--compound-sample-summary'),
                (args.reaction_group_summary is not None, '--reaction-group-summary'),
                (args.compound_group_summary is not None, '--compound-group-summary')
            ) if present_flag]
            if compound_only_flags:
                message = ', '.join(f"'{flag}'" for flag in compound_only_flags)
                raise ConfigError(
                    f"These options were given: {message}. They color a compound layer, or "
                    f"summarize the samples of a draw-kegg-pathways text file, both of which only "
                    f"'--reaction-txt'/'--compound-txt' provide. Database, pangenome and "
                    f"reaction-network-JSON inputs have a reaction layer alone, drawn from one "
                    f"source per database or genome, so none of these can apply."
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

        mapper = Mapper(kegg_dir=args.kegg_dir,
                        overwrite_output=args.overwrite_output_destinations,
                        name_files=args.name_files,
                        categorize_files=args.categorize_files)
        performed = False

        if is_txt_run:
            map_txt_data(args, mapper)
            performed = True

        if args.reaction_network_json is not None:
            map_json_network_ko_data(args, mapper)
            performed = True

        if (
            args.contigs_dbs is not None and
            len(args.contigs_dbs) == 1 and
            args.reaction_colormap is None
        ):
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
