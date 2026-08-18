#!/usr/bin/env python
"""Make KEGG pathway maps incorporating data sourced from external inputs."""

import os
import re
import fitz
import json
import math
import shutil
import numpy as np
import pandas as pd
import matplotlib.colors as mcolors

from argparse import Namespace
from itertools import combinations
from typing import Callable, Dict, Iterable, List, Literal, Set, Tuple, Union
# Colorbars are drawn with Matplotlib's object-oriented API rather than 'pyplot'. Importing
# 'pyplot' selects a backend, and with the pinned Matplotlib that happens at import time: on Linux
# it probes the X display, so merely running this program over an X-forwarded SSH connection (even
# just for '-h') opens an X client connection and fires up the user's X server. The classes below
# are equivalent, never touch a backend, and keep figures out of pyplot's global registry.
from matplotlib import colormaps
from matplotlib.figure import Figure
from matplotlib.cm import ScalarMappable

import anvio.kgml as kgml
import anvio.utils as utils
import anvio.dbinfo as dbinfo
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio import FORCE_OVERWRITE, QUIET, __version__ as VERSION

from anvio.errors import ConfigError
from anvio.metabolism.context import KeggContext
from anvio.dbops import ContigsDatabase, PanSuperclass
from anvio.colorconversions import blend_hexcodes, tint_hexcode
from anvio.metabolism.constants import GLOBAL_MAP_ID_PATTERN, OVERVIEW_MAP_ID_PATTERN


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2026, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = VERSION
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"
__status__ = "Development"


# The colors of qualitative and repeating colormaps are sampled in order, whereas other colormaps,
# including sequential colormaps, are sampled evenly.
qualitative_colormaps: List[str] = [
    'Pastel1',
    'Pastel2',
    'Paired',
    'Accent',
    'Dark2',
    'Set1',
    'Set2',
    'Set3',
    'tab10',
    'tab20',
    'tab20b',
    'tab20c'
]
repeating_colormaps: List[str] = [
    'flag',
    'prism'
]

# A colormap cut down to a fraction of itself is renamed by '_trim_colormap', which is the only
# place that name is written; this reads the name of the colormap it was cut from back out of it, so
# that what a trimmed colormap is made of is still knowable downstream. The two must be changed
# together.
TRIMMED_COLORMAP_PATTERN = re.compile(r'^trunc\((?P<name>.+),[0-9.]+,[0-9.]+\)$')

# Colormaps whose colors run from a neutral middle out to two opposed extremes, which is what makes
# the middle of a color scale mean anything and so what the '--*-value-center' options are for.
# Matplotlib groups these in its documentation but does not publish the grouping programmatically,
# so the names of its diverging class are listed here. A reversed colormap ('RdBu_r') diverges
# exactly as the colormap it reverses does, so only the base names are needed
# ('_is_diverging_colormap').
DIVERGING_COLORMAPS: Tuple[str, ...] = (
    'BrBG',
    'PRGn',
    'PiYG',
    'PuOr',
    'RdBu',
    'RdGy',
    'RdYlBu',
    'RdYlGn',
    'Spectral',
    'bwr',
    'coolwarm',
    'seismic'
)

# Functions for reducing a sequence of numeric values to a single value in quantitative coloring.
# Keys match the recommended choices of the '--*-gene-aggregation'/'--*-accession-aggregation'
# arguments, which reduce the values of a gene's rows into a per-accession value and the values of a
# map element's accessions into a per-element value, both within a sample; they also match the
# recommended value choices of the '--*-sample-summary'/'--*-group-summary' arguments, which pool the
# values of a set of samples or a set of sample groups. Any other pandas aggregation name is also
# accepted, resolved by '_resolve_aggregation'.
#
# These functions are fast paths that MUST compute the same statistic as the pandas aggregation of
# the same name, agreeing to floating-point precision (the two sum in different orders, so the last
# bit can differ). Per-accession values are computed by pandas ('_aggregate_accession_quantities'
# reduces a gene's rows with 'groupby.agg(name)') while every other level of the hierarchy applies
# these functions to a list, so were the two to disagree, one argument would silently mean two
# different statistics at different levels. Hence 'std' is the sample standard deviation (ddof=1), as
# in pandas, which is undefined for a single value: an undefined result means the accession or
# element has no value here, so it is dropped ('_finite_values') and left uncolored.
AGGREGATION_FUNCTIONS = {
    'sum': lambda values: float(np.sum(values)),
    'mean': lambda values: float(np.mean(values)),
    'max': lambda values: float(np.max(values)),
    'min': lambda values: float(np.min(values)),
    'median': lambda values: float(np.median(values)),
    'std': lambda values: float(np.std(values, ddof=1)) if len(values) > 1 else float('nan')
}

# The presence choices of the '--*-sample-summary'/'--*-group-summary' arguments, which summarize a
# set of samples or groups by how many ('count' and 'count_continuous') or exactly which
# ('membership') of them contain an accession, mapped to the colormap schemes that color them. These
# schemes are the choices of '--presence-colormap-scheme' and are unrelated to the sequential
# colormap that colors a value column; a summary named here is a presence summary, and any other is
# an aggregation of values. 'count' and 'count_continuous' color a count from the same colormap —
# identically, for the sequential colormap a count scale calls for — and differ only in how the
# scale is drawn: 'count' gives each count a band of a discrete colorbar, which requires one
# distinguishable color per category, while 'count_continuous' draws a gradient from the first count
# to the last, which any number of categories can share ('_membership_layer_colors').
SUMMARY_PRESENCE_SCHEMES = {
    'count': 'by_count',
    'count_continuous': 'by_count_continuous',
    'membership': 'by_membership'
}

# The presence summary names as a phrase for messages that list them, e.g. "'count',
# 'count_continuous' and 'membership'".
SUMMARY_PRESENCE_PHRASE = ' and '.join(
    part for part in (
        ', '.join(repr(name) for name in list(SUMMARY_PRESENCE_SCHEMES)[:-1]),
        repr(list(SUMMARY_PRESENCE_SCHEMES)[-1])
    ) if part
)

# How to ask for each presence colormap scheme, keyed by the scheme, for the inputs whose scheme
# '--presence-colormap-scheme' chooses: contigs databases, pan genomes, and groups of either. A
# draw-kegg-pathways text layer carries its own version of this map instead ('_build_txt_model'),
# since its scheme comes from the sample or group summary that colors the 'unified' map, whose
# presence names are the keys of 'SUMMARY_PRESENCE_SCHEMES' rather than the schemes themselves. The
# option that applies to one of these inputs is refused for the other, so a message naming the
# option that changes a layer's scheme has to take the wording from the layer.
PRESENCE_SCHEME_OPTIONS = {
    scheme: f'--presence-colormap-scheme {scheme}' for scheme in SUMMARY_PRESENCE_SCHEMES.values()
}

# The name a category colors file gives its color column. Its other column holds category names and
# can be headed anything, exactly as the item column of a groups-txt file can, so that one file can
# describe samples in one run and groups in another without being renamed.
CATEGORY_COLORS_COLUMN = 'color'

# What separates the names of a combination of categories in the first column of a category colors
# file, where a row can override the color that coloring by membership would otherwise derive for
# that combination by category color averaging. It matches how the membership colorbar labels a
# combination, so an override row can be copied straight from the scale it adjusts.
CATEGORY_COMBO_SEPARATOR = ','

# Coloring by membership needs a color for every possible combination of the categories, of which
# there are '2 ** n - 1'. With colors given per category the combinations are blended rather than
# sampled from a colormap, so there is no colormap size to bound them; this is the ceiling instead,
# above which the combinations are refused before they are enumerated. It is set at the size of a
# typical colormap, well past the handful of categories whose combinations any scale can tell apart.
MAX_CATEGORY_COLOR_COMBOS = 256

# The most counts a color scale is drawn in discrete bands of before the count is drawn on a
# continuous scale instead. A band's label is sized to fit the band
# ('ColorbarDrawer.draw_discrete'), so a scale of many counts labels them in type too small to read
# long before it runs out of colors to draw them in: this ceiling is where the labels stop fitting,
# not where the colors do. Past it a bar of bands reads as a gradient anyway, while a gradient labels
# the range rather than every count. Coloring by membership has no such ceiling because it has no
# gradient to fall back to — a color on one says which categories, which a position cannot.
MAX_DISCRETE_COUNT_BANDS = 40

# The '--group-colormap' value asking for each group's individual maps to be colored by a ramp built
# from that group's own color rather than sampled from a named Matplotlib colormap. Every group's
# ramp is built the same way, differing only in the hue it runs to, so the panels of a grid stay
# comparable while each says which group it is.
GROUP_COLORMAP_FROM_CATEGORY = 'category'

# The fractions of the way from white to a group's own color that its ramp runs between, used when
# 'GROUP_COLORMAP_FROM_CATEGORY' is asked for without explicit limits. The pale end stops short of
# white because white is a reserved color of standard and overview maps ('_check_reserved_colors'),
# and a ramp that reached it could not be drawn there at all.
DEFAULT_GROUP_TINT_SPAN = (0.25, 1.0)

# The fraction of a named group colormap used when no limits are given, trimming its darkest and
# lightest tenths.
DEFAULT_GROUP_COLORMAP_LIMITS = (0.1, 0.9)

# How to ask for each scheme of the count scale on individual group maps, keyed by the scheme, for
# messages that name the option ('_group_map_colors'). The schemes are those of a count scale
# elsewhere: 'by_count' gives each count a band of a discrete colorbar, which requires one
# distinguishable color per count, while 'by_count_continuous' draws a gradient from the lowest
# count to the highest, which any number of counts can share. 'by_membership' is not among them,
# since a group's own map counts that group's sources and never shows which of them contain an
# element.
GROUP_SCHEME_OPTIONS = {
    scheme: f'--group-colormap-scheme {scheme}'
    for scheme in ('by_count', 'by_count_continuous')
}

# How many colors a ramp built from a group's own color is sampled at to become a Colormap, which is
# what a continuous colorbar spans. The ramp is perceptual ('tint_hexcode') and Matplotlib
# interpolates in sRGB between the samples, so the count is set high enough that the difference
# between neighbors is invisible and the bar matches the colors drawn on the map.
GROUP_RAMP_COLORMAP_SIZE = 256

# What labels the end of a color scale that a value limit truncated ('_make_quantitative_norm'). The
# color there stands for that value or anything past it, since everything beyond the limit is drawn
# in it, and a bare number would claim the color means that value exactly.
CLAMPED_MIN_PREFIX = '≤ '
CLAMPED_MAX_PREFIX = '≥ '

# How close to a truncated end of a continuous colorbar, as a fraction of the bar, an automatic tick
# may fall before it is dropped in favor of the label marking that end ('draw_continuous'). The two
# would otherwise be set almost on top of each other, and the '≤'/'≥' label is the one a reader
# needs, so it is the one that stays.
MIN_TICK_SEPARATION_FRACTION = 0.08

# Subdirectories of the output directory, one per role a map file can have: the map unifying every
# source, the map of one individual source, and the grid comparing them. A fourth holds those
# individual maps a second time, arranged into a directory per map rather than per source. Every
# name directly in the output directory is therefore anvi'o's own, while the names that come from
# user data — samples, contigs databases, genomes, groups — are confined to 'individual', where
# they cannot collide with anything anvi'o creates for itself.
UNIFIED_SUBDIR = 'unified'
INDIVIDUAL_SUBDIR = 'individual'
GRID_SUBDIR = 'grid'
COLLATED_SUBDIR = 'by_map'

# Where '--categorize-files' nests map files in BRITE subdirectories, this subdirectory of links
# gathers every one of them back into a single place to browse ('_link_map_flat').
FLAT_SUBDIR = 'all_maps'

# The colorbar keying the maps of one individual source, written beside them in its directory.
CATEGORY_COLORBAR_BASENAME = 'colorbar.pdf'

class Mapper:
    """
    Make KEGG pathway maps incorporating data sourced from external inputs.

    Attributes
    ==========
    kegg_context : anvio.metabolism.context.KeggContext
        This contains anvi'o KEGG database attributes, such as filepaths.

    available_pathway_numbers : List[str]
        ID numbers of all pathways set up with PNG and KGML files in the KEGG data directory.

    pathway_names : Dict[str, str]
        The names of all KEGG pathways, including those without files in the KEGG data directory.
        Keys are pathway ID numbers and values are pathway names.

    xml_ops : anvio.kgml.XMLOps
        Used for loading KGML files as pathway objects.

    overwrite_output : bool
        If True, methods in this class overwrite existing output files.

    name_files : bool
        Include the pathway name along with the number in output map file names.

    categorize_files : bool
        Categorize output files by pathway map within subdirectories corresponding to the BRITE
        hierarchy of maps (see https://www.genome.jp/brite/br08901).

    collate_files_by_map : bool
        Alongside the maps drawn for each individual source in a directory of its own, gather those
        same maps into a directory per map, each holding one file per source.

    pathway_categorization : dict[str, list[str]]
        Maps pathway numbers to categorization, with categories listed from general to specific.

    run : anvio.terminal.Run
        This object prints run information to the terminal.

    progress : anvio.terminal.Progress
        This object prints transient progress information to the terminal.

    colorbar_drawer : ColorbarDrawer
        Writes standalone colorbar image files.

    grid_drawer : PDFGridDrawer
        Writes PDF files that are a grid of input PDF files.

    ignore_compound_rectangles : bool
        If True, do not draw KGML compound Entry rectangle Graphics. These are found in a small
        number of KGML files (see 00121, 00621, 01052, 01054), and when rendered by 'anvio.kgml' via
        'Bio.Graphics.KGML_vis.KGMLCanvas' have the effect of obscuring underlying drawings of
        compound structures in the base map image.
    """
    def __init__(
        self,
        kegg_dir: str = None,
        overwrite_output: bool = FORCE_OVERWRITE,
        name_files: bool = False,
        categorize_files: bool = False,
        collate_files_by_map: bool = False,
        run: terminal.Run = terminal.Run(),
        progress: terminal.Progress = terminal.Progress(),
        quiet: bool = QUIET
    ) -> None:
        """
        Parameters
        ==========
        kegg_dir : str, None
            Directory containing an anvi'o KEGG database. The default argument of None expects the
            KEGG database to be set up in the default directory used by the program
            anvi-setup-kegg-data.

        overwrite_output : bool, anvio.FORCE_OVERWRITE
            If True, methods in this class overwrite existing output files.

        name_files : bool, False
            Include the pathway name along with the number in output map file names.

        categorize_files : bool, False
            Categorize output files by pathway map within subdirectories corresponding to the BRITE
            hierarchy of maps (see https://www.genome.jp/brite/br08901).

        collate_files_by_map : bool, False
            Alongside the maps drawn for each individual source in a directory of its own, gather
            those same maps into a directory per map, each holding one file per source.

        run : anvio.terminal.Run, anvio.terminal.Run()
            This object prints run information to the terminal.

        progress : anvio.terminal.Progress, anvio.terminal.Progress()
            This object prints transient progress information to the terminal.

        quiet : bool, anvio.QUIET
            If True, run and progress information is not printed to the terminal.
        """
        args = Namespace()
        args.kegg_data_dir = kegg_dir
        self.kegg_context = KeggContext(args)

        if not os.path.exists(self.kegg_context.kegg_map_image_kgml_file):
            raise ConfigError(
                "One of the key files required by KEGG pathway maps is missing in your active "
                "anvi'o installation. If your KEGG data are not stored at the default KEGG data "
                "location, include that path using the 'kegg_dir' argument. Otherwise, please "
                "consider using the program `anvi-setup-kegg-data` to set up the latest KEGG data "
                "that includes the necessary files for KEGG pathway maps."
            )

        available_pathway_numbers: List[str] = []
        for row in pd.read_csv(
            self.kegg_context.kegg_map_image_kgml_file, sep='\t', index_col=0
        ).itertuples():
            if row.KO + row.EC + row.RN == 0:
                continue
            available_pathway_numbers.append(row.Index[-5:])
        self.available_pathway_numbers = available_pathway_numbers

        pathway_names: Dict[str, str] = {}
        for pathway_number, pathway_name in pd.read_csv(
            self.kegg_context.kegg_pathway_list_file, sep='\t', header=None
        ).itertuples(index=False):
            pathway_names[pathway_number[3:]] = pathway_name
        self.pathway_names = pathway_names

        self.xml_ops = kgml.XMLOps()
        self.drawer = kgml.Drawer(
            kegg_dir=self.kegg_context.kegg_data_dir, overwrite_output=overwrite_output
        )
        self.ignore_compound_rectangles = True
        self.colorbar_drawer = ColorbarDrawer(overwrite_output=overwrite_output)
        self.grid_drawer = PDFGridDrawer(overwrite_output=overwrite_output)

        self.name_files = name_files
        self.categorize_files = categorize_files
        self.collate_files_by_map = collate_files_by_map
        self.pathway_categorization = self._categorize_pathways() if categorize_files else None
        self.overwrite_output = overwrite_output
        self.run = run
        self.progress = progress
        self.quiet = quiet

    def map_contigs_database_kos(
        self,
        contigs_db: str,
        output_dir: str,
        pathway_numbers: Iterable[str] = None,
        reaction_color: str = '#2ca02c',
        draw_maps_lacking_data: bool = False
    ) -> Dict[str, bool]:
        """
        Draw pathway maps, highlighting KOs present in the contigs database.

        Parameters
        ==========
        contigs_db : str
            File path to a contigs database containing KO annotations.

        output_dir : str
            Path to the output directory in which pathway map PDF files are drawn. The directory is
            created if it does not exist.

        pathway_numbers : Iterable[str], None
            Regex patterns to match the ID numbers of the drawn pathway maps. The default of None
            draws all available pathway maps in the KEGG data directory.

        reaction_color : str, '#2ca02c'
            This is the color, by default green, for reaction elements represented by contigs
            database KOs. Alternatively to a color hex code, the string, 'original', can be provided
            to use the original color scheme of the reference map. In global and overview maps, KOs
            are represented in reaction lines. The foreground color of lines is set. In standard
            maps, KOs are represented in boxes, the background color of which is set, or lines.

        draw_maps_lacking_data : bool, False
            If False, by default, only draw maps containing any of the KOs in the contigs database.
            If True, draw maps regardless, meaning that nothing may be colored.

        Returns
        =======
        Dict[str, bool]
            Keys are pathway numbers. Values are True if the map was drawn, False if the map was not
            drawn because it did not contain any of the select KOs and 'draw_maps_lacking_data' was
            False.
        """
        # Retrieve the IDs of all KO annotations in the contigs database.
        self.progress.new("Loading KO data from the contigs database")
        self.progress.update("...")

        self._check_contigs_db(contigs_db)
        self._check_contigs_db_ko_annotation(contigs_db)

        cdb = ContigsDatabase(contigs_db)
        ko_ids = cdb.db.get_single_column_from_table(
            'gene_functions',
            'accession',
            unique=True,
            where_clause='source = "KOfam"'
        )
        self.progress.end()

        drawn = self._map_kos_fixed_colors(
            ko_ids,
            os.path.join(output_dir, UNIFIED_SUBDIR),
            pathway_numbers=pathway_numbers,
            color_hexcode=reaction_color,
            draw_maps_lacking_data=draw_maps_lacking_data
        )
        count = sum(drawn.values()) if drawn else 0
        self.run.info("Number of maps drawn", count)

        return drawn

    def map_reaction_network_json_kos(
        self,
        json_path: str,
        output_dir: str,
        pathway_numbers: Iterable[str] = None,
        reaction_color: str = '#2ca02c',
        draw_maps_lacking_data: bool = False
    ) -> Dict[str, bool]:
        """
        Draw pathway maps highlighting KOs present in a reaction network JSON file.

        The JSON file must be in the format produced by 'anvi-get-metabolic-model-file' or by
        'anvi-reaction-network --enzymes-txt ... --output-json ...'. KO IDs are extracted
        directly from the gene annotations in the JSON, so no reference databases are required.

        Parameters
        ==========
        json_path : str
            Path to an anvi'o reaction network JSON file.

        output_dir : str
            Path to the output directory in which pathway map PDF files are drawn.

        pathway_numbers : Iterable[str], None
            Regex patterns to match the ID numbers of the drawn pathway maps. The default of None
            draws all available pathway maps in the KEGG data directory.

        reaction_color : str, '#2ca02c'
            This is the color, by default green, for reaction elements represented by KOs from the
            network. Instead of a color hex code, the string, 'original', can be provided to use the
            original color scheme of the reference map. In global and overview maps, KOs are
            represented in reaction lines — the foreground color of lines is set. In standard maps,
            KOs are represented in boxes — the background color of which is set — or lines.

        draw_maps_lacking_data : bool, False
            If False, by default, only draw maps containing any of the KOs in the network. If True,
            draw maps regardless, meaning that nothing may be colored.

        Returns
        =======
        Dict[str, bool]
            Keys are pathway numbers. Values are True if the map was drawn, False if the map was not
            drawn because it did not contain any of the select KOs and 'draw_maps_lacking_data' was
            False.
        """
        filesnpaths.is_file_exists(json_path)

        self.progress.new("Loading KO data from the reaction network JSON")
        self.progress.update("...")

        with open(json_path) as f:
            json_dict = json.load(f)

        required_keys = {'genes', 'reactions', 'metabolites'}
        missing_keys = required_keys - set(json_dict)
        if missing_keys:
            self.progress.end()
            raise ConfigError(
                f"The file at '{json_path}' does not appear to be an anvi'o reaction network "
                f"JSON: it is missing the required top-level keys: "
                f"{', '.join(sorted(missing_keys))}."
            )

        ko_ids = set()
        for gene_entry in json_dict.get('genes', []):
            for ko_id in gene_entry.get('annotation', {}).get('ko', {}).keys():
                ko_ids.add(ko_id)

        self.progress.end()

        if not ko_ids:
            raise ConfigError(
                f"No KO annotations were found in the reaction network JSON at '{json_path}'. "
                f"There is nothing to draw."
            )

        self.run.info("KOs found in network JSON", len(ko_ids))

        drawn = self._map_kos_fixed_colors(
            ko_ids,
            os.path.join(output_dir, UNIFIED_SUBDIR),
            pathway_numbers=pathway_numbers,
            color_hexcode=reaction_color,
            draw_maps_lacking_data=draw_maps_lacking_data
        )
        count = sum(drawn.values()) if drawn else 0
        self.run.info("Number of maps drawn", count)

        return drawn

    def _read_element_txt(
        self,
        path: str,
        element_type: Literal['reaction', 'compound']
    ) -> dict:
        """
        Read and validate a per-layer draw-kegg-pathways text file into a map-layer table.

        Each file colors one layer. A 'reaction' file (artifact 'kegg-reaction-txt') colors reaction
        elements, and its accessions must be all KO IDs ('K#####') or all KEGG reaction IDs
        ('R#####'), not a mix; it may carry an optional 'gene_id' label column. A 'compound' file
        (artifact 'kegg-compound-txt') colors compound elements, and its accessions are all KEGG
        compound IDs ('C#####'); it must not carry a 'gene_id' column, since compounds are not gene
        products. In both files, the single column that is not a key column ('accession', 'gene_id'
        for reactions, or 'sample') is auto-detected as a numeric value column for quantitative
        coloring; with no such column the layer is colored by presence. An optional 'sample' column,
        if present, must be filled in every row.

        Parameters
        ==========
        path : str
            Path to the tab-delimited per-layer text file.

        element_type : Literal['reaction', 'compound']
            Which layer the file colors, selecting the accession-typing and key-column rules.

        Returns
        =======
        dict
            Keys: 'element_type' (the argument), 'reaction_source' ('KO'/'Reaction' for a reaction
            file, None for a compound file), 'df' (the validated rows, carrying the normalized
            '__accession' column plus '__sample' when a 'sample' column is present), 'value_column'
            (the auto-detected value column name, or None for presence coloring), and 'sample_names'
            (sorted list, or None if there is no 'sample' column).
        """
        # A single-column presence file (just 'accession') is valid, so 'is_file_tab_delimited'
        # (which rejects a line with no tab) is too strict here; the 'accession'-column check below
        # gives a clear error if the file is mis-delimited.
        filesnpaths.is_file_exists(path)

        self.progress.new(f"Loading the kegg-{element_type}-txt file")
        self.progress.update("...")

        # An empty file has no header row for pandas to find, which it reports as an exception rather
        # than as an empty table.
        try:
            df = pd.read_csv(path, sep='\t', dtype=str)
        except pd.errors.EmptyDataError:
            self.progress.end()
            raise ConfigError(
                f"The kegg-{element_type}-txt file at '{path}' is empty. It needs a header row with "
                f"an 'accession' column, and a row for each accession to color."
            )
        original_columns = list(df.columns)

        if 'accession' not in df.columns:
            self.progress.end()
            raise ConfigError(
                f"The kegg-{element_type}-txt file at '{path}' must have an 'accession' column. It "
                f"has these columns: {', '.join(df.columns)}."
            )

        # Every row needs a non-blank accession.
        accession = df['accession'].fillna('').astype(str).str.strip()
        if (accession == '').any():
            self.progress.end()
            raise ConfigError(
                f"Every row of the kegg-{element_type}-txt file at '{path}' must have a value in "
                f"the 'accession' column."
            )
        df = df.assign(__accession=accession)

        # Type the accessions by their KEGG ID prefix. A reaction file must be all 'K' (KO) or all
        # 'R' (KEGG reaction) accessions, not a mix, since both color reaction elements and mixing
        # them could clash on the same elements; a compound file must be all 'C' accessions.
        reaction_source = None
        if element_type == 'reaction':
            is_ko = df['__accession'].str.match(r'^K\d+$')
            is_reaction = df['__accession'].str.match(r'^R\d+$')
            if is_ko.all():
                reaction_source = 'KO'
            elif is_reaction.all():
                reaction_source = 'Reaction'
            elif is_ko.any() and is_reaction.any():
                self.progress.end()
                raise ConfigError(
                    f"The kegg-reaction-txt file at '{path}' mixes KO accessions ('K' followed by "
                    f"digits) and KEGG reaction accessions ('R' followed by digits), but a "
                    f"reaction file must contain only one type: all KOs or all reactions. Both "
                    f"color the reaction elements of a map, and mixing them could create clashes "
                    f"on the same elements."
                )
            else:
                self.progress.end()
                bad = sorted(set(df.loc[~(is_ko | is_reaction), '__accession']))[:5]
                raise ConfigError(
                    f"The accessions in the kegg-reaction-txt file at '{path}' must all be KO IDs "
                    f"('K' followed by digits, e.g., 'K00844') or all be KEGG reaction IDs ('R' "
                    f"followed by digits, e.g., 'R00200'). These accessions match neither: "
                    f"{', '.join(bad)}."
                )
        else:
            is_compound = df['__accession'].str.match(r'^C\d+$')
            if not is_compound.all():
                self.progress.end()
                bad = sorted(set(df.loc[~is_compound, '__accession']))[:5]
                raise ConfigError(
                    f"The accessions in the kegg-compound-txt file at '{path}' must all be KEGG "
                    f"compound IDs ('C' followed by digits, e.g., 'C00031'). These accessions do "
                    f"not: {', '.join(bad)}."
                )

        # A compound file cannot carry a gene-level identifier, since compounds are not gene
        # products.
        if element_type == 'compound' and 'gene_id' in df.columns:
            self.progress.end()
            raise ConfigError(
                f"The kegg-compound-txt file at '{path}' must not have a 'gene_id' column, since "
                f"compounds are not gene products. A 'gene_id' column labeling the gene a value "
                f"comes from is only meaningful for the reaction layer (kegg-reaction-txt)."
            )

        # Auto-detect the value column: the single column that is not a key column ('accession',
        # 'gene_id' for reactions, or 'sample'). None means presence coloring; one means
        # quantitative coloring by that column (its header labels the colorbar); two or more is
        # ambiguous.
        key_columns = {'accession', 'gene_id', 'sample'}
        value_columns = [column for column in original_columns if column not in key_columns]
        if len(value_columns) > 1:
            self.progress.end()
            raise ConfigError(
                f"The kegg-{element_type}-txt file at '{path}' has more than one candidate value "
                f"column, so it is ambiguous which to color elements by: "
                f"{', '.join(value_columns)}. A file may have at most one value column (any column "
                f"that is not a key column): include zero to color by presence/absence, or one to "
                f"color by that numeric value."
            )
        value_column = value_columns[0] if value_columns else None

        # An optional 'sample' column, if present, must be filled in every row, since the sample is
        # the origin used to color elements across samples.
        sample_names = None
        if 'sample' in df.columns:
            sample = df['sample'].fillna('').astype(str).str.strip()
            if (sample == '').any():
                self.progress.end()
                raise ConfigError(
                    f"When the kegg-{element_type}-txt file at '{path}' has a 'sample' column, "
                    f"every row must have a sample value, since the sample is the origin used to "
                    f"color elements across samples. Some rows have a blank 'sample' value."
                )
            df = df.assign(__sample=sample)
            sample_names = sorted(set(sample))

        # With a value column, each thing the file describes must be given one value. Two rows for
        # the same thing are ambiguous — repeated measurements to be averaged, or separate
        # contributions to be added? — and only whoever produced the data can say, so they are
        # refused rather than combined by a rule the user did not choose. What counts as "the same
        # thing" is the file's own key: a reaction file naming genes describes one gene's value for
        # one accession, so several genes carrying one accession are not repeats and remain the
        # ordinary case.
        if value_column is not None:
            key_columns = [
                column for column in ('__accession', 'gene_id', '__sample') if column in df.columns
            ]
            duplicated = df.duplicated(subset=key_columns, keep=False)
            if duplicated.any():
                described = ' and '.join(
                    column.lstrip('_') for column in key_columns
                ) if len(key_columns) > 1 else 'accession'
                examples = df.loc[duplicated, key_columns].drop_duplicates().head(3)
                shown = '; '.join(
                    ', '.join(str(value) for value in row) for row in examples.itertuples(index=False)
                )
                self.progress.end()
                raise ConfigError(
                    f"Every row of the kegg-{element_type}-txt file at '{path}' must describe a "
                    f"different thing, since each carries a value in the '{value_column}' column, but "
                    f"{int(duplicated.sum())} rows repeat a combination of {described}. The first "
                    f"repeated {'combinations are' if len(key_columns) > 1 else 'accessions are'}: "
                    f"{shown}. Anvi'o will not guess how repeated rows should be combined -- averaging "
                    f"replicate measurements and adding up separate contributions, such as the ions of "
                    f"one metabolite, are both reasonable and give different answers. Please combine "
                    f"them in the way that suits your data before drawing."
                )

        self.progress.end()

        return {
            'element_type': element_type,
            'reaction_source': reaction_source,
            'df': df,
            'value_column': value_column,
            'sample_names': sample_names
        }

    @staticmethod
    def _resolve_aggregation(aggregation: str, flag: str) -> Callable:
        """
        Resolve an aggregation name into a function reducing a sequence of values to one value.

        The recommended names have fast paths in 'AGGREGATION_FUNCTIONS'. Any other pandas
        aggregation name is accepted as well, provided it reduces a series to a single number: the
        name is probed here, on a series, before it is applied to any data, so that a name pandas
        does not recognize, or one that transforms rather than reduces, is reported as a
        configuration error rather than failing deep in a drawing loop. Probing a series is what
        makes the name safe at every level of the reduction hierarchy, since the levels below the
        per-accession values reduce plain lists: an aggregation that only a groupby offers, such as
        'first', is rejected here rather than working for a gene's rows and failing for a map
        element's accessions. A name that does reduce but means something unrelated to the value
        column (pandas 'count', say, which counts rows) is accepted; the colorbar is still labeled
        by the value column, so choosing such a name is the user's business.

        Parameters
        ==========
        aggregation : str
            The requested aggregation name.

        flag : str
            The command-line flag the name comes from, used in an error message.

        Returns
        =======
        Callable
            Reduces a sequence of values to a single value.
        """
        try:
            return AGGREGATION_FUNCTIONS[aggregation]
        except KeyError:
            pass

        def aggregate(values):
            return float(pd.Series(values, dtype=float).agg(aggregation))

        # The name is probed at BOTH levels it will be used at, and the two results are compared.
        # Per accession the name goes to pandas as a string through a groupby
        # ('_aggregate_accession_quantities'), while every level above applies this function to a
        # list, and the two accept different sets of names: 'first' exists only on a groupby, and
        # 'idxmax' means the position within a list at one level but a row label of the whole table
        # at the other. Probing one level alone would let such a name through to fail, or silently
        # disagree, on real data. Anything pandas raises for an unknown name, and the TypeError from
        # a name returning a series rather than a number, become one clear error. The probe table is
        # given an index that does not start at 0, so that a name returning a row label rather than
        # a value ('idxmax') disagrees with the list-level result and is caught; every real
        # reduction ignores the index.
        probe_values = [1.0, 2.0, 4.0]
        try:
            probe = aggregate(probe_values)
            grouped_probe = float(
                pd.DataFrame(
                    {'k': ['x'] * len(probe_values), 'v': probe_values},
                    index=range(10, 10 + len(probe_values))
                ).groupby('k')['v'].agg(aggregation).iloc[0]
            )
        except Exception:
            probe = None
        if probe is None or not np.isfinite(probe) or not np.isclose(
            probe, grouped_probe, rtol=1e-12, atol=0
        ):
            raise ConfigError(
                f"'{flag}' was given as '{aggregation}', which either does not reduce values to a "
                f"single number or does not reduce them the same way at every level. Proven "
                f"acceptable names are "
                f"{', '.join(repr(name) for name in AGGREGATION_FUNCTIONS)}; any other pandas "
                f"aggregation that reduces a series to one number, such as 'var' or 'sem', also "
                f"works. Names that transform rather than reduce, like 'cumsum', cannot be used; "
                f"neither can names offered only by a grouping, like 'first', nor names meaning "
                f"different things for a list of values and for a table, like 'idxmax'. Note that "
                f"the sample and group summary options additionally take "
                f"{SUMMARY_PRESENCE_PHRASE} to summarize presence rather than value."
            )
        return aggregate

    @staticmethod
    def _finite_values(values: Dict[str, float], undefined: Set[str] = None) -> Dict[str, float]:
        """
        Drop the accessions whose aggregated value is not a finite number.

        An aggregation can be undefined for the values it is given — the standard deviation of a
        single value, say — and a non-finite value has no place on a color scale, so such an
        accession is treated as having no value and is left uncolored, exactly as an accession
        absent from the file would be.

        Parameters
        ==========
        values : Dict[str, float]
            Keys are accessions, values are aggregated values.

        undefined : Set[str], None
            Accessions that were dropped are added to this set, so that a caller can report them
            once for the whole layer (see '_warn_undefined_values').

        Returns
        =======
        Dict[str, float]
            The input without the accessions whose value is not a finite number.
        """
        dropped = {
            accession for accession, value in values.items()
            if not (isinstance(value, (int, float, np.number)) and np.isfinite(value))
        }
        if not dropped:
            return values
        if undefined is not None:
            undefined.update(dropped)
        return {
            accession: value for accession, value in values.items() if accession not in dropped
        }

    def _warn_undefined_values(self, undefined: Set[str], aggregation: str, path: str) -> None:
        """
        Report the accessions of one layer whose aggregated value was undefined.

        Parameters
        ==========
        undefined : Set[str]
            Accessions dropped by '_finite_values'. Nothing is reported if this is empty.

        aggregation : str
            The name of the aggregation that was undefined, for the message.

        path : str
            Path to the layer's text file, for the message.
        """
        if not undefined:
            return
        examples = ', '.join(sorted(undefined)[:5])
        self.run.warning(
            f"Reducing the values of the text file at '{path}' with '{aggregation}' was undefined "
            f"for {len(undefined)} accession(s), including these: {examples}. This happens when an "
            f"aggregation needs more values than are available, as the standard deviation does for "
            f"a single value. These accessions are treated as having no value, so the map elements "
            f"that depend on them are left uncolored."
        )

    def _aggregate_accession_quantities(
        self,
        rows_df: pd.DataFrame,
        value_column: str,
        aggregation: str,
        path: str,
        undefined: Set[str] = None
    ) -> Dict[str, float]:
        """
        Validate and aggregate a layer file's auto-detected value column into per-accession values.
        This aggregates each accession across its rows (e.g., a reaction file's per-gene rows).

        Values are coerced to numbers; a blank or non-numeric value in any row is an error. Each
        accession's rows are reduced to a single value by 'aggregation', which is applied by pandas
        here and by the matching function from 'AGGREGATION_FUNCTIONS' at every other level of the
        reduction hierarchy. An accession whose result is undefined is dropped ('_finite_values').

        Parameters
        ==========
        rows_df : pandas.DataFrame
            Rows of one layer, carrying the normalized '__accession' column and 'value_column'.

        value_column : str
            Name of the auto-detected numeric value column.

        aggregation : str
            How to reduce an accession's rows to a single value, as a pandas aggregation name.

        path : str
            Path to the layer's text file, used in error messages.

        undefined : Set[str], None
            Accessions whose aggregated value is undefined are added to this set.

        Returns
        =======
        Dict[str, float]
            Keys are accessions, values are aggregated numeric values.
        """
        numeric_values = pd.to_numeric(rows_df[value_column], errors='coerce')
        # A value must be a finite number: 'coerce' turns blanks/non-numerics into NaN, but leaves
        # 'inf'/'-inf' as infinities that would break the color scale, so reject those too.
        invalid = ~np.isfinite(numeric_values)
        if invalid.any():
            bad_accessions = sorted(set(rows_df.loc[invalid, '__accession']))
            self.progress.end()
            raise ConfigError(
                f"The '{value_column}' value column of the text file at '{path}' must contain a "
                f"finite number in every row for elements to be colored quantitatively. However, "
                f"{int(invalid.sum())} {'rows have' if invalid.sum() > 1 else 'row has'} a blank, "
                f"non-numeric, or infinite value, including these accessions: "
                f"{', '.join(bad_accessions[:5])}."
            )
        valued_df = rows_df.assign(__quantitative_value=numeric_values)
        return self._finite_values(
            valued_df.groupby('__accession')['__quantitative_value'].agg(aggregation).to_dict(),
            undefined
        )

    @staticmethod
    def _summarize_category_values(
        category_values: Dict[str, Dict[str, float]],
        aggregate: Callable,
        undefined: Set[str] = None
    ) -> Dict[str, float]:
        """
        Reduce per-category values to a single value per accession.

        This is the value form of a sample or group summary: for each accession, the values of the
        categories (samples or groups) containing it are reduced by 'aggregate'. A category lacking
        the accession does not contribute, so the summary is over the categories that have it.

        Parameters
        ==========
        category_values : Dict[str, Dict[str, float]]
            Keys are category names, values are {accession: value} dictionaries.

        aggregate : Callable
            Reduces a list of values to a single value (see 'AGGREGATION_FUNCTIONS').

        undefined : Set[str], None
            Accessions whose summarized value is undefined are added to this set.

        Returns
        =======
        Dict[str, float]
            Keys are accessions, values are the summarized values.
        """
        accession_values: Dict[str, List[float]] = {}
        for values in category_values.values():
            for accession, value in values.items():
                accession_values.setdefault(accession, []).append(value)
        return Mapper._finite_values(
            {accession: aggregate(values) for accession, values in accession_values.items()},
            undefined
        )

    def _check_category_names(self, categories: Iterable[str], category_noun: str) -> None:
        """
        Check that category names can serve as output subdirectory names.

        A category (sample, source, or group) drawn on its own maps gets its own subdirectory of the
        output directory, and once map grids are drawn the subdirectories of categories that were
        only needed for a grid are deleted ('_draw_map_grids'). A name that is not a plain directory
        name would therefore write, or delete, outside the output directory: names come from a text
        file's 'sample' column or a groups file, so they cannot be trusted to be safe paths. The
        same names become file names where the maps are gathered by map rather than by category
        ('_collate_maps_by_map'), which the very same requirement makes safe.

        Only the categories actually getting maps are checked, since a category summarized on the
        'unified' map alone contributes color rather than a path, and its name is never joined onto
        one.

        No name is reserved. Categories are drawn into '<output directory>/individual', which anvi'o
        creates for that purpose alone, so a category may be named after anything anvi'o puts in the
        output directory itself — 'unified', 'grid', 'by_map', 'all_maps', or a BRITE category such
        as 'Metabolism' — without the two ever meeting.

        Parameters
        ==========
        categories : Iterable[str]
            The category names that would become subdirectory names.

        category_noun : str
            What a category is, for the error message, e.g. 'sample' or 'group'.
        """
        separators = {os.sep, os.altsep} - {None}
        problems: List[str] = []
        for category in categories:
            if not category or not category.strip() or category in ('.', '..'):
                reason = "is not a usable directory name"
            elif any(separator in category for separator in separators):
                reason = "contains a path separator"
            elif os.path.isabs(category) or os.path.basename(category) != category:
                reason = "is not a plain directory name"
            else:
                continue
            problems.append(f"'{category}' ({reason})")

        if not problems:
            return
        raise ConfigError(
            f"Each {category_noun} drawn on its own maps gets its own subdirectory of the output "
            f"directory, so its name must be usable as a directory name. "
            f"{'These names cannot be used' if len(problems) > 1 else 'This name cannot be used'}: "
            f"{', '.join(problems)}. Please rename {'them' if len(problems) > 1 else 'it'} in the "
            f"input, using single words without path separators, such as 'SAMPLE_1' or "
            f"'HIGH_TEMPERATURE'. Alternatively, leave {'them' if len(problems) > 1 else 'it'} out "
            f"of the {category_noun}s requested with '--draw-individual-files'/'--draw-grid': a "
            f"{category_noun} that is only summarized on the 'unified' map does not need a name "
            f"that works as a directory."
        )

    @staticmethod
    def _relate_accessions_to_samples(
        df: pd.DataFrame,
        all_sample_names: List[str]
    ) -> Tuple[Dict[str, List[str]], Dict[str, Set[str]]]:
        """
        Relate a layer's accessions to the samples containing them, and each sample to its accessions.

        Parameters
        ==========
        df : pandas.DataFrame
            Rows of one layer, carrying the normalized '__accession' and '__sample' columns.

        all_sample_names : List[str]
            Names of all samples across the run's input files, so that a sample contributing no rows
            to this layer still gets an empty entry in 'source_accessions'.

        Returns
        =======
        Tuple[Dict[str, List[str]], Dict[str, Set[str]]]
            membership : maps each accession to the sorted names of the samples containing it.
            source_accessions : maps each sample name to its set of accessions.
        """
        membership: Dict[str, Set[str]] = {}
        source_accessions: Dict[str, Set[str]] = {s: set() for s in all_sample_names}
        for accession, sample_name in zip(df['__accession'], df['__sample']):
            source_accessions[sample_name].add(accession)
            membership.setdefault(accession, set()).add(sample_name)
        return (
            {accession: sorted(samples) for accession, samples in membership.items()},
            source_accessions
        )

    @staticmethod
    def _resolve_summary(
        summary: Union[str, None],
        value_column: Union[str, None],
        flag: str,
        path: str
    ) -> Tuple[Literal['value', 'presence'], Union[Callable, None], Union[str, None]]:
        """
        Resolve a sample or group summary into a coloring kind and its reduction.

        A summary reduces a set of samples, or a set of sample groups, to one statement per
        accession. The names of 'SUMMARY_PRESENCE_SCHEMES' summarize presence: how many
        ('count'/'count_continuous'), or exactly which ('membership'), categories contain the
        accession. Any other name pools the categories' values as an aggregation
        ('_resolve_aggregation'), which requires the layer's file to have a value column. The
        default of None summarizes presence with the colormap scheme left unresolved, so that
        '_membership_layer_colors' picks it from the number of categories (by membership for ≤ 3, by
        count > 3, and by a continuous count scale where a discrete one would run out of
        distinguishable colors or of room to label its bands, the latter past
        'MAX_DISCRETE_COUNT_BANDS' of them).

        Parameters
        ==========
        summary : Union[str, None]
            The requested summary, or None for presence with an unresolved scheme.

        value_column : Union[str, None]
            The layer's value column name, or None if it has none.

        flag : str
            The command-line flag this summary comes from, used in an error message.

        path : str
            Path to the layer's text file, used in an error message.

        Returns
        =======
        Tuple[Literal['value', 'presence'], Union[Callable, None], Union[str, None]]
            The kind of coloring, the value reduction function (None for presence), and the presence
            colormap scheme (None for a value summary or an unresolved presence scheme).
        """
        if summary is None:
            return 'presence', None, None
        if summary in SUMMARY_PRESENCE_SCHEMES:
            return 'presence', None, SUMMARY_PRESENCE_SCHEMES[summary]
        if value_column is None:
            raise ConfigError(
                f"'{flag}' was given as '{summary}', which pools the values of samples or sample "
                f"groups, but the text file at '{path}' has no value column, so there are no "
                f"values to pool. Summarize presence with {SUMMARY_PRESENCE_PHRASE} instead, or "
                f"add a value column to the file."
            )
        return 'value', Mapper._resolve_aggregation(summary, flag), None

    @staticmethod
    def _resolve_value_limits(
        limits: Union[Tuple[Union[float, None], Union[float, None]], None], flag: str
    ) -> Union[Tuple[Union[float, None], Union[float, None]], None]:
        """
        Validate the limits a color scale of values may span.

        A limit truncates a scale only where the values cross it ('_make_quantitative_norm'), so
        either end can be left open with None. A pair leaving both ends open would do nothing at
        all, which is likelier a mistake than an intention, and so is refused rather than accepted
        in silence.

        Parameters
        ==========
        limits : Union[Tuple[Union[float, None], Union[float, None]], None]
            The (minimum, maximum) requested for the scale, either of which can be None to leave
            that end wherever the values put it. None means no limits were requested at all.

        flag : str
            The command-line flag the limits come from, used in error messages.

        Returns
        =======
        Union[Tuple[Union[float, None], Union[float, None]], None]
            The limits as floats, either or both of them None, or None if none were requested.
        """
        if limits is None:
            return None
        try:
            limit_min, limit_max = limits
        except (TypeError, ValueError):
            raise ConfigError(
                f"'{flag}' takes two limits, a minimum and a maximum, either of which can be left "
                f"open. Anvi'o could not read a pair of them out of this: {limits}"
            )
        resolved = []
        for limit, end in ((limit_min, 'minimum'), (limit_max, 'maximum')):
            if limit is None:
                resolved.append(None)
                continue
            try:
                number = float(limit)
            except (TypeError, ValueError):
                raise ConfigError(
                    f"The {end} given to '{flag}' must be a number, or left open, but anvi'o got "
                    f"this instead: '{limit}'."
                )
            if not np.isfinite(number):
                raise ConfigError(
                    f"The {end} given to '{flag}' must be a finite number, since it is a place a "
                    f"color scale stops, but anvi'o got this instead: '{limit}'."
                )
            resolved.append(number)
        limit_min, limit_max = resolved
        if limit_min is None and limit_max is None:
            raise ConfigError(
                f"'{flag}' was given with neither a minimum nor a maximum, so there is nothing for "
                f"it to limit. Please give a number to at least one of its two ends, or drop the "
                f"option."
            )
        if limit_min is not None and limit_max is not None and limit_min >= limit_max:
            raise ConfigError(
                f"The minimum given to '{flag}' ({limit_min:g}) must be less than its maximum "
                f"({limit_max:g}). A color scale runs from its minimum up to its maximum, so the "
                f"two cannot be equal, nor the wrong way around."
            )
        return limit_min, limit_max

    @staticmethod
    def _resolve_value_center(
        center: Union[float, str, None],
        flag: str,
        limits: Union[Tuple[Union[float, None], Union[float, None]], None] = None,
        limits_flag: Union[str, None] = None
    ) -> Union[float, None]:
        """
        Validate the value a color scale is centered on.

        A centered scale runs the same distance either side of this value
        ('_make_quantitative_norm'), which puts it at the middle of the colormap however lopsided
        the values around it are. A center lying outside the limits bounding the same scale is
        refused here, before any value is read: those limits declare it out of the scale's reach, so
        centering that scale on it cannot be what was meant.

        Parameters
        ==========
        center : Union[float, str, None]
            The value requested at the middle of the scale, or None for a scale centered wherever
            its own values leave it.

        flag : str
            The command-line flag the center comes from, used in error messages.

        limits : Union[Tuple[Union[float, None], Union[float, None]], None], None
            The limits bounding the same scale, as '_resolve_value_limits' returns them, or None
            where the scale is not limited.

        limits_flag : Union[str, None], None
            The command-line flag those limits come from, used in an error message.

        Returns
        =======
        Union[float, None]
            The center as a float, or None if none was requested.
        """
        if center is None:
            return None
        try:
            number = float(center)
        except (TypeError, ValueError):
            raise ConfigError(
                f"'{flag}' takes the one value that a color scale is centered on, which must be a "
                f"number, but anvi'o got this instead: '{center}'. Give the option no value at all "
                f"to center the scale on zero."
            )
        if not np.isfinite(number):
            raise ConfigError(
                f"The value given to '{flag}' must be a finite number, since it is a place in the "
                f"middle of a color scale, but anvi'o got this instead: '{center}'."
            )
        if limits is not None:
            limit_min, limit_max = limits
            if limit_min is not None and number < limit_min:
                raise ConfigError(
                    f"'{flag}' centers a color scale on {number:g}, while '{limits_flag}' gives "
                    f"that same scale a minimum of {limit_min:g}, which lies above it. A scale "
                    f"cannot be centered on a value it is not allowed to reach. Please move one of "
                    f"the two."
                )
            if limit_max is not None and number > limit_max:
                raise ConfigError(
                    f"'{flag}' centers a color scale on {number:g}, while '{limits_flag}' gives "
                    f"that same scale a maximum of {limit_max:g}, which lies below it. A scale "
                    f"cannot be centered on a value it is not allowed to reach. Please move one of "
                    f"the two."
                )
        return number

    @staticmethod
    def _base_colormap_name(cmap: mcolors.Colormap) -> str:
        """
        Return the name of the colormap this one was made from.

        A colormap cut down to a fraction of itself is renamed after the colormap it was cut from
        ('_trim_colormap'), and that wrapper is what this unwinds, so that messages can name the
        colormap the user actually asked for.

        Parameters
        ==========
        cmap : matplotlib.colors.Colormap
            The colormap to name.

        Returns
        =======
        str
            The name of the colormap it was cut from, or its own name if it was not cut.
        """
        trimmed = TRIMMED_COLORMAP_PATTERN.match(cmap.name)
        return cmap.name if trimmed is None else trimmed.group('name')

    @staticmethod
    def _is_diverging_colormap(cmap: mcolors.Colormap) -> bool:
        """
        Report whether a colormap runs from a neutral middle out to two opposed extremes.

        Only such a colormap gives the middle of a color scale a meaning of its own, which is what
        the '--*-value-center' options put a value at. The name is matched against
        'DIVERGING_COLORMAPS' with any '_r' reversal suffix dropped, a reversed colormap diverging
        exactly as the colormap it reverses does, and with the wrapper a trimmed colormap carries
        unwound first ('TRIMMED_COLORMAP_PATTERN'), a cut of a diverging colormap being made of the
        same two ramps. Whether the trimming left the middle of what remains neutral is a separate
        question, asked where this is used.

        Parameters
        ==========
        cmap : matplotlib.colors.Colormap
            The colormap to classify.

        Returns
        =======
        bool
            True if the colormap is one of the diverging ones.
        """
        name = Mapper._base_colormap_name(cmap)
        return (name[:-2] if name.endswith('_r') else name) in DIVERGING_COLORMAPS

    def _check_centered_colormap(
        self,
        cmap: mcolors.Colormap,
        colormap_limits: Union[Tuple[float, float], None],
        center: Union[float, None],
        center_flag: str,
        colormap_flag: str,
        subject: str
    ) -> None:
        """
        Warn where centering a color scale asks more of a colormap than it can give.

        Putting a value at the middle of a scale only says something to a reader if the middle of
        the colormap looks like a middle: a diverging colormap's neutral center, with its two
        extremes either side. A sequential colormap has no such landmark, and trimming a diverging
        one off-center ('_trim_colormap') moves the neutral color away from the middle of what is
        actually drawn. Neither is an error, since the scale is still centered where it was asked to
        be, and the colorbar labels the center either way.

        Parameters
        ==========
        cmap : matplotlib.colors.Colormap
            The colormap coloring the centered scale, already trimmed.

        colormap_limits : Union[Tuple[float, float], None]
            The fractions of the colormap that were kept, or None for all of it.

        center : Union[float, None]
            The value put at the middle of the scale, or None where the scale is not centered.

        center_flag : str
            The command-line flag the center comes from.

        colormap_flag : str
            The command-line flag the colormap comes from.

        subject : str
            What the colormap colors, e.g. "reactions on the 'unified' map".
        """
        if center is None:
            return
        colormap_name = self._base_colormap_name(cmap)
        if not self._is_diverging_colormap(cmap):
            self.run.warning(
                f"'{center_flag}' puts {center:g} at the middle of the color scale of {subject}, "
                f"but the colormap coloring that scale, '{colormap_name}', is not a diverging one: "
                f"its colors run from one end to the other rather than out from a neutral middle, "
                f"so there is no color there for the centered value to take, and a reader cannot "
                f"see where the middle of the scale is except from the colorbar. A diverging "
                f"colormap given to '{colormap_flag}' — e.g., 'RdYlGn', 'RdYlBu_r — is what makes "
                f"a centered scale legible."
            )
        elif colormap_limits is not None and abs(
            (colormap_limits[0] + colormap_limits[1]) / 2 - 0.5
        ) > 1e-9:
            self.run.warning(
                f"'{center_flag}' puts {center:g} at the middle of the color scale of {subject}, "
                f"which is drawn in the middle color of '{colormap_name}' as it was trimmed by "
                f"'{colormap_flag}' to the fraction between {colormap_limits[0]:g} and "
                f"{colormap_limits[1]:g}. That trimming is not symmetric, so the neutral color "
                f"this diverging colormap has at its own middle is no longer in the middle of what "
                f"is drawn, and the centered value takes some other color of the two ramps "
                f"instead. Trim the same amount off each end, or none at all, to keep the neutral "
                f"color where the centered value is."
            )

    def _build_txt_model(
        self,
        name: str,
        data: dict,
        path: str,
        gene_aggregation: str,
        accession_aggregation: str,
        color: str,
        colormap: Union[bool, str, mcolors.Colormap, None],
        colormap_limits: Union[Tuple[float, float], None],
        category_colormap: Union[str, mcolors.Colormap, None],
        category_colormap_limits: Union[Tuple[float, float], None],
        category_colors_txt: Union[str, None],
        reverse_overlay: bool,
        all_sample_names: Union[List[str], None],
        group_samples: Union[Dict[str, List[str]], None],
        sample_summary: Union[str, None],
        group_summary: Union[str, None],
        value_limits: Union[Tuple[Union[float, None], Union[float, None]], None] = None,
        category_value_limits: Union[Tuple[Union[float, None], Union[float, None]], None] = None,
        value_center: Union[float, None] = None,
        category_value_center: Union[float, None] = None
    ) -> dict:
        """
        Build one layer's coloring model for '_map_elements' from a per-layer reader result.

        The layer colors elements by mode in each of two map contexts: 'unified_mode' for the
        'unified' map and 'category_mode' for the per-sample or per-group maps. Without a 'sample'
        column there is nothing to summarize, so the two agree: 'quantitative' with a value column
        ('_read_element_txt'), or 'single' (one fixed presence color) without one. An
        '--original-color' layer is the exception, always pairing a 'unified_mode' of 'original' with
        a 'category_mode' of 'membership', for the reason given where it is built.

        With a 'sample' column the two contexts can differ, since the value column only colors a
        single sample's magnitude while the summaries independently choose what the views across
        samples and across groups show. Ungrouped, 'sample_summary' colors the 'unified' map and the
        per-sample maps show each sample on its own; grouped, 'group_summary' colors the 'unified'
        map and 'sample_summary' colors each group's map from its samples. A summary naming an
        aggregation pools values ('quantitative'), while a presence summary or no summary colors
        presence ('membership'), as 'SUMMARY_PRESENCE_SCHEMES' and '_resolve_summary' decide. So a
        value layer can be categorical in the 'unified' map and continuous per sample. A layer
        without a 'sample' column in a run where the other layer has one carries no category values,
        so '_map_elements' holds it constant across the per-sample/group maps.

        Parameters
        ==========
        name : str
            Colorbar filename stem ('reactions'/'compounds').

        data : dict
            The '_read_element_txt' result for this layer.

        path : str
            Path to the layer's text file, used in aggregation error messages.

        Notes
        =====
        'aggregation'/'color'/'colormap'/'colormap_limits'/'category_colormap'/
        'category_colormap_limits'/'category_colors_txt'/'reverse_overlay'/'sample_summary'/
        'group_summary' are this layer's per-layer parameters; 'all_sample_names'/'group_samples'
        carry the shared sample space and grouping. Returns the layer model dict '_map_elements'
        consumes. 'value_limits' and 'category_value_limits' bound the two quantitative scales
        independently — the 'unified' map's and the one its per-sample or per-group maps share —
        since a summary can put the two on quite different scales. Each is refused where its own
        context is not colored by value, and setting one of the two while both contexts are colored
        by value earns a warning, the other scale being an easy one to forget. 'category_colormap'
        colors that same per-sample/per-group scale from a colormap of its own: where the 'unified'
        map shows a summary of a different kind from the values behind it — their spread rather than
        their magnitude, for example — one colormap across both invites reading the two as the same
        quantity. Left unset, both contexts share the layer's 'colormap', and it is refused wherever
        the per-sample/per-group context is not colored by value. 'value_center' and
        'category_value_center' put a value at the middle of those same two scales, each of which
        then runs the same distance either side of it however lopsided its own values are, and are
        refused and warned about on exactly the same footing as the limits.
        """
        element_type = data['element_type']
        use_reaction_attribute = data['reaction_source'] == 'Reaction'
        df = data['df']
        value_column = data['value_column']
        has_sample = data['sample_names'] is not None
        accessions = set(df['__accession'])
        category_colors_flag = f'--{element_type}-category-colors'
        common = {
            'name': name,
            'element_type': element_type,
            'use_reaction_attribute': use_reaction_attribute,
            'accessions': accessions,
            'category_colors_flag': category_colors_flag
        }
        # Colors given per category name are read here, before any of them is checked against the
        # run's categories, so that a malformed file is reported as such rather than as a file full
        # of unrecognized names. '_map_elements' does that check, where the categories are known.
        if category_colors_txt is None:
            common['category_colors'] = None
            common['category_combo_colors'] = {}
        else:
            category_colors, combo_colors = self._read_category_colors_txt(
                category_colors_txt, category_colors_flag
            )
            common['category_colors'] = category_colors
            common['category_combo_colors'] = combo_colors

        # Limits and centers are checked before any value is read, so that a pair or a center anvi'o
        # cannot make sense of is reported as the configuration error it is rather than surfacing
        # later as a strange scale. What each can shape at all is settled here too; which context is
        # actually colored by value is not known until the summaries below have chosen the modes.
        # The centers are resolved against the limits, which is where a center placed out of a
        # limited scale's reach is caught, so the limits are resolved first.
        value_limits_flag = f'--{element_type}-value-limits'
        category_value_limits_flag = f'--{element_type}-category-value-limits'
        value_center_flag = f'--{element_type}-value-center'
        category_value_center_flag = f'--{element_type}-category-value-center'
        value_limits = self._resolve_value_limits(value_limits, value_limits_flag)
        category_value_limits = self._resolve_value_limits(
            category_value_limits, category_value_limits_flag
        )
        value_center = self._resolve_value_center(
            value_center, value_center_flag, value_limits, value_limits_flag
        )
        category_value_center = self._resolve_value_center(
            category_value_center, category_value_center_flag, category_value_limits,
            category_value_limits_flag
        )
        if value_column is None:
            for setting, flag, verb in (
                (value_limits, value_limits_flag, 'bound'),
                (category_value_limits, category_value_limits_flag, 'bound'),
                (value_center, value_center_flag, 'center'),
                (category_value_center, category_value_center_flag, 'center')
            ):
                if setting is not None:
                    raise ConfigError(
                        f"'{flag}' {verb}s the color scale spanning the values of the "
                        f"{element_type} layer, but the file at '{path}' has no value column, so "
                        f"that layer is colored by presence and has no scale of values to {verb}."
                    )
        if not has_sample:
            for setting, flag, single_scale_flag, verb in (
                (category_value_limits, category_value_limits_flag, value_limits_flag, 'bound'),
                (category_value_center, category_value_center_flag, value_center_flag, 'center')
            ):
                if setting is not None:
                    raise ConfigError(
                        f"'{flag}' {verb}s the color scale shared by the maps of the individual "
                        f"samples or groups of the {element_type} layer, but the file at '{path}' "
                        f"has no 'sample' column, so it draws no such maps and its values take a "
                        f"single scale. {verb.capitalize()} that one with '{single_scale_flag}'."
                    )

        # The colormap of the per-sample/per-group scale is settled on the same footing as the limits
        # on it: what it can color at all is checked here, before any value is read, while whether
        # that context is colored by value at all waits on the summaries below.
        category_colormap_flag = f'--{element_type}-category-colormap'
        colormap_flag = f'--{element_type}-colormap'
        if category_colormap is not None:
            if value_column is None:
                raise ConfigError(
                    f"'{category_colormap_flag}' colors the maps of the individual samples or "
                    f"groups of the {element_type} layer by value, but the file at '{path}' has no "
                    f"value column, so that layer is colored by presence rather than along a scale "
                    f"of values. '{colormap_flag}' is what colors a presence scale."
                )
            if not has_sample:
                raise ConfigError(
                    f"'{category_colormap_flag}' colors the maps of the individual samples or "
                    f"groups of the {element_type} layer, but the file at '{path}' has no 'sample' "
                    f"column, so it draws no such maps and its values take a single scale. Color "
                    f"that one with '{colormap_flag}'."
                )

        if color == 'original':
            # A reaction presence layer highlighted in the reference map's colors
            # ('--original-color'), drawn by the reference-color drawer rather than the element
            # engine, since both the colors and their render order come from the reference map, so
            # there is no channel for data. The CLI restricts this to a reaction file with no value
            # column, magnitude having nowhere to go. With a 'sample' column, the 'unified' map is
            # the union of the samples and each sample also gets its own map; without one,
            # 'membership' just carries the accession set the single map colors. 'category_mode' is
            # 'membership' rather than 'original' so that a grouped run takes the same path as the
            # database and pangenome inputs, whose per-group maps show within-group source counts:
            # 'source_accessions' is keyed by sample, so the per-category reference-color branch can
            # only serve samples, not groups. The reference-color drawer finds elements by the KO
            # IDs of a map's ortholog entries, so a file of KEGG reaction IDs would match nothing
            # and draw a blank map, and there is no compound layer to draw at all.
            if element_type != 'reaction':
                raise ConfigError(
                    f"The reference map's own colors highlight reaction elements, so they apply to "
                    f"a reaction layer. The file at '{path}' is a {element_type} file. Use "
                    f"'--original-color' with a reaction file, and color compounds with "
                    f"'--compound-color' or '--compound-colormap'."
                )
            if use_reaction_attribute:
                raise ConfigError(
                    f"The reference map's own colors are applied by matching the KO IDs of a map's "
                    f"reaction elements, but the accessions in the reaction file at '{path}' are "
                    f"KEGG reaction IDs ('R' followed by digits) rather than KO IDs, so nothing "
                    f"would be highlighted. Use a file of KO accessions with '--original-color', "
                    f"or color these reactions with '--reaction-color' or '--reaction-colormap'."
                )
            if common['category_colors'] is not None:
                raise ConfigError(
                    f"'{category_colors_flag}' gives a color to each sample or group, but "
                    f"'--original-color' takes both its colors and its drawing order from the "
                    f"reference map, so there is no color left to choose. Please use only one of "
                    f"them."
                )
            for setting, flag, verb in (
                (value_limits, value_limits_flag, 'bound'),
                (category_value_limits, category_value_limits_flag, 'bound'),
                (value_center, value_center_flag, 'center'),
                (category_value_center, category_value_center_flag, 'center')
            ):
                if setting is not None:
                    raise ConfigError(
                        f"'{flag}' {verb}s a color scale spanning a layer's values, but "
                        f"'--original-color' takes both its colors and its drawing order from the "
                        f"reference map, so nothing is colored by value and there is no scale to "
                        f"{verb}. Please use only one of them."
                    )
            if category_colormap is not None:
                raise ConfigError(
                    f"'{category_colormap_flag}' colors the maps of the individual samples or "
                    f"groups along a scale of values, but '--original-color' takes both its colors "
                    f"and its drawing order from the reference map, so nothing is colored by value "
                    f"and there is no scale to color. Please use only one of them."
                )
            if has_sample:
                membership, source_accessions = self._relate_accessions_to_samples(
                    df, all_sample_names
                )
            else:
                membership = {accession: [] for accession in accessions}
                source_accessions = {}
            return {
                **common,
                'unified_mode': 'original',
                'category_mode': 'membership',
                'membership': membership,
                'source_accessions': source_accessions,
                'color_hexcode': 'original'
            }

        undefined: Set[str] = set()
        # Resolved here, before any values are aggregated, so that an unusable aggregation name is
        # reported as a configuration error rather than reaching pandas: the per-accession values
        # below are computed by passing the name itself to a groupby. 'aggregate' reduces a map
        # element's several accessions, which is the level the drawing colorers work at; the gene
        # aggregation is applied only where the per-accession values are built.
        aggregate = self._resolve_aggregation(
            accession_aggregation, f'--{element_type}-accession-aggregation'
        ) if value_column is not None else None
        if value_column is not None:
            self._resolve_aggregation(gene_aggregation, f'--{element_type}-gene-aggregation')

        if not has_sample:
            # With no samples there is nothing to summarize, so the layer colors the same way in
            # every map context: by its value column if it has one, otherwise a single presence
            # color.
            if value_column is None:
                return {
                    **common,
                    'unified_mode': 'single',
                    'category_mode': 'single',
                    'color_hexcode': color
                }
            unified_values = self._aggregate_accession_quantities(
                df, value_column, gene_aggregation, path, undefined
            )
            self._warn_undefined_values(undefined, gene_aggregation, path)
            cmap = self._resolve_sequential_colormap(
                colormap if colormap is not None else 'plasma_r', colormap_limits,
                subject=f'{element_type}s'
            )
            self._check_centered_colormap(
                cmap, colormap_limits, value_center, value_center_flag, colormap_flag,
                subject=f'{element_type}s'
            )
            return {
                **common,
                'unified_mode': 'quantitative',
                'category_mode': 'quantitative',
                'cmap': cmap,
                # There are no per-sample maps to give a colormap of their own — a category colormap
                # was refused above — so the one map's scale serves both contexts.
                'category_cmap': cmap,
                'reverse_overlay': reverse_overlay,
                'unified_values': unified_values,
                'category_values': None,
                'aggregate': aggregate,
                'colorbar_label': value_column,
                'value_limits': value_limits,
                'category_value_limits': None,
                'value_center': value_center,
                'category_value_center': None
            }

        # Resolve the two summaries into the mode of each map context. Ungrouped, the per-sample
        # maps show one sample each, so no summary applies to them: they are the sample's own
        # magnitude with a value column and its presence without one.
        grouped = group_samples is not None
        sample_kind, sample_aggregate, sample_scheme = self._resolve_summary(
            sample_summary, value_column, f'--{element_type}-sample-summary', path
        )
        if grouped:
            group_kind, group_aggregate, group_scheme = self._resolve_summary(
                group_summary, value_column, f'--{element_type}-group-summary', path
            )
            if group_kind == 'value' and sample_kind != 'value':
                raise ConfigError(
                    f"'--{element_type}-group-summary' was given as '{group_summary}', which pools "
                    f"the values of the sample groups, but '--{element_type}-sample-summary' "
                    f"summarizes each group's samples by presence rather than by value, so the "
                    f"groups have no values to pool. Please also set "
                    f"'--{element_type}-sample-summary' to an aggregation such as 'mean'."
                )
            unified_kind, unified_scheme = group_kind, group_scheme
            category_kind = sample_kind
        else:
            unified_kind, unified_scheme = sample_kind, sample_scheme
            category_kind = 'presence' if value_column is None else 'value'

        membership, source_accessions = self._relate_accessions_to_samples(df, all_sample_names)

        # A static color chosen explicitly overrides comparison across samples on the 'unified' map,
        # which then shows presence in ANY sample in that one color, exactly as it does for multiple
        # contigs databases or a pangenome. The CLI signals the choice by passing a colormap of
        # False, the same way it does for those inputs, and rejects the summary options alongside
        # it. The individual maps are unaffected: each still shows one sample, or, grouped, its
        # group's sample counts.
        static_color = colormap is False
        unified_mode = 'static' if static_color else (
            'quantitative' if unified_kind == 'value' else 'membership'
        )

        # The option that chose the 'unified' map's presence scheme, so that a message about that
        # scheme names the option that can change it rather than the one this input refuses: the
        # group summary colors the 'unified' map of a grouped run, and the sample summary colors it
        # otherwise, exactly as 'unified_scheme' was taken above.
        summary_flag = f"--{element_type}-{'group' if grouped else 'sample'}-summary"

        model = {
            **common,
            'unified_mode': unified_mode,
            'category_mode': 'quantitative' if category_kind == 'value' else 'membership',
            'membership': membership,
            'source_accessions': source_accessions,
            'color_hexcode': color,
            'colormap': True if colormap is None else colormap,
            'colormap_limits': colormap_limits,
            'colormap_scheme': unified_scheme,
            'scheme_options': {
                scheme: f'{summary_flag} {name}'
                for name, scheme in SUMMARY_PRESENCE_SCHEMES.items()
            },
            'reverse_overlay': reverse_overlay,
            'category_values': None,
            'value_limits': value_limits,
            'category_value_limits': category_value_limits,
            'value_center': value_center,
            'category_value_center': category_value_center
        }

        if value_column is None:
            return model

        # A limit bounds, and a center centers, a scale that colors by value, so a context colored
        # by presence, or in one static color, offers nothing for either to act on. Which of the two
        # it is decides what the user should reach for instead, so each case names what put that
        # context where it is.
        for setting, flag, other_flag, verb in (
            (value_limits, value_limits_flag, category_value_limits_flag, 'bound'),
            (value_center, value_center_flag, category_value_center_flag, 'center')
        ):
            if setting is not None and model['unified_mode'] != 'quantitative':
                reason = (
                    f"a single static color was chosen for it with '--{element_type}-color'"
                    if model['unified_mode'] == 'static' else
                    f"'{summary_flag}' summarizes presence rather than pooling values"
                )
                raise ConfigError(
                    f"'{flag}' {verb}s the color scale of values on the 'unified' map, but that "
                    f"map is not colored by value here: {reason}. Color it by value, or {verb} the "
                    f"scale of the maps of the individual samples or groups with '{other_flag}' "
                    f"instead."
                )
        for setting, flag, other_flag, verb in (
            (category_value_limits, category_value_limits_flag, value_limits_flag, 'bound'),
            (category_value_center, category_value_center_flag, value_center_flag, 'center')
        ):
            if setting is not None and model['category_mode'] != 'quantitative':
                raise ConfigError(
                    f"'{flag}' {verb}s the color scale shared by the maps of the individual "
                    f"groups, but those maps are not colored by value here: "
                    f"'--{element_type}-sample-summary' summarizes each group's samples by "
                    f"presence rather than by pooling their values. Set it to an aggregation such "
                    f"as 'mean' to color the group maps by value, or {verb} the 'unified' map's "
                    f"own scale with '{other_flag}'."
                )
        if category_colormap is not None and model['category_mode'] != 'quantitative':
            raise ConfigError(
                f"'{category_colormap_flag}' colors the scale shared by the maps of the individual "
                f"groups, but those maps are not colored along a scale of values here: "
                f"'--{element_type}-sample-summary' summarizes each group's samples by presence "
                f"rather than by pooling their values, and presence there is colored by "
                f"'--group-colormap'. Set the sample summary to an aggregation such as 'mean' to "
                f"color the group maps by value, or color the 'unified' map's own scale with "
                f"'{colormap_flag}'."
            )

        # The two scales are bounded separately because a summary can put them on quite different
        # footings: a mean across samples spans a good deal less than the samples themselves do. The
        # flip side is that a limit on one says nothing about the other, and the scale left alone is
        # an easy one to overlook, so bounding exactly one of two scales that both color by value is
        # worth saying out loud. It is a legitimate thing to ask for, so this is not an error.
        if (
            model['unified_mode'] == 'quantitative' and model['category_mode'] == 'quantitative'
            and (value_limits is None) != (category_value_limits is None)
        ):
            if value_limits is None:
                given_flag, other_flag = category_value_limits_flag, value_limits_flag
                other_subject = "the 'unified' map"
            else:
                given_flag, other_flag = value_limits_flag, category_value_limits_flag
                other_subject = f"the maps of the individual {'groups' if grouped else 'samples'}"
            self.run.warning(
                f"'{given_flag}' bounds the color scale of the {element_type} layer, but "
                f"'{other_flag}' was not given, so the scale of {other_subject} still spans "
                f"whatever its own values happen to reach. The two are bounded separately because "
                f"they can differ a great deal -- a summary such as a mean across samples spans "
                f"less than the samples themselves do -- so this may well be what you intend. If "
                f"it is not, give '{other_flag}' limits of its own."
            )

        # Centering exactly one of two scales that both color by value is worth saying out loud for
        # the same reason, and for one more: unless a colormap was given to the per-sample/per-group
        # scale of its own, the two are drawn from the same colormap, and then the middle color
        # means the centered value on the one map and whatever the values happen to leave in the
        # middle on the other.
        if (
            model['unified_mode'] == 'quantitative' and model['category_mode'] == 'quantitative'
            and (value_center is None) != (category_value_center is None)
        ):
            if value_center is None:
                given_flag, other_flag = category_value_center_flag, value_center_flag
                given_center = category_value_center
                other_subject = "the 'unified' map"
            else:
                given_flag, other_flag = value_center_flag, category_value_center_flag
                given_center = value_center
                other_subject = f"the maps of the individual {'groups' if grouped else 'samples'}"
            self.run.warning(
                f"'{given_flag}' centers the color scale of the {element_type} layer on "
                f"{given_center:g}, but '{other_flag}' was not given, so the scale of "
                f"{other_subject} still sits wherever its own values leave it. Where one colormap "
                f"colors both scales, which is what happens unless "
                f"'{category_colormap_flag}' gives one of them a colormap of its own, the very "
                f"same color then means the centered value on the one map and something else "
                f"entirely on the other. The two are centered separately because a summary can put "
                f"them on quite different footings, so this may well be what you intend. If it is "
                f"not, center the other with '{other_flag}'."
            )

        # Per-sample values: each sample's rows reduced to one value per accession by the
        # within-sample aggregation. These are computed even when no context is colored by value, so
        # that a value column is validated whenever the file has one.
        sample_values: Dict[str, Dict[str, float]] = {s: {} for s in all_sample_names}
        for sample_name, sample_rows in df.groupby('__sample'):
            sample_values[sample_name] = self._aggregate_accession_quantities(
                sample_rows, value_column, gene_aggregation, path, undefined
            )
        self._warn_undefined_values(undefined, gene_aggregation, path)

        if 'quantitative' not in (model['unified_mode'], model['category_mode']):
            # Both summaries color presence, which only happens with groups (the per-sample maps of
            # an ungrouped run always show a value column's magnitude), so nothing is colored by
            # value.
            self.run.warning(
                f"The layer from the text file at '{path}' has a value column, "
                f"'{value_column}', but nothing on the maps is colored by it: with sample groups, "
                f"the per-group maps are colored by '--{element_type}-sample-summary' and the "
                f"'unified' map by '--{element_type}-group-summary', and both of these summarize "
                f"presence rather than value. Set '--{element_type}-sample-summary' to an "
                f"aggregation, such as 'mean', to color the group maps by value."
            )
            return model

        # Each summary is a further reduction that can be undefined where the one below it was not,
        # so the accessions it drops are reported against the summary's own name rather than the
        # within-sample aggregation's.
        if category_kind == 'value':
            if grouped:
                undefined_category: Set[str] = set()
                model['category_values'] = {
                    group_name: self._summarize_category_values(
                        {s: sample_values[s] for s in group_sample_names}, sample_aggregate,
                        undefined_category
                    )
                    for group_name, group_sample_names in group_samples.items()
                }
                self._warn_undefined_values(undefined_category, sample_summary, path)
            else:
                # Each per-sample map shows one sample, so no summary applies to it.
                model['category_values'] = sample_values

        if unified_kind == 'value':
            # The 'unified' map summarizes the groups when there are groups, each of which
            # contributes the summary of its own samples, and otherwise all of the samples.
            undefined_unified: Set[str] = set()
            model['unified_values'] = self._summarize_category_values(
                model['category_values'] if grouped else sample_values,
                group_aggregate if grouped else sample_aggregate,
                undefined_unified
            )
            self._warn_undefined_values(
                undefined_unified, group_summary if grouped else sample_summary, path
            )

        model['cmap'] = self._resolve_sequential_colormap(
            colormap if colormap is not None else 'plasma_r', colormap_limits,
            subject=f'{element_type}s'
        )
        # The maps of the individual samples or groups take a colormap of their own where one was
        # given, so that a 'unified' map showing a summary of another kind than the values behind it
        # — their spread rather than their magnitude — is not read as more of the same quantity.
        # Left unset, the layer's one colormap serves both contexts, and is resolved once so that a
        # warning about it is given once.
        category_subject = (
            f"{element_type}s on the maps of the individual "
            f"{'groups' if grouped else 'samples'}"
        )
        model['category_cmap'] = model['cmap'] if category_colormap is None else (
            self._resolve_sequential_colormap(
                category_colormap, category_colormap_limits, subject=category_subject
            )
        )
        # Whether a centered scale's colormap has a middle worth putting a value at is asked of each
        # scale's own colormap, and only where that scale both colors by value and was centered. The
        # two questions are asked separately even where one colormap serves both, since a scale that
        # was not centered has nothing to ask.
        if model['unified_mode'] == 'quantitative':
            self._check_centered_colormap(
                model['cmap'], colormap_limits, value_center, value_center_flag, colormap_flag,
                subject=f"{element_type}s on the 'unified' map"
            )
        if model['category_mode'] == 'quantitative':
            self._check_centered_colormap(
                model['category_cmap'],
                colormap_limits if category_colormap is None else category_colormap_limits,
                category_value_center, category_value_center_flag,
                colormap_flag if category_colormap is None else category_colormap_flag,
                subject=category_subject
            )
        model['aggregate'] = aggregate
        model['colorbar_label'] = value_column
        return model

    def map_kegg_pathways_txt(
        self,
        output_dir: str,
        reaction_txt: str = None,
        compound_txt: str = None,
        reaction_gene_aggregation: str = 'sum',
        reaction_accession_aggregation: str = 'sum',
        compound_accession_aggregation: str = 'sum',
        reaction_sample_summary: str = None,
        compound_sample_summary: str = None,
        reaction_group_summary: str = None,
        compound_group_summary: str = None,
        groups_txt: str = None,
        group_threshold: float = None,
        pathway_numbers: Iterable[str] = None,
        draw_individual_files: Union[Iterable[str], bool] = False,
        draw_grid: Union[Iterable[str], bool] = False,
        reaction_color: str = '#2ca02c',
        compound_color: str = "#e239af",
        reaction_colormap: Union[bool, str, mcolors.Colormap] = None,
        reaction_colormap_limits: Tuple[float, float] = None,
        reaction_category_colormap: Union[str, mcolors.Colormap] = None,
        reaction_category_colormap_limits: Tuple[float, float] = None,
        reaction_category_colors: str = None,
        reaction_reverse_overlay: bool = False,
        reaction_value_limits: Tuple[Union[float, None], Union[float, None]] = None,
        reaction_category_value_limits: Tuple[Union[float, None], Union[float, None]] = None,
        reaction_value_center: Union[float, None] = None,
        reaction_category_value_center: Union[float, None] = None,
        compound_colormap: Union[bool, str, mcolors.Colormap] = None,
        compound_colormap_limits: Tuple[float, float] = None,
        compound_category_colormap: Union[str, mcolors.Colormap] = None,
        compound_category_colormap_limits: Tuple[float, float] = None,
        compound_category_colors: str = None,
        compound_reverse_overlay: bool = False,
        compound_value_limits: Tuple[Union[float, None], Union[float, None]] = None,
        compound_category_value_limits: Tuple[Union[float, None], Union[float, None]] = None,
        compound_value_center: Union[float, None] = None,
        compound_category_value_center: Union[float, None] = None,
        group_colormap: Union[str, mcolors.Colormap] = 'plasma_r',
        group_colormap_limits: Tuple[float, float] = None,
        group_reverse_overlay: bool = False,
        group_colormap_scheme: Literal['by_count', 'by_count_continuous'] = None,
        count_scale_max: Union[str, int] = 'observed',
        draw_maps_lacking_data: bool = False
    ) -> Dict[Literal['unified', 'individual', 'grid'], Dict]:
        """
        Draw pathway maps from a reaction-layer file and/or a compound-layer file.

        The reaction file (kegg-reaction-txt, '--reaction-txt') colors reaction elements; the
        compound file (kegg-compound-txt, '--compound-txt') colors compound elements. Either or both
        may be given, and each layer is colored independently by what its file provides: a value
        column colors it quantitatively (continuous colorbar), a 'sample' column without a value
        column colors it by sample/group presence, and neither colors it a single presence color.
        The two layers are drawn together on each map, so one map can mix, e.g., presence of KOs
        with quantitative compound values. When any layer has a 'sample' column, a 'unified' map
        summarizes the samples and per-sample maps are added; a 'groups_txt' instead draws per-group
        maps and the 'unified' map summarizes the groups.

        Parameters
        ==========
        output_dir : str
            Path to the output directory in which pathway map and colorbar PDF files are drawn.

        reaction_txt : str, None
            Path to a kegg-reaction-txt file (accessions all KO or all KEGG reaction IDs).

        compound_txt : str, None
            Path to a kegg-compound-txt file (accessions all KEGG compound IDs).

        reaction_gene_aggregation : str, 'sum'
            How to reduce the values of the genes annotated with one accession to that accession's
            value, for a reaction file carrying a 'gene_id' column. Any 'AGGREGATION_FUNCTIONS' name,
            or any other pandas aggregation reducing values to one number (see
            '_resolve_aggregation').

        reaction_accession_aggregation : str, 'sum'
            How to reduce the values of the several accessions a map's reaction element stands to
            that element's value.

        compound_accession_aggregation : str, 'sum'
            The same reduction for the compound layer: the several compounds that one map circle
            stands for. A compound file has no genes, and repeated rows are refused
            ('_read_element_txt'), so it has no reduction below this one.

        reaction_sample_summary : str, None
            How the reaction layer summarizes a set of samples: by presence ('count'/'membership')
            or by pooling their values with an aggregation. This colors the 'unified' map without
            groups and each per-group map with groups. The default of None summarizes presence, by
            membership for 3 or fewer samples and by count above that.

        compound_sample_summary : str, None
            The same summary for the compound layer's samples.

        reaction_group_summary : str, None
            How the reaction layer summarizes the sample groups of 'groups_txt', which colors the
            'unified' map when there are groups. The default of None summarizes presence, by
            membership for 3 or fewer groups and by count above that.

        compound_group_summary : str, None
            The same summary for the compound layer's groups.

        reaction_category_colors : str, None
            Path to a kegg-category-colors-txt file giving the reaction layer a color per category —
            per sample, or per group when the samples are grouped — which colors it by membership in
            place of a colormap, and colors each category's own map. A row naming a combination of
            categories overrides the blend of their colors ('_read_category_colors_txt').

        compound_category_colors : str, None
            The same colors for the compound layer.

        reaction_category_colormap : Union[str, matplotlib.colors.Colormap], None
            The colormap of the single color scale shared by the reaction layer's per-sample or
            per-group maps, which 'reaction_category_colormap_limits' trims as
            'reaction_colormap_limits' trims 'reaction_colormap'. The default of None gives both map
            contexts the layer's own 'reaction_colormap', and a colormap here separates them, so
            that a 'unified' map summarizing the samples by a statistic of another kind than the
            values themselves — their standard deviation, say — is not drawn in the colors that
            stand for magnitude on the per-sample maps.

        compound_category_colormap : Union[str, matplotlib.colors.Colormap], None
            The same colormap for the compound layer's per-sample or per-group scale.

        reaction_value_limits : Tuple[Union[float, None], Union[float, None]], None
            The (minimum, maximum) the reaction layer's color scale may span on the 'unified' map,
            either end of which can be None to leave it wherever the values put it. A limit
            truncates the scale only where the values cross it, and the colorbar then marks that end
            '<=' or '>=' ('_make_quantitative_norm').

        reaction_category_value_limits : Tuple[Union[float, None], Union[float, None]], None
            The same limits for the scale shared by the reaction layer's per-sample or per-group
            maps, which is bounded separately because a summary can put it and the 'unified' map's
            scale on quite different footings.

        compound_value_limits : Tuple[Union[float, None], Union[float, None]], None
            The same limits for the compound layer's 'unified' map scale.

        compound_category_value_limits : Tuple[Union[float, None], Union[float, None]], None
            The same limits for the compound layer's per-sample or per-group scale.

        reaction_value_center : Union[float, None], None
            The value put at the middle of the reaction layer's color scale on the 'unified' map.
            The scale then runs the same distance either side of it, so that the middle color of a
            diverging colormap stands for this value however lopsided the values around it are
            ('_make_quantitative_norm'). The default of None leaves the scale spanning its own
            values from end to end.

        reaction_category_value_center : Union[float, None], None
            The same center for the scale shared by the reaction layer's per-sample or per-group
            maps, which is centered separately because a summary can put it and the 'unified' map's
            scale on quite different footings.

        compound_value_center : Union[float, None], None
            The same center for the compound layer's 'unified' map scale.

        compound_category_value_center : Union[float, None], None
            The same center for the compound layer's per-sample or per-group scale.

        Notes
        =====
        The remaining parameters carry the shared groups, per-layer colors/colormaps, and drawing
        options; see the CLI help and '_map_elements'. 'group_colormap' may be
        'GROUP_COLORMAP_FROM_CATEGORY' to color each group's own maps by a ramp running to that
        group's own color instead of by a named colormap, in which case 'group_colormap_limits' is
        how far from white that ramp runs ('_group_map_colors'). 'group_colormap_scheme' draws the
        count scale of every group's maps in discrete bands ('by_count') or as a gradient
        ('by_count_continuous'), and by default keeps the bands while the colors can be told apart.

        Returns
        =======
        Dict[Literal['unified', 'individual', 'grid'], Dict]
            The record returned by '_map_elements'.
        """
        raw_layers: List[dict] = []
        if reaction_txt is not None:
            raw_layers.append({
                'name': 'reactions',
                'path': reaction_txt,
                'data': self._read_element_txt(reaction_txt, 'reaction'),
                'gene_aggregation': reaction_gene_aggregation,
                'accession_aggregation': reaction_accession_aggregation,
                'color': reaction_color,
                'colormap': reaction_colormap,
                'colormap_limits': reaction_colormap_limits,
                'category_colormap': reaction_category_colormap,
                'category_colormap_limits': reaction_category_colormap_limits,
                'category_colors_txt': reaction_category_colors,
                'reverse_overlay': reaction_reverse_overlay,
                'sample_summary': reaction_sample_summary,
                'group_summary': reaction_group_summary,
                'value_limits': reaction_value_limits,
                'category_value_limits': reaction_category_value_limits,
                'value_center': reaction_value_center,
                'category_value_center': reaction_category_value_center
            })
        if compound_txt is not None:
            raw_layers.append({
                'name': 'compounds',
                'path': compound_txt,
                'data': self._read_element_txt(compound_txt, 'compound'),
                # A compound file has one row per accession per sample, so the gene level is a
                # reduction over a single value and its function cannot matter.
                'gene_aggregation': 'sum',
                'accession_aggregation': compound_accession_aggregation,
                'color': compound_color,
                'colormap': compound_colormap,
                'colormap_limits': compound_colormap_limits,
                'category_colormap': compound_category_colormap,
                'category_colormap_limits': compound_category_colormap_limits,
                'category_colors_txt': compound_category_colors,
                'reverse_overlay': compound_reverse_overlay,
                'sample_summary': compound_sample_summary,
                'group_summary': compound_group_summary,
                'value_limits': compound_value_limits,
                'category_value_limits': compound_category_value_limits,
                'value_center': compound_value_center,
                'category_value_center': compound_category_value_center
            })
        if not raw_layers:
            raise ConfigError(
                "No draw-kegg-pathways input files were provided. Supply a reaction-layer file "
                "('--reaction-txt'), a compound-layer file ('--compound-txt'), or both."
            )

        # One shared sample space across both files: when both carry a 'sample' column they describe
        # the same samples, so the sample names are unioned; a layer with no 'sample' column is held
        # constant across the per-sample/group maps.
        sample_name_sets = [
            layer['data']['sample_names'] for layer in raw_layers
            if layer['data']['sample_names'] is not None
        ]
        all_sample_names = sorted(set().union(*sample_name_sets)) if sample_name_sets else None

        # A threshold outside [0, 1] cannot be met by any proportion of a group's samples, or is met
        # by all of them, so it would silently draw an empty or an undiscriminating map. The
        # database and pangenome methods make the same check.
        if group_threshold is not None and not 0 <= group_threshold <= 1:
            raise ConfigError(
                f"'group_threshold' must be a number between 0 and 1, not {group_threshold}. It is "
                f"the proportion of a group's samples in which an element must occur for the group "
                f"to be considered to contain it."
            )

        group_samples = None
        sample_group = None
        categories = None
        category_noun = None
        if all_sample_names is not None:
            if groups_txt is not None:
                sample_group, group_samples = self._relate_samples_to_groups(
                    groups_txt, all_sample_names
                )
                categories = list(group_samples)
                category_noun = 'group'
            else:
                categories = all_sample_names
                category_noun = 'sample'
            self.run.info("Samples found across the input file(s)", len(all_sample_names))

        models = [
            self._build_txt_model(
                **layer, all_sample_names=all_sample_names, group_samples=group_samples
            )
            for layer in raw_layers
        ]

        # Groups are carried through whenever they are defined: they set which samples each group's
        # map summarizes, and, for a group summary of presence, which groups the 'unified' map shows
        # an element in (per 'group_threshold').
        grouped_membership = None
        if group_samples is not None:
            grouped_membership = {
                'source_group': sample_group,
                'group_sources': group_samples,
                'group_threshold': group_threshold,
                'group_colormap': group_colormap,
                'group_colormap_limits': group_colormap_limits,
                'group_reverse_overlay': group_reverse_overlay,
                'group_colormap_scheme': group_colormap_scheme
            }

        category_subject = f"{category_noun}s" if category_noun is not None else None
        return self._map_elements(
            models,
            output_dir,
            pathway_numbers=pathway_numbers,
            categories=categories,
            category_noun=category_noun,
            grid_source_type='sample',
            colorbar_category_suffix=category_subject,
            subset_subject=category_subject,
            unified_plural='samples',
            membership_count_label='sample count',
            membership_members_label='samples',
            membership_singular='sample',
            grouped_membership=grouped_membership,
            count_scale_max=count_scale_max,
            draw_individual_files=draw_individual_files,
            draw_grid=draw_grid,
            draw_maps_lacking_data=draw_maps_lacking_data
        )

    def map_contigs_databases_kos(
        self,
        contigs_dbs: Iterable[str],
        output_dir: str,
        groups_txt: str = None,
        group_threshold: float = None,
        pathway_numbers: Iterable[str] = None,
        draw_individual_files: Union[Iterable[str], bool] = False,
        draw_grid: Union[Iterable[str], bool] = False,
        reaction_colormap: Union[bool, str, mcolors.Colormap] = True,
        reaction_colormap_limits: Tuple[float, float] = None,
        colormap_scheme: Literal['by_count', 'by_count_continuous', 'by_membership'] = None,
        reaction_category_colors: str = None,
        reaction_reverse_overlay: bool = False,
        reaction_color: str = '#2ca02c',
        group_colormap: Union[str, mcolors.Colormap] = 'plasma_r',
        group_colormap_limits: Tuple[float, float] = None,
        group_reverse_overlay: bool = False,
        group_colormap_scheme: Literal['by_count', 'by_count_continuous'] = None,
        count_scale_max: Union[str, int] = 'observed',
        draw_maps_lacking_data: bool = False
    ) -> Dict[Literal['unified', 'individual', 'grid'], Dict]:
        """
        Draw pathway maps, coloring the reaction layer by KOs across contigs databases
        (representing, for example, genomes or metagenomes) or groups of databases (representing,
        for example, taxonomic groups of genomes or geographical groups of metagenomes).

        A reaction element on a map is represented by one or more KOs, matched to the KO annotations
        of each contigs database; a database "contains" the reaction if it has any of those KOs. The
        reaction elements (lines in global and overview maps, boxes or lines in standard maps) are
        colored by the databases or groups containing them via '_map_element_membership'.

        Parameters
        ==========
        contigs_dbs : Iterable[str]
            File paths to contigs databases containing KO annotations. Databases should have
            different project names, by which they are uniquely identified.

        output_dir : str
            Path to the output directory in which pathway map and colorbar PDF files are drawn. The
            directory is created if it does not exist.

        groups_txt : str, None
            A tab-delimited text file specifying which group each contigs database belongs to. The
            first column, which can have any header, contains the file paths of contigs databases,
            those provided to the 'contigs_dbs' argument. The second column, which must be headed
            'group', contains group names, which are recommended to be single words without fancy
            characters, such as 'HIGH_TEMPERATURE' or 'LOW_FITNESS' rather than 'my group #1' or
            'IS-THIS-OK?'. Each contigs database can only be associated with a single group. The
            'group_threshold' argument must also be used for the groups to take effect, assigning
            colors based on group membership and drawing individual files ('draw_individual_files')
            and map grids ('draw_grid') for groups rather than individual databases.

        group_threshold : float, None
            The proportion of contigs databases in a group containing data of interest for the group
            to be represented in terms of presence/absence in a reaction element. Here is a concrete
            example. Say each contigs database represents a genome, and the 'groups_txt' argument,
            which must be used with this argument, groups these genomes by their species, 'A', 'B',
            and 'C'. You wish to understand the distribution of metabolic capabilities across the 3
            species from KO annotations of genes. Reaction colors are assigned based on the groups
            rather than individual genomic contigs databases containing the reaction. Thresholds
            between 0 and 1 can be set to define group membership: a threshold of 0.0 would mean
            that ANY genome in the group can contain the reaction via KOs for the reaction to be
            considered present in the group; a threshold of 0.75 means at least 75% of the genomes
            in the group must contain the reaction for it to be present; a threshold of 1.0 means
            that ALL genomes in the group must contain the reaction for it to be present. In our
            example, set the threshold to 0.5. Reaction J on a map corresponds to KO X, and Reaction
            K on a map corresponds to KOs Y and Z. 90% of species A genomes, 50% of species B
            genomes, and 10% of species C genomes contain KO X, so Reaction J would be colored to
            indicate that it is represented in species A and B. 0% of species A genomes, 15% of
            species B genomes, and 40% of species C genomes contain KO Y and KO Z, so Reaction K
            would not be colored.

        pathway_numbers : Iterable[str], None
            Regex patterns to match the ID numbers of the drawn pathway maps. The default of None
            draws all available pathway maps in the KEGG data directory.

        reaction_color : str, '#2ca02c'
            The single color, by default green, for reaction elements when dynamic coloring is
            disabled (when 'reaction_colormap' is False). It colors both the unified map (by
            presence/absence in any database) and the individual-database maps. Alternatively, the
            string 'original' uses the reference map's original color scheme.

        draw_maps_lacking_data : bool, False
            If False, by default, only draw maps containing any of the select KOs. If True, draw
            maps regardless, meaning that nothing may be colored.

        Notes
        =====
        The dynamic-coloring and drawing options ('reaction_colormap', 'reaction_colormap_limits',
        'colormap_scheme', 'reaction_category_colors', 'reaction_reverse_overlay',
        'draw_individual_files', 'draw_grid', and the group colormap options) mirror the categorical
        engine; see the CLI help and '_map_element_membership'. 'reaction_category_colors' is the
        path to a kegg-category-colors-txt file giving a color per category — per source, or per
        group when the sources are grouped — which colors the layer by membership in place of a
        colormap, and colors each category's own map. 'group_colormap' may be
        'GROUP_COLORMAP_FROM_CATEGORY' to color each group's own maps by a ramp running to that
        group's own color instead of by a named colormap, in which case 'group_colormap_limits' is
        how far from white that ramp runs ('_group_map_colors'). 'group_colormap_scheme' draws the
        count scale of every group's maps in discrete bands ('by_count') or as a gradient
        ('by_count_continuous'), and by default keeps the bands while the colors can be told apart.

        Returns
        =======
        Dict[Literal['unified', 'individual', 'grid'], Dict]
            The record returned by '_map_element_membership': 'unified' maps show all contigs
            databases or groups, 'individual' maps show single databases or groups, and 'grid'
            images show both. See '_map_element_membership' for the nested structure.
        """
        # This method loads KO membership from contigs databases and hands off the drawing of
        # unified, individual, and grid maps to '_map_element_membership'.

        self.progress.new("Loading metadata from contigs databases")
        self.progress.update("...")

        self._check_contigs_dbs(contigs_dbs)
        self._check_contigs_dbs_ko_annotation(contigs_dbs)

        project_name_contigs_db: Dict[str, str] = {}
        contigs_db_project_name: Dict[str, str] = {}
        for contigs_db in contigs_dbs:
            contigs_db_info = dbinfo.ContigsDBInfo(contigs_db)
            self_table = contigs_db_info.get_self_table()
            project_name = self_table['project_name']
            assert project_name not in project_name_contigs_db
            project_name_contigs_db[project_name] = contigs_db
            contigs_db_project_name[contigs_db] = project_name

        self.progress.end()

        # Load groups.
        if (
            (groups_txt is None and group_threshold is not None) or
            (groups_txt is not None and group_threshold is None)
        ):
            raise ConfigError(
                "To group contigs databases, arguments to both 'groups_txt' and 'group_threshold' "
                "must be provided."
            )

        source_group: Dict[str, str] = None
        group_sources: Dict[str, List[str]] = None
        if groups_txt is not None:
            if not 0 <= group_threshold <= 1:
                raise ConfigError(
                    f"'group_threshold' must be a number between 0 and 1, not {group_threshold}"
                )

            # Match groups-txt entries, which may be relative paths, to the input databases by
            # absolute path, and relate groups to project names.
            groups_txt_source_group, groups_txt_group_sources = utils.get_groups_txt_file_as_dict(
                groups_txt, run=self.run, progress=self.progress
            )
            source_abspath_group = {
                os.path.abspath(source): group
                for source, group in groups_txt_source_group.items()
            }
            source_group = {}
            group_sources = {}
            ungrouped_contigs_dbs: List[str] = []
            for contigs_db in contigs_dbs:
                try:
                    group = source_abspath_group[os.path.abspath(contigs_db)]
                except KeyError:
                    ungrouped_contigs_dbs.append(contigs_db)
                    continue

                project_name = contigs_db_project_name[contigs_db]
                source_group[project_name] = group
                try:
                    group_sources[group].append(project_name)
                except KeyError:
                    group_sources[group] = [project_name]

            if ungrouped_contigs_dbs:
                message = ', '.join([f"'{contigs_db}'" for contigs_db in ungrouped_contigs_dbs])
                raise ConfigError(
                    "The following 'contigs_dbs' were not found in the groups provided by "
                    f"'groups_txt': {message}"
                )

            # Order groups by their appearance in the groups-txt file, which determines the order in
            # which colors are assigned to groups and their combinations in the drawn maps.
            group_sources = {
                group: group_sources[group]
                for group in groups_txt_group_sources
                if group in group_sources
            }

            # Report contigs databases in 'groups_txt' that are not among the input databases.
            contigs_db_abspaths = [os.path.abspath(contigs_db) for contigs_db in contigs_dbs]
            missing_sources: List[str] = []
            for source in groups_txt_source_group:
                if os.path.abspath(source) not in contigs_db_abspaths:
                    missing_sources.append(source)
            if missing_sources:
                message = ', '.join([f"'{source}'" for source in missing_sources])
                self.run.warning(
                    "The following contigs databases were grouped in 'groups_txt' but are not "
                    f"found among input 'contigs_dbs', and so will not factor into maps: {message}"
                )

        # Colors given per category are read before the databases are opened, so that a malformed
        # file is reported before the run spends time loading data it may not draw. They are checked
        # against the categories, which are the databases or their groups, in '_map_elements'.
        category_colors, combo_colors = self._read_category_colors_txt(
            reaction_category_colors, '--reaction-category-colors'
        ) if reaction_category_colors is not None else (None, {})

        self.progress.new("Loading KO data from contigs databases")
        self.progress.update("...")

        # Find which contigs databases contain each KO, and which KOs are in each database.
        ko_membership: Dict[str, List[str]] = {}
        source_kos: Dict[str, List[str]] = {}
        for project_name, contigs_db in project_name_contigs_db.items():
            cdb = ContigsDatabase(contigs_db)
            ko_ids = cdb.db.get_single_column_from_table(
                'gene_functions',
                'accession',
                unique=True,
                where_clause='source = "KOfam"'
            )
            source_kos[project_name] = ko_ids
            for ko_id in ko_ids:
                try:
                    ko_membership[ko_id].append(project_name)
                except KeyError:
                    ko_membership[ko_id] = [project_name]

        self.progress.end()

        layer = {
            'name': 'reactions',
            'element_type': 'reaction',
            'use_reaction_attribute': False,
            'membership': ko_membership,
            'source_accessions': {source: set(kos) for source, kos in source_kos.items()},
            'color_hexcode': reaction_color,
            'colormap': reaction_colormap,
            'colormap_limits': reaction_colormap_limits,
            'colormap_scheme': colormap_scheme,
            'category_colors_flag': '--reaction-category-colors',
            'category_colors': category_colors,
            'category_combo_colors': combo_colors,
            'reverse_overlay': reaction_reverse_overlay
        }
        return self._map_element_membership(
            [layer],
            list(project_name_contigs_db),
            'contigs database',
            source_group=source_group,
            group_sources=group_sources,
            group_threshold=group_threshold,
            pathway_numbers=pathway_numbers,
            draw_individual_files=draw_individual_files,
            draw_grid=draw_grid,
            group_colormap=group_colormap,
            group_colormap_limits=group_colormap_limits,
            group_reverse_overlay=group_reverse_overlay,
            group_colormap_scheme=group_colormap_scheme,
            count_scale_max=count_scale_max,
            output_dir=output_dir,
            draw_maps_lacking_data=draw_maps_lacking_data
        )

    def map_pan_database_kos(
        self,
        pan_db: str,
        genomes_storage_db: str,
        output_dir: str,
        groups_txt: str = None,
        group_threshold: float = None,
        consensus_threshold: float = None,
        discard_ties: bool = None,
        pathway_numbers: Iterable[str] = None,
        draw_individual_files: Union[Iterable[str], bool] = False,
        draw_grid: Union[Iterable[str], bool] = False,
        reaction_colormap: Union[bool, str, mcolors.Colormap] = True,
        reaction_colormap_limits: Tuple[float, float] = None,
        colormap_scheme: Literal['by_count', 'by_count_continuous', 'by_membership'] = None,
        reaction_category_colors: str = None,
        reaction_reverse_overlay: bool = False,
        reaction_color: str = '#2ca02c',
        group_colormap: Union[str, mcolors.Colormap] = 'plasma_r',
        group_colormap_limits: Tuple[float, float] = None,
        group_reverse_overlay: bool = False,
        group_colormap_scheme: Literal['by_count', 'by_count_continuous'] = None,
        count_scale_max: Union[str, int] = 'observed',
        draw_maps_lacking_data: bool = False
    ) -> Dict[Literal['unified', 'individual', 'grid'], Dict]:
        """
        Draw pathway maps, coloring the reaction layer by consensus KOs of gene clusters across
        genomes or groups of genomes (representing, for example, taxa or geographical groups).

        A reaction element on a map is represented by one or more KOs, matched to the consensus KOs
        of pangenomic gene clusters (a consensus KO is imputed to genomes with genes in the
        cluster); a genome "contains" the reaction if it has any of those consensus KOs. The
        reaction elements (lines in global and overview maps, boxes or lines in standard maps) are
        colored by the genomes or groups containing them via '_map_element_membership'.

        Parameters
        ==========
        pan_db : str
            File path to a pangenomic database.

        genomes_storage_db : str
            Path to the genomes storage database associated with the pan database. This must contain
            KO annotations.

        output_dir : str
            Path to the output directory in which pathway map and colorbar PDF files are drawn. The
            directory is created if it does not exist.

        groups_txt : str, None
            A tab-delimited text file specifying which group each genome belongs to. The first
            column, which can have any header, contains the names of genomes in the pan database.
            The second column, which must be headed 'group', contains group names, which are
            recommended to be single words without fancy characters, such as 'HIGH_TEMPERATURE' or
            'LOW_FITNESS' rather than 'my group #1' or 'IS-THIS-OK?'. Each genome can only be
            associated with a single group. The 'group_threshold' argument must also be used for the
            groups to take effect, assigning colors based on group membership and drawing individual
            files ('draw_individual_files') and map grids ('draw_grid') for groups rather than
            individual databases.

        group_threshold : float, None
            The proportion of genomes in a group containing data of interest for the group to be
            represented in terms of presence/absence in a reaction element. Here is a concrete
            example. Say the 'groups_txt' argument, which must be used with this argument, groups
            genomes by their species, 'A', 'B', and 'C'. You wish to understand the distribution of
            metabolic capabilities across the 3 species from KO annotations of genes. Reaction
            colors are assigned based on the groups rather than individual genomes containing the
            reaction. Thresholds between 0 and 1 can be set to define group membership: a threshold
            of 0.0 would mean that ANY genome in the group can contain the reaction via KOs for the
            reaction to be considered present in the group; a threshold of 0.75 means that at least
            75% of the genomes in the group must contain the reaction for it to be present; a
            threshold of 1.0 means that ALL genomes in the group must contain the reaction for it to
            be present. In our example, set the threshold to 0.5. Reaction J on a map corresponds to
            KO X, and Reaction K on a map corresponds to KOs Y and Z. 90% of species A genomes, 50%
            of species B genomes, and 10% of species C genomes contain KO X, so Reaction J would be
            colored to indicate that it is represented in species A and B. 0% of species A genomes,
            15% of species B genomes, and 40% of species C genomes contain KO Y or KO Z, so Reaction
            K would not be colored.

        consensus_threshold : float, None
            If a reaction ntework is stored in the pan database, then by default consensus KOs are
            determined using the 'reaction_network_consensus_threshold' value stored as database
            metadata. If a reaction network is not stored, then by default the consensus threshold
            is set to 0, meaning that the KO annotation most frequent in a gene cluster is assigned
            to the cluster as a whole. Alternatively, a number between 0 and 1 can be provided. At
            least this proportion of genes in the cluster must have the most frequent KO annotation
            for it to be assigned to the cluster as a whole.

        discard_ties : bool, None
            If a reaction network is stored in the pan database, then by default consensus KOs are
            determined using the 'reaction_network_discard_ties' value stored as database metadata.
            If a reaction network is not stored, then by default 'discard_ties' assumes a value of
            False. A value of True means that if multiple KO annotations are most frequent among
            genes in a cluster, then a consensus KO is not assigned to the cluster as a whole,
            whereas a value of False would cause one of the most frequent KOs to be arbitrarily
            chosen.

        pathway_numbers : Iterable[str], None
            Regex patterns to match the ID numbers of the drawn pathway maps. The default of None
            draws all available pathway maps in the KEGG data directory.

        reaction_color : str, '#2ca02c'
            The single color, by default green, for reaction elements when dynamic coloring is
            disabled (when 'reaction_colormap' is False). It colors both the unified map (by
            presence/absence in any genome) and the individual-genome maps. Alternatively, the
            string 'original' uses the reference map's original color scheme.

        draw_maps_lacking_data : bool, False
            If False, by default, only draw maps containing any of the select KOs. If True, draw
            maps regardless, meaning that nothing may be colored.

        Notes
        =====
        The dynamic-coloring and drawing options ('reaction_colormap', 'reaction_colormap_limits',
        'colormap_scheme', 'reaction_category_colors', 'reaction_reverse_overlay',
        'draw_individual_files', 'draw_grid', and the group colormap options) mirror the categorical
        engine; see the CLI help and '_map_element_membership'. 'reaction_category_colors' is the
        path to a kegg-category-colors-txt file giving a color per category — per source, or per
        group when the sources are grouped — which colors the layer by membership in place of a
        colormap, and colors each category's own map. 'group_colormap' may be
        'GROUP_COLORMAP_FROM_CATEGORY' to color each group's own maps by a ramp running to that
        group's own color instead of by a named colormap, in which case 'group_colormap_limits' is
        how far from white that ramp runs ('_group_map_colors'). 'group_colormap_scheme' draws the
        count scale of every group's maps in discrete bands ('by_count') or as a gradient
        ('by_count_continuous'), and by default keeps the bands while the colors can be told apart.

        Returns
        =======
        Dict[Literal['unified', 'individual', 'grid'], Dict]
            The record returned by '_map_element_membership': 'unified' maps show all genomes or
            groups, 'individual' maps show single genomes or groups, and 'grid' images show both.
            See '_map_element_membership' for the nested structure.
        """
        # This method loads consensus-KO membership from a pangenome and hands off the drawing
        # of unified, individual, and grid maps to '_map_element_membership'.

        self.progress.new("Loading metadata from pan database")
        self.progress.update("...")

        self._check_pan_db(pan_db)
        self._check_genomes_storage_db(genomes_storage_db)
        self._check_genomes_storage_ko_annotation(genomes_storage_db)

        # Load pan database metadata.
        pan_db_info = dbinfo.PanDBInfo(pan_db)
        self_table = pan_db_info.get_self_table()
        all_genome_names: List[str] = self_table['external_genome_names'].split(',')

        # Parameterize how consensus KOs are found.
        use_network_consensus_threshold = False
        if consensus_threshold is None:
            consensus_threshold = self_table['reaction_network_consensus_threshold']
            if consensus_threshold is not None:
                consensus_threshold = float(consensus_threshold)
                assert 0 <= consensus_threshold <= 1
                use_network_consensus_threshold = True

        use_network_discard_ties = False
        if discard_ties is None:
            discard_ties = self_table['reaction_network_discard_ties']
            if discard_ties is None:
                discard_ties = False
            else:
                discard_ties = bool(int(discard_ties))
                use_network_discard_ties = True

        self.progress.end()

        if use_network_consensus_threshold:
            self.run.info_single(
                "No consensus threshold was explicitly specified for consensus KO assignment to "
                f"gene clusters, but there was a value of '{consensus_threshold}' stored in the "
                "pan database from reaction network construction, so this was used. (The default "
                "if this were not the case is 0, or no threshold.)"
            )

        if use_network_discard_ties:
            self.run.info_single(
                "It was not explicitly specified whether to discard ties in consensus KO "
                f"assignment to gene clusters, but there was a value of '{discard_ties}' stored in "
                "the pan database from reaction network construction, so this was used. (The "
                "default if this were not the case is False, or do not discard ties.)"
            )

        # Load groups.
        if (
            (groups_txt is None and group_threshold is not None) or
            (groups_txt is not None and group_threshold is None)
        ):
            raise ConfigError(
                "To group genomes, arguments to both 'groups_txt' and 'group_threshold' must be "
                "provided."
            )

        source_group: Dict[str, str] = None
        group_sources: Dict[str, List[str]] = None
        if groups_txt is not None:
            if not 0 <= group_threshold <= 1:
                raise ConfigError(
                    f"'group_threshold' must be a number between 0 and 1, not {group_threshold}"
                )

            groups_txt_source_group, groups_txt_group_sources = utils.get_groups_txt_file_as_dict(
                groups_txt, run=self.run, progress=self.progress
            )

            # Check that groups include the pan genomes. Relate groups and genome names.
            source_group = {}
            group_sources = {}
            ungrouped_genomes: List[str] = []
            for genome_name in all_genome_names:
                try:
                    group = groups_txt_source_group[genome_name]
                except KeyError:
                    ungrouped_genomes.append(genome_name)
                    continue

                source_group[genome_name] = group
                try:
                    group_sources[group].append(genome_name)
                except KeyError:
                    group_sources[group] = [genome_name]

            if ungrouped_genomes:
                message = ', '.join([f"'{genome_name}'" for genome_name in ungrouped_genomes])
                raise ConfigError(
                    f"The following 'pan_db' genomes were not found in the groups provided by "
                    f"'groups_txt': {message}"
                )

            # Order groups by their appearance in the groups-txt file, which determines the order in
            # which colors are assigned to groups and their combinations in the drawn maps.
            group_sources = {
                group: group_sources[group]
                for group in groups_txt_group_sources
                if group in group_sources
            }

            # Report genomes in 'groups_txt' that are not in the pan database.
            missing_sources: List[str] = []
            for source in groups_txt_source_group:
                if source not in all_genome_names:
                    missing_sources.append(source)
            if missing_sources:
                message = ', '.join([f"'{source}'" for source in missing_sources])
                self.run.warning(
                    f"The following genomes were grouped in 'groups_txt' but are not found among "
                    f"'pan_db' genomes, and so will not factor into maps: {message}"
                )

        # Colors given per category are read before the gene clusters are loaded, so that a
        # malformed file is reported before the run spends time on data it may not draw. They are
        # checked against the categories, which are the genomes or their groups, in '_map_elements'.
        category_colors, combo_colors = self._read_category_colors_txt(
            reaction_category_colors, '--reaction-category-colors'
        ) if reaction_category_colors is not None else (None, {})

        self.progress.new("Loading consensus KO data from pan database")
        self.progress.update("...")

        # Load gene cluster data.
        progress = self.progress
        self.progress = terminal.Progress(verbose=False)
        run = self.run
        self.run = terminal.Run(verbose=False)
        args = Namespace()
        args.pan_db = pan_db
        args.genomes_storage = genomes_storage_db
        args.consensus_threshold = consensus_threshold
        args.discard_ties = discard_ties
        pan_super = PanSuperclass(args, r=self.run, p=self.progress)
        pan_super.init_gene_clusters()
        pan_super.init_gene_clusters_functions()
        pan_super.init_gene_clusters_functions_summary_dict()
        gene_clusters: Dict[str, Dict[str, List[int]]] = pan_super.gene_clusters
        gene_clusters_functions_summary_dict: Dict = pan_super.gene_clusters_functions_summary_dict
        self.progress = progress
        self.run = run

        # Find clusters with consensus KO annotations.
        consensus_cluster_kos: Dict[str, str] = {}
        for cluster_id, gene_cluster_functions_data in gene_clusters_functions_summary_dict.items():
            gene_cluster_ko_data = gene_cluster_functions_data['KOfam']
            if gene_cluster_ko_data == {'function': None, 'accession': None}:
                continue
            consensus_cluster_kos[cluster_id] = gene_cluster_ko_data['accession']

        # More than one gene cluster can be represented by the same consensus KO. Find which
        # genomes contribute genes to clusters represented by each consensus KO, and which consensus
        # KOs are in each genome.
        consensus_ko_genomes: Dict[str, List[str]] = {}
        genome_consensus_kos: Dict[str, List[str]] = {}
        for cluster_id, ko_id in consensus_cluster_kos.items():
            for genome_name, gcids in gene_clusters[cluster_id].items():
                if not gcids:
                    continue
                try:
                    consensus_ko_genomes[ko_id].append(genome_name)
                except KeyError:
                    consensus_ko_genomes[ko_id] = [genome_name]
                try:
                    genome_consensus_kos[genome_name].append(ko_id)
                except KeyError:
                    genome_consensus_kos[genome_name] = [ko_id]
        for ko_id, ko_genome_names in consensus_ko_genomes.items():
            consensus_ko_genomes[ko_id] = list(set(ko_genome_names))

        self.progress.end()

        # Give every genome a KO list, so that a genome lacking consensus KOs still gets an (empty)
        # individual map rather than raising a KeyError.
        source_kos: Dict[str, List[str]] = {
            genome_name: genome_consensus_kos.get(genome_name, [])
            for genome_name in all_genome_names
        }

        layer = {
            'name': 'reactions',
            'element_type': 'reaction',
            'use_reaction_attribute': False,
            'membership': consensus_ko_genomes,
            'source_accessions': {source: set(kos) for source, kos in source_kos.items()},
            'color_hexcode': reaction_color,
            'colormap': reaction_colormap,
            'colormap_limits': reaction_colormap_limits,
            'colormap_scheme': colormap_scheme,
            'category_colors_flag': '--reaction-category-colors',
            'category_colors': category_colors,
            'category_combo_colors': combo_colors,
            'reverse_overlay': reaction_reverse_overlay
        }
        return self._map_element_membership(
            [layer],
            all_genome_names,
            'pangenome',
            source_group=source_group,
            group_sources=group_sources,
            group_threshold=group_threshold,
            pathway_numbers=pathway_numbers,
            draw_individual_files=draw_individual_files,
            draw_grid=draw_grid,
            group_colormap=group_colormap,
            group_colormap_limits=group_colormap_limits,
            group_reverse_overlay=group_reverse_overlay,
            group_colormap_scheme=group_colormap_scheme,
            count_scale_max=count_scale_max,
            output_dir=output_dir,
            draw_maps_lacking_data=draw_maps_lacking_data
        )

    def _map_kos_fixed_colors(
        self,
        ko_ids: Iterable[str],
        output_dir: str,
        pathway_numbers: List[str] = None,
        color_hexcode: str = '#2ca02c',
        draw_maps_lacking_data: bool = False
    ) -> Dict[str, bool]:
        """
        Draw pathway maps, highlighting reactions containing select KOs in either a single color
        provided by a hex code or the colors originally used in the reference map.

        Parameters
        ==========
        ko_ids : Iterable[str]
            KO IDs to be highlighted in the maps.

        output_dir : str
            Path to the output directory in which pathway map PDF files are drawn. The directory is
            created if it does not exist.

        pathway_numbers : Iterable[str], None
            Regex patterns to match the ID numbers of the drawn pathway maps. The default of None
            draws all available pathway maps in the KEGG data directory.

        color_hexcode : str, '#2ca02c'
            This is the color, by default green, for reactions containing provided KOs.
            Alternatively to a color hex code, the string, 'original', can be provided to use the
            original color scheme of the reference map. In global maps, KOs are represented in
            reaction lines, and in overview maps, KOs are represented in reaction arrows. The
            foreground color of the lines and arrows is set. In standard maps, KOs are represented
            in boxes, the background color of which is set.

        draw_maps_lacking_data : bool, False
            If False, by default, only draw maps containing any of the select KOs. If True, draw
            maps regardless, meaning that nothing may be colored.

        Returns
        =======
        Dict[str, bool]
            Keys are pathway numbers. Values are True if the map was drawn, False if the map was not
            drawn because it did not contain any of the select KOs and 'draw_maps_lacking_data' was
            False.
        """
        # Find the numeric IDs of the maps to draw.
        pathway_numbers = self._find_maps(output_dir, patterns=pathway_numbers)

        filesnpaths.gen_output_directory(output_dir, progress=self.progress, run=self.run)

        # A single-color highlight is drawn through the generalized element engine as a lone
        # reaction (ortholog) layer. The 'original' scheme preserves each reaction's reference-map
        # colors and their render order, which the element engine does not model, so it keeps its
        # own drawer.
        ko_ids = set(ko_ids)
        presence_layers = [{
            'element_type': 'reaction',
            'use_reaction_attribute': False,
            'accessions': ko_ids,
            'color_hexcode': color_hexcode
        }]

        # Draw maps.
        self.progress.new("Drawing map")
        drawn: Dict[str, bool] = {}
        for pathway_number in pathway_numbers:
            self.progress.update(pathway_number)
            if color_hexcode == 'original':
                drawn[pathway_number] = self._draw_map_kos_original_color(
                    pathway_number,
                    ko_ids,
                    output_dir,
                    draw_map_lacking_data=draw_maps_lacking_data
                )
            else:
                drawn[pathway_number] = self._draw_map_element_presence(
                    pathway_number,
                    presence_layers,
                    output_dir,
                    draw_map_lacking_data=draw_maps_lacking_data
                )
        self.progress.end()

        return drawn

    def _relate_samples_to_groups(
        self,
        groups_txt: str,
        all_sample_names: List[str]
    ) -> Tuple[Dict[str, str], Dict[str, List[str]]]:
        """
        Load a groups file relating the input file's samples to groups.

        Every sample must be assigned to a group. Samples grouped in the file but absent from the
        input file are reported and ignored.

        Parameters
        ==========
        groups_txt : str
            Path to a tab-delimited groups file whose first column holds sample names (those in the
            input file's 'sample' column) and whose 'group' column holds group names.

        all_sample_names : List[str]
            The names of all samples in the input file.

        Returns
        =======
        Tuple[Dict[str, str], Dict[str, List[str]]]
            sample_group : maps each sample name to its group name.
            group_samples : maps each group name to its list of sample names, ordered by the group's
                appearance in the groups file.
        """
        all_sample_names_set = set(all_sample_names)
        source_group, group_sources = utils.get_groups_txt_file_as_dict(
            groups_txt, run=self.run, progress=self.progress
        )

        # Check that groups include all samples. Relate groups and sample names.
        group_samples: Dict[str, List[str]] = {}
        sample_group: Dict[str, str] = {}
        ungrouped_samples: List[str] = []
        for sample_name in all_sample_names:
            try:
                group = source_group[sample_name]
            except KeyError:
                ungrouped_samples.append(sample_name)
                continue

            try:
                group_samples[group].append(sample_name)
            except KeyError:
                group_samples[group] = [sample_name]
            sample_group[sample_name] = group

        if ungrouped_samples:
            message = ', '.join([f"'{sample_name}'" for sample_name in ungrouped_samples])
            raise ConfigError(
                f"The following samples in the draw-kegg-pathways text file were not found in the "
                f"groups provided by 'groups_txt': {message}"
            )

        # Order groups by their appearance in the groups-txt file, which determines the order in
        # which colors are assigned to groups and their combinations in the drawn maps.
        group_samples = {
            group: group_samples[group]
            for group in group_sources
            if group in group_samples
        }

        # Report samples in 'groups_txt' that are not among the draw-kegg-pathways text file
        # samples.
        missing_sources: List[str] = []
        for source in source_group:
            if source not in all_sample_names_set:
                missing_sources.append(source)
        if missing_sources:
            message = ', '.join([f"'{source}'" for source in missing_sources])
            self.run.warning(
                f"The following samples were grouped in 'groups_txt' but are not found among the "
                f"samples in the draw-kegg-pathways text file, and so will not factor into maps: "
                f"{message}"
            )

        return sample_group, group_samples

    def _read_category_colors_txt(
        self,
        path: str,
        flag: str
    ) -> Tuple[Dict[str, str], Dict[Tuple[str, ...], str]]:
        """
        Load a file of per-category colors, and of colors overriding category combinations.

        The file names a color for each category — each sample, contigs database, genome, or group —
        which coloring by membership uses for the elements in that category alone, an individual
        map uses for its own category, and a group's color ramp runs to. A row whose first field
        lists several names separated by 'CATEGORY_COMBO_SEPARATOR' instead gives the color of that
        combination of categories, replacing the blend of their colors that coloring by membership
        would otherwise derive.

        The names are not checked against the run's actual categories here, since a file is read
        before they are all known; '_resolve_category_colors' does that.

        Parameters
        ==========
        path : str
            Path to a tab-delimited file whose first column holds category names, or combinations of
            them, and whose 'CATEGORY_COLORS_COLUMN' column holds color hex codes.

        flag : str
            The command-line flag this file came from, used in error messages.

        Returns
        =======
        Tuple[Dict[str, str], Dict[Tuple[str, ...], str]]
            category_colors : maps each category name to its color hex code.
            combo_colors : maps each combination of category names, sorted and as a tuple, to the
                color hex code overriding the blend of that combination.
        """
        filesnpaths.is_file_tab_delimited(path)
        columns = utils.get_columns_of_TAB_delim_file(path, include_first_column=True)
        if CATEGORY_COLORS_COLUMN not in columns:
            self.progress.end()
            raise ConfigError(
                f"The colors file given to '{flag}', at '{path}', should have a column called "
                f"'{CATEGORY_COLORS_COLUMN}' holding a color hex code for each category. Its "
                f"columns are: {', '.join(repr(column) for column in columns)}."
            )
        if len(columns) < 2:
            self.progress.end()
            raise ConfigError(
                f"The colors file given to '{flag}', at '{path}', should have at least two "
                f"columns: the first one holding the names of the categories — the samples, "
                f"contigs databases, genomes, or groups being colored — and a "
                f"'{CATEGORY_COLORS_COLUMN}' column holding a color hex code for each of them."
            )
        if columns[0] == CATEGORY_COLORS_COLUMN:
            self.progress.end()
            raise ConfigError(
                f"The first column of the colors file given to '{flag}', at '{path}', is the "
                f"'{CATEGORY_COLORS_COLUMN}' column, but anvi'o expects the first column to hold "
                f"category names, so the two columns need to be the other way around."
            )

        table = utils.get_TAB_delimited_file_as_dictionary(path)

        category_colors: Dict[str, str] = {}
        combo_colors: Dict[Tuple[str, ...], str] = {}
        # What each normalized key was written as in the file, for reporting two rows that mean one
        # thing.
        seen_rows: Dict[Tuple[str, ...], str] = {}
        for name, row in table.items():
            color = str(row[CATEGORY_COLORS_COLUMN]).strip()
            # Only hex codes are accepted, even though Matplotlib would also read a color name here:
            # the colors a map reserves for its own unhighlighted elements are compared as hex codes
            # ('_check_reserved_colors'), and a file of hex codes is what the drawn colorbars and
            # the '--reaction-color'/'--compound-color' options already speak.
            if not re.fullmatch(r'#[0-9A-Fa-f]{6}', color):
                self.progress.end()
                raise ConfigError(
                    f"Each color in the file given to '{flag}', at '{path}', must be a six-digit "
                    f"hex code such as '#FFA500', but the color of '{name}' is '{color}'."
                )
            members = tuple(
                sorted(
                    {part.strip() for part in name.split(CATEGORY_COMBO_SEPARATOR) if part.strip()}
                )
            )
            if not members:
                self.progress.end()
                raise ConfigError(
                    f"A row of the colors file given to '{flag}', at '{path}', names no category "
                    f"at all in its first column. Every row should name a category, or a "
                    f"combination of them separated by '{CATEGORY_COMBO_SEPARATOR}'."
                )
            # Two rows can name the same thing without being the same text — stray spaces around a
            # name, or the members of a combination written in another order — and the duplicate
            # check that reading the file has already done compares the text. Left alone, the later
            # row would quietly replace the earlier one and the file would report a count that does
            # not match its rows, so what the text meant is checked here too.
            if members in seen_rows:
                separator = f'{CATEGORY_COMBO_SEPARATOR} '
                self.progress.end()
                raise ConfigError(
                    f"The colors file given to '{flag}', at '{path}', gives two colors to the same "
                    f"{'combination' if len(members) > 1 else 'category'}: "
                    f"'{seen_rows[members]}' and '{name}' both name "
                    f"'{separator.join(members)}'. Please give it one row."
                )
            seen_rows[members] = name
            if len(members) == 1:
                category_colors[members[0]] = kgml.canonical_color(color)
            else:
                combo_colors[members] = kgml.canonical_color(color)

        if not category_colors and not combo_colors:
            self.progress.end()
            raise ConfigError(
                f"The colors file given to '{flag}', at '{path}', has its header but no rows, so "
                f"it gives no color to anything. Each row should name a category in its first "
                f"column and give its color hex code in the '{CATEGORY_COLORS_COLUMN}' column."
            )
        if not category_colors:
            self.progress.end()
            raise ConfigError(
                f"The colors file given to '{flag}', at '{path}', gives colors only for "
                f"combinations of categories, and none for a category on its own. A combination's "
                f"color adjusts the blend of its members' colors, so the members need colors first."
            )

        self.run.info(f"Categories colored by '{flag}'", len(category_colors))
        if combo_colors:
            self.run.info(f"Category combinations recolored by '{flag}'", len(combo_colors))
        return category_colors, combo_colors

    def _resolve_category_colors(
        self,
        layers: List[dict],
        categories: List[str], category_noun: str
    ) -> None:
        """
        Reduce each layer's category colors to those the run will actually use.

        A category with no color of its own cannot be colored at all, so it is an error, reported
        for every such category at once. A color given for a name that is not a category, or for a
        combination naming one, describes something this run does not draw, so it is reported and
        ignored, as a groups file's extra items are ('_relate_samples_to_groups'): one file can then
        cover a set of samples that different runs draw different subsets of. Ignoring such a color
        means dropping it from the layer rather than merely reporting it, which is what makes the
        colors that remain the colors of the run — everything downstream reads them as such.

        Parameters
        ==========
        layers : List[dict]
            The layer models, whose 'category_colors'/'category_combo_colors' are checked against
            the run and reduced to it in place. A layer without colors of its own is skipped.

        categories : List[str]
            The names of every category of the run.

        category_noun : str
            What one category is called, for error messages, e.g. 'sample' or 'group'.
        """
        category_set = set(categories)
        # A combination is written as its names separated by 'CATEGORY_COMBO_SEPARATOR', so a name
        # containing that separator could not be told from a combination of other names.
        containing_separator = sorted(
            name for name in category_set if CATEGORY_COMBO_SEPARATOR in name
        )
        if containing_separator:
            message = ', '.join(f"'{name}'" for name in containing_separator)
            self.progress.end()
            raise ConfigError(
                f"A colors file names a combination of {category_noun}s by separating their names "
                f"with '{CATEGORY_COMBO_SEPARATOR}', so no {category_noun} name may itself contain "
                f"that character. "
                f"{'These names do' if len(containing_separator) > 1 else 'This name does'}: "
                f"{message}. Please rename {'them' if len(containing_separator) > 1 else 'it'} in "
                f"the input, or drop the colors file and color these {category_noun}s from a "
                f"colormap instead."
            )

        for layer in layers:
            category_colors = layer.get('category_colors')
            if category_colors is None:
                continue
            flag = layer['category_colors_flag']

            missing = [category for category in categories if category not in category_colors]
            if missing:
                message = ', '.join(f"'{category}'" for category in missing)
                self.progress.end()
                raise ConfigError(
                    f"The colors file given to '{flag}' gives no color for "
                    f"{'these' if len(missing) > 1 else 'this'} {category_noun}"
                    f"{'s' if len(missing) > 1 else ''} of the run: {message}. Coloring by "
                    f"membership needs a color for every {category_noun}, since an element is "
                    f"colored by exactly which of them contain it, so anvi'o will not fill in the "
                    f"gaps from a colormap: that would put colors chosen here and colors chosen "
                    f"for you on one scale. Please add {'them' if len(missing) > 1 else 'it'} to "
                    f"the file."
                )

            extra = sorted(name for name in category_colors if name not in category_set)
            extra_combos = [
                combo for combo in layer['category_combo_colors']
                if not set(combo) <= category_set
            ]
            if extra:
                message = ', '.join(f"'{name}'" for name in extra)
                self.run.warning(
                    f"The colors file given to '{flag}' colors the following names, which are not "
                    f"{category_noun}s of this run and so will not factor into the maps: {message}"
                )
            if extra_combos:
                separator = f'{CATEGORY_COMBO_SEPARATOR} '
                message = ', '.join(
                    f"'{name}'" for name in sorted(separator.join(combo) for combo in extra_combos)
                )
                self.run.warning(
                    f"The colors file given to '{flag}' recolors the following combinations, which "
                    f"name at least one thing that is not a {category_noun} of this run, and so "
                    f"will not factor into the maps: {message}"
                )

            # Drop them, rather than only report them. Everything downstream reads what is left as
            # the layer's colors: the check that no color is one a map reserves, and the check that
            # the layers agree on each group's ramp. A color kept here for something the run does
            # not draw could fail either check, refusing a run it has no part in.
            for name in extra:
                del category_colors[name]
            for combo in extra_combos:
                del layer['category_combo_colors'][combo]

    @staticmethod
    def _trim_colormap(
        cmap: Union[mcolors.Colormap, None],
        colormap_limits: Tuple[float, float]
    ) -> Union[mcolors.Colormap, None]:
        """
        Trim a colormap to a fraction of its range.

        Parameters
        ==========
        cmap : Union[matplotlib.colors.Colormap, None]
            The colormap to trim, or None.

        colormap_limits : Tuple[float, float], None
            Lower and upper cutoffs on the fraction of the colormap to keep, e.g., (0.2, 0.9) trims
            the bottom 20% and top 10%. If None or (0.0, 1.0), or if 'cmap' is None, the colormap is
            returned unchanged.

        Returns
        =======
        Union[matplotlib.colors.Colormap, None]
            The trimmed colormap, or the input returned unchanged.
        """
        if cmap is None or colormap_limits is None or colormap_limits == (0.0, 1.0):
            return cmap
        lower_limit = colormap_limits[0]
        upper_limit = colormap_limits[1]
        assert 0.0 <= lower_limit <= upper_limit <= 1.0
        return mcolors.LinearSegmentedColormap.from_list(
            f'trunc({cmap.name},{lower_limit:.2f},{upper_limit:.2f})',
            cmap(range(int(lower_limit * cmap.N), math.ceil(upper_limit * cmap.N)))
        )

    def _get_colormap(self, colormap: str) -> mcolors.Colormap:
        """
        Look up a Matplotlib colormap by name.

        Parameters
        ==========
        colormap : str
            The name of a Matplotlib colormap.

        Returns
        =======
        matplotlib.colors.Colormap
            The named colormap.
        """
        try:
            return colormaps[colormap]
        except KeyError:
            self.progress.end()
            raise ConfigError(
                f"'{colormap}' is not the name of a Matplotlib colormap. The names are listed at "
                f"https://matplotlib.org/stable/users/explain/colors/colormaps.html and include "
                f"'plasma', 'viridis' and 'tab10'; adding '_r' to a name, as in 'plasma_r', "
                f"reverses its colors."
            )

    def _resolve_sequential_colormap(
        self,
        colormap: Union[str, mcolors.Colormap],
        colormap_limits: Tuple[float, float],
        subject: str = 'elements'
    ) -> mcolors.Colormap:
        """
        Return a Matplotlib Colormap for quantitative coloring, trimmed to 'colormap_limits'.

        Parameters
        ==========
        colormap : Union[str, matplotlib.colors.Colormap]
            A sequential colormap or its name.

        colormap_limits : Tuple[float, float], None
            Lower and upper cutoffs on the fraction of the colormap to use, or None for the full
            colormap.

        subject : str, 'elements'
            What the colormap colors, for the warning about a qualitative colormap, e.g. 'reactions'
            or 'compounds'.

        Returns
        =======
        matplotlib.colors.Colormap
            The (optionally trimmed) colormap.
        """
        if isinstance(colormap, str):
            cmap = self._get_colormap(colormap)
        elif isinstance(colormap, mcolors.Colormap):
            cmap = colormap
        else:
            self.progress.end()
            raise ConfigError(
                f"A colormap must be given as the name of a Matplotlib colormap or as a Colormap "
                f"object, but this one is neither: {colormap}."
            )

        if cmap.name in qualitative_colormaps + repeating_colormaps:
            self.run.warning(
                f"The colormap, '{cmap.name}', provided to color {subject} by value is qualitative "
                f"rather than sequential, which makes a continuous color scale difficult to "
                f"interpret. We recommend a sequential colormap like 'plasma' instead."
            )

        return self._trim_colormap(cmap, colormap_limits)

    @staticmethod
    def _check_requested_subset(
        request: Union[Iterable[str], bool],
        request_phrase: str,
        valid_names: Set[str],
        subject: str
    ) -> None:
        """
        Raise a ConfigError if a request for a subset of sources names any unrecognized source.

        Individual map files ('draw_individual_files') and map grid panels ('draw_grid') can be
        requested for all sources (True), no sources (False), or a subset of sources (a list of
        names). This checks the subset case, and is a no-op otherwise.

        Parameters
        ==========
        request : Union[Iterable[str], bool]
            A 'draw_individual_files' or 'draw_grid' value.

        request_phrase : str
            How to refer to the request in an error message, e.g., 'Individual maps' or 'Individual
            maps in grids'.

        valid_names : Set[str]
            The recognized source names.

        subject : str
            A plural noun for the sources in an error message, e.g., 'samples', 'contigs databases',
            or 'sample groups'.
        """
        if isinstance(request, bool):
            return
        unrecognized = [name for name in request if name not in valid_names]
        if unrecognized:
            message = ', '.join(f"'{name}'" for name in unrecognized)
            raise ConfigError(
                f"{request_phrase} were requested for a subset of {subject}, but the following "
                f"names were not recognized as any of the {subject}: {message}"
            )

    def _make_quantitative_norm(
        self,
        values: List[float],
        limits: Union[Tuple[Union[float, None], Union[float, None]], None] = None,
        flag: Union[str, None] = None,
        center: Union[float, None] = None,
        center_flag: Union[str, None] = None,
        subject: str = 'these maps'
    ) -> Tuple[Union[mcolors.Normalize, None], Union[float, None], Union[float, None], bool, bool]:
        """
        Make a normalization over reaction values for quantitative coloring.

        'limits' truncate the range the scale spans, but only where the values actually cross them:
        a limit no value passes leaves the scale exactly where the values put it, so that a limit
        set to guard against a long tail does not stretch a scale that turns out not to have one.
        Where a limit does truncate, every value beyond it takes the color of that end of the scale,
        which 'clip=True' on the normalization arranges, and the caller labels that end of the
        colorbar '<=' or '>=' so that its color reads as "this value or past it" rather than as an
        exact value.

        'center' then widens whichever side of the range falls short, until the range runs the same
        distance either side of it. That is what puts the centered value at the middle of the
        colormap — the neutral color of a diverging one — however lopsided the values around it are,
        and it keeps the scale linear, so that the same distance in color goes on meaning the same
        distance in value on either side. Widening comes after truncating, so a limit that was
        holding a tail in check can be pushed past by the other side of the range; that undoes what
        the limit was asked to do, and is refused rather than carried out quietly.

        Parameters
        ==========
        values : List[float]
            All per-reaction values that will be colored on the maps sharing this normalization.

        limits : Union[Tuple[Union[float, None], Union[float, None]], None], None
            The (minimum, maximum) the scale may span, as '_resolve_value_limits' returns them.
            Either end can be None to leave it wherever the values put it, and None means the scale
            is not limited at all.

        flag : Union[str, None], None
            The command-line flag 'limits' came from, used in an error message.

        center : Union[float, None], None
            The value to put at the middle of the scale, as '_resolve_value_center' returns it, or
            None to leave the scale wherever the values and limits leave it.

        center_flag : Union[str, None], None
            The command-line flag 'center' came from, used in messages.

        subject : str, 'these maps'
            What the scale colors, e.g. "the 'unified' map", for messages about the centering.

        Returns
        =======
        Tuple[Union[matplotlib.colors.Normalize, None],
              Union[float, None], Union[float, None], bool, bool]
            The normalization, its (vmin, vmax), and whether values fall below the bottom and above
            the top of the range, which is where the colorbar is marked '<=' or '>='. The first
            three are None if 'values' is empty. The normalization is None but vmin and vmax are set
            (and equal) if the range is a single value, degenerate because all values are equal,
            because a limit landed on the far end of them, or because a centered scale has nothing
            but its own center to span, in which case callers map every reaction to one end of the
            colormap, or to its middle where the scale is centered.
        """
        if not values:
            return None, None, None, False, False
        vmin = data_min = min(values)
        vmax = data_max = max(values)
        limit_min, limit_max = (None, None) if limits is None else limits
        if limit_min is not None and data_min < limit_min:
            vmin = limit_min
        if limit_max is not None and data_max > limit_max:
            vmax = limit_max
        if vmin > vmax:
            # A limit lies past the far end of the values, so nothing is left between the two ends
            # of the scale for it to span. Only one of the two limits can be at fault: a pair whose
            # minimum is below its maximum ('_resolve_value_limits') can miss the values only by
            # sitting wholly to one side of them, so the end that does is the one named.
            if limit_min is not None and limit_min > data_max:
                end_noun, limit, side, remedy = 'minimum', limit_min, 'below', f'below {data_max:g}'
            else:
                end_noun, limit, side, remedy = 'maximum', limit_max, 'above', f'above {data_min:g}'
            raise ConfigError(
                f"'{flag}' was given a {end_noun} of {limit:g}, but every value coloring these "
                f"maps falls {side} it -- they run from {data_min:g} to {data_max:g}. That leaves "
                f"the color scale nothing to span: every element would take the same end color, "
                f"and the map could not tell one from another. Please set the {end_noun} {remedy}, "
                f"or drop the limits."
            )
        if center is not None:
            # The farther of the two ends from the center sets how far the scale reaches on both
            # sides. Only one of the two distances can be negative, which happens when the whole
            # range lies to one side of the center, so the larger of them is never negative.
            low_distance = center - vmin
            high_distance = vmax - center
            half_range = max(low_distance, high_distance)
            centered_min, centered_max = center - half_range, center + half_range
            for (
                limit, distance, centered_end, data_end, end_noun, other_end_noun, comparison,
                loosen_verb
            ) in (
                (
                    limit_min, low_distance, centered_min, data_min, 'minimum', 'maximum', 'below',
                    'lower'
                ),
                (
                    limit_max, high_distance, centered_max, data_max, 'maximum', 'minimum', 'above',
                    'raise'
                )
            ):
                # A limit is undone here only when two things are both true. It has to have been
                # truncating something to begin with: a limit no value crossed is doing nothing
                # already, by design, so centering takes nothing away from it. And this end has to
                # be the nearer of the two to the center, since centering leaves the farther end
                # exactly where it is and stretches only the nearer one out to match. That is what
                # 'half_range > distance' asks, of the distances rather than of the widened end, so
                # that no end has to be recomputed to answer it.
                if limit is None:
                    continue
                truncating = data_end < limit if comparison == 'below' else data_end > limit
                if not (truncating and half_range > distance):
                    continue
                # The two ways out that keep a limit are worth working out for the user rather than
                # describing: where the centered scale ends, which is as far as this limit can be
                # moved, and the mirror of the limit about the center, which is the other limit that
                # would truncate a centered scale evenly.
                mirrored = 2 * center - limit
                other_end = vmax if comparison == 'below' else vmin
                raise ConfigError(
                    f"'{flag}' gives the color scale of {subject} a {end_noun} of {limit:g}, which "
                    f"truncates values reaching {data_end:g}, while '{center_flag}' centers that "
                    f"same scale on {center:g}. The two cannot both be had: with its other end at "
                    f"{other_end:g}, a scale centered on {center:g} has its {end_noun} at "
                    f"{centered_end:g}, which lies {comparison} the limit and undoes the "
                    f"truncation the limit was for. Three things would settle it: {loosen_verb} "
                    f"the {end_noun} to {centered_end:g}, which is where a centered scale ends "
                    f"here; give the scale a {other_end_noun} of {mirrored:g} as well, the same "
                    f"distance from {center:g} as the {end_noun}, which truncates a centered scale "
                    f"at both ends; or drop the {end_noun} and let that end fall where the values "
                    f"put it."
                )
            vmin, vmax = centered_min, centered_max
            # A scale that came out with nothing to span at all is a single band on the center
            # itself, where there is no unused half to speak of.
            if half_range > 0 and (data_min >= center or data_max <= center):
                values_phrase = (
                    f"which are every one of them {data_min:g}" if data_min == data_max else
                    f"which run from {data_min:g} to {data_max:g}"
                )
                self.run.warning(
                    f"'{center_flag}' puts {center:g} at the middle of the color scale of "
                    f"{subject}, but the values it colors, {values_phrase}, all lie on one side of "
                    f"it. Half of the colormap therefore goes unused and the values are squeezed "
                    f"into the other half. This may well be what you intend — it is what keeps one "
                    f"scale comparable across datasets that straddle the center to different "
                    f"degrees — but if the values here are not meant to be read against that "
                    f"center, drop the option and let the scale span the values themselves."
                )
        # Which ends the colorbar marks is asked of the range that ended up being drawn, so that a
        # mark means what it says: values do lie past this end. Centering can widen a truncated end
        # past every value, at which point there is nothing beyond it to mark.
        clamped_low = data_min < vmin
        clamped_high = data_max > vmax
        if vmin == vmax:
            return None, vmin, vmax, clamped_low, clamped_high
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax, clip=True)
        return norm, vmin, vmax, clamped_low, clamped_high

    @staticmethod
    def _get_entry_kegg_ids(entry: kgml.Entry, use_reaction_attribute: bool = False) -> List[str]:
        """
        Return the KEGG accessions an Entry represents.

        By default these come from the Entry's 'name' attribute: KO IDs for ortholog entries (e.g.,
        'ko:K00844 ko:K12407' -> ['K00844', 'K12407']) and compound IDs for compound entries (e.g.,
        'cpd:C00031' -> ['C00031']). With 'use_reaction_attribute', they come instead from the
        'reaction' attribute of ortholog entries: reaction IDs (e.g., 'rn:R00764 rn:R00756' ->
        ['R00764', 'R00756']). The source can name multiple accessions or be absent.

        Parameters
        ==========
        entry : kgml.Entry
            The Entry to read.

        use_reaction_attribute : bool, False
            If True, read reaction IDs from the 'reaction' attribute instead of KO/compound IDs from
            the 'name' attribute.

        Returns
        =======
        List[str]
            The accessions (the part of each token after the colon), or an empty list if the source
            attribute is absent.
        """
        source = entry.reaction if use_reaction_attribute else entry.name
        if not source:
            return []
        return [token.split(':')[1] for token in source.split()]

    @staticmethod
    def _reduce_entry_value(
        entry: kgml.Entry,
        values: Dict[str, float],
        aggregate,
        use_reaction_attribute: bool = False
    ) -> Union[float, None]:
        """
        Reduce the values of an Entry's KEGG accessions to a single element value.

        Parameters
        ==========
        entry : kgml.Entry
            An Entry (ortholog or compound) whose accessions are read by '_get_entry_kegg_ids'.

        values : Dict[str, float]
            Keys are KEGG accessions (KO, reaction, or compound IDs), values are per-accession
            values.

        aggregate : callable
            Reduces a list of values to a single value (see 'AGGREGATION_FUNCTIONS').

        use_reaction_attribute : bool, False
            Passed to '_get_entry_kegg_ids' to read reaction IDs rather than KO/compound IDs.

        Returns
        =======
        Union[float, None]
            The aggregated value, or None if none of the Entry's accessions have a value or the
            aggregation is undefined for the values they do have.
        """
        entry_values = [
            values[kegg_id]
            for kegg_id in Mapper._get_entry_kegg_ids(entry, use_reaction_attribute)
            if kegg_id in values
        ]
        if not entry_values:
            return None
        # An aggregation can be undefined for the values of this Entry even when it is defined
        # elsewhere on the map, as the standard deviation is for an element with a single accession,
        # in which case the Entry has no value and is left uncolored.
        value = aggregate(entry_values)
        return value if np.isfinite(value) else None

    def _draw_quantitative_colorbar(
        self,
        cmap: mcolors.Colormap,
        vmin: float,
        vmax: float,
        out_path: str,
        label: str,
        integer_ticks: bool = False,
        clamped_low: bool = False,
        clamped_high: bool = False,
        center: Union[float, None] = None
    ) -> None:
        """
        Draw a continuous colorbar for quantitative coloring, or a single-value colorbar when the
        value range is degenerate (vmin == vmax).

        Parameters
        ==========
        cmap : matplotlib.colors.Colormap
            The colormap sampled across the value range.

        vmin : float
            Lower bound of the value range.

        vmax : float
            Upper bound of the value range.

        out_path : str
            Path to the PDF output file.

        label : str
            Overall colorbar label.

        integer_ticks : bool, False
            If True, label the bar at whole numbers spanning the range rather than at Matplotlib's
            automatic ticks, for a range that counts things rather than measuring them.

        clamped_low : bool, False
            If True, a value limit truncated the bottom of the range, so values below 'vmin' are
            drawn in the color at that end and its label is marked accordingly.

        clamped_high : bool, False
            The same for the top of the range and values above 'vmax'.

        center : Union[float, None], None
            The value the range was centered on, which the bar is ticked at so that a reader can see
            where the middle of the scale is, and which is the one value a degenerate centered range
            can have left.
        """
        if vmin == vmax:
            # A single band, either because every value is the same or because a limit landed on the
            # far end of them; in the latter case the band stands for everything past the limit too.
            # A centered range that collapsed has collapsed onto its center, so the band takes the
            # middle color rather than the top one.
            prefix = CLAMPED_MIN_PREFIX if clamped_low else (
                CLAMPED_MAX_PREFIX if clamped_high else ''
            )
            self.colorbar_drawer.draw_discrete(
                [mcolors.rgb2hex(cmap(0.5 if center is not None else 1.0))], out_path,
                color_labels=[f'{prefix}{vmin:g}'], label=label
            )
        else:
            self.colorbar_drawer.draw_continuous(
                cmap, vmin, vmax, out_path, label=label, integer_ticks=integer_ticks,
                clamped_low=clamped_low, clamped_high=clamped_high, center=center
            )

    def _map_elements(
        self,
        layers: List[dict],
        output_dir: str,
        pathway_numbers: Iterable[str] = None,
        categories: Union[List[str], None] = None,
        category_noun: Union[str, None] = None,
        grid_source_type: Union[str, None] = None,
        colorbar_category_suffix: Union[str, None] = None,
        subset_subject: Union[str, None] = None,
        unified_plural: Union[str, None] = None,
        membership_count_label: Union[str, None] = None,
        membership_members_label: Union[str, None] = None,
        membership_singular: Union[str, None] = None,
        grouped_membership: Union[dict, None] = None,
        count_scale_max: Union[str, int] = 'observed',
        draw_individual_files: Union[Iterable[str], bool] = False,
        draw_grid: Union[Iterable[str], bool] = False,
        draw_maps_lacking_data: bool = False
    ) -> Dict[Literal['unified', 'individual', 'grid'], Dict]:
        """
        The unified element engine: color one or two map layers, each by its own mode, on one map.

        A layer's mode sets how it colors its elements: 'quantitative' (continuous value),
        'membership' (by the sources/groups containing an element, or their count), 'single' (one
        fixed presence/absence color), 'static' (one fixed color pooled across sources, for the
        db/pan single-color path), or 'original' (preserve the reference map's colors, drawn by a
        separate drawer). The mode can differ between the two map contexts, given as 'unified_mode'
        for the 'unified' map and 'category_mode' for the individual sample/source/group maps: a
        text layer with a value column and a 'sample' column, for instance, colors per-sample
        magnitude continuously while summarizing the samples or groups by presence on the 'unified'
        map. Note that 'original' is a 'unified_mode' only: an individual map can be drawn in the
        reference colors for one SOURCE, since 'source_accessions' is keyed by source, but not for a
        group, so an original layer pairs it with a 'category_mode' of 'membership' and its grouped
        individual maps show within-group source counts.

        Because every layer is reduced to a per-Entry '(color, priority)' by '_draw_map_elements', a
        quantitative layer and a presence/categorical layer can be colored on the same map. A
        'unified' map colors elements summarized across all categories; when categories
        (samples/sources/groups) are drawn individually, each gets its own map, sharing one colorbar
        per quantitative layer so colors stay comparable.

        Parameters
        ==========
        layers : List[dict]
            One or two layer models, ordered so layers drawn beneath others come first (reaction
            before compound). Common keys: 'name' (colorbar filename stem), 'accessions' (every
            accession the layer touches, used to find the entries it colors), 'element_type'
            ('reaction'/'compound'), and 'use_reaction_attribute' (bool). A layer declares its
            per-context modes with 'unified_mode' and 'category_mode', or a single 'mode' used in
            every context (except that a 'static'/'original' layer renders as within-group source
            counts on grouped individual maps). Mode-specific keys: quantitative -> 'cmap',
            'reverse_overlay', 'unified_values', 'category_values' (or None), 'aggregate',
            'colorbar_label', the optional 'value_limits'/'category_value_limits' bounding each
            context's scale and 'value_center'/'category_value_center' putting a value at the middle
            of it ('_make_quantitative_norm'), and the optional 'category_cmap' coloring the
            per-category scale from a colormap of its own rather than from 'cmap';
            membership/static/original -> 'membership', 'source_accessions', 'color_hexcode', and
            (membership) 'colormap'/'colormap_limits'/'colormap_scheme'/'scheme_options'/
            'reverse_overlay'; single -> 'accessions', 'color_hexcode'. 'scheme_options' names the
            option that chooses this layer's presence scheme, which differs by input
            ('PRESENCE_SCHEME_OPTIONS'). A layer that is quantitative in one context and membership
            in the other carries the keys of both, plus 'accessions' for finding the entries whose
            values set the ranges. A membership layer may carry a color per category instead of a
            colormap, as 'category_colors' with the 'category_combo_colors' overriding combinations
            of them and the 'category_colors_flag' they came from; those colors then color the
            membership scale and each category's own map, and are checked against the categories
            here.

        output_dir : str
            Path to the output directory in which pathway map and colorbar PDF files are drawn.

        categories : Union[List[str], None]
            The shared category names (samples, sources, or groups) that get individual maps, in
            color-assignment order, or None for a single map with no category dimension.

        Notes
        =====
        'category_noun'/'grid_source_type'/'colorbar_category_suffix'/'subset_subject'/
        'unified_plural'/'membership_*'/'grouped_membership' carry the terminology and grouping that
        the public methods supply; the remaining parameters mirror those methods. 'grid_source_type'
        labels the per-group grid colorbar (which counts sources) and defaults to 'category_noun'.
        Quantitative layers carry their own colorbar label in 'colorbar_label'.
        """
        has_categories = categories is not None
        grouped = grouped_membership is not None
        source_group = grouped_membership['source_group'] if grouped else None
        group_sources = grouped_membership['group_sources'] if grouped else None
        group_threshold = grouped_membership['group_threshold'] if grouped else None

        # A layer that does not declare per-context modes colors the same way in every context. The
        # exception is a 'static' or 'original' layer, which pools its sources for the 'unified' map
        # but still renders as within-group source counts on grouped individual maps.
        for layer in layers:
            if 'unified_mode' not in layer:
                layer['unified_mode'] = layer['mode']
            if 'category_mode' not in layer:
                layer['category_mode'] = (
                    'membership' if layer['mode'] in ('static', 'original') else layer['mode']
                )
            # A colormap for the category scale alone is optional: without one, the two contexts are
            # colored from the layer's single 'cmap'.
            if layer['category_mode'] == 'quantitative' and 'category_cmap' not in layer:
                layer['category_cmap'] = layer['cmap']

        original_run = any(layer['unified_mode'] == 'original' for layer in layers)
        static = any(layer['unified_mode'] in ('static', 'original') for layer in layers)
        # Groups color individual maps by within-group source counts only for the layers whose
        # per-group context is presence; a group map colored by value gets its colors from the
        # layer's own colormap instead.
        grouped_presence = grouped and any(
            layer['category_mode'] == 'membership' for layer in layers
        )

        def _dedup(items: List[str]) -> List[str]:
            seen: Set[str] = set()
            return [item for item in items if not (item in seen or seen.add(item))]

        if has_categories:
            subset_names = set(categories)
            self._check_requested_subset(
                draw_individual_files, "Individual maps", subset_names, subset_subject
            )
            self._check_requested_subset(
                draw_grid, "Individual maps in grids", subset_names, subset_subject
            )
            draw_files_categories = (
                _dedup(list(categories)) if draw_individual_files is True
                else [] if draw_individual_files is False
                else _dedup(list(draw_individual_files))
            )
            draw_grid_categories = (
                _dedup(list(categories)) if draw_grid is True
                else [] if draw_grid is False
                else _dedup(list(draw_grid))
            )
            draw_categories = _dedup(draw_files_categories + draw_grid_categories)
            # Only a category drawn on its own maps has its name joined onto a path, so only those
            # names have to be usable as directory names. A run that draws no individual maps or
            # grids, or that asks for a subset of them, never writes the other categories' names
            # anywhere: they are summarized on the 'unified' map by color alone. The check comes
            # after the subset checks so that a name that is not a category at all is reported as
            # unrecognized rather than as unusable.
            self._check_category_names(draw_categories, category_noun)
        else:
            draw_files_categories = []
            draw_grid_categories = []
            draw_categories = []
        draw_category_maps = has_categories and (
            draw_individual_files is not False or draw_grid is not False
        )

        # Limits on, and a center for, the scale the individual maps share have nothing to act on
        # when this run draws no individual map, just as the group-map coloring options have nothing
        # to act on then. Either one accepted and quietly dropped would look, from the output,
        # exactly like a limit that did nothing because no value crossed it, so say which it is.
        if not draw_category_maps:
            for model_key, flag_suffix, subject_phrase, remedy_phrase in (
                (
                    'category_value_limits', 'category-value-limits',
                    'Limits were given for', 'bound'
                ),
                (
                    'category_value_center', 'category-value-center',
                    'A center was given for', 'center'
                )
            ):
                affected_layers = [
                    layer for layer in layers if layer.get(model_key) is not None
                ]
                if not affected_layers:
                    continue
                flags = ', '.join(
                    f"'--{layer['element_type']}-{flag_suffix}'" for layer in affected_layers
                )
                raise ConfigError(
                    f"{subject_phrase} the color scale shared by the maps of the individual "
                    f"samples or groups ({flags}), but this run draws no such maps, so that scale "
                    f"is never drawn and this would go nowhere. Ask for those maps with "
                    f"'--draw-individual-files' and/or '--draw-grid', or {remedy_phrase} the scale "
                    f"of the 'unified' map instead, with "
                    f"'--reaction-{flag_suffix.replace('category-', '')}'/"
                    f"'--compound-{flag_suffix.replace('category-', '')}'."
                )

        # Colors given per category name are checked against the categories the run has, and against
        # what this run would actually do with them, before anything is drawn. A layer whose contexts
        # are all colored some other way would leave them unused, which is worth an error rather than
        # a map that quietly ignores the file it was handed.
        colored_layers = [layer for layer in layers if layer.get('category_colors') is not None]
        if colored_layers:
            if not has_categories:
                flags = ', '.join(f"'{layer['category_colors_flag']}'" for layer in colored_layers)
                raise ConfigError(
                    f"A color was given per category ({flags}), but this run has no categories to "
                    f"color: there is a single source of data, so there is one map rather than one "
                    f"per sample, database, genome, or group. Set the color of the single map with "
                    f"'--reaction-color'/'--compound-color' instead."
                )
            self._resolve_category_colors(colored_layers, categories, category_noun)
            group_ramps = grouped and (
                grouped_membership['group_colormap'] == GROUP_COLORMAP_FROM_CATEGORY
            )
            for layer in colored_layers:
                # Where the colors land: the combinations of the 'unified' map's membership scale,
                # the one color of each category's own map when the categories are not grouped, and
                # the ramp of each group's own maps when they are and that ramp was asked for.
                if layer['unified_mode'] == 'membership':
                    continue
                if layer['category_mode'] == 'membership' and draw_category_maps and (
                    group_ramps if grouped else True
                ):
                    continue
                grouped_clause = (
                    f" Grouped, a group's own maps are colored by the number of the group's "
                    f"{membership_singular}s containing an element, which takes its colors from "
                    f"'--group-colormap': ask for '--group-colormap "
                    f"{GROUP_COLORMAP_FROM_CATEGORY}' to build each group's scale from that group's "
                    f"own color."
                ) if grouped else ""
                raise ConfigError(
                    f"The colors file given to '{layer['category_colors_flag']}' gives a color to "
                    f"each {category_noun}, which colors an element by exactly which "
                    f"{category_noun}s contain it, and colors each {category_noun}'s own map. "
                    f"Nothing in this run is colored either way for the {layer['element_type']} "
                    f"layer, so the file would go unused.{grouped_clause} Please remove it, or "
                    f"color this layer by presence."
                )

        if static and grouped:
            # The individual group maps are the exception to "static overrides dynamic": a single
            # color cannot distinguish a group's sources, so those maps fall back to within-group
            # source counts. Only say so when such maps are actually requested.
            if draw_individual_files is not False or draw_grid is not False:
                group_map_clause = (
                    f" The individual group maps are an exception: since one color cannot "
                    f"distinguish a group's {membership_singular}s, they are colored by the number "
                    f"of {membership_singular}s in each group containing an element, styled by "
                    f"'--group-colormap'/'--group-reverse-overlay' rather than by the static color."
                )
            else:
                group_map_clause = ""
            self.run.warning(
                f"Groups were provided, but these will be ignored for the 'unified' map, since a "
                f"static color (or the reference map's own colors) was set: dynamic coloring based "
                f"on membership in groups is overridden by static coloring based on "
                f"presence/absence in any {membership_singular}.{group_map_clause}"
            )

        unified_dir = os.path.join(output_dir, UNIFIED_SUBDIR)
        pathway_numbers = self._find_maps(unified_dir, patterns=pathway_numbers)
        filesnpaths.gen_output_directory(output_dir, progress=self.progress, run=self.run)

        drawn: Dict[Literal['unified', 'individual', 'grid'], Dict] = {
            'unified': {}, 'individual': {}, 'grid': {}
        }

        # Finalize the coloring model of each layer. A context colored quantitatively needs a value
        # range; a single parse of each map collects, per layer, the unified values and (shared
        # across categories) the per-category values of whichever contexts are quantitative, using
        # the same extraction and aggregation as drawing.
        norm_layers = [
            layer for layer in layers
            if 'quantitative' in (layer['unified_mode'], layer['category_mode'])
        ]
        if norm_layers:
            self.progress.new("Computing the range of values across maps")
            for layer in norm_layers:
                layer['_unified_vals'] = []
                layer['_category_vals'] = []
            for pathway_number in pathway_numbers:
                self.progress.update(pathway_number)
                pathway = self._get_pathway(pathway_number)
                for layer in norm_layers:
                    use_reaction = layer['use_reaction_attribute']
                    unified_valued = layer['unified_mode'] == 'quantitative'
                    per_category = (
                        has_categories and layer['category_mode'] == 'quantitative'
                        and layer['category_values'] is not None
                    )
                    for entry in self._find_element_entries(
                        pathway, use_reaction, layer['accessions']
                    ):
                        if unified_valued:
                            value = self._reduce_entry_value(
                                entry, layer['unified_values'], layer['aggregate'],
                                use_reaction_attribute=use_reaction
                            )
                            if value is not None:
                                layer['_unified_vals'].append(value)
                        if per_category:
                            for category_name in categories:
                                category_value = self._reduce_entry_value(
                                    entry, layer['category_values'][category_name],
                                    layer['aggregate'], use_reaction_attribute=use_reaction
                                )
                                if category_value is not None:
                                    layer['_category_vals'].append(category_value)
            self.progress.end()
            for layer in norm_layers:
                # No values at all means no element of any drawn map has one, so the layer colors
                # nothing and gets no colorbar: either its accessions are absent from these maps, or
                # its aggregation was undefined everywhere (the standard deviation of a single
                # value, say). Say so rather than leaving a blank map to be puzzled over.
                if not layer['_unified_vals'] and not layer['_category_vals']:
                    self.run.warning(
                        f"Nothing on the maps could be colored by the values of the "
                        f"'{layer['colorbar_label']}' column of the {layer['element_type']} layer, "
                        f"so no color scale was drawn for it. Either none of the drawn maps "
                        f"contains its accessions, or the aggregation reducing them is undefined "
                        f"for every map element — the standard deviation of a single value, for "
                        f"instance."
                    )
                if layer['unified_mode'] == 'quantitative':
                    norm, vmin, vmax, clamped_low, clamped_high = self._make_quantitative_norm(
                        layer['_unified_vals'], layer.get('value_limits'),
                        f"--{layer['element_type']}-value-limits",
                        center=layer.get('value_center'),
                        center_flag=f"--{layer['element_type']}-value-center",
                        subject=(
                            f"the {layer['element_type']} layer of the 'unified' map"
                            if has_categories else f"the {layer['element_type']} layer"
                        )
                    )
                    layer['_unified_norm'] = norm
                    layer['_unified_range'] = (vmin, vmax)
                    layer['_unified_clamped'] = (clamped_low, clamped_high)
                    layer['_unified_center'] = layer.get('value_center')
                if layer['category_mode'] != 'quantitative':
                    continue
                if has_categories and layer['category_values'] is not None:
                    norm, vmin, vmax, clamped_low, clamped_high = self._make_quantitative_norm(
                        layer['_category_vals'], layer.get('category_value_limits'),
                        f"--{layer['element_type']}-category-value-limits",
                        center=layer.get('category_value_center'),
                        center_flag=f"--{layer['element_type']}-category-value-center",
                        subject=(
                            f"the {layer['element_type']} layer of the maps of the individual "
                            f"{category_noun}s"
                        )
                    )
                    layer['_category_norm'] = norm
                    layer['_category_range'] = (vmin, vmax)
                    layer['_category_clamped'] = (clamped_low, clamped_high)
                    layer['_category_center'] = layer.get('category_value_center')
                else:
                    # A layer without a category dimension is constant across the category maps.
                    layer['_category_norm'] = layer['_unified_norm']
                    layer['_category_range'] = layer['_unified_range']
                    layer['_category_clamped'] = layer['_unified_clamped']
                    layer['_category_center'] = layer['_unified_center']

        # A context colored by membership needs its by-count/by-membership color scheme over the
        # categories. Only the 'unified' map draws such a scale (and its colorbar) from the
        # categories themselves; the presence coloring of grouped individual maps counts each
        # group's own sources, precomputed below.
        membership_layers = [layer for layer in layers if layer['unified_mode'] == 'membership']
        group_membership_layers = [
            layer for layer in layers if layer['category_mode'] == 'membership'
        ]
        # Where a count scale is asked to stop at the highest count in the data, that count has to
        # be found first, which takes one parse of the maps about to be drawn.
        group_observed_counts: Dict[str, int] = {}
        if count_scale_max == 'observed':
            scaled_group_layers = group_membership_layers if (
                draw_category_maps and grouped_presence
            ) else []
            if membership_layers or scaled_group_layers:
                group_observed_counts = self._find_observed_counts(
                    membership_layers, scaled_group_layers, pathway_numbers, group_sources,
                    group_threshold
                )
        if membership_layers:
            self.progress.new("Setting map colors")
            self.progress.update("...")
            for layer in membership_layers:
                layer['_count_scale_top'] = self._resolve_count_scale_top(
                    count_scale_max, layer.get('_observed_count', 0), len(categories)
                )
                layer['_colors'] = self._membership_layer_colors(
                    layer, categories, layer['_count_scale_top']
                )
            self.progress.end()

        # Each group's individual maps are colored by how many of that group's own sources contain
        # an element, on a scale of its own. The colors are resolved here, before anything is
        # checked or drawn, so that they can be checked against the colors a map reserves along with
        # every other color of the run; the per-layer specs that use them are built further below,
        # once the within-group membership they need has been narrowed to each group.
        group_color_priorities: Dict[str, List[Tuple[str, float]]] = {}
        group_scale_tops: Dict[str, int] = {}
        group_cmaps: Dict[str, mcolors.Colormap] = {}
        group_colormap_scheme = None
        if grouped_presence:
            group_color_priorities, group_scale_tops, group_cmaps, group_colormap_scheme = (
                self._group_map_colors(
                    grouped_membership, layers, draw_categories, group_observed_counts,
                    count_scale_max, unified_plural, category_noun
                )
            )

        def _reaction_derived(layer, mode, cmap_key='cmap'):
            # How a reaction layer derives compound colors on a reaction-only global/overview map.
            # The derived colors come from the same scale as the reactions they are derived from, so
            # the caller names the context's colormap: the two can differ ('category_cmap').
            if layer['element_type'] != 'reaction':
                return None
            if mode == 'quantitative':
                cmap = layer[cmap_key]
                return ('average', cmap.reversed() if layer['reverse_overlay'] else cmap)
            return ('high', None)

        def _unified_spec(layer):
            mode = layer['unified_mode']
            if mode == 'quantitative':
                return {
                    'element_type': layer['element_type'],
                    'use_reaction_attribute': layer['use_reaction_attribute'],
                    'entry_keys': layer['unified_values'],
                    'colorer': self._quantitative_colorer(
                        layer['unified_values'], layer['_unified_norm'], layer['cmap'],
                        layer['reverse_overlay'], layer['aggregate'],
                        layer['use_reaction_attribute'], center=layer['_unified_center']
                    ),
                    'derived_compound': _reaction_derived(layer, mode)
                }
            if mode == 'membership':
                _, color_priorities, category_combos, _ = layer['_colors']
                return {
                    'element_type': layer['element_type'],
                    'use_reaction_attribute': layer['use_reaction_attribute'],
                    'entry_keys': layer['membership'],
                    'colorer': self._membership_colorer(
                        layer['membership'], color_priorities, category_combos, group_sources,
                        group_threshold, layer['use_reaction_attribute']
                    ),
                    'derived_compound': _reaction_derived(layer, mode)
                }
            # 'static' pools accessions across all sources; 'single' has its own accession set.
            entry_keys = set(layer['membership']) if mode == 'static' else layer['accessions']
            return {
                'element_type': layer['element_type'],
                'use_reaction_attribute': layer['use_reaction_attribute'],
                'entry_keys': entry_keys,
                'colorer': self._single_color_colorer(layer['color_hexcode']),
                'derived_compound': _reaction_derived(layer, mode)
            }

        def _category_spec(layer, category):
            mode = layer['category_mode']
            if mode == 'quantitative':
                values = (
                    layer['unified_values'] if layer['category_values'] is None
                    else layer['category_values'][category]
                )
                return {
                    'element_type': layer['element_type'],
                    'use_reaction_attribute': layer['use_reaction_attribute'],
                    'entry_keys': values,
                    'colorer': self._quantitative_colorer(
                        values, layer['_category_norm'], layer['category_cmap'],
                        layer['reverse_overlay'], layer['aggregate'],
                        layer['use_reaction_attribute'], center=layer['_category_center']
                    ),
                    'derived_compound': _reaction_derived(layer, mode, 'category_cmap')
                }
            if mode == 'single':
                return {
                    'element_type': layer['element_type'],
                    'use_reaction_attribute': layer['use_reaction_attribute'],
                    'entry_keys': layer['accessions'],
                    'colorer': self._single_color_colorer(layer['color_hexcode']),
                    'derived_compound': _reaction_derived(layer, mode, 'category_cmap')
                }
            # Ungrouped membership/static: an individual source's map colors that source's elements
            # in one fixed color (the grouped case is precomputed below). With a color per category,
            # that color is the category's own, so a panel of a grid says which category it is and
            # matches the band it takes on the 'unified' map's membership scale.
            accessions = layer['source_accessions'].get(category, set())
            category_colors = layer.get('category_colors')
            color_hexcode = layer['color_hexcode'] if category_colors is None else (
                category_colors[category]
            )
            return {
                'element_type': layer['element_type'],
                'use_reaction_attribute': layer['use_reaction_attribute'],
                'entry_keys': accessions,
                'colorer': self._single_color_colorer(color_hexcode),
                'derived_compound': _reaction_derived(layer, mode, 'category_cmap')
            }

        self._check_reserved_colors(
            layers, pathway_numbers, group_color_priorities=group_color_priorities,
            group_element_types={layer['element_type'] for layer in group_membership_layers}
        )

        # Per-layer colorbars for the unified map (layer-prefixed so two layers do not collide),
        # plus a shared colorbar for the category maps of each layer colored quantitatively there.
        # Each context is keyed by its own mode, so a layer summarized by presence in the unified
        # map gets a presence colorbar there and a continuous one for its category maps.
        for layer in layers:
            if layer['unified_mode'] == 'quantitative':
                vmin, vmax = layer['_unified_range']
                if vmin is not None:
                    clamped_low, clamped_high = layer['_unified_clamped']
                    self._draw_quantitative_colorbar(
                        layer['cmap'], vmin, vmax,
                        os.path.join(output_dir, f"colorbar_{layer['name']}.pdf"),
                        layer['colorbar_label'],
                        clamped_low=clamped_low, clamped_high=clamped_high,
                        center=layer['_unified_center']
                    )
            elif layer['unified_mode'] == 'membership':
                scheme, color_priorities, category_combos, presence_cmap = layer['_colors']
                count_scale_top = layer['_count_scale_top']
                colorbar_path = os.path.join(output_dir, f"colorbar_{layer['name']}.pdf")
                count_label = 'group count' if grouped else membership_count_label
                if scheme == 'by_count_continuous':
                    # The same color per count that 'by_count' assigns, but shown as a gradient
                    # across the colormap from the lowest count to the highest rather than as one
                    # labeled band per count, which is what frees it from needing a distinct color
                    # for each.
                    self._draw_quantitative_colorbar(
                        presence_cmap, 1, count_scale_top, colorbar_path, count_label,
                        integer_ticks=True
                    )
                else:
                    if scheme == 'by_count':
                        labels = range(1, count_scale_top + 1)
                        label = count_label
                    else:
                        labels = [', '.join(combo) for combo in category_combos]
                        label = 'groups' if grouped else membership_members_label
                    self.colorbar_drawer.draw_discrete(
                        [color for color, _ in color_priorities],
                        colorbar_path,
                        color_labels=labels,
                        label=label
                    )
            if (
                draw_category_maps and layer['category_mode'] == 'quantitative'
                and layer['category_values'] is not None
            ):
                vmin, vmax = layer['_category_range']
                if vmin is not None:
                    clamped_low, clamped_high = layer['_category_clamped']
                    self._draw_quantitative_colorbar(
                        layer['category_cmap'], vmin, vmax,
                        os.path.join(
                            output_dir, f"colorbar_{layer['name']}_{colorbar_category_suffix}.pdf"
                        ),
                        layer['colorbar_label'],
                        clamped_low=clamped_low, clamped_high=clamped_high,
                        center=layer['_category_center']
                    )

        # Draw the unified map (the single map when there are no categories).
        self.progress.new(
            f"Drawing 'unified' map incorporating data from all {unified_plural}"
            if has_categories else "Drawing map"
        )
        unified_specs = None if original_run else [_unified_spec(layer) for layer in layers]
        for pathway_number in pathway_numbers:
            self.progress.update(pathway_number)
            if original_run:
                drawn['unified'][pathway_number] = self._draw_map_kos_original_color(
                    pathway_number, set(layers[0]['membership']), unified_dir,
                    draw_map_lacking_data=draw_maps_lacking_data
                )
            else:
                drawn['unified'][pathway_number] = self._draw_map_elements(
                    pathway_number, unified_specs, unified_dir,
                    draw_map_lacking_data=draw_maps_lacking_data
                )
        self.progress.end()

        if not draw_categories:
            count = sum(drawn['unified'].values()) if drawn['unified'] else 0
            self.run.info("Number of maps drawn", count)
            return drawn

        # For grouped membership individual maps, precompute each group's within-group membership,
        # narrowing the layer's membership to the group's own sources so that a group's map counts
        # only them. The colors these counts take were resolved above, with every other color of the
        # run.
        group_layer_membership: Dict[str, Tuple] = {}
        if grouped_presence:
            for group in draw_categories:
                specs = []
                for layer in group_membership_layers:
                    inner_membership = {}
                    for accession, sources in layer['membership'].items():
                        in_group = [s for s in sources if source_group.get(s) == group]
                        if in_group:
                            inner_membership[accession] = in_group
                    specs.append({
                        'element_type': layer['element_type'],
                        'use_reaction_attribute': layer['use_reaction_attribute'],
                        'entry_keys': inner_membership,
                        'colorer': self._membership_colorer(
                            inner_membership, group_color_priorities[group], None, None, None,
                            layer['use_reaction_attribute']
                        ),
                        'derived_compound': _reaction_derived(layer, 'membership')
                    })
                group_layer_membership[group] = (
                    specs, group_color_priorities[group], group_scale_tops[group]
                )

        for category in draw_categories:
            drawn_category: Dict[str, bool] = {}
            self.progress.new(f"Drawing maps for {category_noun} '{category}'")
            self.progress.update("...")
            progress = self.progress
            self.progress = terminal.Progress(verbose=False)
            run = self.run
            self.run = terminal.Run(verbose=False)
            category_output_dir = os.path.join(output_dir, INDIVIDUAL_SUBDIR, category)
            filesnpaths.gen_output_directory(
                category_output_dir, progress=self.progress, run=self.run
            )

            if grouped_presence:
                group_specs, category_color_priorities, category_scale_top = (
                    group_layer_membership[category]
                )
                colorbar_path = os.path.join(category_output_dir, CATEGORY_COLORBAR_BASENAME)
                if group_colormap_scheme == 'by_count_continuous':
                    # The same color per count that the discrete bands assign, shown as a gradient
                    # across the ramp from the lowest count to the highest rather than as one labeled
                    # band per count, which is what frees it from needing a distinct color for each.
                    self._draw_quantitative_colorbar(
                        group_cmaps[category], 1, category_scale_top, colorbar_path,
                        membership_count_label, integer_ticks=True
                    )
                else:
                    self.colorbar_drawer.draw_discrete(
                        [color for color, _ in category_color_priorities],
                        colorbar_path,
                        color_labels=range(1, category_scale_top + 1),
                        label=membership_count_label
                    )
                # Grouped membership specs are precomputed; a layer colored by value or a single
                # color in its per-group context colors its own map and rides along on the same map.
                specs = group_specs + [
                    _category_spec(layer, category) for layer in layers
                    if layer['category_mode'] != 'membership'
                ]
                for pathway_number in pathway_numbers:
                    drawn_category[pathway_number] = self._draw_map_elements(
                        pathway_number, specs, category_output_dir,
                        draw_map_lacking_data=draw_maps_lacking_data
                    )
            elif original_run:
                for pathway_number in pathway_numbers:
                    drawn_category[pathway_number] = self._draw_map_kos_original_color(
                        pathway_number, layers[0]['source_accessions'].get(category, set()),
                        category_output_dir, draw_map_lacking_data=draw_maps_lacking_data
                    )
            else:
                specs = [_category_spec(layer, category) for layer in layers]
                for pathway_number in pathway_numbers:
                    drawn_category[pathway_number] = self._draw_map_elements(
                        pathway_number, specs, category_output_dir,
                        draw_map_lacking_data=draw_maps_lacking_data
                    )

            self.progress = progress
            self.run = run
            self.progress.end()
            drawn['individual'][category] = drawn_category

        if draw_grid is not False:
            grid_group_color_priorities = None
            grid_group_scale_tops = None
            if grouped_presence:
                grid_group_color_priorities = {
                    group: priorities
                    for group, (_, priorities, _) in group_layer_membership.items()
                }
                grid_group_scale_tops = {
                    group: scale_top
                    for group, (_, _, scale_top) in group_layer_membership.items()
                }
            self._draw_map_grids(
                pathway_numbers,
                draw_categories,
                draw_grid_categories,
                draw_files_categories,
                output_dir,
                drawn,
                group_scale_tops=grid_group_scale_tops,
                group_color_priorities=grid_group_color_priorities,
                group_cmaps=group_cmaps,
                group_colormap_scheme=group_colormap_scheme,
                check_maps_lacking_kos=not draw_maps_lacking_data,
                source_type=grid_source_type if grid_source_type is not None else category_noun
            )

        # The individual maps are gathered only once nothing else stands to change them: drawing the
        # grids deletes the maps that were drawn as grid panels alone, along with the directories of
        # the categories that were never asked for individually.
        collated_count = 0
        if self.collate_files_by_map and draw_files_categories:
            self.progress.new(f"Gathering the maps of each {category_noun}")
            self.progress.update("...")
            collated_count = self._collate_maps_by_map(output_dir, draw_files_categories)
            self.progress.end()

        count = sum(drawn['unified'].values()) if drawn['unified'] else 0
        self.run.info(
            f"Number of 'unified' maps drawn incorporating data from all {unified_plural}", count
        )
        if draw_individual_files:
            count = sum(
                sum(d.values()) if d else 0 for d in drawn['individual'].values()
            ) if drawn['individual'] else 0
            self.run.info(f"Number of maps drawn for individual {category_noun}s", count)
            if self.collate_files_by_map:
                self.run.info(
                    f"Number of maps gathered across individual {category_noun}s", collated_count
                )
        count = sum(drawn['grid'].values()) if drawn['grid'] else 0
        self.run.info("Number of map grids drawn", count)

        return drawn

    def _check_reserved_colors(
        self,
        layers: List[dict],
        pathway_numbers: List[str],
        group_color_priorities: Dict[str, List[Tuple[str, float]]] = None,
        group_element_types: Set[str] = None
    ) -> None:
        """
        Check that no layer color is one that a map reserves for its unidentified elements.

        The reactions and compounds a map does not highlight are recolored so that highlighted ones
        stand out ('kgml.Pathway.set_color_priority'), which puts those colors out of bounds for a
        layer: an element colored one of them could not be told apart from the background, so 'kgml'
        refuses to draw such a map. It refuses per map, once drawing is under way, which would leave
        the output half written -- so the same clash is caught here first, across every class of map
        about to be drawn, before any file is created.

        Parameters
        ==========
        layers : List[dict]
            The layer models, whose colors are checked. Layers drawn in the reference map's own
            colors are skipped, having no colors of their own.

        pathway_numbers : List[str]
            The maps about to be drawn, whose classes decide which colors are reserved: global maps
            recolor to gray, overview maps to black reactions and white compounds, and standard maps
            to white.

        group_color_priorities : Dict[str, List[Tuple[str, float]]], None
            The colors of each group's individual maps ('_group_map_colors'), which belong to no one
            layer: one scale colors the within-group source counts of every layer whose per-group
            context is presence. A ramp built from a group's own color runs towards white by
            construction, so it can reach a reserved color even where a layer's own colors do not.

        group_element_types : Set[str], None
            The element types the group scale colors, deciding which reserved colors apply to it.
        """
        reserved: Dict[str, Set[str]] = {'reaction': set(), 'compound': set()}
        for pathway_number in pathway_numbers:
            is_global = re.match(GLOBAL_MAP_ID_PATTERN, pathway_number) is not None
            is_overview = re.match(OVERVIEW_MAP_ID_PATTERN, pathway_number) is not None
            recolor_colors = kgml.reserved_recolor_colors(
                'g' if is_global else 'w', is_overview
            )
            reserved['reaction'].add(kgml.canonical_color(recolor_colors['ortholog']))
            reserved['compound'].add(kgml.canonical_color(recolor_colors['compound']))

        # One scale colors the within-group source counts of every layer whose per-group context is
        # presence, so it is out of bounds for the reserved colors of all of their element types.
        group_reserved = set().union(
            *(reserved[element_type] for element_type in group_element_types)
        ) if group_element_types else set()
        for group, color_priorities in (group_color_priorities or {}).items():
            clashing = sorted(
                {
                    color for color in map(
                        kgml.canonical_color, (color for color, _ in color_priorities)
                    ) if color in group_reserved
                }
            )
            if not clashing:
                continue
            self.progress.end()
            raise ConfigError(
                f"The individual maps of group '{group}' would be colored "
                f"{'colors' if len(clashing) > 1 else 'a color'} that the pathway maps keep for "
                f"their own unidentified elements: {', '.join(clashing)}. Anvi'o recolors the "
                f"reactions and compounds a map does not highlight — gray on global maps, and "
                f"black reactions with white compounds elsewhere — so a highlighted element in one "
                f"of those colors would be invisible. These maps color the number of the group's "
                f"own sources containing an element, styled by '--group-colormap': choose a "
                f"colormap that does not reach these colors, or, where the scale is a ramp running "
                f"from a pale tint to the group's own color ('--group-colormap "
                f"{GROUP_COLORMAP_FROM_CATEGORY}'), start that ramp further from white with the "
                f"two limits that option also takes."
            )

        for layer in layers:
            if layer['unified_mode'] == 'original':
                continue
            # Every color the layer can stage: sampled from its colormap at the values it will color
            # by, taken from the scale it colors presence by, or its one fixed color.
            staged: Set[str] = set()
            for norm_key, values_key, cmap_key in (
                ('_unified_norm', '_unified_vals', 'cmap'),
                ('_category_norm', '_category_vals', 'category_cmap')
            ):
                if norm_key not in layer:
                    continue
                norm = layer[norm_key]
                cmap = layer[cmap_key]
                for value in layer[values_key]:
                    fraction = 1.0 if norm is None else float(norm(value))
                    staged.add(mcolors.rgb2hex(cmap(fraction)))
            if '_colors' in layer:
                staged.update(color for color, _ in layer['_colors'][1])
            if layer.get('category_colors') is not None:
                staged.update(layer['category_colors'].values())
            if layer.get('color_hexcode') is not None:
                staged.add(layer['color_hexcode'])

            clashing = sorted(
                {
                    color for color in map(kgml.canonical_color, staged)
                    if color in reserved[layer['element_type']]
                }
            )
            if not clashing:
                continue
            self.progress.end()
            raise ConfigError(
                f"The {layer['element_type']} layer would be colored "
                f"{'colors' if len(clashing) > 1 else 'a color'} that the pathway maps keep for "
                f"their own unidentified elements: {', '.join(clashing)}. Anvi'o recolors the "
                f"reactions and compounds a map does not highlight -- gray on global maps, and "
                f"black reactions with white compounds elsewhere -- so a highlighted element in "
                f"one of those colors would be invisible. Choose a different color, or a colormap "
                f"that does not reach these colors: grayscale colormaps and those running to pure "
                f"white or black, such as 'Greys', 'hot' and 'bone', all do."
            )

    @staticmethod
    def _resolve_count_scale_top(count_scale_max: Union[str, int], observed: int, total: int) -> int:
        """
        Resolve the count a color scale runs up to from the '--count-scale-max' choice.

        There are three choices. 'observed' stops the scale at the highest count anything on the
        drawn maps actually has, so that the colors spread over the counts that occur. 'total' runs
        it to every category there is, which keeps the scale the same however few of them the data
        reaches, and so is comparable between runs. A number pins the scale, which is how separate
        figures are given one scale when their data differs.

        Parameters
        ==========
        count_scale_max : Union[str, int]
            'observed', 'total', or the count to stop at.

        observed : int
            The highest count anything on the drawn maps has, 0 if nothing was counted.

        total : int
            How many categories there are in all.

        Returns
        =======
        int
            The count the scale runs up to.
        """
        if count_scale_max == 'observed':
            # An observed maximum of 0 means nothing on the drawn maps was counted at all, and a
            # scale of no counts cannot be drawn, so it runs to every category instead. The layer
            # colors nothing either way, so which of the two it is never shows on a map.
            top = observed if observed > 0 else total
        elif count_scale_max == 'total':
            top = total
        else:
            top = int(count_scale_max)

        # Whatever the choice, a scale reaches at least 1, since an element colored by its count is
        # in at least one category.
        return max(top, 1)

    def _find_observed_counts(
        self,
        unified_layers: List[dict],
        group_layers: List[dict],
        pathway_numbers: Iterable[str],
        group_sources: Union[Dict[str, List[str]], None],
        group_threshold: Union[float, None]
    ) -> Dict[str, int]:
        """
        Find the highest count that presence coloring reaches, in one parse of the drawn maps.

        Where a count scale stops at the highest count in the data ('--count-scale-max observed'),
        that count is not known until every element of every map about to be drawn has been counted.
        The pass that computes a value scale's range works the same way and for the same reason: one
        scale has to serve every map for the colors on them to be comparable.

        Parameters
        ==========
        unified_layers : List[dict]
            Layers colored by count or membership on the 'unified' map. Each has the highest count
            it reaches there recorded as '_observed_count'.

        group_layers : List[dict]
            Layers colored by within-group source counts on grouped individual maps. What they reach
            is returned rather than recorded, since one scale serves every layer of a group's map.

        pathway_numbers : Iterable[str]
            Numeric IDs of the maps that will be drawn.

        Returns
        =======
        Dict[str, int]
            The highest within-group source count of each group, empty when there are no groups or no
            layer is colored by them.

        Notes
        =====
        'group_sources'/'group_threshold' group the sources, as elsewhere.
        """
        for layer in unified_layers:
            layer['_observed_count'] = 0
        categorizers = [
            self._membership_categorizer(
                layer['membership'], group_sources, group_threshold,
                layer['use_reaction_attribute']
            )
            for layer in unified_layers
        ]

        source_group: Dict[str, str] = {}
        if group_sources is not None:
            for group, sources in group_sources.items():
                for source in sources:
                    source_group[source] = group
        group_counts: Dict[str, int] = (
            {group: 0 for group in group_sources} if group_layers and source_group else {}
        )

        self.progress.new("Finding the highest count across maps")
        # Groups are counted inside the pass over the maps rather than by looping over the groups
        # themselves, so that a run with many groups reads each map once.
        for pathway_number in pathway_numbers:
            self.progress.update(pathway_number)
            pathway = self._get_pathway(pathway_number)
            for layer, categorize in zip(unified_layers, categorizers):
                for entry in self._find_element_entries(
                    pathway, layer['use_reaction_attribute'], layer['membership']
                ):
                    categories = categorize(entry)
                    if categories is not None and len(categories) > layer['_observed_count']:
                        layer['_observed_count'] = len(categories)
            for layer in group_layers if group_counts else ():
                for entry in self._find_element_entries(
                    pathway, layer['use_reaction_attribute'], layer['membership']
                ):
                    # A group's own map counts the group's sources containing an element, whatever
                    # the group threshold, which only decides whether the GROUP counts as containing
                    # it on the 'unified' map.
                    within: Dict[str, int] = {}
                    for source in self._entry_sources(
                        entry, layer['membership'], layer['use_reaction_attribute']
                    ):
                        group = source_group.get(source)
                        if group is not None:
                            within[group] = within.get(group, 0) + 1
                    for group, count in within.items():
                        if count > group_counts[group]:
                            group_counts[group] = count
        self.progress.end()

        return group_counts

    def _group_map_colors(
        self,
        grouped_membership: dict,
        layers: List[dict],
        draw_categories: List[str],
        group_observed_counts: Dict[str, int],
        count_scale_max: Union[str, int],
        unified_plural: str,
        category_noun: str
    ) -> Tuple[
        Dict[str, List[Tuple[str, float]]],
        Dict[str, int],
        Dict[str, mcolors.Colormap],
        Literal['by_count', 'by_count_continuous']
    ]:
        """
        Resolve the colors of each group's individual maps, which count the group's own sources.

        Each group's map colors elements by how many of that group's own sources contain them. Which
        color a count takes never depends on the scheme resolved here: the colors are always sampled
        at even fractions of the group's ramp, so the scheme decides only how the scale is DRAWN, in
        discrete bands of one color per count ('by_count') or as a gradient from the lowest count to
        the highest ('by_count_continuous'). Unlike the schemes of the 'unified' map, which each
        layer's summary chooses for itself ('_membership_layer_colors'), this one is chosen for the
        whole run by the group style, and rightly so: both layers of a group's map share one set of
        colors and one colorbar, so there is no per-layer scale for a per-layer scheme to describe.

        Where the colors run short of the counts, the gradient is what a discrete band per count
        cannot be, so the scheme falls back to it and says so, exactly as the 'unified' map's count
        scale does. Asking for 'by_count' outright keeps the bands, and neighboring counts then share
        a color, as the warning below reports; every element is still colored by its own count.

        The colors come from a named Matplotlib colormap, shared by every group, or, when
        'GROUP_COLORMAP_FROM_CATEGORY' is asked for, from a ramp per group running from a pale tint
        to that group's own color ('tint_hexcode'). The ramp binds the identity channel to the
        magnitude channel: every panel of a grid is built the same way and so stays comparable,
        while its hue says which group it is. A named colormap is the default because one engineered
        for magnitude reads better as magnitude than a color chosen to identify a group does.

        Parameters
        ==========
        grouped_membership : dict
            The grouping record, carrying 'group_sources' and the group style ('group_colormap'/
            'group_colormap_limits'/'group_reverse_overlay'/'group_colormap_scheme').

        layers : List[dict]
            Every layer of the run, of which those given a color per category say what color a ramp
            runs to. Since one style covers every layer's group maps, two layers coloring one group
            differently is refused.

        draw_categories : List[str]
            The groups whose own maps are drawn.

        group_observed_counts : Dict[str, int]
            The highest within-group source count found on the drawn maps, per group, or empty when
            the scale does not stop at what was observed.

        count_scale_max : Union[str, int]
            Where each group's count scale stops ('_resolve_count_scale_top').

        unified_plural : str
            The plural of what a group's sources are, for messages, e.g. 'samples'.

        category_noun : str
            What one category is called, for messages, e.g. 'sample group'.

        Returns
        =======
        Tuple[Dict[str, List[Tuple[str, float]]], Dict[str, int], Dict[str, mcolors.Colormap], str]
            Per group: the (color_hexcode, priority) pairs in ascending order of count, as
            '_membership_colorer' looks them up; the count its scale runs up to; and the colormap
            its colors were sampled from, which its colorbar spans, left empty for a scale drawn in
            discrete bands, which spans no colormap. Then the scheme ('by_count'/
            'by_count_continuous') that every group's colorbar is drawn by.
        """
        group_sources = grouped_membership['group_sources']
        group_colormap = grouped_membership['group_colormap']
        colormap_limits = grouped_membership['group_colormap_limits']
        group_reverse_overlay = grouped_membership['group_reverse_overlay']
        requested_scheme = grouped_membership['group_colormap_scheme']

        from_category = group_colormap == GROUP_COLORMAP_FROM_CATEGORY
        group_colors: Dict[str, str] = {}
        group_cmap = None
        if from_category:
            # One ramp per group needs one color per group, and one style covers every layer's group
            # maps, so the layers that have colors must agree about them.
            colored_layers = [layer for layer in layers if layer.get('category_colors')]
            if not colored_layers:
                self.progress.end()
                raise ConfigError(
                    f"The group colormap was given as '{GROUP_COLORMAP_FROM_CATEGORY}', which "
                    f"colors each group's own maps by a ramp running from a pale tint to that "
                    f"group's own color, but no color was given for any group. Give one per group "
                    f"with the category colors option of a layer colored by presence, or name a "
                    f"Matplotlib colormap for the group maps instead."
                )
            # A ramp is built only for the groups whose own maps are drawn, so they are the only ones
            # the layers must agree about: a group the run has but draws no map for has no ramp to
            # disagree over. Every one of them is guaranteed a color by '_resolve_category_colors'.
            for layer in colored_layers:
                for group in draw_categories:
                    color = layer['category_colors'][group]
                    if group_colors.setdefault(group, color) != color:
                        self.progress.end()
                        raise ConfigError(
                            f"The group colormap was given as '{GROUP_COLORMAP_FROM_CATEGORY}', so "
                            f"each group's own maps are colored by a ramp running to that group's "
                            f"color. One ramp per group styles every layer's group maps at once, "
                            f"but the reaction and compound layers were given different colors for "
                            f"the {category_noun} '{group}': '{group_colors[group]}' and "
                            f"'{color}'. Please give that {category_noun} one color in both files "
                            f"— pointing both options at a single file is the simplest way — or "
                            f"name a Matplotlib colormap for the group maps instead."
                        )
            if colormap_limits is None:
                colormap_limits = DEFAULT_GROUP_TINT_SPAN
            if not 0.0 <= colormap_limits[0] <= colormap_limits[1] <= 1.0:
                self.progress.end()
                raise ConfigError(
                    f"The two limits of a group colormap of '{GROUP_COLORMAP_FROM_CATEGORY}' are "
                    f"how far from white towards a group's own color its ramp starts and stops, so "
                    f"they must lie between 0.0 and 1.0 with the smaller one first. These do not: "
                    f"{colormap_limits[0]}, {colormap_limits[1]}."
                )
        else:
            if isinstance(group_colormap, str):
                group_cmap = self._get_colormap(group_colormap)
            else:
                group_cmap = group_colormap
            if group_cmap.name in qualitative_colormaps + repeating_colormaps:
                self.run.warning(
                    f"The group colormap, '{group_cmap.name}', that was provided to color "
                    f"individual group maps is not especially useful for displaying the count of "
                    f"{unified_plural}. We recommend a sequential colormap like 'plasma' instead."
                )
            if colormap_limits is None:
                colormap_limits = DEFAULT_GROUP_COLORMAP_LIMITS
            group_cmap = self._trim_colormap(group_cmap, colormap_limits)

        # A group's scale stops where '--count-scale-max' says, exactly as the 'unified' map's does,
        # so that a group of many sources whose elements are in only a few of them is not drawn in
        # one shade at the bottom of the scale. Every group's top is resolved before any colors are
        # built, because the scheme below is decided for the whole run and so needs all of them.
        group_scale_tops: Dict[str, int] = {
            group: self._resolve_count_scale_top(
                count_scale_max, group_observed_counts.get(group, 0), len(group_sources[group])
            )
            for group in draw_categories
        }

        group_color_priorities: Dict[str, List[Tuple[str, float]]] = {}
        group_distinct_colors: Dict[str, int] = {}
        for group in draw_categories:
            group_scale_top = group_scale_tops[group]
            if group_scale_top == 1:
                sample_points = np.linspace(1, 1, 1)
            else:
                sample_points = np.linspace(0, 1, group_scale_top)
            if from_category:
                lower_limit, upper_limit = colormap_limits
                # A group of one source is drawn in that group's color exactly, the top of its ramp,
                # since there is no range of counts for a ramp to spread over.
                colors = [
                    tint_hexcode(
                        group_colors[group],
                        lower_limit + (upper_limit - lower_limit) * sample_point
                    )
                    for sample_point in sample_points
                ]
            else:
                colors = [
                    mcolors.rgb2hex(group_cmap(sample_point)) for sample_point in sample_points
                ]
            # In ascending order of count, as '_membership_colorer' looks them up.
            group_color_priorities[group] = [
                (color, 1 - sample_point if group_reverse_overlay else sample_point)
                for color, sample_point in zip(colors, sample_points)
            ]
            # Rounding to 8-bit color means the supply can fall short of what the ramp holds, so
            # what matters is how many distinct colors actually came out of it.
            group_distinct_colors[group] = len(set(colors))

        # Two things stop a band per count from being drawn: a ramp with too few DISTINCT colors to
        # tell the bands apart, and more bands than a colorbar can label in type large enough to
        # read ('MAX_DISCRETE_COUNT_BANDS'). The scheme is settled once for the whole run rather
        # than per group, so that every panel of a grid carries the same kind of bar; the group that
        # runs furthest past a limit decides it, since a scheme good for that one is good for the
        # rest. A group with nothing of its own on the drawn maps colors nothing, and under the
        # default its scale top is not a count that occurs but a stand-in for the scale it cannot
        # otherwise have ('_resolve_count_scale_top'). Letting that stand-in decide would hand every
        # other group a gradient on account of a group with nothing to show, and name the empty one
        # as the reason, so it sits the decision out. A top that was asked for outright still
        # counts: its bands are really drawn to it, however little lands on them.
        deciding_groups = [
            group for group in draw_categories
            if count_scale_max != 'observed' or group_observed_counts.get(group, 0) > 0
        ] or list(draw_categories)
        short_groups = [
            group for group in deciding_groups
            if group_distinct_colors[group] < group_scale_tops[group]
        ]
        crowded_groups = [
            group for group in deciding_groups
            if group_scale_tops[group] > MAX_DISCRETE_COUNT_BANDS
        ]
        scheme = requested_scheme
        if scheme is None:
            # Nobody asked for either scheme, so the bands are kept while they can be drawn and read,
            # and the gradient takes over when they cannot.
            scheme = 'by_count_continuous' if short_groups or crowded_groups else 'by_count'

        def _short_group_clause(group: str) -> str:
            return (
                f"The ramp running to the color of {category_noun} '{group}' could supply only "
                f"{group_distinct_colors[group]} colors that can be told apart, fewer than the "
                f"{group_scale_tops[group]} counts its scale runs over. Widening the ramp's span "
                f"with the two limits of the group colormap option gives it more room."
            ) if from_category else (
                f"The group colormap could supply only {group_distinct_colors[group]} colors that "
                f"can be told apart, fewer than the {group_scale_tops[group]} counts the scale of "
                f"{category_noun} '{group}' runs over."
            )

        def _crowded_group_clause(group: str) -> str:
            return (
                f"The scale of {category_noun} '{group}' runs over {group_scale_tops[group]} "
                f"counts, and a discrete colorbar labels one band per count, which is more labels "
                f"than it can set in type large enough to read past approximately "
                f"{MAX_DISCRETE_COUNT_BANDS} of them."
            )

        # Too few colors is the harder limit of the two, so it is the one reported when both apply:
        # a ramp that cannot fill the bands cannot fill fewer of them either.
        def _worst_group_clause() -> str:
            if short_groups:
                return _short_group_clause(max(
                    short_groups,
                    key=lambda group: group_scale_tops[group] - group_distinct_colors[group]
                ))
            return _crowded_group_clause(max(crowded_groups, key=group_scale_tops.get))

        if scheme == 'by_count_continuous':
            # Running past either limit is what the gradient is for, so it is worth reporting only
            # when the gradient was not asked for: the bands were the default, and this is why they
            # were not drawn. Only the colorbar changes, since the colors a count takes are sampled
            # from fractions of the ramp either way.
            if (short_groups or crowded_groups) and requested_scheme is None:
                self.run.warning(
                    f"{_worst_group_clause()} The count on individual {category_noun} maps is "
                    f"therefore drawn on a CONTINUOUS color scale rather than in discrete bands of "
                    f"one color per count. Each colorbar is a gradient running from a count of 1 "
                    f"to the top of that {category_noun}'s scale, on which a color reads as a "
                    f"position along that range rather than as an exact count. Ask for this scale "
                    f"explicitly with '{GROUP_SCHEME_OPTIONS['by_count_continuous']}', or ask for "
                    f"'{GROUP_SCHEME_OPTIONS['by_count']}' to insist on the discrete bands and be "
                    f"told when they cannot be drawn. One scheme covers every {category_noun}'s "
                    f"maps, so the {category_noun} that runs furthest past a limit is the one "
                    f"reported here."
                )
        else:
            # The bands were asked for outright, so each group that runs past a limit is reported
            # once, by the limit that says the most about it.
            for group in short_groups:
                self.run.warning(
                    f"{_short_group_clause(group)} Neighboring counts therefore share a color on "
                    f"that {category_noun}'s individual maps, and its colorbar labels more bands "
                    f"than it has distinct colors. Every element is still colored by the count of "
                    f"the {unified_plural} containing it. Drawing the count on a continuous color "
                    f"scale instead, with '{GROUP_SCHEME_OPTIONS['by_count_continuous']}', needs "
                    f"no distinct color per count."
                )
            for group in crowded_groups:
                if group in short_groups:
                    continue
                self.run.warning(
                    f"{_crowded_group_clause(group)} The bands are drawn as asked, and every "
                    f"element is still colored by the count of the {unified_plural} containing it, "
                    f"but that {category_noun}'s colorbar labels will be very small. Drawing the "
                    f"count on a continuous color scale instead, with "
                    f"'{GROUP_SCHEME_OPTIONS['by_count_continuous']}', labels the range rather "
                    f"than every count."
                )

        # A continuous colorbar spans the colormap its colors were sampled from. A named colormap,
        # trimmed to its limits above, is that colormap already and is shared by every group; a ramp
        # built from a group's own color has to be made into one. Neither is built for a scale drawn
        # in discrete bands, which spans no colormap: its bar is the colors themselves.
        group_cmaps: Dict[str, mcolors.Colormap] = {}
        if scheme == 'by_count_continuous':
            for group in draw_categories:
                if from_category:
                    lower_limit, upper_limit = colormap_limits
                    group_cmaps[group] = mcolors.LinearSegmentedColormap.from_list(
                        f'tint({group_colors[group]},{lower_limit:.2f},{upper_limit:.2f})',
                        [
                            tint_hexcode(
                                group_colors[group],
                                lower_limit + (upper_limit - lower_limit) * fraction
                            )
                            for fraction in np.linspace(0, 1, GROUP_RAMP_COLORMAP_SIZE)
                        ]
                    )
                else:
                    group_cmaps[group] = group_cmap

        return group_color_priorities, group_scale_tops, group_cmaps, scheme

    def _membership_layer_category_colors(
        self,
        layer: dict,
        categories: List[str],
        count_scale_top: int
    ) -> Tuple[str, List[Tuple[str, float]], List[Tuple[str]], None]:
        """
        Resolve a membership layer's colors from a color given per category.

        This is the 'category_colors' branch of '_membership_layer_colors', which see: it returns
        the same record, and the drawing and colorbar code cannot tell where the colors came from.
        Each category alone is drawn in its own color, and a combination of categories in the blend
        of their colors ('blend_hexcodes') unless the colors file overrides that combination.
        Drawing priority follows the order of the combinations exactly as it does when they are
        sampled from a colormap, so an element in more categories is still drawn over one in fewer.

        Parameters
        ==========
        layer : dict
            The layer model, carrying 'category_colors', 'category_combo_colors', the flag they came
            from, and 'reverse_overlay'.

        categories : List[str]
            The categories whose membership colors the layer, in color-assignment order.

        count_scale_top : int
            The count a color scale would run up to, named in the message that points at coloring by
            count instead.

        Returns
        =======
        Tuple[str, List[Tuple[str, float]], List[Tuple[str]], None]
            'by_membership', the (color_hexcode, priority) pairs in combination order, the
            combinations themselves, and None in place of the colormap the colors did not come from.
        """
        category_colors: Dict[str, str] = layer['category_colors']
        combo_colors: Dict[Tuple[str, ...], str] = layer['category_combo_colors']
        flag = layer['category_colors_flag']
        scheme_options = layer.get('scheme_options', PRESENCE_SCHEME_OPTIONS)
        reverse_overlay = layer.get('reverse_overlay', False)

        # A count says how many categories contain an element, not which, so there is no category
        # whose color it could take. Refused rather than ignored, since the two were asked for
        # together and only one of them can be honored.
        colormap_scheme = layer.get('colormap_scheme')
        if colormap_scheme is not None and colormap_scheme != 'by_membership':
            self.progress.end()
            raise ConfigError(
                f"The {layer['element_type']} layer was given a color for each category by "
                f"'{flag}', and asked to be colored by count with "
                f"'{scheme_options[colormap_scheme]}'. A count is how many categories contain an "
                f"element rather than which ones, so it cannot take a category's own color: the "
                f"two requests cannot both be honored. Either drop '{flag}' and color the counts "
                f"from a colormap, or ask for '{scheme_options['by_membership']}' so that the "
                f"colors given per category are what an element is colored by."
            )

        # Coloring by membership needs a color per combination of the categories, so the count of
        # them doubles with each category. The ceiling is checked before the combinations are
        # enumerated, so that a large number of categories cannot spend gigabytes on its way to the
        # same refusal.
        combo_count = 2 ** len(categories) - 1
        if combo_count > MAX_CATEGORY_COLOR_COMBOS:
            self.progress.end()
            combos_shown = f'{combo_count}' if combo_count <= 10 ** 6 else f'{combo_count:.1e}'
            raise ConfigError(
                f"Coloring the {layer['element_type']} layer by membership needs a distinct color "
                f"for every combination of the {len(categories)} categories, of which there are "
                f"{combos_shown}. Anvi'o blends the colors given by '{flag}' to derive them, and "
                f"refuses past {MAX_CATEGORY_COLOR_COMBOS} combinations, well beyond the handful "
                f"that any color scale could tell apart. Color by count instead with "
                f"'{scheme_options['by_count']}' and drop '{flag}', since a count takes its colors "
                f"from a colormap and the two cannot both be honored; it needs just "
                f"{count_scale_top} colors. Alternatively, reduce the number of categories, for "
                f"example by grouping them."
            )

        category_combos: List[Tuple[str]] = []
        for category_count in range(1, len(categories) + 1):
            category_combos += list(combinations(categories, category_count))

        color_priorities: List[Tuple[str, float]] = []
        blended: List[Tuple[str]] = []
        for combo_index, combo in enumerate(category_combos):
            if len(combo) == 1:
                color = category_colors[combo[0]]
            else:
                override = combo_colors.get(tuple(sorted(combo)))
                if override is None:
                    color = blend_hexcodes(category_colors[category] for category in combo)
                    blended.append(combo)
                else:
                    color = override
            # The position in the combination order is both the color's place on the scale and its
            # drawing priority, which 'reverse_overlay' inverts. Counting from the last combination
            # rather than subtracting from 1 keeps every priority non-negative, as a Pathway
            # requires.
            color_priorities.append((
                color,
                1 - combo_index / len(category_combos) if reverse_overlay
                else (combo_index + 1) / len(category_combos)
            ))

        # A discrete colorbar labels one band per combination, so two combinations sharing a color
        # would leave the bar with more labels than colors a reader can tell apart. Blending is what
        # usually causes this — two blends can land on one color, and a blend can land on a
        # category's own color — so the message points at the rows that would fix it.
        distinct = len({color for color, _ in color_priorities})
        if distinct != len(category_combos):
            self.progress.end()
            blend_clause = (
                f" {len(blended)} of the {len(category_combos)} combinations took a color blended "
                f"from their members' colors; a row of the file naming a combination, with its "
                f"names separated by '{CATEGORY_COMBO_SEPARATOR}', sets that combination's color "
                f"directly instead."
            ) if blended else ""
            raise ConfigError(
                f"The colors of the {layer['element_type']} layer of this map could not be "
                f"assigned. Coloring by membership needs a distinct color for every combination of "
                f"the {len(categories)} categories, which is {len(category_combos)} of them, and "
                f"the colors from '{flag}' came to only {distinct}.{blend_clause} Alternatively, "
                f"color by count with '{scheme_options['by_count']}' and drop '{flag}', since a "
                f"count takes its colors from a colormap and the two cannot both be honored; it "
                f"needs just {count_scale_top} colors."
            )

        return 'by_membership', color_priorities, category_combos, None

    def _membership_layer_colors(
        self,
        layer: dict,
        categories: List[str],
        count_scale_top: int = None
    ) -> Tuple[str, List[Tuple[str, float]], Union[List[Tuple[str]], None], mcolors.Colormap]:
        """
        Resolve a membership layer's coloring scheme, colors and priorities, and category combos.

        Resolves the by-count/by-membership colormap logic for one layer. 'categories' are the
        sources (or groups) whose count/membership colors the layer, and 'count_scale_top' is the
        count the color scale runs up to, defaulting to all of them ('_resolve_count_scale_top').
        Stopping the scale where the data does spreads the colors over the counts that occur instead
        of over counts nothing reaches, which for sparse data is the difference between a readable
        map and one colored in a single shade; the cost is that the scale then depends on what was
        drawn.

        The colors come back in the order the scheme assigns them — by ascending count, or by the
        order of the category combinations — rather than keyed by color, because a continuous count
        scale ('by_count_continuous') deliberately gives the same color to neighboring counts
        wherever the colormap has no distinguishable color left for each one. The two count schemes
        sample the same colormap and differ only in how the scale is drawn, and therefore in whether
        the colors they assign must be distinguishable: a discrete colorbar labels one band per
        count, so a count without its own color would leave the bar mislabeled, whereas a gradient
        from the first count to the last stays honest however many counts share a color. For the
        sequential colormap a count scale calls for, the colors the two assign are identical; only a
        qualitative colormap makes them differ, since 'by_count' samples it at whole positions while
        a gradient has to span a range.

        A layer given a color per category ('category_colors') is colored by membership from those
        colors instead of from a colormap: each category alone takes its own color, and a
        combination of them takes the blend of their colors unless a row of the colors file
        overrides it. Since a count identifies no category, the count schemes have nothing to take
        from such a file and are refused for it.

        Returns
        =======
        Tuple[str, List[Tuple[str, float]], Union[List[Tuple[str]], None], matplotlib.colors.Colormap]
            The scheme ('by_count'/'by_count_continuous'/'by_membership'), the (color_hexcode,
            priority) pairs in assignment order, the list of category combinations (for
            by-membership) or None (for the count schemes), and the trimmed colormap the colors were
            sampled from, which a continuous count scale's colorbar spans (None when the colors came
            from a colors file rather than a colormap, which only by-membership allows).
        """
        colormap = layer.get('colormap', True)
        colormap_scheme = layer.get('colormap_scheme')
        category_colors = layer.get('category_colors')
        # A message that tells the reader how to ask for a scheme has to name the option this
        # layer's input actually takes, since the option the other input takes is refused here.
        scheme_options = layer.get('scheme_options', PRESENCE_SCHEME_OPTIONS)
        reverse_overlay = layer.get('reverse_overlay', False)
        # Only the count schemes have a scale that can stop early: coloring by membership needs a
        # color for every combination of the categories however few of them the data reaches.
        if count_scale_top is None:
            count_scale_top = len(categories)

        if category_colors is not None:
            return self._membership_layer_category_colors(layer, categories, count_scale_top)

        if colormap_scheme is not None:
            scheme = colormap_scheme
        else:
            scheme = 'by_membership' if len(categories) < 4 else 'by_count'

        colormap_limits = layer.get('colormap_limits')
        if colormap is True:
            if scheme == 'by_membership':
                cmap = self._get_colormap('tab10')
                colormap_limits = (0.0, 1.0) if colormap_limits is None else colormap_limits
            else:
                cmap = self._get_colormap('plasma_r')
                colormap_limits = (0.1, 0.9) if colormap_limits is None else colormap_limits
        elif isinstance(colormap, str):
            cmap = self._get_colormap(colormap)
            colormap_limits = (0.0, 1.0) if colormap_limits is None else colormap_limits
        elif isinstance(colormap, mcolors.Colormap):
            cmap = colormap
            colormap_limits = (0.0, 1.0) if colormap_limits is None else colormap_limits
        else:
            raise AssertionError

        # A qualitative colormap is sampled at whole positions rather than at fractions of its range.
        # A continuous count scale is the exception: its colorbar is a gradient across the fraction of
        # the colormap in use, so its colors have to be sampled from that same fraction to be the
        # ones the colorbar shows.
        colormap_name = cmap.name
        qualitative = colormap_name in qualitative_colormaps + repeating_colormaps
        cmap = self._trim_colormap(cmap, colormap_limits)

        # Coloring by membership needs a color per combination of the categories, so the count of
        # them doubles with each category. A colormap holds at most 'cmap.N' colors, so the check
        # below would fail anyway; making it before the combinations are enumerated keeps a large
        # number of categories from spending gigabytes on its way to the same error.
        if scheme == 'by_membership' and 2 ** len(categories) - 1 > cmap.N:
            self.progress.end()
            # Written out in full up to a million and in scientific notation above it, since past a
            # few dozen categories the exact number runs to hundreds of digits that say nothing.
            combos = 2 ** len(categories) - 1
            combos_shown = f'{combos}' if combos <= 10 ** 6 else f'{combos:.1e}'
            raise ConfigError(
                f"Coloring the {layer['element_type']} layer by membership needs a distinct color "
                f"for every combination of the {len(categories)} categories, of which there are "
                f"{combos_shown}, and its colormap holds only {cmap.N}. Color by count instead "
                f"with '{scheme_options['by_count']}', which needs just {count_scale_top} colors, "
                f"or give a colormap with more colors. Note that no color scale can distinguish "
                f"combinations of more than a handful of categories."
            )

        def _count_colors(in_order: bool) -> List[Tuple[str, float]]:
            # The color and drawing priority of each count, in ascending order of count, from a
            # count of 1 up to the top of the scale.
            if count_scale_top == 1:
                sample_points = range(1, 2) if in_order else np.linspace(1, 1, 1)
            else:
                sample_points = range(count_scale_top) if in_order else np.linspace(
                    0, 1, count_scale_top
                )
            # A sample point is both a position in the colormap and a drawing priority, which
            # 'reverse_overlay' inverts. Inverting it as '1 - point' suits the fractions but goes
            # negative at the whole positions an 'in_order' colormap is sampled at (1 - 2 = -1), and
            # a Pathway requires non-negative priorities. Counting back from the last point descends
            # without going negative, and for fractions that point IS 1, so one expression does
            # both.
            last_point = max(sample_points)
            return [
                (
                    mcolors.rgb2hex(cmap(sample_point)),
                    (last_point - sample_point) if reverse_overlay else sample_point
                )
                for sample_point in sample_points
            ]

        category_combos = None
        if scheme == 'by_membership':
            category_combos = []
            for category_count in range(1, len(categories) + 1):
                category_combos += list(combinations(categories, category_count))
            if qualitative:
                sample_points = range(len(category_combos))
            else:
                sample_points = np.linspace(0, 1, len(category_combos))
            color_priorities = [
                (
                    mcolors.rgb2hex(cmap(sample_point)),
                    1 - sample_point / cmap.N if reverse_overlay else (sample_point + 1) / cmap.N
                )
                for sample_point in sample_points
            ]
        else:
            color_priorities = _count_colors(qualitative and scheme == 'by_count')

        # A discrete colorbar labels one band per count or per membership combination, so a colormap
        # that cannot supply a DISTINCT color for each of them leaves the bar with more labels than
        # colors a reader can tell apart. Rounding to 8-bit color means the supply can fall short of
        # 'cmap.N' too, so what matters is how many distinct colors actually came out, not how many
        # the colormap claims to hold.
        needed = len(category_combos) if scheme == 'by_membership' else count_scale_top
        distinct = len({color for color, _ in color_priorities})
        # Two things stop a band per count from being drawn: a colormap with too few distinct colors
        # to tell the bands apart, and more bands than a colorbar can label in type large enough to
        # read. Either is reason for the gradient, and the message says which it was.
        short_colors = scheme == 'by_count' and distinct != needed
        too_many_bands = scheme == 'by_count' and count_scale_top > MAX_DISCRETE_COUNT_BANDS
        if (short_colors or too_many_bands) and colormap_scheme is None:
            # Nobody asked for the discrete bands: the scheme was chosen from the number of
            # categories, and there are more of them than bands can carry, so the count is drawn as
            # a gradient instead. The colors are resampled from fractions of the colormap's range,
            # which is what that gradient spans.
            scheme = 'by_count_continuous'
            color_priorities = _count_colors(False)
            # Too few colors is the harder limit of the two, so it is the one reported when both
            # apply: a colormap that cannot fill the bands cannot fill fewer of them either.
            reason = (
                f"its colormap could supply only {distinct} colors that can be told apart"
            ) if short_colors else (
                f"a discrete colorbar labels one band per count, which is more labels than it can "
                f"set in type large enough to read past approximately {MAX_DISCRETE_COUNT_BANDS} "
                f"of them"
            )
            self.run.warning(
                f"The {layer['element_type']} layer is colored by counts running up to "
                f"{count_scale_top}, and {reason}, so that count is drawn on a CONTINUOUS color "
                f"scale rather than in discrete bands of one color per count. The colorbar is a "
                f"gradient running from a count of 1 to a count of {count_scale_top}, on which a "
                f"color reads as a position along that range rather than as an exact count. Ask "
                f"for this scale explicitly with '{scheme_options['by_count_continuous']}', or ask "
                f"for '{scheme_options['by_count']}' to insist on the discrete bands and be told "
                f"when they cannot be drawn. Each layer decides this for itself, so with a "
                f"reaction layer and a compound layer you may see this twice, once per layer.",
                progress=self.progress
            )
        elif too_many_bands and not short_colors:
            # The bands were asked for outright. Unlike a colormap short of colors, which cannot
            # draw them at all, this only makes them hard to read, so they are drawn as asked — so
            # a colormap that is ALSO short of colors is left to the refusal below, which is what
            # actually happens next, rather than being told here that the bands are drawn.
            self.run.warning(
                f"The {layer['element_type']} layer was asked to color counts running up to "
                f"{count_scale_top} in discrete bands, which is more bands than a colorbar can "
                f"label in type large enough to read: it holds about {MAX_DISCRETE_COUNT_BANDS}. "
                f"The bands are drawn as asked, and every element is colored by its own count, but "
                f"the colorbar's labels will be very small. Drawing the count on a continuous "
                f"scale instead, with '{scheme_options['by_count_continuous']}', labels the range "
                f"rather than every count.",
                progress=self.progress
            )
        if scheme == 'by_count_continuous' and qualitative:
            self.run.warning(
                f"The colormap, '{colormap_name}', that colors the {layer['element_type']} layer "
                f"by count is qualitative rather than sequential, which makes a continuous color "
                f"scale difficult to interpret. We recommend a sequential colormap like 'plasma' "
                f"instead.",
                progress=self.progress
            )
        if scheme != 'by_count_continuous' and distinct != needed:
            self.progress.end()
            if scheme == 'by_membership':
                advice = (
                    f"Coloring by membership needs a distinct color for every combination of the "
                    f"{len(categories)} categories, which is {needed} of them, and the colormap "
                    f"supplied only {distinct}. Color by count instead with "
                    f"'{scheme_options['by_count']}', which needs just {count_scale_top} colors, "
                    f"or give a colormap with more distinct colors."
                )
            else:
                advice = (
                    f"Coloring by count in discrete bands needs a distinct color for each count "
                    f"up to {needed}, and the colormap supplied only {distinct}. Color by "
                    f"count on a continuous scale instead, which needs no distinct color per count "
                    f"and is asked for with '{scheme_options['by_count_continuous']}'; "
                    f"alternatively, reduce the number of categories, for example by grouping "
                    f"them, or give a colormap with more distinct colors."
                )
            raise ConfigError(
                f"The colors of the {layer['element_type']} layer of this map could not be "
                f"assigned: {advice} Note that a color scale can hold at most a few hundred "
                f"distinguishable colors in any case, so a very large number of categories cannot "
                f"be told apart by color even where it can be drawn."
            )

        return scheme, color_priorities, category_combos, cmap

    def _map_element_membership(
        self,
        layers: List[dict],
        all_sources: List[str],
        source_type: Literal['contigs database', 'pangenome', 'sample'],
        source_group: Dict[str, str] = None,
        group_sources: Dict[str, List[str]] = None,
        group_threshold: float = None,
        pathway_numbers: Iterable[str] = None,
        draw_individual_files: Union[Iterable[str], bool] = False,
        draw_grid: Union[Iterable[str], bool] = False,
        group_colormap: Union[str, mcolors.Colormap] = 'plasma_r',
        group_colormap_limits: Tuple[float, float] = None,
        group_reverse_overlay: bool = False,
        group_colormap_scheme: Literal['by_count', 'by_count_continuous'] = None,
        count_scale_max: Union[str, int] = 'observed',
        output_dir: str = None,
        draw_maps_lacking_data: bool = False
    ) -> Dict[Literal['unified', 'individual', 'grid'], Dict]:
        """
        Adapt presence/absence membership layers to the unified '_map_elements' engine.

        Each layer colors reaction and/or compound elements by presence/absence across sources
        (contigs databases, pangenome genomes, or samples) or groups of sources, by source/group
        count or membership, or a single static color (or the reference map's original colors) when
        a layer's colormap is False. This method classifies each layer's mode and translates the
        source-type terminology and grouping into a call to '_map_elements'. The single-layer
        contigs-db and pan-db paths route through here (reaction layer only).

        Parameters
        ==========
        layers : List[dict]
            One or two layer descriptors, reaction before compound. Each has: 'name'
            ('reactions'/'compounds'), 'element_type', 'use_reaction_attribute', 'membership'
            ({accession: [sources]}), 'source_accessions' ({source: set of accessions}, for
            individual ungrouped maps), 'color_hexcode' (single color for individual ungrouped
            maps), the colormap options 'colormap'/'colormap_limits'/'colormap_scheme'/
            'reverse_overlay', and the colors given per category, if any, as 'category_colors'/
            'category_combo_colors' with the 'category_colors_flag' they came from.

        all_sources : List[str]
            Names of all sources (samples, contigs databases, or genomes), in color-assignment
            order.

        source_type : Literal['contigs database', 'pangenome', 'sample']
            The kind of source, selecting terminology for messages and colorbar labels.

        Notes
        =====
        'source_group'/'group_sources'/'group_threshold' group the sources; the remaining parameters
        mirror the other engines.
        """
        grouped = group_sources is not None

        # Terminology used in messages and colorbar labels, keyed by source type.
        singular, plural, count_label, members_label, group_phrase = {
            'contigs database': (
                'contigs database', 'contigs databases', 'database count', 'databases',
                'contigs database group'
            ),
            'pangenome': ('genome', 'genomes', 'genome count', 'genomes', 'group'),
            'sample': ('sample', 'samples', 'sample count', 'samples', 'sample group')
        }[source_type]

        # A layer whose 'colormap' is False is colored statically (a single fixed color, or the
        # reference map's original colors) by presence/absence in any source, rather than by
        # source/group count or membership. Only the single-layer db/pan paths use it; the
        # draw-kegg-pathways layers are never static.
        models = []
        for layer in layers:
            if layer.get('colormap') is False:
                mode = 'original' if layer['color_hexcode'] == 'original' else 'static'
            else:
                mode = 'membership'
            # These inputs take their presence scheme from '--presence-colormap-scheme', which is
            # what a message about the scheme should name for them.
            models.append({**layer, 'mode': mode, 'scheme_options': PRESENCE_SCHEME_OPTIONS})

        grouped_membership = None
        if grouped:
            grouped_membership = {
                'source_group': source_group,
                'group_sources': group_sources,
                'group_threshold': group_threshold,
                'group_colormap': group_colormap,
                'group_colormap_limits': group_colormap_limits,
                'group_reverse_overlay': group_reverse_overlay,
                'group_colormap_scheme': group_colormap_scheme
            }

        return self._map_elements(
            models,
            output_dir,
            pathway_numbers=pathway_numbers,
            categories=list(group_sources) if grouped else list(all_sources),
            # The categories are groups when the sources are grouped, so the noun that names one has
            # to follow: it labels the per-category messages and the subdirectory each category is
            # drawn into.
            category_noun=group_phrase if grouped else singular,
            subset_subject=f"{singular} groups" if grouped else plural,
            unified_plural=plural,
            membership_count_label=count_label,
            membership_members_label=members_label,
            membership_singular=singular,
            grouped_membership=grouped_membership,
            count_scale_max=count_scale_max,
            draw_individual_files=draw_individual_files,
            draw_grid=draw_grid,
            draw_maps_lacking_data=draw_maps_lacking_data
        )

    def _find_element_entries(
        self,
        pathway: kgml.Pathway,
        use_reaction_attribute: bool,
        values: Dict[str, float]
    ) -> List[kgml.Entry]:
        """
        Find the entries a layer colors: those with accessions among 'values'.

        For KO and compound layers, accessions are matched against 'Entry.name' via
        'get_entries(kegg_ids=...)'. For a reaction-by-R-number layer, ortholog entries are matched
        by the reaction IDs in 'Entry.reaction', which the KEGG-ID index does not cover.

        Parameters
        ==========
        pathway : kgml.Pathway
            The pathway to search.

        use_reaction_attribute : bool
            If True, match reaction IDs from 'Entry.reaction'; otherwise match KO/compound IDs from
            'Entry.name'.

        values : Dict[str, float]
            Keys are the accessions of interest.

        Returns
        =======
        List[kgml.Entry]
            The matching entries.
        """
        if use_reaction_attribute:
            return [
                entry for entry in pathway.get_entries(entry_type='ortholog')
                if any(
                    reaction_id in values
                    for reaction_id in self._get_entry_kegg_ids(entry, use_reaction_attribute=True)
                )
            ]
        return pathway.get_entries(kegg_ids=values)

    def _stage_element_color(
        self,
        pathway: kgml.Pathway,
        entry: kgml.Entry,
        element_type: Literal['reaction', 'compound'],
        color_hexcode: str,
        priority: float,
        color_priority: dict
    ) -> None:
        """
        Set an entry's graphics colors for a layer and register them in 'color_priority'.

        Reaction (ortholog) entries are lines in global and overview maps and boxes or lines in
        standard maps; compound entries are circles. The registered '(fgcolor, bgcolor) -> priority'
        keeps the entry from being treated as unprioritized and recolored to the background.

        Parameters
        ==========
        pathway : kgml.Pathway
            The pathway being colored.

        entry : kgml.Entry
            The entry to color.

        element_type : Literal['reaction', 'compound']
            Which kind of element (and thus which graphics/coloring convention) this is.

        color_hexcode : str
            The color for the element.

        priority : float
            Drawing-order priority (higher renders on top).

        color_priority : dict
            The accumulating '{entry_type: {graphics_type: {(fg, bg): priority}}}' dictionary,
            shared across a map's layers and passed once to 'set_color_priority'.
        """
        if element_type == 'compound':
            for uuid in entry.children['graphics']:
                graphics: kgml.Graphics = pathway.uuid_element_lookup[uuid]
                # Compounds are circles. On a few maps some compounds are rectangles that are zeroed
                # out of the base image; those cannot be colored and are skipped (the caller warns).
                if graphics.type != 'circle':
                    continue
                if pathway.is_global_map:
                    graphics.fgcolor = color_hexcode
                    graphics.bgcolor = color_hexcode
                    colors = (color_hexcode, color_hexcode)
                else:
                    graphics.fgcolor = '#000000'
                    graphics.bgcolor = color_hexcode
                    colors = ('#000000', color_hexcode)
                color_priority.setdefault(
                    'compound', {}
                ).setdefault('circle', {})[colors] = priority
            return

        for uuid in entry.children['graphics']:
            graphics: kgml.Graphics = pathway.uuid_element_lookup[uuid]
            if pathway.is_global_map:
                assert graphics.type == 'line'
                graphics.fgcolor = color_hexcode
                graphics.bgcolor = '#FFFFFF'
                graphics_type = 'line'
                colors = (color_hexcode, '#FFFFFF')
            elif pathway.is_overview_map:
                assert graphics.type == 'line'
                graphics.fgcolor = color_hexcode
                graphics.bgcolor = '#FFFFFF'
                graphics.width = 5.0
                graphics_type = 'line'
                colors = (color_hexcode, '#FFFFFF')
            else:
                if graphics.type == 'rectangle':
                    graphics.fgcolor = '#000000'
                    graphics.bgcolor = color_hexcode
                    graphics_type = 'rectangle'
                    colors = ('#000000', color_hexcode)
                elif graphics.type == 'line':
                    graphics.fgcolor = color_hexcode
                    graphics.bgcolor = '#FFFFFF'
                    graphics.width = 5.0
                    graphics_type = 'line'
                    colors = (color_hexcode, '#FFFFFF')
                else:
                    self.progress.end()
                    raise ConfigError(
                        f"Reaction elements are expected to be drawn as a rectangle or a line, but "
                        f"an ortholog entry of KEGG pathway map {pathway.number} has a graphics "
                        f"element of type '{graphics.type}', which anvi'o cannot color."
                    )
            color_priority.setdefault(
                'ortholog', {}
            ).setdefault(graphics_type, {})[colors] = priority

    def _warn_unrenderable_compounds(
        self,
        pathway: kgml.Pathway,
        compound_accessions: Iterable[str]
    ) -> None:
        """
        Warn about supplied compounds that cannot be colored on a map.

        Compounds are colored via their circle Graphics. On a handful of maps (00121, 00621,
        01052, 01054), some compounds are drawn only as rectangles, which are zeroed out of the
        base image (see '_zero_out_compound_rectangles') so they do not obscure the chemical
        structure drawings there. Such rectangles cannot be colored, so a supplied compound present
        on the map only as rectangles would be invisible; this warns rather than dropping it
        silently.

        Parameters
        ==========
        pathway : kgml.Pathway
            The map being drawn.

        compound_accessions : Iterable[str]
            Compound accessions supplied for the compound layer of this map.
        """
        unrenderable: List[str] = []
        for accession in compound_accessions:
            entries = pathway.get_entries(kegg_ids=[accession])
            if not entries:
                continue
            if not any(
                pathway.uuid_element_lookup[uuid].type == 'circle'
                for entry in entries
                for uuid in entry.children['graphics']
            ):
                unrenderable.append(accession)

        if not unrenderable:
            return

        self.run.warning(
            f"On KEGG pathway map {pathway.number}, the following supplied compound(s) could not "
            f"be colored because they are represented on this map only as rectangles rather than "
            f"circles: {', '.join(sorted(unrenderable))}. Anvi'o zeroes out these compound "
            f"rectangles so they do not obscure the chemical structure drawings on the few maps "
            f"where they occur (such as 00121, 00621, 01052, and 01054), so these particular "
            f"compounds cannot be shown in color. Everything else was drawn as usual.",
            progress=self.progress
        )

    def _quantitative_colorer(
        self,
        values: Dict[str, float],
        norm: Union[mcolors.Normalize, None],
        cmap: mcolors.Colormap,
        reverse_overlay: bool,
        aggregate,
        use_reaction_attribute: bool,
        center: Union[float, None] = None
    ):
        """
        Build a colorer that colors an Entry by the continuous value of its accessions.

        See '_draw_map_elements' for the colorer contract. An Entry's value is the 'aggregate' of
        its accessions' values; None (no accession has a value) leaves the Entry uncolored. The
        color is 'cmap' sampled at the normalized value, and the priority is that fraction, or its
        complement under 'reverse_overlay', which 'clip=True' on the norm keeps in [0, 1]. A
        degenerate range (no norm) leaves every element at the top of the colormap, except on a
        centered scale, where the one value the range collapsed to is the center itself and so takes
        the middle color, the very color the centering was asked for.
        """
        def colorer(entry: kgml.Entry) -> Union[Tuple[str, float], None]:
            value = self._reduce_entry_value(
                entry, values, aggregate, use_reaction_attribute=use_reaction_attribute
            )
            if value is None:
                return None
            if norm is None:
                fraction = 0.5 if center is not None else 1.0
            else:
                fraction = float(norm(value))
            priority = (1.0 - fraction) if reverse_overlay else fraction
            return mcolors.rgb2hex(cmap(fraction)), priority
        return colorer

    def _single_color_colorer(self, color_hexcode: str):
        """
        Build a colorer that colors every matching Entry a single fixed color at priority 1.0.

        See '_draw_map_elements' for the colorer contract. Used wherever a layer's elements all take
        one fixed color: a 'single' layer, the pooled 'unified' map of a 'static' layer, and each
        individual source or sample map of a layer colored by membership.
        """
        def colorer(entry: kgml.Entry) -> Tuple[str, float]:
            return color_hexcode, 1.0
        return colorer

    def _entry_sources(
        self,
        entry: kgml.Entry,
        membership: Dict[str, List[str]],
        use_reaction_attribute: bool
    ) -> Set[str]:
        """
        The sources containing a map element, pooled across the accessions it stands for.

        A single line, box or circle can stand for several KOs, reactions or compounds, so the
        sources it is present in are the union of its accessions' sources. This one definition serves
        everything that counts an element's sources: the colorer that colors by that count, the pass
        that finds the highest count on the drawn maps, and the within-group counts of grouped
        individual maps.

        Parameters
        ==========
        entry : anvio.kgml.Entry
            The map element.

        membership : Dict[str, List[str]]
            Maps each accession to the sources containing it.

        use_reaction_attribute : bool
            Read the element's KEGG reaction IDs rather than its KO IDs.

        Returns
        =======
        Set[str]
            The sources containing the element, empty if it is in none of them.
        """
        sources: Set[str] = set()
        for accession in self._get_entry_kegg_ids(entry, use_reaction_attribute):
            if accession in membership:
                sources.update(membership[accession])
        return sources

    def _membership_categorizer(
        self,
        membership: Dict[str, List[str]],
        group_sources: Union[Dict[str, List[str]], None],
        group_threshold: Union[float, None],
        use_reaction_attribute: bool
    ):
        """
        Build a function giving the categories a map element counts as being in, or None for none.

        Ungrouped, an element's categories are the sources containing it ('_entry_sources'). With
        'group_sources' they are the qualifying groups instead: those where the proportion of the
        group's own sources containing the element meets 'group_threshold'. Both the colorer that
        colors an element by how many categories it is in ('_membership_colorer') and the pass that
        finds the highest such count across the drawn maps ('_find_observed_counts') go through
        here, so the count a color stands for and the count the scale was built for cannot disagree.
        """
        grouped = group_sources is not None
        group_source_count: Dict[str, int] = {}
        source_group: Dict[str, str] = {}
        if grouped:
            for group, sources in group_sources.items():
                group_source_count[group] = len(sources)
                for source in sources:
                    source_group[source] = group

        def categorize(entry: kgml.Entry) -> Union[Set[str], None]:
            sources = self._entry_sources(entry, membership, use_reaction_attribute)
            if not sources:
                return None
            if not grouped:
                return sources
            group_counts = {group: 0 for group in group_source_count}
            for source in sources:
                if source in source_group:
                    group_counts[source_group[source]] += 1
            categories = set()
            for group, count in group_counts.items():
                proportion = count / group_source_count[group]
                if (proportion > 0) if group_threshold == 0 else (proportion >= group_threshold):
                    categories.add(group)
            return categories or None
        return categorize

    def _membership_colorer(
        self,
        membership: Dict[str, List[str]],
        color_priorities: List[Tuple[str, float]],
        category_combos: Union[List[Tuple[str]], None],
        group_sources: Union[Dict[str, List[str]], None],
        group_threshold: Union[float, None],
        use_reaction_attribute: bool
    ):
        """
        Build a colorer that colors an Entry by the sources (or groups) containing its accessions.

        See '_draw_map_elements' for the colorer contract. An Entry's containing sources are pooled
        across its accessions via 'membership'; with 'group_sources', the qualifying groups (those
        meeting 'group_threshold' among their sources) are used instead. The color and its priority
        are looked up in 'color_priorities' by the count of categories ('category_combos' None) or
        by the position of their exact combination (by membership), so a color shared by more than
        one count — which a continuous count scale allows — still carries that count's own priority.
        An Entry in no source, or, when grouped, in no qualifying group, is left uncolored.
        """
        combo_lookup: Dict[Tuple[str], Tuple[str]] = {}
        if category_combos is not None:
            for combo in category_combos:
                combo_lookup[tuple(sorted(combo))] = combo
        categorize = self._membership_categorizer(
            membership, group_sources, group_threshold, use_reaction_attribute
        )
        count_scale_top = len(color_priorities)

        def colorer(entry: kgml.Entry) -> Union[Tuple[str, float], None]:
            categories = categorize(entry)
            if categories is None:
                return None
            if category_combos is None:
                # A count scale can stop below the highest count there is, since '--count-scale-max'
                # takes a ceiling, so a count above the top of the scale takes the top color, as a
                # value above the top of a value scale does.
                return color_priorities[min(len(categories), count_scale_top) - 1]
            return color_priorities[
                category_combos.index(combo_lookup[tuple(sorted(categories))])
            ]
        return colorer

    def _draw_map_elements(
        self,
        pathway_number: str,
        layer_specs: List[dict],
        output_dir: str,
        draw_map_lacking_data: bool = False
    ) -> bool:
        """
        Draw one pathway map, coloring each layer's elements via that layer's colorer.

        This is the single mode-agnostic draw path of the element engine. Each layer contributes a
        'colorer' mapping a matching Entry to a '(color_hexcode, priority)' pair, or None to leave it
        uncolored, computed however that layer's mode requires (continuous value, sample/group
        membership, or a single fixed color). All layers are staged onto one 'Pathway' and applied in
        a single 'set_color_priority' call (list order sets draw order: earlier layers render beneath
        later ones), so one map can mix a quantitative layer with a presence/categorical layer. The
        map is drawn once.

        Parameters
        ==========
        pathway_number : str
            Numeric ID of the map to draw.

        layer_specs : List[dict]
            Per-layer specs, ordered so earlier layers render beneath later ones. Each has:
            'element_type' ('reaction'/'compound'), 'use_reaction_attribute' (bool), 'entry_keys'
            (the accessions this layer touches, for entry lookup and the compound warning), 'colorer'
            (callable mapping an Entry to '(color_hexcode, priority)' or None), and 'derived_compound'
            (how a reaction layer derives compound colors when no compound layer is present on a
            global/overview map: a '(mode, colormap)' pair such as ('average', cmap) or ('high',
            None); None on a compound layer, where it is never read).

        output_dir : str
            Path to the output directory in which the map PDF is drawn.

        draw_map_lacking_data : bool, False
            If False, only draw the map if some layer matched any of its accessions.

        Returns
        =======
        bool
            True if the map was drawn, False if it was skipped for lacking data.
        """
        pathway = self._get_pathway(pathway_number)

        color_priority: dict = {}
        found_entries = False
        for spec in layer_specs:
            entries = self._find_element_entries(
                pathway, spec['use_reaction_attribute'], spec['entry_keys']
            )
            if entries:
                found_entries = True
            for entry in entries:
                colored = spec['colorer'](entry)
                if colored is None:
                    continue
                color_hexcode, priority = colored
                self._stage_element_color(
                    pathway, entry, spec['element_type'], color_hexcode, priority, color_priority
                )

        if not found_entries and not draw_map_lacking_data:
            return False

        compound_layer = next(
            (spec for spec in layer_specs if spec['element_type'] == 'compound'), None
        )
        if compound_layer is not None:
            self._warn_unrenderable_compounds(pathway, compound_layer['entry_keys'])

        # When a compound layer supplies compound colors, or the map is neither global nor overview,
        # colors are applied directly. Otherwise, on a reaction-only global/overview map, compounds
        # are derived from their associated reactions per the reaction layer's 'derived_compound'.
        recolor = 'g' if pathway.is_global_map else 'w'
        if compound_layer is not None or not (pathway.is_global_map or pathway.is_overview_map):
            pathway.set_color_priority(color_priority, recolor_unprioritized_entries=recolor)
        else:
            reaction_layer = next(
                spec for spec in layer_specs if spec['element_type'] == 'reaction'
            )
            mode, colormap = reaction_layer['derived_compound']
            if colormap is None:
                pathway.set_color_priority(
                    color_priority,
                    recolor_unprioritized_entries=recolor,
                    color_associated_compounds=mode
                )
            else:
                pathway.set_color_priority(
                    color_priority,
                    recolor_unprioritized_entries=recolor,
                    color_associated_compounds=mode,
                    colormap=colormap
                )

        self._draw_map(pathway, output_dir)
        return True

    def _draw_map_element_presence(
        self,
        pathway_number: str,
        layers: List[dict],
        output_dir: str,
        draw_map_lacking_data: bool = False
    ) -> bool:
        """
        Draw one pathway map, coloring each layer's present elements a single fixed color.

        Every present element of a layer is colored that layer's 'color_hexcode', with no value or
        source driving the choice. Both layers are staged onto one 'Pathway' and applied in a single
        'set_color_priority' call (earlier layers render beneath later ones), then the map is drawn
        once. This is the single-color path used by '_map_kos_fixed_colors', which serves the
        single-source database and reaction-network-JSON inputs; layers whose color varies by value
        or by source go through '_map_elements' instead.

        Parameters
        ==========
        pathway_number : str
            Numeric ID of the map to draw.

        layers : List[dict]
            Per-layer specs, each with keys 'element_type' ('reaction'/'compound'),
            'use_reaction_attribute' (bool), 'accessions' (a set/dict of the accessions to color),
            and 'color_hexcode'. Ordered so earlier layers render beneath later ones.

        output_dir : str
            Path to the output directory in which the map PDF is drawn.

        draw_map_lacking_data : bool, False
            If False, only draw the map if some layer contains any of its accessions.

        Returns
        =======
        bool
            True if the map was drawn, False if it was skipped for lacking data.
        """
        specs = []
        for layer in layers:
            # On a reaction-only global/overview map, compounds take the color of the highest-
            # priority reaction they touch (as the single-source database path does).
            derived_compound = ('high', None) if layer['element_type'] == 'reaction' else None
            specs.append({
                'element_type': layer['element_type'],
                'use_reaction_attribute': layer['use_reaction_attribute'],
                'entry_keys': layer['accessions'],
                'colorer': self._single_color_colorer(layer['color_hexcode']),
                'derived_compound': derived_compound
            })
        return self._draw_map_elements(
            pathway_number, specs, output_dir, draw_map_lacking_data=draw_map_lacking_data
        )

    @staticmethod
    def _check_contigs_db(contigs_db: str) -> None:
        """
        Check the validity of an expected contigs database.

        Parameters
        ==========
        contigs_db : str
            File path to an expected contigs database.
        """
        if not os.path.exists(contigs_db):
            raise ConfigError(
                f"There was no file at the following expected contigs database path: '{contigs_db}'"
            )

        contigs_db_info = dbinfo.ContigsDBInfo(contigs_db, dont_raise=True, expecting='contigs')
        if contigs_db_info is None:
            raise ConfigError(
                "The file at the following expected contigs database path is not a contigs "
                f"database: '{contigs_db}'"
            )

    @staticmethod
    def _check_contigs_db_ko_annotation(contigs_db: str) -> None:
        """
        Check that a contigs database was annotated with KOs.

        Parameters
        ==========
        contigs_db : str
            File path to a contigs database.
        """
        contigs_db_info = dbinfo.ContigsDBInfo(contigs_db, expecting='contigs')
        if 'KOfam' not in contigs_db_info.get_functional_annotation_sources():
            raise ConfigError(
                f"The contigs database, '{contigs_db}', was never annotated with KOs. This can be "
                "rectified by running `anvi-run-kegg-kofams` on the database."
            )

    @staticmethod
    def _check_genomes_storage_db(genomes_storage_db: str) -> None:
        """
        Check the validity of an expected genomes storage database.

        Parameters
        ==========
        genomes_storage_db : str
            File path to an expected genomes storage database.
        """
        if not os.path.exists(genomes_storage_db):
            raise ConfigError(
                "There was no file at the following expected genomes storage database path: "
                f"'{genomes_storage_db}'"
            )

        gsdb_info = dbinfo.GenomeStorageDBInfo(
            genomes_storage_db, dont_raise=True, expecting='genomestorage'
        )
        if gsdb_info is None:
            raise ConfigError(
                "The file at the following expected genomes storage database path is not a genomes "
                f"storage database: '{genomes_storage_db}'"
            )

    @staticmethod
    def _check_genomes_storage_ko_annotation(genomes_storage_db: str) -> None:
        """
        Check that a genomes storage database was annotated with KOs.

        Parameters
        ==========
        genomes_storage_db : str
            File path to a genomes storage database.
        """
        gsdb_info = dbinfo.GenomeStorageDBInfo(genomes_storage_db, expecting='genomestorage')
        if 'KOfam' not in gsdb_info.get_functional_annotation_sources():
            raise ConfigError(
                f"The genomes storage database, '{genomes_storage_db}', was never annotated with "
                "KOs. The genomes storage should be remade with annotated genomes, which can be "
                "rectified by running `anvi-run-kegg-kofams` on the genome databases."
            )

    @staticmethod
    def _check_contigs_dbs(contigs_dbs: Iterable[str]) -> None:
        """
        Check the validity of expected contigs databases.

        Parameters
        ==========
        contigs_dbs : Iterable[str]
            File paths to expected contigs databases.
        """
        invalid_paths: List[str] = []
        invalid_filetypes: List[str] = []
        for contigs_db in contigs_dbs:
            if not os.path.exists(contigs_db):
                invalid_paths.append(contigs_db)
            if invalid_paths:
                continue

            contigs_db_info = dbinfo.ContigsDBInfo(contigs_db, dont_raise=True, expecting='contigs')
            if contigs_db_info is None:
                invalid_filetypes.append(contigs_db)
            if invalid_filetypes:
                continue

        if invalid_paths:
            paths = ', '.join([f'{path}' for path in invalid_paths])
            raise ConfigError(
                f"There were no files at the following expected contigs database paths: {paths}"
            )

        if invalid_filetypes:
            paths = ', '.join([f'{path}' for path in invalid_filetypes])
            raise ConfigError(
                "The files at the following expected contigs database paths are not contigs "
                f"databases: {paths}"
            )

    @staticmethod
    def _check_contigs_dbs_ko_annotation(contigs_dbs: Iterable[str]) -> None:
        unannotated: List[str] = []
        for contigs_db in contigs_dbs:
            contigs_db_info = dbinfo.ContigsDBInfo(contigs_db, expecting='contigs')
            if 'KOfam' not in contigs_db_info.get_functional_annotation_sources():
                unannotated.append(contigs_db)
            if unannotated:
                continue

        if unannotated:
            paths = ', '.join([f'{path}' for path in unannotated])
            raise ConfigError(
                "The following contigs databases were never annotated with KOs, but this can be "
                f"rectified by running `anvi-run-kegg-kofams` on them: {paths}"
            )

    @staticmethod
    def _check_pan_db(pan_db: str) -> None:
        """
        Check the validity of an expected pan database.

        Parameters
        ==========
        pan_db : str
            File path to an expected pan database.
        """
        if not os.path.exists(pan_db):
            raise ConfigError(
                f"There was no file at the following expected pan database path: '{pan_db}'"
            )

        pan_db_info = dbinfo.PanDBInfo(pan_db, dont_raise=True, expecting='pan')
        if pan_db_info is None:
            raise ConfigError(
                "The file at the following expected pan database path is not a pan database: "
                f"'{pan_db}'"
            )

    def _find_maps(self, output_dir: str, patterns: List[str] = None) -> List[str]:
        """
        Find the numeric IDs of maps to draw, checking that the map can be drawn in the target
        output directory.

        Parameters
        ==========
        output_dir : str
            Path to the output directory in which pathway map PDF files are drawn.

        patterns : List[str], None
            Regex patterns of pathway numbers, which are five digits.
        """
        if patterns is None:
            pathway_numbers = self.available_pathway_numbers
        else:
            pathway_numbers = self._get_pathway_numbers_from_patterns(patterns)

        if self.pathway_categorization is not None:
            missing_pathway_numbers: list[str] = []
            for pathway_number in pathway_numbers:
                if pathway_number not in self.pathway_categorization:
                    missing_pathway_numbers.append(pathway_number)
            if missing_pathway_numbers:
                message = ', '.join(f"'{p}'" for p in missing_pathway_numbers)
                raise AssertionError(
                    "The KEGG BRITE hierarchy of pathway maps, 'br08901', did not contain all of "
                    "the pathway numbers requested to be drawn. This prevents output files from "
                    "being categorized in a subdirectory structure corresponding to the hierarchy. "
                    "The option to categorize files cannot be used. It would be worthwhile to make "
                    "the developers aware of this error so they can hopefully figure out a "
                    "solution. Here is the list of pathway numbers missing from the hierarchy: "
                    f"{message}"
                )

        if not self.overwrite_output:
            for pathway_number in pathway_numbers:
                out_basename = self._map_basename(pathway_number)
                if self.pathway_categorization is None:
                    out_path = os.path.join(output_dir, out_basename)
                else:
                    out_path = os.path.join(
                        output_dir, *self.pathway_categorization[pathway_number], out_basename
                    )
                if os.path.exists(out_path):
                    raise ConfigError(
                        f"Output files would be overwritten in the output directory, {output_dir}. "
                        "Either delete the contents of the directory, or use the option to "
                        "overwrite output destinations."
                    )

        return pathway_numbers

    def _get_pathway_numbers_from_patterns(self, patterns: Iterable[str]) -> List[str]:
        """
        Among pathways available in the KEGG data directory, get those with ID numbers matching the
        given regex patterns.

        Parameters
        ==========
        patterns : Iterable[str]
            Regex patterns of pathway numbers, which are five digits.

        Returns
        =======
        List[str]
            Pathway numbers matching the regex patterns.
        """
        pathway_numbers: List[str] = []
        for pattern in patterns:
            for available_pathway_number in self.available_pathway_numbers:
                if re.match(pattern, available_pathway_number):
                    pathway_numbers.append(available_pathway_number)

        # Maintain the order of pathway numbers recovered from patterns.
        seen = set()
        return [
            pathway_number for pathway_number in pathway_numbers
            if not (pathway_number in seen or seen.add(pathway_number))
        ]

    def _draw_map_kos_original_color(
        self,
        pathway_number: str,
        ko_ids: Iterable[str],
        output_dir: str,
        draw_map_lacking_data: bool = False
    ) -> bool:
        """
        Draw a pathway map, highlighting reactions containing select KOs in the color or colors
        originally used in the reference map.

        Parameters
        ==========
        pathway_number : str, None
            Numeric ID of the map to draw.

        ko_ids : Iterable[str]
            Select KOs, any of which in the map are colored.

        output_dir : str
            Path to the output directory in which map PDF files are drawn, created if it doesn't
            already exist.

        draw_map_lacking_data : bool, False
            If False, by default, only draw the map if it contains any of the select KOs. If True,
            draw the map regardless, meaning that nothing may be highlighted.

        Returns
        =======
        bool
            True if the map was drawn, False if the map was not drawn because it did not contain any
            of the select KOs and 'draw_map_lacking_data' was False.
        """
        pathway = self._get_pathway(pathway_number)

        select_entries = pathway.get_entries(kegg_ids=ko_ids)
        if not select_entries and not draw_map_lacking_data:
            return False

        # Set "secondary" colors of ortholog Graphics elements for reactions containing select KOs:
        # white background color of lines or black foreground text of rectangles. For other Graphics
        # elements, change the 'fgcolor' attribute to a nonsense value to ensure that the elements
        # with prioritized colors can be distinguished from other elements. Also, in overview and
        # standard maps, widen lines from the base map default of 1.0.
        all_entries = pathway.get_entries(entry_type='ortholog')
        select_uuids = [entry.uuid for entry in select_entries]
        prioritized_colors: Dict[str, List[Tuple[str, str]]] = {}
        for entry in all_entries:
            if entry.uuid in select_uuids:
                for uuid in entry.children['graphics']:
                    graphics: kgml.Graphics = pathway.uuid_element_lookup[uuid]
                    if pathway.is_global_map:
                        assert graphics.type == 'line'
                        graphics.bgcolor = '#FFFFFF'
                    elif pathway.is_overview_map:
                        assert graphics.type == 'line'
                        graphics.bgcolor = '#FFFFFF'
                        graphics.width = 5.0
                    else:
                        if graphics.type == 'rectangle':
                            graphics.fgcolor = '#000000'
                        elif graphics.type == 'line':
                            graphics.bgcolor = '#FFFFFF'
                            graphics.width = 5.0
                        else:
                            self.progress.end()
                            raise ConfigError(
                                f"Reaction elements are expected to be drawn as a rectangle or a "
                                f"line, but an ortholog entry of KEGG pathway map "
                                f"{pathway.number} has a graphics element of type "
                                f"'{graphics.type}', which anvi'o cannot color."
                            )
                    try:
                        graphics_type_prioritized_colors = prioritized_colors[graphics.type]
                    except:
                        prioritized_colors[graphics.type] = graphics_type_prioritized_colors = []
                    graphics_type_prioritized_colors.append((graphics.fgcolor, graphics.bgcolor))
            else:
                for uuid in entry.children['graphics']:
                    graphics: kgml.Graphics = pathway.uuid_element_lookup[uuid]
                    graphics.fgcolor = '0'

        # By default, global maps but not overview and standard maps display reaction graphics in
        # more than one color. Give higher priority to reaction entries that are encountered later
        # (occur further down in the KGML file), and would thus be rendered above earlier reactions.
        color_priority: Dict[str, Dict[str, Dict[Tuple[str, str], float]]] = {'ortholog': {}}
        for graphics_type, graphics_type_prioritized_colors in prioritized_colors.items():
            seen = set()
            unique_prioritized_colors = [
                colors for colors in graphics_type_prioritized_colors
                if not (colors in seen or seen.add(colors))
            ]
            priorities = np.linspace(0, 1, len(unique_prioritized_colors) + 1)[1: ]
            graphics_type_color_priority = {
                colors: priority for colors, priority in zip(unique_prioritized_colors, priorities)
            }
            color_priority['ortholog'][graphics_type] = graphics_type_color_priority

        # Recolor "unprioritized" reactions to a background color. In global and overview maps,
        # recolor circles to reflect the colors of prioritized reactions involving the compounds.
        if pathway.is_global_map:
            recolor_unprioritized_entries = 'g'
            color_associated_compounds = 'high'
        elif pathway.is_overview_map:
            recolor_unprioritized_entries = 'w'
            color_associated_compounds = 'high'
        else:
            recolor_unprioritized_entries = 'w'
            color_associated_compounds = None
        pathway.set_color_priority(
            color_priority,
            recolor_unprioritized_entries=recolor_unprioritized_entries,
            color_associated_compounds=color_associated_compounds
        )

        self._draw_map(pathway, output_dir)

        return True

    def _draw_map(self, pathway: kgml.Pathway, output_dir: str) -> None:
        """
        Draw a map given the KGML pathway data.

        Parameters
        ==========
        pathway : anvio.kgml.Pathway
            KGML pathway element object.

        output_dir : str
            Path to the output directory in which the pathway map PDF file is drawn. The directory
            is created if it does not exist.
        """
        out_basename = self._map_basename(pathway.number)

        if self.pathway_categorization is None:
            out_dir = output_dir
            out_path = os.path.join(output_dir, out_basename)
        else:
            out_dir = os.path.join(output_dir, *self.pathway_categorization[pathway.number])
            out_path = os.path.join(out_dir, out_basename)
        os.makedirs(out_dir, exist_ok=True)

        self.drawer.draw_map(pathway, out_path)

        if self.pathway_categorization is None:
            return

        self._link_map_flat(output_dir, out_path)

    def _link_map_flat(self, output_dir: str, map_path: str) -> None:
        """
        Link to a map file from the flat directory, which gathers all maps in a single directory.

        Where '--categorize-files' nests map files in BRITE subdirectories, the flat directory
        contains all of the files for easier browsing.

        Parameters
        ==========
        output_dir : str
            Path to the output directory the map file was drawn in, whether directly or in a
            subdirectory of it. The flat directory is made there if it does not already exist.

        map_path : str
            Path of the map file to link to.
        """
        self._link_map(map_path, os.path.join(output_dir, FLAT_SUBDIR, os.path.basename(map_path)))

    def _collate_maps_by_map(self, output_dir: str, categories: List[str]) -> int:
        """
        Gather the maps drawn for individual categories into a directory per map.

        Each map gets a directory holding one file per category, named after the category, so that a
        single map can be compared across categories by stepping through one directory — the way a
        file browser's preview moves from one file to the next — rather than by opening a file in
        each category's own directory.

        Every file gathered is a link to a map already drawn under 'individual', so this second
        arrangement costs no disk space of its own.

        Which maps to gather is read off the drawn files themselves rather than rebuilt from pathway
        numbers, so it can neither disagree with what was drawn nor lose track of where
        '--categorize-files' put it: a map's place within a category's directory becomes its place
        here, BRITE subdirectories and all.

        Parameters
        ==========
        output_dir : str
            Path to the output directory holding the 'individual' subdirectory.

        categories : List[str]
            Names of the categories whose maps are gathered, each a subdirectory of 'individual'.

        Returns
        =======
        int
            The number of maps gathered, which is the number of directories written.
        """
        collated_dir = os.path.join(output_dir, COLLATED_SUBDIR)
        map_dirs: Set[str] = set()
        for category in categories:
            self.progress.update(category)
            category_dir = os.path.join(output_dir, INDIVIDUAL_SUBDIR, category)
            for dir_path, subdir_names, basenames in os.walk(category_dir):
                # Everything a category's directory holds is one of its maps, save for the colorbar
                # keying them and the flat directory of links to those very same maps.
                if FLAT_SUBDIR in subdir_names:
                    subdir_names.remove(FLAT_SUBDIR)
                relative_dir = os.path.relpath(dir_path, category_dir)
                for basename in basenames:
                    map_name, extension = os.path.splitext(basename)
                    if extension != '.pdf' or basename == CATEGORY_COLORBAR_BASENAME:
                        continue
                    map_dir = os.path.normpath(os.path.join(collated_dir, relative_dir, map_name))
                    map_dirs.add(map_dir)
                    self._link_map(
                        os.path.join(dir_path, basename),
                        os.path.join(map_dir, f'{category}{extension}')
                    )
        return len(map_dirs)

    def _link_map(self, map_path: str, link_path: str) -> None:
        """
        Link from elsewhere in the output to a map file that has already been drawn.

        A hard link is made, which a file browser cannot tell from the map file itself and which
        survives the output directory being moved or renamed. Where the filesystem refuses a hard
        link, a relative symbolic link stands in: a file browser follows it just as readily, and it
        holds for as long as the two files keep their places within the output directory.

        Parameters
        ==========
        map_path : str
            Path to the map file that was drawn.

        link_path : str
            Path of the link to make. Its directory is created if it does not already exist, and
            whatever the path may already hold is replaced.
        """
        link_dir = os.path.dirname(link_path)
        os.makedirs(link_dir, exist_ok=True)
        # 'lexists' rather than 'exists', so that a link an earlier run left broken is replaced
        # rather than left in place for the call below to trip over.
        if os.path.lexists(link_path):
            os.remove(link_path)
        try:
            os.link(map_path, link_path)
        except OSError:
            os.symlink(os.path.relpath(map_path, link_dir), link_path)

    def _draw_map_grids(
        self,
        pathway_numbers: List[str],
        draw_categories: List[str],
        draw_grid_categories: List[str],
        draw_files_categories: List[str],
        output_dir: str,
        drawn: Dict[Literal['unified', 'individual', 'grid'], Dict],
        group_scale_tops: Dict[str, int] = None,
        group_color_priorities: Dict[str, List[Tuple[str, float]]] = None,
        group_cmaps: Dict[str, mcolors.Colormap] = None,
        group_colormap_scheme: str = None,
        check_maps_lacking_kos: bool = True,
        source_type: str = 'unknown'
    ) -> None:
        """
        Make map grids from arbitrary categories of data sources or groups of data sources, where
        each category corresponds an individual map cell in the grid, e.g., categories of the
        ungrouped 'pangenome' data source type are genomes, categories of the ungrouped 'contigs
        database' type are contigs databases, and categories of grouped 'pangenome' or 'contigs
        database' types are groups.

        This method picks up in map methods for multiple data sources after unified and individual
        map files are drawn.

        Parameters
        ==========
        pathway_numbers : List[str]
            IDs of pathway maps to draw.

        draw_categories : List[str]
            All categories for which map files are attempted to be drawn.

        draw_grid_categories : List[str]
            All categories for which map grids are attempted to be drawn.

        draw_files_categories : List[str]
            All categories for which individual map files were drawn.

        output_dir : str
            Path to the output directory, which should already exist, in which pathway map PDF files
            are drawn.

        drawn : Dict[Literal['unified', 'individual', 'grid'], Dict]
            Record of drawn map files.

        group_scale_tops : Dict[str, int], None
            Used to draw a per-group colorbar of source counts. Keys are group names; values are the
            count each group's scale runs up to ('_resolve_count_scale_top'). Left None when no
            layer colors its individual group maps by within-group source counts, as when a group
            map is colored by value instead, in which case no per-group colorbars are drawn.

        group_color_priorities : Dict[str, List[Tuple[str, float]]], None
            Used together with 'group_scale_tops' to draw the per-group colorbars. Keys are group
            names; values are lists of (color hex code, priority) pairs in ascending order of
            within-group source count. Reactions assigned higher priority colors are drawn over
            reactions assigned lower priority colors. Left None whenever 'group_scale_tops' is.

        group_cmaps : Dict[str, matplotlib.colors.Colormap], None
            The colormap each group's colors were sampled from, which its colorbar spans when
            'group_colormap_scheme' draws that bar as a gradient ('_group_map_colors'). Empty or
            None for a scheme that draws discrete bands, which span no colormap.

        group_colormap_scheme : str, None
            How the per-group colorbars are drawn ('by_count' in discrete bands of one color per
            count, 'by_count_continuous' as a gradient from the lowest count to the highest),
            settled for the whole run by '_group_map_colors'. Left None whenever 'group_scale_tops'
            is.

        check_maps_lacking_kos : bool, True
            If True, check for "empty" individual map files that are needed to complete the map grid
            and draw them as temporary files, deleting them at after map grids are drawn.

        source_type : Literal['contigs database', 'pangenome', 'sample'], 'unknown'
            The kind of source a group is made of, which labels the per-group colorbars of source
            counts. The text engine passes 'sample', including for a grouped run, since a group is
            made of samples.
        """
        self.progress.new("Drawing map grid")
        self.progress.update("...")

        # Draw empty maps needed to fill in grids.
        paths_to_remove: List[str] = []
        if check_maps_lacking_kos:
            # Make a new dictionary with outer keys being pathway numbers, inner dictionaries
            # indicating which maps were drawn per category (e.g., database, pan genome, or group).
            drawn_pathway_number: Dict[str, Dict[str, bool]] = {}
            for category, drawn_category in drawn['individual'].items():
                for pathway_number, drawn_map in drawn_category.items():
                    try:
                        drawn_pathway_number[pathway_number][category] = drawn_map
                    except KeyError:
                        drawn_pathway_number[pathway_number] = {category: drawn_map}

            # Draw empty maps as needed, for pathways with some but not all maps drawn.
            progress = self.progress
            self.progress = terminal.Progress(verbose=False)
            run = self.run
            self.run = terminal.Run(verbose=False)

            for pathway_number, drawn_category in drawn_pathway_number.items():
                if set(drawn_category.values()) != set([True, False]):
                    continue

                pathway_basename = self._map_basename(pathway_number)
                for category, drawn_map in drawn_category.items():
                    if drawn_map:
                        continue

                    out_dir = os.path.join(output_dir, INDIVIDUAL_SUBDIR, category)

                    self._map_kos_fixed_colors(
                        [], out_dir, [pathway_number], draw_maps_lacking_data=True
                    )

                    if self.pathway_categorization is None:
                        out_path = os.path.join(out_dir, pathway_basename)
                    else:
                        out_path = os.path.join(
                            out_dir, *self.pathway_categorization[pathway_number], pathway_basename
                        )
                    paths_to_remove.append(out_path)
                    if self.categorize_files:
                        paths_to_remove.append(os.path.join(out_dir, FLAT_SUBDIR, pathway_basename))

            self.progress = progress
            self.run = run

        # Draw map grids.
        grid_dir = os.path.join(output_dir, GRID_SUBDIR)
        filesnpaths.gen_output_directory(grid_dir, progress=self.progress, run=self.run)

        if group_scale_tops is not None:
            # Draw colorbars for each group.
            if source_type == 'pangenome':
                label = 'genome count'
            elif source_type == 'contigs database':
                label = 'database count'
            elif source_type == 'sample':
                label = 'sample count'
            else:
                label = 'source count'
            for group in draw_categories:
                colorbar_path = os.path.join(grid_dir, f'colorbar_{group}.pdf')
                if group_colormap_scheme == 'by_count_continuous':
                    # A gradient from the lowest count to the highest, as the same group's own map
                    # directory gets, so a grid and the maps it is made of are keyed the same way.
                    self._draw_quantitative_colorbar(
                        group_cmaps[group], 1, group_scale_tops[group], colorbar_path, label,
                        integer_ticks=True
                    )
                else:
                    self.colorbar_drawer.draw_discrete(
                        [color for color, _ in group_color_priorities[group]],
                        colorbar_path,
                        color_labels=range(1, group_scale_tops[group] + 1),
                        label=label
                    )

        for pathway_number in pathway_numbers:
            self.progress.update(pathway_number)
            pathway_basename = self._map_basename(pathway_number)
            unified_dir = os.path.join(output_dir, UNIFIED_SUBDIR)
            if self.pathway_categorization is None:
                unified_map_path = os.path.join(unified_dir, pathway_basename)
            else:
                out_dir = os.path.join(unified_dir, *self.pathway_categorization[pathway_number])
                unified_map_path = os.path.join(out_dir, pathway_basename)
            if not os.path.exists(unified_map_path):
                continue
            in_paths = [unified_map_path]
            if source_type == 'pangenome':
                labels = ['pangenome']
            else:
                labels = ['all']

            for category in draw_grid_categories:
                if self.pathway_categorization is None:
                    individual_map_path = os.path.join(
                        output_dir, INDIVIDUAL_SUBDIR, category, pathway_basename
                    )
                else:
                    out_dir = os.path.join(
                        output_dir, INDIVIDUAL_SUBDIR, category,
                        *self.pathway_categorization[pathway_number]
                    )
                    individual_map_path = os.path.join(out_dir, pathway_basename)
                if not os.path.exists(individual_map_path):
                    break
                in_paths.append(individual_map_path)
                labels.append(category)
            else:
                if self.pathway_categorization is None:
                    out_path = os.path.join(grid_dir, pathway_basename)
                else:
                    out_dir = os.path.join(grid_dir, *self.pathway_categorization[pathway_number])
                    os.makedirs(out_dir, exist_ok=True)
                    out_path = os.path.join(out_dir, pathway_basename)
                self.grid_drawer.draw(in_paths, out_path, labels=labels)
                if self.pathway_categorization is not None:
                    self._link_map_flat(grid_dir, out_path)
                drawn['grid'][pathway_number] = True

        self.progress.end()

        # Remove individual maps and their links that were only needed for map grids.
        for path in paths_to_remove:
            os.remove(path)
        for category in set(draw_categories).difference(set(draw_files_categories)):
            shutil.rmtree(os.path.join(output_dir, INDIVIDUAL_SUBDIR, category))
            drawn['individual'].pop(category)

        # The directory holding the individual maps exists only for them, so it should not outlive
        # them: a run that asked for grids alone drew those maps as grid panels and has just removed
        # them again, which would otherwise leave an empty directory behind.
        individual_dir = os.path.join(output_dir, INDIVIDUAL_SUBDIR)
        if os.path.isdir(individual_dir) and not os.listdir(individual_dir):
            os.rmdir(individual_dir)

        # If map files were categorized in a subdirectory structure, remove subdirectories that no
        # longer contain files.
        if not self.categorize_files:
            return
        for category in draw_files_categories:
            category_dir = os.path.join(output_dir, INDIVIDUAL_SUBDIR, category)
            for dir_path, subdir_names, filenames in os.walk(category_dir, topdown=False):
                if dir_path != category_dir and not os.listdir(dir_path):
                    os.rmdir(dir_path)

    def _get_pathway(self, pathway_number: str) -> kgml.Pathway:
        """
        Get a Pathway object for the KGML file used in drawing a pathway map.

        Parameters
        ==========
        pathway_number : str
            Numeric ID of the map to draw.

        Returns
        =======
        kgml.Pathway
            Representation of the KGML file as an object.
        """
        # KOs correspond to arrows rather than boxes in global and overview maps.
        is_global_map = False
        is_overview_map = False
        if re.match(GLOBAL_MAP_ID_PATTERN, pathway_number):
            is_global_map = True
        elif re.match(OVERVIEW_MAP_ID_PATTERN, pathway_number):
            is_overview_map = True

        # A 1x resolution global 'KO' image is used as the base of the drawing, whereas a 2x
        # overview or standard 'map' image is used as the base. The global 'KO' image grays out
        # all reaction arrows that are not annotated by KO ID. Select the KGML file accordingly.
        if is_global_map:
            kgml_path = os.path.join(
                self.kegg_context.kgml_1x_ko_dir, f'ko{pathway_number}.xml'
            )
        else:
            kgml_path = os.path.join(
                self.kegg_context.kgml_2x_ko_dir, f'ko{pathway_number}.xml'
            )
        pathway = self.xml_ops.load(kgml_path)
        if self.ignore_compound_rectangles:
            self._zero_out_compound_rectangles(pathway)

        return pathway

    def _map_basename(self, pathway_number: str) -> str:
        """
        Name the output file of a pathway map.

        The name is the map's KEGG accession — its number under the 'ko' prefix of the reference
        KGML files the maps are drawn from — and, with 'name_files', the pathway name as well.

        Every context that draws a map or goes looking for one names it here, so the same map is
        named the same under 'unified', 'individual', and 'grid'. That is how a grid finds the
        panels it is made of.

        Parameters
        ==========
        pathway_number : str
            ID number of the pathway, which is five digits.

        Returns
        =======
        str
            The file name, e.g. 'ko00010.pdf', or 'ko00010_Glycolysis_Gluconeogenesis.pdf' where
            file names carry pathway names.
        """
        pathway_name = f'_{self._name_pathway(pathway_number)}' if self.name_files else ''
        return f'ko{pathway_number}{pathway_name}.pdf'

    def _name_pathway(self, pathway_number: str) -> str:
        """
        Format the pathway name corresponding to the number for suitability in file paths.

        Replace all non-alphanumeric characters except parentheses, brackets, and curly braces with
        underscores. Replace multiple consecutive underscores with a single underscore. Strip
        leading and trailing underscores.

        Parameters
        ==========
        pathway_number : str
            Numeric ID of a pathway map.

        Returns
        =======
        str
            Altered version of the pathway name.
        """
        try:
            pathway_name = self.pathway_names[pathway_number]
        except KeyError:
            raise ConfigError(
                f"The pathway number, '{pathway_number}', is not recognized in the table of KEGG "
                "pathway names set up in the KEGG data directory, which can be found here: "
                f"'{self.kegg_context.kegg_pathway_list_file}'."
            )

        altered = re.sub(r'[^a-zA-Z0-9()\[\]\{\}]', '_', pathway_name)
        altered = re.sub(r'_+', '_', altered)
        altered = altered.strip('_')

        return altered

    def _categorize_pathways(self) -> dict[str, list[str]]:
        """
        Categorize pathways in the BRITE hierarchy, 'br08901'.

        Alter category names to make suitable for directory paths. Replace all non-alphanumeric
        characters except parentheses, brackets, and curly braces with underscores. Replace multiple
        consecutive underscores with a single underscore. Strip leading and trailing underscores.

        Returns
        =======
        dict[str, list[str]]
            Keys are pathway numbers. Values are lists of the categories from general to specific.
            For example, '00010': ['Metabolism', 'Carbohydrate metabolism']
        """
        with open(self.kegg_context.kegg_brite_pathways_file) as f:
            hierarchy = json.load(f)
        pathway_categorizations: Dict[str, list[list[str]]] = (
            self.kegg_context.invert_brite_json_dict(hierarchy)
        )

        assert len(set(pathway_categorizations)) == len(pathway_categorizations)

        if sum(set(
            [len(categorizations) for categorizations in pathway_categorizations.values()]
        )) != 1:
            raise AssertionError(
                "The KEGG BRITE hierarchy of pathway maps, 'br08901', did not meet the expectation "
                "that each pathway be categorized in exactly one place in the hierarchy. This "
                "prevents output files from being categorized in a subdirectory structure "
                "corresponding to the hierarchy. The option to categorize files cannot be used. It "
                "would be worthwhile to make the developers aware of this error so they can "
                "hopefully figure out a solution."
            )

        pathway_categorization: Dict[str, list[str]] = {}
        for pathway, categorizations in pathway_categorizations.items():
            pathway_number = pathway[:5]
            assert pathway_number.isdigit()
            assert pathway_number not in pathway_categorization

            categorization = categorizations[0]
            assert categorization[0] == 'br08901'

            altered_categorization: list[str] = []
            for category in categorization[1:]:
                altered = re.sub(r'[^a-zA-Z0-9()\[\]\{\}]', '_', category)
                altered = re.sub(r'_+', '_', altered)
                altered = altered.strip('_')
                altered_categorization.append(altered)

            pathway_categorization[pathway_number] = altered_categorization

        return pathway_categorization

    @staticmethod
    def _zero_out_compound_rectangles(pathway: kgml.Pathway) -> int:
        """
        Zero out the size of KGML compound Entry rectangle Graphics in the Pathway.

        These are found in a small number of KGML files (see 00121, 00621, 01052, 01054), and when
        rendered by 'anvio.kgml' via 'Bio.Graphics.KGML_vis.KGMLCanvas' have the effect of obscuring
        underlying drawings of compound structures in the base map image.

        Parameters
        ==========
        pathway : anvio.kgml.Pathway
            KGML pathway element object.

        Returns
        =======
        int
            Count of zeroed out Graphics.
        """
        rectangle_count = 0
        for compound_entry in pathway.get_entries(entry_type='compound'):
            compound_rectangle_uuids: List[str] = []
            for uuid in compound_entry.children['graphics']:
                graphics: kgml.Graphics = pathway.uuid_element_lookup[uuid]
                if graphics.type == 'rectangle':
                    graphics.width = 0.0
                    graphics.height = 0.0
                    rectangle_count += 1
        return rectangle_count

class ColorbarDrawer:
    """
    Writes standalone colorbar image files.

    Parameters
    ==========
    overwrite_output : bool
        If True, methods in this class overwrite existing output files.

    figsize : Tuple[int, int]
        Dimensions of the figure in inches.

    orientation : Literal['horizontal', 'vertical']
        Orientation of the colobar.

    tick_fontsize : Union[int, None]
        If None, tick labels are sized to fit the bar: to one color segment on a discrete bar
        ('draw_discrete'), and to 'max_tick_fontsize' on a continuous one, whose few labels always
        have the room. Otherwise set to the provided value.

    tick_fontsize_segment_fraction : float
        How much of one color segment a dynamically sized tick label is allowed to fill, since
        labels one segment apart need a gap between them ('draw_discrete').

    max_tick_fontsize : int
        The largest a dynamically sized tick label gets, however much room a color segment has.

    label_rotation : Union[int, None]
        If None, the colorbar label is rotated 270° if 'orientation' is vertical or 0° if
        horizontal. Otherwise rotate the label by the provided value.

    label_fontsize : int
        Font size of colorbar label.

    labelpad : int
        Spacing of colorbar label from tick labels in points.

    max_integer_ticks : int
        The most ticks a continuous colorbar is given by 'draw_continuous' with 'integer_ticks'.
    """
    def __init__(self, overwrite_output: bool = FORCE_OVERWRITE) -> None:
        """
        Parameters
        ==========
        overwrite_output : bool, FORCE_OVERWRITE
            If True, methods in this class overwrite existing output files.
        """
        self.overwrite_output = overwrite_output

        self.figsize: Tuple[int, int] = (1, 6)
        self.orientation: Literal['horizontal', 'vertical'] = 'vertical'
        self.tick_fontsize: Union[int, None] = None
        self.tick_fontsize_segment_fraction: float = 0.8
        self.max_tick_fontsize: int = 24
        self.label_rotation: int = None
        self.label_fontsize: int = 24
        self.labelpad: int = 30
        self.max_integer_ticks: int = 6

    def draw_discrete(
        self,
        colors: Iterable,
        out_path: str,
        color_labels: Iterable[str] = None,
        label: str = None
    ) -> None:
        """
        Save a standalone discrete (segmented) colorbar to a file, with one color band per color.

        Parameters
        ==========
        colors : Iterable
            Sequence of Matplotlib color specifications for 'matplotlib.colors.ListedColormap' color
            parameter.

        out_path : str
            Path to PDF output file.

        color_labels : Iterable[str], None
            Color segment labels.

        label : str, None
            Overall colorbar label.
        """
        if color_labels is not None:
            assert len(colors) == len(color_labels)

        fig = Figure(figsize=self.figsize)
        ax = fig.subplots()

        cmap = mcolors.ListedColormap(colors)
        norm = mcolors.BoundaryNorm(boundaries=range(len(colors) + 1), ncolors=len(colors))

        cb = fig.colorbar(
            ScalarMappable(norm=norm, cmap=cmap),
            cax=ax,
            orientation=self.orientation
        )

        # Don't show tick marks.
        cb.ax.tick_params(size=0)

        if color_labels:
            if self.tick_fontsize is None:
                # Fit the tick labels to one color segment, which is 1 in the data coordinates of a
                # colorbar built on integer boundaries. Labels sit one segment apart, so type as
                # tall as a segment would leave neighbors touching. 'transData' lands in display
                # coordinates, which are pixels at the figure's dpi, while a font size is in points,
                # so the segment is converted rather than handed over as it stands: without that the
                # type comes out dpi/72 too large, and changes size with a Matplotlib setting that
                # has nothing to do with the geometry of the PDF.
                origin_in_pixels = ax.transData.transform((0, 0))
                if self.orientation == 'vertical':
                    segment_in_pixels = (ax.transData.transform((0, 1)) - origin_in_pixels)[1]
                elif self.orientation == 'horizontal':
                    segment_in_pixels = (ax.transData.transform((1, 0)) - origin_in_pixels)[0]
                else:
                    raise AssertionError
                segment_in_points = segment_in_pixels * 72 / fig.dpi
                tick_fontsize = min(
                    segment_in_points * self.tick_fontsize_segment_fraction,
                    self.max_tick_fontsize
                )
            else:
                tick_fontsize = self.tick_fontsize

            cb.set_ticks(np.arange(len(colors)) + 0.5)
            cb.set_ticklabels(color_labels, fontsize=tick_fontsize)

        if label:
            if self.label_rotation is None:
                if self.orientation == 'vertical':
                    label_rotation = 270
                elif self.orientation == 'horizontal':
                    label_rotation = 0
                else:
                    raise AssertionError
            else:
                label_rotation = self.label_rotation
            cb.set_label(
                label,
                rotation=label_rotation,
                labelpad=self.labelpad,
                fontsize=self.label_fontsize
            )

        filesnpaths.is_output_file_writable(out_path, ok_if_exists=self.overwrite_output)
        fig.savefig(out_path, format='pdf', bbox_inches='tight')

    def draw_continuous(
        self,
        colormap: mcolors.Colormap,
        vmin: float,
        vmax: float,
        out_path: str,
        label: str = None,
        integer_ticks: bool = False,
        clamped_low: bool = False,
        clamped_high: bool = False,
        center: Union[float, None] = None
    ) -> None:
        """
        Save a standalone continuous colorbar to a file.

        Parameters
        ==========
        colormap : matplotlib.colors.Colormap
            Colormap sampled across the value range.

        vmin : float
            Lower bound of the value range.

        vmax : float
            Upper bound of the value range.

        out_path : str
            Path to PDF output file.

        label : str, None
            Overall colorbar label.

        integer_ticks : bool, False
            If True, label the bar at up to 'max_integer_ticks' whole numbers evenly spanning the
            range, both ends included, rather than at Matplotlib's automatic ticks. Pass this for a
            range that counts things — how many samples contain an element, say — where automatic
            ticks can fall between whole numbers and leave the ends of the range unlabeled, which
            are the two values a reader most needs in order to tell what a color stands for. The
            range is assumed to run between whole numbers, as a count does.

        clamped_low : bool, False
            If True, a value limit truncated the bottom of the range, so 'vmin' is labeled with
            'CLAMPED_MIN_PREFIX': its color is what everything below it is drawn in, and a bare
            number there would claim the color stands for that value alone.

        clamped_high : bool, False
            The same for the top of the range, labeling 'vmax' with 'CLAMPED_MAX_PREFIX'.

        center : Union[float, None], None
            The value the range was centered on, which is ticked and labeled wherever it falls.
            Without it a reader has only the two ends to go by, and a center that is not a round
            number, or a colormap whose middle color is not obvious, leaves the middle of the scale
            impossible to place. A centered range is symmetric about this value, which therefore
            always lies strictly inside it.
        """
        fig = Figure(figsize=self.figsize)
        ax = fig.subplots()

        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
        cb = fig.colorbar(
            ScalarMappable(norm=norm, cmap=colormap),
            cax=ax,
            orientation=self.orientation
        )

        if integer_ticks:
            # Every whole number where they all fit, and otherwise a whole-number stride, so that
            # the gaps between labels are even: spacing the ticks evenly and rounding each to a
            # whole number afterwards would leave gaps of different sizes, which reads as a missing
            # label. The last tick is the top of the range rather than the last stride, so both ends
            # are labeled; only the final gap can come out shorter, as it does on any axis.
            lower = int(vmin)
            upper = int(vmax)
            if upper - lower + 1 <= self.max_integer_ticks:
                ticks = list(range(lower, upper + 1))
            else:
                stride = math.ceil((upper - lower) / (self.max_integer_ticks - 1))
                ticks = list(range(lower, upper, stride)) + [upper]
            cb.set_ticks(ticks)

        if clamped_low or clamped_high or center is not None:
            # A truncated end is labeled with the value it stops at, marked '≤' or '≥' because its
            # color is what everything past the limit is drawn in, and a centered range is labeled
            # at its center, which is otherwise the one place on the bar a reader cannot find. The
            # rest of the bar keeps Matplotlib's own ticks, run through the bar's own formatter, so
            # that such a bar is typeset exactly like a plain one; any of them falling all but on
            # top of one of these labels is dropped, since of the two labels in that spot, this is
            # the one a reader needs. A tick landing on an end that was NOT marked is kept, which is
            # what labels that end of the bar.
            span = vmax - vmin
            marked = (
                ([vmin] if clamped_low else [])
                + ([center] if center is not None else [])
                + ([vmax] if clamped_high else [])
            )
            kept = [
                tick for tick in cb.get_ticks()
                if vmin <= tick <= vmax
                and all(
                    abs(tick - value) >= MIN_TICK_SEPARATION_FRACTION * span for value in marked
                )
            ]
            ticks = sorted(kept + marked)
            formatter = cb.formatter
            formatter.set_locs(ticks)
            if formatter.get_offset():
                # The formatter wants to set a shared offset or multiplier beside the bar, which the
                # fixed labels below would drop, taking the meaning of every tick with it. Spelling
                # the numbers out in full is longer but says what it means without it.
                formatter.set_useOffset(False)
                formatter.set_scientific(False)
                formatter.set_locs(ticks)
            tick_labels = [formatter(tick) for tick in ticks]
            if clamped_low:
                tick_labels[0] = f'{CLAMPED_MIN_PREFIX}{tick_labels[0]}'
            if clamped_high:
                tick_labels[-1] = f'{CLAMPED_MAX_PREFIX}{tick_labels[-1]}'
            cb.set_ticks(ticks)
            cb.set_ticklabels(tick_labels)

        # A continuous bar has no segments to fit labels between, and never more than
        # 'max_integer_ticks' of them, so they take the largest a tick label is allowed — which is
        # what 'draw_discrete' settles on for a bar with as few bands. Left to Matplotlib's default
        # they would be a third the size of the label beside them, and a third the size of the ticks
        # on the discrete bars a single run writes into the same directory.
        cb.ax.tick_params(
            labelsize=self.max_tick_fontsize if self.tick_fontsize is None else self.tick_fontsize
        )

        if label:
            if self.label_rotation is None:
                if self.orientation == 'vertical':
                    label_rotation = 270
                elif self.orientation == 'horizontal':
                    label_rotation = 0
                else:
                    raise AssertionError
            else:
                label_rotation = self.label_rotation
            cb.set_label(
                label,
                rotation=label_rotation,
                labelpad=self.labelpad,
                fontsize=self.label_fontsize
            )

        filesnpaths.is_output_file_writable(out_path, ok_if_exists=self.overwrite_output)
        fig.savefig(out_path, format='pdf', bbox_inches='tight')

class PDFGridDrawer:
    """
    Writes PDF files that are a grid of input PDF files.

    Attributes
    ==========
    overwrite_output : bool
        If True, methods in this class overwrite existing output files.

    paper_format : Union[str, None]
        If None, automatically use the first input PDF file, placed in the upper left of the grid,
        to set the page layout to 'letter-L' (landscape) if aspect ratio > 1, 'letter' if aspect
        ratio ≤ 1. Alternatively, a paper format string can be provided, e.g., 'letter', 'letter-l',
        'A4', 'A4-L'.

    margin : float
        Minimum space between grid cells.

    label_fontsize_scale : Union[float, None]
        If None, the font size of labels over grid cells is set to 80% of the 'margin' argument.
        Alternatively, provide a float on (0.0, 1.0] for the proportion of the minimum margin to use
        as the font size.
    """
    def __init__(self, overwrite_output: bool = FORCE_OVERWRITE) -> None:
        """
        Parameters
        ==========
        overwrite_output : bool, FORCE_OVERWRITE
            If True, methods in this class overwrite existing output files.
        """
        self.overwrite_output = overwrite_output

        self.paper_format: Union[str, None] = None
        self.margin: float = 10.0
        self.label_fontsize_scale: float = 0.8

    def draw(
        self,
        in_paths: Iterable[str],
        out_path: str,
        labels: Iterable[str] = None
    ) -> None:
        """
        Write a PDF containing a grid of input PDF images.

        Parameters
        ==========
        in_paths : Iterable[str]
            Paths to input PDFs.

        out_path : str
            Path to output PDF.

        labels : Iterable[str], None
            Labels displayed over grid cells corresponding to input files.
        """
        assert len(in_paths) > 0
        if labels:
            assert len(in_paths) == len(labels)
        assert 0 < self.label_fontsize_scale <= 1

        # Lay out the page.
        if self.paper_format is None:
            pdf_doc = fitz.open(in_paths[0])
            page = pdf_doc.load_page(0)
            first_input_aspect_ratio = page.rect.width / page.rect.height
            paper_format = 'letter-L' if first_input_aspect_ratio > 1 else 'letter'
        else:
            paper_format = self.paper_format

        # Find the number of rows and columns in the grid.
        cols = math.ceil(math.sqrt(len(in_paths)))
        rows = math.ceil(len(in_paths) / cols)

        # Find the width and height of each cell.
        width, height = fitz.paper_size(paper_format)
        cell_width = (width - (cols + 1) * self.margin) / cols
        cell_height = (height - (rows + 1) * self.margin) / rows

        fontsize = self.margin * self.label_fontsize_scale

        # Create a new PDF document.
        output_doc = fitz.open()
        output_page = output_doc.new_page(width=width, height=height)

        # Loop through input PDF files, placing them in the grid.
        for i, pdf_path in enumerate(in_paths):
            pdf_doc = fitz.open(pdf_path)
            page = pdf_doc.load_page(0)

            # Calculate position in the grid.
            row = i // cols
            col = i % cols
            x = self.margin + col * (cell_width + self.margin)
            y = self.margin + row * (cell_height + self.margin)

            # Resize the input PDF to the cell by the longest dimension, maintaining aspect ratio.
            input_aspect_ratio = page.rect.width / page.rect.height
            if input_aspect_ratio > 1:
                draw_width = cell_width
                draw_height = cell_width / input_aspect_ratio
            else:
                draw_height = cell_height
                draw_width = cell_height * input_aspect_ratio

            # If the resized shorter side still exceeds the cell size, resize by the shorter side.
            if draw_width > cell_width:
                draw_width = cell_width
                draw_height = cell_width / input_aspect_ratio
            if draw_height > cell_height:
                draw_height = cell_height
                draw_width = cell_height * input_aspect_ratio

            # Find upper left drawing coordinates.
            draw_x = x + (cell_width - draw_width) / 2
            draw_y = y + (cell_height - draw_height) / 2

            # Place the input PDF.
            rect = fitz.Rect(draw_x, draw_y, draw_x + draw_width, draw_y + draw_height)
            output_page.show_pdf_page(rect, pdf_doc, 0)

            if labels:
                # Draw labels above each image.
                label = labels[i]
                label_x = draw_x
                label_y = draw_y
                output_page.insert_text((label_x, label_y), label, fontsize=fontsize)

        filesnpaths.is_output_file_writable(out_path, ok_if_exists=self.overwrite_output)
        output_doc.save(out_path)
