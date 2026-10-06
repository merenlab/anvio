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

# Colormaps that anvi'o defines itself, keyed by name. Each is a series of colors spaced evenly from
# one end of the colormap to the other. '_get_colormap' looks for a name here before it looks among
# Matplotlib's colormaps. Every colormap option therefore accepts these names. An '_r' suffix
# reverses one of them. The suffix works the same way for Matplotlib's colormaps.
#
# 'clocktime' colors clock times in hours, from 0 to 24. Its colors are those of hours 0 through 24,
# six hours to a row. The last color repeats the first. Midnight then has the same color at both
# ends of the scale. Noon is yellow. Midnight is purple. The morning runs through blue, teal and
# green. The afternoon and evening run through orange, rose and plum. Lightness rises at an even
# rate from midnight to noon. It falls at the same rate after noon. So two hours equally far from
# noon, such as 6 h and 18 h, are equally light but differ in hue. The colormap is modeled on
# Matplotlib's 'twilight_shifted'. That colormap is nearly white at noon, which is hard to see on a
# white map. It is also nearly black at midnight, which makes black box labels hard to read. Black
# text has a contrast of at least 3:1 on every color of 'clocktime'.
ANVIO_COLORMAPS: Dict[str, Tuple[str, ...]] = {
    'clocktime': (
        '#6e4c80', '#73559c', '#7060b7', '#646ecc', '#577dd5', '#4f8dd7',
        '#519bd3', '#59a8cb', '#5fb6c1', '#74c3ab', '#95cd8a', '#c1d35e',
        '#efd53b', '#f8c13b', '#f9b050', '#f5a15e', '#e9956a', '#df8872',
        '#dd7771', '#d66877', '#c85d81', '#b55686', '#a05286', '#874f85',
        '#6e4c80'
    )
}

# Colormaps for values that repeat, such as clock times. Each has the same color at both ends. On a
# global or overview map drawn from reactions alone, each compound takes a color from the reaction
# lines that touch it. With one of these colormaps, the values of those lines are averaged on a
# circle. A layer with a period ('--reaction-value-period') averages them on a circle with any
# colormap. '_is_cyclic_colormap' matches colormap names against this list. Matplotlib's cyclic
# colormaps are not listed. People often use one of them, 'hsv', for values that do not repeat.
CYCLIC_COLORMAPS: Tuple[str, ...] = (
    'clocktime',
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

# Aggregations of values that repeat after a period, such as clock times in hours, which repeat
# every 24. A layer's period is given by '--reaction-value-period' or '--compound-value-period'.
# Every other aggregation treats values just below the period and just above 0 as far apart
# ('_resolve_aggregation'). Pandas has none of these names for circular aggregations. At the gene
# level, a pandas groupby reduces the rows of each accession ('_aggregate_accession_quantities'); it
# usually takes the aggregation's name, but for these names, it takes anvi'o's own function instead.
# Every level of reduction beyond genes then also uses the same function.
CIRCULAR_AGGREGATIONS: Tuple[str, ...] = (
    'circular_mean',
)

# Summaries of how far values that repeat after a period are spread apart. An example is how far
# apart the peak times of a map element are in different samples. They need a period, as
# 'CIRCULAR_AGGREGATIONS' do. A spread is a plain number, not a point in the period. Only the
# summary that colors the 'unified' map ('--*-sample-summary' without groups or '--*-group-summary'
# with groups) takes these names. The two kinds of map then differ ('_build_txt_model'). Each
# individual map shows a point in the period for each element, its circular mean in one sample
# or group. Its scale runs from 0 to the period, with the cyclic 'clocktime' by default, and its
# colorbar is labeled with the value column's name. The 'unified' map shows the spread of those
# points. Its scale runs from 0 to the largest spread there can be, 5.40 with a period of 24. It
# takes a colormap that is not cyclic, 'plasma_r' by default, since the smallest and the largest
# spread are not the same point. Its colorbar is labeled '{value} angular deviation'. On global and
# overview maps, its compound colors are averaged as plain numbers, not on a circle.
CIRCULAR_SPREADS: Tuple[str, ...] = (
    'angular_deviation',
)

# The fewest values whose angular deviation is defined ('_angular_deviation'). The spread of two
# values is only the gap between them. Fewer values leave the element uncolored.
MIN_ANGULAR_DEVIATION_VALUES = 3

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

# The statistic a normalization compares each element's value against, computed over the categories
# in which that element has a value at all. These reduce a list like
# 'AGGREGATION_FUNCTIONS' do: they are applied to one element's values across the samples or groups
# rather than through pandas.
ELEMENT_NORMALIZATION_REFERENCES = {
    'mean': lambda values: float(np.mean(values)),
    'median': lambda values: float(np.median(values)),
    'max': lambda values: float(np.max(values)),
    'total': lambda values: float(np.sum(values))
}

# The normalizations of '--*-element-normalization', which rescale a map element's value in each
# sample or group against the same element's values across all of them. This transformation is not a
# reduction like an aggregation or a summary. Normalized values are displayed on sample or group
# maps. 'form' is how a value is compared to the reference ('_element_normalization_function'), and
# 'reference' is the statistic it is compared against, one of 'ELEMENT_NORMALIZATION_REFERENCES',
# 'circular_mean' for 'CIRCULAR_ELEMENT_NORMALIZATIONS', or None for a form needing none. 'centered'
# says whether zero is the neutral middle of the resulting quantity, which triggers a diverging
# colormap and a scale centered on zero by default. 'label' names the quantity on the colorbar, with
# '{value}' standing in for the value column's own name.
ELEMENT_NORMALIZATIONS = {
    'relative_to_mean': {
        'form': 'relative', 'reference': 'mean', 'centered': True,
        'label': '{value} relative to mean'
    },
    'relative_to_median': {
        'form': 'relative', 'reference': 'median', 'centered': True,
        'label': '{value} relative to median'
    },
    'difference_from_mean': {
        'form': 'difference', 'reference': 'mean', 'centered': True,
        'label': '{value} - mean'
    },
    'difference_from_median': {
        'form': 'difference', 'reference': 'median', 'centered': True,
        'label': '{value} - median'
    },
    'log2_ratio_to_mean': {
        'form': 'log2_ratio', 'reference': 'mean', 'centered': True,
        'label': 'log2({value} / mean)'
    },
    'log2_ratio_to_median': {
        'form': 'log2_ratio', 'reference': 'median', 'centered': True,
        'label': 'log2({value} / median)'
    },
    'z_score': {
        'form': 'z_score', 'reference': 'mean', 'centered': True,
        'label': '{value} z-score'
    },
    'rank': {
        'form': 'rank', 'reference': None, 'centered': False,
        'label': '{value} rank'
    },
    'fraction_of_max': {
        'form': 'fraction', 'reference': 'max', 'centered': False,
        'label': '{value} / max'
    },
    'fraction_of_total': {
        'form': 'fraction', 'reference': 'total', 'centered': False,
        'label': '{value} / total'
    },
    'difference_from_circular_mean': {
        'form': 'circular_difference', 'reference': 'circular_mean', 'centered': True,
        'label': '{value} - circular mean'
    }
}

# The normalization names as a phrase for messages that list them.
ELEMENT_NORMALIZATION_PHRASE = ', '.join(repr(name) for name in ELEMENT_NORMALIZATIONS)

# Normalizations of values that repeat after a period. Each value becomes its signed offset from the
# element's circular mean, wrapped to [-period / 2, period / 2). The two ends of that range are the
# same point. The circular mean comes from '_circular_mean', not from
# 'ELEMENT_NORMALIZATION_REFERENCES'; the other normalizations take their reference from that table,
# such as the plain mean, and compare values on a line. Circular normalization names need a period,
# and a layer with a period takes only these names ('_resolve_element_normalization').
CIRCULAR_ELEMENT_NORMALIZATIONS: Tuple[str, ...] = (
    'difference_from_circular_mean',
)

# The forms whose reference has to be positive for the comparison to mean anything: dividing by a
# reference of zero has no scale to measure against, and dividing by a negative one inverts the sign
# of every deviation, so that a value above the reference would be drawn as though it were below it.
# An element whose reference is not positive gets no value under these forms and is left uncolored.
ELEMENT_NORMALIZATION_POSITIVE_FORMS = ('relative', 'log2_ratio', 'fraction')

# The colormap the per-sample or per-group scale takes when a normalization centered on zero colors
# it and no colormap was named for it. A normalized value says which side of the element's own
# reference a sample falls on, which a diverging colormap shows, unlike the sequential default.
DEFAULT_CENTERED_COLORMAP = 'RdYlGn'

# The colormap of a layer whose values repeat after a period, when no colormap was named for it. 0
# and the period are the same point, and a cyclic colormap gives them the same color. 'clocktime' is
# laid out for hours, but it is cyclic for any period.
DEFAULT_PERIOD_COLORMAP = 'clocktime'

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

# What marks the label of a value limit that values lie past ('_make_quantitative_norm'). Everything
# beyond the limit is drawn in the color at that end of the scale. The color there therefore stands
# for that value or anything past it. A bare number would claim it means that value exactly.
CLAMPED_MIN_PREFIX = '≤ '
CLAMPED_MAX_PREFIX = '≥ '

# How close to a limit or center labeled on a continuous colorbar, as a fraction of the bar, an
# automatic tick may fall before it is dropped ('draw_continuous'). The two labels would otherwise
# be set almost on top of each other. The label of the limit or center is the one a reader needs, so
# it is the one that stays.
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
        directly from the gene annotations in the JSON, so no reference databases are required. In
        a JSON from a pangenome, the genes are gene clusters, and each one has a single KO.

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
            ko_annotation = gene_entry.get('annotation', {}).get('ko', {})
            if 'id' in ko_annotation:
                # A JSON from a pangenome lists gene clusters. Each has one KO, given by 'id'.
                ko_ids.add(ko_annotation['id'])
            else:
                # A JSON from a contigs database lists genes. Their KOs are the keys.
                ko_ids.update(ko_annotation)

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
    def _resolve_aggregation(
        aggregation: str,
        flag: str,
        period: Union[float, None] = None,
        period_flag: Union[str, None] = None
    ) -> Callable:
        """
        Resolve an aggregation name into a function reducing a sequence of values to one value.

        Values that repeat after a period are reduced only by the names in 'CIRCULAR_AGGREGATIONS'.
        These names need the period. Both rules are checked here.

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

        period : Union[float, None], None
            The period after which the layer's values repeat, or None if they do not repeat.

        period_flag : Union[str, None], None
            The command-line flag that gives the layer's period, used in an error message.

        Returns
        =======
        Callable
            Reduces a sequence of values to a single value.
        """
        if aggregation in CIRCULAR_AGGREGATIONS + CIRCULAR_SPREADS:
            if period is None:
                raise ConfigError(
                    f"'{flag}' was given as '{aggregation}', which is for values that repeat after "
                    f"a period, such as clock times. No period was given for these values. Give "
                    f"one with '{period_flag}', such as 24 for clock times in hours."
                )
            circular_function = (
                Mapper._angular_deviation if aggregation in CIRCULAR_SPREADS
                else Mapper._circular_mean
            )

            def circular_aggregate(values):
                return circular_function(values, period)
            return circular_aggregate
        if period is not None:
            raise ConfigError(
                f"'{flag}' was given as '{aggregation}', but '{period_flag}' says that the values "
                f"repeat every {period:g}. Values just below {period:g} and just above 0 are then "
                f"close together. '{aggregation}' treats them as far apart. Use "
                f"{', '.join(repr(name) for name in CIRCULAR_AGGREGATIONS)} instead."
            )

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
    def _circular_mean(values: Iterable[float], period: float) -> float:
        """
        Average values that repeat after a period, such as clock times, on a circle.

        Each value is a point on a circle whose circumference is the period. 0 and the period are
        the same point. The mean is the direction of the average of these points. With a period of
        24, the mean of 23.5 and 0.5 is 0, not 12. The mean lies in [0, period).

        Values that cancel out have no mean direction. Examples are 6 and 18, or 0, 8 and 16, with a
        period of 24. The mean is then undefined, and the result is NaN. The value is then dropped
        as any undefined aggregation is ('_finite_values', '_reduce_entry_value').

        Parameters
        ==========
        values : Iterable[float]
            The values to average. There is at least one. They can be a pandas Series, as a groupby
            passes, so they are read by position.

        period : float
            The positive period after which the values repeat.

        Returns
        =======
        float
            The mean, in [0, period), or NaN where it is undefined.
        """
        # Each value is first taken modulo the period. The angles are then exact, however far from
        # 0 the values lie.
        values = np.mod(np.asarray(values, dtype=float), period)
        # Equal values return that value. A single value then comes back exactly, with no rounding
        # from the trigonometry.
        if np.all(values == values[0]):
            mean = float(values[0])
        else:
            angles = 2 * np.pi * values / period
            sine = np.mean(np.sin(angles))
            cosine = np.mean(np.cos(angles))
            # The tolerance is the one kgml.Pathway uses to find compound colors that cancel out.
            if np.hypot(sine, cosine) < 1e-9:
                return float('nan')
            # The fraction of the circle is rounded as kgml.Pathway rounds it. This removes
            # floating-point error. Without it, a mean of 0 could come out a hair below the period.
            # Those are the same point, but the two ends of a scale that is not cyclic.
            mean = (round(float(np.arctan2(sine, cosine) / (2 * np.pi)), 12) % 1.0) * period
        # A value a hair below 0, taken modulo the period, rounds to the period itself. That is the
        # same point as 0.
        return 0.0 if mean >= period else mean

    @staticmethod
    def _angular_deviation(values: Iterable[float], period: float) -> float:
        """
        Measure how far values that repeat after a period are spread apart on a circle.

        Each value is a point on a circle whose circumference is the period, as in
        '_circular_mean'. R is the length of the average of these points. R is 1 where the values
        are equal, and 0 where they cancel out. The angular deviation is sqrt(2 * (1 - R)), in units
        of the period divided by 2 * pi. It runs from 0, where the values are equal, to
        '_largest_angular_deviation', where they cancel out. With a period of 24, that is 5.40. For
        values close together, it is close to their standard deviation.

        Fewer than 'MIN_ANGULAR_DEVIATION_VALUES' values give NaN. The value is then dropped as any
        undefined aggregation is ('_finite_values', '_reduce_entry_value').

        Parameters
        ==========
        values : Iterable[float]
            The values whose spread is measured.

        period : float
            The positive period after which the values repeat.

        Returns
        =======
        float
            The angular deviation, or NaN where there are too few values.
        """
        values = np.mod(np.asarray(values, dtype=float), period)
        if len(values) < MIN_ANGULAR_DEVIATION_VALUES:
            return float('nan')
        angles = 2 * np.pi * values / period
        resultant = np.hypot(np.mean(np.sin(angles)), np.mean(np.cos(angles)))
        # Rounding can put R a hair above 1 where the values are equal.
        return float(np.sqrt(max(0.0, 2 * (1 - resultant))) * period / (2 * np.pi))

    @staticmethod
    def _largest_angular_deviation(period: float) -> float:
        """
        Return the angular deviation of values that cancel out, the largest there can be.

        Parameters
        ==========
        period : float
            The positive period after which the values repeat.

        Returns
        =======
        float
            period * sqrt(2) / (2 * pi), which is 5.40 for a period of 24.
        """
        return period * np.sqrt(2) / (2 * np.pi)

    @staticmethod
    def _element_normalization_function(
        form: str,
        reference: Union[str, None],
        period: Union[float, None] = None
    ) -> Callable:
        """
        Build the function behind one of 'ELEMENT_NORMALIZATIONS'.

        The function takes one map element's values in the categories that have one, in order, and
        returns a new value for each of them, in the same order. A value it cannot define is
        returned as 'nan' or as an infinity, either of which leaves that element uncolored on that
        category's map, exactly as an undefined aggregation leaves it uncolored
        ('_reduce_entry_value').

        Parameters
        ==========
        form : str
            How a value is compared to the reference: the 'form' of an 'ELEMENT_NORMALIZATIONS'
            entry.

        reference : Union[str, None]
            The statistic to compare against, named among 'ELEMENT_NORMALIZATION_REFERENCES', or
            'circular_mean', or None for a form that needs no reference.

        period : Union[float, None], None
            The period after which the values repeat. A circular reference and the
            'circular_difference' form need it.

        Returns
        =======
        Callable
            Maps a sequence of values to a NumPy array of normalized values of the same length.
        """
        needs_positive = form in ELEMENT_NORMALIZATION_POSITIVE_FORMS
        if reference in CIRCULAR_AGGREGATIONS:
            def statistic(values):
                return Mapper._circular_mean(values, period)
        else:
            statistic = None if reference is None else ELEMENT_NORMALIZATION_REFERENCES[reference]

        def normalize(values) -> np.ndarray:
            # NumPy is asked not to raise or warn on the divisions and logarithms below, so that a
            # reference or a value that leaves one element undefined is returned as 'nan' or an
            # infinity for that element alone rather than interrupting a run over thousands of them.
            with np.errstate(all='ignore'):
                array = np.asarray(values, dtype=float)
                if form == 'rank':
                    # Ties share the average of the ranks they span, as pandas ranks by default.
                    return pd.Series(array).rank().to_numpy()
                undefined = np.full(array.shape, np.nan)
                center = statistic(array)
                if not np.isfinite(center) or (needs_positive and center <= 0):
                    return undefined
                if form == 'relative':
                    return (array - center) / center
                if form == 'difference':
                    return array - center
                if form == 'circular_difference':
                    # The offset is the shorter way around the circle, wrapped to
                    # [-period / 2, period / 2). It is rounded as '_circular_mean' rounds, so that
                    # an offset of half a period always lands on the same end.
                    fraction = np.round((array - center) / period + 0.5, 12) % 1.0 - 0.5
                    return fraction * period
                if form == 'log2_ratio':
                    return np.log2(array / center)
                if form == 'fraction':
                    return array / center
                if form == 'z_score':
                    # The spread of a single value is undefined, and a spread of zero would divide
                    # deviations that are themselves all zero by it.
                    if array.size < 2:
                        return undefined
                    spread = float(np.std(array, ddof=1))
                    if not np.isfinite(spread) or spread == 0:
                        return undefined
                    return (array - center) / spread
                raise ConfigError(
                    f"_element_normalization_function :: '{form}' is not a form of comparison "
                    f"anvi'o knows how to make. This is worth reporting as a bug."
                )

        return normalize

    @staticmethod
    def _resolve_element_normalization(
        normalization: str,
        flag: str,
        element_type: str,
        period: Union[float, None] = None,
        period_flag: Union[str, None] = None
    ) -> Tuple[Callable, str, bool]:
        """
        Resolve a normalization name into a function rescaling one element's per-category values.

        The recommended names are the keys of 'ELEMENT_NORMALIZATIONS'. Any other name is taken to
        be a pandas Series method that transforms the values it is given into one new value each,
        such as 'abs', and is probed here, before any data is touched, so that a name pandas does
        not recognize, or one that does something other than transform, is reported as a
        configuration error rather than failing deep in a drawing loop. This mirrors
        '_resolve_aggregation' one rung up, and the two are opposites: an aggregation must reduce a
        set of values to one number and so refuses a name that transforms, while a normalization
        must return one value per category and so refuses a name that reduces.

        Probes ask four things of a name. It must return a Series of one value per category, with
        the same index in the same order — a name that reorders or drops categories would hand a
        sample another sample's value undetected. It must return numbers, which alone can go on a
        color scale. It must give the same answer twice for the same values, since the range of the
        scale and the colors on the maps are worked out in separate passes. And it must not depend
        on the order the categories come in, which is an artifact of how the samples happen to be
        named and listed and says nothing about the data: 'cumsum' and 'ffill' answer differently
        for a sample depending on what precedes it, so are refused. What a name MAKES of the values
        is its own business, as it is for an aggregation: a name that transforms each value without
        reference to the others, such as 'abs', is accepted, and the colorbar says which name was
        used.

        Parameters
        ==========
        normalization : str
            The requested normalization name.

        flag : str
            The command-line flag the name comes from, used in error messages.

        element_type : str
            'reaction' or 'compound', naming the layer in error messages.

        period : Union[float, None], None
            The period after which the layer's values repeat, or None if they do not repeat. Values
            with a period take only the names of 'CIRCULAR_ELEMENT_NORMALIZATIONS', and these names
            need a period. Both rules are checked here.

        period_flag : Union[str, None], None
            The command-line flag that gives the layer's period, used in error messages.

        Returns
        =======
        Tuple[Callable, str, bool]
            The normalizing function, the colorbar label with '{value}' standing in for the value
            column's name, and whether zero is the neutral middle of the quantity it makes.
        """
        if normalization in CIRCULAR_ELEMENT_NORMALIZATIONS and period is None:
            raise ConfigError(
                f"'{flag}' was given as '{normalization}', which is for values that repeat after a "
                f"period, such as clock times. No period was given for these values. Give one "
                f"with '{period_flag}', such as 24 for clock times in hours."
            )
        if period is not None and normalization not in CIRCULAR_ELEMENT_NORMALIZATIONS:
            circular_names = ', '.join(repr(name) for name in CIRCULAR_ELEMENT_NORMALIZATIONS)
            raise ConfigError(
                f"'{flag}' was given as '{normalization}', but '{period_flag}' says that the "
                f"values repeat every {period:g}. Values just below {period:g} and just above 0 "
                f"are then close together. '{normalization}' compares values on a line, so it "
                f"treats them as far apart. Use {circular_names} instead."
            )
        try:
            preset = ELEMENT_NORMALIZATIONS[normalization]
        except KeyError:
            pass
        else:
            return (
                Mapper._element_normalization_function(
                    preset['form'], preset['reference'], period
                ),
                preset['label'],
                preset['centered']
            )

        def refuse(reason: str) -> None:
            raise ConfigError(
                f"'{flag}' was given as '{normalization}', which {reason}. It takes one of "
                f"anvi'o's normalizations ({ELEMENT_NORMALIZATION_PHRASE}), or the name of any "
                f"pandas Series method that transforms the values it is given into one new value "
                f"each, such as 'abs'. A name that reduces a set of values to a single number, "
                f"such as 'mean', belongs to '--{element_type}-sample-summary' instead, which is "
                f"what summarizes samples."
            )

        if not hasattr(pd.Series(dtype=float), normalization):
            refuse("is neither one of anvi'o's normalizations nor a method pandas offers a Series")

        def apply(series: pd.Series) -> pd.Series:
            return getattr(series, normalization)()

        def as_floats(result: pd.Series) -> np.ndarray:
            # A nullable dtype ('convert_dtypes') carries a missing value that no array of plain
            # floats can hold, which is a refusal rather than an error raised from the probes.
            try:
                return result.to_numpy(dtype=float)
            except (TypeError, ValueError):
                refuse(
                    "returns values of a kind that cannot be made into the plain numbers a color "
                    "scale is drawn from"
                )

        # The probes give the categories names rather than positions, and names whose alphabetical
        # order differs from the order they are in, so that a name sorting by either the values or
        # the index ('sort_values', 'sort_index') disagrees with the order it was handed. Ties in
        # the values unmask a name that quietly drops repeated values ('drop_duplicates'). A missing
        # value unmasks a name that drops categories ('dropna') or fills one from its neighbors
        # ('ffill'), which are no-ops on a complete set of values. The last two probes are a pair and
        # a single value: two categories are the fewest a normalization across them can say
        # anything about, and one is what an element found in a single sample gives it, which is
        # where a name that collapses to a scalar ('squeeze') would otherwise fail mid-drawing.
        probe_index = ['SAMPLE_2', 'SAMPLE_1', 'SAMPLE_10', 'SAMPLE_3', 'SAMPLE_20', 'SAMPLE_4']
        probe_values = [3.5, 11.25, 0.75, 7.0, 2.25, 19.5]
        tied_values = [3.0, 3.0, 7.0, 7.0, 1.0, 9.0]
        missing_values = [3.5, float('nan'), 0.75, 7.0, 2.25, 19.5]
        for values in (
            probe_values, tied_values, missing_values, probe_values[:2], probe_values[:1]
        ):
            index = probe_index[:len(values)]

            try:
                result = apply(pd.Series(values, index=index, dtype=float))
            except Exception:
                refuse("pandas could not apply to a set of category values")

            if isinstance(result, (float, int, np.number)) or np.isscalar(result):
                refuse(
                    "reduces a set of values to a single number rather than giving each of them a "
                    "new value"
                )

            if not isinstance(result, pd.Series):
                refuse("does not give each category a value of its own")

            if list(result.index) != index:
                refuse(
                    "does not return the values of the same categories in the same order, so "
                    "anvi'o cannot tell to which of them each value belongs"
                )

            if not pd.api.types.is_numeric_dtype(result):
                refuse("does not return numbers, so its results cannot be put on a color scale")

            # The range of the scale and the colors of the maps are worked out in two passes over
            # the maps, so a name that answers differently each time would color the maps on a scale
            # that does not span them.
            again = apply(pd.Series(values, index=index, dtype=float))
            if not np.allclose(
                as_floats(result), as_floats(again), rtol=1e-12, atol=0, equal_nan=True
            ):
                refuse("does not give the same answer twice for the same values")

        # Reversing the categories, rotating them by one, and swapping the first two catch a
        # dependence on their order that any single rearrangement might leave looking the same. The
        # tied values are probed as well, since a name can answer the same for values that all
        # differ and still turn on which of two equal ones comes first ('duplicated'), and so is
        # the missing value, since a name that fills a category from its neighbors only reveals
        # that it does so when there is something to fill.
        for values in (probe_values, tied_values, missing_values):
            index = probe_index[:len(values)]
            probe = pd.Series(values, index=index, dtype=float)
            given = as_floats(apply(probe))
            positions = list(range(len(values)))
            for order in (
                positions[::-1], positions[1:] + positions[:1], [1, 0] + positions[2:]
            ):
                try:
                    rearranged = apply(probe.iloc[order]).reindex(index)
                except Exception:
                    refuse("pandas could not apply to a set of category values")

                if not np.allclose(
                    given, as_floats(rearranged), rtol=1e-12, atol=0, equal_nan=True
                ):
                    refuse(
                        "answers differently when the categories are given in a different order, "
                        "and their order says nothing about the data — it is an artifact of how "
                        "the samples happen to be named — so its results would depend on those "
                        "names"
                    )

        def normalize(values) -> np.ndarray:
            # The categories are numbered here rather than named, since a name that survived the
            # probes above answers the same whatever their labels are, and numbering is cheaper.
            result = apply(pd.Series(np.asarray(values, dtype=float)))
            if len(result) != len(values):
                raise ConfigError(
                    f"'{flag}' was given as '{normalization}', which passed anvi'o's checks but "
                    f"has now returned {len(result)} values for a map element with "
                    f"{len(values)} of them. This is worth reporting as a bug."
                )
            return result.to_numpy(dtype=float)

        return normalize, '{value} (' + normalization + ')', False

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
        # Example accessions are given only for a circular mean, as in '_warn_undefined_summaries'.
        examples = ''
        if aggregation in CIRCULAR_AGGREGATIONS:
            examples = f", including these: {', '.join(sorted(undefined)[:5])}"
        self.run.warning(
            f"Reducing the values of the text file at '{path}' with '{aggregation}' was undefined "
            f"for {len(undefined)} accession(s){examples}. {Mapper._undefined_cause(aggregation)} "
            f"These accessions are treated as having no value, so the map elements that depend on "
            f"them are left uncolored."
        )

    @staticmethod
    def _undefined_cause(aggregation: str, noun: str = 'an aggregation') -> str:
        """
        Say why an aggregation or summary can be undefined, for a warning that names it.

        Parameters
        ==========
        aggregation : str
            The name of the aggregation or summary.

        noun : str, 'an aggregation'
            What the warning calls it, 'an aggregation' or 'a summary'.

        Returns
        =======
        str
            A sentence giving the cause.
        """
        if aggregation in CIRCULAR_AGGREGATIONS:
            return (
                "A circular mean is undefined where the values cancel out, such as 6 and 18 with a "
                "period of 24."
            )
        return (
            f"This happens when {noun} needs more values than are available, as the standard "
            f"deviation does for a single value."
        )

    def _warn_undefined_summaries(
        self,
        layer: dict,
        draw_unified_maps: bool,
        draw_category_maps: bool
    ) -> None:
        """
        Report the map elements of one layer whose sample or group summary was undefined.

        '_summarize_entry_values' collects these elements while the range of values is found. A
        summary is reported only where it changes a map that is drawn. Ungrouped, the sample summary
        colors the 'unified' map. Grouped, the sample summary colors the map of each group. The
        group summary then pools the groups' values on the 'unified' map. So a sample summary that
        is undefined in a group also changes the 'unified' map where the group summary pools values.

        Parameters
        ==========
        layer : dict
            A 'quantitative' layer model with samples.

        draw_unified_maps : bool
            True if the 'unified' map is drawn.

        draw_category_maps : bool
            True if the maps of the individual samples or groups are drawn.
        """
        element_type = layer['element_type']
        grouped = layer['group_samples'] is not None
        reports = []
        if grouped:
            consequences = []
            if draw_category_maps:
                consequences.append("It is left uncolored on that group's map.")
            if draw_unified_maps and layer['group_aggregate'] is not None:
                consequences.append(
                    "The 'unified' map pools only the groups where the element has a value. An "
                    "element with no value in any group is left uncolored there."
                )
            if consequences:
                reports.append((
                    layer['_undefined_category'], f'--{element_type}-sample-summary',
                    layer['sample_summary'], ' in at least one group',
                    ' '.join(
                        ["An element has no value in a group where its summary is undefined."]
                        + consequences
                    )
                ))
        if draw_unified_maps:
            reports.append((
                layer['_undefined_unified'],
                f"--{element_type}-{'group' if grouped else 'sample'}-summary",
                layer['group_summary'] if grouped else layer['sample_summary'], '',
                "These elements are left uncolored on the 'unified' map."
            ))
        for undefined, flag, summary, where, consequence in reports:
            if not undefined:
                continue
            # Example accessions are given only for a circular mean. It is undefined where the
            # values of an element cancel out, so the examples name elements worth checking. Other
            # summaries, such as 'std', are undefined where an element has too few values. Their
            # examples would only be the first accessions in alphabetical order. An angular
            # deviation needs a fixed number of values, so its warning states that number.
            if summary in CIRCULAR_SPREADS:
                self.run.warning(
                    f"The summary '{flag} {summary}' needs values from at least "
                    f"{MIN_ANGULAR_DEVIATION_VALUES} {'groups' if grouped else 'samples'}. "
                    f"{len(undefined)} distinct map element(s) have fewer, so they are left "
                    f"uncolored on the 'unified' map. The spread of two values is only the gap "
                    f"between them."
                )
                continue
            examples = ''
            if summary in CIRCULAR_AGGREGATIONS:
                examples = ', '.join(sorted({
                    kegg_id for key in undefined for kegg_id in key
                    if kegg_id in layer['accessions']
                })[:5])
                examples = f"Their accessions in the file include these: {examples}. "
            self.run.warning(
                f"The summary '{flag} {summary}' was undefined{where} for {len(undefined)} "
                f"distinct map element(s). {examples}"
                f"{Mapper._undefined_cause(summary, 'a summary')} {consequence}"
            )

    def _aggregate_accession_quantities(
        self,
        rows_df: pd.DataFrame,
        value_column: str,
        aggregation: Union[str, Callable],
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
        A name of 'CIRCULAR_AGGREGATIONS' is not a pandas name. It comes here as the function that
        '_resolve_aggregation' made for it, which every other level applies too.

        Parameters
        ==========
        rows_df : pandas.DataFrame
            Rows of one layer, carrying the normalized '__accession' column and 'value_column'.

        value_column : str
            Name of the auto-detected numeric value column.

        aggregation : Union[str, Callable]
            How to reduce an accession's rows to a single value, as a pandas aggregation name or as
            the function for a circular aggregation.

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

    def _check_category_names(self, categories: Iterable[str], category_noun: str) -> None:
        """
        Check that category names can serve as output subdirectory names.

        A category (sample, source, or group) drawn on its own maps gets its own subdirectory of the
        output directory, and once map grids are drawn the subdirectories of categories that were
        only needed for a grid are deleted ('_draw_map_grids'). Category names must be checked to
        ensure that they form proper paths that lie within the output directory: names come from a
        text file's 'sample' column or a groups file, so they cannot be trusted to be safe paths
        (for example, a sample named "Station 5/Depth 10" is not safe). The same names become file
        names where the maps are gathered by map rather than by category ('_collate_maps_by_map'),
        which the same requirement makes safe. Only categories getting individual maps as
        independent files or as part of a grid are checked, not those that are only summarized on a
        'unified' map.

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
        plural = len(problems) > 1
        raise ConfigError(
            f"Each {category_noun} drawn on its own maps gets its own subdirectory of the output "
            f"directory, so its name must be usable as a directory name. "
            f"{'These names' if plural else 'This name'} cannot be used: {', '.join(problems)}. "
            f"Please rename {'them' if plural else 'it'} in the input, using single words without "
            f"path separators, such as 'SAMPLE_1' or 'HIGH_TEMPERATURE'."
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
        path: str,
        period: Union[float, None] = None,
        period_flag: Union[str, None] = None
    ) -> Tuple[Literal['value', 'presence'], Union[Callable, None], Union[str, None]]:
        """
        Resolve a sample or group summary into a coloring kind and its reduction.

        A summary reduces a set of samples, or a set of sample groups, to one statement per map
        element. The names of 'SUMMARY_PRESENCE_SCHEMES' summarize presence: how many
        ('count'/'count_continuous'), or exactly which ('membership'), categories contain the
        element. Any other name pools the element's values in the categories as an aggregation
        ('_resolve_aggregation', '_summarize_entry_values'). This requires the layer's file to have
        a value column. The default of None summarizes presence with the colormap scheme left
        unresolved, so that '_membership_layer_colors' picks it from the number of categories (by
        membership for ≤ 3, by count > 3, and by a continuous count scale where a discrete one would
        run out of distinguishable colors or of room to label its bands, the latter past
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

        period : Union[float, None], None
            The period after which the layer's values repeat, or None ('_resolve_aggregation').

        period_flag : Union[str, None], None
            The command-line flag that gives the layer's period, used in an error message.

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
        return 'value', Mapper._resolve_aggregation(summary, flag, period, period_flag), None

    @staticmethod
    def _resolve_value_limits(
        limits: Union[Tuple[Union[float, None], Union[float, None]], None], flag: str
    ) -> Union[Tuple[Union[float, None], Union[float, None]], None]:
        """
        Validate the limits a color scale of values may span.

        A limit sets that end of a scale, whether or not any value reaches it
        ('_make_quantitative_norm'). Either end can be left open with None. The values then set it.
        A pair leaving both ends open would do nothing at all, which is likelier a mistake than an
        intention, and so is refused rather than accepted in silence.

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
    def _resolve_value_period(period: Union[float, str, None], flag: str) -> Union[float, None]:
        """
        Resolve the period after which a layer's values repeat, as typed or as a caller gave it.

        Parameters
        ==========
        period : Union[float, str, None]
            The period, or None if the values do not repeat.

        flag : str
            The command-line flag the period comes from, used in an error message.

        Returns
        =======
        Union[float, None]
            The period as a positive finite number, or None.
        """
        if period is None:
            return None
        try:
            number = float(period)
        except (TypeError, ValueError):
            number = None
        if number is None or not np.isfinite(number) or number <= 0:
            raise ConfigError(
                f"'{flag}' must be a positive number, such as 24 for clock times in hours. It is "
                f"how far the values go before they repeat. Anvi'o got this instead: '{period}'."
            )
        return number

    @staticmethod
    def _format_accepted(number: float, accepted: Callable[[float], bool]) -> str:
        """
        Format a number that a message suggests, in as few digits as the check on it accepts.

        A user may type a suggested number back in as it was printed. Six significant digits, as
        ':g' prints, can round it to a value that the same check refuses. More digits are then used,
        up to the full precision of a float.

        Parameters
        ==========
        number : float
            The number to suggest.

        accepted : Callable[[float], bool]
            Whether the check the number is suggested for accepts a value.

        Returns
        =======
        str
            The number in the fewest significant digits, from six up, that the check accepts.
        """
        for precision in range(6, 18):
            text = f'{number:.{precision}g}'
            if accepted(float(text)):
                return text
        return repr(number)

    @staticmethod
    def _resolve_value_center(
        center: Union[float, str, None],
        flag: str,
        limits: Union[Tuple[Union[float, None], Union[float, None]], None] = None,
        limits_flag: Union[str, None] = None,
        from_normalization: bool = False
    ) -> Union[float, None]:
        """
        Validate the value a color scale is centered on.

        A centered scale runs the same distance either side of this value
        ('_make_quantitative_norm'), which puts it at the middle of the colormap however lopsided
        the values around it are. A center lying outside the limits bounding the same scale is
        refused here, before any value is read: those limits declare it out of the scale's reach, so
        centering that scale on it cannot be what was meant. A center on a limit is refused too.
        Neither end of a centered scale can lie on its center. A center between two limits is
        refused unless it is their midpoint. Each limit sets its end of the scale, so the two leave
        nothing for centering to widen.

        Parameters
        ==========
        center : Union[float, str, None]
            The value requested at the middle of the scale, or None for a scale left where its
            values and any limits set it.

        flag : str
            The command-line flag the center comes from, used in error messages.

        limits : Union[Tuple[Union[float, None], Union[float, None]], None], None
            The limits bounding the same scale, as '_resolve_value_limits' returns them, or None
            where the scale is not limited.

        limits_flag : Union[str, None], None
            The command-line flag those limits come from, used in an error message.

        from_normalization : bool, False
            If True, a normalization supplied the center rather than the user. A message then does
            not suggest moving the center. The center is the normalization's neutral value. Nobody
            gave it by hand.

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

            def advice(end_noun: str) -> str:
                # A center from a normalization was not given by hand, so only the limit is named.
                return (
                    f"Please move the {end_noun}, or leave that end open." if from_normalization
                    else "Please move one of the two."
                )

            if limit_min is not None and number < limit_min:
                raise ConfigError(
                    f"'{flag}' centers a color scale on {number:g}, while '{limits_flag}' gives "
                    f"that same scale a minimum of {limit_min:g}, which lies above it. A scale "
                    f"cannot be centered on a value it is not allowed to reach. "
                    f"{advice('minimum')}"
                )
            if limit_max is not None and number > limit_max:
                raise ConfigError(
                    f"'{flag}' centers a color scale on {number:g}, while '{limits_flag}' gives "
                    f"that same scale a maximum of {limit_max:g}, which lies below it. A scale "
                    f"cannot be centered on a value it is not allowed to reach. "
                    f"{advice('maximum')}"
                )
            for limit, end_noun in ((limit_min, 'minimum'), (limit_max, 'maximum')):
                if limit is not None and number == limit:
                    raise ConfigError(
                        f"'{flag}' centers a color scale on {number:g}, while '{limits_flag}' sets "
                        f"the {end_noun} of that same scale at {limit:g} as well. A centered scale "
                        f"runs the same distance either side of its center. Neither of its ends "
                        f"can lie on the center. {advice(end_noun)}"
                    )
            # The midpoint is compared within a tolerance scaled to the span of the limits. A center
            # typed as the midpoint of two decimals can differ from their computed midpoint in the
            # last digit.
            if limit_min is not None and limit_max is not None:
                midpoint = (limit_min + limit_max) / 2
                tolerance = 1e-9 * (limit_max - limit_min)
                if abs(number - midpoint) > tolerance:
                    # The given numbers are printed exactly, so that none of them prints the same
                    # as the midpoint without being it.
                    number_text, min_text, max_text = (
                        Mapper._format_accepted(given, lambda value, given=given: value == given)
                        for given in (number, limit_min, limit_max)
                    )
                    midpoint_text = Mapper._format_accepted(
                        midpoint, lambda value: abs(value - midpoint) <= tolerance
                    )
                    move_clause = (
                        '' if from_normalization else f"move the center to {midpoint_text}, "
                    )
                    raise ConfigError(
                        f"'{flag}' centers a color scale on {number_text}, while '{limits_flag}' "
                        f"sets both ends of that same scale, at {min_text} and {max_text}. A "
                        f"centered scale runs the same distance either side of its center. With "
                        f"both ends set, it can only be centered on their midpoint, "
                        f"{midpoint_text}. Please {move_clause}give limits that lie the same "
                        f"distance either side of {number_text}, or leave one of the two ends open."
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

    @staticmethod
    def _is_cyclic_colormap(cmap: mcolors.Colormap) -> bool:
        """
        Report whether a colormap is for values that repeat, with the same color at both ends.

        The name is matched against 'CYCLIC_COLORMAPS' with one '_r' reversal suffix dropped. A
        reversed cyclic colormap is still cyclic. The name of a trimmed colormap keeps its
        'trunc(...)' wrapper ('TRIMMED_COLORMAP_PATTERN'). That name is not in the list. A trimmed
        colormap is therefore never treated as cyclic. A trim that drops colors leaves the two ends
        with different colors.

        Parameters
        ==========
        cmap : matplotlib.colors.Colormap
            The colormap to classify.

        Returns
        =======
        bool
            True if the colormap is one of the cyclic ones.
        """
        name = cmap.name
        return (name[:-2] if name.endswith('_r') else name) in CYCLIC_COLORMAPS

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
                f"colormap given to '{colormap_flag}' — e.g., 'RdYlGn', 'RdYlBu_r' — is what makes "
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
        gene_aggregation: Union[str, None],
        accession_aggregation: Union[str, None],
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
        category_value_center: Union[float, None] = None,
        element_normalization: Union[str, None] = None,
        element_normalization_label: Union[str, None] = None,
        value_period: Union[float, str, None] = None
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

        'element_normalization' rescales each sample's or group's value for a map element against
        the same element's values across all samples or groups ('_resolve_element_normalization').
        This allows the enrichment or depletion of an element in a sample to be displayed rather
        than how much of it there is. Therefore, normalization requires a value column and samples.
        Normalization applies to the per-sample/per-group context alone, while the 'unified' map
        summarizes the unnormalized category values. A normalization whose neutral value is zero
        uses a centered scale on the diverging colormap, 'DEFAULT_CENTERED_COLORMAP', unless
        overridden by 'category_value_center' and 'category_colormap'. 'element_normalization_label'
        names the rescaled quantity on that scale's colorbar in place of the label the normalization
        derives for itself.

        'value_period' says that the values repeat after that period, as clock times in hours repeat
        every 24. Every reduction of the values is then circular ('CIRCULAR_AGGREGATIONS'). An
        aggregation left as None is 'circular_mean' with a period and 'sum' without one. A period
        refuses every other aggregation and summary of values. Its only normalization is
        'difference_from_circular_mean' ('CIRCULAR_ELEMENT_NORMALIZATIONS'). Presence summaries are
        unaffected. Each scale of values runs from 0 to the period ('_map_elements'). The scale of
        offsets from a normalization runs from -period / 2 to period / 2 instead. Other limits are
        refused, and so are centers. A colormap left as None is 'DEFAULT_PERIOD_COLORMAP'.
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
        # A period says that the values repeat, as clock times do every 24 hours. It is resolved
        # before the normalization, which it rules out.
        value_period_flag = f'--{element_type}-value-period'
        value_period = self._resolve_value_period(value_period, value_period_flag)

        # The normalization is settled before the centers are, since one whose neutral value is zero
        # centers the per-sample/per-group scale on zero unless asked for another center, and that
        # default has to be validated against the limits on that scale like any other center. It
        # rescales the value of each sample or group against the element's values across all samples
        # or groups, so it needs both a value column to rescale and samples.
        element_normalization_flag = f'--{element_type}-element-normalization'
        element_normalize = None
        element_normalization_label_template = None
        # Whether the center of the per-sample/per-group scale came from a normalization rather than
        # from the user, which the warning about centering one scale and not the other has to know:
        # nobody asked for this center, so there is nothing inconsistent about the other scale
        # lacking one.
        center_from_normalization = False
        if element_normalization is None and element_normalization_label is not None:
            raise ConfigError(
                f"A label was given for what '{element_normalization_flag}' would put on the "
                f"colorbar of the {element_type} layer, but no normalization was given for it to "
                f"label. That colorbar is labeled by the value column of the file at '{path}' when "
                f"nothing rescales it."
            )
        if element_normalization_label is not None and not element_normalization_label.strip():
            raise ConfigError(
                f"The label given to '{element_normalization_flag}' for the colorbar of the "
                f"{element_type} layer's rescaled scale is blank, so nothing would say what that "
                f"scale shows. Give the option a label with something in it, or give it the "
                f"normalization alone and let anvi'o compose one."
            )
        if element_normalization is not None:
            if value_column is None:
                raise ConfigError(
                    f"'{element_normalization_flag}' rescales the values of the {element_type} "
                    f"layer, but the file at '{path}' has no value column, so that layer is "
                    f"colored by presence and has no values to rescale."
                )
            if not has_sample:
                raise ConfigError(
                    f"'{element_normalization_flag}' rescales the value of each sample or group "
                    f"against the values of all of them, but the file at '{path}' has no 'sample' "
                    f"column, so there is a single set of values with nothing to compare them "
                    f"against. Add a 'sample' column to compare samples."
                )
            element_normalize, element_normalization_label_template, centered = (
                self._resolve_element_normalization(
                    element_normalization, element_normalization_flag, element_type, value_period,
                    value_period_flag
                )
            )
            # With a period, the offsets run from -period / 2 to period / 2. That scale is fixed
            # in '_map_elements', with 0 at its middle, so no center is set here. Its two ends are
            # the same point, so its default colormap is the cyclic one. A diverging colormap
            # would give the two ends opposite colors.
            if value_period is not None:
                if category_colormap is None:
                    category_colormap = DEFAULT_PERIOD_COLORMAP
            else:
                if centered and category_value_center is None:
                    category_value_center = 0.0
                    center_from_normalization = True
                if category_colormap is None and centered:
                    category_colormap = DEFAULT_CENTERED_COLORMAP

        value_center = self._resolve_value_center(
            value_center, value_center_flag, value_limits, value_limits_flag
        )
        # A center that came from a centered normalization rather than from the user is reported
        # under the normalization's own flag: nobody gave a center, so naming the center's flag
        # would name a flag that was never used.
        category_center_source_flag = (
            element_normalization_flag if center_from_normalization else category_value_center_flag
        )
        category_value_center = self._resolve_value_center(
            category_value_center, category_center_source_flag, category_value_limits,
            category_value_limits_flag, from_normalization=center_from_normalization
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
            if value_period is not None:
                raise ConfigError(
                    f"'{value_period_flag}' says how far the values of the {element_type} layer "
                    f"go before they repeat, but the file at '{path}' has no value column, so that "
                    f"layer is colored by presence and has no values."
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
        # The summary that colors the 'unified' map is the group summary with groups, and the sample
        # summary without them. A spread there ('CIRCULAR_SPREADS') gives that map a scale, a
        # colormap, and a label of its own.
        unified_summary_flag = (
            f"--{element_type}-{'group' if group_samples is not None else 'sample'}-summary"
        )
        unified_summary = group_summary if group_samples is not None else sample_summary
        unified_spread = has_sample and unified_summary in CIRCULAR_SPREADS

        # A period fixes each scale of values to run from 0 to the period ('_map_elements'). A color
        # then always means the same point in the period. Limits of exactly 0 and the period ask for
        # the same scale, so they are accepted. Any other limit, or a center, is refused. A spread
        # on the 'unified' map is fixed to run from 0 to the largest spread there can be, so no
        # limits are accepted for that map.
        if value_period is not None:
            if unified_spread and value_limits is not None:
                raise ConfigError(
                    f"'{value_limits_flag}' bounds the color scale of the 'unified' map, but "
                    f"'{unified_summary_flag} {unified_summary}' colors that map by how far the "
                    f"values are spread apart. That scale runs from 0 to "
                    f"{self._largest_angular_deviation(value_period):.3g}, the largest spread "
                    f"there can be with a period of {value_period:g}. Please drop "
                    f"'{value_limits_flag}'."
                )
            # The offsets of a normalization run from -period / 2 to period / 2 instead, with 0 at
            # the middle of the scale.
            offsets = element_normalization is not None
            category_limits = (
                (-value_period / 2, value_period / 2) if offsets else (0.0, value_period)
            )
            for setting, flag, limits in (
                (value_limits, value_limits_flag, (0.0, value_period)),
                (category_value_limits, category_value_limits_flag, category_limits)
            ):
                if setting is not None and setting != limits:
                    raise ConfigError(
                        f"'{value_period_flag}' fixes the color scale that '{flag}' bounds. It "
                        f"runs from {limits[0]:g} to {limits[1]:g}, so that a color always means "
                        f"the same point in the period. Other limits were given. Please drop "
                        f"'{flag}', or give it '{limits[0]:g} {limits[1]:g}'."
                    )
            for setting, flag in (
                (value_center, value_center_flag),
                (category_value_center, category_value_center_flag)
            ):
                if setting is not None and offsets and flag == category_value_center_flag:
                    raise ConfigError(
                        f"'{flag}' centers the color scale of the maps of the individual samples "
                        f"or groups of the {element_type} layer, but "
                        f"'{element_normalization_flag} {element_normalization}' already puts 0 "
                        f"at the middle of that scale. It runs from {category_limits[0]:g} to "
                        f"{category_limits[1]:g}. Please drop '{flag}'."
                    )
                if setting is not None:
                    raise ConfigError(
                        f"'{flag}' centers a color scale of the {element_type} layer, but "
                        f"'{value_period_flag}' fixes each such scale to run from 0 to "
                        f"{value_period:g}. Values that repeat have no middle to read the others "
                        f"against. Please use only one of the two options."
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
        # A spread is a plain number, not a point in the period. So only the summary that colors the
        # 'unified' map takes one. Within a sample, the KOs of an element must still give a point in
        # the period. With groups, the sample summary colors each group's map, and the group summary
        # then pools those values.
        for name, flag in (
            (gene_aggregation, f'--{element_type}-gene-aggregation'),
            (accession_aggregation, f'--{element_type}-accession-aggregation'),
            (sample_summary if group_samples is not None else None,
             f'--{element_type}-sample-summary')
        ):
            if name in CIRCULAR_SPREADS:
                raise ConfigError(
                    f"'{flag}' was given as '{name}', which measures how far values are spread "
                    f"apart. A spread is a plain number, not a point in the period. So it can only "
                    f"summarize the {'groups' if group_samples is not None else 'samples'} for the "
                    f"'unified' map, with '{unified_summary_flag}'."
                )
        # An aggregation left unset is 'circular_mean' for values with a period, and 'sum'
        # otherwise. A colormap of values left unset is cyclic for values with a period.
        default_aggregation = 'circular_mean' if value_period is not None else 'sum'
        # Messages that suggest a summary of values name one that the period accepts.
        suggested_summary = 'circular_mean' if value_period is not None else 'mean'
        default_colormap = DEFAULT_PERIOD_COLORMAP if value_period is not None else 'plasma_r'
        if gene_aggregation is None:
            gene_aggregation = default_aggregation
        if accession_aggregation is None:
            accession_aggregation = default_aggregation
        # Resolved here, before any values are aggregated, so that an unusable aggregation name is
        # reported as a configuration error rather than reaching pandas: the per-accession values
        # below are computed by passing the name itself to a groupby. 'aggregate' reduces a map
        # element's several accessions, which is the level the drawing colorers work at; the gene
        # aggregation is applied only where the per-accession values are built.
        aggregate = self._resolve_aggregation(
            accession_aggregation, f'--{element_type}-accession-aggregation', value_period,
            value_period_flag
        ) if value_column is not None else None
        # A circular aggregation is not a pandas name, so its function reduces a gene's rows. A
        # pandas name is passed as the name itself, which pandas applies faster.
        gene_reduction = gene_aggregation
        if value_column is not None:
            gene_aggregate = self._resolve_aggregation(
                gene_aggregation, f'--{element_type}-gene-aggregation', value_period,
                value_period_flag
            )
            if gene_aggregation in CIRCULAR_AGGREGATIONS:
                gene_reduction = gene_aggregate

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
                df, value_column, gene_reduction, path, undefined
            )
            self._warn_undefined_values(undefined, gene_aggregation, path)
            cmap = self._resolve_sequential_colormap(
                colormap if colormap is not None else default_colormap, colormap_limits,
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
                'sample_values': None,
                'aggregate': aggregate,
                'colorbar_label': value_column,
                # A normalization was refused above, there being no samples to rescale against, so
                # the one scale of this layer carries the value column's own name.
                'category_colorbar_label': value_column,
                'element_normalize': None,
                'value_limits': value_limits,
                'category_value_limits': None,
                'value_center': value_center,
                'category_value_center': None,
                'value_period': value_period
            }

        # Resolve the two summaries into the mode of each map context. Ungrouped, the per-sample
        # maps show one sample each, so no summary applies to them: they are the sample's own
        # magnitude with a value column and its presence without one.
        grouped = group_samples is not None
        sample_kind, sample_aggregate, sample_scheme = self._resolve_summary(
            sample_summary, value_column, f'--{element_type}-sample-summary', path, value_period,
            value_period_flag
        )
        if grouped:
            group_kind, group_aggregate, group_scheme = self._resolve_summary(
                group_summary, value_column, f'--{element_type}-group-summary', path,
                value_period, value_period_flag
            )
            if group_kind == 'value' and sample_kind != 'value':
                raise ConfigError(
                    f"'--{element_type}-group-summary' was given as '{group_summary}', which pools "
                    f"the values of the sample groups, but '--{element_type}-sample-summary' "
                    f"summarizes each group's samples by presence rather than by value, so the "
                    f"groups have no values to pool. Please also set "
                    f"'--{element_type}-sample-summary' to an aggregation such as "
                    f"'{suggested_summary}'."
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
            'sample_values': None,
            'element_normalize': element_normalize,
            'value_limits': value_limits,
            'category_value_limits': category_value_limits,
            'value_center': value_center,
            'category_value_center': category_value_center,
            # A center that a normalization supplied is reported under the normalization's flag.
            'category_value_center_flag': category_center_source_flag,
            'value_period': value_period
        }

        if value_column is None:
            return model

        # A normalization rescales the values that color the maps of the individual samples or
        # groups, so it has nothing to act on where those maps are not colored by value. Grouped,
        # that is what a sample summary of presence leaves them: each map shows how many of a
        # group's samples contain an element rather than how much of it there is. This is checked
        # before the limits and the center below, since a centered normalization supplies a center
        # of its own and their messages would name that center rather than what put it there.
        if element_normalization is not None and model['category_mode'] != 'quantitative':
            raise ConfigError(
                f"'{element_normalization_flag}' rescales the values that color the maps of the "
                f"individual groups, but those maps are not colored by value here: "
                f"'--{element_type}-sample-summary' summarizes each group's samples by presence "
                f"rather than by pooling their values, so each map shows how many of a group's "
                f"samples contain an element rather than how much of it there is. Set the sample "
                f"summary to an aggregation such as '{suggested_summary}' to color the group maps "
                f"by value."
            )

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
                    f"as '{suggested_summary}' to color the group maps by value, or {verb} the "
                    f"'unified' map's own scale with '{other_flag}'."
                )
        if category_colormap is not None and model['category_mode'] != 'quantitative':
            raise ConfigError(
                f"'{category_colormap_flag}' colors the scale shared by the maps of the individual "
                f"groups, but those maps are not colored along a scale of values here: "
                f"'--{element_type}-sample-summary' summarizes each group's samples by presence "
                f"rather than by pooling their values, and presence there is colored by "
                f"'--group-colormap'. Set the sample summary to an aggregation such as "
                f"'{suggested_summary}' to color the group maps by value, or color the 'unified' "
                f"map's own scale with '{colormap_flag}'."
            )

        # The 'unified' scale and the scale of the individual maps take separate limits. A summary
        # can span much less than the values it summarizes. For example, a mean across samples spans
        # less than the samples do. So limits given for one scale do not apply to the other. Where
        # both scales color by value and only one has limits, the user may not have meant to leave
        # the other without limits. This warns about it. It is not refused, since limiting only one
        # scale can be intended. With a period, both scales have fixed limits, so there is nothing
        # to warn about.
        if (
            model['unified_mode'] == 'quantitative' and model['category_mode'] == 'quantitative'
            and (value_limits is None) != (category_value_limits is None)
            and value_period is None
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

        # Centering only one of the two scales also warns, for the same reason as the limits above.
        # There is a second reason. Without '--*-category-colormap', both scales use the same
        # colormap. The middle color then means the center on one map, and the midpoint of that
        # map's values on the other. A center that a normalization supplied does not warn. The user
        # did not ask for it, and a normalization colors its scale with a diverging colormap that the
        # summary does not share.
        if (
            model['unified_mode'] == 'quantitative' and model['category_mode'] == 'quantitative'
            and (value_center is None) != (category_value_center is None)
            and not center_from_normalization
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
                f"{given_center:g}, but '{other_flag}' was not given. So the scale of "
                f"{other_subject} is not centered. Unless '{category_colormap_flag}' gives the two "
                f"scales different colormaps, they use the same colormap. Its middle color then "
                f"means {given_center:g} on one map, and a different value on the other. The two "
                f"scales are centered separately, because a summary can span much less than the "
                f"values it summarizes. So this may be what you intend. If it is not, center the "
                f"other scale with '{other_flag}'."
            )

        # Per-sample values: each sample's rows reduced to one value per accession by the
        # within-sample aggregation. These are computed even when no context is colored by value, so
        # that a value column is validated whenever the file has one.
        sample_values: Dict[str, Dict[str, float]] = {s: {} for s in all_sample_names}
        for sample_name, sample_rows in df.groupby('__sample'):
            sample_values[sample_name] = self._aggregate_accession_quantities(
                sample_rows, value_column, gene_reduction, path, undefined
            )
        self._warn_undefined_values(undefined, gene_aggregation, path)

        if 'quantitative' not in (model['unified_mode'], model['category_mode']):
            # Both summaries color presence, which only happens with groups (the per-sample maps of
            # an ungrouped run always show a value column's magnitude), so nothing is colored by
            # value.
            self.run.warning(
                f"The layer from the text file at '{path}' has a value column, '{value_column}', "
                f"but nothing on the maps is colored by it: with sample groups, the per-group maps "
                f"are colored by '--{element_type}-sample-summary' and the 'unified' map by "
                f"'--{element_type}-group-summary', and both of these summarize presence rather "
                f"than value. Set '--{element_type}-sample-summary' to an aggregation, such as "
                f"'{suggested_summary}', to color the group maps by value."
            )
            return model

        # A summary pools the values of a map element, not of an accession. An element's value in a
        # sample is the aggregate of its accessions in that sample. It is what that sample's own map
        # draws. One accession can belong to different elements on different maps. So these values
        # are known only once a map is read. The summaries are therefore taken map by map
        # ('_summarize_entry_values'). Here the layer keeps what they need. Ungrouped, the sample
        # summary colors the 'unified' map. Grouped, it colors each group's map from that group's
        # samples. The group summary then colors the 'unified' map from the groups.
        model['sample_values'] = sample_values
        model['group_samples'] = group_samples
        model['sample_summary'] = sample_summary if sample_kind == 'value' else None
        model['sample_aggregate'] = sample_aggregate
        model['group_summary'] = group_summary if grouped and group_kind == 'value' else None
        model['group_aggregate'] = group_aggregate if grouped else None

        # A spread on the 'unified' map does not repeat. Its scale takes a sequential colormap. A
        # cyclic one would give the smallest and the largest spread the same color.
        model['unified_spread'] = unified_spread
        unified_default_colormap = 'plasma_r' if unified_spread else default_colormap
        model['cmap'] = self._resolve_sequential_colormap(
            colormap if colormap is not None else unified_default_colormap, colormap_limits,
            subject=f'{element_type}s'
        )
        if unified_spread and self._is_cyclic_colormap(model['cmap']):
            raise ConfigError(
                f"'{colormap_flag}' was given as '{colormap}', a cyclic colormap, but "
                f"'{unified_summary_flag} {unified_summary}' colors the 'unified' map by how far "
                f"the values are spread apart. A spread does not repeat. A cyclic colormap would "
                f"give the smallest and the largest spread the same color. Please give a "
                f"sequential colormap, such as 'plasma_r'."
            )
        # The maps of the individual samples or groups take a colormap of their own where one was
        # given, so that a 'unified' map showing a summary of another kind than the values behind it
        # — their spread rather than their magnitude — is not read as more of the same quantity.
        # Left unset, the layer's one colormap serves both contexts, and is resolved once so that a
        # warning about it is given once. The exception is a spread on the 'unified' map. The maps
        # of the individual samples or groups then show points in the period. Left unset, their
        # colormap is the cyclic default.
        category_subject = (
            f"{element_type}s on the maps of the individual "
            f"{'groups' if grouped else 'samples'}"
        )
        if category_colormap is not None:
            model['category_cmap'] = self._resolve_sequential_colormap(
                category_colormap, category_colormap_limits, subject=category_subject
            )
        elif unified_spread:
            model['category_cmap'] = self._resolve_sequential_colormap(
                default_colormap, None, subject=category_subject
            )
        else:
            model['category_cmap'] = model['cmap']
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
        # A spread on the 'unified' map shows a quantity of its own. Its colorbar says what it is.
        if unified_spread:
            model['unified_colorbar_label'] = f'{value_column} angular deviation'
        # The 'unified' map is derived from the values themselves, so only the scale that the
        # per-sample or per-group maps share is relabeled by a normalization. Every normalization
        # composes its label from the value column's name: one anvi'o knows names the quantity it
        # makes ('{value} - mean'), and any other is named by the pandas method it is
        # ('{value} (abs)').
        model['category_colorbar_label'] = (
            value_column if element_normalization_label_template is None
            else element_normalization_label_template.format(value=value_column)
        )
        if element_normalization_label is not None:
            model['category_colorbar_label'] = element_normalization_label
        return model

    def map_kegg_pathways_txt(
        self,
        output_dir: str,
        reaction_txt: str = None,
        compound_txt: str = None,
        reaction_gene_aggregation: str = None,
        reaction_accession_aggregation: str = None,
        compound_accession_aggregation: str = None,
        reaction_value_period: Union[float, str, None] = None,
        compound_value_period: Union[float, str, None] = None,
        reaction_sample_summary: str = None,
        compound_sample_summary: str = None,
        reaction_group_summary: str = None,
        compound_group_summary: str = None,
        groups_txt: str = None,
        group_threshold: float = None,
        pathway_numbers: Iterable[str] = None,
        draw_unified_maps: bool = True,
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
        reaction_element_normalization: Union[str, None] = None,
        reaction_element_normalization_label: Union[str, None] = None,
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
        compound_element_normalization: Union[str, None] = None,
        compound_element_normalization_label: Union[str, None] = None,
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

        reaction_gene_aggregation : str, None
            How to reduce the values of the genes annotated with one accession to that accession's
            value, for a reaction file carrying a 'gene_id' column. Any 'AGGREGATION_FUNCTIONS' name,
            or any other pandas aggregation reducing values to one number (see
            '_resolve_aggregation'). The default of None is 'sum', or 'circular_mean' when
            'reaction_value_period' is given.

        reaction_accession_aggregation : str, None
            How to reduce the values of the several accessions a map's reaction element stands to
            that element's value. The default is the same as for 'reaction_gene_aggregation'.

        compound_accession_aggregation : str, None
            The same reduction for the compound layer: the several compounds that one map circle
            stands for. A compound file has no genes, and repeated rows are refused
            ('_read_element_txt'), so it has no reduction below this one. The default of None is
            'sum', or 'circular_mean' when 'compound_value_period' is given.

        reaction_value_period : Union[float, str, None], None
            The period after which the values of the reaction layer repeat, such as 24 for clock
            times in hours. Every reduction of the values is then circular, and 'circular_mean' is
            the only aggregation and value summary accepted ('CIRCULAR_AGGREGATIONS'). The only
            normalization accepted is 'difference_from_circular_mean'. Each scale of values runs
            from 0 to the period. The scale of its offsets runs from -period / 2 to period / 2.
            Other limits are refused, and so are centers. Without a colormap, the scales take
            'DEFAULT_PERIOD_COLORMAP'. None means that the values do not repeat.

        compound_value_period : Union[float, str, None], None
            The same period for the values of the compound layer.

        reaction_sample_summary : str, None
            How the reaction layer summarizes a set of samples: by presence ('count'/'membership')
            or by pooling their values with an aggregation. An aggregation pools the values of each
            map element. These are the values that the maps of the individual samples draw. This
            colors the 'unified' map without groups and each per-group map with groups. The default
            of None summarizes presence, by membership for 3 or fewer samples and by count above
            that.

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
            either end of which can be None to leave it wherever the values put it. A limit sets
            that end of the scale. The colorbar labels the limit. It marks the limit '<=' or '>='
            where values lie past it ('_make_quantitative_norm').

        reaction_category_value_limits : Tuple[Union[float, None], Union[float, None]], None
            The same limits for the scale shared by the reaction layer's per-sample or per-group
            maps. It is bounded separately from the 'unified' map's scale, because a summary can
            span much less than the values it summarizes.

        compound_value_limits : Tuple[Union[float, None], Union[float, None]], None
            The same limits for the compound layer's 'unified' map scale.

        compound_category_value_limits : Tuple[Union[float, None], Union[float, None]], None
            The same limits for the compound layer's per-sample or per-group scale.

        reaction_value_center : Union[float, None], None
            The value put at the middle of the reaction layer's color scale on the 'unified' map.
            The scale then runs the same distance either side of it, so that the middle color of a
            diverging colormap stands for this value however lopsided the values around it are
            ('_make_quantitative_norm'). The default of None leaves the scale where its values and
            any limits set it.

        reaction_category_value_center : Union[float, None], None
            The same center for the scale shared by the reaction layer's per-sample or per-group
            maps. It is centered separately from the 'unified' map's scale, because a summary can
            span much less than the values it summarizes.

        compound_value_center : Union[float, None], None
            The same center for the compound layer's 'unified' map scale.

        compound_category_value_center : Union[float, None], None
            The same center for the compound layer's per-sample or per-group scale.

        reaction_element_normalization : Union[str, None], None
            How to rescale each sample's or group's value for a reaction element against the same
            element's values across all samples or groups. This allows the enrichment or depletion
            of an element in a sample to be displayed rather than how much of it there is. The
            argument must be one of 'ELEMENT_NORMALIZATIONS' or the name of a pandas Series method
            that transforms each value into another value ('_resolve_element_normalization'). It
            rescales the maps of the individual samples or groups, while the 'unified' map
            summarizes the unnormalized values of the categories. The default of None draws each
            sample's own values.

        reaction_element_normalization_label : Union[str, None], None
            What the colorbar of the rescaled reaction scale is labeled, in place of the label the
            normalization derives from the value column's name.

        compound_element_normalization : Union[str, None], None
            The same normalization for the compound layer.

        compound_element_normalization_label : Union[str, None], None
            The same colorbar label for the compound layer's rescaled scale.

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
            The record returned by '_map_elements'. All three keys are always present. The
            dictionary under 'unified' is empty when 'draw_unified_maps' is False, since no
            'unified' map is drawn.
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
                'category_value_center': reaction_category_value_center,
                'element_normalization': reaction_element_normalization,
                'element_normalization_label': reaction_element_normalization_label,
                'value_period': reaction_value_period
            })
        if compound_txt is not None:
            raw_layers.append({
                'name': 'compounds',
                'path': compound_txt,
                'data': self._read_element_txt(compound_txt, 'compound'),
                # A compound file has one row per accession per sample, so the gene level is a
                # reduction over a single value. Left unset, it takes the layer's default. A period
                # makes that default circular, and a period would refuse a 'sum' here.
                'gene_aggregation': None,
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
                'category_value_center': compound_category_value_center,
                'element_normalization': compound_element_normalization,
                'element_normalization_label': compound_element_normalization_label,
                'value_period': compound_value_period
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
            draw_unified_maps=draw_unified_maps,
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
        draw_unified_maps: bool = True,
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
        'draw_unified_maps', 'draw_individual_files', 'draw_grid', and the group colormap options)
        mirror the categorical engine; see the CLI help and '_map_element_membership'.
        'reaction_category_colors' is the path to a kegg-category-colors-txt file giving a color per
        category — per source, or per group when the sources are grouped — which colors the layer by
        membership in place of a colormap, and colors each category's own map. 'group_colormap' may
        be 'GROUP_COLORMAP_FROM_CATEGORY' to color each group's own maps by a ramp running to that
        group's own color instead of by a named colormap, in which case 'group_colormap_limits' is
        how far from white that ramp runs ('_group_map_colors'). 'group_colormap_scheme' draws the
        count scale of every group's maps in discrete bands ('by_count') or as a gradient
        ('by_count_continuous'), and by default keeps the bands while the colors can be told apart.

        Returns
        =======
        Dict[Literal['unified', 'individual', 'grid'], Dict]
            The record returned by '_map_element_membership': 'unified' maps show all contigs
            databases or groups, 'individual' maps show single databases or groups, and 'grid'
            images show both. See '_map_element_membership' for the nested structure. All three
            keys are always present. The dictionary under 'unified' is empty when
            'draw_unified_maps' is False, since no 'unified' map is drawn.
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

        # Load groups. A threshold decides when a group counts as containing an element. Only the
        # 'unified' map uses that. A group's own map counts the group's own databases whatever the
        # threshold is. A run not drawing the 'unified' map therefore refuses a threshold rather
        # than accepting one that decides nothing. Such a run may give 'groups_txt' on its own.
        if group_threshold is not None and not draw_unified_maps:
            raise ConfigError(
                "A 'group_threshold' was given, which sets the proportion of a group's contigs "
                "databases that must contain an element for the group to count as containing it. "
                "Only the 'unified' map uses that threshold. A group's own map counts the group's "
                "own databases whatever the threshold is. This run does not draw the 'unified' "
                "map. Please pass 'group_threshold' or set 'draw_unified_maps' to False, not both."
            )
        if (groups_txt is None and group_threshold is not None) or (
            groups_txt is not None and group_threshold is None and draw_unified_maps
        ):
            raise ConfigError(
                "To group contigs databases, arguments to both 'groups_txt' and 'group_threshold' "
                "must be provided."
            )

        source_group: Dict[str, str] = None
        group_sources: Dict[str, List[str]] = None
        if groups_txt is not None:
            if group_threshold is not None and not 0 <= group_threshold <= 1:
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
            draw_unified_maps=draw_unified_maps,
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
        draw_unified_maps: bool = True,
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
        'draw_unified_maps', 'draw_individual_files', 'draw_grid', and the group colormap options)
        mirror the categorical engine; see the CLI help and '_map_element_membership'.
        'reaction_category_colors' is the path to a kegg-category-colors-txt file giving a color per
        category — per source, or per group when the sources are grouped — which colors the layer by
        membership in place of a colormap, and colors each category's own map. 'group_colormap' may
        be 'GROUP_COLORMAP_FROM_CATEGORY' to color each group's own maps by a ramp running to that
        group's own color instead of by a named colormap, in which case 'group_colormap_limits' is
        how far from white that ramp runs ('_group_map_colors'). 'group_colormap_scheme' draws the
        count scale of every group's maps in discrete bands ('by_count') or as a gradient
        ('by_count_continuous'), and by default keeps the bands while the colors can be told apart.

        Returns
        =======
        Dict[Literal['unified', 'individual', 'grid'], Dict]
            The record returned by '_map_element_membership': 'unified' maps show all genomes or
            groups, 'individual' maps show single genomes or groups, and 'grid' images show both.
            See '_map_element_membership' for the nested structure. All three keys are always
            present. The dictionary under 'unified' is empty when 'draw_unified_maps' is False,
            since no 'unified' map is drawn.
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

        # Load groups. A threshold decides when a group counts as containing an element. Only the
        # 'unified' map uses that. A group's own map counts the group's own genomes whatever the
        # threshold is. A run not drawing the 'unified' map therefore refuses a threshold rather
        # than accepting one that decides nothing. Such a run may give 'groups_txt' on its own.
        if group_threshold is not None and not draw_unified_maps:
            raise ConfigError(
                "A 'group_threshold' was given, which sets the proportion of a group's genomes "
                "that must contain an element for the group to count as containing it. Only the "
                "'unified' map uses that threshold. A group's own map counts the group's own "
                "genomes whatever the threshold is. This run does not draw the 'unified' map. "
                "Please pass 'group_threshold' or set 'draw_unified_maps' to False, not both."
            )
        if (groups_txt is None and group_threshold is not None) or (
            groups_txt is not None and group_threshold is None and draw_unified_maps
        ):
            raise ConfigError(
                "To group genomes, arguments to both 'groups_txt' and 'group_threshold' must be "
                "provided."
            )

        source_group: Dict[str, str] = None
        group_sources: Dict[str, List[str]] = None
        if groups_txt is not None:
            if group_threshold is not None and not 0 <= group_threshold <= 1:
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
            draw_unified_maps=draw_unified_maps,
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
        pathway_numbers = self._find_maps([output_dir], patterns=pathway_numbers)

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

        # The color is checked before any map is drawn, as '_map_elements' checks its colors. It
        # colors the reactions, and, on global and overview maps, the compounds they touch.
        if color_hexcode != 'original':
            check_layers = [
                {**layer, 'unified_mode': 'single', 'category_mode': 'single'}
                for layer in presence_layers
            ]
            self._check_reserved_colors(check_layers, pathway_numbers, drawn_categories=[])
            self._check_derived_compound_colors(
                check_layers, pathway_numbers,
                lambda: [(None, True, self._element_presence_specs(presence_layers))]
            )

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
        categories: List[str],
        category_noun: str,
        draw_unified_maps: bool = True
    ) -> None:
        """
        Reduce each layer's category colors to those the run will actually use.

        A category with no color of its own cannot be colored at all, so it is an error, reported
        for every such category at once. A color given for a name that is not a category, or for a
        combination naming one, describes something this run does not draw, so it is reported and
        ignored, as a groups file's extra items are ('_relate_samples_to_groups'): one file can then
        cover a set of samples that different runs draw different subsets of. Ignoring such a color
        means dropping it from the layer rather than merely reporting it, which is what makes the
        colors that remain the colors of the run — everything downstream reads them as such. A
        color naming a combination of categories is dropped the same way where the 'unified' map is
        not drawn: only that map's scale gives a combination a color of its own.

        Parameters
        ==========
        layers : List[dict]
            The layer models, whose 'category_colors'/'category_combo_colors' are checked against
            the run and reduced to it in place. A layer without colors of its own is skipped.

        categories : List[str]
            The names of every category of the run.

        category_noun : str
            What one category is called, for error messages, e.g. 'sample' or 'group'.

        draw_unified_maps : bool, True
            If True, the 'unified' map is drawn. If False, it is not. Only that map's color scale
            blends a combination of categories. When it is not drawn, a color given for a
            combination is reported and dropped.
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

            # Only the 'unified' map's scale blends a combination of categories, so a color
            # overriding one of those blends has nothing to act on where that map is not drawn.
            # The rows go with a word rather than refusing the run: the colors given per category
            # still draw the individual maps, so one file can serve a run that draws the 'unified'
            # map and a run that does not.
            combo_colors = layer['category_combo_colors']
            if not draw_unified_maps and combo_colors:
                separator = f'{CATEGORY_COMBO_SEPARATOR} '
                message = ', '.join(
                    f"'{name}'" for name in
                    sorted(separator.join(combo) for combo in combo_colors)
                )
                self.run.warning(
                    f"The colors file given to '{flag}' recolors the following combinations of "
                    f"{category_noun}s: {message}. A combination takes a color only on the "
                    f"'unified' map, where an element is colored by exactly which {category_noun}s "
                    f"contain it. That map is not drawn here, so these rows are ignored. The "
                    f"colors given to single {category_noun}s are unaffected. Each one still draws "
                    f"that {category_noun}'s own map."
                )
                combo_colors.clear()

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
        Look up a colormap by name.

        The name is looked for first among the colormaps anvi'o defines ('ANVIO_COLORMAPS') and then
        among Matplotlib's. An '_r' suffix reverses a colormap of either kind.

        Parameters
        ==========
        colormap : str
            The name of an anvi'o or Matplotlib colormap.

        Returns
        =======
        matplotlib.colors.Colormap
            The named colormap.
        """
        reverse = colormap.endswith('_r')
        anvio_name = colormap[:-2] if reverse else colormap
        if anvio_name in ANVIO_COLORMAPS:
            cmap = mcolors.LinearSegmentedColormap.from_list(
                anvio_name, ANVIO_COLORMAPS[anvio_name]
            )
            return cmap.reversed() if reverse else cmap
        try:
            return colormaps[colormap]
        except KeyError:
            self.progress.end()
            anvio_names = ', '.join(f"'{name}'" for name in ANVIO_COLORMAPS)
            raise ConfigError(
                f"'{colormap}' is not the name of a colormap. Matplotlib's colormaps include "
                f"'plasma', 'viridis' and 'tab10'. Anvi'o also defines colormaps of its own: "
                f"{anvio_names}. Adding '_r' to any of these names reverses the colormap. For "
                f"example, 'plasma_r' is 'plasma' reversed. Matplotlib's colormaps are listed at "
                f"https://matplotlib.org/stable/users/explain/colors/colormaps.html"
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
                f"A colormap must be given as a name or as a Colormap object. This one is neither. "
                f"It was given as {colormap}."
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

        Each of the 'limits' sets that end of the scale, whether or not any value reaches it. An end
        left as None is set by the values. A limit keeps its end the same from run to run. With both
        ends set, a color means the same value on maps drawn from different data. Every value past a
        limit takes the color of that end of the scale, which 'clip=True' on the normalization
        arranges. The caller then labels that end of the colorbar '<=' or '>=', so that its color
        reads as "this value or past it" rather than as an exact value.

        'center' then widens whichever side of the range falls short, until the range runs the same
        distance either side of it. That is what puts the centered value at the middle of the
        colormap — the neutral color of a diverging one — however lopsided the values around it are,
        and it keeps the scale linear, so that the same distance in color goes on meaning the same
        distance in value on either side. Widening can only move an end that the values set. Moving
        a limit would undo what the limit was asked to do, so it is refused rather than carried out
        quietly. Two limits leave no end to move. '_resolve_value_center' has already refused a
        center that is not their midpoint, so such a scale is used as the limits set it.

        Parameters
        ==========
        values : List[float]
            All per-reaction values that will be colored on the maps sharing this normalization.

        limits : Union[Tuple[Union[float, None], Union[float, None]], None], None
            The (minimum, maximum) the scale may span, as '_resolve_value_limits' returns them.
            Either end can be None to leave it wherever the values put it. None means the scale is
            not limited at all.

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
            Five values:
            - The normalization.
            - vmin, the bottom of the scale.
            - vmax, the top of the scale.
            - Whether some value lies below the minimum limit. The colorbar then marks vmin '<='.
            - Whether some value lies above the maximum limit. The colorbar then marks vmax '>='.

            With no values, the first three are None and both flags are False.

            The normalization is also None where the scale is a single value. vmin then equals vmax.
            This happens in three cases:
            - No limit or center is given, and every value is the same.
            - The one limit given lies on the far end of the values.
            - A center is given, and every value equals it.
            Callers then give every element the top color of the colormap. Where the scale is
            centered, they give it the middle color instead.
        """
        if not values:
            return None, None, None, False, False
        data_min = min(values)
        data_max = max(values)
        limit_min, limit_max = (None, None) if limits is None else limits
        vmin = data_min if limit_min is None else limit_min
        vmax = data_max if limit_max is None else limit_max
        if (limit_min is not None and limit_min > data_max) or (
            limit_max is not None and limit_max < data_min
        ):
            # Every value lies past one limit, so every element would take the color at that end of
            # the scale. Only one of the two limits can be at fault: a pair whose minimum is below
            # its maximum ('_resolve_value_limits') can miss the values only by sitting wholly to
            # one side of them, so the end that does is the one named.
            if limit_min is not None and limit_min > data_max:
                end_noun, limit, side, remedy = 'minimum', limit_min, 'below', f'below {data_max:g}'
            else:
                end_noun, limit, side, remedy = 'maximum', limit_max, 'above', f'above {data_min:g}'
            raise ConfigError(
                f"'{flag}' was given a {end_noun} of {limit:g}, but every value coloring these "
                f"maps falls {side} it -- they run from {data_min:g} to {data_max:g}. Every "
                f"element would therefore take the color at that end of the scale. The map could "
                f"not tell one from another. Please set the {end_noun} {remedy}, or drop the "
                f"limits."
            )
        # Two limits set both ends. The center is their midpoint ('_resolve_value_center'), so there
        # is nothing for centering to do.
        if center is not None and (limit_min is None or limit_max is None):
            # The farther of the two ends from the center sets how far the scale reaches on both
            # sides. Only one of the two distances can be negative, which happens when the whole
            # range lies to one side of the center, so the larger of them is never negative.
            low_distance = center - vmin
            high_distance = vmax - center
            half_range = max(low_distance, high_distance)
            centered_min, centered_max = center - half_range, center + half_range

            def moves_end(distance: float, far_distance: float) -> bool:
                # Centering leaves the farther end exactly where it is and stretches only the nearer
                # one out to match. An end is moved if it is the nearer one by more than rounding
                # error. A limit typed at the centered end can differ from it in the last digit.
                reach = max(distance, far_distance)
                return reach - distance > 1e-9 * reach

            for (
                limit, distance, far_distance, centered_end, end_noun, other_end_noun, comparison,
                loosen_verb
            ) in (
                (
                    limit_min, low_distance, high_distance, centered_min, 'minimum', 'maximum',
                    'below', 'lower'
                ),
                (
                    limit_max, high_distance, low_distance, centered_max, 'maximum', 'minimum',
                    'above', 'raise'
                )
            ):
                if limit is None or not moves_end(distance, far_distance):
                    continue
                # The changes that keep a limit are worth working out for the user rather than
                # describing. One is where the centered scale ends, which is as far as this limit
                # can be moved. The other is the mirror of the limit about the center. That is the
                # other limit that would set a centered scale at both ends. It is offered only where
                # some value lies short of it, since otherwise every element would take the color
                # at that end. Each suggested number is printed in enough digits that typing it back
                # in works. The given numbers are printed exactly. No two numbers then print the
                # same without being the same.
                def short_of_data(value: float) -> bool:
                    return value > data_min if comparison == 'below' else value < data_max

                end_text = self._format_accepted(
                    centered_end,
                    lambda value: not moves_end(
                        center - value if comparison == 'below' else value - center, far_distance
                    )
                )
                mirrored = 2 * center - limit
                mirror_fits = short_of_data(mirrored)
                mirror_text = self._format_accepted(
                    mirrored,
                    lambda value: short_of_data(value)
                    and abs(center - (limit + value) / 2) <= 1e-9 * abs(value - limit)
                )
                limit_text, center_text = (
                    self._format_accepted(given, lambda value, given=given: value == given)
                    for given in (limit, center)
                )
                changes = [
                    f"{loosen_verb.capitalize()} the {end_noun} to {end_text}, which is where a "
                    f"centered scale ends here."
                ]
                if mirror_fits:
                    changes.append(
                        f"Or give the scale a {other_end_noun} of {mirror_text} as well, the same "
                        f"distance from {center_text} as the {end_noun}."
                    )
                changes.append(f"Or drop the {end_noun} and let the values set that end.")
                other_end = vmax if comparison == 'below' else vmin
                raise ConfigError(
                    f"'{flag}' sets the {end_noun} of the color scale of {subject} at "
                    f"{limit_text}, while '{center_flag}' centers that same scale on "
                    f"{center_text}. The two conflict. The values put the other end of the scale "
                    f"at {other_end:g}. A scale centered on {center_text} would then have its "
                    f"{end_noun} at {end_text}, which lies {comparison} the limit. "
                    f"{'Three' if mirror_fits else 'Two'} changes would settle it. "
                    f"{' '.join(changes)}"
                )
            # A limit stays exactly where it was given. The centered end can differ from it by
            # rounding error.
            vmin = centered_min if limit_min is None else limit_min
            vmax = centered_max if limit_max is None else limit_max
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
                    f"center, drop the option. The values and any limits then set the scale."
                )
        # Only a limit can be marked, and only where values do lie past it. Centering can leave an
        # end the values set a rounding error short of them. Each limit is therefore compared with
        # the values, rather than the end of the range.
        clamped_low = limit_min is not None and data_min < limit_min
        clamped_high = limit_max is not None and data_max > limit_max
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

    @staticmethod
    def _summarize_entry_values(
        entry: kgml.Entry,
        layer: dict
    ) -> Tuple[Dict[str, float], Union[float, None]]:
        """
        Find an Entry's value in each sample or group of a layer, and on the 'unified' map.

        An element's value in a sample is the aggregate of its accessions in that sample
        ('_reduce_entry_value'). This is the value the sample's own map draws. A summary then pools
        these element values. So the 'unified' map summarizes what the maps of the individual
        samples or groups draw. Grouped, each group's value is the sample summary of its samples'
        element values. The 'unified' value is then the group summary of the groups' values.

        A summary is over the samples or groups in which the element has a value. A sample with no
        row for any of the element's accessions does not count as a zero. A summary can be
        undefined, as the standard deviation of a single value is. The element then has no value in
        that group, or on the 'unified' map. The tuple of its accessions is then added to the
        layer's '_undefined_category' or '_undefined_unified' set. '_warn_undefined_summaries'
        reports these sets once for the layer.

        Parameters
        ==========
        entry : kgml.Entry
            An Entry (ortholog or compound) whose accessions are read by '_get_entry_kegg_ids'.

        layer : dict
            A 'quantitative' layer model from '_build_txt_model'. A layer without samples has
            'sample_values' of None. Its one value per element comes from 'unified_values'. The
            range-finding pass of '_map_elements' gives the layer an '_element_cache' and the two
            undefined sets. Results are kept in the cache, keyed by an Entry's accessions. These are
            all an answer depends on. Entries on different maps, with different IDs, can stand for
            the same accessions.

        Returns
        =======
        Tuple[Dict[str, float], Union[float, None]]
            The element's value in each sample (ungrouped) or group (grouped) that has one, and its
            value on the 'unified' map, or None where it has none. The first is empty for a layer
            without samples.
        """
        use_reaction_attribute = layer['use_reaction_attribute']
        # One call works out the element's value in every sample or group and on the 'unified' map.
        # The range-finding pass makes the first call. The result is then reused by the 'unified'
        # map and the map of each sample or group. The result is keyed by the accessions the Entry
        # stands for, not by its ID, which is unique to one map. An element's value depends only on
        # its accessions. So the result is also reused for an Entry on another KEGG map, or a second
        # Entry on one KEGG map, that stands for the same accessions.
        key = tuple(Mapper._get_entry_kegg_ids(entry, use_reaction_attribute))
        cache = layer['_element_cache']
        if key in cache:
            return cache[key]

        # A layer without samples has one value per accession. The element's value is the aggregate
        # of its accessions, and every map shows it.
        aggregate = layer['aggregate']
        if layer['sample_values'] is None:
            result = {}, Mapper._reduce_entry_value(
                entry, layer['unified_values'], aggregate,
                use_reaction_attribute=use_reaction_attribute
            )
            cache[key] = result
            return result

        # The element's value in each sample is the aggregate of the accessions the sample has
        # values for. A sample with none of them has no value. It is left out of every summary
        # below.
        sample_element_values: Dict[str, float] = {}
        for sample, values in layer['sample_values'].items():
            value = Mapper._reduce_entry_value(
                entry, values, aggregate, use_reaction_attribute=use_reaction_attribute
            )
            if value is not None:
                sample_element_values[sample] = value

        # A summary can be undefined for the values it is given, as a standard deviation is for a
        # single value. The element then has no value at that level. Its key is recorded for
        # '_warn_undefined_summaries'.
        def _summarize(values: List[float], summary: Callable, undefined_key: str):
            value = summary(values)
            if np.isfinite(value):
                return value
            layer[undefined_key].add(key)
            return None

        if layer['group_samples'] is None:
            # Ungrouped, the map of each sample shows the sample's own value. The sample summary
            # pools these values on the 'unified' map.
            category_values = sample_element_values
            unified_aggregate = layer['sample_aggregate']
        else:
            # Grouped, the map of each group shows the sample summary of its own samples' values. A
            # group none of whose samples has a value has no value either. A sample summary of
            # presence colors the group maps by count instead ('_group_map_colors'), so the groups
            # get no values here. The group summary pools the groups' values on the 'unified' map.
            category_values = {}
            if layer['sample_aggregate'] is not None:
                for group, samples in layer['group_samples'].items():
                    group_values = [
                        sample_element_values[sample] for sample in samples
                        if sample in sample_element_values
                    ]
                    if not group_values:
                        continue
                    value = _summarize(
                        group_values, layer['sample_aggregate'], '_undefined_category'
                    )
                    if value is not None:
                        category_values[group] = value
            unified_aggregate = layer['group_aggregate']

        # The 'unified' value pools the samples or groups that have a value. The summary is None
        # where the 'unified' map shows presence, and then there is no value to work out.
        unified_value = None
        if unified_aggregate is not None and category_values:
            unified_value = _summarize(
                list(category_values.values()), unified_aggregate, '_undefined_unified'
            )
        result = category_values, unified_value
        cache[key] = result
        return result

    @staticmethod
    def _normalize_entry_value(
        entry: kgml.Entry,
        layer: dict,
        categories: Iterable[str],
        normalize,
        cache: Union[Dict[Tuple[str, ...], Dict[str, float]], None] = None
    ) -> Dict[str, float]:
        """
        Rescale an Entry's value in each category against its values across all of them.

        An element's value in each sample or group comes from '_summarize_entry_values'. It is built
        from the aggregate of the element's accessions. Which accessions an element stands for is
        known only from the KEGG map itself. So a ratio against the element's mean, for instance, is
        a ratio of aggregates.

        The values the normalization sees are those of the categories in which the element has a
        value at all, so a sample with no row for any of the element's accessions is not counted as
        a zero — the same rule a sample or group summary follows.

        Parameters
        ==========
        entry : kgml.Entry
            An Entry (ortholog or compound) whose accessions are read by '_get_entry_kegg_ids'.

        layer : dict
            A 'quantitative' layer model with samples. Its element values are the ones rescaled.

        categories : Iterable[str]
            The categories to normalize across, in the order they are colored in.

        normalize : callable
            Rescales an element's values across the categories that have one
            ('_resolve_element_normalization').

        cache : Union[Dict[Tuple[str, ...], Dict[str, float]], None], None
            Results already worked out for this layer, keyed by an Entry's accessions, which is all
            an answer depends on once the layer's values, aggregation and normalization are fixed. A
            normalized value needs every category's value for the element, so the map of each
            category would otherwise work out every other category's values again; the
            range-finding pass fills this in beforehand, leaving the drawing passes nothing to
            recompute. Pass None to rescale every element afresh. The element values themselves stay
            cached on the layer ('_summarize_entry_values').

        Returns
        =======
        Dict[str, float]
            The normalized value per category, holding only the categories the normalization
            defined one for. Empty if the element has no value in any category.
        """
        key = None
        if cache is not None:
            key = tuple(Mapper._get_entry_kegg_ids(entry, layer['use_reaction_attribute']))
            if key in cache:
                return cache[key]
        element_values = Mapper._summarize_entry_values(entry, layer)[0]
        values: Dict[str, float] = {
            category: element_values[category] for category in categories
            if category in element_values
        }
        if not values:
            normalized_values: Dict[str, float] = {}
        else:
            # A normalization can be undefined for the values of this element even where it is
            # defined elsewhere on the map (a z-score needs more than one value, and a ratio needs a
            # positive reference), in which case the element has no value and is left uncolored.
            normalized = normalize(list(values.values()))
            normalized_values = {
                category: float(value)
                for category, value in zip(values, normalized) if np.isfinite(value)
            }
        if cache is not None:
            cache[key] = normalized_values
        return normalized_values

    def _draw_quantitative_colorbar(
        self,
        cmap: mcolors.Colormap,
        vmin: float,
        vmax: float,
        out_path: str,
        label: str,
        integer_ticks: bool = False,
        limited_low: bool = False,
        limited_high: bool = False,
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

        limited_low : bool, False
            If True, a value limit set the bottom of the range, so 'vmin' is labeled.

        limited_high : bool, False
            The same for the top of the range and 'vmax'.

        clamped_low : bool, False
            If True, values lie below the limit at the bottom of the range. They are drawn in the
            color at that end, and its label is marked accordingly.

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
                limited_low=limited_low, limited_high=limited_high, clamped_low=clamped_low,
                clamped_high=clamped_high, center=center
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
        draw_unified_maps: bool = True,
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
            One or two layer models. They are ordered from the bottom layer to the top one, so a
            reaction layer comes before a compound layer. The keys a model needs depend on how it
            colors its elements in each context.

            Every layer has these keys:
            - 'name': the stem of its colorbar file names.
            - 'element_type': 'reaction' or 'compound'.
            - 'use_reaction_attribute': True if its accessions are KEGG reaction IDs rather than KO
              IDs.

            A text file layer gives its mode for the 'unified' map as 'unified_mode', and its mode
            for the individual maps as 'category_mode'. A database or pangenome layer gives one
            'mode' for both. The exception is a 'static' or 'original' layer. Its individual maps
            are colored by 'membership'.

            A 'quantitative' layer needs these keys:
            - 'accessions': every accession the layer touches. They find the map entries whose
              values set the range of each scale.
            - 'sample_values': the value of each accession in each sample, or None for a layer
              without samples.
            - 'unified_values': the value of each accession, for a layer without samples. All of its
              maps show these values.
            - 'aggregate': reduces the values of an element's accessions to one value.
            - 'cmap': the colormap.
            - 'reverse_overlay': True to draw lower values on top of higher ones.
            - 'colorbar_label': the label of the colorbar.

            A 'quantitative' layer can also have these keys:
            - 'category_cmap': the colormap of the individual maps. It defaults to 'cmap'.
            - 'value_limits' and 'category_value_limits': the limits of the 'unified' map's scale
              and of the individual maps' scale ('_make_quantitative_norm').
            - 'value_center' and 'category_value_center': the centers of the same two scales.
            - 'category_value_center_flag': the option that the individual maps' center came from.
            - 'value_period': the period after which the values repeat, or None. Each scale of
              values then runs from 0 to the period, and derived compound colors are averaged on a
              circle. With 'element_normalize', the scale of the individual maps shows offsets. It
              runs from -period / 2 to period / 2, with a center tick at 0.
            - 'unified_spread': True where the 'unified' map shows how far values that repeat are
              spread apart ('CIRCULAR_SPREADS'). That scale then runs from 0 to the largest spread
              there can be, and derived compound colors are averaged as plain numbers.
            - 'group_samples': the samples of each group, for a layer with samples. It is None in an
              ungrouped run.
            - 'sample_aggregate' and 'group_aggregate': the summaries that pool an element's values
              across samples and across groups ('_summarize_entry_values'). Each is None where it
              summarizes presence. 'group_aggregate' is also None in an ungrouped run.
              'sample_summary' and 'group_summary' hold their names for messages.
            - 'element_normalize': rescales each category's value for an element against the
              values of all categories ('_normalize_entry_value').
            - 'category_colorbar_label': the label of the individual maps' colorbar. It defaults to
              'colorbar_label'. A normalization makes the scale show a different quantity, so it
              needs its own label.
            - 'unified_colorbar_label': the label of the 'unified' map's colorbar. It defaults to
              'colorbar_label'. A spread makes that scale show a different quantity, so it needs
              its own label.

            A 'membership' layer needs these keys:
            - 'membership': the sources that contain each accession.
            - 'source_accessions': the accessions of each source.
            - 'color_hexcode': the color of an individual source's map.

            A 'membership' layer can also have these keys:
            - 'colormap', 'colormap_limits', 'colormap_scheme' and 'reverse_overlay', which set
              how presence is colored.
            - 'scheme_options': the option that chooses the presence scheme, which differs by input.
              It defaults to 'PRESENCE_SCHEME_OPTIONS'.
            - 'category_colors': a color for each category, in place of a colormap. These colors
              color the membership scale and each category's own map. They are checked against the
              categories here.
            - 'category_combo_colors': colors for combinations of categories, which override the
              mix of 'category_colors'. It is needed with 'category_colors'.
            - 'category_colors_flag': the option that 'category_colors' came from. It is needed
              with 'category_colors'.

            A 'static' or 'original' layer needs 'membership', 'source_accessions' and
            'color_hexcode'. 'color_hexcode' is 'original' for an 'original' layer.

            A 'single' layer needs 'accessions' and 'color_hexcode'.

            A layer that is 'quantitative' in one context and 'membership' in the other has the keys
            of both.

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

        def _unified_scale_drawn(layer: dict) -> bool:
            # A layer with no category dimension of its own colors the individual maps from the
            # 'unified' scale. That scale is therefore still worked out where the 'unified' map is
            # not drawn. The options bounding and centering it still apply.
            return draw_unified_maps or (
                layer['category_mode'] == 'quantitative' and layer['sample_values'] is None
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
            # Every category reaches an ordinary run's output, if only as color on the 'unified'
            # map. Without that map, the ones left out of a subset request are drawn nowhere at
            # all.
            if not draw_unified_maps and len(draw_categories) < len(set(categories)):
                self.run.warning(
                    f"Maps were asked for {len(draw_categories)} of the {len(set(categories))} "
                    f"{category_noun}s. The 'unified' map is not drawn in this run, so the other "
                    f"{category_noun}s appear nowhere in the output."
                )
            # Without the 'unified' map's panel, a grid of one category holds a single map. That
            # is not an error, but it is worth saying what the file will look like.
            if not draw_unified_maps and len(draw_grid_categories) == 1:
                self.run.warning(
                    f"The 'unified' map is not drawn in this run, so no grid leads with its "
                    f"panel. Only one {category_noun} was named for the grid. Each grid will hold "
                    f"a single map on a page. It is the same map that '--draw-individual-files' "
                    f"writes, with the {category_noun}'s name over it."
                )
        else:
            draw_files_categories = []
            draw_grid_categories = []
            draw_categories = []
        draw_category_maps = has_categories and (
            draw_individual_files is not False or draw_grid is not False
        )

        # The 'unified' map is what a run draws when it is asked for nothing else. Leaving it out
        # means the individual maps are the whole of the output. Without them there is no output at
        # all. Say so before an empty directory is created.
        if not draw_unified_maps and not draw_categories:
            if has_categories:
                raise ConfigError(
                    f"No maps would be drawn. The 'unified' map summarizing every {category_noun} "
                    f"is not drawn here. Neither is a map for any individual {category_noun}. Ask "
                    f"for those maps with '--draw-individual-files' and/or '--draw-grid', or let "
                    f"the 'unified' map be drawn."
                )
            raise ConfigError(
                "No maps would be drawn. The 'unified' map is not drawn here. This run has a "
                "single source of data, so that map is the only one it has to draw. Please let it "
                "be drawn."
            )

        # Limits on, and a center for, the scale the individual maps share have nothing to act on
        # when this run draws no individual map, just as the group-map coloring options have nothing
        # to act on then. Either one accepted and quietly dropped would leave no sign in the output
        # that it was ignored, so say so instead.
        if not draw_category_maps:
            # A normalization is checked before the limits and the center are, since one centered on
            # zero supplies a center of its own.
            normalized_layers = [
                layer for layer in layers if layer.get('element_normalize') is not None
            ]
            if normalized_layers:
                flags = ', '.join(
                    f"'--{layer['element_type']}-element-normalization'"
                    for layer in normalized_layers
                )
                raise ConfigError(
                    f"A normalization was given for the values of the individual samples or groups "
                    f"({flags}), which rescales each of their values against the values of all of "
                    f"them, but this run draws no map for any individual sample or group, so there "
                    f"is nowhere for a rescaled value to be drawn. Ask for those maps with "
                    f"'--draw-individual-files' and/or '--draw-grid'. Note that the 'unified' map "
                    f"summarizes the unnormalized sample or group values, so it is never "
                    f"rescaled."
                )
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

        # The options bounding and centering the 'unified' map's own scale are the mirror case. A
        # run leaving out that map never draws that scale. The exception is a layer with no category
        # dimension. Its individual maps take that very scale, so the options still apply.
        if not draw_unified_maps:
            for model_key, flag_suffix, subject_phrase, remedy_phrase in (
                ('value_limits', 'value-limits', 'Limits were given for', 'Bound'),
                ('value_center', 'value-center', 'A center was given for', 'Center')
            ):
                affected_layers = [
                    layer for layer in layers
                    if layer.get(model_key) is not None and not _unified_scale_drawn(layer)
                ]
                if not affected_layers:
                    continue
                flags = ', '.join(
                    f"'--{layer['element_type']}-{flag_suffix}'" for layer in affected_layers
                )
                raise ConfigError(
                    f"{subject_phrase} the color scale of the 'unified' map ({flags}). That map is "
                    f"not drawn here, so its scale is never drawn either. {remedy_phrase} the "
                    f"scale shared by the maps of the individual {category_noun}s instead, with "
                    f"'--reaction-category-{flag_suffix}'/'--compound-category-{flag_suffix}', or "
                    f"let the 'unified' map be drawn."
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
            self._resolve_category_colors(
                colored_layers, categories, category_noun, draw_unified_maps=draw_unified_maps
            )
            group_ramps = grouped and (
                grouped_membership['group_colormap'] == GROUP_COLORMAP_FROM_CATEGORY
            )
            for layer in colored_layers:
                # These colors are used in two places. The first is the 'unified' map's membership
                # scale. That scale gives a band to each category, and to each combination of them.
                # The second is the individual maps. Without groups, each category's own map is
                # drawn in that category's color. With groups, each group's own map is drawn in a
                # ramp running to the group's color, but only when '--group-colormap category' asked
                # for that ramp. The two checks below skip a layer whose colors reach either place.
                # A layer reaching neither would leave its colors file unused.
                if draw_unified_maps and layer['unified_mode'] == 'membership':
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

        # With groups, a static color or reference map colors only apply to the 'unified' map. A
        # group map colored by source count cannot utilize a static color or the reference map color
        # scheme, so leaving out the 'unified' map refuses these options.
        static_layers = [
            layer for layer in layers if layer['unified_mode'] in ('static', 'original')
        ]
        if static_layers and grouped and not draw_unified_maps:
            flags = ', '.join(_dedup([
                "'--original-color'" if layer['unified_mode'] == 'original'
                else f"'--{layer['element_type']}-color'" for layer in static_layers
            ]))
            raise ConfigError(
                f"A static color or reference map colors were requested ({flags}). The only map "
                f"these can color is the 'unified' map, and that map is not drawn here. One color "
                f"cannot distinguish a group's {membership_singular}s. Each group's own map is "
                f"instead colored by the number of the group's {membership_singular}s containing "
                f"an element, styled by '--group-colormap'/'--group-reverse-overlay'. Please drop "
                f"the color, or let the 'unified' map be drawn."
            )

        if static and grouped and draw_unified_maps:
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

        # The overwrite check looks at every directory this run writes map files into. Checking
        # the 'unified' directory alone is not enough. A run skipping the 'unified' map writes no
        # such directory, so that check would pass while individual maps and grids were
        # overwritten.
        unified_dir = os.path.join(output_dir, UNIFIED_SUBDIR)
        check_dirs = [unified_dir] if draw_unified_maps else []
        check_dirs += [
            os.path.join(output_dir, INDIVIDUAL_SUBDIR, category) for category in draw_categories
        ]
        if draw_grid is not False:
            check_dirs.append(os.path.join(output_dir, GRID_SUBDIR))
        pathway_numbers = self._find_maps(check_dirs, patterns=pathway_numbers)
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
                # Shared with the drawing passes below, which then have nothing left to work out:
                # this pass asks for the normalized values of every element of every drawn map, and
                # a normalized value covers all of the categories at once.
                layer['_normalized_cache'] = {}
                # The same holds for each element's values, and for the elements whose summary is
                # undefined ('_summarize_entry_values').
                layer['_element_cache'] = {}
                layer['_undefined_unified'] = set()
                layer['_undefined_category'] = set()
                # The 'unified' map and the maps of the groups show summaries. A summary can be
                # undefined for every element of a KEGG map. An example is a standard deviation
                # where each element has a value in a single sample. The layer then has nothing to
                # color there, so the maps where it has a summarized value are recorded.
                layer['_unified_maps'] = set()
                layer['_category_maps'] = {
                    group: set() for group in (layer.get('group_samples') or {})
                }
            for pathway_number in pathway_numbers:
                self.progress.update(pathway_number)
                pathway = self._get_pathway(pathway_number)
                for layer in norm_layers:
                    use_reaction = layer['use_reaction_attribute']
                    unified_valued = (
                        layer['unified_mode'] == 'quantitative' and _unified_scale_drawn(layer)
                    )
                    per_category = (
                        has_categories and layer['category_mode'] == 'quantitative'
                        and layer['sample_values'] is not None
                    )
                    normalize = layer.get('element_normalize')
                    for entry in self._find_element_entries(
                        pathway, use_reaction, layer['accessions']
                    ):
                        category_values, unified_value = self._summarize_entry_values(entry, layer)
                        if unified_value is not None:
                            layer['_unified_maps'].add(pathway_number)
                            if unified_valued:
                                layer['_unified_vals'].append(unified_value)
                        for group, group_maps in layer['_category_maps'].items():
                            if group in category_values:
                                group_maps.add(pathway_number)
                        if per_category and normalize is not None:
                            # A normalization sets each category's value from all of them at once,
                            # so the scale must span the calculated normalized values.
                            layer['_category_vals'].extend(
                                self._normalize_entry_value(
                                    entry, layer, categories, normalize,
                                    cache=layer['_normalized_cache']
                                ).values()
                            )
                        elif per_category:
                            layer['_category_vals'].extend(
                                category_values[category] for category in categories
                                if category in category_values
                            )
            self.progress.end()
            unaffected_clause = (
                " The 'unified' map is unaffected, being drawn from the unnormalized values."
            ) if draw_unified_maps else ""
            for layer in norm_layers:
                if layer['sample_values'] is not None:
                    self._warn_undefined_summaries(layer, draw_unified_maps, draw_category_maps)
                # No values at all means no element of any drawn map has one, so the layer colors
                # nothing and gets no colorbar: either its accessions are absent from these maps, or
                # its aggregation was undefined everywhere (the standard deviation of a single
                # value, say). Say so rather than leaving a blank map to be puzzled over. Where a
                # normalization colors those maps it is named first, and is asked about on its own:
                # the maps it colors come out blank whether or not the 'unified' map has values,
                # that map being drawn from the unnormalized values, so a test of both together
                # would let a blank set of individual maps pass unremarked.
                offsets = (
                    layer.get('value_period') is not None
                    and layer.get('element_normalize') is not None
                )
                if layer.get('element_normalize') is not None and not layer['_category_vals']:
                    cause = (
                        "An offset from a circular mean is undefined where the values cancel out, "
                        "such as 6 and 18 with a period of 24." if offsets else
                        f"A ratio is undefined wherever the value it is measured against is not a "
                        f"positive number, and a z-score wherever an element is found in a single "
                        f"{category_noun}."
                    )
                    self.run.warning(
                        f"Nothing on the maps of the individual {category_noun}s could be colored "
                        f"by '--{layer['element_type']}-element-normalization', so those maps "
                        f"carry no colors from the '{layer['colorbar_label']}' column of the "
                        f"{layer['element_type']} layer and no scale was drawn for them. Either "
                        f"none of the drawn maps contains its accessions, or the normalization is "
                        f"undefined for every map element. {cause}{unaffected_clause}"
                    )
                elif offsets:
                    # An element whose values cancel out has no circular mean to be offset from. It
                    # is left uncolored on the map of every sample or group that has it. That would
                    # look like an element the samples lack, so the user is told how many there
                    # are. The element values and the offsets are cached under the same key.
                    undefined = {
                        key for key, normalized in layer['_normalized_cache'].items()
                        if not normalized and layer['_element_cache'][key][0]
                    }
                    if undefined:
                        examples = ', '.join(sorted({
                            kegg_id for key in undefined for kegg_id in key
                            if kegg_id in layer['accessions']
                        })[:5])
                        self.run.warning(
                            f"'--{layer['element_type']}-element-normalization' was undefined for "
                            f"{len(undefined)} distinct map element(s). Their accessions in the "
                            f"file include these: {examples}. "
                            f"{self._undefined_cause('circular_mean')} These elements have no "
                            f"circular mean to be offset from. They are left uncolored on the maps "
                            f"of the individual {category_noun}s."
                        )
                elif not layer['_unified_vals'] and not layer['_category_vals']:
                    self.run.warning(
                        f"Nothing on the maps could be colored by the values of the "
                        f"'{layer['colorbar_label']}' column of the {layer['element_type']} layer, "
                        f"so no color scale was drawn for it. Either none of the drawn maps "
                        f"contains its accessions, or the aggregation reducing them is undefined "
                        f"for every map element — the standard deviation of a single value, for "
                        f"instance."
                    )
                # A period fixes each scale of values to run from 0 to the period. The limits are
                # set here, not in the layer model. The checks on limits read a limit in the model
                # as one that was given, and these were not. A spread on the 'unified' map runs from
                # 0 to the largest spread there can be.
                period = layer.get('value_period')
                if period is None:
                    value_limits = layer.get('value_limits')
                elif layer.get('unified_spread'):
                    value_limits = (0.0, self._largest_angular_deviation(period))
                else:
                    value_limits = (0.0, period)
                # The offsets of a normalization run from -period / 2 to period / 2, with a center
                # tick at 0. The two ends are the same point.
                category_value_center = layer.get('category_value_center')
                if period is None:
                    category_value_limits = layer.get('category_value_limits')
                elif layer.get('element_normalize') is not None:
                    category_value_limits = (-period / 2, period / 2)
                    category_value_center = 0.0
                else:
                    category_value_limits = (0.0, period)
                if layer['unified_mode'] == 'quantitative' and _unified_scale_drawn(layer):
                    norm, vmin, vmax, clamped_low, clamped_high = self._make_quantitative_norm(
                        layer['_unified_vals'], value_limits,
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
                    # A limit sets its end of the scale, so the colorbar labels every end a limit
                    # set, and marks those that values lie past. The top of a spread's scale is the
                    # largest spread there can be, 5.40 with a period of 24. A label there would
                    # print every tick to three decimals, so the top is left to the usual ticks.
                    layer['_unified_limited'] = tuple(
                        limit is not None for limit in (value_limits or (None, None))
                    )
                    if layer.get('unified_spread'):
                        layer['_unified_limited'] = (True, False)
                    layer['_unified_clamped'] = (clamped_low, clamped_high)
                    layer['_unified_center'] = layer.get('value_center')
                if layer['category_mode'] != 'quantitative':
                    continue
                if has_categories and layer['sample_values'] is not None:
                    norm, vmin, vmax, clamped_low, clamped_high = self._make_quantitative_norm(
                        layer['_category_vals'], category_value_limits,
                        f"--{layer['element_type']}-category-value-limits",
                        center=category_value_center,
                        center_flag=layer.get(
                            'category_value_center_flag',
                            f"--{layer['element_type']}-category-value-center"
                        ),
                        subject=(
                            f"the {layer['element_type']} layer of the maps of the individual "
                            f"{category_noun}s"
                        )
                    )
                    layer['_category_norm'] = norm
                    layer['_category_range'] = (vmin, vmax)
                    layer['_category_limited'] = tuple(
                        limit is not None for limit in (category_value_limits or (None, None))
                    )
                    layer['_category_clamped'] = (clamped_low, clamped_high)
                    layer['_category_center'] = category_value_center
                else:
                    # A layer without a category dimension is constant across the category maps.
                    layer['_category_norm'] = layer['_unified_norm']
                    layer['_category_range'] = layer['_unified_range']
                    layer['_category_limited'] = layer['_unified_limited']
                    layer['_category_clamped'] = layer['_unified_clamped']
                    layer['_category_center'] = layer['_unified_center']

        # A context colored by membership needs its by-count/by-membership color scheme over the
        # categories. Only the 'unified' map draws such a scale (and its colorbar) from the
        # categories themselves; the presence coloring of grouped individual maps counts each
        # group's own sources, precomputed below. A run leaving that map out works out no such
        # scale. The limits that keep such a scale legible do not apply to it.
        membership_layers = [
            layer for layer in layers
            if draw_unified_maps and layer['unified_mode'] == 'membership'
        ]
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

        def _reaction_derived(layer, mode, cmap_key='cmap', spread=False):
            # How a reaction layer derives compound colors on a reaction-only global/overview map.
            # The derived colors come from the same scale as the reactions they are derived from, so
            # the caller names the context's colormap: the two can differ ('category_cmap'). The
            # caller also says whether the context shows a spread, which is a plain number and so is
            # averaged as one.
            if layer['element_type'] != 'reaction':
                return None
            if mode == 'quantitative' and spread:
                cmap = layer[cmap_key]
                return ('average', cmap.reversed() if layer['reverse_overlay'] else cmap)
            if mode == 'quantitative':
                cmap = layer[cmap_key]
                # Values that repeat after a period are averaged on a circle, whatever the colormap.
                # Their scale runs from 0 to the period, so a circle of the scale is one period.
                # Clock times of 0.1 h and 23.9 h then give 0 h, not 12 h. Without a period, a
                # cyclic colormap also averages on a circle, since it has the same color at both
                # ends of the scale. The colormap is checked before it is reversed. Reversing
                # 'clocktime_r' would name it 'clocktime_r_r'.
                circular = (
                    layer.get('value_period') is not None or self._is_cyclic_colormap(cmap)
                )
                transfer = 'circular_average' if circular else 'average'
                return (transfer, cmap.reversed() if layer['reverse_overlay'] else cmap)
            return ('high', None)

        def _sample_accessions(layer, samples):
            # Return the accessions that the samples have values for. The samples are taken in turn,
            # and the accessions of each sample in turn. This order must be the same in every run.
            # Map entries are colored in this order. Each color is stored with a priority. The
            # priority is the position of the element's value on the color scale. Two close values
            # can get the same color. The color then keeps the priority of the element colored last.
            # Priorities set which elements are drawn on top. On global and overview maps, they also
            # set the colors of compounds. Python orders a set of strings differently in each run,
            # so a set would not work here.
            return dict.fromkeys(
                accession for sample in samples for accession in layer['sample_values'][sample]
            )

        def _unified_spec(layer):
            mode = layer['unified_mode']
            if mode == 'quantitative':
                def entry_value(entry):
                    return self._summarize_entry_values(entry, layer)[1]
                if layer['sample_values'] is None:
                    entry_keys, maps = layer['unified_values'], None
                else:
                    # A summary colors only the maps where some element has a summarized value.
                    # Grouped, the samples are taken group by group. The 'unified' value is built
                    # from the groups in this order too.
                    group_samples = layer['group_samples']
                    samples = layer['sample_values'] if group_samples is None else [
                        sample for members in group_samples.values() for sample in members
                    ]
                    entry_keys = _sample_accessions(layer, samples)
                    maps = layer['_unified_maps']
                return {
                    'element_type': layer['element_type'],
                    'use_reaction_attribute': layer['use_reaction_attribute'],
                    'entry_keys': entry_keys,
                    'pathway_numbers': maps,
                    'colorer': self._quantitative_colorer(
                        entry_value, layer['_unified_norm'], layer['cmap'],
                        layer['reverse_overlay'], center=layer['_unified_center']
                    ),
                    'derived_compound': _reaction_derived(
                        layer, mode, spread=layer.get('unified_spread', False)
                    )
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
                    'keep_highest_priority': True,
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
                normalize = layer.get('element_normalize')
                maps = None
                if layer['sample_values'] is None:
                    # A layer without samples draws the same values on every map.
                    entry_keys = layer['unified_values']

                    def entry_value(entry):
                        return self._summarize_entry_values(entry, layer)[1]
                else:
                    # A sample's map finds the elements of the accessions the sample has values for.
                    # A group's map finds those of the accessions its samples have values for. It
                    # colors only the maps where some element has a value for the group. A
                    # normalization can still give a category no value for an element, where it had
                    # nothing to rescale or where the rescaling itself was undefined.
                    if layer['group_samples'] is None:
                        entry_keys = layer['sample_values'][category]
                    else:
                        entry_keys = _sample_accessions(layer, layer['group_samples'][category])
                        maps = layer['_category_maps'][category]
                    if normalize is None:
                        def entry_value(entry):
                            return self._summarize_entry_values(entry, layer)[0].get(category)
                    else:
                        def entry_value(entry):
                            return self._normalize_entry_value(
                                entry, layer, categories, normalize,
                                cache=layer['_normalized_cache']
                            ).get(category)
                return {
                    'element_type': layer['element_type'],
                    'use_reaction_attribute': layer['use_reaction_attribute'],
                    'entry_keys': entry_keys,
                    'pathway_numbers': maps,
                    'colorer': self._quantitative_colorer(
                        entry_value, layer['_category_norm'], layer['category_cmap'],
                        layer['reverse_overlay'], center=layer['_category_center']
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
                        'keep_highest_priority': True,
                        'derived_compound': _reaction_derived(layer, 'membership')
                    })
                group_layer_membership[group] = (
                    specs, group_color_priorities[group], group_scale_tops[group]
                )

        def _category_specs(category):
            # Grouped membership specs are precomputed; a layer colored by value or a single color
            # in its per-group context colors its own map and rides along on the same map.
            if grouped_presence:
                return group_layer_membership[category][0] + [
                    _category_spec(layer, category) for layer in layers
                    if layer['category_mode'] != 'membership'
                ]
            return [_category_spec(layer, category) for layer in layers]

        def _element_map_specs():
            # The specs of each kind of map that '_draw_map_elements' draws, with a phrase naming
            # it. A run in the reference map's own colors draws its other maps with
            # '_draw_map_kos_original_color'.
            if draw_unified_maps and not original_run:
                yield (
                    "the 'unified' map" if has_categories else None, True,
                    [_unified_spec(layer) for layer in layers]
                )
            for category in draw_categories:
                if grouped_presence or not original_run:
                    yield f"{category_noun} '{category}'", False, _category_specs(category)

        self._check_reserved_colors(
            layers, pathway_numbers, group_color_priorities=group_color_priorities,
            group_element_types={layer['element_type'] for layer in group_membership_layers},
            drawn_categories=draw_categories
        )
        self._check_derived_compound_colors(
            layers, pathway_numbers, _element_map_specs,
            group_color_priorities=group_color_priorities
        )

        # Per-layer colorbars for the unified map (layer-prefixed so two layers do not collide),
        # plus a shared colorbar for the category maps of each layer colored quantitatively there.
        # Each context is keyed by its own mode, so a layer summarized by presence in the unified
        # map gets a presence colorbar there and a continuous one for its category maps. A run
        # leaving out the unified map draws neither of its bars. The exception is a layer with no
        # category dimension of its own. Its individual maps take the unified scale, so that
        # colorbar is drawn as the only key to their colors.
        for layer in layers:
            if layer['unified_mode'] == 'quantitative' and _unified_scale_drawn(layer):
                vmin, vmax = layer['_unified_range']
                if vmin is not None:
                    limited_low, limited_high = layer['_unified_limited']
                    clamped_low, clamped_high = layer['_unified_clamped']
                    self._draw_quantitative_colorbar(
                        layer['cmap'], vmin, vmax,
                        os.path.join(output_dir, f"colorbar_{layer['name']}.pdf"),
                        layer.get('unified_colorbar_label', layer['colorbar_label']),
                        limited_low=limited_low, limited_high=limited_high,
                        clamped_low=clamped_low, clamped_high=clamped_high,
                        center=layer['_unified_center']
                    )
            elif layer['unified_mode'] == 'membership' and draw_unified_maps:
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
                and layer['sample_values'] is not None
            ):
                vmin, vmax = layer['_category_range']
                if vmin is not None:
                    limited_low, limited_high = layer['_category_limited']
                    clamped_low, clamped_high = layer['_category_clamped']
                    self._draw_quantitative_colorbar(
                        layer['category_cmap'], vmin, vmax,
                        os.path.join(
                            output_dir, f"colorbar_{layer['name']}_{colorbar_category_suffix}.pdf"
                        ),
                        layer.get('category_colorbar_label', layer['colorbar_label']),
                        limited_low=limited_low, limited_high=limited_high,
                        clamped_low=clamped_low, clamped_high=clamped_high,
                        center=layer['_category_center']
                    )

        # Draw the unified map (the single map when there are no categories).
        if draw_unified_maps:
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
                _, category_color_priorities, category_scale_top = (
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
                specs = _category_specs(category)
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
                specs = _category_specs(category)
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
                source_type=grid_source_type if grid_source_type is not None else category_noun,
                include_unified=draw_unified_maps
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

        if draw_unified_maps:
            count = sum(drawn['unified'].values()) if drawn['unified'] else 0
            self.run.info(
                f"Number of 'unified' maps drawn incorporating data from all {unified_plural}",
                count
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
        drawn_categories: List[str],
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

        drawn_categories : List[str]
            The categories drawn on their own maps. A category color is the color of that
            category's own map, so only these colors are checked here. Where the 'unified' map is
            drawn, its own scale is checked as well, and that scale carries every category's color.

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

        # The value scale shared by the maps of the individual categories colors only those maps.
        value_scales = [('_unified_norm', '_unified_center', '_unified_vals', 'cmap')]
        if drawn_categories:
            value_scales.append(
                ('_category_norm', '_category_center', '_category_vals', 'category_cmap')
            )
        for layer in layers:
            if layer['unified_mode'] == 'original':
                continue
            # Every color the layer can stage: sampled from its colormap at the values it will color
            # by, taken from the scale it colors presence by, or its one fixed color.
            staged: Set[str] = set()
            for norm_key, center_key, values_key, cmap_key in value_scales:
                if norm_key not in layer or not layer[values_key]:
                    continue
                norm = layer[norm_key]
                cmap = layer[cmap_key]
                # A scale with no norm spans a single value. The colorers give every element the
                # middle color where the scale is centered. Otherwise they give it the top color.
                if norm is None:
                    fraction = 0.5 if layer[center_key] is not None else 1.0
                    staged.add(mcolors.rgb2hex(cmap(fraction)))
                    continue
                # A run with many samples can have millions of values here. A colormap has only a
                # few hundred colors. The values are therefore colored all at once. Each distinct
                # color is then converted to a hex code once.
                rgba = cmap(norm(np.asarray(layer[values_key])))
                staged.update(mcolors.rgb2hex(color) for color in np.unique(rgba, axis=0))
            if '_colors' in layer:
                staged.update(color for color, _ in layer['_colors'][1])
            category_colors = layer.get('category_colors')
            if category_colors is not None:
                staged.update(category_colors[category] for category in drawn_categories)
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

    def _check_derived_compound_colors(
        self,
        layers: List[dict],
        pathway_numbers: List[str],
        element_map_specs: Callable[[], Iterable[Tuple[Union[str, None], bool, List[dict]]]],
        group_color_priorities: Dict[str, List[Tuple[str, float]]] = None
    ) -> None:
        """
        Check that no compound is given a reserved color derived from the reactions it touches.

        A global or overview map drawn without a compound layer colors each compound from the
        reactions it touches ('kgml.Pathway.color_associated_compounds'). Where the reactions are
        colored by value, a compound takes the color at the average position of its reactions on
        the color scale. Otherwise it takes the color of its reaction drawn on top. That color can
        be the one the map keeps for the compounds it does not highlight. This can happen even where
        no reaction has that color. 'kgml' then refuses to draw the map. This check finds such a
        clash before any map or colorbar is drawn. '_check_reserved_colors' does the same for the
        colors of the layers themselves.

        The derived colors depend on which reactions touch each compound. Each map is therefore
        colored here without being drawn. That costs about as much as drawing the map, less the
        rendering. So it is done only on the maps whose reserved color a derived color can reach at
        all. A derived color is always a color of a reaction colormap or a fixed color of the
        reactions.

        Parameters
        ==========
        layers : List[dict]
            The layer models. Compound colors are derived only where no layer is a compound layer.

        pathway_numbers : List[str]
            The maps about to be drawn. Only global and overview maps derive compound colors.

        element_map_specs : Callable[[], Iterable[Tuple[Union[str, None], bool, List[dict]]]]
            Gives the layer specs of each kind of map that '_draw_map_elements' is about to draw,
            such as the 'unified' map or the map of one sample. Each comes with a phrase naming that
            kind of map, and with True if it is the 'unified' map. The phrase is None where the run
            has no categories and so draws a single map. The specs are built only if the maps have
            to be colored.

        group_color_priorities : Dict[str, List[Tuple[str, float]]], None
            The colors of each group's individual maps ('_group_map_colors'). A compound on a
            group's map can take one of them.
        """
        # A compound layer colors the compounds itself, so no compound color is derived.
        if any(layer['element_type'] == 'compound' for layer in layers):
            return
        reserved: Dict[str, str] = {}
        for pathway_number in pathway_numbers:
            is_global = re.match(GLOBAL_MAP_ID_PATTERN, pathway_number) is not None
            is_overview = re.match(OVERVIEW_MAP_ID_PATTERN, pathway_number) is not None
            if is_global or is_overview:
                reserved[pathway_number] = kgml.canonical_color(
                    kgml.reserved_recolor_colors('g' if is_global else 'w', is_overview)['compound']
                )
        if not reserved:
            return

        # Every color a compound can be derived in. A reversed overlay derives compound colors from
        # the reversed colormap ('_reaction_derived'). Each colormap is therefore read both ways.
        derivable: Set[str] = set()
        for layer in layers:
            for mode_key, cmap_key in (
                ('unified_mode', 'cmap'), ('category_mode', 'category_cmap')
            ):
                if layer[mode_key] != 'quantitative':
                    continue
                for cmap in (layer[cmap_key], layer[cmap_key].reversed()):
                    derivable.update(mcolors.rgb2hex(color) for color in cmap(np.arange(cmap.N)))
            if '_colors' in layer:
                derivable.update(color for color, _ in layer['_colors'][1])
            if layer.get('category_colors') is not None:
                derivable.update(layer['category_colors'].values())
            if layer.get('color_hexcode') is not None:
                derivable.add(layer['color_hexcode'])
        for color_priorities in (group_color_priorities or {}).values():
            derivable.update(color for color, _ in color_priorities)
        derivable = set(map(kgml.canonical_color, derivable))
        # Only the maps whose reserved color can be derived are colored.
        reserved = {
            pathway_number: color for pathway_number, color in reserved.items()
            if color in derivable
        }
        if not reserved:
            return

        self.progress.new("Checking the compound colors derived from reactions")
        for drawing, unified, specs in element_map_specs():
            reaction_spec = next(spec for spec in specs if spec['element_type'] == 'reaction')
            mode, colormap = reaction_spec['derived_compound']
            for pathway_number, reserved_color in reserved.items():
                self.progress.update(
                    pathway_number if drawing is None else f"{pathway_number}, {drawing}"
                )
                pathway = self._get_pathway(pathway_number)
                color_priority, _ = self._stage_map_elements(pathway_number, pathway, specs)
                # These are the steps 'kgml.Pathway.set_color_priority' takes to derive compound
                # colors, short of the last one. That step recolors the compounds left without a
                # color. It is where 'kgml' refuses a derived color that is reserved.
                pathway.set_color_priority(
                    color_priority,
                    recolor_unprioritized_entries='g' if pathway.is_global_map else 'w'
                )
                pathway.color_associated_compounds(mode, colormap=colormap)
                derived = {
                    kgml.canonical_color(bgcolor)
                    for _, bgcolor in pathway.color_priority.get('compound', {}).get('circle', {})
                }
                if reserved_color not in derived:
                    continue
                self.progress.end()
                drawn_for = '' if drawing is None else f" ({drawing})"
                remedy = self._derived_compound_remedy(
                    layers, unified, group_color_priorities, reserved_color
                )
                raise ConfigError(
                    f"Some compounds would be colored {reserved_color} on pathway map "
                    f"{pathway_number}{drawn_for}. The map keeps this color for the compounds it "
                    f"does not highlight. They could not be told apart from those compounds. No "
                    f"compound layer is drawn. The map therefore colors each compound from the "
                    f"reactions it touches. Where the reactions are colored by value, a compound "
                    f"takes the color at the average position of its reactions on the color scale. "
                    f"Otherwise it takes the color of its reaction drawn on top. {remedy}"
                )
        self.progress.end()

    @staticmethod
    def _derived_compound_remedy(
        layers: List[dict],
        unified: bool,
        group_color_priorities: Union[Dict[str, List[Tuple[str, float]]], None],
        reserved_color: str
    ) -> str:
        """
        Say which option to change where a compound would be derived in a reserved color.

        The option is the one that set the colors of the reactions on that kind of map
        ('_check_derived_compound_colors').

        Parameters
        ==========
        layers : List[dict]
            The layer models. The run has no compound layer, so its one reaction layer is here.

        unified : bool
            True for the 'unified' map, False for the map of one category.

        group_color_priorities : Union[Dict[str, List[Tuple[str, float]]], None]
            The colors of each group's individual maps. A grouped run colors those maps by counts.

        reserved_color : str
            The reserved color that a compound would take.

        Returns
        =======
        str
            The sentences that name the option and what to do with it.
        """
        layer = next(layer for layer in layers if layer['element_type'] == 'reaction')
        mode = layer['unified_mode'] if unified else layer['category_mode']
        if reserved_color == '#ffffff':
            colormaps = (
                "Grayscale colormaps and those running to pure white, such as 'Greys', 'hot' and "
                "'bone', reach it."
            )
        else:
            colormaps = "Colormaps with grays in them, such as 'Greys' and 'RdGy', reach it."
        if not unified and group_color_priorities and mode == 'membership':
            return (
                f"These maps color the number of the group's own sources containing each "
                f"reaction, styled by '--group-colormap'. Choose a colormap that does not reach "
                f"{reserved_color}. {colormaps} A ramp from a pale tint to the group's own color "
                f"('--group-colormap {GROUP_COLORMAP_FROM_CATEGORY}') can start further from "
                f"white instead. Raise the first of the two limits that option also takes."
            )
        if mode == 'quantitative' and not unified:
            return (
                f"Choose a colormap for '--reaction-category-colormap' that does not reach "
                f"{reserved_color}. Without that option, these maps use the colormap given to "
                f"'--reaction-colormap'. {colormaps}"
            )
        # A color per category colors the presence of the categories in place of a colormap.
        if mode == 'membership' and layer.get('category_colors') is not None:
            return (
                f"Remove {reserved_color} from the file given to "
                f"'{layer['category_colors_flag']}'."
            )
        if mode == 'quantitative' or (mode == 'membership' and unified):
            return (
                f"Choose a colormap for '--reaction-colormap' that does not reach "
                f"{reserved_color}. {colormaps}"
            )
        return f"Choose a color other than {reserved_color} for '--reaction-color'."

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
        draw_unified_maps: bool = True,
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
            draw_unified_maps=draw_unified_maps,
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
        color_priority: dict,
        keep_highest_priority: bool = False
    ) -> None:
        """
        Set an entry's graphics colors for a layer and register them in 'color_priority'.

        Reaction (ortholog) entries are lines in global and overview maps and boxes or lines in
        standard maps; compound entries are circles. The registered '(fgcolor, bgcolor) -> priority'
        keeps the entry from being treated as unprioritized and recolored to the background. A color
        has one priority on a map. Two entries given one color can bring different priorities.

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

        keep_highest_priority : bool, False
            If True, a color that is already registered keeps the higher of its priority and this
            one. If False, this priority replaces it. The last entry given a color then sets its
            priority.
        """
        def _register(entry_type: str, graphics_type: str, colors: Tuple[str, str]) -> None:
            priorities = color_priority.setdefault(entry_type, {}).setdefault(graphics_type, {})
            if keep_highest_priority and colors in priorities:
                priorities[colors] = max(priorities[colors], priority)
            else:
                priorities[colors] = priority

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
                _register('compound', 'circle', colors)
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
            _register('ortholog', graphics_type, colors)

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
        entry_value: Callable[[kgml.Entry], Union[float, None]],
        norm: Union[mcolors.Normalize, None],
        cmap: mcolors.Colormap,
        reverse_overlay: bool,
        center: Union[float, None] = None
    ):
        """
        Build a colorer that colors an Entry by a continuous value.

        See '_draw_map_elements' for the colorer contract. 'entry_value' gives an Entry's value on
        this map. It is the element's value in one sample, a summary of samples or groups
        ('_summarize_entry_values'), or a value rescaled across them ('_normalize_entry_value').
        None leaves the Entry uncolored. The color is 'cmap' sampled at the normalized value, and
        the priority is that fraction, or its complement under 'reverse_overlay', which 'clip=True'
        on the norm keeps in [0, 1]. A degenerate range (no norm) leaves every element at the top of
        the colormap, except on a centered scale, where the one value the range collapsed to is the
        center itself and so takes the middle color, the very color the centering was asked for.
        """
        def colorer(entry: kgml.Entry) -> Union[Tuple[str, float], None]:
            value = entry_value(entry)
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
        by the position of their exact combination (by membership). Two counts can share a color, on
        a continuous count scale or with a cyclic colormap. A map keeps one priority per color. So
        the spec of a presence layer asks '_draw_map_elements' to keep the highest priority given to
        a color on the map ('keep_highest_priority'). An Entry in no source, or, when grouped, in no
        qualifying group, is left uncolored.
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
            There is one spec per layer. Earlier specs are drawn beneath later ones.

            Each spec has these keys:
            - 'element_type': 'reaction' or 'compound'.
            - 'use_reaction_attribute': True if the layer's accessions are KEGG reaction IDs rather
              than KO IDs.
            - 'entry_keys': the accessions the layer touches. They find the map entries to color. A
              compound layer also uses them to warn about compounds that are drawn only as
              rectangles. Those cannot be colored. Entries are colored in the order of these
              accessions. A color keeps the priority of the last entry given it. So if the colorer
              can give a color more than one priority, this order must not change between runs.
            - 'colorer': a function that takes an Entry. It returns a '(color_hexcode, priority)'
              pair, or None to leave the Entry uncolored.
            - 'derived_compound': how a reaction layer colors compounds by the reactions that touch
              them. It is used only on a global or overview map without a compound layer. It is a
              '(mode, colormap)' pair, such as ('average', cmap), ('circular_average', cmap) or
              ('high', None). The mode is passed to 'kgml.Pathway.set_color_priority' as
              'color_associated_compounds'. A compound layer has None here; it is not read.

            A spec can also have these keys:
            - 'pathway_numbers': the maps on which the layer has a value to color. On any other map
              the layer colors nothing. It then does not count as matching the map's accessions. A
              summary of samples or groups can be undefined for every element of a map. This key
              lets such a map be skipped.
            - 'keep_highest_priority': True to give each color the highest priority of the entries
              given that color on the map. Otherwise a color keeps the priority of the last entry
              given it. A presence layer sets this. Two of its counts can share a color. Its entries
              are colored in the order of its input, such as the rows of a text file. Without this
              key, that order could change the map.

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

        color_priority, found_entries = self._stage_map_elements(
            pathway_number, pathway, layer_specs
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

    def _stage_map_elements(
        self,
        pathway_number: str,
        pathway: kgml.Pathway,
        layer_specs: List[dict]
    ) -> Tuple[dict, bool]:
        """
        Stage the colors of every layer on one map, without applying them.

        Each layer's colorer colors the entries the layer matches. Drawing a map
        ('_draw_map_elements') and checking the compound colors a map derives
        ('_check_derived_compound_colors') both stage colors here. The check therefore sees the
        colors the drawing will use.

        Parameters
        ==========
        pathway_number : str
            Numeric ID of the map.

        pathway : kgml.Pathway
            The map, freshly loaded. Its matched entries are given their colors.

        layer_specs : List[dict]
            One spec per layer, with the keys '_draw_map_elements' describes. Earlier specs are
            drawn beneath later ones.

        Returns
        =======
        Tuple[dict, bool]
            Two values:
            - The '{entry_type: {graphics_type: {(fg, bg): priority}}}' dictionary, to be passed to
              'kgml.Pathway.set_color_priority'.
            - Whether any layer matched an entry of the map.
        """
        color_priority: dict = {}
        found_entries = False
        for spec in layer_specs:
            if spec.get('pathway_numbers') is not None and (
                pathway_number not in spec['pathway_numbers']
            ):
                continue
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
                    pathway, entry, spec['element_type'], color_hexcode, priority, color_priority,
                    keep_highest_priority=spec.get('keep_highest_priority', False)
                )
        return color_priority, found_entries

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
        return self._draw_map_elements(
            pathway_number, self._element_presence_specs(layers), output_dir,
            draw_map_lacking_data=draw_map_lacking_data
        )

    def _element_presence_specs(self, layers: List[dict]) -> List[dict]:
        """
        Build the layer specs of '_draw_map_element_presence' from its layers.

        Drawing the maps and checking their colors first ('_map_kos_fixed_colors') both build the
        specs here, so the check sees the colors the drawing will use.
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
        return specs

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

    def _find_maps(self, check_dirs: List[str], patterns: List[str] = None) -> List[str]:
        """
        Find the numeric IDs of maps to draw, checking that the maps can be drawn without
        overwriting output already there.

        Parameters
        ==========
        check_dirs : List[str]
            The directories this run writes map files into. Each is checked for files that would
            be overwritten. A run gives every one of them, so that the check covers everything it
            is about to write.

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
            for check_dir in check_dirs:
                for pathway_number in pathway_numbers:
                    out_basename = self._map_basename(pathway_number)
                    if self.pathway_categorization is None:
                        out_path = os.path.join(check_dir, out_basename)
                    else:
                        out_path = os.path.join(
                            check_dir, *self.pathway_categorization[pathway_number], out_basename
                        )
                    if os.path.exists(out_path):
                        raise ConfigError(
                            f"Output files would be overwritten in {check_dir}. Either delete the "
                            f"contents of that directory, or use the option to overwrite output "
                            f"destinations."
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
        source_type: str = 'unknown',
        include_unified: bool = True
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

        include_unified : bool, True
            If True, each grid leads with the 'unified' map summarizing every category, labeled
            'all' or 'pangenome'. If False, a grid holds the individual maps alone. The run has not
            drawn the 'unified' map.
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
            # With the 'unified' map left out, a grid whose every panel would be a blank filler has
            # nothing to show. The pathway's data is in a category this grid does not hold. Fillers
            # even out a grid that has content. They do not make a grid out of nothing.
            if not include_unified and not any(
                drawn['individual'].get(category, {}).get(pathway_number, False)
                for category in draw_grid_categories
            ):
                continue
            pathway_basename = self._map_basename(pathway_number)
            in_paths: List[str] = []
            labels: List[str] = []
            if include_unified:
                unified_dir = os.path.join(output_dir, UNIFIED_SUBDIR)
                if self.pathway_categorization is None:
                    unified_map_path = os.path.join(unified_dir, pathway_basename)
                else:
                    out_dir = os.path.join(
                        unified_dir, *self.pathway_categorization[pathway_number]
                    )
                    unified_map_path = os.path.join(out_dir, pathway_basename)
                if not os.path.exists(unified_map_path):
                    continue
                in_paths.append(unified_map_path)
                labels.append('pangenome' if source_type == 'pangenome' else 'all')

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
        limited_low: bool = False,
        limited_high: bool = False,
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

        limited_low : bool, False
            If True, a value limit set the bottom of the range, so 'vmin' is labeled. A limit is a
            value someone chose, and an automatic tick rarely falls on it.

        limited_high : bool, False
            The same for the top of the range, labeling 'vmax'.

        clamped_low : bool, False
            If True, values lie below the limit at the bottom of the range, so 'vmin' is labeled
            with 'CLAMPED_MIN_PREFIX'. Its color is what everything below it is drawn in, and a bare
            number there would claim the color stands for that value alone. The end is labeled
            whether or not 'limited_low' is also given.

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

        labeled_low = limited_low or clamped_low
        labeled_high = limited_high or clamped_high
        if labeled_low or labeled_high or center is not None:
            # An end set by a limit is labeled with the limit. It is marked '≤' or '≥' where values
            # lie past it, because its color is what everything past the limit is drawn in. A
            # centered range is labeled at its center, which is otherwise the one place on the bar a
            # reader cannot find. The rest of the bar keeps Matplotlib's own ticks, run through the
            # bar's own formatter, so that such a bar is typeset exactly like a plain one; any of
            # them falling all but on top of one of these labels is dropped, since of the two labels
            # in that spot, this is the one a reader needs. A tick landing on an end that was NOT
            # labeled is kept, which is what labels that end of the bar.
            span = vmax - vmin
            marked = (
                ([vmin] if labeled_low else [])
                + ([center] if center is not None else [])
                + ([vmax] if labeled_high else [])
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
            # The formatter's own call prints any value below 1e-8 as zero. Every label of a bar of
            # tiny values would then read zero. Its format string is applied directly instead. The
            # offset and multiplier are both zero here, so nothing else differs. Only a value within
            # rounding error of zero, as tick arithmetic leaves it, is printed as zero.
            tick_labels = [
                formatter.fix_minus(formatter.format % (0 if abs(tick) < 1e-9 * span else tick))
                for tick in ticks
            ]
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
