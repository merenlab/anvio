#!/usr/bin/env python
"""
Make KEGG pathway maps incorporating data sourced from external inputs.

Most maps are drawn by the element engine. '_map_elements' builds a spec for each layer of each
drawn map from the layer's models, as '_element_presence_specs' does for a single fixed color, and
'_draw_map_elements' stages the specs on each KEGG map and draws it. The terms below are defined as
the code uses them, those of the data before those of drawing.

Terms
=====
map
    A KEGG pathway map, a 'kgml.Pathway', whose number the code passes as 'pathway_number'. Its kind
    ('MAP_KINDS') is standard, overview or global, and reactions are drawn as boxes only on standard
    maps.

accession
    A KEGG ID of a KO, a KEGG reaction or a compound, given in a layer's input or named by an Entry
    ('_get_entry_kegg_ids'). Reaction data hold only KO IDs or, with 'use_reaction_attribute', only
    KEGG reaction IDs.

map element
    A box, line or circle that stands for a reaction or a compound on a map. It is a Graphics
    ('kgml.Graphics') of an Entry ('kgml.Entry'), the KGML node that names its accessions in 'name',
    or in 'reaction' with 'use_reaction_attribute'. The code sometimes calls the Entry itself the
    element, while 'kgml.Element' and "Graphics element" mean nodes of KGML. A layer's
    'element_type' is 'reaction' for ortholog Entries and 'compound' for compound Entries.

layer
    Everything drawn from the input of one element type, such as the text file of a text layer. A
    run has at most two layers, which are drawn from the bottom up, reaction before compound. The
    key columns of a text file are 'accession', 'gene_id' and 'sample', and options name its border,
    bar and weights columns. The one other column, if any, is its value column, and its values are
    'numeric' or 'categorical' ('VALUE_TYPES').

individual, source
    A sample, contigs database, or genome that gets a map of its own, any of which is also called a
    source. In a grouped run the individuals are instead the groups of sources from a groups-txt
    file. 'individuals' lists the names in color-assignment order, and it is None in a run with a
    single source of data.

'unified' map, context
    The 'unified' map summarizes every individual, and it is the only map of a run with a single
    source of data. A drawn map is the 'unified' map or the map of one individual, as the
    'drawn_map' column of the key table names it beside the KEGG 'map'. Keys such as 'unified_mode'
    and 'individual_mode' set apart the two contexts, the 'unified' map and the individual maps, and
    "kind of map" sometimes means a context rather than a map's kind. The 'context' argument of the
    slicing and ring code is instead a sample or group, or None for the 'unified' map or a layer
    without samples.

aggregation, summary
    The reductions of many values to one. Within a sample, the gene aggregation reduces an
    accession's rows, and the accession aggregation an element's accessions. The sample summary
    pools samples into the map of their group, or into the 'unified' map in an ungrouped run. The
    group summary pools groups into the 'unified' map. Each reduction of categorical values is a
    rule ('CATEGORICAL_RULES'). A model holds the accession aggregation in 'aggregate' and the
    summaries in 'sample_aggregate' and 'group_aggregate'. An element normalization rescales the
    values of individual maps and is not a reduction.

presence
    Whether an individual contains an element. A source contains it if it has any of the element's
    accessions, and a group if enough of its sources do ('group_threshold', a fraction). A text
    layer without a value column shows presence, and so does a fill whose summary is a presence name
    ('count', 'count_continuous', or 'membership') or is left at its default. How presence is
    colored is its scheme ('by_count', 'by_count_continuous', or 'by_membership').

channel, role
    What one column of a text layer drives, named by its role. The value column colors the 'fill' of
    an element (by presence where there is none), a border column colors its 'border', and a bar
    column sets the length of its 'bar'. The fill is the inside of a box or circle, or a line
    itself. '_txt_channel' builds a channel from 'TXT_CHANNEL_ROLES', which gives each role's option
    stem ('flag_stem'), message words, key-file 'name' and defaults.

border, ring
    A border is the channel of a border column, and a layer's 'border' holds its model. It is drawn
    as a ring ('kgml.Graphics.border'), whole or in bands or segments, inside the outline of a box
    on a standard map or of a circle on a standard or overview map.

bar
    The channel of the reaction bar column, and a layer's 'bars' is a list holding its model. A bar
    is drawn inside a box on a standard map, and it runs from 0 to its value. Its style
    ('BAR_STYLES') draws it in the fill's color over a gray track, or as a dark bar on a white strip
    under the whole fill. On a scale with values on both sides of 0, bars grow both ways from a zero
    line, and the boxes of the 'fill' style drop their labels. "Bar" can also mean a colorbar.

model
    The dict that says how one channel of a layer is colored, or measured for bars, in each context,
    with the keys that the '_map_elements' docstring lists. '_build_txt_model' builds the models of
    a text layer, and '_map_element_membership' those of databases and pangenomes, which have no
    channel and count as fills. A layer is passed as the model of its fill, and the variable 'layer'
    may hold any model. Keys that start with '_' hold working state.

mode
    How a model colors its elements in one context ('unified_mode', 'individual_mode'):
    'quantitative' (numeric values), 'categorical', 'membership' (presence, by count or combination
    on a map that summarizes sources, or in one color on a source's own map), 'single' or 'static'
    (one presence color), 'original' (KEGG's reference colors) or 'none' (no fill, only rings).
    'static' and 'original' are modes of the 'unified' map only. A border or bar model's mode is
    None, not 'none', in a context whose maps do not show its values. "Membership" also names a
    presence summary, a scheme, a model key and a colors file.

scale
    The mapping of values, or of presence counts or combinations, to colors by a 'norm' and a
    colormap, or of values to the lengths of bars. A value model has one scale for the 'unified' map
    and one that the individual maps share, which is the same scale in a layer without samples. Each
    group's presence map has a count scale of its own. A 'norm' is a 'matplotlib.colors.Normalize',
    not an element normalization.

slice
    One unit's part of a reaction element's fill, in the color of its own value: a column of a box
    ('kgml.Graphics.glyph') or a dash along a line ('kgml.Graphics.slices'). The unit is an
    accession or a gene ('SLICE_UNITS') and takes a share of the element. An element with too many
    slices is drawn whole in the color of a stand-in value. A ring can be sliced as well, with the
    fill or on its own ('slice_layout').

spec
    A drawing instruction for one layer on one drawn map, used on each KEGG map and drawn by
    '_draw_map_elements', whose docstring lists its keys. '_layer_spec' builds most specs from
    models, and '_element_presence_specs' those of a single fixed color. Its 'entry_keys' are the
    accessions the layer touches, which find the Entries to color. Specs are staged in list order,
    and reactions are drawn beneath compounds whatever the order.

colorer, priority
    A colorer takes a 'kgml.Entry' and returns '(color_hexcode, priority)', or None to leave it
    uncolored. Colors of higher priority are drawn on top, and the priority of a quantitative color
    is the fraction of the scale at its value, reversed by 'reverse_overlay'.

stage
    To write a map's colors, slices, rings and bars onto a freshly loaded 'kgml.Pathway' and gather
    their priorities ('_stage_map_elements'). 'kgml.Pathway.set_color_priority' then orders the
    Entries, derives compound colors, and recolors each Entry left uncolored in a reserved color of
    its kind of map ('kgml.reserved_recolor_colors').

derived compound colors
    The colors that compound circles take from the reaction lines that end at them, on a global or
    overview map where a reaction layer stages fills and no compound layer does
    ('_derives_compounds'). The spec key 'derived_compound' says how.

key (file)
    A colorbar, a legend of categorical values, or a bar key: a file named from its model's 'name'
    ('ColorbarDrawer'), or from its group for a group's count colorbar. The key table is the text
    file that lists the slices of fills and rings ('_write_slice_keys'). The key table of unlabeled
    boxes lists the boxes that bars drew without their labels ('_write_label_keys'). Elsewhere a key
    is usually a dict key, or one of the key columns of a text file.
"""

import os
import re
import csv
import html
import fitz
import json
import math
import base64
import shutil
import hashlib
import numpy as np
import pandas as pd
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker

from io import BytesIO
from argparse import Namespace
from collections import Counter
from PIL import Image, features
from itertools import combinations
# The text that pandas reads as a missing value by default, such as 'NA' and 'NaN'.
from pandas._libs.parsers import STR_NA_VALUES
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
# rate from midnight to noon. It falls at the same rate after noon, so two hours equally far from
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
# global or overview map where no compound layer colors the fills, each compound takes a color from
# the reaction lines that touch it. With one of these colormaps, the values of those lines are
# averaged on a circle. A layer with a period ('--reaction-value-period') averages them on a circle
# with any colormap. '_is_cyclic_colormap' matches colormap names against this list. Matplotlib's
# cyclic colormaps are not listed. People often use one of them, 'hsv', for values that do not
# repeat.
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
# individual map shows a point in the period for each element, its circular mean in one sample or
# group. Its scale runs from 0 to the period, with the cyclic 'clocktime' by default, and its
# colorbar is labeled with the value column's name. The 'unified' map shows the spread of those
# points. Its scale runs from 0 to the largest spread there can be, 5.40 with a period of 24. It
# takes a colormap that is not cyclic, since the smallest and the largest spread are not the same
# point. By default this is 'plasma_r' for the fill and 'viridis' for a border
# ('TXT_CHANNEL_ROLES'). Its colorbar is labeled '{value} angular deviation'. On global and overview
# maps, its compound colors are averaged as plain numbers, not on a circle.
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
# an aggregation of values. 'count' and 'count_continuous' color a count from the same colormap. A
# count scale calls for a sequential colormap. With such a colormap, the two give a count the same
# color. They differ only in how the scale is drawn. 'count' gives each count a band of a discrete
# colorbar. This requires one distinguishable color for each count. 'count_continuous' draws a
# gradient from the first count to the last. A gradient works for any number of samples or groups
# ('_membership_layer_colors').
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

# The statistic a normalization compares each element's value against. It is computed over the
# samples or groups in which that element has a value at all. These reduce a list like
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

# The kinds of values the value column of a text file can hold ('--reaction-value-type',
# '--compound-value-type'). 'numeric' values are numbers. They are colored along a scale.
# 'categorical' values are names, such as 'up' and 'down'. Each name takes a color of its own. The
# kind is declared, not detected from the file. Names written as numbers, such as 1 and 2, then stay
# names. A typo in a column of numbers still raises an error.
VALUE_TYPES: Tuple[str, ...] = ('numeric', 'categorical')

# The colors that categorical values take when no colors file gives them ('_resolve_value_colors').
# The values take the colors in order. The first nine colors are those of Matplotlib's 'tab10'
# without its gray. The next nine are the light colors of 'tab20' without its light gray, in the
# same order of hues. The first value is blue, for instance, and the tenth value is light blue. Gray
# is left out because global maps draw the elements they do not highlight in gray.
DEFAULT_VALUE_COLORS: Tuple[str, ...] = (
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#bcbd22',
    '#17becf',
    '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', '#dbdb8d',
    '#9edae5'
)

# The colors that categorical border values take when no colors file gives them
# ('TXT_CHANNEL_ROLES'). These are the colors of Matplotlib's 'Dark2' without its gray. They differ
# from the colors of the fill, so that a color means a value of one column only. The last two are
# close to the default compound and reaction colors, so they come last.
DEFAULT_BORDER_VALUE_COLORS: Tuple[str, ...] = (
    '#d95f02', '#7570b3', '#e6ab02', '#a6761d', '#1b9e77', '#e7298a', '#66a61e'
)

# The colormap of numeric border values when no colormap is named for them, and the part of it that
# is used. A border lies between two black lines. 'viridis' is nearly black below 0.2, where a
# border would merge with those lines. The fill takes 'plasma_r' instead, so that a color means a
# value of one column only.
DEFAULT_BORDER_COLORMAP = 'viridis'
DEFAULT_BORDER_COLORMAP_LIMITS: Tuple[float, float] = (0.2, 1.0)

# The channels of a layer from a text file, keyed by role ('_txt_channel'). The value column of a
# file colors the fill of each map element. A border column colors the border of the element instead
# ('--reaction-border-column', '--compound-border-column'). A bar column sets the length of a bar
# inside the fill, drawn in the color of the fill ('--reaction-bar-column'). Every channel is built
# by the same methods ('_build_txt_model'). A channel names the options of its values in messages,
# and it gives their defaults. The channels have different defaults, so that a color means a value
# of one column only. A bar takes numbers only, with no colormap, and its scale always includes 0
# ('baseline'). In the text of a role, '{element_type}' stands for 'reaction' or 'compound'.
TXT_CHANNEL_ROLES: Dict[str, dict] = {
    'fill': {
        'flag_stem': '{element_type}',
        'noun': '{element_type}',
        'plural': '{element_type}s',
        'name': '{element_type}s',
        'column': 'value column',
        'values': 'values',
        'scale': 'color scale',
        'presence': True,
        'categorical': True,
        'period': True,
        'baseline': None,
        'default_aggregation': 'sum',
        'default_colormap': 'plasma_r',
        'default_colormap_limits': None,
        'default_centered_colormap': DEFAULT_CENTERED_COLORMAP,
        'default_value_colors': DEFAULT_VALUE_COLORS,
        'choice_salt': (),
        'consequence': 'is left uncolored',
        'consequence_plural': 'are left uncolored'
    },
    'border': {
        'flag_stem': '{element_type}-border',
        'noun': '{element_type} border',
        'plural': '{element_type} borders',
        'name': '{element_type}s_border',
        'column': 'border column',
        'values': 'border values',
        'scale': 'color scale',
        'presence': False,
        'categorical': True,
        'period': True,
        'baseline': None,
        'default_aggregation': 'mean',
        'default_colormap': DEFAULT_BORDER_COLORMAP,
        'default_colormap_limits': DEFAULT_BORDER_COLORMAP_LIMITS,
        'default_centered_colormap': DEFAULT_CENTERED_COLORMAP,
        'default_value_colors': DEFAULT_BORDER_VALUE_COLORS,
        'choice_salt': ('border',),
        'consequence': 'shows no border',
        'consequence_plural': 'show no border'
    },
    'bar': {
        'flag_stem': '{element_type}-bar',
        'noun': '{element_type} bar',
        'plural': '{element_type} bars',
        'name': '{element_type}s',
        'column': 'bar column',
        'values': 'bar values',
        'scale': 'bar scale',
        'presence': False,
        'categorical': False,
        'period': False,
        'baseline': 0.0,
        'default_aggregation': 'mean',
        'default_colormap': None,
        'default_colormap_limits': None,
        'default_centered_colormap': None,
        'default_value_colors': None,
        'choice_salt': ('bar',),
        'consequence': 'is left uncolored',
        'consequence_plural': 'are left uncolored'
    }
}

# The rules that choose one categorical value where several meet ('_choose_categorical_value').
# Values meet among the genes of one accession ('--reaction-gene-aggregation'), among the accessions
# of one map element ('--*-accession-aggregation'), and among the samples on the 'unified' map
# ('--*-sample-summary'). With groups, they meet among the samples of each group
# ('--*-sample-summary'), and then among the groups on the 'unified' map ('--*-group-summary'). On
# global and overview maps where no compound layer colors the fills, they meet among the reactions
# that touch a compound, where the accession rule chooses. Each gene, accession, sample, group, or
# reaction with a value gives one vote. These are the rules:
# - 'plurality_order' takes the most common value, even with less than half of the votes. A tie
#   goes to the value that comes first in the order of the values. That order is the row order of
#   the colors file, or else sorted order ('_resolve_value_colors').
# - 'plurality_random' takes the most common value. A tie is broken at random.
# - 'plurality' takes the most common value. A tie gives no value.
# - 'majority' takes the value with more than half of the votes. Without one, there is no value.
# - 'order' takes the value that comes first in the order of the values, whatever the votes.
# - 'random' takes one of the values at random, whatever the votes.
# Where a rule gives no value, the thing it chooses for has no value. It does not vote at the next
# level, as an undefined numeric summary does not ('_finite_values').
CATEGORICAL_RULES: Tuple[str, ...] = (
    'plurality_order', 'plurality_random', 'plurality', 'majority', 'order', 'random'
)

# The categorical rules as a phrase for messages that list them.
CATEGORICAL_RULES_PHRASE = ', '.join(repr(name) for name in CATEGORICAL_RULES)

# The rule that chooses a categorical value where no rule is given.
DEFAULT_CATEGORICAL_RULE = 'plurality_order'

# The rules that make random choices. A random choice depends on a seed ('--seed').
RANDOM_CATEGORICAL_RULES: Tuple[str, ...] = ('plurality_random', 'random')

# The seed of random choices where no seed is given.
DEFAULT_CATEGORICAL_SEED = 42

# What the reaction elements of a map can be sliced by ('--reaction-slice-by'). A slice is one
# accession of the reaction file's own type, a KO or a KEGG reaction. A slice can instead be one
# gene of the file's 'gene_id' column. Each slice of a box or line takes the color of its own value
# and a share of the element's width or length ('_slice_result'). Only maps that show the values of
# the layer are sliced. Compounds are never sliced.
SLICE_UNITS: Tuple[str, ...] = ('accession', 'gene')

# How the shares of slices are pooled on a map that summarizes several samples, which is the map of
# a group or the 'unified' map ('--reaction-slice-share-summary'). They apply only with a column of
# weights. 'share_mean', the default, averages the shares that the maps of the samples draw.
# 'pooled_sum' adds up the weights of each accession over the samples and divides by the total.
SLICE_SHARE_SUMMARIES: Tuple[str, ...] = ('share_mean', 'pooled_sum')

# How the weights of the genes of an accession are reduced to the weight of the accession in a
# sample ('--reaction-slice-weight-aggregation'). 'sum' is the default.
SLICE_WEIGHT_AGGREGATIONS: Tuple[str, ...] = ('sum', 'mean')

# The options of the border of a layer ('map_kegg_pathways_txt'). 'column' names the border column
# of the layer's text file ('--reaction-border-column', '--compound-border-column'). Every other key
# is the border version of an option of the fill, such as 'colormap' for
# '--reaction-border-colormap'. 'colormap_limits' and 'individual_colormap_limits' come with their
# colormaps, and 'element_normalization_label' comes with its normalization. 'slice_by' divides the
# rings of boxes among accessions or genes ('--reaction-border-slice-by'). A compound file has no
# genes, so a compound border takes no 'gene_aggregation'. Compounds are never sliced, so it takes
# no 'slice_by' either.
BORDER_OPTIONS: Tuple[str, ...] = (
    'column', 'value_type', 'value_colors', 'value_period', 'gene_aggregation',
    'accession_aggregation', 'sample_summary', 'group_summary', 'colormap', 'colormap_limits',
    'individual_colormap', 'individual_colormap_limits', 'value_limits', 'individual_value_limits',
    'value_center', 'individual_value_center', 'element_normalization',
    'element_normalization_label', 'slice_by'
)

# The options of a border that only the reaction layer has.
REACTION_BORDER_OPTIONS: Tuple[str, ...] = ('gene_aggregation', 'slice_by')

# The options of the bars of the reaction layer ('map_kegg_pathways_txt'). 'column' names the bar
# column of the reaction file ('--reaction-bar-column'). 'style' is how a bar is drawn in its box
# ('BAR_STYLES'). Every other key is the bar version of an option of the fill, such as
# 'sample_summary' for '--reaction-bar-sample-summary'. Bars take numbers only, so they take no
# value type, period, colormap, center, or normalization.
BAR_OPTIONS: Tuple[str, ...] = (
    'column', 'style', 'gene_aggregation', 'accession_aggregation', 'sample_summary',
    'group_summary', 'value_limits', 'individual_value_limits'
)

# How a bar is drawn in its box ('--reaction-bar-style'). A 'fill' bar is drawn in the color of the
# fill over a gray track, so the fill shows only as far as the bar reaches. A 'strip' bar is a thin
# dark bar on a white strip along the bottom of the box, and the fill is drawn whole above it. The
# strip keeps both values readable, while a 'fill' bar weights the color of the fill by its length.
# 'fill' is the default.
BAR_STYLES: Tuple[str, ...] = ('fill', 'strip')

# The track of a bar, which is the part of the fill that a bar does not cover. A gray track tells a
# bar of 0 from a box without data, which is white. The end line is drawn centered on the free end
# of every bar, at 0 and at the full length too. A light bar differs from the track by little more
# than its hue, and the line marks its end. Widths are in points.
BAR_TRACK_COLOR = '#E5E5E5'
BAR_END_LINE_COLOR = '#4D4D4D'
BAR_END_LINE_WIDTH = 0.75

# The zero line of a bar scale with values on both sides of 0. Its bars grow both ways from 0, so a
# thin gray line marks 0 inside each box, under the end line. The label of a box of the 'fill' style
# on such a scale is dropped, since the label would hide the zero line and the ends of short bars
# ('LABEL_KEY_COLUMNS'). A strip lies below the label and keeps it. The width is in points.
BAR_ZERO_LINE_COLOR = '#808080'
BAR_ZERO_LINE_WIDTH = 0.5

# The strip of the 'strip' style, along the bottom of the fill of a box, below its label. Its white
# track means 0 only inside a colored fill, and a box without data is white throughout. The bar is
# dark on every fill, and it takes the same end line. The height is in points.
BAR_STRIP_HEIGHT = 5
BAR_STRIP_COLOR = '#333333'
BAR_STRIP_TRACK_COLOR = '#FFFFFF'

# The narrowest that a slice of a box may be, in points, where every slice takes the same share. A
# box with more slices than that allows is drawn whole in one color ('_slice_cap'). Boxes on
# standard maps are 92 points wide, so they hold up to 11 slices. In a layer with a border column,
# slices fit inside the ring, in 85 points, so a box holds up to 10 slices.
MIN_BOX_SLICE_WIDTH = 8

# The most slices that the ring of a box can be divided into where the ring is sliced on its own
# ('--reaction-border-slice-by'). A ring with more slices is drawn whole in one color. The middle
# line of a ring around a box of 92 by 34 points is about 238 points long, so each of 11 equal
# slices is about 21.6 points long.
MAX_RING_SLICES = 11

# The columns of the key table that lists the slices of a layer ('_write_slice_keys'). A slice of a
# gene lists the gene's accessions in the element, joined by commas. 'gene_id' is empty for a slice
# of an accession.
SLICE_KEY_COLUMNS: Tuple[str, ...] = (
    'map', 'drawn_map', 'entry_id', 'x', 'y', 'label', 'slice', 'accession', 'gene_id', 'value',
    'share', 'color', 'drawn_whole'
)

# The columns of the key table that lists the boxes that are drawn without their labels
# ('_write_label_keys'). A box of the 'fill' style on a bar scale with values on both sides of 0
# drops its label. The table names the box by its accessions in the layer's file, joined by commas
# in KGML order, and gives its bar value.
LABEL_KEY_COLUMNS: Tuple[str, ...] = (
    'map', 'drawn_map', 'entry_id', 'x', 'y', 'label', 'accession', 'bar'
)

# The header of the color column of a kegg-membership-colors-txt file and of a kegg-value-colors-txt
# file. The first column holds the names of samples, contigs databases, genomes, or groups in the
# first file, and categorical values in the second. That column can be headed anything, as the item
# column of a groups-txt file can. One membership colors file can then describe samples in one run
# and groups in another without being renamed.
COLORS_COLUMN = 'color'

# What separates the names of a combination in the first column of a kegg-membership-colors-txt
# file. A combination is a set of samples, contigs databases, genomes, or groups. Coloring by
# membership derives the color of a combination by averaging the colors of its members. A row that
# names the combination overrides that color. The separator matches how the membership colorbar
# labels a combination, so an override row can be copied straight from the scale it adjusts.
MEMBERSHIP_COMBO_SEPARATOR = ','

# Coloring by membership needs a color for every possible combination of the samples, contigs
# databases, genomes, or groups. With n of them, there are '2 ** n - 1' combinations. When a
# kegg-membership-colors-txt file gives each of them a color, the combinations are blended rather
# than sampled from a colormap, so no colormap size bounds them. This is the ceiling instead. Above
# it, the combinations are refused before they are enumerated. It is set at the size of a typical
# colormap. No scale can tell apart the combinations of more than a handful of samples, contigs
# databases, genomes, or groups. The ceiling is well past that.
MAX_MEMBERSHIP_COLOR_COMBOS = 256

# The most counts a color scale is drawn in discrete bands of before the count is drawn on a
# continuous scale instead. A band's label is sized to fit the band
# ('ColorbarDrawer.draw_discrete'), so a scale of many counts labels them in type too small to read
# long before it runs out of colors to draw them in: this ceiling is where the labels stop fitting,
# not where the colors do. Past it a bar of bands reads as a gradient anyway, while a gradient
# labels the range rather than every count. Coloring by membership has no such ceiling. It has no
# gradient to fall back to. A color on its scale says which samples, contigs databases, genomes, or
# groups contain an element. A position on a gradient cannot say that.
MAX_DISCRETE_COUNT_BANDS = 40

# The '--group-colormap' value asking for each group's individual maps to be colored by a ramp built
# from that group's own color rather than sampled from a named Matplotlib colormap. Every group's
# ramp is built the same way, differing only in the hue it runs to, so the panels of a grid stay
# comparable while each says which group it is.
GROUP_COLORMAP_TINT = 'tint'

# The fractions of the way from white to a group's own color that its ramp runs between, used when
# 'GROUP_COLORMAP_TINT' is asked for without explicit limits. The pale end stops short of white
# because white is a reserved color of standard and overview maps ('_check_reserved_colors'), and a
# ramp that reached it could not be drawn there at all.
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

# The same for the ruler of a key of bar lengths ('draw_length'). Its labels sit side by side along
# the ruler rather than one above another, so a label of a limit, such as '≥ 35', needs about twice
# the room.
MIN_LENGTH_TICK_SEPARATION_FRACTION = 0.16

# Subdirectories of the output directory, one per role a map file can have: the map unifying every
# source, the map of one individual source, and the grid comparing them. A fourth holds those
# individual maps a second time, arranged into a directory per map rather than per source. Beside
# each of those directories, it can also hold a player of the same maps. Every name directly in the
# output directory is therefore anvi'o's own, while the names that come from user data — samples,
# contigs databases, genomes, groups — are confined to 'individual', where they cannot collide with
# anything anvi'o creates for itself.
UNIFIED_SUBDIR = 'unified'
INDIVIDUAL_SUBDIR = 'individual'
GRID_SUBDIR = 'grid'
COLLATED_SUBDIR = 'by_map'

# Where '--categorize-files' nests map files in BRITE subdirectories, this subdirectory of links
# gathers every one of them back into a single place to browse ('_link_map_flat').
FLAT_SUBDIR = 'all_maps'

# The colorbar keying the maps of one individual source, written beside them in its directory.
INDIVIDUAL_COLORBAR_BASENAME = 'colorbar.pdf'

# The player of a map gathered by map ('_draw_map_players'). It is an HTML file written beside the
# map's directory under 'by_map', and it shows the same maps one after another. The page is filled
# in from a template that ships with anvi'o.
PLAYER_EXTENSION = '.html'
PLAYER_TEMPLATE_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), 'data', 'static', 'template',
    'kegg-map-player.tmpl'
)
# The number of seconds a player first shows each frame. The viewer can change it on the page.
PLAYER_SECONDS_PER_FRAME = 1
# The note on a frame of a player whose source has no map drawn, so the frame shows the bare map.
PLAYER_BARE_NOTE = 'No data on this map'

# The kinds of map whose reaction elements and compound circles can be labeled ('--reaction-labels'
# and '--compound-labels'). Global maps match 'GLOBAL_MAP_ID_PATTERN' (011XX), overview maps match
# 'OVERVIEW_MAP_ID_PATTERN' (012XX and 013XX), and all other maps are standard. By default, only the
# reaction elements of standard maps are labeled.
MAP_KINDS: Tuple[str, ...] = ('standard', 'overview', 'global')
DEFAULT_REACTION_LABEL_MAP_KINDS: Tuple[str, ...] = ('standard',)
DEFAULT_COMPOUND_LABEL_MAP_KINDS: Tuple[str, ...] = ()

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

    animate_files_by_map : bool
        Beside each directory of maps gathered by 'collate_files_by_map', write an HTML player that
        shows the maps one after another.

    reaction_label_map_kinds : Set[str]
        The kinds of map whose reaction elements are labeled, among those in 'MAP_KINDS'.

    compound_label_map_kinds : Set[str]
        The kinds of map whose compound circles are labeled, among those in 'MAP_KINDS'.

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

    player_drawer : HTMLPlayerDrawer
        Writes HTML files that show input PDF files one after another.

    ignore_compound_rectangles : bool
        If True, do not draw KGML compound Entry rectangle Graphics. These are found in a small
        number of KGML files (see 00121, 00621, 01052, 01054), and when rendered by 'anvio.kgml' via
        'Bio.Graphics.KGML_vis.KGMLCanvas' have the effect of obscuring underlying drawings of
        compound structures in the base map image.

    visible_line_lengths : Dict[Tuple[str, str, int], float]
        The length of the part of each line that compound circles leave visible, as the cap on
        slices finds it ('_get_visible_line_length'). Keys are the map number, the entry ID, and the
        position of the line among the Graphics elements of the entry.

    compound_circle_grids : Dict[str, anvio.kgml.CircleGrid]
        The compound circles of each map whose lines the cap on slices has checked, sorted into a
        grid. Keys are map numbers.
    """
    def __init__(
        self,
        kegg_dir: str = None,
        overwrite_output: bool = FORCE_OVERWRITE,
        name_files: bool = False,
        categorize_files: bool = False,
        collate_files_by_map: bool = False,
        animate_files_by_map: bool = False,
        reaction_label_map_kinds: Iterable[str] = None,
        compound_label_map_kinds: Iterable[str] = None,
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

        animate_files_by_map : bool, False
            Beside each directory of maps gathered by 'collate_files_by_map', write an HTML player
            that shows the maps one after another. Without 'collate_files_by_map', no player is
            written.

        reaction_label_map_kinds : Iterable[str], None
            The kinds of map whose reaction elements are labeled, among 'standard', 'overview', and
            'global'. A label is KEGG's name for the element, which is its first KO. None labels
            standard maps only ('DEFAULT_REACTION_LABEL_MAP_KINDS'). Their boxes are blank without
            labels, since the fills hide the text of the base map.

        compound_label_map_kinds : Iterable[str], None
            The kinds of map whose compound circles are labeled with their KEGG IDs, among
            'standard', 'overview', and 'global'. None labels no map
            ('DEFAULT_COMPOUND_LABEL_MAP_KINDS').

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

        # A player holds each frame as a WebP image. Pillow writes WebP only if it was built with
        # libwebp, which is the case for the usual builds, but not for every build.
        if animate_files_by_map and not features.check('webp'):
            raise ConfigError(
                "A player written by '--animate-files-by-map' holds its frames as WebP images, but "
                "the installed Python package 'pillow' cannot write WebP. Please install a build "
                "of 'pillow' with WebP support, or drop this flag."
            )

        if reaction_label_map_kinds is None:
            reaction_label_map_kinds = DEFAULT_REACTION_LABEL_MAP_KINDS
        if compound_label_map_kinds is None:
            compound_label_map_kinds = DEFAULT_COMPOUND_LABEL_MAP_KINDS
        # The program offers only the known kinds of map, so this check is for other callers.
        reaction_label_map_kinds = set(reaction_label_map_kinds)
        compound_label_map_kinds = set(compound_label_map_kinds)
        for parameter_name, label_map_kinds in (
            ('reaction_label_map_kinds', reaction_label_map_kinds),
            ('compound_label_map_kinds', compound_label_map_kinds)
        ):
            unknown_map_kinds = sorted(label_map_kinds.difference(MAP_KINDS), key=str)
            if unknown_map_kinds:
                unknown = ', '.join(f"'{map_kind}'" for map_kind in unknown_map_kinds)
                known = ', '.join(f"'{map_kind}'" for map_kind in MAP_KINDS[:-1])
                known += f", and '{MAP_KINDS[-1]}'"
                raise ConfigError(
                    f"'{parameter_name}' names kinds of map that anvi'o does not know: {unknown}. "
                    f"The kinds of map are {known}."
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
        self.visible_line_lengths: Dict[Tuple[str, str, int], float] = {}
        self.compound_circle_grids: Dict[str, kgml.CircleGrid] = {}
        self.colorbar_drawer = ColorbarDrawer(overwrite_output=overwrite_output)
        self.grid_drawer = PDFGridDrawer(overwrite_output=overwrite_output)
        self.player_drawer = HTMLPlayerDrawer(overwrite_output=overwrite_output)

        self.name_files = name_files
        self.categorize_files = categorize_files
        self.collate_files_by_map = collate_files_by_map
        self.animate_files_by_map = animate_files_by_map
        self.reaction_label_map_kinds = reaction_label_map_kinds
        self.compound_label_map_kinds = compound_label_map_kinds
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
        element_type: Literal['reaction', 'compound'],
        value_type: Literal['numeric', 'categorical'] = 'numeric',
        slice_weights: Union[str, None] = None,
        slice_by: Union[str, None] = None,
        border_column: Union[str, None] = None,
        border_value_type: Literal['numeric', 'categorical'] = 'numeric',
        border_slice_by: Union[str, None] = None,
        bar_column: Union[str, None] = None
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

        A value column can hold categorical values instead of numbers ('VALUE_TYPES'). Such values
        are names. They are read as text, so 'NA' is a name like any other. Every row needs one.

        Every file is read with quotes left as text. A pair of double quotes enclosing a whole field
        or column name is dropped, as R writes them ('_unquote'). Other quotes are kept, so a number
        or a weight that holds one is refused, and a name keeps it. Each column name may appear only
        once in the header line ('_read_txt_header'). A name cannot start with '__', which marks the
        columns that this method adds.

        A column of slice weights is named, not detected ('--reaction-slice-weights'). It is left
        out of the search for the value column. It cannot be a key column, and the file still needs
        a value column besides it, unless the rings of a border column are sliced
        ('--reaction-border-slice-by'). Every row needs a weight, a finite number of 0 or more. The
        weights are read as text first, as every column is, and then converted to numbers.

        Slices of genes ('--reaction-slice-by gene') need a 'gene_id' column. Every row then needs a
        gene ID. The spaces around each ID are dropped before repeated rows are looked for. Without
        slices of genes, the 'gene_id' column is read as it is. A file without a value column is not
        checked for gene IDs, since slices of presence are refused ('_build_slice_model'). Rings
        sliced by gene ('--reaction-border-slice-by gene') need gene IDs in the same way, even in a
        file without a value column.

        A border column is named, not detected ('--reaction-border-column',
        '--compound-border-column'). Its values color the borders of map elements, and the value
        column still colors their fills. It is left out of the search for the value column, so a
        file with a border column and no other value column colors its fills by presence. It cannot
        be a key column or the column of slice weights. A row can leave its border cell empty. That
        row has no border value, but its other columns are read as usual. A column with no border
        value at all is refused. A file with a border column is read as text, so 'NA' is a name in
        every column, as in a file of categorical values. A numeric border value must be a finite
        number.

        A bar column is named, not detected ('--reaction-bar-column'). Its values set the lengths of
        bars drawn inside the fills of map elements. It is left out of the search for the value
        column, as a border column is, so a file with a bar column and no other value column colors
        its fills by presence. It cannot be a key column, the column of slice weights, or the border
        column. A row can leave its bar cell empty, and then it has no bar value. A column with no
        bar value at all is refused. A bar value must be a finite number, which can be negative. A
        file with a bar column is read as text, as a file with a border column is.

        Parameters
        ==========
        path : str
            Path to the tab-delimited per-layer text file.

        element_type : Literal['reaction', 'compound']
            Which layer the file colors, selecting the accession-typing and key-column rules.

        value_type : Literal['numeric', 'categorical'], 'numeric'
            What the value column holds. 'categorical' requires a value column.

        slice_weights : Union[str, None], None
            The name of the column of slice weights, or None for a file without one.

        slice_by : Union[str, None], None
            What the elements are sliced by ('SLICE_UNITS'), or None where they are not sliced.

        border_column : Union[str, None], None
            The name of the border column, or None for a file without one.

        border_value_type : Literal['numeric', 'categorical'], 'numeric'
            What the border column holds.

        border_slice_by : Union[str, None], None
            What the rings of the border are sliced by ('SLICE_UNITS'), or None where they are not
            sliced.

        bar_column : Union[str, None], None
            The name of the bar column, or None for a file without one.

        Returns
        =======
        dict
            The layer table, with these keys:
            - 'element_type': the 'element_type' argument.
            - 'reaction_source': 'KO' for a reaction file of KO IDs, and 'Reaction' for a reaction
              file of KEGG reaction IDs. It is None for a compound file.
            - 'df': a DataFrame of every row of the file. It keeps the file's own columns. It adds
              an '__accession' column, which holds each row's accession without the spaces around
              it. A file with a 'sample' column also gets a '__sample' column, made the same way. A
              file of categorical values also gets a '__value' column, made the same way from the
              value column. Every column is read as text. Numbers are converted later
              ('_aggregate_accession_quantities'). With slices of genes, a file with a value column
              also gets a '__gene_id' column, which holds each row's gene ID, made the same way. A
              file whose rings are sliced by gene gets it too.
            - 'value_column': the name of the value column. It is None for a file without one,
              whose layer is colored by presence.
            - 'value_type': the 'value_type' argument.
            - 'sample_names': the names in the 'sample' column, sorted, each listed once. It is
              None for a file without a 'sample' column.
            - 'slice_weights': the 'slice_weights' argument. A file with a column of weights also
              gets a '__slice_weight' column in 'df', which holds each row's weight as a number.
            - 'border_column' and 'border_value_type': the arguments. A file with a border column
              also gets a '__border_value' column in 'df'. It holds each row's border value, as a
              number or as text without the spaces around it, and is missing where the row has none
              ('_channel_view').
            - 'border_empty': the number of rows without a border value.
            - 'bar_column': the 'bar_column' argument. A file with a bar column also gets a
              '__bar_value' column in 'df'. It holds each row's bar value as a number, and is
              missing where the row has none ('_channel_view').
            - 'bar_empty': the number of rows without a bar value.
        """
        # A single-column presence file (just 'accession') is valid, so 'is_file_tab_delimited'
        # (which rejects a line with no tab) is too strict here; the 'accession'-column check below
        # gives a clear error if the file is mis-delimited.
        filesnpaths.is_file_exists(path)

        self.progress.new(f"Loading the kegg-{element_type}-txt file")
        self.progress.update("...")

        # An empty file has no header row for pandas to find, which it reports as an exception
        # rather than as an empty table. Pandas reads text such as 'NA' as a missing value by
        # default. Categorical values keep such text as a name. Pandas also takes a double quote to
        # start a quoted field, which can hold tabs and line breaks, so a stray quote would merge the
        # rows after it into one field. Every file is read with quotes left as text instead. Then
        # only a pair of quotes enclosing a whole field is dropped, as R writes it and as in a colors
        # file ('_read_value_colors_txt'). Any other quote stays in its field, where a number or a
        # weight that holds it is refused. A missing value in a file of numbers is left missing. A
        # file with a border or bar column is read as text too, so that a value such as 'NA' in that
        # column is kept as it is written.
        categorical = value_type == 'categorical'
        try:
            df = pd.read_csv(
                path, sep='\t', dtype=str,
                keep_default_na=not categorical and border_column is None and bar_column is None,
                quoting=csv.QUOTE_NONE
            )
        except pd.errors.EmptyDataError:
            self.progress.end()
            raise ConfigError(
                f"The kegg-{element_type}-txt file at '{path}' is empty. It needs a header row with "
                f"an 'accession' column, and a row for each accession to color."
            )
        df.columns = [self._unquote(column) for column in df.columns]
        df = df.applymap(lambda cell: cell if pd.isna(cell) else self._unquote(cell))
        # Quotes around nothing, as pandas writes a missing value when it quotes every field, leave
        # empty text. Outside a file of categorical values, such a cell is missing.
        if not categorical:
            df = df.replace('', np.nan)

        # Pandas renames a repeated column name, so a second 'coverage' column would be read as
        # 'coverage.1', and it could be taken for the value column without a word. The header line
        # is read as it is written to find such names. Blank names are left out, since pandas names
        # each blank column itself, such as 'Unnamed: 2'. A name that starts with '__' could clash
        # with a column that this method adds, such as '__accession'.
        names = self._read_txt_header(path)
        repeated = [name for name in dict.fromkeys(names) if name and names.count(name) > 1]
        if repeated:
            self.progress.end()
            raise ConfigError(
                f"The header line of the kegg-{element_type}-txt file at '{path}' names these "
                f"columns more than once: {', '.join(repr(name) for name in repeated)}. A name may "
                f"appear only once."
            )
        reserved = [name for name in names if name.startswith('__')]
        if reserved:
            self.progress.end()
            raise ConfigError(
                f"The header line of the kegg-{element_type}-txt file at '{path}' has these names "
                f"that start with '__': {', '.join(repr(name) for name in reserved)}. Names that "
                f"start with '__' are kept for anvi'o's own columns, such as '__accession'. Please "
                f"rename these columns."
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

        # A column of slice weights is named by an option. It must be a column of the file, and it
        # cannot be a key column.
        weights_flag = f'--{element_type}-slice-weights'
        if slice_weights is not None:
            if slice_weights in ('accession', 'gene_id', 'sample'):
                self.progress.end()
                raise ConfigError(
                    f"'{weights_flag}' names the '{slice_weights}' column of the "
                    f"kegg-{element_type}-txt file at '{path}' as the column of slice weights. "
                    f"That is a key column, which cannot hold weights. Name a column of numbers "
                    f"instead."
                )
            if slice_weights not in original_columns:
                self.progress.end()
                raise ConfigError(
                    f"'{weights_flag}' names the '{slice_weights}' column as the column of slice "
                    f"weights, but the kegg-{element_type}-txt file at '{path}' has no such "
                    f"column. Its columns are: {', '.join(original_columns)}."
                )

        # A border column is named by an option, as the column of weights is. It must be a column of
        # the file. It cannot be a key column, and it cannot hold the weights too.
        border_flag = f'--{element_type}-border-column'
        if border_column is not None:
            if border_column in ('accession', 'gene_id', 'sample'):
                self.progress.end()
                raise ConfigError(
                    f"'{border_flag}' names the '{border_column}' column of the "
                    f"kegg-{element_type}-txt file at '{path}' as the border column. That is a key "
                    f"column, which cannot hold border values. Name a column of values instead."
                )
            if border_column not in original_columns:
                self.progress.end()
                raise ConfigError(
                    f"'{border_flag}' names the '{border_column}' column as the border column, but "
                    f"the kegg-{element_type}-txt file at '{path}' has no such column. Its columns "
                    f"are: {', '.join(original_columns)}."
                )
            if border_column == slice_weights:
                self.progress.end()
                raise ConfigError(
                    f"'{border_flag}' and '{weights_flag}' both name the '{border_column}' column "
                    f"of the kegg-{element_type}-txt file at '{path}'. A column can either color "
                    f"borders or hold the weights of slices, not both. Please name another column "
                    f"for one of them."
                )

        # A bar column is named by an option, as a border column is. It must be a column of the
        # file. It cannot be a key column, and it cannot hold the weights or the border values too.
        bar_flag = f'--{element_type}-bar-column'
        if bar_column is not None:
            if bar_column in ('accession', 'gene_id', 'sample'):
                self.progress.end()
                raise ConfigError(
                    f"'{bar_flag}' names the '{bar_column}' column of the kegg-{element_type}-txt "
                    f"file at '{path}' as the bar column. That is a key column, which cannot hold "
                    f"bar values. Name a column of numbers instead."
                )
            if bar_column not in original_columns:
                self.progress.end()
                raise ConfigError(
                    f"'{bar_flag}' names the '{bar_column}' column as the bar column, but the "
                    f"kegg-{element_type}-txt file at '{path}' has no such column. Its columns "
                    f"are: {', '.join(original_columns)}."
                )
            for other_flag, other_column, holds in (
                (weights_flag, slice_weights, 'hold the weights of slices'),
                (border_flag, border_column, 'color borders')
            ):
                if bar_column == other_column:
                    self.progress.end()
                    raise ConfigError(
                        f"'{bar_flag}' and '{other_flag}' both name the '{bar_column}' column of "
                        f"the kegg-{element_type}-txt file at '{path}'. A column can either set "
                        f"the lengths of bars or {holds}, not both. Please name another column for "
                        f"one of them."
                    )

        # Auto-detect the value column: the single column that is not a key column ('accession',
        # 'gene_id' for reactions, or 'sample'). None means presence coloring; one means
        # quantitative coloring by that column (its header labels the colorbar); two or more is
        # ambiguous. A column of slice weights, a border column, and a bar column are left out of
        # the search.
        key_columns = {
            'accession', 'gene_id', 'sample', slice_weights, border_column, bar_column
        } - {None}
        value_columns = [column for column in original_columns if column not in key_columns]
        if len(value_columns) > 1:
            self.progress.end()
            # A column of weights or a border column not named by its option is a likely cause.
            hint = (
                f" If one of these columns holds the weights of slices, name it with "
                f"'{weights_flag}'." if element_type == 'reaction' and slice_weights is None else ""
            )
            if border_column is None:
                hint += f" If one of them holds border values, name it with '{border_flag}'."
            if element_type == 'reaction' and bar_column is None:
                hint += f" If one of them sets the lengths of bars, name it with '{bar_flag}'."
            raise ConfigError(
                f"The kegg-{element_type}-txt file at '{path}' has more than one candidate value "
                f"column, so it is ambiguous which to color elements by: "
                f"{', '.join(value_columns)}. A file may have at most one value column (any column "
                f"that is not a key column): include zero to color by presence/absence, or one to "
                f"color by that value.{hint}"
            )
        value_column = value_columns[0] if value_columns else None
        # Slices are colored by their values, so a file with weights needs a value column too. The
        # weights can instead size the slices of the rings of a border column.
        sliced_border = border_column is not None and border_slice_by is not None
        if slice_weights is not None and value_column is None and not sliced_border:
            others = 'its key columns, its border column,' if border_column is not None else (
                'its key columns'
            )
            remedy = f"Add a value column, or drop '{weights_flag}'" if border_column is None else (
                f"Add a value column, slice only the rings of the border column with "
                f"'--{element_type}-border-slice-by' in place of '--{element_type}-slice-by', or "
                f"drop '{weights_flag}'"
            )
            self.progress.end()
            raise ConfigError(
                f"'{weights_flag}' names the '{slice_weights}' column of the "
                f"kegg-{element_type}-txt file at '{path}' as the column of slice weights. The "
                f"file has no column besides {others} and this one, so it has no values to color "
                f"the slices by. {remedy} to color equal slices by '{slice_weights}'."
            )

        # Slices of genes are named by the 'gene_id' column, so every row needs a gene ID. Blanks
        # around an ID are dropped. Without slices of genes, the column is read as it is. A file
        # without a value column is not checked here. Its layer shows presence, and slicing it is
        # refused for that reason first ('_build_slice_model'). Rings sliced by gene need gene IDs
        # too, even in a file without a value column. The messages name the option of the rings only
        # where the fill is not sliced by gene.
        fill_genes = slice_by == 'gene' and value_column is not None
        ring_genes = sliced_border and border_slice_by == 'gene'
        if fill_genes or ring_genes:
            if fill_genes:
                gene_flag = f'--{element_type}-slice-by'
                gives = "gives each gene a slice of its own"
            else:
                gene_flag = f'--{element_type}-border-slice-by'
                gives = "gives each gene its own slice of a ring"
            if 'gene_id' not in df.columns:
                self.progress.end()
                raise ConfigError(
                    f"'{gene_flag} gene' {gives}, but the kegg-{element_type}-txt file at '{path}' "
                    f"has no 'gene_id' column. Its columns are: {', '.join(original_columns)}. Add "
                    f"a 'gene_id' column that names the gene of each row, or give '{gene_flag} "
                    f"accession' instead."
                )
            gene_id = df['gene_id'].fillna('').astype(str).str.strip()
            blank = gene_id == ''
            if blank.any():
                bad = sorted(set(df.loc[blank, '__accession']))[:5]
                self.progress.end()
                raise ConfigError(
                    f"'{gene_flag} gene' {gives}, so every row of the kegg-{element_type}-txt file "
                    f"at '{path}' needs a gene ID in its 'gene_id' column. {int(blank.sum())} "
                    f"{'rows have' if blank.sum() > 1 else 'row has'} none. The accessions with "
                    f"such rows include these: {', '.join(bad)}."
                )
            df = df.assign(__gene_id=gene_id)

        # Categorical values are names, so each row needs one. Blanks around a name are dropped.
        if value_type == 'categorical':
            if value_column is None:
                self.progress.end()
                raise ConfigError(
                    f"'--{element_type}-value-type categorical' says that the value column of the "
                    f"kegg-{element_type}-txt file at '{path}' holds categorical values. The file "
                    f"has no value column. Its columns are: {', '.join(original_columns)}. Add a "
                    f"column of values, or drop the option to color the layer by presence."
                )
            value = df[value_column].fillna('').astype(str).str.strip()
            blank = value == ''
            if blank.any():
                bad = sorted(set(df.loc[blank, '__accession']))[:5]
                self.progress.end()
                raise ConfigError(
                    f"Every row of the kegg-{element_type}-txt file at '{path}' needs a "
                    f"categorical value in its '{value_column}' column. {int(blank.sum())} "
                    f"{'rows have' if blank.sum() > 1 else 'row has'} none. The accessions with "
                    f"such rows include these: {', '.join(bad)}."
                )
            df = df.assign(__value=value)

        # Every row needs a weight, a finite number of 0 or more. Every column is read as text, so
        # the weights are converted here. A blank or a word such as 'NA' is not a number.
        if slice_weights is not None:
            weights = pd.to_numeric(
                df[slice_weights].fillna('').astype(str).str.strip(), errors='coerce'
            )
            invalid = ~(np.isfinite(weights) & (weights >= 0))
            if invalid.any():
                bad = sorted(set(df.loc[invalid, '__accession']))[:5]
                self.progress.end()
                raise ConfigError(
                    f"Every row of the kegg-{element_type}-txt file at '{path}' needs a slice "
                    f"weight in its '{slice_weights}' column ('{weights_flag}'). A weight is a "
                    f"finite number of 0 or more. {int(invalid.sum())} "
                    f"{'rows have' if invalid.sum() > 1 else 'row has'} a blank, a negative "
                    f"number, or something other than a finite number. Their accessions include "
                    f"these: {', '.join(bad)}."
                )
            df = df.assign(__slice_weight=weights.astype(float))

        # A row can leave its border cell empty. It then has no border value, and its other columns
        # are read as usual. An empty cell is read as missing, or as empty text in a file of
        # categorical values. A numeric border value must be a finite number. A border column with
        # no value at all colors nothing, so it is refused.
        border_empty = 0
        if border_column is not None:
            border_text = df[border_column].fillna('').astype(str).str.strip()
            has_border = border_text != ''
            if not has_border.any():
                self.progress.end()
                raise ConfigError(
                    f"The border column '{border_column}' of the kegg-{element_type}-txt file at "
                    f"'{path}' ('{border_flag}') has no values. Every cell of it is empty, so it "
                    f"would color no border. Name a column that holds border values."
                )
            border_empty = int((~has_border).sum())
            if border_value_type == 'categorical':
                border_values = border_text.where(has_border, np.nan)
            else:
                border_values = pd.to_numeric(border_text.where(has_border, ''), errors='coerce')
                invalid = has_border & ~np.isfinite(border_values)
                if invalid.any():
                    bad = sorted(set(df.loc[invalid, '__accession']))[:5]
                    self.progress.end()
                    raise ConfigError(
                        f"The border column '{border_column}' of the kegg-{element_type}-txt file "
                        f"at '{path}' holds numbers, so each of its values must be a finite "
                        f"number. {int(invalid.sum())} "
                        f"{'rows have' if invalid.sum() > 1 else 'row has'} something else, such "
                        f"as text or 'NA'. Their accessions include these: {', '.join(bad)}. A row "
                        f"without a border value can leave its cell empty. If the column holds "
                        f"categorical values, such as names, read it as such with "
                        f"'--{element_type}-border-value-type categorical'."
                    )
            df = df.assign(__border_value=border_values)

        # A row can leave its bar cell empty. It then has no bar value, and its other columns are
        # read as usual. A bar value must be a finite number. A negative value grows its bar the
        # other way from 0. A bar column with no value at all sets no bar, so it is refused.
        bar_empty = 0
        if bar_column is not None:
            bar_text = df[bar_column].fillna('').astype(str).str.strip()
            has_bar = bar_text != ''
            if not has_bar.any():
                self.progress.end()
                raise ConfigError(
                    f"The bar column '{bar_column}' of the kegg-{element_type}-txt file at "
                    f"'{path}' ('{bar_flag}') has no values. Every cell of it is empty, so it "
                    f"would set no bar. Name a column that holds bar values."
                )
            bar_empty = int((~has_bar).sum())
            bar_values = pd.to_numeric(bar_text.where(has_bar, ''), errors='coerce')
            rows = has_bar & ~np.isfinite(bar_values)
            if rows.any():
                bad = sorted(set(df.loc[rows, '__accession']))[:5]
                self.progress.end()
                raise ConfigError(
                    f"The bar column '{bar_column}' of the kegg-{element_type}-txt file at "
                    f"'{path}' sets the lengths of bars, so each of its values must be a finite "
                    f"number. {int(rows.sum())} {'rows have' if rows.sum() > 1 else 'row has'} "
                    f"something other than a finite number, such as text or 'NA'. Their accessions "
                    f"include these: {', '.join(bad)}. A row without a bar value can leave its "
                    f"cell empty."
                )
            df = df.assign(__bar_value=bar_values)

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
        # ordinary case. Slices of genes read the gene IDs without the spaces around them. The rule
        # also applies to a file with a border or bar column, even one that colors its fills by
        # presence.
        if value_column is not None or border_column is not None or bar_column is not None:
            gene_column = '__gene_id' if '__gene_id' in df.columns else 'gene_id'
            key_columns = [
                column for column in ('__accession', gene_column, '__sample')
                if column in df.columns
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
                # Rows that differ only in the spaces around their gene IDs look different in the
                # file, so the message says why they repeat.
                stripped = (
                    f" With '{gene_flag} gene', the spaces around gene IDs are dropped before rows "
                    f"are compared." if gene_column == '__gene_id' else ''
                )
                if value_column is None:
                    # The values of such a file are those of its border and bar columns.
                    held = [
                        (noun, column) for noun, column in (
                            ('border value', border_column), ('bar value', bar_column)
                        ) if column is not None
                    ]
                    can_have = ' and '.join(
                        f"a {noun} in the '{column}' column" for noun, column in held
                    )
                    raise ConfigError(
                        f"Each row of the kegg-{element_type}-txt file at '{path}' can have "
                        f"{can_have}, so every row must describe a different thing. "
                        f"{int(duplicated.sum())} rows repeat a combination of {described}. The "
                        f"first repeated "
                        f"{'combinations are' if len(key_columns) > 1 else 'accessions are'}: "
                        f"{shown}. Each of these can have only one "
                        f"{' and one '.join(noun for noun, _ in held)}. Please keep one row for "
                        f"each.{stripped}"
                    )
                if value_type == 'categorical':
                    raise ConfigError(
                        f"Every row of the kegg-{element_type}-txt file at '{path}' must describe "
                        f"a different thing, but {int(duplicated.sum())} rows repeat a combination "
                        f"of {described}. The first repeated "
                        f"{'combinations are' if len(key_columns) > 1 else 'accessions are'}: "
                        f"{shown}. Each of these can have only one categorical value in the "
                        f"'{value_column}' column. Please keep one row for each.{stripped}"
                    )
                raise ConfigError(
                    f"Every row of the kegg-{element_type}-txt file at '{path}' must describe a "
                    f"different thing, since each carries a value in the '{value_column}' column, but "
                    f"{int(duplicated.sum())} rows repeat a combination of {described}. The first "
                    f"repeated {'combinations are' if len(key_columns) > 1 else 'accessions are'}: "
                    f"{shown}. Anvi'o will not guess how repeated rows should be combined -- averaging "
                    f"replicate measurements and adding up separate contributions, such as the ions of "
                    f"one metabolite, are both reasonable and give different answers. Please combine "
                    f"them in the way that suits your data before drawing.{stripped}"
                )

        self.progress.end()

        return {
            'element_type': element_type,
            'reaction_source': reaction_source,
            'df': df,
            'value_column': value_column,
            'value_type': value_type,
            'sample_names': sample_names,
            'slice_weights': slice_weights,
            'border_column': border_column,
            'border_value_type': border_value_type,
            'border_empty': border_empty,
            'bar_column': bar_column,
            'bar_empty': bar_empty
        }

    @staticmethod
    def _channel_view(data: dict, role: Literal['border', 'bar']) -> dict:
        """
        Make a reader result whose values are those of the border or bar column of a file.

        The model of the border or the bars is built from this view by the same method as the model
        of the fill ('_build_txt_model'). The view holds only the rows with a value in the column.
        In these rows, the column holds that value, and so does '__value', as text for categorical
        values. The view names the column as its value column. A bar column holds numbers only.

        Parameters
        ==========
        data : dict
            The '_read_element_txt' result of a file with a border or bar column.

        role : Literal['border', 'bar']
            The channel of the column ('TXT_CHANNEL_ROLES').

        Returns
        =======
        dict
            A reader result of the same form, with the values of the column in place of the fill
            values. The names of the samples are those of the whole file.
        """
        column = data[f'{role}_column']
        df = data['df']
        rows = df[df[f'__{role}_value'].notna()]
        values = rows[f'__{role}_value']
        rows = rows.assign(**{column: values, '__value': values.astype(str)})
        return {
            **data,
            'df': rows,
            'value_column': column,
            'value_type': data['border_value_type'] if role == 'border' else 'numeric'
        }

    @staticmethod
    def _resolve_aggregation(
        aggregation: str,
        flag: str,
        channel: dict,
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

        channel : dict
            The channel of the values ('_txt_channel'). It names the column and the flag of the
            value type in an error message.

        period : Union[float, None], None
            The period after which the layer's values repeat, or None if they do not repeat.

        period_flag : Union[str, None], None
            The command-line flag that gives the layer's period, used in an error message.

        Returns
        =======
        Callable
            Reduces a sequence of values to a single value.
        """
        # A categorical rule chooses among names rather than reducing numbers. A layer of
        # categorical values resolves its rules elsewhere ('_resolve_categorical_rule'). A bar takes
        # numbers only, so its message names no value type.
        if aggregation in CATEGORICAL_RULES:
            hint = (
                f" To read the {channel['column']} as categorical values, give "
                f"'--{channel['flag_stem']}-value-type categorical'."
            ) if channel['categorical'] else ""
            raise ConfigError(
                f"'{flag}' was given as '{aggregation}', which chooses one categorical value among "
                f"several. The values here are numbers, so they are reduced with an aggregation, "
                f"such as 'sum' or 'mean'.{hint}"
            )
        if aggregation in CIRCULAR_AGGREGATIONS + CIRCULAR_SPREADS:
            if not channel['period']:
                raise ConfigError(
                    f"'{flag}' was given as '{aggregation}', which is for values that repeat after "
                    f"a period, such as clock times. The {channel['values']} take no period, since "
                    f"the length of a bar starts at 0. Give an aggregation, such as 'mean', "
                    f"instead."
                )
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
            # A border takes no presence summaries ('_resolve_summary'), so its message names none.
            presence_note = (
                f" Note that the sample and group summary options additionally take "
                f"{SUMMARY_PRESENCE_PHRASE} to summarize presence rather than value."
            ) if channel['presence'] else ""
            raise ConfigError(
                f"'{flag}' was given as '{aggregation}', which either does not reduce values to a "
                f"single number or does not reduce them the same way at every level. Proven "
                f"acceptable names are "
                f"{', '.join(repr(name) for name in AGGREGATION_FUNCTIONS)}; any other pandas "
                f"aggregation that reduces a series to one number, such as 'var' or 'sem', also "
                f"works. Names that transform rather than reduce, like 'cumsum', cannot be used; "
                f"neither can names offered only by a grouping, like 'first', nor names meaning "
                f"different things for a list of values and for a table, like 'idxmax'."
                f"{presence_note}"
            )
        return aggregate

    @staticmethod
    def _circular_mean(
        values: Iterable[float],
        period: float,
        weights: Union[Iterable[float], None] = None
    ) -> float:
        """
        Average values that repeat after a period, such as clock times, on a circle.

        Each value is a point on a circle whose circumference is the period. 0 and the period are
        the same point. The mean is the direction of the average of these points. With a period of
        24, the mean of 23.5 and 0.5 is 0, not 12. The mean lies in [0, period).

        Values that cancel out have no mean direction. Examples are 6 and 18, or 0, 8 and 16, with a
        period of 24. The mean is then undefined, and the result is NaN. The value is then dropped
        as any undefined aggregation is ('_finite_values', '_reduce_element_value').

        Parameters
        ==========
        values : Iterable[float]
            The values to average. There is at least one. They can be a pandas Series, as a groupby
            passes, so they are read by position.

        period : float
            The positive period after which the values repeat.

        weights : Union[Iterable[float], None], None
            A weight for each value, such as the shares of the slices of a map element
            ('_slice_standin_value'). The average of the points is then weighted. None weighs every
            value the same.

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
            if weights is None:
                sine = np.mean(np.sin(angles))
                cosine = np.mean(np.cos(angles))
            else:
                sine = np.average(np.sin(angles), weights=weights)
                cosine = np.average(np.cos(angles), weights=weights)
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
        undefined aggregation is ('_finite_values', '_reduce_element_value').

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

        The function takes one map element's values in the samples or groups that have one, in
        order. It returns a new value for each of them, in the same order. A value it cannot define
        is returned as 'nan' or as an infinity. Either leaves that element uncolored on the map of
        that sample or group, exactly as an undefined aggregation leaves it uncolored
        ('_reduce_element_value').

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
        channel: dict,
        period: Union[float, None] = None,
        period_flag: Union[str, None] = None
    ) -> Tuple[Callable, str, bool]:
        """
        Resolve a normalization name into a function rescaling an element's sample or group values.

        The recommended names are the keys of 'ELEMENT_NORMALIZATIONS'. Any other name is taken to
        be a pandas Series method that transforms the values it is given into one new value each,
        such as 'abs', and is probed here, before any data is touched, so that a name pandas does
        not recognize, or one that does something other than transform, is reported as a
        configuration error rather than failing deep in a drawing loop. This mirrors
        '_resolve_aggregation' at the next level up. The two are opposites. An aggregation must
        reduce a set of values to one number, so it refuses a name that transforms. A normalization
        must return one value per sample or group, so it refuses a name that reduces.

        Probes ask four things of a name. It must return a Series of one value per sample or group,
        with the same index in the same order. A name that reorders or drops samples or groups would
        give a sample another sample's value without any error. It must return numbers, which alone
        can go on a color scale. It must give the same answer twice for the same values, since the
        range of the scale and the colors on the maps are worked out in separate passes. And it must
        not depend on the order the samples or groups come in. That order comes only from how the
        samples happen to be named and listed. It says nothing about the data. 'cumsum' and 'ffill'
        answer differently for a sample depending on what precedes it, so they are refused. What a
        name MAKES of the values is its own business, as it is for an aggregation: a name that
        transforms each value without reference to the others, such as 'abs', is accepted, and the
        colorbar says which name was used.

        Parameters
        ==========
        normalization : str
            The requested normalization name.

        flag : str
            The command-line flag the name comes from, used in error messages.

        channel : dict
            The channel of the values ('_txt_channel'). It names the flag of the sample summary in
            an error message.

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
                f"such as 'mean', belongs to '--{channel['flag_stem']}-sample-summary' instead, "
                f"which is what summarizes samples."
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

        # The probes label the values with sample names rather than positions. The names are not in
        # alphabetical order, so a name that sorts by either the values or the index ('sort_values',
        # 'sort_index') changes their order. Ties in the values unmask a name that quietly drops
        # repeated values ('drop_duplicates'). A missing value reveals a name that drops a sample
        # ('dropna') or fills its value from its neighbors ('ffill'). Neither changes a complete set
        # of values. The last two probes are a pair of values and a single value. Two samples or
        # groups are the fewest that a normalization across them can say anything about. One value
        # is what an element found in a single sample gives it. A name that collapses to a scalar
        # ('squeeze') would otherwise fail there, in the middle of drawing.
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
                refuse("pandas could not apply to the values of a set of samples or groups")

            if isinstance(result, (float, int, np.number)) or np.isscalar(result):
                refuse(
                    "reduces a set of values to a single number rather than giving each of them a "
                    "new value"
                )

            if not isinstance(result, pd.Series):
                refuse("does not give each sample or group a value of its own")

            if list(result.index) != index:
                refuse(
                    "does not return the values of the same samples or groups in the same order, "
                    "so anvi'o cannot tell which of them each value belongs to"
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

        # Reversing the samples, rotating them by one, and swapping the first two catch a dependence
        # on their order. Any single rearrangement alone might give the same answer and miss it. The
        # tied values are probed as well. A name can answer the same for values that all differ and
        # still depend on which of two equal values comes first ('duplicated'). The missing value is
        # probed too. A name that fills a sample's value from its neighbors only shows that it does
        # so when there is something to fill.
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
                    refuse("pandas could not apply to the values of a set of samples or groups")

                if not np.allclose(
                    given, as_floats(rearranged), rtol=1e-12, atol=0, equal_nan=True
                ):
                    refuse(
                        "answers differently when the samples or groups are given in a different "
                        "order. Their order says nothing about the data. It comes only from how "
                        "the samples happen to be named, so the results of this normalization "
                        "would depend on those names"
                    )

        def normalize(values) -> np.ndarray:
            # The samples or groups are numbered here rather than named. A name that passed the
            # probes above answers the same whatever their labels are. Numbering is cheaper.
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

    def _warn_undefined_values(
        self,
        undefined: Set[str],
        aggregation: str,
        path: str,
        channel: dict
    ) -> None:
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

        channel : dict
            The channel of the values ('_txt_channel'), for the message.
        """
        if not undefined:
            return
        examples = ', '.join(sorted(undefined)[:5])
        self.run.warning(
            f"Reducing the {channel['values']} of the text file at '{path}' with '{aggregation}' "
            f"was undefined for {len(undefined)} accession(s), including these: {examples}. "
            f"{Mapper._undefined_cause(aggregation)} These accessions are treated as having no "
            f"value, so the map elements that depend on them {channel['consequence_plural']}."
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
        if aggregation in CIRCULAR_SPREADS:
            return (
                f"An angular deviation needs at least {MIN_ANGULAR_DEVIATION_VALUES} values. The "
                f"spread of two values is only the gap between them."
            )
        return (
            f"This happens when {noun} needs more values than are available, as the standard "
            f"deviation does for a single value."
        )

    def _warn_undefined_summaries(
        self,
        layer: dict,
        draw_unified_maps: bool,
        draw_individual_maps: bool
    ) -> None:
        """
        Report the map elements of one layer whose sample or group summary was undefined.

        '_summarize_entry_values' collects these elements while the range of values is found. A
        summary is reported only where it changes a map that is drawn. Ungrouped, the sample summary
        colors the 'unified' map. Grouped, the sample summary colors the map of each group. The
        group summary then pools the groups' values on the 'unified' map, so a sample summary that
        is undefined in a group also changes the 'unified' map where the group summary pools values.
        A layer sliced by accession is summarized one accession at a time, so its report counts
        accessions. A layer sliced by gene is summarized one slice of a gene at a time. Its report
        counts these slices, and a gene with the same rows in two elements counts once. A border
        with sliced rings is reported in the same way, in terms of its rings.

        Parameters
        ==========
        layer : dict
            A 'quantitative' layer model with samples.

        draw_unified_maps : bool
            True if the 'unified' map is drawn.

        draw_individual_maps : bool
            True if the maps of the individual samples or groups are drawn.
        """
        channel = layer['channel']
        flag_stem = channel['flag_stem']
        grouped = layer['group_samples'] is not None
        # A sliced layer summarizes each accession or gene of an element on its own
        # ('_slice_result'). A slice whose summary is undefined is not drawn. The other slices of
        # the element then take its whole width or length. A slice of a ring without a value takes
        # no bands where the rings are aligned with the fill. The other slices of a ring sliced on
        # its own take its length.
        sliced = layer.get('slice_by') is not None
        slice_noun, value_noun = 'slice', 'value'
        group_effect = "It is not drawn on that group's map."
        unified_effect = "These slices are not drawn on the 'unified' map."
        others = " The other slices of an element take its whole width or length."
        if channel['role'] == 'border' and layer.get('slice_layout') == 'aligned':
            slice_noun, value_noun = 'slice of a ring', 'border value'
            group_effect = (
                "On that group's map, its column takes no bands, or it does not count toward a "
                "ring drawn whole."
            )
            unified_effect = (
                "On the 'unified' map, their columns take no bands, or they do not count toward "
                "rings drawn whole."
            )
            others = ""
        elif channel['role'] == 'border':
            slice_noun, value_noun = 'slice of a ring', 'border value'
            unified_effect = "These slices of rings are not drawn on the 'unified' map."
            others = " The other slices of a ring take its length."
        reports = []
        if grouped:
            consequences = []
            if draw_individual_maps:
                consequences.append(
                    group_effect if sliced
                    else f"It {channel['consequence']} on that group's map."
                )
            if draw_unified_maps and layer['group_aggregate'] is not None:
                consequences.append(
                    f"The 'unified' map pools only the groups where the {slice_noun} has a "
                    f"{value_noun}." if sliced else (
                        f"The 'unified' map pools only the groups where the element has a value. "
                        f"An element with no value in any group {channel['consequence']} there."
                    )
                )
            if consequences:
                reports.append((
                    layer['_undefined_individual'], f'--{flag_stem}-sample-summary',
                    layer['sample_summary'], ' in at least one group',
                    ' '.join(
                        [
                            f"A {slice_noun} has no {value_noun} in a group where its summary is "
                            f"undefined." if sliced else
                            "An element has no value in a group where its summary is undefined."
                        ]
                        + consequences
                    )
                ))
        if draw_unified_maps:
            reports.append((
                layer['_undefined_unified'],
                f"--{flag_stem}-{'group' if grouped else 'sample'}-summary",
                layer['group_summary'] if grouped else layer['sample_summary'], '',
                unified_effect if sliced
                else f"These elements {channel['consequence_plural']} on the 'unified' map."
            ))
        for undefined, flag, summary, where, consequence in reports:
            if not undefined:
                continue
            if layer.get('slice_by') == 'gene':
                # Each key names the rows of one gene, as the gene ID and the accession joined by a
                # tab ('_gene_members').
                genes = self._sort_values({
                    member.rsplit('\t', 1)[0] for key in undefined for member in key
                })
                self.run.warning(
                    f"The summary '{flag} {summary}' was undefined{where} for {len(undefined)} "
                    f"slice(s) of genes. Their genes include these: {', '.join(genes[:5])}. "
                    f"{Mapper._undefined_cause(summary, 'a summary')} {consequence}{others}"
                )
                continue
            examples = ', '.join(sorted({
                kegg_id for key in undefined for kegg_id in key if kegg_id in layer['accessions']
            })[:5])
            if sliced:
                self.run.warning(
                    f"The summary '{flag} {summary}' was undefined{where} for {len(undefined)} "
                    f"accession(s), which include these: {examples}. "
                    f"{Mapper._undefined_cause(summary, 'a summary')} {consequence}{others}"
                )
                continue
            self.run.warning(
                f"The summary '{flag} {summary}' was undefined{where} for {len(undefined)} map "
                f"element(s). Elements with the same accessions count once. Their accessions in "
                f"the file include these: {examples}. "
                f"{Mapper._undefined_cause(summary, 'a summary')} {consequence}"
            )

    def _report_categorical_choices(
        self,
        layers: List[dict],
        draw_unified_maps: bool,
        drawn_individuals: List[str],
        draw_maps_lacking_data: bool = False,
        draw_grid: bool = False
    ) -> None:
        """
        Report how often a rule chose one categorical value among different values.

        Values meet among
        - the genes of an accession
        - the accessions of a map element
        - the samples on the 'unified' map or in a group
        - the groups on the 'unified' map
        - the reactions that touch a compound.
        Where they differ, a rule chooses one of them ('_choose_categorical_value'). The choices
        among genes are counted for every gene of the file, when the file is read
        ('_select_accession_values'). The choices among accessions, samples, and groups are counted
        where a drawn map uses them ('_summarize_entry_individuals'). The choices for compounds are
        counted once per compound and set of votes ('_categorical_compound_chooser'). They are
        reported here once per layer and level. Each count says how the choices were made, so the
        effect of the rule can be seen: by plurality, by majority, by the order of the values, at
        random, or without a value. A choice without a value can leave an element uncolored, so it
        looks absent. With slices, the accession rule chooses only the stand-in values of elements
        ('_slicer'), and a choice without a value can leave a slice undrawn. With slices of genes,
        it also chooses the value of a gene among its rows for the accessions of an element. The
        genes of an accession are then never reduced, so the choices among them are not reported.
        The same holds for the slices of rings of a border ('_ring_result'). These choices are also
        reported in a warning that says what they change in this run. A layer whose values color no
        drawn map is not reported.

        Parameters
        ==========
        layers : List[dict]
            The layer models. Only layers of categorical values are reported.

        draw_unified_maps : bool
            True if the 'unified' map is drawn.

        drawn_individuals : List[str]
            The samples or groups drawn on their own maps.

        draw_maps_lacking_data : bool, False
            True if maps with nothing to color are drawn.

        draw_grid : bool, False
            True if map grids are drawn.
        """
        ways = (
            ('plurality', 'by plurality'), ('majority', 'by majority'), ('order', 'by value order'),
            ('random', 'at random'), ('none', 'without a value')
        )
        draw_individual_maps = bool(drawn_individuals)
        # A map with nothing to color is not drawn unless such maps are asked for. A grid shows a
        # blank panel in place of a missing 'unified' map ('_draw_map_grids').
        unified_skip = ''
        if not draw_maps_lacking_data:
            unified_skip = " A 'unified' map on which no element has a value is not drawn."
            if draw_grid:
                unified_skip += " The map grid shows a blank panel in its place."
        for layer in layers:
            if layer.get('value_type') != 'categorical':
                continue
            if not (
                (draw_individual_maps and layer['individual_mode'] == 'categorical')
                or (draw_unified_maps and layer['unified_mode'] == 'categorical')
            ):
                continue
            channel = layer['channel']
            has_samples = layer['sample_values'] is not None
            grouped = layer['group_samples'] is not None
            in_sample = ' in a sample' if has_samples else ''

            # What a choice without a value changes depends on what the run draws, and on whether a
            # rule at the next level counts its votes. The 'unified' map has such a rule only where
            # one was given for it.
            unified_votes = draw_unified_maps and (
                layer['group_rule' if grouped else 'sample_rule'] is not None
            )
            # A choice that the 'unified' map uses is counted even where the map of its sample or
            # group is not drawn.
            if_drawn = ''
            if has_samples and unified_votes and draw_individual_maps:
                individuals = layer['group_samples'] if grouped else layer['sample_values']
                if not set(individuals) <= set(drawn_individuals):
                    if_drawn = ' if that map is drawn'
            # An accession without a value can leave an element with no votes. The element then has
            # no value either.
            gene_consequence = (
                f"Such an accession has no value{' in that sample' if has_samples else ''}. It "
                f"does not vote in its map elements."
            )
            if not has_samples:
                gene_consequence += (
                    " An element is left uncolored on its map if none of its accessions has a "
                    "value. It looks like an element the file lacks."
                )
            elif draw_individual_maps and not grouped:
                gene_consequence += (
                    " An element is left uncolored on a sample's map if none of its accessions has "
                    "a value in that sample. It looks like an element the file lacks."
                )
            elif draw_individual_maps:
                gene_consequence += (
                    " An element is left uncolored on a group's map if none of its accessions has "
                    "a value in any of the group's samples. It looks like an element the group's "
                    "samples lack."
                )
            if has_samples and unified_votes:
                gene_consequence += (
                    " An element is left uncolored on the 'unified' map if none of its accessions "
                    "has a value in any sample."
                )
            # With slices, a choice without a value leaves a slice undrawn, or an element without a
            # stand-in value ('_slicer'). Such an element casts no vote for its compounds, and it is
            # left uncolored only where it would be drawn whole.
            sliced = layer.get('slice_by') is not None
            standin_effect = (
                "casts no vote for the color of its compounds, and it is left uncolored if it has "
                "more slices than fit"
            )
            # With slices of genes, the rule also chooses the value of a gene among its rows for the
            # accessions of an element. A choice without a value leaves the gene without a slice.
            gene_sliced = layer.get('slice_by') == 'gene'
            if gene_sliced:
                element_consequence = (
                    "With slices of genes, this rule chooses the value of a gene with rows for two "
                    "or more accessions of an element, and the stand-in value of an element."
                )
                if not has_samples:
                    element_consequence += (
                        f" Such a gene gets no slice. Such an element {standin_effect}."
                    )
                elif grouped:
                    element_consequence += (
                        " Such a gene or element has no value in that sample. It does not vote "
                        "among the samples of its group."
                    )
                    if draw_individual_maps or unified_votes:
                        element_consequence += (
                            f" A gene with no value on a drawn map gets no slice there. An element "
                            f"with no stand-in value on a drawn map {standin_effect}."
                        )
                else:
                    element_consequence += " Such a gene or element has no value in that sample."
                    if draw_individual_maps:
                        element_consequence += (
                            f" On the map of that sample{if_drawn}, such a gene gets no slice. "
                            f"Such an element there {standin_effect}."
                        )
                    if unified_votes:
                        element_consequence += (
                            " Such a gene or element does not vote across the samples for the "
                            "'unified' map."
                        )
            elif sliced:
                element_consequence = (
                    "With slices, this rule chooses only the stand-in value of an element."
                )
                if not has_samples:
                    element_consequence += f" Such an element {standin_effect}."
                elif grouped:
                    element_consequence += (
                        " Such an element has no stand-in value in that sample. It does not vote "
                        "among the samples of its group."
                    )
                    if draw_individual_maps or unified_votes:
                        element_consequence += (
                            f" An element with no stand-in value on a drawn map {standin_effect}."
                        )
                else:
                    element_consequence += " Such an element has no stand-in value in that sample."
                    if draw_individual_maps:
                        element_consequence += (
                            f" On the map of that sample{if_drawn}, it {standin_effect}."
                        )
                    if unified_votes:
                        element_consequence += (
                            " It does not vote across the samples for the stand-in value on the "
                            "'unified' map."
                        )
            elif not has_samples:
                element_consequence = (
                    "Such an element is left uncolored on its map. It looks like an element the "
                    "file lacks."
                )
            elif grouped:
                element_consequence = (
                    "Such an element has no value in that sample. It does not vote among the "
                    "samples of its group."
                )
                if draw_individual_maps:
                    element_consequence += (
                        " An element is left uncolored on a group's map if it has no value in any "
                        "of the group's samples. It looks like an element the group's samples lack."
                    )
                if unified_votes:
                    element_consequence += (
                        " An element is left uncolored on the 'unified' map if it has no value in "
                        "any sample."
                    )
            else:
                element_consequence = "Such an element has no value in that sample."
                if draw_individual_maps:
                    element_consequence += (
                        f" It is left uncolored on the map of that sample{if_drawn}. It looks like "
                        f"an element the file lacks."
                    )
                if unified_votes:
                    element_consequence += (
                        " It does not vote across the samples. An element is left uncolored on the "
                        "'unified' map if it has no value in any sample."
                    )
            if channel['role'] == 'border' and sliced:
                # With sliced rings, the rules choose the border values of the slices of rings, and
                # the stand-in values of rings drawn whole ('_ring_result'). A slice without a value
                # shows no border. A ring drawn whole without a stand-in value shows none either.
                in_that_sample = ' in that sample' if has_samples else ''
                if layer.get('slice_layout') == 'aligned':
                    slice_effect = "the column of such a slice shows no bands"
                else:
                    slice_effect = (
                        "such a slice of a ring is not drawn, and the other slices of its ring "
                        "take its length"
                    )
                ring_effect = "shows no border where its ring is drawn whole"
                gene_consequence = (
                    f"Such an accession has no border value{in_that_sample}. Its slice of a ring "
                    f"shows no border, and it does not vote for the stand-in value of a ring."
                )
                if gene_sliced:
                    element_consequence = (
                        "With rings sliced by gene, this rule chooses the border value of a gene "
                        "with rows for two or more accessions of an element, and the stand-in "
                        "value of a ring."
                    )
                    subject = 'a gene or element'
                    gene_effect = "the slice of such a gene shows no border"
                else:
                    element_consequence = (
                        "With sliced rings, this rule chooses only the stand-in value of a ring."
                    )
                    subject = 'an element'
                    gene_effect = None
                if not has_samples:
                    element_consequence += (
                        f" {gene_effect.capitalize()}." if gene_effect else ""
                    ) + f" Such an element {ring_effect}."
                elif grouped:
                    element_consequence += (
                        f" Such {subject} has no border value in that sample. It does not vote "
                        f"among the samples of its group."
                    )
                else:
                    element_consequence += f" Such {subject} has no border value in that sample."
                    if draw_individual_maps:
                        element_consequence += (
                            f" On the map of that sample{if_drawn}, "
                            f"{gene_effect + ', and such an element' if gene_effect else 'it'} "
                            f"{ring_effect}."
                        )
                    if unified_votes:
                        element_consequence += " It does not vote across the samples."
                group_consequence = (
                    "Such a slice of a ring or map element has no border value in that group."
                )
                if draw_individual_maps:
                    group_consequence += (
                        f" On the map of that group{if_drawn}, {slice_effect}. Such an element "
                        f"{ring_effect}."
                    )
                if unified_votes:
                    group_consequence += " It does not vote across the groups."
                unified_consequence = f"On that map, {slice_effect}. Such an element {ring_effect}."
                pooled_what = 'slice(s) of rings or map element(s)'
            elif channel['role'] == 'border':
                # A border shows values only, so an element without a border value shows no border.
                # Its fill does not change.
                in_that_sample = ' in that sample' if has_samples else ''
                gene_consequence = (
                    f"Such an accession has no border value{in_that_sample}. It does not vote in "
                    f"its map elements. An element without votes shows no border."
                )
                if not has_samples:
                    element_consequence = "Such an element shows no border on its map."
                else:
                    element_consequence = "Such an element has no border value in that sample."
                    if grouped:
                        element_consequence += " It does not vote among the samples of its group."
                    elif draw_individual_maps:
                        element_consequence += (
                            f" It shows no border on the map of that sample{if_drawn}."
                        )
                    if not grouped and unified_votes:
                        element_consequence += " It does not vote across the samples."
                group_consequence = "Such an element has no border value in that group."
                if draw_individual_maps:
                    group_consequence += f" It shows no border on the map of that group{if_drawn}."
                if unified_votes:
                    group_consequence += " It does not vote across the groups."
                unified_consequence = "Such an element shows no border on the 'unified' map."
                pooled_what = 'map element(s)'
            elif sliced:
                # The samples and groups choose the values of slices, and the stand-in values of
                # elements.
                group_consequence = "Such a slice or map element has no value in that group."
                if draw_individual_maps:
                    group_consequence += (
                        f" On the map of that group{if_drawn}, such a slice is not drawn, and the "
                        f"other slices of its element take its whole width or length. Such an "
                        f"element {standin_effect}."
                    )
                if unified_votes:
                    group_consequence += " It does not vote across the groups."
                unified_consequence = (
                    f"Such a slice is not drawn on that map, and the other slices of its element "
                    f"take its whole width or length. Such an element {standin_effect}."
                    f"{unified_skip}"
                )
                pooled_what = 'slice(s) or map element(s)'
            else:
                group_consequence = "Such an element has no value in that group."
                if draw_individual_maps:
                    group_consequence += (
                        f" It is left uncolored on the map of that group{if_drawn}. It looks like "
                        f"an element the group's samples lack."
                    )
                if unified_votes:
                    group_consequence += (
                        " It does not vote across the groups. An element is left uncolored on the "
                        "'unified' map if no group has a value for it."
                    )
                unified_consequence = None
                pooled_what = 'map element(s)'

            for level, voters, rule_key, flag_suffix, what, consequence in (
                (
                    'gene', 'genes', 'gene_rule', 'gene-aggregation', f'accession(s){in_sample}',
                    gene_consequence
                ),
                (
                    'element',
                    f"accessions{' or genes' if gene_sliced else ''} of a map element",
                    'accession_rule', 'accession-aggregation',
                    f"{'gene(s) or ' if gene_sliced else ''}map element(s){in_sample}",
                    element_consequence
                ),
                (
                    'sample', 'samples of a group', 'sample_rule', 'sample-summary',
                    f'{pooled_what} in a group', group_consequence
                ) if grouped else (
                    'sample', 'samples', 'sample_rule', 'sample-summary',
                    f"{pooled_what} on the 'unified' map",
                    unified_consequence or (
                        f"Such an element is left uncolored on that map. It looks like an element "
                        f"that no sample has.{unified_skip}"
                    )
                ),
                (
                    'group', 'groups', 'group_rule', 'group-summary',
                    f"{pooled_what} on the 'unified' map",
                    unified_consequence or (
                        f"Such an element is left uncolored on that map. It looks like an element "
                        f"that no group has.{unified_skip}"
                    )
                ),
                (
                    'compound', 'reactions touching a compound', 'accession_rule',
                    'accession-aggregation', 'compound(s) on global or overview maps',
                    "Such a compound keeps the color of a compound that no colored reaction "
                    "touches."
                )
            ):
                # Slices of genes never reduce the genes of an accession, so the choices made among
                # them when the file was read color no map.
                if level == 'gene' and gene_sliced:
                    continue
                counts = layer['_choices'][level]
                total = sum(counts.values())
                if not total:
                    continue
                rule = layer[rule_key]
                made = ', '.join(f"{counts[how]} {phrase}" for how, phrase in ways if counts[how])
                self.run.info(
                    f"Choices of a {channel['noun']} value among {voters}",
                    f"{total} with '{rule}' ({made})"
                )
                if not counts['none']:
                    continue
                reason = (
                    "no value had more than half of the votes" if rule == 'majority' else
                    "two or more values tied for the most votes"
                )
                self.run.warning(
                    f"'--{channel['flag_stem']}-{flag_suffix} {rule}' gave no value to "
                    f"{counts['none']} {what}, since {reason}. {consequence} To give every one of "
                    f"them a value, choose 'plurality_order' or 'plurality_random' instead."
                )

    def _aggregate_accession_quantities(
        self,
        rows_df: pd.DataFrame,
        value_column: str,
        aggregation: Union[str, Callable],
        path: str,
        channel: dict,
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

        channel : dict
            The channel of the values ('_txt_channel'). It names the column and the flag of the
            value type in an error message.

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
            # A column of names is likelier a column of categorical values than a column of numbers
            # with typos in it. Pandas reads a blank cell, and text such as 'NA', as a missing
            # value, not as text. A file with a border column is read as text ('_read_element_txt').
            # Such text then names no category either.
            text = rows_df.loc[invalid & numeric_values.isna(), value_column].dropna()
            if '__border_value' in rows_df.columns:
                text = text[~text.isin(STR_NA_VALUES)]
            categorical_clause = (
                f" If the column holds categorical values, such as names, read it as such with "
                f"'--{channel['flag_stem']}-value-type categorical'."
            ) if not text.empty else ""
            self.progress.end()
            raise ConfigError(
                f"The '{value_column}' {channel['column']} of the text file at '{path}' must "
                f"contain a finite number in every row for elements to be colored quantitatively. "
                f"However, {int(invalid.sum())} {'rows have' if invalid.sum() > 1 else 'row has'} "
                f"a blank, non-numeric, or infinite value, including these accessions: "
                f"{', '.join(bad_accessions[:5])}.{categorical_clause}"
            )
        valued_df = rows_df.assign(__quantitative_value=numeric_values)
        return self._finite_values(
            valued_df.groupby('__accession')['__quantitative_value'].agg(aggregation).to_dict(),
            undefined
        )

    def _check_individual_names(self, individuals: Iterable[str], individual_noun: str) -> None:
        """
        Check that sample, contigs database, genome, or group names can serve as subdirectory names.

        Each sample, contigs database, genome, or group drawn on its own maps gets its own
        subdirectory of the output directory. Once map grids are drawn, a subdirectory that was only
        needed for a grid is deleted ('_draw_map_grids'). The names must be checked to ensure that
        they form proper paths that lie within the output directory. Names come from a text file's
        'sample' column or a groups file, so they cannot be trusted to be safe paths. For example, a
        sample named "Station 5/Depth 10" is not safe. The same names become file names where the
        maps are gathered into a directory per map ('_collate_maps_by_map'). The same check makes
        those file names safe. Only the samples, contigs databases, genomes, or groups that get
        individual maps are checked, as independent files or as part of a grid. Those that are only
        summarized on a 'unified' map are not checked.

        No name is reserved. The individual maps are drawn into '<output directory>/individual'.
        anvi'o creates that directory for this purpose alone, so a sample, contigs database, genome,
        or group may be named after anything anvi'o puts in the output directory itself. Examples
        are 'unified', 'grid', 'by_map', 'all_maps', and a BRITE category such as 'Metabolism'. A
        name from the input never shares a directory with a name that anvi'o uses.

        Parameters
        ==========
        individuals : Iterable[str]
            The names of the samples, contigs databases, genomes, or groups that would become
            subdirectory names.

        individual_noun : str
            What each name belongs to, e.g. 'sample' or 'group'. It is used in the error message.
        """
        separators = {os.sep, os.altsep} - {None}
        problems: List[str] = []
        for individual in individuals:
            if not individual or not individual.strip() or individual in ('.', '..'):
                reason = "is not a usable directory name"
            elif any(separator in individual for separator in separators):
                reason = "contains a path separator"
            elif os.path.isabs(individual) or os.path.basename(individual) != individual:
                reason = "is not a plain directory name"
            else:
                continue
            problems.append(f"'{individual}' ({reason})")

        if not problems:
            return
        plural = len(problems) > 1
        raise ConfigError(
            f"Each {individual_noun} drawn on its own maps gets its own subdirectory of the output "
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
        channel: dict,
        period: Union[float, None] = None,
        period_flag: Union[str, None] = None
    ) -> Tuple[Literal['value', 'presence'], Union[Callable, None], Union[str, None]]:
        """
        Resolve a sample or group summary into a coloring kind and its reduction.

        A summary reduces a set of samples, or a set of sample groups, to one statement per map
        element. The names of 'SUMMARY_PRESENCE_SCHEMES' summarize presence. They say how many
        ('count'/'count_continuous') or exactly which ('membership') samples or groups contain the
        element. Any other name pools the element's values in the samples or groups as an
        aggregation ('_resolve_aggregation', '_summarize_entry_values'). This requires the layer's
        file to have a value column. The default of None summarizes presence with the colormap
        scheme left unresolved. Then '_membership_layer_colors' picks the scheme from the number of
        samples or groups. It colors by membership for ≤ 3 and by count for > 3. It uses a
        continuous count scale where a discrete one would run out of distinguishable colors or of
        room to label its bands. Room runs out past 'MAX_DISCRETE_COUNT_BANDS' bands.

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

        channel : dict
            The channel of the values ('_txt_channel'). It names the flag of the value type in an
            error message.

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
        voters = 'groups' if flag.endswith('-group-summary') else 'samples'
        # A border or a bar shows values only. Its fill already shows presence.
        if summary in SUMMARY_PRESENCE_SCHEMES and not channel['presence']:
            raise ConfigError(
                f"'{flag}' was given as '{summary}', which summarizes presence across the "
                f"{voters}. A {channel['role']} shows values only, and the fill already shows "
                f"presence. Give an aggregation, such as "
                f"'{'mean' if period is None else 'circular_mean'}', to pool the "
                f"{channel['values']}, or drop '{flag}'."
            )
        if summary in SUMMARY_PRESENCE_SCHEMES:
            return 'presence', None, SUMMARY_PRESENCE_SCHEMES[summary]
        # A layer of categorical values resolves its summary elsewhere
        # ('_build_categorical_txt_model'), so a categorical rule here is given to a layer without
        # categorical values. The message names what the user can do, which depends on the file.
        if summary in CATEGORICAL_RULES:
            if value_column is None:
                raise ConfigError(
                    f"'{flag}' was given as '{summary}', which chooses one categorical value "
                    f"across the {voters}. The text file at '{path}' has no value column, so there "
                    f"are no values to choose among. Summarize presence with "
                    f"{SUMMARY_PRESENCE_PHRASE} instead."
                )
            hint = (
                f"To read them as categorical values, give '--{channel['flag_stem']}-value-type "
                f"categorical'."
            ) if channel['categorical'] else "Pool them with an aggregation, such as 'mean'."
            raise ConfigError(
                f"'{flag}' was given as '{summary}', which chooses one categorical value across "
                f"the {voters}. The {channel['values']} of the text file at '{path}' are numbers. "
                f"{hint}"
            )
        if value_column is None:
            raise ConfigError(
                f"'{flag}' was given as '{summary}', which pools the values of samples or sample "
                f"groups, but the text file at '{path}' has no value column, so there are no "
                f"values to pool. Summarize presence with {SUMMARY_PRESENCE_PHRASE} instead, or "
                f"add a value column to the file."
            )
        aggregate = Mapper._resolve_aggregation(summary, flag, channel, period, period_flag)
        return 'value', aggregate, None

    @staticmethod
    def _sort_values(values: Iterable[str]) -> List[str]:
        """
        Sort categorical values, with the numbers in them compared as numbers.

        Plain sorting puts '10' before '2', and '-1' before '-2'. Where every value is a finite
        number, such as '-1', '0' and '1.5', the values are sorted as numbers. Otherwise each value
        is split into runs of digits and runs of other characters, and the runs of digits are
        compared as numbers, so 'cluster 2' comes before 'cluster 10'. A minus sign or a decimal
        point in such a value is compared as text. In both cases, the text itself breaks a tie, such
        as between '01' and '1'.

        Parameters
        ==========
        values : Iterable[str]
            The values to sort.

        Returns
        =======
        List[str]
            The values in sorted order.
        """
        values = list(values)
        numbers: Dict[str, float] = {}
        for value in values:
            try:
                number = float(value)
            except ValueError:
                break
            if not math.isfinite(number):
                break
            numbers[value] = number
        else:
            return sorted(values, key=lambda value: (numbers[value], value))

        def natural_key(text: str) -> Tuple[Tuple[Union[str, int], ...], str]:
            # Splitting on a captured group puts the runs of digits at the odd positions.
            parts = re.split(r'(\d+)', text)
            return tuple(
                int(part) if position % 2 else part for position, part in enumerate(parts)
            ), text
        return sorted(values, key=natural_key)

    @staticmethod
    def _unquote(text: str) -> str:
        """
        Drop a pair of double quotes that encloses a whole field, as R writes it.

        A doubled quote inside such a field stands for one quote, as R's 'write.csv' writes it. A
        quote anywhere else is part of the text, so a stray quote stays in the value it was typed
        in.

        Parameters
        ==========
        text : str
            The text of one field.

        Returns
        =======
        str
            The text without its enclosing quotes.
        """
        if len(text) >= 2 and text[0] == text[-1] == '"':
            return text[1:-1].replace('""', '"')
        return text

    @staticmethod
    def _read_txt_header(path: str) -> List[str]:
        """
        Read the column names of a tab-delimited text file as they are written.

        A normal read cannot keep a repeated name, since pandas renames the second copy, so
        'coverage' becomes 'coverage.1'. The header line is therefore read as a row of text. Quotes
        are read as text, as the file itself is read ('_read_element_txt'), and a pair of quotes
        enclosing a whole name is dropped ('_unquote').

        Parameters
        ==========
        path : str
            The path to the file.

        Returns
        =======
        List[str]
            The names in the order of the file, or an empty list for an empty file.
        """
        try:
            header = pd.read_csv(
                path, sep='\t', header=None, nrows=1, dtype=str, keep_default_na=False,
                quoting=csv.QUOTE_NONE
            )
        except pd.errors.EmptyDataError:
            return []
        return [Mapper._unquote(name) for name in header.iloc[0]]

    def _read_value_colors_txt(self, path: str, flag: str) -> Dict[str, str]:
        """
        Load a file that gives a color to each categorical value (artifact 'kegg-value-colors-txt').

        The file has the format of a kegg-membership-colors-txt file
        ('_read_membership_colors_txt'). Its first column holds the values and can be headed
        anything. Its 'COLORS_COLUMN' column holds six-digit color hex codes. A value is a whole
        name, so a comma in it does not name a combination of values. The order of the rows is the
        order of the values ('_resolve_value_colors').

        Parameters
        ==========
        path : str
            Path to the tab-delimited colors file.

        flag : str
            The command-line flag the file came from, used in error messages.

        Returns
        =======
        Dict[str, str]
            Maps each value to its color hex code, in the order of the rows.
        """
        filesnpaths.is_file_tab_delimited(path)
        # The header is read as the rows are, so a header written with quotes is accepted.
        raw_columns = utils.get_columns_of_TAB_delim_file(path, include_first_column=True)
        columns = [self._unquote(column) for column in raw_columns]
        if COLORS_COLUMN not in columns or len(columns) < 2 or columns[0] == COLORS_COLUMN:
            self.progress.end()
            raise ConfigError(
                f"The colors file given to '{flag}', at '{path}', should have two columns. The "
                f"first holds the categorical values and can be headed anything. The second is "
                f"headed '{COLORS_COLUMN}' and holds a color hex code for each value. Its columns "
                f"are: {', '.join(repr(column) for column in columns)}."
            )
        # The rows are keyed by line, not by value. The utility then does not stop at a value named
        # twice, and the check below reports it. A blank cell comes back as None.
        table = utils.get_TAB_delimited_file_as_dictionary(path, indexing_field=-1)
        value_colors: Dict[str, str] = {}
        for row in table.values():
            value = self._unquote(row[raw_columns[0]].strip())
            color = self._unquote((row[raw_columns[columns.index(COLORS_COLUMN)]] or '').strip())
            # Only hex codes are accepted, as in a kegg-membership-colors-txt file. Reserved colors
            # are compared as hex codes ('_check_reserved_colors').
            if not re.fullmatch(r'#[0-9A-Fa-f]{6}', color):
                shown = f"'{color}'" if color else 'blank'
                self.progress.end()
                raise ConfigError(
                    f"Each color in the file given to '{flag}', at '{path}', must be a six-digit "
                    f"hex code such as '#FFA500'. The color of '{value}' is {shown}."
                )
            # Two rows can name one value, exactly or with different spaces around it.
            if value in value_colors:
                self.progress.end()
                raise ConfigError(
                    f"The colors file given to '{flag}', at '{path}', names the value '{value}' in "
                    f"more than one row. Please give each value one row."
                )
            value_colors[value] = kgml.canonical_color(color)
        if not value_colors:
            self.progress.end()
            raise ConfigError(
                f"The colors file given to '{flag}', at '{path}', has its header but no rows. Each "
                f"row should name a value in its first column and give its color in the "
                f"'{COLORS_COLUMN}' column."
            )
        return value_colors

    def _resolve_value_colors(
        self,
        values: Iterable[str],
        value_colors_txt: Union[str, None],
        channel: dict,
        path: str
    ) -> Dict[str, str]:
        """
        Put the categorical values of a layer in order, and give each value a color.

        The order sets where each value sits in the legend. It breaks ties under 'plurality_order',
        and it decides every choice under 'order' ('CATEGORICAL_RULES'). For the fill, it also sets
        which value is drawn on top where map elements overlap. A border does not change which
        element is drawn on top.

        With a colors file, the order is the order of its rows. Every value of the layer needs a
        row. Rows that name values the layer does not have are reported and ignored, so one file can
        serve several runs. Two values cannot share a color, since the maps and the legend could not
        tell them apart.

        Without a colors file, the values are sorted, with the numbers in them compared as numbers
        ('_sort_values'). They then take the default colors of the channel in order. These are the
        18 colors of 'DEFAULT_VALUE_COLORS' for the fill, and the 7 colors of
        'DEFAULT_BORDER_VALUE_COLORS' for the border. More values than that need a colors file.

        Parameters
        ==========
        values : Iterable[str]
            The categorical values of the layer.

        value_colors_txt : Union[str, None]
            Path to a colors file ('_read_value_colors_txt'), or None for the default colors.

        channel : dict
            The channel of the values ('_txt_channel'). It gives the default colors and names the
            flag of the colors file.

        path : str
            Path to the layer's text file, for messages.

        Returns
        =======
        Dict[str, str]
            Maps each value to its color hex code, in the order of the values.
        """
        flag = f"--{channel['flag_stem']}-value-colors"
        element_type = channel['element_type']
        # A message about border values names the border column.
        file_values = f"the {element_type} file at '{path}'"
        if channel['role'] == 'border':
            file_values = f"the border column of {file_values}"
        default_colors = channel['default_value_colors']
        values = set(values)
        if value_colors_txt is None:
            ordered = self._sort_values(values)
            if len(ordered) > len(default_colors):
                raise ConfigError(
                    f"The {element_type} file at '{path}' has {len(ordered)} different categorical "
                    f"{channel['values']}. Anvi'o has {len(default_colors)} default colors for "
                    f"them. Please give each value a color of your choosing with '{flag}'."
                )
            return dict(zip(ordered, default_colors))

        def _listed(names: List[str]) -> str:
            # A long list is cut at 10 values.
            shown = ', '.join(repr(name) for name in names[:10])
            if len(names) > 10:
                return f"The first 10 are: {shown}."
            return f"The value is {shown}." if len(names) == 1 else f"The values are {shown}."

        file_colors = self._read_value_colors_txt(value_colors_txt, flag)
        missing = self._sort_values(values - set(file_colors))
        if missing:
            raise ConfigError(
                f"The colors file given to '{flag}', at '{value_colors_txt}', has no color for "
                f"{len(missing)} of the categorical {channel['values']} of the {element_type} file "
                f"at '{path}'. {_listed(missing)} Every value needs a color. Anvi'o does not fill "
                f"in missing colors, since chosen colors and default colors would then share one "
                f"legend."
            )
        extra = [value for value in file_colors if value not in values]
        if extra:
            self.run.warning(
                f"The colors file given to '{flag}', at '{value_colors_txt}', names {len(extra)} "
                f"{'value' if len(extra) == 1 else 'values'} that {file_values} does not have. "
                f"{_listed(extra)} {'It is' if len(extra) == 1 else 'They are'} ignored."
            )
        value_colors = {value: color for value, color in file_colors.items() if value in values}
        values_by_color: Dict[str, List[str]] = {}
        for value, color in value_colors.items():
            values_by_color.setdefault(color, []).append(value)
        shared = [names for names in values_by_color.values() if len(names) > 1]
        if shared:
            raise ConfigError(
                f"The colors file given to '{flag}', at '{value_colors_txt}', gives one color to "
                f"more than one value: "
                f"{'; '.join(', '.join(repr(name) for name in names) for names in shared)}. The "
                f"maps and the legend could not tell these values apart. Please give each value a "
                f"color of its own."
            )
        return value_colors

    @staticmethod
    def _resolve_categorical_rule(rule: Union[str, None], flag: str, channel: dict) -> str:
        """
        Check the rule that chooses one categorical value where several meet.

        Parameters
        ==========
        rule : Union[str, None]
            The rule given, or None for 'DEFAULT_CATEGORICAL_RULE'.

        flag : str
            The command-line flag the rule comes from, used in an error message.

        channel : dict
            The channel of the values ('_txt_channel'), for the error message.

        Returns
        =======
        str
            One of 'CATEGORICAL_RULES'.
        """
        if rule is None:
            return DEFAULT_CATEGORICAL_RULE
        if rule in CATEGORICAL_RULES:
            return rule
        raise ConfigError(
            f"'{flag}' was given as '{rule}'. The values of the {channel['noun']} layer are "
            f"categorical ('--{channel['flag_stem']}-value-type categorical'). Categorical values "
            f"cannot be added up or averaged. One of them is chosen instead. These are the rules "
            f"that choose: {CATEGORICAL_RULES_PHRASE}."
        )

    @staticmethod
    def _choose_categorical_value(
        votes: List[str],
        rule: str,
        order: Dict[str, int],
        seed: int,
        key: Tuple[str, ...]
    ) -> Tuple[Union[str, None], Union[str, None]]:
        """
        Choose one categorical value from the votes of genes, accessions, samples, or groups.

        The rules are those of 'CATEGORICAL_RULES'. 'plurality' and 'majority' can give no value. A
        lone vote is always a majority, and so are votes that all agree. A random choice is not
        drawn from a random number generator. A generator would give the same element different
        values on two maps, or in two runs that draw the maps in a different order. Each random
        choice is taken instead from a hash of the seed and of 'key', which names what the value is
        chosen for. The same element then gets the same value on every map. The hash is SHA-256,
        since Python's own hash of a string changes from one run to the next.

        Parameters
        ==========
        votes : List[str]
            One value per voter. A value can occur more than once.

        rule : str
            One of 'CATEGORICAL_RULES'.

        order : Dict[str, int]
            The position of each value in the order of the values ('_resolve_value_colors').

        seed : int
            The seed of random choices.

        key : Tuple[str, ...]
            What the value is chosen for, such as the accessions of a map element in one sample.

        Returns
        =======
        Tuple[Union[str, None], Union[str, None]]
            The chosen value, or None where the rule gives no value. Then how it was chosen. How is
            None where all of the votes are for one value. It is otherwise 'plurality', 'majority',
            'order', 'random', or 'none' where the rule gave no value.
        """
        candidates = sorted(set(votes), key=order.__getitem__)
        if len(candidates) == 1:
            return candidates[0], None
        counts = Counter(votes)
        if rule == 'majority':
            # More than half of the votes. At most one value can have that many.
            winners = [value for value in candidates if 2 * counts[value] > len(votes)]
            return (winners[0], 'majority') if winners else (None, 'none')
        if rule in ('plurality_order', 'plurality_random', 'plurality'):
            top = max(counts.values())
            candidates = [value for value in candidates if counts[value] == top]
            if len(candidates) == 1:
                return candidates[0], 'plurality'
            if rule == 'plurality':
                return None, 'none'
        if rule in ('plurality_order', 'order'):
            return candidates[0], 'order'
        digest = hashlib.sha256('\t'.join((str(seed),) + key).encode()).digest()
        return candidates[int.from_bytes(digest[:8], 'big') % len(candidates)], 'random'

    def _select_accession_values(
        self,
        rows_df: pd.DataFrame,
        rule: str,
        order: Dict[str, int],
        seed: int,
        sample: str,
        choices: Counter,
        salt: Tuple[str, ...] = ()
    ) -> Dict[str, str]:
        """
        Give each accession one categorical value, chosen among the values of its genes.

        This is the categorical counterpart of '_aggregate_accession_quantities'. Each row of a
        reaction file with a 'gene_id' column is one gene. Each gene gives one vote, and the rule
        chooses among the votes of an accession ('_choose_categorical_value'). Without genes,
        repeated rows are refused ('_read_element_txt'), so each accession has one row and there is
        nothing to choose.

        Parameters
        ==========
        rows_df : pandas.DataFrame
            Rows of one layer in one sample, carrying the '__accession' and '__value' columns.

        rule : str
            One of 'CATEGORICAL_RULES'.

        order : Dict[str, int]
            The position of each value in the order of the values.

        seed : int
            The seed of random choices.

        sample : str
            The sample of the rows, or '' for a file without samples. It is part of the key of a
            random choice.

        choices : collections.Counter
            Counts how each choice among different values was made. It is updated here.

        salt : Tuple[str, ...], ()
            What starts the key of a random choice ('_txt_channel'). The fill has none.

        Returns
        =======
        Dict[str, str]
            Keys are accessions, values are their categorical values. An accession for which the
            rule gave no value is left out.
        """
        values: Dict[str, str] = {}
        for accession, votes in rows_df.groupby('__accession')['__value'].agg(list).items():
            value, how = self._choose_categorical_value(
                votes, rule, order, seed, salt + ('gene', sample, accession)
            )
            if how is not None:
                choices[how] += 1
            # A rule can give no value. The accession then has no value in this sample.
            if value is not None:
                values[accession] = value
        return values

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

    @staticmethod
    def _txt_channel(
        element_type: Literal['reaction', 'compound'],
        role: Literal['fill', 'border', 'bar'] = 'fill'
    ) -> dict:
        """
        Describe one channel of a layer from a text file: its fill, its border, or its bars.

        A channel is the row of its role in 'TXT_CHANNEL_ROLES', with the layer named in its text.

        Parameters
        ==========
        element_type : Literal['reaction', 'compound']
            The layer of the channel.

        role : Literal['fill', 'border', 'bar'], 'fill'
            The channel.

        Returns
        =======
        dict
            The channel, with these keys:
            - 'role': the 'role' argument.
            - 'element_type': the 'element_type' argument.
            - 'flag_stem': what the options of the channel's values start with after '--', such as
              'reaction' in '--reaction-colormap' or 'reaction-border' in
              '--reaction-border-colormap'.
            - 'noun': what messages call the values, as in "a reaction value" or "the reaction
              border layer".
            - 'plural': what the channel colors, for messages about its colormaps.
            - 'name': the stem of the channel's key files, such as 'reactions' in
              'colorbar_reactions.pdf' or 'reactions_border' in 'colorbar_reactions_border.pdf'.
            - 'column': what messages call the column of the values.
            - 'values': what messages call the values of that column.
            - 'scale': what messages call the scale of the values, a color scale or a bar scale.
            - 'presence': whether the channel can show presence. A border and a bar show only
              values.
            - 'categorical' and 'period': whether the values can be categorical
              ('--*-value-type'), and whether they can repeat after a period ('--*-value-period').
              Bars take neither.
            - 'baseline': the value that every scale of the channel includes, or None. A bar grows
              from 0, so its scale always includes 0.
            - 'default_aggregation': how numbers are reduced where no aggregation is given and the
              values have no period. A border value given once per accession and repeated in the
              rows of its genes keeps its size under 'mean', but 'sum' would multiply it.
            - 'default_colormap': the colormap of numbers where none is given and the values have no
              period. It is None for a channel without colormaps.
            - 'default_colormap_limits': the part of that colormap that is used, or None for all of
              it.
            - 'default_centered_colormap': the colormap of a scale that a normalization centers.
            - 'default_value_colors': the colors of categorical values where no file gives them.
            - 'choice_salt': what starts the key of each random choice among categorical values
              ('_choose_categorical_value'). The fill and the border then choose independently, and
              the keys of the fill stay as they were before borders existed.
            - 'consequence' and 'consequence_plural': what happens to an element without a value, as
              in "an element is left uncolored" or "elements show no border".
        """
        return {
            'role': role,
            'element_type': element_type,
            **{
                key: value.format(element_type=element_type) if isinstance(value, str) else value
                for key, value in TXT_CHANNEL_ROLES[role].items()
            }
        }

    def _txt_model_common(
        self,
        name: str,
        data: dict,
        path: str,
        membership_colors_txt: Union[str, None],
        channel: dict
    ) -> dict:
        """
        Build the keys that every layer model from a text file has ('_map_elements').

        Parameters
        ==========
        name : str
            Colorbar filename stem ('reactions'/'compounds').

        data : dict
            The '_read_element_txt' result for this layer.

        path : str
            Path to the layer's text file. Messages name it.

        membership_colors_txt : Union[str, None]
            Path to a kegg-membership-colors-txt file for the layer, or None.

        channel : dict
            The channel of the values ('_txt_channel').

        Returns
        =======
        dict
            The keys shared by every layer model. These include 'channel' and 'value_column', the
            header of the column of the values or None. Messages name the column by its header.
        """
        element_type = data['element_type']
        membership_colors_flag = f'--{element_type}-membership-colors'
        common = {
            'name': name,
            'path': path,
            'element_type': element_type,
            'channel': channel,
            'value_column': data['value_column'],
            'use_reaction_attribute': data['reaction_source'] == 'Reaction',
            'accessions': set(data['df']['__accession']),
            'membership_colors_flag': membership_colors_flag
        }
        # The membership colors file is read here. Its names are checked against the run's samples
        # or groups later, in '_map_elements', where those are known. Reading the file first means a
        # malformed file is reported as malformed, not as a file full of unrecognized names.
        if membership_colors_txt is None:
            common['membership_colors'] = None
            common['membership_combo_colors'] = {}
        else:
            membership_colors, combo_colors = self._read_membership_colors_txt(
                membership_colors_txt, membership_colors_flag
            )
            common['membership_colors'] = membership_colors
            common['membership_combo_colors'] = combo_colors
        return common

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
        individual_colormap: Union[str, mcolors.Colormap, None],
        individual_colormap_limits: Union[Tuple[float, float], None],
        membership_colors_txt: Union[str, None],
        reverse_overlay: bool,
        all_sample_names: Union[List[str], None],
        group_samples: Union[Dict[str, List[str]], None],
        sample_summary: Union[str, None],
        group_summary: Union[str, None],
        value_limits: Union[Tuple[Union[float, None], Union[float, None]], None] = None,
        individual_value_limits: Union[Tuple[Union[float, None], Union[float, None]], None] = None,
        value_center: Union[float, None] = None,
        individual_value_center: Union[float, None] = None,
        element_normalization: Union[str, None] = None,
        element_normalization_label: Union[str, None] = None,
        value_period: Union[float, str, None] = None,
        value_colors_txt: Union[str, None] = None,
        seed: int = DEFAULT_CATEGORICAL_SEED,
        slice_by: Union[str, None] = None,
        channel: Union[dict, None] = None
    ) -> dict:
        """
        Build one layer's coloring model for '_map_elements' from a per-layer reader result.

        The layer colors elements by mode in each of two map contexts: 'unified_mode' for the
        'unified' map and 'individual_mode' for the per-sample or per-group maps. Without a 'sample'
        column there is nothing to summarize, so the two agree: 'quantitative' with a value column
        ('_read_element_txt'), or 'single' (one fixed presence color) without one. An
        '--original-color' layer is the exception. It always pairs a 'unified_mode' of 'original'
        with an 'individual_mode' of 'membership'. The reason is given where it is built. A 'color'
        of 'none' gives a layer without values the mode 'none' in both contexts, with or without
        samples. Its fills are left uncolored, and it draws only the rings of its border column.

        With a 'sample' column the two contexts can differ, since the value column only colors a
        single sample's magnitude while the summaries independently choose what the views across
        samples and across groups show. Ungrouped, 'sample_summary' colors the 'unified' map and the
        per-sample maps show each sample on its own; grouped, 'group_summary' colors the 'unified'
        map and 'sample_summary' colors each group's map from its samples. A summary naming an
        aggregation pools values ('quantitative'), while a presence summary or no summary colors
        presence ('membership'), as 'SUMMARY_PRESENCE_SCHEMES' and '_resolve_summary' decide. This
        lets a value layer show presence in the 'unified' map and values per sample. A layer may
        lack a 'sample' column in a run where the other layer has one. Such a layer has no values
        for each sample or group, so '_map_elements' holds it constant across the per-sample/group
        maps.

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
        'aggregation'/'color'/'colormap'/'colormap_limits'/'individual_colormap'/
        'individual_colormap_limits'/'membership_colors_txt'/'reverse_overlay'/'sample_summary'/
        'group_summary' are this layer's per-layer parameters; 'all_sample_names'/'group_samples'
        carry the shared sample space and grouping. Returns the layer model dict '_map_elements'
        consumes. 'value_limits' and 'individual_value_limits' bound the two quantitative scales
        independently — the 'unified' map's and the one its per-sample or per-group maps share —
        since a summary can put the two on quite different scales. Each is refused where its own
        context is not colored by value, and setting one of the two while both contexts are colored
        by value earns a warning, the other scale being an easy one to forget. 'individual_colormap'
        colors that same per-sample/per-group scale from a colormap of its own: where the 'unified'
        map shows a summary of a different kind from the values behind it — their spread rather than
        their magnitude, for example — one colormap across both invites reading the two as the same
        quantity. Left unset, both contexts share the layer's 'colormap', and it is refused wherever
        the per-sample/per-group context is not colored by value. 'value_center' and
        'individual_value_center' put a value at the middle of those same two scales, each of which
        then runs the same distance either side of it however lopsided its own values are, and are
        refused and warned about on exactly the same footing as the limits.

        'element_normalization' rescales each sample's or group's value for a map element against
        the same element's values across all samples or groups ('_resolve_element_normalization').
        This allows the enrichment or depletion of an element in a sample to be displayed rather
        than how much of it there is. Therefore, normalization requires a value column and samples.
        Normalization applies to the per-sample/per-group context alone. The 'unified' map
        summarizes the unnormalized values of the samples or groups. A normalization whose neutral
        value is zero uses a centered scale on the diverging colormap, 'DEFAULT_CENTERED_COLORMAP',
        unless overridden by 'individual_value_center' and 'individual_colormap'.
        'element_normalization_label' names the rescaled quantity on that scale's colorbar in place
        of the label the normalization derives for itself.

        'value_period' says that the values repeat after that period, as clock times in hours repeat
        every 24. Every reduction of the values is then circular ('CIRCULAR_AGGREGATIONS'). An
        aggregation left as None is 'circular_mean' with a period and 'sum' without one. A period
        refuses every other aggregation and summary of values. Its only normalization is
        'difference_from_circular_mean' ('CIRCULAR_ELEMENT_NORMALIZATIONS'). Presence summaries are
        unaffected. Each scale of values runs from 0 to the period ('_map_elements'). The scale of
        offsets from a normalization runs from -period / 2 to period / 2 instead. Other limits are
        refused, and so are centers. A colormap left as None is 'DEFAULT_PERIOD_COLORMAP'.

        A layer of categorical values ('VALUE_TYPES') is built by '_build_categorical_txt_model'
        instead. 'value_colors_txt' and 'seed' apply only to such a layer.

        'slice_by' says what the elements of the layer are sliced by ('SLICE_UNITS'), or None.
        Slices of genes take their values from the rows of each gene ('_build_slice_model'). The
        values of the accessions are then not drawn, so an accession whose genes reduce to an
        undefined value is not reported.

        'channel' is the channel that the values color ('_txt_channel'), or None for the fill.
        Messages name the options of the channel, which start with its 'flag_stem'. An option left
        as None takes the channel's default. An aggregation left as None is then 'mean' for a border
        without a period.
        """
        if channel is None:
            channel = self._txt_channel(data['element_type'])
        flag_stem = channel['flag_stem']
        # A color of 'none' leaves the fills of a layer without values uncolored, so that the layer
        # draws only the rings of its border column. The layer takes the fill mode 'none' on every
        # map. Every option of the fill is refused, since it would have nothing to act on. The CLI
        # makes these checks too.
        if color == 'none':
            element_type = data['element_type']
            none_flag = f"--{element_type}-color none"
            unfilled = (
                f"'{none_flag}' leaves the fills of the {element_type} layer uncolored, so that "
                f"only the rings of its border column are drawn."
            )
            if data['border_column'] is None:
                raise ConfigError(
                    f"{unfilled} No border column was named with '--{element_type}-border-column', "
                    f"so the layer would draw nothing. Name the column whose values color the "
                    f"rings, or give a color hex code."
                )
            if data['value_column'] is not None:
                raise ConfigError(
                    f"{unfilled} The file at '{path}' has a fill value column, "
                    f"'{data['value_column']}'. Please drop '{none_flag}' to color the fills by it."
                )
            fill_flags = list(dict.fromkeys(f"'--{element_type}-{option}'" for option, given in (
                ('gene-aggregation', gene_aggregation is not None),
                ('accession-aggregation', accession_aggregation is not None),
                ('colormap', colormap is not None and colormap is not False),
                ('colormap', colormap_limits is not None),
                ('individual-colormap', individual_colormap is not None),
                ('individual-colormap', individual_colormap_limits is not None),
                ('membership-colors', membership_colors_txt is not None),
                ('reverse-overlay', reverse_overlay),
                ('sample-summary', sample_summary is not None),
                ('group-summary', group_summary is not None),
                ('value-limits', value_limits is not None),
                ('individual-value-limits', individual_value_limits is not None),
                ('value-center', value_center is not None),
                ('individual-value-center', individual_value_center is not None),
                ('element-normalization', element_normalization is not None),
                ('element-normalization', element_normalization_label is not None),
                ('value-period', value_period is not None),
                ('slice-by', slice_by is not None)
            ) if given))
            if fill_flags:
                raise ConfigError(
                    f"These options color the fills of the {element_type} layer: "
                    f"{', '.join(fill_flags)}. {unfilled} Please drop these options. The options "
                    f"of the rings start with '--{element_type}-border-'."
                )
            return {
                **self._txt_model_common(name, data, path, None, channel),
                'unified_mode': 'none',
                'individual_mode': 'none'
            }
        if data['value_type'] == 'categorical':
            return self._build_categorical_txt_model(
                name, data, path, gene_aggregation, accession_aggregation, color, colormap,
                colormap_limits, membership_colors_txt, reverse_overlay, all_sample_names,
                group_samples, sample_summary, group_summary, value_colors_txt, seed, {
                    f"--{flag_stem}-value-limits": value_limits,
                    f"--{flag_stem}-individual-value-limits": individual_value_limits,
                    f"--{flag_stem}-value-center": value_center,
                    f"--{flag_stem}-individual-value-center": individual_value_center,
                    f"--{flag_stem}-element-normalization": (
                        element_normalization if element_normalization is not None
                        else element_normalization_label
                    ),
                    f"--{flag_stem}-value-period": value_period,
                    f"--{flag_stem}-individual-colormap": individual_colormap
                }, channel
            )
        element_type = data['element_type']
        noun = channel['noun']
        use_reaction_attribute = data['reaction_source'] == 'Reaction'
        df = data['df']
        value_column = data['value_column']
        has_sample = data['sample_names'] is not None
        accessions = set(df['__accession'])
        membership_colors_flag = f'--{element_type}-membership-colors'
        common = self._txt_model_common(name, data, path, membership_colors_txt, channel)

        # Limits and centers are checked before any value is read, so that a pair or a center anvi'o
        # cannot make sense of is reported as the configuration error it is rather than surfacing
        # later as a strange scale. What each can shape at all is settled here too; which context is
        # actually colored by value is not known until the summaries below have chosen the modes.
        # The centers are resolved against the limits, which is where a center placed out of a
        # limited scale's reach is caught, so the limits are resolved first.
        value_limits_flag = f'--{flag_stem}-value-limits'
        individual_value_limits_flag = f'--{flag_stem}-individual-value-limits'
        value_center_flag = f'--{flag_stem}-value-center'
        individual_value_center_flag = f'--{flag_stem}-individual-value-center'
        value_limits = self._resolve_value_limits(value_limits, value_limits_flag)
        individual_value_limits = self._resolve_value_limits(
            individual_value_limits, individual_value_limits_flag
        )
        # A bar grows from 0, so each of its scales always includes 0 ('TXT_CHANNEL_ROLES'). Limits
        # that leave 0 out would draw bars that start elsewhere, and their lengths would misstate
        # the values. A scale can reach below 0. Where it reaches above 0 too, its bars grow both
        # ways from 0.
        baseline = channel['baseline']
        for limits, flag in (
            (value_limits, value_limits_flag),
            (individual_value_limits, individual_value_limits_flag)
        ):
            if baseline is None or limits is None:
                continue
            limit_min, limit_max = limits
            shown = ' to '.join('open' if limit is None else f'{limit:g}' for limit in limits)
            if (limit_min is not None and limit_min > baseline) or (
                limit_max is not None and limit_max < baseline
            ):
                raise ConfigError(
                    f"'{flag}' was given limits that leave out {baseline:g}: {shown}. A bar grows "
                    f"from {baseline:g}, so its scale always includes {baseline:g}, and a bar that "
                    f"started elsewhere would misstate its value. Please give limits that include "
                    f"{baseline:g}, or leave an end open."
                )
        # A period says that the values repeat, as clock times do every 24 hours. It is resolved
        # before the normalization, which it rules out.
        value_period_flag = f'--{flag_stem}-value-period'
        value_period = self._resolve_value_period(value_period, value_period_flag)

        # The normalization is settled before the centers are, since one whose neutral value is zero
        # centers the per-sample/per-group scale on zero unless asked for another center, and that
        # default has to be validated against the limits on that scale like any other center. It
        # rescales the value of each sample or group against the element's values across all samples
        # or groups, so it needs both a value column to rescale and samples.
        element_normalization_flag = f'--{flag_stem}-element-normalization'
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
                f"colorbar of the {noun} layer, but no normalization was given for it to label. "
                f"That colorbar is labeled by the {channel['column']} of the file at '{path}' when "
                f"nothing rescales it."
            )
        if element_normalization_label is not None and not element_normalization_label.strip():
            raise ConfigError(
                f"The label given to '{element_normalization_flag}' for the colorbar of the {noun} "
                f"layer's rescaled scale is blank, so nothing would say what that scale shows. "
                f"Give the option a label with something in it, or give it the normalization alone "
                f"and let anvi'o compose one."
            )
        if element_normalization is not None:
            if value_column is None:
                raise ConfigError(
                    f"'{element_normalization_flag}' rescales the values of the {noun} layer, but "
                    f"the file at '{path}' has no value column, so that layer is colored by "
                    f"presence and has no values to rescale."
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
                    element_normalization, element_normalization_flag, channel, value_period,
                    value_period_flag
                )
            )
            # With a period, the offsets run from -period / 2 to period / 2. That scale is fixed
            # in '_map_elements', with 0 at its middle, so no center is set here. Its two ends are
            # the same point, so its default colormap is the cyclic one. A diverging colormap
            # would give the two ends opposite colors.
            if value_period is not None:
                if individual_colormap is None:
                    individual_colormap = DEFAULT_PERIOD_COLORMAP
            else:
                if centered and individual_value_center is None:
                    individual_value_center = 0.0
                    center_from_normalization = True
                if individual_colormap is None and centered:
                    individual_colormap = channel['default_centered_colormap']

        value_center = self._resolve_value_center(
            value_center, value_center_flag, value_limits, value_limits_flag
        )
        # A center that came from a centered normalization rather than from the user is reported
        # under the normalization's own flag: nobody gave a center, so naming the center's flag
        # would name a flag that was never used.
        individual_center_source_flag = (
            element_normalization_flag if center_from_normalization
            else individual_value_center_flag
        )
        individual_value_center = self._resolve_value_center(
            individual_value_center, individual_center_source_flag, individual_value_limits,
            individual_value_limits_flag, from_normalization=center_from_normalization
        )
        if value_column is None:
            for setting, flag, verb in (
                (value_limits, value_limits_flag, 'bound'),
                (individual_value_limits, individual_value_limits_flag, 'bound'),
                (value_center, value_center_flag, 'center'),
                (individual_value_center, individual_value_center_flag, 'center')
            ):
                if setting is not None:
                    raise ConfigError(
                        f"'{flag}' {verb}s the color scale spanning the values of the {noun} "
                        f"layer, but the file at '{path}' has no value column, so that layer is "
                        f"colored by presence and has no scale of values to {verb}."
                    )
            if value_period is not None:
                raise ConfigError(
                    f"'{value_period_flag}' says how far the values of the {noun} layer go before "
                    f"they repeat, but the file at '{path}' has no value column, so that layer is "
                    f"colored by presence and has no values."
                )
        if not has_sample:
            for setting, flag, single_scale_flag, verb in (
                (individual_value_limits, individual_value_limits_flag, value_limits_flag, 'bound'),
                (individual_value_center, individual_value_center_flag, value_center_flag, 'center')
            ):
                if setting is not None:
                    raise ConfigError(
                        f"'{flag}' {verb}s the color scale shared by the maps of the individual "
                        f"samples or groups of the {noun} layer, but the file at '{path}' has no "
                        f"'sample' column, so it draws no such maps and its values take a single "
                        f"scale. {verb.capitalize()} that one with '{single_scale_flag}'."
                    )
        # The summary that colors the 'unified' map is the group summary with groups, and the sample
        # summary without them. A spread there ('CIRCULAR_SPREADS') gives that map a scale, a
        # colormap, and a label of its own.
        unified_summary_flag = (
            f"--{flag_stem}-{'group' if group_samples is not None else 'sample'}-summary"
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
            individual_limits = (
                (-value_period / 2, value_period / 2) if offsets else (0.0, value_period)
            )
            for setting, flag, limits in (
                (value_limits, value_limits_flag, (0.0, value_period)),
                (individual_value_limits, individual_value_limits_flag, individual_limits)
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
                (individual_value_center, individual_value_center_flag)
            ):
                if setting is not None and offsets and flag == individual_value_center_flag:
                    raise ConfigError(
                        f"'{flag}' centers the color scale of the maps of the individual samples "
                        f"or groups of the {noun} layer, but '{element_normalization_flag} "
                        f"{element_normalization}' already puts 0 at the middle of that scale. It "
                        f"runs from {individual_limits[0]:g} to {individual_limits[1]:g}. Please "
                        f"drop '{flag}'."
                    )
                if setting is not None:
                    raise ConfigError(
                        f"'{flag}' centers a color scale of the {noun} layer, but "
                        f"'{value_period_flag}' fixes each such scale to run from 0 to "
                        f"{value_period:g}. Values that repeat have no middle to read the others "
                        f"against. Please use only one of the two options."
                    )

        # The colormap of the per-sample/per-group scale is settled on the same footing as the limits
        # on it: what it can color at all is checked here, before any value is read, while whether
        # that context is colored by value at all waits on the summaries below.
        individual_colormap_flag = f'--{flag_stem}-individual-colormap'
        colormap_flag = f'--{flag_stem}-colormap'
        if individual_colormap is not None:
            if value_column is None:
                raise ConfigError(
                    f"'{individual_colormap_flag}' colors the maps of the individual samples or "
                    f"groups of the {noun} layer by value, but the file at '{path}' has no value "
                    f"column, so that layer is colored by presence rather than along a scale of "
                    f"values. '{colormap_flag}' is what colors a presence scale."
                )
            if not has_sample:
                raise ConfigError(
                    f"'{individual_colormap_flag}' colors the maps of the individual samples or "
                    f"groups of the {noun} layer, but the file at '{path}' has no 'sample' column, "
                    f"so it draws no such maps and its values take a single scale. Color that one "
                    f"with '{colormap_flag}'."
                )

        if color == 'original':
            # A reaction presence layer highlighted in the reference map's colors
            # ('--original-color'), drawn by the reference-color drawer rather than the element
            # engine, since both the colors and their render order come from the reference map, so
            # there is no channel for data. The CLI restricts this to a reaction file with no value
            # column, magnitude having nowhere to go. With a 'sample' column, the 'unified' map is
            # the union of the samples and each sample also gets its own map; without one,
            # 'membership' just carries the accession set the single map colors. 'individual_mode'
            # is 'membership' rather than 'original', so that a grouped run takes the same path as
            # the database and pangenome inputs. Their per-group maps show within-group source
            # counts. 'source_accessions' is keyed by sample, so the reference-color branch for the
            # individual maps can serve only samples, not groups. The reference-color drawer finds
            # elements by the KO IDs of a map's ortholog entries, so a file of KEGG reaction IDs
            # would match nothing and draw a blank map, and there is no compound layer to draw at
            # all.
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
            if common['membership_colors'] is not None:
                raise ConfigError(
                    f"'{membership_colors_flag}' gives a color to each sample or group, but "
                    f"'--original-color' takes both its colors and its drawing order from the "
                    f"reference map, so there is no color left to choose. Please use only one of "
                    f"them."
                )
            for setting, flag, verb in (
                (value_limits, value_limits_flag, 'bound'),
                (individual_value_limits, individual_value_limits_flag, 'bound'),
                (value_center, value_center_flag, 'center'),
                (individual_value_center, individual_value_center_flag, 'center')
            ):
                if setting is not None:
                    raise ConfigError(
                        f"'{flag}' {verb}s a color scale spanning a layer's values, but "
                        f"'--original-color' takes both its colors and its drawing order from the "
                        f"reference map, so nothing is colored by value and there is no scale to "
                        f"{verb}. Please use only one of them."
                    )
            if individual_colormap is not None:
                raise ConfigError(
                    f"'{individual_colormap_flag}' colors the maps of the individual samples or "
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
                'individual_mode': 'membership',
                'membership': membership,
                'source_accessions': source_accessions,
                'color_hexcode': 'original'
            }

        undefined: Set[str] = set()
        # A spread is a plain number, not a point in the period, so only the summary that colors the
        # 'unified' map takes one. Within a sample, the KOs of an element must still give a point in
        # the period. With groups, the sample summary colors each group's map, and the group summary
        # then pools those values.
        for name, flag in (
            (gene_aggregation, f'--{flag_stem}-gene-aggregation'),
            (accession_aggregation, f'--{flag_stem}-accession-aggregation'),
            (sample_summary if group_samples is not None else None,
             f'--{flag_stem}-sample-summary')
        ):
            if name in CIRCULAR_SPREADS:
                raise ConfigError(
                    f"'{flag}' was given as '{name}', which measures how far values are spread "
                    f"apart. A spread is a plain number, not a point in the period, so it can only "
                    f"summarize the {'groups' if group_samples is not None else 'samples'} for the "
                    f"'unified' map, with '{unified_summary_flag}'."
                )
        # An aggregation left unset is 'circular_mean' for values with a period. Otherwise it is
        # 'sum' for the fill and 'mean' for the border ('_txt_channel'). A colormap of values left
        # unset is cyclic for values with a period. Otherwise it is the channel's default, of which
        # only a part may be used.
        default_aggregation = (
            'circular_mean' if value_period is not None else channel['default_aggregation']
        )
        # Messages that suggest a summary of values name one that the period accepts.
        suggested_summary = 'circular_mean' if value_period is not None else 'mean'
        default_colormap = (
            DEFAULT_PERIOD_COLORMAP if value_period is not None else channel['default_colormap']
        )
        default_colormap_limits = (
            None if value_period is not None else channel['default_colormap_limits']
        )

        def _limits(default_limits):
            # Limits given without a colormap apply to the default colormap. Where neither is given,
            # the default colormap takes its default limits.
            if colormap is None and colormap_limits is None:
                return default_limits
            return colormap_limits
        if gene_aggregation is None:
            gene_aggregation = default_aggregation
        if accession_aggregation is None:
            accession_aggregation = default_aggregation
        # The run reports the aggregations of a layer with a border ('map_kegg_pathways_txt').
        common['aggregations'] = {'gene': gene_aggregation, 'accession': accession_aggregation}
        # Resolved here, before any values are aggregated, so that an unusable aggregation name is
        # reported as a configuration error rather than reaching pandas: the per-accession values
        # below are computed by passing the name itself to a groupby. 'aggregate' reduces a map
        # element's several accessions, which is the level the drawing colorers work at; the gene
        # aggregation is applied only where the per-accession values are built.
        aggregate = self._resolve_aggregation(
            accession_aggregation, f'--{flag_stem}-accession-aggregation', channel, value_period,
            value_period_flag
        ) if value_column is not None else None
        # A circular aggregation is not a pandas name, so its function reduces a gene's rows. A
        # pandas name is passed as the name itself, which pandas applies faster.
        gene_reduction = gene_aggregation
        if value_column is not None:
            gene_aggregate = self._resolve_aggregation(
                gene_aggregation, f'--{flag_stem}-gene-aggregation', channel, value_period,
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
                    'individual_mode': 'single',
                    'color_hexcode': color
                }
            unified_values = self._aggregate_accession_quantities(
                df, value_column, gene_reduction, path, channel, undefined
            )
            if slice_by != 'gene':
                self._warn_undefined_values(undefined, gene_aggregation, path, channel)
            # A channel without colormaps, such as a bar, sets lengths rather than colors.
            cmap = None
            if channel['default_colormap'] is not None:
                limits = _limits(default_colormap_limits)
                cmap = self._resolve_sequential_colormap(
                    colormap if colormap is not None else default_colormap, limits,
                    subject=channel['plural']
                )
                self._check_centered_colormap(
                    cmap, limits, value_center, value_center_flag, colormap_flag,
                    subject=channel['plural']
                )
            return {
                **common,
                'unified_mode': 'quantitative',
                'individual_mode': 'quantitative',
                'cmap': cmap,
                # There are no per-sample maps to give a colormap of their own, so
                # 'individual_colormap' was refused above. The one map's scale serves both contexts.
                'individual_cmap': cmap,
                'reverse_overlay': reverse_overlay,
                'unified_values': unified_values,
                'sample_values': None,
                'aggregate': aggregate,
                'colorbar_label': value_column,
                # A normalization was refused above, there being no samples to rescale against, so
                # the one scale of this layer carries the value column's own name.
                'individual_colorbar_label': value_column,
                'element_normalize': None,
                'value_limits': value_limits,
                'individual_value_limits': None,
                'value_center': value_center,
                'individual_value_center': None,
                'value_period': value_period
            }

        # Resolve the two summaries into the mode of each map context. Ungrouped, the per-sample
        # maps show one sample each, so no summary applies to them: they are the sample's own
        # magnitude with a value column and its presence without one.
        grouped = group_samples is not None
        sample_kind, sample_aggregate, sample_scheme = self._resolve_summary(
            sample_summary, value_column, f'--{flag_stem}-sample-summary', path, channel,
            value_period, value_period_flag
        )
        if grouped:
            group_kind, group_aggregate, group_scheme = self._resolve_summary(
                group_summary, value_column, f'--{flag_stem}-group-summary', path, channel,
                value_period, value_period_flag
            )
            if group_kind == 'value' and sample_kind != 'value':
                raise ConfigError(
                    f"'--{flag_stem}-group-summary' was given as '{group_summary}', which pools "
                    f"the values of the sample groups, but '--{flag_stem}-sample-summary' "
                    f"summarizes each group's samples by presence rather than by value, so the "
                    f"groups have no values to pool. Please also set "
                    f"'--{flag_stem}-sample-summary' to an aggregation such as "
                    f"'{suggested_summary}'."
                )
            unified_kind, unified_scheme = group_kind, group_scheme
            individual_kind = sample_kind
        else:
            unified_kind, unified_scheme = sample_kind, sample_scheme
            individual_kind = 'presence' if value_column is None else 'value'

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
        # A border or a bar shows values only. A map that a summary would color by presence shows no
        # border or bar, which the mode None says ('_txt_channel'). Its summaries take no presence
        # names ('_resolve_summary'), so such a map is one without a summary.
        values_only = not channel['presence']
        individual_mode = 'quantitative' if individual_kind == 'value' else 'membership'
        if values_only:
            unified_mode = 'quantitative' if unified_kind == 'value' else None
            individual_mode = 'quantitative' if individual_kind == 'value' else None

        # The option that chose the 'unified' map's presence scheme, so that a message about that
        # scheme names the option that can change it rather than the one this input refuses: the
        # group summary colors the 'unified' map of a grouped run, and the sample summary colors it
        # otherwise, exactly as 'unified_scheme' was taken above.
        summary_flag = f"--{flag_stem}-{'group' if grouped else 'sample'}-summary"

        model = {
            **common,
            'unified_mode': unified_mode,
            'individual_mode': individual_mode,
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
            'individual_value_limits': individual_value_limits,
            'value_center': value_center,
            'individual_value_center': individual_value_center,
            # A center that a normalization supplied is reported under the normalization's flag.
            'individual_value_center_flag': individual_center_source_flag,
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
        # A border shows no ring, and a bar no bar, on a map without a summary, as their summaries
        # take no presence names. Their messages say so.
        absent = 'no border' if channel['role'] == 'border' else 'no bars'
        scale_of = (
            f"{channel['scale']} of {channel['values']}" if channel['role'] == 'border'
            else channel['scale']
        )
        no_values = (
            f"Those maps show {absent} here, since '--{flag_stem}-sample-summary' was not given. "
            f"Give it an aggregation, such as '{suggested_summary}', to show {channel['values']} "
            f"on the group maps"
        )
        if element_normalization is not None and model['individual_mode'] != 'quantitative':
            if values_only:
                raise ConfigError(
                    f"'{element_normalization_flag}' rescales the {channel['values']} of the maps "
                    f"of the individual groups. {no_values}."
                )
            raise ConfigError(
                f"'{element_normalization_flag}' rescales the values that color the maps of the "
                f"individual groups, but those maps are not colored by value here: "
                f"'--{flag_stem}-sample-summary' summarizes each group's samples by presence "
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
            (value_limits, value_limits_flag, individual_value_limits_flag, 'bound'),
            (value_center, value_center_flag, individual_value_center_flag, 'center')
        ):
            if setting is not None and model['unified_mode'] != 'quantitative' and values_only:
                raise ConfigError(
                    f"'{flag}' {verb}s the {scale_of} on the 'unified' map. That map shows "
                    f"{absent} here, since '{summary_flag}' was not given. Give it an aggregation, "
                    f"such as '{suggested_summary}', or {verb} the scale of the maps of the "
                    f"individual samples or groups with '{other_flag}' instead."
                )
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
            (individual_value_limits, individual_value_limits_flag, value_limits_flag, 'bound'),
            (individual_value_center, individual_value_center_flag, value_center_flag, 'center')
        ):
            if setting is not None and model['individual_mode'] != 'quantitative' and values_only:
                raise ConfigError(
                    f"'{flag}' {verb}s the {channel['scale']} shared by the maps of the individual "
                    f"groups. {no_values}, or {verb} the 'unified' map's own scale with "
                    f"'{other_flag}'."
                )
            if setting is not None and model['individual_mode'] != 'quantitative':
                raise ConfigError(
                    f"'{flag}' {verb}s the color scale shared by the maps of the individual "
                    f"groups, but those maps are not colored by value here: "
                    f"'--{flag_stem}-sample-summary' summarizes each group's samples by presence "
                    f"rather than by pooling their values. Set it to an aggregation such as "
                    f"'{suggested_summary}' to color the group maps by value, or {verb} the "
                    f"'unified' map's own scale with '{other_flag}'."
                )
        if (
            individual_colormap is not None and model['individual_mode'] != 'quantitative'
            and values_only
        ):
            raise ConfigError(
                f"'{individual_colormap_flag}' colors the scale shared by the maps of the "
                f"individual groups. {no_values}, or color the 'unified' map's own scale with "
                f"'{colormap_flag}'."
            )
        if individual_colormap is not None and model['individual_mode'] != 'quantitative':
            raise ConfigError(
                f"'{individual_colormap_flag}' colors the scale shared by the maps of the "
                f"individual groups, but those maps are not colored along a scale of values here: "
                f"'--{flag_stem}-sample-summary' summarizes each group's samples by presence "
                f"rather than by pooling their values, and presence there is colored by "
                f"'--group-colormap'. Set the sample summary to an aggregation such as "
                f"'{suggested_summary}' to color the group maps by value, or color the 'unified' "
                f"map's own scale with '{colormap_flag}'."
            )

        # The 'unified' scale and the scale of the individual maps take separate limits. A summary
        # can span much less than the values it summarizes. For example, a mean across samples spans
        # less than the samples do, so limits given for one scale do not apply to the other. Where
        # both scales color by value and only one has limits, the user may not have meant to leave
        # the other without limits. This warns about it. It is not refused, since limiting only one
        # scale can be intended. With a period, both scales have fixed limits, so there is nothing
        # to warn about.
        if (
            model['unified_mode'] == 'quantitative' and model['individual_mode'] == 'quantitative'
            and (value_limits is None) != (individual_value_limits is None)
            and value_period is None
        ):
            if value_limits is None:
                given_flag, other_flag = individual_value_limits_flag, value_limits_flag
                other_subject = "the 'unified' map"
            else:
                given_flag, other_flag = value_limits_flag, individual_value_limits_flag
                other_subject = f"the maps of the individual {'groups' if grouped else 'samples'}"
            self.run.warning(
                f"'{given_flag}' bounds the {channel['scale']} of the {noun} layer, but "
                f"'{other_flag}' was not given, so the scale of {other_subject} still spans "
                f"whatever its own values happen to reach. The two are bounded separately because "
                f"they can differ a great deal -- a summary such as a mean across samples spans "
                f"less than the samples themselves do -- so this may well be what you intend. If "
                f"it is not, give '{other_flag}' limits of its own."
            )

        # Centering only one of the two scales also warns, for the same reason as the limits above.
        # There is a second reason. Without '--*-individual-colormap', both scales use the same
        # colormap. The middle color then means the center on one map, and the midpoint of that
        # map's values on the other. A center that a normalization supplied does not warn. The user
        # did not ask for it, and a normalization colors its scale with a diverging colormap that
        # the summary does not share.
        if (
            model['unified_mode'] == 'quantitative' and model['individual_mode'] == 'quantitative'
            and (value_center is None) != (individual_value_center is None)
            and not center_from_normalization
        ):
            if value_center is None:
                given_flag, other_flag = individual_value_center_flag, value_center_flag
                given_center = individual_value_center
                other_subject = "the 'unified' map"
            else:
                given_flag, other_flag = value_center_flag, individual_value_center_flag
                given_center = value_center
                other_subject = f"the maps of the individual {'groups' if grouped else 'samples'}"
            self.run.warning(
                f"'{given_flag}' centers the color scale of the {noun} layer on {given_center:g}, "
                f"but '{other_flag}' was not given, so the scale of {other_subject} is not "
                f"centered. Unless '{individual_colormap_flag}' gives the two scales different "
                f"colormaps, they use the same colormap. Its middle color then means "
                f"{given_center:g} on one map, and a different value on the other. The two scales "
                f"are centered separately because a summary can span much less than the values it "
                f"summarizes. Leaving one of them uncentered may be what you intend. If it is not, "
                f"center the other scale with '{other_flag}'."
            )

        # Per-sample values: each sample's rows reduced to one value per accession by the
        # within-sample aggregation. These are computed even when no context is colored by value, so
        # that a value column is validated whenever the file has one.
        sample_values: Dict[str, Dict[str, float]] = {s: {} for s in all_sample_names}
        for sample_name, sample_rows in df.groupby('__sample'):
            sample_values[sample_name] = self._aggregate_accession_quantities(
                sample_rows, value_column, gene_reduction, path, channel, undefined
            )
        if slice_by != 'gene':
            self._warn_undefined_values(undefined, gene_aggregation, path, channel)

        if 'quantitative' not in (model['unified_mode'], model['individual_mode']) and values_only:
            # A border or a bar that no map shows is refused once the maps of the run are known
            # ('_check_channel_maps'). The refusal says whether the samples are grouped.
            model['group_samples'] = group_samples
            return model
        if 'quantitative' not in (model['unified_mode'], model['individual_mode']):
            # Both summaries color presence, which only happens with groups (the per-sample maps of
            # an ungrouped run always show a value column's magnitude), so nothing is colored by
            # value.
            self.run.warning(
                f"The layer from the text file at '{path}' has a value column, '{value_column}', "
                f"but nothing on the maps is colored by it: with sample groups, the per-group maps "
                f"are colored by '--{flag_stem}-sample-summary' and the 'unified' map by "
                f"'--{flag_stem}-group-summary', and both of these summarize presence rather than "
                f"value. Set '--{flag_stem}-sample-summary' to an aggregation, such as "
                f"'{suggested_summary}', to color the group maps by value."
            )
            return model

        # A summary pools the values of a map element, not of an accession. An element's value in a
        # sample is the aggregate of its accessions in that sample. It is what that sample's own map
        # draws. One accession can belong to different elements on different maps, so these values
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
        # A channel without colormaps, such as a bar, sets lengths rather than colors. It has no
        # normalization to relabel its scale.
        if channel['default_colormap'] is None:
            model['cmap'] = None
            model['individual_cmap'] = None
            model['aggregate'] = aggregate
            model['colorbar_label'] = value_column
            model['individual_colorbar_label'] = value_column
            return model
        unified_default_colormap = (
            channel['default_colormap'] if unified_spread else default_colormap
        )
        unified_limits = _limits(
            channel['default_colormap_limits'] if unified_spread else default_colormap_limits
        )
        model['cmap'] = self._resolve_sequential_colormap(
            colormap if colormap is not None else unified_default_colormap, unified_limits,
            subject=channel['plural']
        )
        if unified_spread and self._is_cyclic_colormap(model['cmap']):
            raise ConfigError(
                f"'{colormap_flag}' was given as '{colormap}', a cyclic colormap, but "
                f"'{unified_summary_flag} {unified_summary}' colors the 'unified' map by how far "
                f"the values are spread apart. A spread does not repeat. A cyclic colormap would "
                f"give the smallest and the largest spread the same color. Please give a "
                f"sequential colormap, such as '{channel['default_colormap']}'."
            )
        # The maps of the individual samples or groups take a colormap of their own where one was
        # given, so that a 'unified' map showing a summary of another kind than the values behind it
        # — their spread rather than their magnitude — is not read as more of the same quantity.
        # Left unset, the layer's one colormap serves both contexts, and is resolved once so that a
        # warning about it is given once. The exception is a spread on the 'unified' map. The maps
        # of the individual samples or groups then show points in the period. Left unset, their
        # colormap is the cyclic default.
        individual_subject = (
            f"{channel['plural']} on the maps of the individual "
            f"{'groups' if grouped else 'samples'}"
        )
        if individual_colormap is not None:
            model['individual_cmap'] = self._resolve_sequential_colormap(
                individual_colormap, individual_colormap_limits, subject=individual_subject
            )
        elif unified_spread:
            model['individual_cmap'] = self._resolve_sequential_colormap(
                default_colormap, None, subject=individual_subject
            )
        else:
            model['individual_cmap'] = model['cmap']
        # Whether a centered scale's colormap has a middle worth putting a value at is asked of each
        # scale's own colormap, and only where that scale both colors by value and was centered. The
        # two questions are asked separately even where one colormap serves both, since a scale that
        # was not centered has nothing to ask.
        if model['unified_mode'] == 'quantitative':
            self._check_centered_colormap(
                model['cmap'], unified_limits, value_center, value_center_flag, colormap_flag,
                subject=f"{channel['plural']} on the 'unified' map"
            )
        if model['individual_mode'] == 'quantitative':
            self._check_centered_colormap(
                model['individual_cmap'],
                unified_limits if individual_colormap is None else individual_colormap_limits,
                individual_value_center, individual_value_center_flag,
                colormap_flag if individual_colormap is None else individual_colormap_flag,
                subject=individual_subject
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
        model['individual_colorbar_label'] = (
            value_column if element_normalization_label_template is None
            else element_normalization_label_template.format(value=value_column)
        )
        if element_normalization_label is not None:
            model['individual_colorbar_label'] = element_normalization_label
            # A border's keys add ' (border)' to their labels, but not to a label given here.
            model['individual_colorbar_label_given'] = True
        return model

    def _build_categorical_txt_model(
        self,
        name: str,
        data: dict,
        path: str,
        gene_aggregation: Union[str, None],
        accession_aggregation: Union[str, None],
        color: str,
        colormap: Union[bool, str, mcolors.Colormap, None],
        colormap_limits: Union[Tuple[float, float], None],
        membership_colors_txt: Union[str, None],
        reverse_overlay: bool,
        all_sample_names: Union[List[str], None],
        group_samples: Union[Dict[str, List[str]], None],
        sample_summary: Union[str, None],
        group_summary: Union[str, None],
        value_colors_txt: Union[str, None],
        seed: int,
        numeric_settings: Dict[str, object],
        channel: dict
    ) -> dict:
        """
        Build the coloring model of a layer whose values are categorical ('VALUE_TYPES').

        Each categorical value takes a color of its own ('_resolve_value_colors'). The order of the
        values sets their drawing priority. The first value is drawn on top, unless
        'reverse_overlay' puts the last value on top. Where several values meet, a rule chooses one
        ('CATEGORICAL_RULES'). 'gene_aggregation' chooses among the genes of an accession, and
        'accession_aggregation' among the accessions of a map element. Each is
        'DEFAULT_CATEGORICAL_RULE' where it is None.

        Without a 'sample' column, the layer colors every map by its values, so 'categorical' is
        both its 'unified_mode' and its 'individual_mode'. With a 'sample' column, the map of each
        sample shows that sample's values. A value in one sample need not mean the same thing in
        another sample, so by default the 'unified' map shows presence across the samples, as a
        layer without values does ('_resolve_summary'). It counts the samples, however few there
        are. A categorical rule given as 'sample_summary' says that the values mean the same thing
        in every sample. The rule then chooses one value per element across the samples, and the
        'unified' map shows it.

        With groups of samples, a rule given as 'sample_summary' chooses one value per element among
        the samples of each group. Each group's map shows it. Without that rule, each group's map
        counts the group's samples. A rule given as 'group_summary' then chooses one value per
        element across the groups for the 'unified' map. It needs a rule for the samples, since the
        groups otherwise have no values. Without it, the 'unified' map counts the groups. Without a
        rule for the samples, no map shows the values. They then take no colors, however many there
        are.

        Parameters
        ==========
        name : str
            Colorbar filename stem ('reactions'/'compounds').

        data : dict
            The '_read_element_txt' result for this layer.

        path : str
            Path to the layer's text file, used in messages.

        numeric_settings : Dict[str, object]
            The options that apply only to numeric values, keyed by flag. Each must be None.

        channel : dict
            The channel that the values color ('_txt_channel').

        Notes
        =====
        The other parameters are those of '_build_txt_model'. 'value_colors_txt' is a
        kegg-value-colors-txt file, or None for the default colors of the channel. 'seed' is the
        seed of random choices ('_choose_categorical_value'). Their keys start with the salt of the
        channel.

        Returns
        =======
        dict
            The layer model that '_map_elements' consumes.
        """
        element_type = data['element_type']
        flag_stem = channel['flag_stem']
        noun = channel['noun']
        df = data['df']
        has_sample = data['sample_names'] is not None
        type_flag = f'--{flag_stem}-value-type'
        colormap_flag = f'--{flag_stem}-colormap'
        sample_summary_flag = f'--{flag_stem}-sample-summary'
        common = self._txt_model_common(name, data, path, membership_colors_txt, channel)

        # Categorical values take colors of their own, not colors from a scale, so the options that
        # shape a scale of numbers have nothing to act on.
        for flag, setting in numeric_settings.items():
            if setting is not None:
                raise ConfigError(
                    f"'{flag}' works on numeric values. The values of the {noun} layer are "
                    f"categorical ('{type_flag} categorical'). They are names, and each takes a "
                    f"color of its own. Please drop '{flag}'."
                )
        # Ungrouped, the sample summary colors the 'unified' map. Grouped, it colors the map of each
        # group from the group's samples, and the group summary colors the 'unified' map from the
        # groups. Each is a rule, a presence name, or None for presence. A rule chooses one value
        # across the samples or groups. A rule says that a value means the same thing in each of
        # them. Presence on a group's map is always the count of its samples ('_group_map_colors').
        grouped = has_sample and group_samples is not None
        group_summary_flag = f'--{flag_stem}-group-summary'

        border = channel['role'] == 'border'

        def _resolve_categorical_summary(summary, flag, voters):
            if summary is None:
                return None, None
            if summary in CATEGORICAL_RULES:
                return summary, None
            # A border shows values only. Its fill already shows presence.
            if summary in SUMMARY_PRESENCE_SCHEMES and border:
                raise ConfigError(
                    f"'{flag}' was given as '{summary}', which summarizes presence across the "
                    f"{voters}. A border shows values only, and the fill already shows presence. "
                    f"Choose one border value across the {voters} with one of "
                    f"{CATEGORICAL_RULES_PHRASE}, or drop '{flag}'."
                )
            if summary in SUMMARY_PRESENCE_SCHEMES:
                return None, SUMMARY_PRESENCE_SCHEMES[summary]
            # A border takes no presence names, so its message suggests none.
            presence = "" if border else (
                f"Summarize presence across the {voters} with {SUMMARY_PRESENCE_PHRASE}. Or "
            )
            raise ConfigError(
                f"'{flag}' was given as '{summary}'. The values of the {noun} layer are "
                f"categorical ('{type_flag} categorical'). They cannot be added up or averaged "
                f"across the {voters}. {presence}{'choose' if presence else 'Choose'} one value "
                f"across the {voters} with one of {CATEGORICAL_RULES_PHRASE}. Choosing assumes "
                f"that a value means the same thing in every one of them."
            )

        sample_rule, sample_scheme = _resolve_categorical_summary(
            sample_summary, sample_summary_flag, 'samples'
        ) if has_sample else (None, None)
        group_rule, group_scheme = _resolve_categorical_summary(
            group_summary, group_summary_flag, 'groups'
        ) if grouped else (None, None)
        if group_rule is not None and sample_rule is None:
            raise ConfigError(
                f"'{group_summary_flag}' was given as '{group_summary}', which chooses one value "
                f"among the values of the groups. Each group's value comes from "
                f"'{sample_summary_flag}'. Without a rule, that option counts each group's "
                f"samples, so the groups have no values to choose among. Please also give "
                f"'{sample_summary_flag}' a rule, such as 'plurality_order'."
            )
        # Grouped, the values reach a map only through a rule among each group's samples. Without
        # one, a colors file of values would color nothing.
        if grouped and sample_rule is None and value_colors_txt is not None and border:
            raise ConfigError(
                f"'--{flag_stem}-value-colors' gives a color to each categorical border value. "
                f"Here the samples are grouped, and '{sample_summary_flag}' gives no rule. The "
                f"maps of the groups and the 'unified' map then show no border. Give "
                f"'{sample_summary_flag}' a rule, such as 'plurality_order', or drop the colors "
                f"file."
            )
        if grouped and sample_rule is None and value_colors_txt is not None:
            raise ConfigError(
                f"'--{flag_stem}-value-colors' gives a color to each categorical value. Here the "
                f"samples are grouped, and '{sample_summary_flag}' gives no rule. Each group's map "
                f"therefore counts the group's samples, and no map shows the values. Give "
                f"'{sample_summary_flag}' a rule, such as 'plurality_order', to color each group's "
                f"map by a value chosen among its samples. Or drop the colors file."
            )
        unified_rule = group_rule if grouped else sample_rule
        unified_scheme = group_scheme if grouped else sample_scheme
        unified_summary_flag = group_summary_flag if grouped else sample_summary_flag
        presence_unified = has_sample and unified_rule is None
        # The colormap of a border colors numbers only. A border shows no presence for it to color.
        if border and colormap is not None:
            raise ConfigError(
                f"'{colormap_flag}' colors numeric border values. The border values of the "
                f"{element_type} layer are categorical ('{type_flag} categorical'). Each takes a "
                f"color of its own, from '--{flag_stem}-value-colors' or from anvi'o's default "
                f"border colors. Please drop '{colormap_flag}'."
            )
        if not presence_unified and colormap is not None:
            reason = (
                f"'{unified_summary_flag} {unified_rule}' colors the 'unified' map by the values "
                f"too."
                if has_sample else
                f"The file at '{path}' has no 'sample' column, so its one map shows the values."
            )
            option = colormap_flag if colormap is not False else f'--{element_type}-color'
            raise ConfigError(
                f"'{option}' colors presence on the 'unified' map. Here no map of the {noun} layer "
                f"is colored by presence. {reason} Categorical values take their colors from "
                f"'--{flag_stem}-value-colors', or from anvi'o's default colors."
            )

        gene_rule = self._resolve_categorical_rule(
            gene_aggregation, f'--{flag_stem}-gene-aggregation', channel
        )
        accession_rule = self._resolve_categorical_rule(
            accession_aggregation, f'--{flag_stem}-accession-aggregation', channel
        )

        # The order of the values sets their colors, their place in the legend, and their drawing
        # priority. The first value is drawn on top. Its priority is 1. The priorities of the others
        # fall evenly towards 0, which no value reaches. 'reverse_overlay' reverses the priorities,
        # so the last value is drawn on top. Where no map shows the values, they take no colors. Any
        # number of them can then be read without a colors file. The rule among the genes of an
        # accession still needs their order ('_select_accession_values').
        if grouped and sample_rule is None:
            value_colors = None
            value_color_priorities = None
            value_order = {
                value: position
                for position, value in enumerate(self._sort_values(set(df['__value'])))
            }
        else:
            value_colors = self._resolve_value_colors(
                df['__value'], value_colors_txt, channel, path
            )
            value_order = {value: position for position, value in enumerate(value_colors)}
            value_count = len(value_colors)
            value_color_priorities = {
                value: (
                    color,
                    (position + 1) / value_count if reverse_overlay
                    else (value_count - position) / value_count
                )
                for position, (value, color) in enumerate(value_colors.items())
            }

        # The genes of an accession can have different values. A rule chooses one value for the
        # accession in each sample. The choices are counted for '_report_categorical_choices'.
        choices = {
            'gene': Counter(), 'element': Counter(), 'sample': Counter(), 'group': Counter(),
            'compound': Counter()
        }
        if has_sample:
            unified_values = None
            sample_values: Dict[str, Dict[str, str]] = {sample: {} for sample in all_sample_names}
            for sample, sample_rows in df.groupby('__sample'):
                sample_values[sample] = self._select_accession_values(
                    sample_rows, gene_rule, value_order, seed, sample, choices['gene'],
                    channel['choice_salt']
                )
            membership, source_accessions = self._relate_accessions_to_samples(
                df, all_sample_names
            )
        else:
            sample_values = None
            unified_values = self._select_accession_values(
                df, gene_rule, value_order, seed, '', choices['gene'], channel['choice_salt']
            )
            membership, source_accessions = {}, {}

        if not presence_unified:
            unified_mode = 'categorical'
        elif colormap is False:
            unified_mode = 'static'
        else:
            unified_mode = 'membership'
        # A group's map shows a value only where a rule chooses one among the group's samples.
        # Otherwise it counts the group's samples. A border shows no presence, so a map that the
        # fill would color by presence shows no border, which the mode None says.
        individual_mode = 'membership' if grouped and sample_rule is None else 'categorical'
        if border:
            unified_mode = None if presence_unified else unified_mode
            individual_mode = None if individual_mode == 'membership' else individual_mode
        model = {
            **common,
            'aggregations': {'gene': gene_rule, 'accession': accession_rule},
            'value_type': 'categorical',
            'unified_mode': unified_mode,
            'individual_mode': individual_mode,
            'membership': membership,
            'source_accessions': source_accessions,
            'colormap': True if colormap is None else colormap,
            'colormap_limits': colormap_limits,
            'colormap_scheme': unified_scheme,
            # Presence across the samples is counted by default, however few samples there are.
            # Coloring by membership would take its colors from 'tab10', as the values do. A color
            # would then mean a sample on the 'unified' map and a value on the maps of the samples.
            'default_colormap_scheme': 'by_count',
            'scheme_options': {
                scheme: f'{unified_summary_flag} {summary}'
                for summary, scheme in SUMMARY_PRESENCE_SCHEMES.items()
            },
            'reverse_overlay': reverse_overlay,
            'unified_values': unified_values,
            'sample_values': sample_values,
            'group_samples': group_samples if grouped else None,
            'value_colors': value_colors,
            'value_order': value_order,
            'value_color_priorities': value_color_priorities,
            'gene_rule': gene_rule,
            'accession_rule': accession_rule,
            'sample_rule': sample_rule,
            'group_rule': group_rule,
            'seed': seed,
            'choice_salt': channel['choice_salt'],
            'colorbar_label': data['value_column'],
            '_element_cache': {},
            '_choices': choices
        }
        # The presence color of a static 'unified' map. A layer without samples has no presence
        # color, as a layer of numeric values without samples has none.
        if has_sample:
            model['color_hexcode'] = color
        return model

    @staticmethod
    def _border_flag(element_type: str, option: str) -> str:
        """
        Name the command-line flag of an option of the border of a layer ('BORDER_OPTIONS').

        Parameters
        ==========
        element_type : str
            'reaction' or 'compound'.

        option : str
            A key of 'BORDER_OPTIONS'.

        Returns
        =======
        str
            The flag, such as '--reaction-border-colormap'. Limits of a colormap and the label of a
            normalization are given with the flag of their colormap or normalization.
        """
        option = {
            'colormap_limits': 'colormap',
            'individual_colormap_limits': 'individual_colormap',
            'element_normalization_label': 'element_normalization'
        }.get(option, option)
        return f"--{element_type}-border-{option.replace('_', '-')}"

    @staticmethod
    def _check_border_options(
        reaction_txt: Union[str, None],
        compound_txt: Union[str, None],
        reaction_border: Union[Dict[str, object], None],
        compound_border: Union[Dict[str, object], None],
        reaction_color: Union[str, None]
    ) -> None:
        """
        Check the options of the borders of the two layers, as far as that can be done without the
        files.

        Every option of a border needs the border column of its layer, and that column needs the
        layer's file. The builder of the border's model checks the rest ('_build_txt_model'). The
        command line makes these checks too, and these repeat them for programmatic callers. Where a
        border column is named, the private methods of Biopython that draw rings are checked too
        ('kgml.SlicedKGMLCanvas.check_biopython').

        Parameters
        ==========
        reaction_txt, compound_txt : Union[str, None]
            Paths to the two files, or None.

        reaction_border, compound_border : Union[Dict[str, object], None]
            The options of the two borders, keyed by 'BORDER_OPTIONS', or None.

        reaction_color : Union[str, None]
            The presence color of the reaction layer. It is 'original' for '--original-color'.
        """
        bordered = False
        for element_type, txt, border in (
            ('reaction', reaction_txt, reaction_border), ('compound', compound_txt, compound_border)
        ):
            if not border:
                continue
            allowed = [
                option for option in BORDER_OPTIONS
                if element_type == 'reaction' or option not in REACTION_BORDER_OPTIONS
            ]
            unknown = [option for option in border if option not in allowed]
            if unknown:
                raise ConfigError(
                    f"The options of the {element_type} border include "
                    f"{', '.join(repr(option) for option in unknown)}, which "
                    f"{'is not an option' if len(unknown) == 1 else 'are not options'} of a "
                    f"{element_type} border. It takes these: "
                    f"{', '.join(repr(option) for option in allowed)}."
                )
            column_flag = Mapper._border_flag(element_type, 'column')
            # A colormap and its limits, or a normalization and its label, come with one flag.
            given = list(dict.fromkeys(
                Mapper._border_flag(element_type, option) for option in BORDER_OPTIONS
                if option != 'column' and border.get(option) is not None
            ))
            if border.get('column') is None:
                if given:
                    listed = ', '.join(f"'{flag}'" for flag in given)
                    raise ConfigError(
                        f"{listed} {'apply' if len(given) > 1 else 'applies'} to the borders of "
                        f"the {element_type} layer, but no border column was named with "
                        f"'{column_flag}'. Name the column whose values color the borders."
                    )
                continue
            bordered = True
            if txt is None:
                raise ConfigError(
                    f"'{column_flag}' names the border column of the {element_type} file, but no "
                    f"{element_type} file ('--{element_type}-txt') was provided."
                )
            value_type = border.get('value_type')
            type_flag = Mapper._border_flag(element_type, 'value_type')
            if value_type is not None and value_type not in VALUE_TYPES:
                raise ConfigError(
                    f"'{type_flag}' was given as '{value_type}'. It can be "
                    f"{' or '.join(repr(name) for name in VALUE_TYPES)}."
                )
            if border.get('value_colors') is not None and value_type != 'categorical':
                raise ConfigError(
                    f"'{Mapper._border_flag(element_type, 'value_colors')}' gives a color to each "
                    f"categorical border value of the {element_type} layer. The border values of "
                    f"that layer are numeric. To read them as categorical values, give "
                    f"'{type_flag} categorical'."
                )
            # A ring sliced by gene gives each gene a slice of its own, so no genes are reduced to
            # an accession. Slices of rings are not normalized, as slices of fills are not.
            slice_by = border.get('slice_by')
            slice_flag = Mapper._border_flag(element_type, 'slice_by')
            if slice_by is not None and slice_by not in SLICE_UNITS:
                raise ConfigError(
                    f"'{slice_flag}' was given as '{slice_by}'. It can be "
                    f"{' or '.join(repr(unit) for unit in SLICE_UNITS)}."
                )
            if slice_by == 'gene' and border.get('gene_aggregation') is not None:
                gene_flag = Mapper._border_flag(element_type, 'gene_aggregation')
                raise ConfigError(
                    f"'{gene_flag}' reduces the border values of the genes of an accession to one "
                    f"value. '{slice_flag} gene' gives each gene its own slice of a ring, so the "
                    f"genes of an accession are never reduced. Please drop '{gene_flag}'."
                )
            if slice_by is not None and border.get('element_normalization') is not None:
                raise ConfigError(
                    f"'{slice_flag}' cannot yet be used with "
                    f"'{Mapper._border_flag(element_type, 'element_normalization')}'. The border "
                    f"values of slices cannot be normalized in this version. Please use only one "
                    f"of the two options."
                )
            # The reference map's colors are drawn without rings.
            if reaction_color == 'original':
                raise ConfigError(
                    f"'{column_flag}' colors the borders of {element_type} elements by a column of "
                    f"values, but '--original-color' takes both the colors and the drawing order "
                    f"of the reactions from the reference map, and it draws no borders. Please use "
                    f"only one of them."
                )
        # Rings are drawn with private methods of Biopython, which a release can rename. This is
        # checked before anything is written.
        if bordered:
            kgml.SlicedKGMLCanvas.check_biopython()

    @staticmethod
    def _check_bar_options(
        reaction_txt: Union[str, None],
        reaction_bar: Union[Dict[str, object], None],
        reaction_color: Union[str, None],
        slice_by: Union[str, None]
    ) -> None:
        """
        Check the options of the bars of the reaction layer, as far as that can be done without the
        file.

        Every option of the bars needs the bar column, and that column needs the reaction file. A
        bar is drawn in a box whose fill is colored, in the color of the fill unless it is a strip
        ('BAR_STYLES'), so fills left uncolored are refused. The reference map's colors are refused
        too, since a renderer of their own draws them, and it draws no bars. Bars in the slices of
        boxes are not supported yet. The builder of the model of the bars checks the rest
        ('_build_txt_model'). The command line makes these checks too, and these repeat them for
        programmatic callers. Where a bar column is named, the private methods of Biopython that
        draw bars are checked too ('kgml.SlicedKGMLCanvas.check_biopython').

        Parameters
        ==========
        reaction_txt : Union[str, None]
            Path to the reaction file, or None.

        reaction_bar : Union[Dict[str, object], None]
            The options of the bars, keyed by 'BAR_OPTIONS', or None.

        reaction_color : Union[str, None]
            The presence color of the reaction layer. It is 'original' for '--original-color', and
            'none' where the fills are left uncolored.

        slice_by : Union[str, None]
            What the reaction elements are sliced by ('SLICE_UNITS'), or None.
        """
        if not reaction_bar:
            return
        unknown = [option for option in reaction_bar if option not in BAR_OPTIONS]
        if unknown:
            raise ConfigError(
                f"The options of the reaction bars include "
                f"{', '.join(repr(option) for option in unknown)}, which "
                f"{'is not an option' if len(unknown) == 1 else 'are not options'} of reaction "
                f"bars. They take these: {', '.join(repr(option) for option in BAR_OPTIONS)}."
            )
        column_flag = '--reaction-bar-column'
        given = [
            f"'--reaction-bar-{option.replace('_', '-')}'" for option in BAR_OPTIONS
            if option != 'column' and reaction_bar.get(option) is not None
        ]
        if reaction_bar.get('column') is None:
            if given:
                raise ConfigError(
                    f"{', '.join(given)} {'apply' if len(given) > 1 else 'applies'} to the bars of "
                    f"the reaction layer, but no bar column was named with '{column_flag}'. Name "
                    f"the column whose values set the lengths of the bars."
                )
            return
        if reaction_txt is None:
            raise ConfigError(
                f"'{column_flag}' names the bar column of the reaction file, but no reaction file "
                f"('--reaction-txt') was provided."
            )
        style = reaction_bar.get('style')
        if style is not None and style not in BAR_STYLES:
            raise ConfigError(
                f"'--reaction-bar-style' takes {' or '.join(repr(name) for name in BAR_STYLES)}, "
                f"not '{style}'."
            )
        # A strip is dark on any fill, but it is drawn only in a colored box, as a bar of the fill
        # style is.
        if style == 'strip':
            draws = "draws bars as strips in colored reaction boxes"
            uncolored = "so no box would show a bar"
        else:
            draws = "draws bars in the colors of the fills of reaction boxes"
            uncolored = "so the bars would have no color"
        if reaction_color == 'original':
            raise ConfigError(
                f"'{column_flag}' {draws}, but '--original-color' takes both the colors and the "
                f"drawing order of the reactions from the reference map, and it draws no bars. "
                f"Please use only one of them."
            )
        if reaction_color == 'none':
            raise ConfigError(
                f"'{column_flag}' {draws}, but '--reaction-color none' leaves the fills uncolored, "
                f"{uncolored}. Please give the fills a color, or drop one of the two options."
            )
        if slice_by is not None:
            raise ConfigError(
                f"'--reaction-slice-by' cannot be combined with '{column_flag}' yet. Bars in the "
                f"slices of boxes are not supported. Please use only one of them."
            )
        # Bars are drawn with private methods of Biopython, which a release can rename. This is
        # checked before anything is written.
        kgml.SlicedKGMLCanvas.check_biopython()

    @staticmethod
    def _check_slice_options(
        reaction_txt: Union[str, None],
        slice_by: Union[str, None],
        weights: Union[str, None],
        weight_aggregation: Union[str, None],
        share_summary: Union[str, None],
        element_normalization: Union[str, None],
        color: str,
        colormap: Union[bool, str, mcolors.Colormap, None],
        gene_aggregation: Union[str, None] = None,
        border_column: Union[str, None] = None,
        border_slice_by: Union[str, None] = None
    ) -> None:
        """
        Check the options that slice reaction elements, as far as that can be done without the file.

        '_build_slice_model' checks the rest once the file is read and the layer model is built. The
        reader checks the 'gene_id' column that slices of genes need ('_read_element_txt'). The
        command line makes the same checks, and these repeat them for programmatic callers. Where
        slicing is on, the private methods of Biopython that draw slices are checked too
        ('kgml.SlicedKGMLCanvas.check_biopython').

        The weights size the slices of elements and the slices of rings that are sliced on their own
        ('--reaction-border-slice-by'). The options of weights therefore need either kind of
        slicing. Their aggregation needs a kind that slices by accession.

        Parameters
        ==========
        reaction_txt : Union[str, None]
            Path to the reaction file, or None.

        border_column : Union[str, None], None
            The border column of the reaction file, or None.

        border_slice_by : Union[str, None], None
            What the rings of the reaction border are sliced by, or None.

        Notes
        =====
        The other parameters are the arguments of 'map_kegg_pathways_txt' of the same names, without
        the 'reaction_' prefix. 'color' is 'original' for '--original-color', and 'none' where the
        fills are left uncolored. 'colormap' is False where a color was chosen with
        '--reaction-color'.
        """
        given = [
            flag for flag, setting in (
                ('--reaction-slice-by', slice_by),
                ('--reaction-slice-weights', weights),
                ('--reaction-slice-weight-aggregation', weight_aggregation),
                ('--reaction-slice-share-summary', share_summary)
            ) if setting is not None
        ]
        if not given:
            return
        listed = ', '.join(f"'{flag}'" for flag in given)
        plural = len(given) > 1
        if reaction_txt is None:
            raise ConfigError(
                f"{listed} {'apply' if plural else 'applies'} to the slices of reaction elements, "
                f"but no reaction file ('--reaction-txt') was provided."
            )
        if slice_by is None and border_slice_by is None:
            # Weights of genes are reduced to an accession only when slicing by accession.
            units = ('accession',) if weight_aggregation is not None else SLICE_UNITS
            suggestion = ' or '.join(f"'--reaction-slice-by {unit}'" for unit in units)
            drop = 'these options' if plural else 'it'
            if border_column is None:
                raise ConfigError(
                    f"{listed} {'apply' if plural else 'applies'} to the slices of reaction "
                    f"elements. No element is sliced here, because '--reaction-slice-by' was not "
                    f"given. To slice elements, give {suggestion}. Otherwise, drop {drop}."
                )
            ring_suggestion = ' or '.join(f"'--reaction-border-slice-by {unit}'" for unit in units)
            raise ConfigError(
                f"{listed} {'apply' if plural else 'applies'} to the slices of reaction elements "
                f"and of their rings. Nothing is sliced here, because neither "
                f"'--reaction-slice-by' nor '--reaction-border-slice-by' was given. To slice "
                f"elements, give {suggestion}. To slice the rings of the border column "
                f"'{border_column}', give {ring_suggestion}. Otherwise, drop {drop}."
            )
        for flag, setting, choices in (
            ('--reaction-slice-by', slice_by, SLICE_UNITS),
            ('--reaction-slice-weight-aggregation', weight_aggregation, SLICE_WEIGHT_AGGREGATIONS),
            ('--reaction-slice-share-summary', share_summary, SLICE_SHARE_SUMMARIES)
        ):
            if setting is not None and setting not in choices:
                raise ConfigError(
                    f"'{flag}' was given as '{setting}'. It can be "
                    f"{' or '.join(repr(choice) for choice in choices)}."
                )
        # Each gene is a slice of its own, so no genes are reduced to an accession. The weights of
        # genes are reduced where the elements or the rings are sliced by accession.
        if slice_by == 'gene' and gene_aggregation is not None:
            raise ConfigError(
                "'--reaction-gene-aggregation' reduces the values of the genes of an accession to "
                "one value. '--reaction-slice-by gene' gives each gene a slice of its own, so the "
                "genes of an accession are never reduced. Please drop "
                "'--reaction-gene-aggregation'."
            )
        if weight_aggregation is not None and 'accession' not in (slice_by, border_slice_by):
            by_gene = [
                f"'{flag} gene'" for flag, unit in (
                    ('--reaction-slice-by', slice_by),
                    ('--reaction-border-slice-by', border_slice_by)
                ) if unit == 'gene'
            ]
            raise ConfigError(
                f"'--reaction-slice-weight-aggregation' reduces the weights of the genes of an "
                f"accession to one weight. {' and '.join(by_gene)} "
                f"{'give' if len(by_gene) > 1 else 'gives'} each gene a slice of its own, so the "
                f"genes of an accession are never reduced. Please drop "
                f"'--reaction-slice-weight-aggregation'."
            )
        if weights is None:
            for flag, setting, action in (
                (
                    '--reaction-slice-weight-aggregation', weight_aggregation,
                    "reduces the weights of the genes of an accession to one weight"
                ),
                (
                    '--reaction-slice-share-summary', share_summary,
                    "pools the shares of slices across samples or groups"
                )
            ):
                if setting is not None:
                    raise ConfigError(
                        f"'{flag}' {action}. No column of weights was named with "
                        f"'--reaction-slice-weights', so every slice on a map takes the same "
                        f"share. To size slices by weights, name their column. Otherwise, drop "
                        f"'{flag}'."
                    )
        # The rest concerns the slices of elements. Rings sliced on their own need none of it.
        if slice_by is None:
            kgml.SlicedKGMLCanvas.check_biopython()
            return
        if element_normalization is not None:
            raise ConfigError(
                "'--reaction-slice-by' cannot yet be used with '--reaction-element-normalization'. "
                "The values of slices cannot be normalized in this version. Please use only one of "
                "the two options."
            )
        if color == 'original':
            raise ConfigError(
                "'--reaction-slice-by' colors each slice of a reaction element by its own value, "
                "but '--original-color' takes the colors of the reactions from the reference map. "
                "Please use only one of them."
            )
        if color == 'none':
            raise ConfigError(
                "'--reaction-slice-by' colors each slice of a reaction element by its own value, "
                "but '--reaction-color none' leaves the fills of the reaction layer uncolored. To "
                "slice the rings of the border column, give '--reaction-border-slice-by' in place "
                "of '--reaction-slice-by'."
            )
        if colormap is False:
            raise ConfigError(
                "'--reaction-slice-by' colors each slice of a reaction element by its own value, "
                "but '--reaction-color' gives the reaction layer a single color. Please use only "
                "one of them."
            )
        # Slices are drawn with private methods of Biopython, which a release can rename. This is
        # checked before anything is written.
        kgml.SlicedKGMLCanvas.check_biopython()

    def _build_slice_model(
        self,
        model: dict,
        data: dict,
        slice_by: str,
        weight_aggregation: Union[str, None],
        share_summary: Union[str, None],
        accession_aggregation: Union[str, None],
        all_sample_names: Union[List[str], None],
        group_samples: Union[Dict[str, List[str]], None],
        draw_unified_maps: bool,
        draw_individual_maps: bool,
        fill: Union[dict, None] = None
    ) -> dict:
        """
        Check that a reaction layer can be sliced, and build the keys of its model for slicing.

        Slices are drawn only on the maps that show the values of the layer, which are those whose
        mode is 'quantitative' or 'categorical'. Slicing is refused where no drawn map shows the
        values. An aggregation of the accessions of an element is refused for numeric values sliced
        by accession. For numeric values sliced by gene, it is refused if it does not return a
        single value unchanged. A summary of shares needs a drawn map that pools the samples and
        shows slices. Such maps are the maps of groups and the 'unified' map of a layer with
        samples. Whether one exists is checked once the fill and the border are both built
        ('map_kegg_pathways_txt').

        With slices of accessions, the weights of the genes of an accession are reduced to one
        weight per accession and sample, or per accession for a layer without samples. Slices of
        genes keep the values and weights of each row ('_build_gene_slice_keys').

        The rings of the border of a layer are sliced by this method too
        ('--reaction-border-slice-by'). The model of the border is built from the rows with a border
        value ('_channel_view'), and the messages name the options of the border. A border that no
        drawn map shows is refused with the maps of the run ('_check_channel_maps'). Where the fill
        is sliced by the same unit, the rings are aligned with its columns. The fill then gives the
        slices and their shares ('_ring_result'), so the border takes no weights. With slices of
        genes, it takes the genes of the fill, so that both find the same slices. Otherwise, the
        rings are sliced on their own. Their shares then come from the weights of the rows with a
        border value.

        Parameters
        ==========
        model : dict
            The layer model built from the reaction file, or the model of its border
            ('_build_txt_model').

        data : dict
            The '_read_element_txt' result for the reaction file, or its border view
            ('_channel_view').

        slice_by : str
            One of 'SLICE_UNITS'.

        weight_aggregation : Union[str, None]
            One of 'SLICE_WEIGHT_AGGREGATIONS', or None for 'sum'.

        share_summary : Union[str, None]
            One of 'SLICE_SHARE_SUMMARIES', or None for 'share_mean'.

        accession_aggregation : Union[str, None]
            The aggregation or rule given for the accessions of an element, or None.

        all_sample_names : Union[List[str], None]
            The samples of the run, or None for a run without samples.

        group_samples : Union[Dict[str, List[str]], None]
            The samples of each group, or None for a run without groups.

        draw_unified_maps : bool
            True if the 'unified' map is drawn.

        draw_individual_maps : bool
            True if the maps of individual samples or groups are drawn.

        fill : Union[dict, None], None
            The sliced model of the fill, for a border whose rings are aligned with it. It is None
            for the fill, and for a border whose rings are sliced on their own.

        Returns
        =======
        dict
            The keys that '_map_elements' reads for slicing:
            - 'slice_by': the 'slice_by' argument.
            - 'slice_weights': True if the shares of slices come from a column of weights.
            - 'slice_share_summary': how shares are pooled, or None without weights.
            - 'unified_weights': the weight of each accession, for a layer without samples. It is
              None otherwise.
            - 'sample_weights': the weight of each accession in each sample, for a layer with
              samples. It is None otherwise.
            - 'slice_pooled': True if a drawn map pools samples and shows the values.
            - '_slice_cache': the numeric values of slices ('_slice_values').
            - '_slice_results': the slices of each element on each map ('_slice_result').
            - '_slice_key' and '_slice_uncolored': what the drawn maps record for the key table
              and the warnings ('_write_slice_keys').
            A border adds 'slice_layout', which is 'aligned' or 'independent', and '_slice_unshown',
            where the drawn maps record the units with border values that no column of the fill
            shows. Slices of genes add the keys of '_build_gene_slice_keys'.
        """
        path = model['path']
        df = data['df']
        channel = model['channel']
        border = channel['role'] == 'border'
        flag_stem = channel['flag_stem']
        slice_flag = f'--{flag_stem}-slice-by'
        aggregation_flag = f'--{flag_stem}-accession-aggregation'
        # The view of a border always names its border column as its value column.
        if data['value_column'] is None:
            raise ConfigError(
                f"'--reaction-slice-by' colors each slice of a reaction element by its own value, "
                f"but the reaction file at '{path}' has no value column. Its layer is colored by "
                f"presence, and slices of presence are not supported yet. Add a value column to "
                f"the file, or drop '--reaction-slice-by'."
            )
        if border:
            each_slice, value, whole = 'each slice of a ring', 'border value', 'A ring drawn whole'
        else:
            each_slice, value, whole = 'each slice', 'value', 'An element drawn whole'
        # A slice of a gene can stand for rows of several accessions of an element. The accession
        # aggregation reduces these rows, so it is accepted there.
        if (
            accession_aggregation is not None and data['value_type'] != 'categorical'
            and slice_by == 'accession'
        ):
            raise ConfigError(
                f"'{aggregation_flag}' reduces the {channel['values']} of the accessions of a map "
                f"element to one value. With '{slice_flag} accession', {each_slice} shows the "
                f"{value} of one accession. {whole} takes the mean of its slice values, weighted "
                f"by their shares, in place of this reduction. Please drop '{aggregation_flag}'."
            )
        # A gene with one row for an element shows the value of that row ('_build_gene_slice_keys').
        # An aggregation that does not return a single value unchanged, such as 'std' or 'count',
        # would give genes with several rows values of another kind on the same scale. With a
        # period, only 'circular_mean' is accepted, and it keeps a single value.
        if (
            accession_aggregation is not None and data['value_type'] != 'categorical'
            and slice_by == 'gene' and model['value_period'] is None
        ):
            try:
                single = float(model['aggregate']([3.0]))
            except (TypeError, ValueError):
                single = float('nan')
            if not (np.isfinite(single) and abs(single - 3.0) < 1e-9):
                raise ConfigError(
                    f"'{aggregation_flag}' was given as '{accession_aggregation}'. With "
                    f"'{slice_flag} gene', it reduces the rows of a gene for the accessions of an "
                    f"element to one value. A gene with one such row shows the {value} of that "
                    f"row. '{accession_aggregation}' does not return a single value unchanged, so "
                    f"genes with several rows would show another kind of value on the same color "
                    f"scale. Please give 'mean', 'sum', 'max', 'min', or 'median' instead."
                )
        # Rings aligned with the fill take the shares of the fill, so only the fill reads weights
        # then. The weights of genes are reduced only for slices of accessions.
        weighted = data['slice_weights'] is not None and fill is None
        if (
            weight_aggregation is not None and weighted and slice_by == 'accession'
            and 'gene_id' not in df.columns
        ):
            raise ConfigError(
                f"'--reaction-slice-weight-aggregation' reduces the weights of the genes of an "
                f"accession to one weight. The reaction file at '{path}' has no 'gene_id' column, "
                f"so there are no weights of genes to reduce. Please drop "
                f"'--reaction-slice-weight-aggregation'."
            )

        value_modes = ('quantitative', 'categorical')
        has_sample = data['sample_names'] is not None
        grouped = has_sample and group_samples is not None
        unified_values = draw_unified_maps and model['unified_mode'] in value_modes
        individual_values = draw_individual_maps and model['individual_mode'] in value_modes
        if not border and not (unified_values or individual_values):
            noun = 'groups' if grouped else 'samples'
            suggestion = (
                "a rule, such as 'plurality_order'" if model.get('value_type') == 'categorical'
                else "an aggregation, such as "
                f"'{'mean' if model.get('value_period') is None else 'circular_mean'}'"
            )
            parts = []
            if model['unified_mode'] in value_modes:
                parts.append("The 'unified' map would show them, but it is not drawn here.")
            else:
                # Grouped, the 'unified' map pools the values of the groups, which come from the
                # sample summary.
                flags = (
                    "'--reaction-sample-summary' and '--reaction-group-summary' both"
                    if grouped and model['individual_mode'] not in value_modes
                    else f"'--reaction-{'group' if grouped else 'sample'}-summary'"
                )
                parts.append(
                    f"The 'unified' map shows presence across the {noun}. To show the values "
                    f"there, set {flags} to {suggestion}."
                )
            if has_sample and model['individual_mode'] in value_modes:
                parts.append(
                    f"The maps of the {noun} would show them, but none are drawn here. Draw them "
                    f"with '--draw-individual-files' or '--draw-grid'."
                )
            elif grouped:
                parts.append(
                    f"The maps of the groups count the samples of each group. To show the values "
                    f"there, set '--reaction-sample-summary' to {suggestion}."
                )
            raise ConfigError(
                f"'--reaction-slice-by' slices the reaction elements of the maps that show the "
                f"values of the reaction layer. No drawn map shows them here. {' '.join(parts)}"
            )

        keys = {
            'slice_by': slice_by,
            'slice_weights': weighted,
            'slice_share_summary': None,
            'unified_weights': None,
            'sample_weights': None,
            'slice_pooled': has_sample and (unified_values or (grouped and individual_values)),
            '_slice_cache': {},
            '_slice_results': {},
            '_slice_key': {},
            '_slice_uncolored': {}
        }
        if border:
            keys['slice_layout'] = 'independent' if fill is None else 'aligned'
            keys['_slice_unshown'] = {}
        if weighted:
            keys['slice_share_summary'] = (
                SLICE_SHARE_SUMMARIES[0] if share_summary is None else share_summary
            )
        if slice_by == 'gene':
            keys.update(self._build_gene_slice_keys(
                model, data if weighted else {**data, 'slice_weights': None},
                accession_aggregation, all_sample_names
            ))
            # Aligned rings find the genes of each element as the fill does. Each slice of a gene
            # then takes its border value from the rows of the gene that have one.
            if fill is not None:
                keys['accession_genes'] = fill['accession_genes']
            return keys
        if not weighted:
            return keys
        aggregation = (
            SLICE_WEIGHT_AGGREGATIONS[0] if weight_aggregation is None else weight_aggregation
        )

        def _accession_weights(rows: pd.DataFrame) -> Dict[str, float]:
            return {
                accession: float(weight) for accession, weight
                in rows.groupby('__accession')['__slice_weight'].agg(aggregation).items()
            }

        if has_sample:
            sample_weights: Dict[str, Dict[str, float]] = {
                sample: {} for sample in all_sample_names
            }
            for sample, sample_rows in df.groupby('__sample'):
                sample_weights[sample] = _accession_weights(sample_rows)
            keys['sample_weights'] = sample_weights
        else:
            keys['unified_weights'] = _accession_weights(df)
        return keys

    def _build_gene_slice_keys(
        self,
        model: dict,
        data: dict,
        accession_aggregation: Union[str, None],
        all_sample_names: Union[List[str], None]
    ) -> dict:
        """
        Build the keys of a layer model for slices of genes.

        One gene ID names the same gene in every sample. A slice of a gene stands for the rows of
        the gene for the accessions of an element. Each row is a member, named by its gene and its
        accession ('_gene_members'). Members keep the values and weights of their rows. A view of
        the layer puts them in place of accessions ('_gene_view').

        A numeric slice reduces the rows of its gene in a sample by the accession aggregation. Where
        none was given, the rows are reduced by 'mean', or by 'circular_mean' with a period. Rows of
        one gene usually repeat one measurement, so a sum would count it more than once. A gene with
        one row takes the value of that row. An aggregation that does not return a single value
        unchanged, such as 'std', is refused ('_build_slice_model'). A categorical slice takes the
        value that the accession rule chooses among the rows.

        Parameters
        ==========
        model : dict
            The layer model built from the reaction file ('_build_txt_model').

        data : dict
            The '_read_element_txt' result for the reaction file, read for slices of genes.

        accession_aggregation : Union[str, None]
            The aggregation or rule given for the accessions of an element, or None.

        all_sample_names : Union[List[str], None]
            The samples of the run, or None for a run without samples.

        Returns
        =======
        dict
            These keys:
            - 'accession_genes': the genes that have a row for each accession, in natural sort order
              ('_sort_values').
            - 'gene_sample_values': the value of each member in each sample, for a layer with
              samples. It is None otherwise.
            - 'gene_unified_values': the value of each member, for a layer without samples. It is
              None otherwise.
            - 'gene_sample_weights' and 'gene_unified_weights': the weights of the members in the
              same form. They are None without a column of weights.
            - 'gene_sample_accessions': the accessions of the rows of each gene in each sample in
              which it has rows ('_warn_gene_ids'). It is None without samples.
            - 'slice_accessions': the accessions with rows, sorted, in each sample, or in the file
              for a layer without samples. They find the entries to slice ('_map_elements').
            - 'gene_aggregate': reduces the values of the rows of a gene in a sample to one value.
              It is None for categorical values.
            - '_gene_elements', '_gene_cache', and '_gene_standin_cache': the slices that each
              element can have ('_gene_element'), the values of slices ('_slice_values'), and the
              categorical stand-ins of elements ('_gene_standin_individuals').
        """
        df = data['df']
        categorical = data['value_type'] == 'categorical'
        genes = df['__gene_id']
        members = (genes + '\t' + df['__accession']).tolist()
        if categorical:
            values = df['__value'].tolist()
        else:
            # The values were checked when the model was built ('_aggregate_accession_quantities').
            values = pd.to_numeric(df[data['value_column']], errors='coerce').astype(float).tolist()
        weighted = data['slice_weights'] is not None
        weights = df['__slice_weight'].tolist() if weighted else [None] * len(df)

        # The genes of each accession are kept in natural sort order. The slices of an element
        # follow the KGML order of their genes' first accessions, and then this order
        # ('_gene_element').
        ranks = {gene: rank for rank, gene in enumerate(self._sort_values(set(genes)))}
        accession_genes: Dict[str, Set[str]] = {}
        for accession, gene in zip(df['__accession'], genes):
            accession_genes.setdefault(accession, set()).add(gene)
        keys = {
            'accession_genes': {
                accession: tuple(sorted(gene_set, key=ranks.__getitem__))
                for accession, gene_set in accession_genes.items()
            },
            'gene_sample_values': None,
            'gene_unified_values': None,
            'gene_sample_weights': None,
            'gene_unified_weights': None,
            'gene_sample_accessions': None,
            '_gene_elements': {},
            '_gene_cache': {},
            '_gene_standin_cache': {}
        }
        if data['sample_names'] is None:
            keys['gene_unified_values'] = dict(zip(members, values))
            if weighted:
                keys['gene_unified_weights'] = dict(zip(members, weights))
            keys['slice_accessions'] = dict.fromkeys(sorted(accession_genes))
        else:
            sample_values: Dict[str, Dict[str, Union[float, str]]] = {
                sample: {} for sample in all_sample_names
            }
            sample_weights: Dict[str, Dict[str, float]] = {
                sample: {} for sample in all_sample_names
            }
            sample_accessions: Dict[str, Set[str]] = {sample: set() for sample in all_sample_names}
            gene_sample_accessions: Dict[str, Dict[str, Set[str]]] = {}
            for sample, gene, accession, member, value, weight in zip(
                df['__sample'], genes, df['__accession'], members, values, weights
            ):
                sample_values[sample][member] = value
                sample_weights[sample][member] = weight
                sample_accessions[sample].add(accession)
                gene_sample_accessions.setdefault(gene, {}).setdefault(sample, set()).add(accession)
            keys['gene_sample_values'] = sample_values
            if weighted:
                keys['gene_sample_weights'] = sample_weights
            keys['gene_sample_accessions'] = gene_sample_accessions
            keys['slice_accessions'] = {
                sample: dict.fromkeys(sorted(sample_accessions[sample]))
                for sample in all_sample_names
            }

        if categorical:
            keys['gene_aggregate'] = None
            return keys
        period = model['value_period']
        flag_stem = model['channel']['flag_stem']
        aggregate = model['aggregate'] if accession_aggregation is not None else (
            self._resolve_aggregation(
                'circular_mean' if period is not None else 'mean',
                f'--{flag_stem}-accession-aggregation', model['channel'], period,
                f'--{flag_stem}-value-period'
            )
        )

        def reduce_rows(row_values: List[float]) -> float:
            # A gene with one row takes the value of that row. With a period, the circular mean
            # gives that value as a point in the period, from 0 up to the period.
            if period is None and len(row_values) == 1:
                return float(row_values[0])
            return aggregate(row_values)
        keys['gene_aggregate'] = reduce_rows
        return keys

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
        reaction_value_type: Literal['numeric', 'categorical'] = 'numeric',
        compound_value_type: Literal['numeric', 'categorical'] = 'numeric',
        reaction_value_colors: str = None,
        compound_value_colors: str = None,
        seed: int = None,
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
        reaction_individual_colormap: Union[str, mcolors.Colormap] = None,
        reaction_individual_colormap_limits: Tuple[float, float] = None,
        reaction_membership_colors: str = None,
        reaction_reverse_overlay: bool = False,
        reaction_value_limits: Tuple[Union[float, None], Union[float, None]] = None,
        reaction_individual_value_limits: Tuple[Union[float, None], Union[float, None]] = None,
        reaction_value_center: Union[float, None] = None,
        reaction_individual_value_center: Union[float, None] = None,
        reaction_element_normalization: Union[str, None] = None,
        reaction_element_normalization_label: Union[str, None] = None,
        reaction_slice_by: str = None,
        reaction_slice_weights: str = None,
        reaction_slice_weight_aggregation: str = None,
        reaction_slice_share_summary: str = None,
        compound_colormap: Union[bool, str, mcolors.Colormap] = None,
        compound_colormap_limits: Tuple[float, float] = None,
        compound_individual_colormap: Union[str, mcolors.Colormap] = None,
        compound_individual_colormap_limits: Tuple[float, float] = None,
        compound_membership_colors: str = None,
        compound_reverse_overlay: bool = False,
        compound_value_limits: Tuple[Union[float, None], Union[float, None]] = None,
        compound_individual_value_limits: Tuple[Union[float, None], Union[float, None]] = None,
        compound_value_center: Union[float, None] = None,
        compound_individual_value_center: Union[float, None] = None,
        compound_element_normalization: Union[str, None] = None,
        compound_element_normalization_label: Union[str, None] = None,
        reaction_border: Union[Dict[str, object], None] = None,
        compound_border: Union[Dict[str, object], None] = None,
        reaction_bar: Union[Dict[str, object], None] = None,
        group_colormap: Union[str, mcolors.Colormap, None] = None,
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
        column colors it quantitatively (continuous colorbar), or by categorical values with a
        'reaction_value_type' or 'compound_value_type' of 'categorical' (a color per value). A
        'sample' column without a value column colors it by sample/group presence, and neither
        colors it a single presence color. The two layers are drawn together on each map, so one map
        can mix, e.g., presence of KOs with quantitative compound values. When any layer has a
        'sample' column, a 'unified' map summarizes the samples and per-sample maps are added; a
        'groups_txt' instead draws per-group maps and the 'unified' map summarizes the groups.

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
            that element's value. The default is the same as for 'reaction_gene_aggregation'. With
            slices of genes, it reduces the rows of a gene for the accessions of an element. The
            default is then 'mean', or 'circular_mean' when 'reaction_value_period' is given.

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

        reaction_value_type : Literal['numeric', 'categorical'], 'numeric'
            What the value column of the reaction file holds ('VALUE_TYPES'). Categorical values are
            names. Each takes a color of its own ('_build_categorical_txt_model'). The gene and
            accession aggregations and the sample and group summaries then take the rules of
            'CATEGORICAL_RULES'. A rule chooses one value where several meet. Without a rule, the
            sample summary shows presence on the 'unified' map. With groups and without a rule, it
            counts each group's samples on that group's map.

        compound_value_type : Literal['numeric', 'categorical'], 'numeric'
            The same for the values of the compound file.

        reaction_value_colors : str, None
            Path to a kegg-value-colors-txt file giving a color to each categorical value of the
            reaction layer ('_resolve_value_colors'). The order of its rows is the order of the
            values. The default of None takes 'DEFAULT_VALUE_COLORS' in sorted order.

        compound_value_colors : str, None
            The same colors for the categorical values of the compound layer.

        seed : int, None
            The seed of the random choices among categorical values ('RANDOM_CATEGORICAL_RULES').
            The default of None is 'DEFAULT_CATEGORICAL_SEED'. A seed is refused where no rule
            chooses at random.

        reaction_sample_summary : str, None
            How the reaction layer summarizes a set of samples: by presence ('count'/'membership')
            or by pooling their values with an aggregation. An aggregation pools the values of each
            map element. These are the values that the maps of the individual samples draw. This
            colors the 'unified' map without groups and each per-group map with groups. The default
            of None summarizes presence, by membership for 3 or fewer samples and by count above
            that. For categorical values, the default counts the samples, however few there are.

        compound_sample_summary : str, None
            The same summary for the compound layer's samples.

        reaction_group_summary : str, None
            How the reaction layer summarizes the sample groups of 'groups_txt', which colors the
            'unified' map when there are groups. The default of None summarizes presence, by
            membership for 3 or fewer groups and by count above that. For categorical values, the
            default counts the groups, however few there are. A rule of 'CATEGORICAL_RULES' chooses
            one value among the groups. It needs a rule in 'reaction_sample_summary'.

        compound_group_summary : str, None
            The same summary for the compound layer's groups.

        reaction_membership_colors : str, None
            Path to a kegg-membership-colors-txt file. It gives the color of each sample, or of each
            group when the samples are grouped. The reaction layer is then colored by membership in
            place of a colormap. Without groups, the file also colors the individual map of each
            sample. With groups, a group's maps take the group's color only with
            'GROUP_COLORMAP_TINT'. A row naming a combination of samples or groups overrides the
            blend of their colors ('_read_membership_colors_txt').

        compound_membership_colors : str, None
            The same colors for the compound layer.

        reaction_individual_colormap : Union[str, matplotlib.colors.Colormap], None
            The colormap of the single color scale shared by the reaction layer's per-sample or
            per-group maps, which 'reaction_individual_colormap_limits' trims as
            'reaction_colormap_limits' trims 'reaction_colormap'. The default of None gives both map
            contexts the layer's own 'reaction_colormap', and a colormap here separates them, so
            that a 'unified' map summarizing the samples by a statistic of another kind than the
            values themselves — their standard deviation, say — is not drawn in the colors that
            stand for magnitude on the per-sample maps.

        compound_individual_colormap : Union[str, matplotlib.colors.Colormap], None
            The same colormap for the compound layer's per-sample or per-group scale.

        reaction_value_limits : Tuple[Union[float, None], Union[float, None]], None
            The (minimum, maximum) the reaction layer's color scale may span on the 'unified' map,
            either end of which can be None to leave it wherever the values put it. A limit sets
            that end of the scale. The colorbar labels the limit. It marks the limit '<=' or '>='
            where values lie past it ('_make_quantitative_norm').

        reaction_individual_value_limits : Tuple[Union[float, None], Union[float, None]], None
            The same limits for the scale shared by the reaction layer's per-sample or per-group
            maps. It is bounded separately from the 'unified' map's scale, because a summary can
            span much less than the values it summarizes.

        compound_value_limits : Tuple[Union[float, None], Union[float, None]], None
            The same limits for the compound layer's 'unified' map scale.

        compound_individual_value_limits : Tuple[Union[float, None], Union[float, None]], None
            The same limits for the compound layer's per-sample or per-group scale.

        reaction_value_center : Union[float, None], None
            The value put at the middle of the reaction layer's color scale on the 'unified' map.
            The scale then runs the same distance either side of it, so that the middle color of a
            diverging colormap stands for this value however lopsided the values around it are
            ('_make_quantitative_norm'). The default of None leaves the scale where its values and
            any limits set it.

        reaction_individual_value_center : Union[float, None], None
            The same center for the scale shared by the reaction layer's per-sample or per-group
            maps. It is centered separately from the 'unified' map's scale, because a summary can
            span much less than the values it summarizes.

        compound_value_center : Union[float, None], None
            The same center for the compound layer's 'unified' map scale.

        compound_individual_value_center : Union[float, None], None
            The same center for the compound layer's per-sample or per-group scale.

        reaction_element_normalization : Union[str, None], None
            How to rescale each sample's or group's value for a reaction element against the same
            element's values across all samples or groups. This allows the enrichment or depletion
            of an element in a sample to be displayed rather than how much of it there is. The
            argument must be one of 'ELEMENT_NORMALIZATIONS' or the name of a pandas Series method
            that transforms each value into another value ('_resolve_element_normalization'). It
            rescales the maps of the individual samples or groups. The 'unified' map summarizes the
            unnormalized values of the samples or groups. The default of None draws each sample's
            own values.

        reaction_element_normalization_label : Union[str, None], None
            What the colorbar of the rescaled reaction scale is labeled, in place of the label the
            normalization derives from the value column's name.

        reaction_slice_by : str, None
            Slice each reaction element by its accessions ('accession') or by its genes ('gene'), as
            'SLICE_UNITS' lists. Each accession or gene of an element with a value on a map takes a
            slice in the color of its own value. Only the maps that show the values of the reaction
            layer are sliced. The reaction file needs a value column. Slices of accessions are
            summarized one accession at a time, so 'reaction_accession_aggregation' is refused for
            numeric values. Slices of genes need a 'gene_id' column with an ID in every row. One
            gene ID names the same gene in every sample. For numeric values sliced by gene,
            'reaction_accession_aggregation' reduces the rows of a gene for the accessions of an
            element. It is refused if it does not return a single value unchanged.
            'reaction_gene_aggregation' is refused, since no genes are reduced to an accession.
            'reaction_slice_weight_aggregation' is refused too, unless the rings of the border are
            sliced by accession. 'reaction_element_normalization' is not supported with slices yet.
            The slices drawn are listed in a key table, 'slice_key_reactions.txt', in the output
            directory. The default of None draws every element whole.

        reaction_slice_weights : str, None
            The column of the reaction file that weighs the slices. The shares of the slices of an
            element then follow the weights. A slice of a gene weighs the mean of the weights of its
            rows for the element's accessions. The weights also size the slices of rings that are
            sliced on their own ('reaction_border'). It needs 'reaction_slice_by' or the 'slice_by'
            of the reaction border. The default of None gives every slice on a map the same share.

        reaction_slice_weight_aggregation : str, None
            How the weights of the genes of an accession are reduced to the accession's weight, for
            a reaction file with a 'gene_id' column ('SLICE_WEIGHT_AGGREGATIONS'). It needs elements
            or rings sliced by accession. The default of None is 'sum'.

        reaction_slice_share_summary : str, None
            How the shares of slices are pooled on the maps of groups and on the 'unified' map,
            which summarize several samples ('SLICE_SHARE_SUMMARIES'). It needs a column of weights,
            and such a drawn map that shows slices of elements or of rings sliced on their own. The
            default of None is 'share_mean'.

        compound_element_normalization : Union[str, None], None
            The same normalization for the compound layer.

        compound_element_normalization_label : Union[str, None], None
            The same colorbar label for the compound layer's rescaled scale.

        reaction_border : Union[Dict[str, object], None], None
            The options of the border of the reaction layer, keyed by 'BORDER_OPTIONS', or None for
            no border. 'column' names the border column of the reaction file. Its values color a
            ring inside the outline of each box on standard maps ('kgml.SlicedKGMLCanvas'). The
            other keys are the border versions of the options of the fill. Each takes the values
            that its fill option takes, except these. A border shows values only, so its summaries
            take no presence names. Without a summary, the 'unified' map of a file with samples
            shows no border. The aggregations default to 'mean' for numbers. The colormap of numbers
            defaults to 'viridis' from 0.2 to 1. The categorical values take the colors of
            'DEFAULT_BORDER_VALUE_COLORS'. The border is built by the builder of the fill
            ('_build_txt_model', '_txt_channel'). Options left out are None. 'slice_by' divides the
            ring of each box among its accessions or genes with border values ('SLICE_UNITS'). Where
            'reaction_slice_by' names the same unit, the rings are aligned with the slices of the
            fill. Each column of the fill then takes the border value of its own unit, as bands
            above and below it. Otherwise, each ring is divided on its own, clockwise from the top
            left corner, in equal shares or in shares from 'reaction_slice_weights'
            ('_ring_result'). The slices of rings are listed in 'slice_key_reactions_border.txt'.

        compound_border : Union[Dict[str, object], None], None
            The same for the compound layer, whose rings are drawn inside compound circles on
            standard and overview maps. It takes no 'gene_aggregation' and no 'slice_by'.

        reaction_bar : Union[Dict[str, object], None], None
            The options of the bars of the reaction layer, keyed by 'BAR_OPTIONS', or None for no
            bars. 'column' names the bar column of the reaction file. Its values set the length of a
            bar in each reaction box on standard maps, drawn in the color of the fill of the box
            over a gray track ('BAR_TRACK_COLOR'). With a 'style' of 'strip', the fill is drawn
            whole, and the bar is a dark bar on a white strip along the bottom of the box
            ('BAR_STYLES'). The other keys are the bar versions of the options of the fill. A bar
            shows values only, so its summaries take no presence names. Without a summary, the
            'unified' map of a file with samples shows no bars, and its fills are drawn whole. The
            aggregations default to 'mean'. Each scale of bars includes 0, and limits that leave 0
            out are refused. On a scale with values on both sides of 0, bars grow both ways from a
            zero line, and boxes of the 'fill' style drop their labels ('label_key_reactions.txt').
            On a map that shows bars, a box with a fill and no bar value is left uncolored, or keeps
            its whole fill without a strip. The bars are built by the builder of the fill
            ('_build_txt_model', '_txt_channel'). Options left out are None. '--original-color', a
            'reaction_color' of 'none', and 'reaction_slice_by' are refused.

        Notes
        =====
        The remaining parameters carry the shared groups, per-layer colors/colormaps, and drawing
        options; see the CLI help and '_map_elements'. 'group_colormap' may be 'GROUP_COLORMAP_TINT'
        to color each group's own maps by a ramp running to that group's own color instead of by a
        named colormap, in which case 'group_colormap_limits' is how far from white that ramp runs
        ('_group_map_colors'). A 'group_colormap' of None is 'plasma_r'. 'group_colormap_scheme'
        draws the count scale of every group's maps in discrete bands ('by_count') or as a gradient
        ('by_count_continuous'), and by default keeps the bands while the colors can be told apart.
        A layer whose 'reaction_color' or 'compound_color' is 'none' uses none of the group options.
        They are refused where no other layer uses them.

        Returns
        =======
        Dict[Literal['unified', 'individual', 'grid'], Dict]
            The record returned by '_map_elements'. All three keys are always present. The
            dictionary under 'unified' is empty when 'draw_unified_maps' is False, since no
            'unified' map is drawn.
        """
        # A value type and a colors file of values are checked before any file is read.
        for element_type, value_type, value_colors in (
            ('reaction', reaction_value_type, reaction_value_colors),
            ('compound', compound_value_type, compound_value_colors)
        ):
            if value_type not in VALUE_TYPES:
                raise ConfigError(
                    f"'--{element_type}-value-type' was given as '{value_type}'. It can be "
                    f"{' or '.join(repr(name) for name in VALUE_TYPES)}."
                )
            if value_colors is not None and value_type != 'categorical':
                raise ConfigError(
                    f"'--{element_type}-value-colors' gives a color to each categorical value of "
                    f"the {element_type} layer. The values of that layer are numeric. To read them "
                    f"as categorical values, give '--{element_type}-value-type categorical'."
                )
        # The options of borders are checked as far as they can be before the files are read.
        reaction_border = dict(reaction_border or {})
        compound_border = dict(compound_border or {})
        self._check_border_options(
            reaction_txt, compound_txt, reaction_border, compound_border, reaction_color
        )
        # The options of bars are checked the same way.
        reaction_bar = dict(reaction_bar or {})
        self._check_bar_options(reaction_txt, reaction_bar, reaction_color, reaction_slice_by)
        # The options that slice reaction elements are checked as far as they can be before the file
        # is read.
        self._check_slice_options(
            reaction_txt, reaction_slice_by, reaction_slice_weights,
            reaction_slice_weight_aggregation, reaction_slice_share_summary,
            reaction_element_normalization, reaction_color, reaction_colormap,
            gene_aggregation=reaction_gene_aggregation,
            border_column=reaction_border.get('column'),
            border_slice_by=reaction_border.get('slice_by')
        )

        raw_layers: List[dict] = []
        if reaction_txt is not None:
            raw_layers.append({
                'name': 'reactions',
                'path': reaction_txt,
                'data': self._read_element_txt(
                    reaction_txt, 'reaction', reaction_value_type,
                    slice_weights=reaction_slice_weights, slice_by=reaction_slice_by,
                    border_column=reaction_border.get('column'),
                    border_value_type=reaction_border.get('value_type') or 'numeric',
                    border_slice_by=reaction_border.get('slice_by'),
                    bar_column=reaction_bar.get('column')
                ),
                'border': reaction_border,
                'bar': reaction_bar,
                'slice_by': reaction_slice_by,
                'gene_aggregation': reaction_gene_aggregation,
                'accession_aggregation': reaction_accession_aggregation,
                'color': reaction_color,
                'colormap': reaction_colormap,
                'colormap_limits': reaction_colormap_limits,
                'individual_colormap': reaction_individual_colormap,
                'individual_colormap_limits': reaction_individual_colormap_limits,
                'membership_colors_txt': reaction_membership_colors,
                'reverse_overlay': reaction_reverse_overlay,
                'sample_summary': reaction_sample_summary,
                'group_summary': reaction_group_summary,
                'value_limits': reaction_value_limits,
                'individual_value_limits': reaction_individual_value_limits,
                'value_center': reaction_value_center,
                'individual_value_center': reaction_individual_value_center,
                'element_normalization': reaction_element_normalization,
                'element_normalization_label': reaction_element_normalization_label,
                'value_period': reaction_value_period,
                'value_colors_txt': reaction_value_colors
            })
        if compound_txt is not None:
            raw_layers.append({
                'name': 'compounds',
                'path': compound_txt,
                'data': self._read_element_txt(
                    compound_txt, 'compound', compound_value_type,
                    border_column=compound_border.get('column'),
                    border_value_type=compound_border.get('value_type') or 'numeric'
                ),
                'border': compound_border,
                'bar': {},
                # A compound file has one row per accession per sample, so the gene level is a
                # reduction over a single value. Left unset, it takes the layer's default. A period
                # makes that default circular, and a period would refuse a 'sum' here.
                'gene_aggregation': None,
                'accession_aggregation': compound_accession_aggregation,
                'color': compound_color,
                'colormap': compound_colormap,
                'colormap_limits': compound_colormap_limits,
                'individual_colormap': compound_individual_colormap,
                'individual_colormap_limits': compound_individual_colormap_limits,
                'membership_colors_txt': compound_membership_colors,
                'reverse_overlay': compound_reverse_overlay,
                'sample_summary': compound_sample_summary,
                'group_summary': compound_group_summary,
                'value_limits': compound_value_limits,
                'individual_value_limits': compound_individual_value_limits,
                'value_center': compound_value_center,
                'individual_value_center': compound_individual_value_center,
                'element_normalization': compound_element_normalization,
                'element_normalization_label': compound_element_normalization_label,
                'value_period': compound_value_period,
                'value_colors_txt': compound_value_colors
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
        individuals = None
        individual_noun = None
        if all_sample_names is not None:
            if groups_txt is not None:
                sample_group, group_samples = self._relate_samples_to_groups(
                    groups_txt, all_sample_names
                )
                individuals = list(group_samples)
                individual_noun = 'group'
            else:
                individuals = all_sample_names
                individual_noun = 'sample'
            self.run.info("Samples found across the input file(s)", len(all_sample_names))

        # A seed is used only where a rule chooses at random among categorical values. A seed with
        # nothing to choose would change nothing, so it is refused. A random rule given to a layer
        # of numbers counts here, so that the refusal of that rule names the real problem
        # ('_resolve_aggregation', '_resolve_summary').
        random_rules = [
            rule for layer in raw_layers
            for rule in (
                layer['gene_aggregation'], layer['accession_aggregation'], layer['sample_summary'],
                layer['group_summary'], *(
                    layer[role].get(option) for role in ('border', 'bar') for option in (
                        'gene_aggregation', 'accession_aggregation', 'sample_summary',
                        'group_summary'
                    )
                )
            ) if rule in RANDOM_CATEGORICAL_RULES
        ]
        if seed is not None and not random_rules:
            raise ConfigError(
                f"'--seed' sets the seed of random choices among categorical values. No rule here "
                f"chooses at random. The rules that do are "
                f"{' and '.join(repr(name) for name in RANDOM_CATEGORICAL_RULES)}. Please drop "
                f"'--seed', or choose with one of those rules."
            )

        seed = DEFAULT_CATEGORICAL_SEED if seed is None else seed
        models = [
            self._build_txt_model(
                **{key: value for key, value in layer.items() if key not in ('border', 'bar')},
                all_sample_names=all_sample_names, group_samples=group_samples, seed=seed
            )
            for layer in raw_layers
        ]

        # The aggregations that reduce the values of a channel, for the run's report. A channel
        # sliced by accession shows the value of each accession, so the report leaves out its
        # accession aggregation. A channel sliced by gene leaves out its gene aggregation too. These
        # options can still choose stand-in values, or reduce the rows of a gene.
        def _used(aggregations: Union[dict, None], slice_by: Union[str, None]) -> dict:
            aggregations = dict(aggregations or {})
            if slice_by is not None:
                aggregations.pop('accession', None)
            if slice_by == 'gene':
                aggregations.pop('gene', None)
            return aggregations

        # The border of a layer is built from its border column by the builder of the fill. Its
        # channel names its options and gives its defaults. Its options left out are None. The bars
        # of a layer are built from its bar column in the same way. They are kept as a list, which
        # holds one model.
        for layer, model in zip(raw_layers, models):
            data = layer['data']
            element_type = data['element_type']
            model['border'] = None
            model['bars'] = []
            if data['border_column'] is not None:
                border = {option: layer['border'].get(option) for option in BORDER_OPTIONS}
                channel = self._txt_channel(element_type, 'border')
                model['border'] = self._build_txt_model(
                    channel['name'], self._channel_view(data, 'border'), layer['path'],
                    border['gene_aggregation'], border['accession_aggregation'], None,
                    border['colormap'], border['colormap_limits'], border['individual_colormap'],
                    border['individual_colormap_limits'], None, False, all_sample_names,
                    group_samples, border['sample_summary'], border['group_summary'],
                    value_limits=border['value_limits'],
                    individual_value_limits=border['individual_value_limits'],
                    value_center=border['value_center'],
                    individual_value_center=border['individual_value_center'],
                    element_normalization=border['element_normalization'],
                    element_normalization_label=border['element_normalization_label'],
                    value_period=border['value_period'], value_colors_txt=border['value_colors'],
                    seed=seed, slice_by=border['slice_by'], channel=channel
                )
            if data['bar_column'] is not None:
                bar = {option: layer['bar'].get(option) for option in BAR_OPTIONS}
                channel = self._txt_channel(element_type, 'bar')
                # A box without a bar value keeps its whole fill in the strip style, so the messages
                # about undefined values say that it shows no strip.
                style = bar['style'] or 'fill'
                if style == 'strip':
                    channel['consequence'] = 'shows no strip'
                    channel['consequence_plural'] = 'show no strip'
                bar_model = self._build_txt_model(
                    channel['name'], self._channel_view(data, 'bar'), layer['path'],
                    bar['gene_aggregation'], bar['accession_aggregation'], None, None, None, None,
                    None, None, False, all_sample_names, group_samples, bar['sample_summary'],
                    bar['group_summary'], value_limits=bar['value_limits'],
                    individual_value_limits=bar['individual_value_limits'], seed=seed,
                    channel=channel
                )
                bar_model['style'] = style
                model['bars'] = [bar_model]
            if model['border'] is None and not model['bars']:
                continue
            # Each column and the aggregations that reduce its values are reported. The channels
            # reduce values by different defaults, so the same column can give an element different
            # values as its fill and as its border. A layer in the fill mode 'none' colors no fills.
            # Its rows without a border value then show nothing.
            label = element_type.capitalize()
            no_column = (
                f"none, and '--{element_type}-color none' leaves the fills uncolored"
                if model['unified_mode'] == 'none' else "none, so the fills show presence"
            )
            columns = [(
                'value', data['value_column'], data['value_type'],
                _used(model.get('aggregations'), layer.get('slice_by'))
            )]
            if model['border'] is not None:
                columns.append((
                    'border', data['border_column'], data['border_value_type'],
                    _used(model['border'].get('aggregations'), border['slice_by'])
                ))
            for bar_model in model['bars']:
                columns.append(('bar', data['bar_column'], 'numeric', bar_model['aggregations']))
            for name, column, value_type, aggregations in columns:
                self.run.info(
                    f"{label} {name} column",
                    f"{column} ({value_type})" if column is not None else no_column
                )
                for level, aggregation in (aggregations or {}).items():
                    if column is None or level == 'gene' and 'gene_id' not in data['df'].columns:
                        continue
                    self.run.info(f"{label} {name} {level} aggregation", aggregation)
            if model['border'] is not None:
                self.run.info(f"{label} rows without a border value", data['border_empty'])
            if model['bars']:
                self.run.info(f"{label} bar style", model['bars'][0]['style'])
                self.run.info(f"{label} rows without a bar value", data['bar_empty'])
            if model['border'] is not None and border['slice_by'] is not None:
                self.run.info(
                    f"{label} border slices",
                    f"by {border['slice_by']}, " + (
                        "aligned with the fill" if border['slice_by'] == layer.get('slice_by')
                        else "independent of the fill"
                    )
                )

        # The reaction layer comes first. Slicing it needs a drawn map that shows its values.
        draw_individual_maps = all_sample_names is not None and (
            draw_individual_files is not False or draw_grid is not False
        )
        if reaction_slice_by is not None:
            models[0].update(self._build_slice_model(
                models[0], raw_layers[0]['data'], reaction_slice_by,
                reaction_slice_weight_aggregation, reaction_slice_share_summary,
                reaction_accession_aggregation, all_sample_names, group_samples, draw_unified_maps,
                draw_individual_maps
            ))
        # The rings of the reaction border are sliced from the rows with a border value. Where the
        # fill is sliced by the same unit, the rings are aligned with its columns.
        border_slice_by = reaction_border.get('slice_by')
        if border_slice_by is not None:
            fill = models[0]
            border_slice_model = self._build_slice_model(
                fill['border'], self._channel_view(raw_layers[0]['data'], 'border'),
                border_slice_by, reaction_slice_weight_aggregation, reaction_slice_share_summary,
                reaction_border.get('accession_aggregation'), all_sample_names, group_samples,
                draw_unified_maps, draw_individual_maps,
                fill=fill if border_slice_by == reaction_slice_by else None
            )
            fill['border'].update(border_slice_model)
        # Shares are pooled on the maps that summarize several samples. A summary of shares needs
        # such a map that shows the slices of the fill, or rings that are sliced on their own.
        if reaction_slice_share_summary is not None:
            fill = models[0]
            ring = fill['border'] if border_slice_by is not None else None
            if not (fill.get('slice_pooled') or (
                ring is not None and ring['slice_layout'] == 'independent' and ring['slice_pooled']
            )):
                shown = "the values of the reaction layer" if ring is None else (
                    "the values of the reaction layer or the slices of its rings"
                )
                raise ConfigError(
                    f"'--reaction-slice-share-summary' pools the shares of slices on the maps that "
                    f"summarize several samples, which are the maps of groups and the 'unified' "
                    f"map. Here no drawn map of this kind shows {shown}, so no shares are pooled. "
                    f"Please drop '--reaction-slice-share-summary'."
                )

        # A layer in the fill mode 'none' colors no fills, so it neither summarizes the groups by
        # presence nor counts the samples of a group on the group's map. The group options that do
        # these things are refused where no other layer uses them. The CLI makes this check too.
        unfilled = [model for model in models if model['unified_mode'] == 'none']
        if unfilled and group_samples is not None:
            others = [model for model in models if model['unified_mode'] != 'none']
            threshold_used = any(model['unified_mode'] == 'membership' for model in others)
            counts_used = any(model['individual_mode'] == 'membership' for model in others)
            idle = [f"'{flag}'" for flag, given, used in (
                ('--group-threshold', group_threshold is not None, threshold_used),
                (
                    '--group-colormap',
                    group_colormap is not None or group_colormap_limits is not None, counts_used
                ),
                ('--group-colormap-scheme', group_colormap_scheme is not None, counts_used),
                ('--group-reverse-overlay', group_reverse_overlay, counts_used)
            ) if given and not used]
            if idle:
                unfilled_clauses = ' '.join(
                    f"'--{model['element_type']}-color none' leaves the fills of the "
                    f"{model['element_type']} layer uncolored." for model in unfilled
                )
                raise ConfigError(
                    f"These options shape how the fills of a layer show groups of samples: "
                    f"{', '.join(idle)}. No layer here uses them. {unfilled_clauses} Such a layer "
                    f"draws only the rings of its border column. Please drop these options."
                )

        # Groups are carried through whenever they are defined: they set which samples each group's
        # map summarizes, and, for a group summary of presence, which groups the 'unified' map shows
        # an element in (per 'group_threshold').
        grouped_membership = None
        if group_samples is not None:
            grouped_membership = {
                'source_group': sample_group,
                'group_sources': group_samples,
                'group_threshold': group_threshold,
                'group_colormap': 'plasma_r' if group_colormap is None else group_colormap,
                'group_colormap_limits': group_colormap_limits,
                'group_reverse_overlay': group_reverse_overlay,
                'group_colormap_scheme': group_colormap_scheme
            }

        individual_subject = f"{individual_noun}s" if individual_noun is not None else None
        return self._map_elements(
            models,
            output_dir,
            pathway_numbers=pathway_numbers,
            individuals=individuals,
            individual_noun=individual_noun,
            grid_source_type='sample',
            colorbar_individual_suffix=individual_subject,
            subset_subject=individual_subject,
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
        reaction_membership_colors: str = None,
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
        'colormap_scheme', 'reaction_membership_colors', 'reaction_reverse_overlay',
        'draw_unified_maps', 'draw_individual_files', 'draw_grid', and the group colormap options)
        mirror the presence engine; see the CLI help and '_map_element_membership'.
        'reaction_membership_colors' is the path to a kegg-membership-colors-txt file. It gives the
        color of each contigs database, or of each group when the contigs databases are grouped. The
        layer is then colored by membership in place of a colormap. Without groups, the file also
        colors the individual map of each contigs database. With groups, a group's own maps take the
        group's color only where 'group_colormap' is 'GROUP_COLORMAP_TINT'. That value colors each
        group's own maps by a ramp running to the group's own color instead of by a named colormap.
        'group_colormap_limits' is then how far from white the ramp runs ('_group_map_colors').
        'group_colormap_scheme' draws the count scale of every group's maps in discrete bands
        ('by_count') or as a gradient ('by_count_continuous'), and by default keeps the bands while
        the colors can be told apart.

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

        # The membership colors file is read before the databases are opened, so a malformed file is
        # reported before the run spends time loading data it may not draw. Its names are checked
        # against the databases or their groups in '_map_elements'.
        membership_colors, combo_colors = self._read_membership_colors_txt(
            reaction_membership_colors, '--reaction-membership-colors'
        ) if reaction_membership_colors is not None else (None, {})

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
            'membership_colors_flag': '--reaction-membership-colors',
            'membership_colors': membership_colors,
            'membership_combo_colors': combo_colors,
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
        reaction_membership_colors: str = None,
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
        'colormap_scheme', 'reaction_membership_colors', 'reaction_reverse_overlay',
        'draw_unified_maps', 'draw_individual_files', 'draw_grid', and the group colormap options)
        mirror the presence engine; see the CLI help and '_map_element_membership'.
        'reaction_membership_colors' is the path to a kegg-membership-colors-txt file. It gives the
        color of each genome, or of each group when the genomes are grouped. The layer is then
        colored by membership in place of a colormap. Without groups, the file also colors the
        individual map of each genome. With groups, a group's own maps take the group's color only
        where 'group_colormap' is 'GROUP_COLORMAP_TINT'. That value colors each group's own maps by
        a ramp running to the group's own color instead of by a named colormap.
        'group_colormap_limits' is then how far from white the ramp runs ('_group_map_colors').
        'group_colormap_scheme' draws the count scale of every group's maps in discrete bands
        ('by_count') or as a gradient ('by_count_continuous'), and by default keeps the bands while
        the colors can be told apart.

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

        # The membership colors file is read before the gene clusters are loaded, so a malformed
        # file is reported before the run spends time on data it may not draw. Its names are checked
        # against the genomes or their groups in '_map_elements'.
        membership_colors, combo_colors = self._read_membership_colors_txt(
            reaction_membership_colors, '--reaction-membership-colors'
        ) if reaction_membership_colors is not None else (None, {})

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
            'membership_colors_flag': '--reaction-membership-colors',
            'membership_colors': membership_colors,
            'membership_combo_colors': combo_colors,
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
                {**layer, 'unified_mode': 'single', 'individual_mode': 'single'}
                for layer in presence_layers
            ]
            self._check_reserved_colors(check_layers, pathway_numbers, drawn_individuals=[])
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

    def _read_membership_colors_txt(
        self,
        path: str,
        flag: str
    ) -> Tuple[Dict[str, str], Dict[Tuple[str, ...], str]]:
        """
        Load the colors of a kegg-membership-colors-txt file.

        The file names a color for each sample, contigs database, genome, or group. Coloring by
        membership uses that color for the elements in that sample, database, genome, or group
        alone. The individual map of a sample, database, or genome is drawn in its color. A group's
        ramp of tints ('GROUP_COLORMAP_TINT') is made from the group's color. A row whose first
        field lists several names separated by 'MEMBERSHIP_COMBO_SEPARATOR' instead gives the color
        of that combination. That color replaces the blend of the members' colors that coloring by
        membership would otherwise derive.

        The names are not checked against the run's samples, contigs databases, genomes, or groups
        here. The file is read before they are all known. '_resolve_membership_colors' checks them.

        Parameters
        ==========
        path : str
            Path to a tab-delimited file. Its first column holds the names of samples, contigs
            databases, genomes, or groups, or combinations of them. Its 'COLORS_COLUMN' column holds
            color hex codes.

        flag : str
            The command-line flag this file came from, used in error messages.

        Returns
        =======
        Tuple[Dict[str, str], Dict[Tuple[str, ...], str]]
            membership_colors : maps each name of a sample, contigs database, genome, or group
                to its color hex code.
            combo_colors : maps each combination of those names, sorted and as a tuple, to the
                color hex code overriding the blend of that combination.
        """
        filesnpaths.is_file_tab_delimited(path)
        columns = utils.get_columns_of_TAB_delim_file(path, include_first_column=True)
        if COLORS_COLUMN not in columns:
            self.progress.end()
            raise ConfigError(
                f"The colors file given to '{flag}', at '{path}', should have a column called "
                f"'{COLORS_COLUMN}' holding a color hex code for each sample, contigs database, "
                f"genome, or group. Its columns are: "
                f"{', '.join(repr(column) for column in columns)}."
            )
        if len(columns) < 2:
            self.progress.end()
            raise ConfigError(
                f"The colors file given to '{flag}', at '{path}', should have at least two "
                f"columns. The first column should hold the names of the samples, contigs "
                f"databases, genomes, or groups being colored. A '{COLORS_COLUMN}' column should "
                f"hold a color hex code for each of them."
            )
        if columns[0] == COLORS_COLUMN:
            self.progress.end()
            raise ConfigError(
                f"The first column of the colors file given to '{flag}', at '{path}', is the "
                f"'{COLORS_COLUMN}' column. Anvi'o expects the first column to hold the names of "
                f"samples, contigs databases, genomes, or groups, so the two columns need to be "
                f"the other way around."
            )

        table = utils.get_TAB_delimited_file_as_dictionary(path)

        membership_colors: Dict[str, str] = {}
        combo_colors: Dict[Tuple[str, ...], str] = {}
        # What each normalized key was written as in the file, for reporting two rows that mean one
        # thing.
        seen_rows: Dict[Tuple[str, ...], str] = {}
        for name, row in table.items():
            color = str(row[COLORS_COLUMN]).strip()
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
            members = tuple(sorted(
                {part.strip() for part in name.split(MEMBERSHIP_COMBO_SEPARATOR) if part.strip()}
            ))
            if not members:
                self.progress.end()
                raise ConfigError(
                    f"A row of the colors file given to '{flag}', at '{path}', names nothing at "
                    f"all in its first column. Every row should name a sample, contigs database, "
                    f"genome, or group. A row can also name a combination of them, separated by "
                    f"'{MEMBERSHIP_COMBO_SEPARATOR}'."
                )
            # Two rows can name the same thing without being the same text — stray spaces around a
            # name, or the members of a combination written in another order — and the duplicate
            # check that reading the file has already done compares the text. Left alone, the later
            # row would quietly replace the earlier one and the file would report a count that does
            # not match its rows, so what the text meant is checked here too.
            if members in seen_rows:
                separator = f'{MEMBERSHIP_COMBO_SEPARATOR} '
                self.progress.end()
                raise ConfigError(
                    f"The colors file given to '{flag}', at '{path}', gives two colors to "
                    f"{'the combination ' if len(members) > 1 else ''}'{separator.join(members)}'. "
                    f"The rows '{seen_rows[members]}' and '{name}' both name it. Please give it "
                    f"one row."
                )
            seen_rows[members] = name
            if len(members) == 1:
                membership_colors[members[0]] = kgml.canonical_color(color)
            else:
                combo_colors[members] = kgml.canonical_color(color)

        if not membership_colors and not combo_colors:
            self.progress.end()
            raise ConfigError(
                f"The colors file given to '{flag}', at '{path}', has its header but no rows, so "
                f"it gives no color to anything. Each row should name a sample, contigs database, "
                f"genome, or group in its first column. Each row should also give a color hex code "
                f"in the '{COLORS_COLUMN}' column."
            )
        if not membership_colors:
            self.progress.end()
            raise ConfigError(
                f"The colors file given to '{flag}', at '{path}', gives colors only for "
                f"combinations. It gives no color for a single sample, contigs database, genome, "
                f"or group. A combination's color adjusts the blend of its members' colors, so the "
                f"members need colors first."
            )

        self.run.info(f"Names colored by '{flag}'", len(membership_colors))
        if combo_colors:
            self.run.info(f"Combinations recolored by '{flag}'", len(combo_colors))
        return membership_colors, combo_colors

    def _resolve_membership_colors(
        self,
        layers: List[dict],
        individuals: List[str],
        individual_noun: str,
        draw_unified_maps: bool = True
    ) -> None:
        """
        Reduce each layer's membership colors to those the run will actually use.

        A sample, contigs database, genome, or group with no color of its own cannot be colored at
        all, so it is an error. The error lists all of them at once. A color may be given for a name
        that is not a sample, contigs database, genome, or group of the run, or for a combination
        that names one. Such a color describes something this run does not draw, so it is reported
        and ignored, as a groups file's extra items are ('_relate_samples_to_groups'). One file can
        then cover a set of samples that different runs draw different subsets of. Ignoring such a
        color means dropping it from the layer rather than merely reporting it, which is what makes
        the colors that remain the colors of the run — everything downstream reads them as such. A
        color naming a combination is dropped the same way where the 'unified' map is not drawn.
        Only that map's scale gives a combination a color of its own.

        Parameters
        ==========
        layers : List[dict]
            The layer models, whose 'membership_colors'/'membership_combo_colors' are checked
            against the run and reduced to it in place. A layer without colors of its own is
            skipped.

        individuals : List[str]
            The names of every sample, contigs database, genome, or group of the run.

        individual_noun : str
            What one of 'individuals' is called in error messages, e.g. 'sample' or 'group'.

        draw_unified_maps : bool, True
            If True, the 'unified' map is drawn. If False, it is not. Only that map's color scale
            blends the colors of a combination. When it is not drawn, a color given for a
            combination is reported and dropped.
        """
        individual_set = set(individuals)
        # A combination is written as its names separated by 'MEMBERSHIP_COMBO_SEPARATOR', so a name
        # containing that separator could not be told from a combination of other names.
        containing_separator = sorted(
            name for name in individual_set if MEMBERSHIP_COMBO_SEPARATOR in name
        )
        if containing_separator:
            message = ', '.join(f"'{name}'" for name in containing_separator)
            self.progress.end()
            raise ConfigError(
                f"A colors file names a combination of {individual_noun}s by separating their "
                f"names with '{MEMBERSHIP_COMBO_SEPARATOR}', so no {individual_noun} name may "
                f"itself contain that character. "
                f"{'These names do' if len(containing_separator) > 1 else 'This name does'}: "
                f"{message}. Please rename {'them' if len(containing_separator) > 1 else 'it'} in "
                f"the input, or drop the colors file and color these {individual_noun}s from a "
                f"colormap instead."
            )

        for layer in layers:
            membership_colors = layer.get('membership_colors')
            if membership_colors is None:
                continue
            flag = layer['membership_colors_flag']

            missing = [
                individual for individual in individuals if individual not in membership_colors
            ]
            if missing:
                message = ', '.join(f"'{individual}'" for individual in missing)
                self.progress.end()
                raise ConfigError(
                    f"The colors file given to '{flag}' gives no color for "
                    f"{'these' if len(missing) > 1 else 'this'} {individual_noun}"
                    f"{'s' if len(missing) > 1 else ''} of the run: {message}. Coloring by "
                    f"membership needs a color for every {individual_noun}, since an element is "
                    f"colored by exactly which of them contain it, so anvi'o will not fill in the "
                    f"gaps from a colormap: that would put colors chosen here and colors chosen "
                    f"for you on one scale. Please add {'them' if len(missing) > 1 else 'it'} to "
                    f"the file."
                )

            extra = sorted(name for name in membership_colors if name not in individual_set)
            extra_combos = [
                combo for combo in layer['membership_combo_colors']
                if not set(combo) <= individual_set
            ]
            if extra:
                message = ', '.join(f"'{name}'" for name in extra)
                self.run.warning(
                    f"The colors file given to '{flag}' colors the following names, which are not "
                    f"{individual_noun}s of this run and so will not factor into the maps: "
                    f"{message}"
                )
            if extra_combos:
                separator = f'{MEMBERSHIP_COMBO_SEPARATOR} '
                message = ', '.join(
                    f"'{name}'" for name in sorted(separator.join(combo) for combo in extra_combos)
                )
                self.run.warning(
                    f"The colors file given to '{flag}' recolors the following combinations, which "
                    f"name at least one thing that is not a {individual_noun} of this run, and so "
                    f"will not factor into the maps: {message}"
                )

            # Drop them, rather than only report them. Everything downstream reads what is left as
            # the layer's colors: the check that no color is one a map reserves, and the check that
            # the layers agree on each group's ramp. A color kept here for something the run does
            # not draw could fail either check, refusing a run it has no part in.
            for name in extra:
                del membership_colors[name]
            for combo in extra_combos:
                del layer['membership_combo_colors'][combo]

            # Only the 'unified' map's scale blends the colors of a combination, so a color that
            # overrides one of those blends has no effect where that map is not drawn. These rows
            # are dropped with a warning. The run is not refused. The colors of single samples,
            # contigs databases, genomes, or groups still draw the individual maps, so one file can
            # serve a run that draws the 'unified' map and a run that does not.
            combo_colors = layer['membership_combo_colors']
            if not draw_unified_maps and combo_colors:
                separator = f'{MEMBERSHIP_COMBO_SEPARATOR} '
                message = ', '.join(
                    f"'{name}'" for name in
                    sorted(separator.join(combo) for combo in combo_colors)
                )
                self.run.warning(
                    f"The colors file given to '{flag}' recolors the following combinations of "
                    f"{individual_noun}s: {message}. A combination takes a color only on the "
                    f"'unified' map, where an element is colored by exactly which "
                    f"{individual_noun}s contain it. That map is not drawn here, so these rows are "
                    f"ignored. The colors given to single {individual_noun}s are unaffected. Each "
                    f"one still draws that {individual_noun}'s own map."
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
        subject: str = 'these maps',
        baseline: Union[float, None] = None
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

        baseline : Union[float, None], None
            A value that the scale always includes, or None. A bar grows from 0, so the scale of
            bars includes 0 wherever the values would put its ends ('TXT_CHANNEL_ROLES'). Its limits
            include it too ('_build_txt_model'). Such a scale takes no center.

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
            centered, they give it the middle color instead. With a baseline, the scale is a single
            value only where every value is the baseline, and callers draw bars of no length
            ('_bar_extent').
        """
        if not values:
            return None, None, None, False, False
        data_min = min(values)
        data_max = max(values)
        limit_min, limit_max = (None, None) if limits is None else limits
        vmin = data_min if limit_min is None else limit_min
        vmax = data_max if limit_max is None else limit_max
        if baseline is not None:
            vmin = min(vmin, baseline) if limit_min is None else vmin
            vmax = max(vmax, baseline) if limit_max is None else vmax
            # A value a rounding error past the baseline, such as the mean of values that cancel
            # out, does not carry the scale past it. Its bar has no length, and the scale is not
            # taken to have values on both sides of the baseline ('_stage_entry_bar').
            tolerance = 1e-9 * (vmax - vmin)
            if limit_min is None and 0 < baseline - vmin <= tolerance:
                vmin = baseline
            if limit_max is None and 0 < vmax - baseline <= tolerance:
                vmax = baseline
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
        ['R00764', 'R00756']). The source can name multiple accessions or be absent. An accession
        that it names twice, as entry 179 of map 00670 names R08557, is returned once, so it counts
        once in the value of the element.

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
            The accessions (the part of each token after the colon), each once and in the order they
            first appear, or an empty list if the source attribute is absent.
        """
        source = entry.reaction if use_reaction_attribute else entry.name
        if not source:
            return []
        return list(dict.fromkeys(token.split(':')[1] for token in source.split()))

    @staticmethod
    def _reduce_element_value(
        accessions: Iterable[str],
        values: Dict[str, float],
        aggregate: Callable
    ) -> Union[float, None]:
        """
        Reduce the values of a map element's KEGG accessions to a single element value.

        Parameters
        ==========
        accessions : Iterable[str]
            The accessions the element stands for, each once ('_get_entry_kegg_ids').

        values : Dict[str, float]
            Keys are KEGG accessions (KO, reaction, or compound IDs), values are per-accession
            values.

        aggregate : Callable
            Reduces a list of values to a single value (see 'AGGREGATION_FUNCTIONS').

        Returns
        =======
        Union[float, None]
            The aggregated value, or None if none of the accessions have a value or the aggregation
            is undefined for the values they do have.
        """
        element_values = [values[accession] for accession in accessions if accession in values]
        if not element_values:
            return None
        # An aggregation can be undefined for the values of this element even when it is defined
        # elsewhere on the map, as the standard deviation is for an element with a single accession,
        # in which case the element has no value and is left uncolored.
        value = aggregate(element_values)
        return value if np.isfinite(value) else None

    @staticmethod
    def _summarize_entry_values(
        entry: kgml.Entry,
        layer: dict
    ) -> Tuple[Dict[str, float], Union[float, None]]:
        """
        Find an Entry's value in each sample or group of a layer, and on the 'unified' map.

        The Entry's accessions name its map element ('_get_entry_kegg_ids').
        '_summarize_element_values' finds the element's values with the layer's accession
        aggregation and keeps them in the layer's '_element_cache'.

        Parameters
        ==========
        entry : kgml.Entry
            An Entry (ortholog or compound).

        layer : dict
            A 'quantitative' layer model from '_build_txt_model'. The range-finding pass of
            '_map_elements' gives it an '_element_cache' ('_summarize_element_values').

        Returns
        =======
        Tuple[Dict[str, float], Union[float, None]]
            The element's value in each sample (ungrouped) or group (grouped) that has one, and its
            value on the 'unified' map, or None where it has none. The first is empty for a layer
            without samples.
        """
        accessions = tuple(Mapper._get_entry_kegg_ids(entry, layer['use_reaction_attribute']))
        return Mapper._summarize_element_values(
            accessions, layer, layer['aggregate'], layer['_element_cache']
        )

    @staticmethod
    def _summarize_element_values(
        accessions: Tuple[str, ...],
        layer: dict,
        aggregate: Callable,
        cache: Dict[Tuple[str, ...], Tuple[Dict[str, float], Union[float, None]]]
    ) -> Tuple[Dict[str, float], Union[float, None]]:
        """
        Find a map element's value in each sample or group of a layer, and on the 'unified' map.

        An element's value in a sample is the aggregate of its accessions in that sample
        ('_reduce_element_value'). This is the value the sample's own map draws. A summary then
        pools these element values, so the 'unified' map summarizes what the maps of the individual
        samples or groups draw. Grouped, each group's value is the sample summary of its samples'
        element values. The 'unified' value is then the group summary of the groups' values.

        A summary is over the samples or groups in which the element has a value. A sample with no
        row for any of the element's accessions does not count as a zero. A summary can be
        undefined, as the standard deviation of a single value is. The element then has no value in
        that group, or on the 'unified' map. The tuple of its accessions is then added to the
        layer's '_undefined_individual' or '_undefined_unified' set. '_warn_undefined_summaries'
        reports these sets once for the layer.

        Parameters
        ==========
        accessions : Tuple[str, ...]
            The accessions the element stands for, each once ('_get_entry_kegg_ids'). They name the
            element in the cache and in the layer's undefined sets.

        layer : dict
            A 'quantitative' layer model from '_build_txt_model'. A layer without samples has
            'sample_values' of None. Its one value per element comes from 'unified_values'. The
            range-finding pass of '_map_elements' gives the layer the two undefined sets.

        aggregate : Callable
            Reduces the values of the element's accessions to one value ('_reduce_element_value').
            Each sample gets its own value, and a layer without samples gets one from
            'unified_values'. '_summarize_entry_values' passes the layer's accession aggregation.

        cache : Dict[Tuple[str, ...], Tuple[Dict[str, float], Union[float, None]]]
            Results already worked out for this layer and aggregation, keyed by an element's
            accessions. With the layer and aggregation fixed, a result depends only on the
            accessions. Entries on different maps, with different IDs, can stand for the same
            accessions and share one result. '_summarize_entry_values' passes the layer's
            '_element_cache'.

        Returns
        =======
        Tuple[Dict[str, float], Union[float, None]]
            The element's value in each sample (ungrouped) or group (grouped) that has one, and its
            value on the 'unified' map, or None where it has none. The first is empty for a layer
            without samples.
        """
        # One call works out the element's value in every sample or group and on the 'unified' map.
        # The range-finding pass makes the first call. The result is then reused by the 'unified'
        # map and the map of each sample or group. The result is keyed by the element's accessions,
        # not by the ID of an Entry, which is unique to one map. An element's value depends only on
        # its accessions, so the result is also reused for an Entry on another KEGG map, or a second
        # Entry on one KEGG map, that stands for the same accessions.
        if accessions in cache:
            return cache[accessions]

        # A layer without samples has one value per accession. The element's value is the aggregate
        # of its accessions, and every map shows it.
        if layer['sample_values'] is None:
            result = {}, Mapper._reduce_element_value(
                accessions, layer['unified_values'], aggregate
            )
            cache[accessions] = result
            return result

        # The element's value in each sample is the aggregate of the accessions the sample has
        # values for. A sample with none of them has no value. It is left out of every summary
        # below.
        sample_element_values: Dict[str, float] = {}
        for sample, values in layer['sample_values'].items():
            value = Mapper._reduce_element_value(accessions, values, aggregate)
            if value is not None:
                sample_element_values[sample] = value

        # A summary can be undefined for the values it is given, as a standard deviation is for a
        # single value. The element then has no value at that level. Its accessions are recorded for
        # '_warn_undefined_summaries'.
        def _summarize(values: List[float], summary: Callable, undefined_key: str):
            value = summary(values)
            if np.isfinite(value):
                return value
            layer[undefined_key].add(accessions)
            return None

        if layer['group_samples'] is None:
            # Ungrouped, the map of each sample shows the sample's own value. The sample summary
            # pools these values on the 'unified' map.
            individual_values = sample_element_values
            unified_aggregate = layer['sample_aggregate']
        else:
            # Grouped, the map of each group shows the sample summary of its own samples' values. A
            # group none of whose samples has a value has no value either. A sample summary of
            # presence colors the group maps by count instead ('_group_map_colors'), so the groups
            # get no values here. The group summary pools the groups' values on the 'unified' map.
            individual_values = {}
            if layer['sample_aggregate'] is not None:
                for group, samples in layer['group_samples'].items():
                    group_values = [
                        sample_element_values[sample] for sample in samples
                        if sample in sample_element_values
                    ]
                    if not group_values:
                        continue
                    value = _summarize(
                        group_values, layer['sample_aggregate'], '_undefined_individual'
                    )
                    if value is not None:
                        individual_values[group] = value
            unified_aggregate = layer['group_aggregate']

        # The 'unified' value pools the samples or groups that have a value. The summary is None
        # where the 'unified' map shows presence, and then there is no value to work out.
        unified_value = None
        if unified_aggregate is not None and individual_values:
            unified_value = _summarize(
                list(individual_values.values()), unified_aggregate, '_undefined_unified'
            )
        result = individual_values, unified_value
        cache[accessions] = result
        return result

    @staticmethod
    def _summarize_entry_individuals(
        entry: kgml.Entry,
        layer: dict,
        individual: Union[str, None] = None,
        count_choices: bool = True
    ) -> Tuple[Dict[str, str], Union[str, None]]:
        """
        Find an Entry's categorical value in each sample or group, and on the 'unified' map.

        The Entry's accessions name its map element ('_get_entry_kegg_ids').
        '_summarize_element_individuals' finds the element's values.

        Parameters
        ==========
        entry : kgml.Entry
            An Entry (ortholog or compound).

        layer : dict
            A layer model of categorical values ('_build_categorical_txt_model').

        individual : Union[str, None], None
            The sample or group whose map asks for the values, or None for the 'unified' map
            ('_summarize_element_individuals').

        count_choices : bool, True
            False to find the values without counting the choices made for them
            ('_summarize_element_individuals').

        Returns
        =======
        Tuple[Dict[str, str], Union[str, None]]
            The element's value in each sample (ungrouped) or group (grouped) that has one, and its
            value on the 'unified' map, or None where it has none. The first is empty for a layer
            without samples, and for a grouped layer without a rule for its samples.
        """
        accessions = tuple(Mapper._get_entry_kegg_ids(entry, layer['use_reaction_attribute']))
        return Mapper._summarize_element_individuals(accessions, layer, individual, count_choices)

    @staticmethod
    def _summarize_element_individuals(
        accessions: Tuple[str, ...],
        layer: dict,
        individual: Union[str, None] = None,
        count_choices: bool = True
    ) -> Tuple[Dict[str, str], Union[str, None]]:
        """
        Find a map element's categorical value in each sample or group, and on the 'unified' map.

        This is the categorical counterpart of '_summarize_element_values'. Each accession of the
        element that has a value gives one vote. The accession rule chooses among the votes in each
        sample ('_choose_categorical_value'). Where the layer has no samples, it chooses among the
        votes for the one map. Where the sample summary is a rule, each sample with a value for the
        element gives one vote, and that rule chooses among them. Ungrouped, it chooses for the
        'unified' map. Grouped, it chooses for each group among the group's samples. The group rule
        then chooses for the 'unified' map among the groups. A sample or group without a value does
        not vote. Without a rule, the 'unified' map shows presence, and no value is chosen for it.

        Results are kept in the layer's '_element_cache', keyed by the element's sorted accessions.
        The values of an element are worked out for every sample at once. The layer's '_choices'
        counts how each choice among different values was made ('_report_categorical_choices'). A
        choice is counted the first time a map that uses it asks for the element. It is therefore
        counted once, and only if a drawn map uses it. The map of a sample uses the choice made in
        that sample. The map of a group uses the choices made in its samples and the choice among
        them. The 'unified' map uses every choice that leads to its value. A layer without samples
        has one choice per element, which every map uses.

        Parameters
        ==========
        accessions : Tuple[str, ...]
            The accessions the element stands for, each once ('_get_entry_kegg_ids').

        layer : dict
            A layer model of categorical values ('_build_categorical_txt_model').

        individual : Union[str, None], None
            The sample or group whose map asks for the values, or None for the 'unified' map. It
            sets which choices are counted. The values returned do not depend on it.

        count_choices : bool, True
            False to count no choice. The range of a border asks for values this way, since no map
            is drawn then ('_map_elements').

        Returns
        =======
        Tuple[Dict[str, str], Union[str, None]]
            The element's value in each sample (ungrouped) or group (grouped) that has one, and its
            value on the 'unified' map, or None where it has none. The first is empty for a layer
            without samples, and for a grouped layer without a rule for its samples.
        """
        # The element is named by its sorted accessions, in the cache and in the key of a random
        # choice. Two maps can list the same accessions in different orders. The element then gets
        # the same value on both maps, and its choices are counted once.
        element_key = tuple(sorted(accessions))
        cache = layer['_element_cache']
        choices = layer['_choices']

        def _count_used(result):
            # Count the choices that the asking map uses. A counted choice is dropped, so no choice
            # is counted twice.
            individual_values, unified_value, uncounted = result
            if not count_choices:
                return individual_values, unified_value
            if individual is None or layer['sample_values'] is None:
                used = list(uncounted)
            elif layer['group_samples'] is None:
                used = [('element', individual)]
            else:
                used = [('element', sample) for sample in layer['group_samples'][individual]]
                used.append(('sample', individual))
            for level, voter in used:
                how = uncounted.pop((level, voter), None)
                if how is not None:
                    choices[level][how] += 1
            return individual_values, unified_value

        if element_key in cache:
            return _count_used(cache[element_key])
        order = layer['value_order']
        seed = layer['seed']
        # The key of a random choice starts with the salt of the channel ('_txt_channel').
        salt = layer['choice_salt']
        # How each choice among different values was made, keyed by its level and by the sample or
        # group it was made for. The key holds '' instead for a choice on the 'unified' map, and for
        # a choice in a layer without samples.
        uncounted: Dict[Tuple[str, str], str] = {}

        def _choose(values: Dict[str, str], sample: str) -> Union[str, None]:
            votes = [values[accession] for accession in accessions if accession in values]
            if not votes:
                return None
            value, how = Mapper._choose_categorical_value(
                votes, layer['accession_rule'], order, seed,
                salt + ('accession', sample) + element_key
            )
            if how is not None:
                uncounted[('element', sample)] = how
            return value

        if layer['sample_values'] is None:
            cache[element_key] = {}, _choose(layer['unified_values'], ''), uncounted
            return _count_used(cache[element_key])

        sample_element_values: Dict[str, str] = {}
        for sample, values in layer['sample_values'].items():
            value = _choose(values, sample)
            if value is not None:
                sample_element_values[sample] = value

        def _choose_across(
            votes: List[str], rule: str, level: str, voter: str, key: Tuple[str, ...]
        ) -> Union[str, None]:
            value, how = Mapper._choose_categorical_value(votes, rule, order, seed, salt + key)
            if how is not None:
                uncounted[(level, voter)] = how
            return value

        group_samples = layer['group_samples']
        if group_samples is None:
            # Ungrouped, each sample's map shows the sample's own value. The sample rule chooses
            # among them for the 'unified' map.
            individual_values = sample_element_values
            unified_rule, unified_level = layer['sample_rule'], 'sample'
        else:
            # Grouped, the sample rule chooses each group's value among the group's samples. Each
            # sample with a value gives one vote. A group with no such sample has no value, and
            # neither does a group whose rule gave none. The group rule then chooses among the
            # groups for the 'unified' map. Each group with a value gives one vote, whatever its
            # size.
            individual_values = {}
            if layer['sample_rule'] is not None:
                for group, samples in group_samples.items():
                    votes = [
                        sample_element_values[sample] for sample in samples
                        if sample in sample_element_values
                    ]
                    if not votes:
                        continue
                    value = _choose_across(
                        votes, layer['sample_rule'], 'sample', group,
                        ('sample', group) + element_key
                    )
                    if value is not None:
                        individual_values[group] = value
            unified_rule, unified_level = layer['group_rule'], 'group'
        unified_value = None
        if unified_rule is not None and individual_values:
            unified_value = _choose_across(
                list(individual_values.values()), unified_rule, unified_level, '',
                (unified_level,) + element_key
            )
        cache[element_key] = individual_values, unified_value, uncounted
        return _count_used(cache[element_key])

    @staticmethod
    def _normalize_entry_value(
        entry: kgml.Entry,
        layer: dict,
        individuals: Iterable[str],
        normalize,
        cache: Union[Dict[Tuple[str, ...], Dict[str, float]], None] = None
    ) -> Dict[str, float]:
        """
        Rescale an Entry's value in each sample or group against its values across all of them.

        An element's value in each sample or group comes from '_summarize_entry_values', and
        '_normalize_element_values' rescales these values. Each value is built from the aggregate of
        the element's accessions. Which accessions an element stands for is known only from the KEGG
        map itself, so a ratio against the element's mean, for instance, is a ratio of aggregates.

        The normalization sees the values of only the samples or groups in which the element has a
        value at all. A sample with no row for any of the element's accessions is not counted as a
        zero. A sample or group summary follows the same rule.

        Parameters
        ==========
        entry : kgml.Entry
            An Entry (ortholog or compound) whose accessions are read by '_get_entry_kegg_ids'.

        layer : dict
            A 'quantitative' layer model with samples. Its element values are the ones rescaled.

        individuals : Iterable[str]
            The samples or groups to normalize across, in the order they are colored in.

        normalize : callable
            Rescales an element's values across the samples or groups that have one
            ('_resolve_element_normalization').

        cache : Union[Dict[Tuple[str, ...], Dict[str, float]], None], None
            Results already worked out for this layer, keyed by an Entry's accessions, which is all
            an answer depends on once the layer's values, aggregation, and normalization are fixed.
            A normalized value needs the element's value in every sample or group. Without the
            cache, the map of each sample or group would work out the values of all the others
            again. The range-finding pass fills the cache beforehand. The drawing passes then have
            nothing to recompute. Pass None to rescale every element afresh. The element values
            themselves stay cached on the layer ('_summarize_entry_values').

        Returns
        =======
        Dict[str, float]
            The normalized value in each sample or group. It holds only the samples or groups that
            the normalization defined a value for. It is empty if the element has no value in any
            sample or group.
        """
        key = None
        if cache is not None:
            key = tuple(Mapper._get_entry_kegg_ids(entry, layer['use_reaction_attribute']))
            if key in cache:
                return cache[key]
        normalized_values = Mapper._normalize_element_values(
            Mapper._summarize_entry_values(entry, layer)[0], individuals, normalize
        )
        if cache is not None:
            cache[key] = normalized_values
        return normalized_values

    @staticmethod
    def _normalize_element_values(
        element_values: Dict[str, float],
        individuals: Iterable[str],
        normalize
    ) -> Dict[str, float]:
        """
        Rescale a map element's values in the samples or groups against each other.

        Parameters
        ==========
        element_values : Dict[str, float]
            The element's value in each sample or group that has one ('_summarize_element_values').

        individuals : Iterable[str]
            The samples or groups to normalize across, in the order they are colored in.

        normalize : callable
            Rescales an element's values across the samples or groups that have one
            ('_resolve_element_normalization').

        Returns
        =======
        Dict[str, float]
            The normalized value in each sample or group. It holds only the samples or groups that
            the normalization defined a value for. It is empty if the element has no value in any of
            the samples or groups.
        """
        values: Dict[str, float] = {
            individual: element_values[individual] for individual in individuals
            if individual in element_values
        }
        if not values:
            return {}
        # A normalization can be undefined for the values of this element even where it is defined
        # elsewhere on the map (a z-score needs more than one value, and a ratio needs a positive
        # reference), in which case the element has no value and is left uncolored.
        normalized = normalize(list(values.values()))
        return {
            individual: float(value)
            for individual, value in zip(values, normalized) if np.isfinite(value)
        }

    @staticmethod
    def _take_single_value(values: List[float]) -> float:
        """
        Return the one value of a slice, which stands for a single accession ('_slice_values').

        Parameters
        ==========
        values : List[float]
            The values of the slice's one accession in a sample. There is one.

        Returns
        =======
        float
            That value.
        """
        return float(values[0])

    @staticmethod
    def _gene_members(unit: Tuple[str, Tuple[str, ...]]) -> Tuple[str, ...]:
        """
        Name the rows of the file that a slice of a gene stands for.

        A member is the row of a gene for one accession. It is named by the gene ID and the
        accession, joined by a tab ('_build_gene_slice_keys'). The members are sorted, so the same
        rows get the same name on every map. The value of a slice depends on which rows of its gene
        fall in the element, so the members name the slice in the caches.

        Parameters
        ==========
        unit : Tuple[str, Tuple[str, ...]]
            A gene and its accessions in an element ('_gene_element').

        Returns
        =======
        Tuple[str, ...]
            The members.
        """
        gene, accessions = unit
        return tuple(sorted(f'{gene}\t{accession}' for accession in accessions))

    @staticmethod
    def _gene_view(layer: dict, cache: dict) -> dict:
        """
        Make a view of a layer sliced by gene, in which members take the place of accessions.

        '_summarize_element_values' and '_summarize_element_individuals' read the values of
        accessions from a layer. The view gives them the values of members instead
        ('_gene_members'). A slice of a gene is then summarized as an element whose accessions are
        its members. The view shares the rules, summaries, seed, salt, choices, and undefined sets
        of the layer.

        Parameters
        ==========
        layer : dict
            A layer model sliced by gene ('_build_gene_slice_keys').

        cache : dict
            The cache of the view, which takes the place of the layer's '_element_cache'.

        Returns
        =======
        dict
            The view.
        """
        view = {
            key: layer[key] for key in (
                'group_samples', 'sample_aggregate', 'group_aggregate', 'value_order', 'seed',
                'choice_salt', 'accession_rule', 'sample_rule', 'group_rule', '_choices',
                '_undefined_individual', '_undefined_unified'
            ) if key in layer
        }
        view['sample_values'] = layer['gene_sample_values']
        view['unified_values'] = layer['gene_unified_values']
        view['_element_cache'] = cache
        return view

    @staticmethod
    def _gene_element(accessions: Tuple[str, ...], layer: dict) -> dict:
        """
        Find the slices of genes that a map element can have, and their weights.

        A gene can have a slice where it has a row for any of the element's accessions. A gene with
        rows for two or more of them has one slice. The slices are ordered by the position of each
        gene's first accession in the element, then by gene ID in natural sort order. The weight of
        a gene in a sample is the mean weight of its rows for the element's accessions in that
        sample, so a gene with two accessions does not take twice the width. Results are kept in the
        layer's '_gene_elements', keyed by the accessions.

        Parameters
        ==========
        accessions : Tuple[str, ...]
            The accessions of the element, each once, in KGML order ('_get_entry_kegg_ids').

        layer : dict
            A layer model sliced by gene ('_build_gene_slice_keys').

        Returns
        =======
        dict
            These keys:
            - 'units': the slices in order, each as a gene and its accessions in the element in KGML
              order.
            - 'sample_weights': the weight of each slice in each sample, for a layer with samples
              and weights. It is None otherwise.
            - 'unified_weights': the weight of each slice, for a layer without samples and with
              weights. It is None otherwise.
        """
        cache = layer['_gene_elements']
        if accessions in cache:
            return cache[accessions]
        genes: Dict[str, List[str]] = {}
        for accession in accessions:
            for gene in layer['accession_genes'].get(accession, ()):
                genes.setdefault(gene, []).append(accession)
        units = tuple((gene, tuple(gene_accessions)) for gene, gene_accessions in genes.items())
        members = [Mapper._gene_members(unit) for unit in units]

        def _weights(row_weights: Dict[str, float]) -> Dict[Tuple[str, Tuple[str, ...]], float]:
            unit_weights = {}
            for unit, unit_members in zip(units, members):
                found = [row_weights[member] for member in unit_members if member in row_weights]
                if found:
                    unit_weights[unit] = sum(found) / len(found)
            return unit_weights

        element = {'units': units, 'sample_weights': None, 'unified_weights': None}
        if layer['gene_sample_weights'] is not None:
            element['sample_weights'] = {
                sample: _weights(row_weights)
                for sample, row_weights in layer['gene_sample_weights'].items()
            }
        elif layer['gene_unified_weights'] is not None:
            element['unified_weights'] = _weights(layer['gene_unified_weights'])
        cache[accessions] = element
        return element

    @staticmethod
    def _gene_standin_individuals(
        units: Tuple[Tuple[str, Tuple[str, ...]], ...],
        layer: dict,
        individual: Union[str, None],
        count_choices: bool = True
    ) -> Tuple[Dict[str, str], Union[str, None]]:
        """
        Find the categorical stand-in value of an element sliced by gene.

        The element stands for its drawn slices. In each sample, each gene has the value that the
        accession rule chose among its rows ('_slice_values'). The accession rule then chooses among
        the values of the genes. Each gene gives one vote. The sample and group rules choose across
        the samples and groups ('_summarize_element_individuals'). Results are kept in the layer's
        '_gene_standin_cache'. A gene is named there by its members, since its value depends on
        which of its rows fall in the element.

        Parameters
        ==========
        units : Tuple[Tuple[str, Tuple[str, ...]], ...]
            The drawn slices, each as a gene and its accessions in the element ('_gene_element').

        layer : dict
            A layer model of categorical values sliced by gene ('_build_gene_slice_keys').

        individual : Union[str, None]
            The sample or group whose map asks for the value, or None for the 'unified' map.

        count_choices : bool, True
            False to count no choice ('_summarize_element_individuals').

        Returns
        =======
        Tuple[Dict[str, str], Union[str, None]]
            The stand-in value in each sample or group, and on the 'unified' map, in the form of
            '_summarize_element_individuals'.
        """
        names = {'\n'.join(members): members for members in map(Mapper._gene_members, units)}
        cache = layer['_gene_standin_cache']
        view = Mapper._gene_view(layer, cache)
        # The votes of the genes are needed only where the stand-in is not yet in the cache, which
        # is keyed by the sorted names. The value of each gene repeats the choice made for its own
        # slice. That choice was counted when the map asked for the slice.
        if tuple(sorted(names)) not in cache:
            def _gene_values(values: Dict[str, str], sample: str) -> Dict[str, str]:
                chosen = {}
                for name, members in names.items():
                    votes = [values[member] for member in members if member in values]
                    if not votes:
                        continue
                    value, _ = Mapper._choose_categorical_value(
                        votes, layer['accession_rule'], layer['value_order'], layer['seed'],
                        layer['choice_salt'] + ('accession', sample) + members
                    )
                    if value is not None:
                        chosen[name] = value
                return chosen
            if view['sample_values'] is None:
                view['unified_values'] = _gene_values(view['unified_values'], '')
            else:
                view['sample_values'] = {
                    sample: _gene_values(values, sample)
                    for sample, values in view['sample_values'].items()
                }
        return Mapper._summarize_element_individuals(tuple(names), view, individual, count_choices)

    @staticmethod
    def _slice_values(
        accessions: Tuple[str, ...],
        layer: dict,
        context: Union[str, None],
        count_choices: bool = True
    ) -> Dict[Union[str, Tuple[str, Tuple[str, ...]]], Union[float, str]]:
        """
        Find the value of each slice of a map element on one map.

        A slice of an accession takes the value of the accession, summarized as if the accession
        were a map element of its own. The accession aggregation is never applied to it. Within a
        sample, the value is the accession's value, to which the gene aggregation has already
        reduced its genes. The sample and group summaries then pool these values, as they pool the
        values of elements. A numeric slice is summarized by '_summarize_element_values' with
        '_take_single_value', and its result is kept in the layer's '_slice_cache'. The element
        cache cannot hold it, since a one-accession element and a slice share a key, and an
        element's value comes from the accession aggregation. A categorical slice is summarized by
        '_summarize_element_individuals'. It shares the layer's '_element_cache' with the elements.
        A slice and a one-accession element then get the same value, and their choices are counted
        once.

        A slice of a gene is summarized as an element whose accessions are its members, the rows of
        the gene for the accessions of the element ('_gene_members', '_gene_view'). Within a sample,
        a numeric slice takes the value that 'gene_aggregate' reduces these rows to. A categorical
        slice takes the value that the accession rule chooses among them. The sample and group
        summaries then pool these values. A sample in which the gene has no row for the element is
        left out. The results are kept in the layer's '_gene_cache', keyed by the members, so a
        slice never takes a result found for other rows.

        Parameters
        ==========
        accessions : Tuple[str, ...]
            The accessions of the element, each once, in KGML order ('_get_entry_kegg_ids').

        layer : dict
            A sliced layer model of numeric or categorical values ('_build_slice_model').

        context : Union[str, None]
            The sample or group whose map is drawn, or None for the 'unified' map. A layer without
            samples shows the same values on every map. It takes None for every map.

        count_choices : bool, True
            False to count no choice among categorical values ('_summarize_element_individuals').

        Returns
        =======
        Dict[Union[str, Tuple[str, Tuple[str, ...]]], Union[float, str]]
            The value of each slice that has one on the map, in the order of the slices. A slice of
            an accession is keyed by the accession. A slice of a gene is keyed by the gene and its
            accessions in the element ('_gene_element'). A slice whose summary is undefined has no
            value.
        """
        categorical = layer.get('value_type') == 'categorical'
        values: Dict[Union[str, Tuple[str, Tuple[str, ...]]], Union[float, str]] = {}
        if layer['slice_by'] == 'gene':
            view = Mapper._gene_view(layer, layer['_gene_cache'])
            for unit in Mapper._gene_element(accessions, layer)['units']:
                members = Mapper._gene_members(unit)
                if categorical:
                    individual_values, unified_value = Mapper._summarize_element_individuals(
                        members, view, context, count_choices
                    )
                else:
                    individual_values, unified_value = Mapper._summarize_element_values(
                        members, view, layer['gene_aggregate'], layer['_gene_cache']
                    )
                value = unified_value if context is None else individual_values.get(context)
                if value is not None:
                    values[unit] = value
            return values
        for accession in accessions:
            if accession not in layer['accessions']:
                continue
            if categorical:
                individual_values, unified_value = Mapper._summarize_element_individuals(
                    (accession,), layer, context, count_choices
                )
            else:
                individual_values, unified_value = Mapper._summarize_element_values(
                    (accession,), layer, Mapper._take_single_value, layer['_slice_cache']
                )
            value = unified_value if context is None else individual_values.get(context)
            if value is not None:
                values[accession] = value
        return values

    @staticmethod
    def _slice_shares(
        accessions: Tuple[str, ...],
        layer: dict,
        context: Union[str, None],
        drawn: List[str]
    ) -> Dict[str, float]:
        """
        Find the share of each slice of a map element on one map, before the shares are scaled.

        Without weights, every slice takes the same share. With weights, w(s, a) is the weight of
        accession a in sample s ('_build_slice_model'). A layer without samples gives each accession
        its one weight. The map of a sample gives each accession its weight in that sample. A map
        that pools samples, which is the map of a group or the 'unified' map, pools the weights as
        'slice_share_summary' says ('SLICE_SHARE_SUMMARIES').

        'share_mean' takes the shares of each sample, p(s, a) = w(s, a) / T(s). T(s) is the total
        weight of all of the element's accessions in that sample. An accession absent from a sample
        has a share of 0 there. The map of a group takes the mean of these shares over the group's
        samples where T(s) > 0. Ungrouped, the 'unified' map takes the mean over the samples where
        T(s) > 0. Grouped, it takes the mean of the groups' shares over the groups that have shares.

        'pooled_sum' adds up the weights of each accession over the samples of the group, or over
        every sample for the 'unified' map.

        Slices of genes take the place of accessions in all of this. '_slice_result' then passes the
        genes of the element and their weights in this element ('_gene_element'). T(s) is the total
        weight of all of the element's genes in s.

        Parameters
        ==========
        accessions : Tuple[str, ...]
            The accessions of the element, each once ('_get_entry_kegg_ids'). For slices of genes,
            these are the genes of the element ('_gene_element').

        layer : dict
            A sliced layer model ('_build_slice_model'). For slices of genes, its weights are those
            of the genes in the element.

        context : Union[str, None]
            The sample or group whose map is drawn, or None for the 'unified' map, as in
            '_slice_values'.

        drawn : List[str]
            The accessions with a value on the map ('_slice_values'). Only these get shares.

        Returns
        =======
        Dict[str, float]
            The share of each accession of 'drawn'. Each share is 0 or more. The shares need not sum
            to 1 ('_slice_result').
        """
        if not layer['slice_weights']:
            return {accession: 1.0 for accession in drawn}
        if layer['sample_values'] is None:
            weights = layer['unified_weights']
            return {accession: weights.get(accession, 0.0) for accession in drawn}
        sample_weights = layer['sample_weights']
        group_samples = layer['group_samples']
        if context is not None and group_samples is None:
            weights = sample_weights[context]
            return {accession: weights.get(accession, 0.0) for accession in drawn}

        def _pool(samples: Iterable[str]) -> Union[Dict[str, float], None]:
            # Pool the weights of a set of samples. None means that no sample has weight here.
            if layer['slice_share_summary'] == 'pooled_sum':
                return {
                    accession: sum(sample_weights[sample].get(accession, 0.0) for sample in samples)
                    for accession in accessions
                }
            sample_shares = []
            for sample in samples:
                weights = sample_weights[sample]
                total = sum(weights.get(accession, 0.0) for accession in accessions)
                if total > 0:
                    sample_shares.append({
                        accession: weights.get(accession, 0.0) / total for accession in accessions
                    })
            if not sample_shares:
                return None
            return {
                accession: sum(shares[accession] for shares in sample_shares) / len(sample_shares)
                for accession in accessions
            }

        if context is not None:
            pooled = _pool(group_samples[context])
        elif group_samples is None or layer['slice_share_summary'] == 'pooled_sum':
            pooled = _pool(sample_weights)
        else:
            group_shares = [
                shares for shares in (_pool(samples) for samples in group_samples.values())
                if shares is not None
            ]
            pooled = {
                accession: sum(shares[accession] for shares in group_shares) / len(group_shares)
                for accession in accessions
            } if group_shares else None
        if pooled is None:
            return {accession: 0.0 for accession in drawn}
        return {accession: pooled[accession] for accession in drawn}

    @staticmethod
    def _slice_result(
        accessions: Tuple[str, ...],
        layer: dict,
        context: Union[str, None],
        count_choices: bool = True
    ) -> Union[List[Tuple[str, Union[float, str], float]], None]:
        """
        Find the slices of a map element on one map.

        An accession or gene gets a slice where it has a value on the map ('_slice_values') and a
        share above 0 ('_slice_shares'). The shares of these slices are then scaled to sum to 1.
        Slices of accessions follow the order of the accessions in the KGML entry. Slices of genes
        follow the order of '_gene_element'. Results are kept in the layer's '_slice_results'. The
        range of values is found from them, and the maps are drawn from them. A result found without
        counting the choices among categorical values is not kept, so that the map that draws it
        still counts them.

        Parameters
        ==========
        accessions : Tuple[str, ...]
            The accessions of the element, each once, in KGML order ('_get_entry_kegg_ids').

        layer : dict
            A sliced layer model ('_build_slice_model').

        context : Union[str, None]
            The sample or group whose map is drawn, or None for the 'unified' map, as in
            '_slice_values'.

        count_choices : bool, True
            False to count no choice among categorical values, as in '_slice_values'.

        Returns
        =======
        Union[List[Tuple[str, Union[float, str], float]], None]
            None where no accession or gene of the element has a value on the map. Otherwise the
            slices, as (unit, value, share) tuples. The unit is the accession, or the gene and its
            accessions in the element ('_slice_values'). The list is empty where every slice with a
            value has a share of 0. The element is then left uncolored.
        """
        key = (accessions, context)
        cache = layer['_slice_results']
        if key in cache:
            return cache[key]
        if not count_choices and layer.get('value_type') == 'categorical':
            cache = {}
        values = Mapper._slice_values(accessions, layer, context, count_choices)
        if not values:
            cache[key] = None
            return None
        units, shares_layer = accessions, layer
        if layer['slice_by'] == 'gene':
            element = Mapper._gene_element(accessions, layer)
            units = element['units']
            if layer['slice_weights']:
                shares_layer = {
                    'slice_weights': True,
                    'sample_values': layer['sample_values'],
                    'group_samples': layer.get('group_samples'),
                    'slice_share_summary': layer['slice_share_summary'],
                    'sample_weights': element['sample_weights'],
                    'unified_weights': element['unified_weights']
                }
        shares = Mapper._slice_shares(units, shares_layer, context, list(values))
        total = sum(share for share in shares.values() if share > 0)
        cache[key] = [
            (unit, value, shares[unit] / total)
            for unit, value in values.items() if shares[unit] > 0
        ]
        return cache[key]

    @staticmethod
    def _slice_standin_value(
        values: List[float],
        shares: List[float],
        period: Union[float, None] = None
    ) -> Union[float, None]:
        """
        Find the numeric value that stands in for the slices of a map element.

        An element with more slices than fit is drawn whole in the color of this value
        ('_stage_entry_slices'). The value also gives the element's one vote for the color of the
        compounds its reactions touch. It is the mean of the slice values, weighted by their shares.
        Values that repeat after a period are averaged on a circle ('_circular_mean').

        Parameters
        ==========
        values : List[float]
            The values of the slices.

        shares : List[float]
            The shares of the slices, which sum to 1.

        period : Union[float, None], None
            The period after which the values repeat, or None.

        Returns
        =======
        Union[float, None]
            The value, or None where a circular mean is undefined.
        """
        if period is None:
            # The mean lies between the smallest and the largest value. Rounding could put it
            # slightly outside them, which would widen a scale found from the values of the slices.
            mean = float(np.average(values, weights=shares))
            return min(max(mean, min(values)), max(values))
        value = Mapper._circular_mean(values, period, weights=shares)
        return value if np.isfinite(value) else None

    @staticmethod
    def _slice_period(layer: dict, unified: bool) -> Union[float, None]:
        """
        Return the period over which the stand-in values of a sliced layer are averaged, or None.

        Values that repeat after a period are averaged on a circle ('_slice_standin_value'). A
        spread on the 'unified' map is a plain number, so it is averaged as one. Categorical values
        are not averaged.

        Parameters
        ==========
        layer : dict
            A sliced layer model ('_build_slice_model').

        unified : bool
            True for the 'unified' map.

        Returns
        =======
        Union[float, None]
            The period, or None where values are averaged as plain numbers or not at all.
        """
        if layer.get('value_type') == 'categorical' or (unified and layer.get('unified_spread')):
            return None
        return layer.get('value_period')

    def _slicer(
        self,
        layer: dict,
        individual: Union[str, None] = None
    ) -> Callable[[kgml.Entry], Union[dict, None]]:
        """
        Build the function that slices the Entries of one layer on one kind of map.

        See '_draw_map_elements' for the slicer contract. A numeric slice takes the color of its
        value on the scale of the map ('_quantitative_color_priority'). A categorical slice takes
        the color of its value ('value_color_priorities').

        The stand-in of an element is worked out only where it is asked for, since a categorical
        stand-in counts the choices that make it ('_report_categorical_choices'). A numeric stand-in
        is the share-weighted mean of the slice values ('_slice_standin_value'). A categorical
        stand-in is the value that '_summarize_element_individuals' gives the accessions of the
        drawn slices on the map, as if they formed an element of their own. With slices of genes,
        each drawn gene gives one vote in each sample instead ('_gene_standin_individuals'). An
        accession or gene without a slice on the map does not vote in either case.

        Parameters
        ==========
        layer : dict
            A sliced layer model ('_build_slice_model').

        individual : Union[str, None], None
            The sample or group whose map is drawn, or None for the 'unified' map.

        Returns
        =======
        Callable[[kgml.Entry], Union[dict, None]]
            The slicer.
        """
        # A layer without samples shows the same slices on every map.
        context = None if layer['sample_values'] is None else individual
        categorical = layer.get('value_type') == 'categorical'
        if categorical:
            value_color_priorities = layer['value_color_priorities']

            def color_priority(value) -> Tuple[str, float]:
                return value_color_priorities[value]
        else:
            norm, cmap, center = (
                (layer['_unified_norm'], layer['cmap'], layer['_unified_center'])
                if individual is None else
                (layer['_individual_norm'], layer['individual_cmap'], layer['_individual_center'])
            )

            def color_priority(value) -> Tuple[str, float]:
                return Mapper._quantitative_color_priority(
                    value, norm, cmap, layer['reverse_overlay'], center
                )
        period = Mapper._slice_period(layer, individual is None)

        def slicer(entry: kgml.Entry) -> Union[dict, None]:
            accessions = tuple(
                Mapper._get_entry_kegg_ids(entry, layer['use_reaction_attribute'])
            )
            result = Mapper._slice_result(accessions, layer, context)
            if result is None:
                return None
            slices = [
                (unit, value, share) + tuple(color_priority(value))
                for unit, value, share in result
            ]

            def standin() -> Union[Tuple[Union[float, str], str, float], None]:
                if categorical:
                    # Only the accessions or genes of the drawn slices vote, as in the numeric mean.
                    drawn = tuple(unit for unit, _, _ in result)
                    individual_values, unified_value = (
                        Mapper._gene_standin_individuals(drawn, layer, context)
                        if layer['slice_by'] == 'gene' else
                        Mapper._summarize_element_individuals(drawn, layer, context)
                    )
                    value = unified_value if context is None else individual_values.get(context)
                else:
                    value = Mapper._slice_standin_value(
                        [value for _, value, _ in result], [share for _, _, share in result],
                        period
                    )
                if value is None:
                    return None
                return (value,) + tuple(color_priority(value))
            return {'accessions': accessions, 'slices': slices, 'standin': standin}
        return slicer

    def _ring_result(
        self,
        pathway: kgml.Pathway,
        entry: kgml.Entry,
        layer: dict,
        individual: Union[str, None],
        count_choices: bool = True
    ) -> Union[dict, None]:
        """
        Find what the sliced ring of a box shows on one map ('--reaction-border-slice-by').

        The range of the border scale and the staging of rings both ask this, so that the scale
        spans the values that the rings show. The rings of a layer are aligned with the fill or
        sliced on their own ('slice_layout').

        Aligned with the fill, a ring follows the fill of its box on the map. Where the fill is
        drawn in slices, each column takes the border value of its own accession or gene, as bands
        above and below it ('kgml.SlicedKGMLCanvas'). A column whose unit has no border value takes
        no bands. Where the fill has one slice, the ring is whole in the border value of that slice.
        Where the fill has more slices than fit ('_slice_cap'), the ring is whole in the color of a
        stand-in value. That value stands for the units with both a slice of the fill and a border
        value, with the shares of the fill. Where the fill is not sliced, since the map shows
        presence or the box has no fill value there, the ring is whole in the color of a stand-in
        value. That value stands for every unit with a border value, with equal shares. A unit with
        a border value but no slice of the fill is not shown.

        Sliced on its own, a ring is divided among the units with a border value and a share above
        0, in the order of their slices ('_slice_result'). Each unit takes one segment, clockwise
        from the top left corner. A ring with one slice is whole in the color of that slice. A ring
        with more than 'MAX_RING_SLICES' slices is whole in the color of the stand-in value of all
        of them.

        A stand-in value is found as for a fill drawn whole ('_slice_standin_value',
        '_summarize_element_individuals', '_gene_standin_individuals'). A ring whose stand-in value
        is undefined is not drawn. The choices among categorical values are counted for the slices
        that a ring shows or lists in the key table, and for its stand-in. A unit that the ring
        leaves out, by a share of 0 or for want of a column of the fill, counts no choice.

        Parameters
        ==========
        pathway : kgml.Pathway
            The map.

        entry : kgml.Entry
            An ortholog Entry of the map with a box that shows rings ('_shows_border').

        layer : dict
            The model of the fill of a reaction layer. Its 'border' holds the sliced model of its
            border ('_build_slice_model').

        individual : Union[str, None]
            The sample or group whose map is drawn, or None for the 'unified' map. The border shows
            its values on that map.

        count_choices : bool, True
            False to count no choice among categorical values. The range of a border asks this way,
            since no map is drawn then.

        Returns
        =======
        Union[dict, None]
            None where the box has no value to show in its ring. Otherwise a dictionary with these
            keys:
            - 'kind': 'bands', 'whole', or 'segments'. It is 'uncolored' where every slice with a
              value has a share of 0, and 'none' where only units that are not shown have border
              values. The ring is not drawn in these two cases.
            - 'slices': the rows of the key table, as (position, unit, value, share) tuples. The
              position counts from 1. It is None in a ring drawn whole for a reason.
            - 'values': the border value of each column of the fill, or None for a column without
              bands. Only bands have this key.
            - 'value': the value of a whole ring. Only a whole ring has this key.
            - 'reason': why a ring is drawn whole, such as 'over the cap' or 'fill not sliced', or
              None.
            - 'unshown': the units with a border value but no slice of the fill.
            - 'accessions': the accessions of the Entry.
            A result of the kind 'uncolored' or 'none' has only 'kind', 'unshown', and
            'accessions'.
        """
        border = layer['border']
        context = None if border['sample_values'] is None else individual
        accessions = tuple(self._get_entry_kegg_ids(entry, border['use_reaction_attribute']))
        categorical = border.get('value_type') == 'categorical'
        period = self._slice_period(border, individual is None)

        def _standin(slices) -> Union[float, str, None]:
            # The value that stands in for the slices of a ring drawn whole. Only the units of these
            # slices vote, as for a fill drawn whole ('_slicer').
            if categorical:
                units = tuple(unit for unit, _, _ in slices)
                individual_values, unified_value = (
                    self._gene_standin_individuals(units, border, context, count_choices)
                    if border['slice_by'] == 'gene' else
                    self._summarize_element_individuals(units, border, context, count_choices)
                )
                return unified_value if context is None else individual_values.get(context)
            return self._slice_standin_value(
                [value for _, value, _ in slices], [share for _, _, share in slices], period
            )

        def _whole(slices, reason, unshown=()) -> Union[dict, None]:
            # A ring of one slice takes the value of that slice. Any other ring takes a stand-in.
            value = slices[0][1] if len(slices) == 1 else _standin(slices)
            if value is None:
                return None
            return {
                'kind': 'whole', 'value': value, 'reason': reason, 'unshown': list(unshown),
                'accessions': accessions,
                'slices': [
                    (None if reason else position, unit, slice_value, share)
                    for position, (unit, slice_value, share) in enumerate(slices, 1)
                ]
            }

        # The choices among categorical values are counted only for the slices that a ring shows,
        # once they are known ('_count_slice_choices').
        if border['slice_layout'] == 'independent':
            slices = self._slice_result(accessions, border, context, count_choices=False)
            if slices is None:
                return None
            if not slices:
                return {'kind': 'uncolored', 'unshown': [], 'accessions': accessions}
            if count_choices:
                self._count_slice_choices([unit for unit, _, _ in slices], border, context)
            if len(slices) > MAX_RING_SLICES:
                return _whole(slices, 'over the cap')
            if len(slices) == 1:
                return _whole(slices, None)
            return {
                'kind': 'segments', 'reason': None, 'unshown': [], 'accessions': accessions,
                'slices': [
                    (position, unit, value, share)
                    for position, (unit, value, share) in enumerate(slices, 1)
                ]
            }

        # Aligned rings follow the fill of the box on this map. The fill is sliced only on the maps
        # that show its values ('_sliced').
        mode = layer['unified_mode'] if individual is None else layer['individual_mode']
        fill_slices = None
        if mode in ('quantitative', 'categorical'):
            fill_slices = self._slice_result(accessions, layer, context, count_choices)
        border_values = self._slice_values(accessions, border, context, count_choices=False)
        if not border_values:
            return None
        if fill_slices is not None and not fill_slices:
            return {'kind': 'uncolored', 'unshown': [], 'accessions': accessions}
        if fill_slices is None:
            if count_choices:
                self._count_slice_choices(list(border_values), border, context)
            share = 1.0 / len(border_values)
            return _whole(
                [(unit, value, share) for unit, value in border_values.items()], 'fill not sliced'
            )
        units = [unit for unit, _, _ in fill_slices]
        unshown = [unit for unit in border_values if unit not in units]
        shown = [
            (position, unit, border_values[unit], share)
            for position, (unit, _, share) in enumerate(fill_slices, 1) if unit in border_values
        ]
        if not shown:
            return {'kind': 'none', 'unshown': unshown, 'accessions': accessions}
        if count_choices:
            self._count_slice_choices([unit for _, unit, _, _ in shown], border, context)
        reason = self._slice_cap(
            pathway, entry, len(fill_slices), ringed=True
        ) if len(fill_slices) > 1 else None
        if len(fill_slices) > 1 and reason is None:
            return {
                'kind': 'bands', 'reason': None, 'unshown': unshown, 'accessions': accessions,
                'values': [border_values.get(unit) for unit in units], 'slices': shown
            }
        return _whole([(unit, value, share) for _, unit, value, share in shown], reason, unshown)

    @staticmethod
    def _count_slice_choices(
        units: List[Union[str, Tuple[str, Tuple[str, ...]]]],
        layer: dict,
        context: Union[str, None]
    ) -> None:
        """
        Count the choices among categorical values that some slices use on one map.

        The values of the slices were found without counting their choices ('_slice_values'). They
        are kept in the caches of the layer, so asking for them again only counts the choices
        ('_summarize_element_individuals'). Nothing is counted for numeric values.

        Parameters
        ==========
        units : List[Union[str, Tuple[str, Tuple[str, ...]]]]
            The units of the slices, as in '_slice_values'.

        layer : dict
            A sliced layer model ('_build_slice_model').

        context : Union[str, None]
            The sample or group whose map is drawn, or None for the 'unified' map, as in
            '_slice_values'.
        """
        if layer.get('value_type') != 'categorical':
            return
        if layer['slice_by'] == 'gene':
            view = Mapper._gene_view(layer, layer['_gene_cache'])
            for unit in units:
                Mapper._summarize_element_individuals(Mapper._gene_members(unit), view, context)
            return
        for unit in units:
            Mapper._summarize_element_individuals((unit,), layer, context)

    def _ring_values(
        self,
        pathway: kgml.Pathway,
        entry: kgml.Entry,
        layer: dict,
        individual: Union[str, None]
    ) -> List[Union[float, str]]:
        """
        List the values that the sliced ring of a box shows on one map, without counting choices.

        Parameters
        ==========
        pathway, entry, layer, individual
            As in '_ring_result'.

        Returns
        =======
        List[Union[float, str]]
            The value of each band or segment, or the one value of a whole ring. The list is empty
            where the box shows no ring.
        """
        ring = self._ring_result(pathway, entry, layer, individual, count_choices=False)
        if ring is None:
            return []
        if ring['kind'] == 'bands':
            return [value for value in ring['values'] if value is not None]
        if ring['kind'] == 'segments':
            return [value for _, _, value, _ in ring['slices']]
        if ring['kind'] == 'whole':
            return [ring['value']]
        return []

    def _ringer(
        self,
        layer: dict,
        individual: Union[str, None] = None
    ) -> Callable[[kgml.Pathway, kgml.Entry], Union[dict, None]]:
        """
        Build the function that finds the sliced rings of one layer on one kind of map.

        A numeric value takes the color of its value on the border scale of the map
        ('_quantitative_color_priority'). A categorical value takes the color of its value
        ('value_color_priorities').

        Parameters
        ==========
        layer : dict
            The model of the fill of a reaction layer whose border is sliced ('_ring_result').

        individual : Union[str, None], None
            The sample or group whose map is drawn, or None for the 'unified' map.

        Returns
        =======
        Callable[[kgml.Pathway, kgml.Entry], Union[dict, None]]
            The function. It takes the map and an ortholog Entry with a box that shows rings. It
            returns what '_ring_result' returns. A ring that is drawn also gets these keys:
            - 'border': the ring, as 'Graphics.border' takes it.
            - 'colors': the color of each row of 'slices'. The rows of a ring drawn whole all take
              the color of the ring.
        """
        border = layer['border']
        if border.get('value_type') == 'categorical':
            value_color_priorities = border['value_color_priorities']

            def color(value) -> str:
                return value_color_priorities[value][0]
        else:
            norm, cmap, center = (
                (border['_unified_norm'], border['cmap'], border['_unified_center'])
                if individual is None else (
                    border['_individual_norm'], border['individual_cmap'],
                    border['_individual_center']
                )
            )

            def color(value) -> str:
                return Mapper._quantitative_color_priority(
                    value, norm, cmap, border['reverse_overlay'], center
                )[0]

        def ringer(pathway: kgml.Pathway, entry: kgml.Entry) -> Union[dict, None]:
            ring = self._ring_result(pathway, entry, layer, individual)
            if ring is None or ring['kind'] in ('uncolored', 'none'):
                return ring
            if ring['kind'] == 'whole':
                ring['border'] = color(ring['value'])
                ring['colors'] = [ring['border']] * len(ring['slices'])
                return ring
            ring['colors'] = [color(value) for _, _, value, _ in ring['slices']]
            if ring['kind'] == 'bands':
                ring['border'] = (
                    'bands', [None if value is None else color(value) for value in ring['values']]
                )
            else:
                ring['border'] = ('segments', [
                    (hexcode, share)
                    for hexcode, (_, _, _, share) in zip(ring['colors'], ring['slices'])
                ])
            return ring
        return ringer

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
        individuals: Union[List[str], None] = None,
        individual_noun: Union[str, None] = None,
        grid_source_type: str = 'unknown',
        colorbar_individual_suffix: Union[str, None] = None,
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
        'categorical' (a color per categorical value), 'membership' (by the sources/groups
        containing an element, or their count), 'single' (one fixed presence/absence color),
        'static' (one fixed color pooled across sources, for the db/pan single-color path),
        'original' (preserve the reference map's colors, drawn by a separate drawer), or 'none' (no
        fill, so that a text layer draws only the rings of its border column). The mode can differ
        between the two map contexts, given as 'unified_mode' for the 'unified' map and
        'individual_mode' for the individual sample/source/group maps: a text layer with a value
        column and a 'sample' column, for instance, colors per-sample magnitude continuously while
        summarizing the samples or groups by presence on the 'unified' map. Note that 'original' is
        a 'unified_mode' only. An individual map can be drawn in the reference colors for one
        SOURCE, since 'source_accessions' is keyed by source. It cannot be drawn that way for a
        group, so an original layer has an 'individual_mode' of 'membership'. Its grouped individual
        maps show within-group source counts.

        Because every layer is reduced to a per-Entry '(color, priority)' by '_draw_map_elements', a
        quantitative layer and a presence layer can be colored on the same map. A 'unified' map
        colors elements summarized across all samples, contigs databases, genomes, or groups. When
        they are drawn individually, each gets its own map. These individual maps share one colorbar
        per quantitative layer, so colors stay comparable.

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
            for the individual maps as 'individual_mode'. A database or pangenome layer gives one
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
            - 'value_column': the header of the column of the values, which messages name.
            - 'channel': the channel that the values color, the fill or the border ('_txt_channel').
              Messages name its options.

            A 'quantitative' layer can also have these keys:
            - 'individual_cmap': the colormap of the individual maps. It defaults to 'cmap'.
            - 'value_limits' and 'individual_value_limits': the limits of the 'unified' map's scale
              and of the individual maps' scale ('_make_quantitative_norm').
            - 'value_center' and 'individual_value_center': the centers of the same two scales.
            - 'individual_value_center_flag': the option that the individual maps' center came from.
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
            - 'element_normalize': rescales an element's value in each sample or group against its
              values in all of them ('_normalize_entry_value').
            - 'individual_colorbar_label': the label of the individual maps' colorbar. It defaults
              to 'colorbar_label'. A normalization makes the scale show a different quantity, so it
              needs its own label.
            - 'unified_colorbar_label': the label of the 'unified' map's colorbar. It defaults to
              'colorbar_label'. A spread makes that scale show a different quantity, so it needs its
              own label.

            A 'categorical' layer comes from '_build_categorical_txt_model'. It needs these keys:
            - 'value_type': 'categorical'.
            - 'accessions', 'sample_values', 'unified_values', 'group_samples', 'colorbar_label',
              'value_column', and 'channel', as for a 'quantitative' layer. The values are names.
            - 'value_color_priorities': the color and priority of each value, in the order of the
              values. It is None where no map shows the values.
            - 'value_order', 'gene_rule', 'accession_rule', 'sample_rule', 'group_rule', 'seed', and
              'choice_salt', which choose one value where several meet
              ('_summarize_entry_individuals'). 'sample_rule' and 'group_rule' are None where no
              rule was given.
            - '_element_cache' and '_choices', which keep and count those choices.

            A 'membership' layer needs these keys:
            - 'membership': the sources that contain each accession.
            - 'source_accessions': the accessions of each source.
            - 'color_hexcode': the color of an individual source's map.

            A 'membership' layer can also have these keys:
            - 'colormap', 'colormap_limits', 'colormap_scheme' and 'reverse_overlay', which set how
              presence is colored.
            - 'default_colormap_scheme': the scheme used where 'colormap_scheme' is None. Without
              it, the scheme depends on the number of samples, contigs databases, genomes, or groups
              ('_membership_layer_colors').
            - 'scheme_options': the option that chooses the presence scheme, which differs by input.
              It defaults to 'PRESENCE_SCHEME_OPTIONS'.
            - 'membership_colors': the color of each sample, contigs database, genome, or group, in
              place of a colormap. These colors color the membership scale. Without groups, each one
              also colors the individual map of its sample, contigs database, or genome. With
              groups, a group's maps take the group's color only with 'GROUP_COLORMAP_TINT'. They
              are checked against 'individuals' here.
            - 'membership_combo_colors': the colors of combinations of samples, contigs databases,
              genomes, or groups. They override the mix of 'membership_colors'. This key is needed
              with 'membership_colors'.
            - 'membership_colors_flag': the option that 'membership_colors' came from. It is needed
              with 'membership_colors'.

            A 'static' or 'original' layer needs 'membership', 'source_accessions' and
            'color_hexcode'. 'color_hexcode' is 'original' for an 'original' layer.

            A 'single' layer needs 'accessions' and 'color_hexcode'.

            A 'none' layer needs 'accessions' and a 'border'. It is 'none' in both contexts. Its
            matched entries do not count as data, so a map is drawn for it only where a ring shows.

            A layer that is 'quantitative' in one context and 'membership' in the other has the keys
            of both.

            A sliced reaction layer has the keys that '_build_slice_model' returns. Every drawn map
            that shows its values ('quantitative' or 'categorical') is sliced ('_slicer'). The range
            of each scale is found from the values of the slices drawn. The slices drawn are listed
            in a key table ('_write_slice_keys'). A border with sliced rings has these keys too. Its
            rings are found from the fill of its layer ('_ring_result').

            A text layer can also have 'border', the model of its border column, and 'bars', a list
            of the models of its bar columns, which holds one model. These have the keys of a
            'quantitative' layer, except that a bar has no colormap. Their mode is None on a map
            that shows none of their values. A bar also has 'style', how it is drawn in its box
            ('BAR_STYLES').

        output_dir : str
            Path to the output directory in which pathway map and colorbar PDF files are drawn.

        individuals : Union[List[str], None]
            The shared names of the samples, contigs databases, genomes, or groups that get
            individual maps, in color-assignment order. It is None in a run with a single source of
            data, which draws a single map.

        Notes
        =====
        'individual_noun'/'grid_source_type'/'colorbar_individual_suffix'/'subset_subject'/
        'unified_plural'/'membership_*'/'grouped_membership' carry the terminology and grouping that
        the public methods supply; the remaining parameters mirror those methods. 'grid_source_type'
        is the kind of source: 'contigs database', 'pangenome', or 'sample'. It labels the per-group
        grid colorbars, which count sources, and the grid's 'unified' panel. Its default, 'unknown',
        labels the colorbars 'source count' and the panel 'all'. Quantitative layers carry their own
        colorbar label in 'colorbar_label'.
        """
        has_individuals = individuals is not None
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
            if 'individual_mode' not in layer:
                layer['individual_mode'] = (
                    'membership' if layer['mode'] in ('static', 'original') else layer['mode']
                )
            # A colormap for the scale of the individual maps alone is optional. Without one, the
            # two contexts are colored from the layer's single 'cmap'.
            if layer['individual_mode'] == 'quantitative' and 'individual_cmap' not in layer:
                layer['individual_cmap'] = layer['cmap']

        # The border of a text layer is a model of its own, built from its border column
        # ('map_kegg_pathways_txt'). It has the keys of a 'quantitative' or 'categorical' layer. Its
        # mode is None on a map that shows no border. The bars of a text layer are models too, built
        # from its bar column. The steps that find, check, and key a scale take the borders and bars
        # with the layers. The checks of reserved colors and of derived compound colors do not,
        # since a ring is never part of the colors of an element, and a bar takes the color of its
        # fill.
        borders = [layer['border'] for layer in layers if layer.get('border') is not None]
        bars = [bar for layer in layers for bar in layer.get('bars', [])]
        models = layers + borders + bars

        def _stem(layer: dict) -> str:
            # The options of a layer from a text file start with the stem of its channel.
            return layer['channel']['flag_stem'] if 'channel' in layer else layer['element_type']

        def _noun(layer: dict) -> str:
            return layer['channel']['noun'] if 'channel' in layer else layer['element_type']

        def _is_border(layer: dict) -> bool:
            return 'channel' in layer and layer['channel']['role'] == 'border'

        def _role(layer: dict) -> str:
            # The channel of a model: the fill of its layer, or its border or bars.
            return layer['channel']['role'] if 'channel' in layer else 'fill'

        original_run = any(layer['unified_mode'] == 'original' for layer in layers)
        static = any(layer['unified_mode'] in ('static', 'original') for layer in layers)
        # Groups color individual maps by within-group source counts only for the layers whose
        # per-group context is presence; a group map colored by value gets its colors from the
        # layer's own colormap instead.
        grouped_presence = grouped and any(
            layer['individual_mode'] == 'membership' for layer in layers
        )

        def _unified_scale_drawn(layer: dict) -> bool:
            # A layer without samples of its own colors the individual maps from the 'unified'
            # scale. That scale is therefore still worked out where the 'unified' map is not drawn.
            # The options bounding and centering it still apply.
            return draw_unified_maps or (
                layer['individual_mode'] == 'quantitative' and layer['sample_values'] is None
            )

        def _dedup(items: List[str]) -> List[str]:
            seen: Set[str] = set()
            return [item for item in items if not (item in seen or seen.add(item))]

        if has_individuals:
            subset_names = set(individuals)
            self._check_requested_subset(
                draw_individual_files, "Individual maps", subset_names, subset_subject
            )
            self._check_requested_subset(
                draw_grid, "Individual maps in grids", subset_names, subset_subject
            )
            draw_files_individuals = (
                _dedup(list(individuals)) if draw_individual_files is True
                else [] if draw_individual_files is False
                else _dedup(list(draw_individual_files))
            )
            draw_grid_individuals = (
                _dedup(list(individuals)) if draw_grid is True
                else [] if draw_grid is False
                else _dedup(list(draw_grid))
            )
            draw_individuals = _dedup(draw_files_individuals + draw_grid_individuals)
            # Only the names of the samples, contigs databases, genomes, or groups drawn on their
            # own maps are joined onto a path, so only those names have to be usable as directory
            # names. A run that draws no individual maps or grids never writes the other names
            # anywhere. Neither does a run that asks for a subset of them. The 'unified' map
            # summarizes the others by color alone. The check comes after the subset checks, so a
            # name that is not a sample, contigs database, genome, or group of the run is reported
            # as unrecognized rather than as unusable.
            self._check_individual_names(draw_individuals, individual_noun)
            # A player of one frame would show the same map that the directory of the map holds, so
            # the flag would have no effect.
            if (
                self.collate_files_by_map and self.animate_files_by_map
                and len(draw_files_individuals) == 1
            ):
                raise ConfigError(
                    f"'--animate-files-by-map' writes a player that shows the map of each "
                    f"{individual_noun} in turn. Maps were asked for only one {individual_noun}, "
                    f"'{draw_files_individuals[0]}', so each player would hold a single frame. "
                    f"Please ask for the maps of more {individual_noun}s, or drop this flag."
                )
            # In an ordinary run, every sample, contigs database, genome, or group appears in the
            # output. At the least, it appears as color on the 'unified' map. Without that map, the
            # ones left out of a subset request are drawn nowhere at all.
            if not draw_unified_maps and len(draw_individuals) < len(set(individuals)):
                self.run.warning(
                    f"Maps were asked for {len(draw_individuals)} of the {len(set(individuals))} "
                    f"{individual_noun}s. The 'unified' map is not drawn in this run, so the other "
                    f"{individual_noun}s appear nowhere in the output."
                )
            # Without the 'unified' map's panel, a grid of one sample, contigs database, genome, or
            # group holds a single map. That is not an error, but it is worth saying what the file
            # will look like.
            if not draw_unified_maps and len(draw_grid_individuals) == 1:
                self.run.warning(
                    f"The 'unified' map is not drawn in this run, so no grid leads with its panel. "
                    f"Only one {individual_noun} was named for the grid. Each grid will hold a "
                    f"single map on a page. It is the same map that '--draw-individual-files' "
                    f"writes, with the {individual_noun}'s name over it."
                )
        else:
            draw_files_individuals = []
            draw_grid_individuals = []
            draw_individuals = []
        draw_individual_maps = has_individuals and (
            draw_individual_files is not False or draw_grid is not False
        )

        # The 'unified' map is what a run draws when it is asked for nothing else. Leaving it out
        # means the individual maps are the whole of the output. Without them there is no output at
        # all. Say so before an empty directory is created.
        if not draw_unified_maps and not draw_individuals:
            if has_individuals:
                raise ConfigError(
                    f"No maps would be drawn. The 'unified' map summarizing every "
                    f"{individual_noun} is not drawn here. Neither is a map for any individual "
                    f"{individual_noun}. Ask for those maps with '--draw-individual-files' and/or "
                    f"'--draw-grid', or let the 'unified' map be drawn."
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
        if not draw_individual_maps:
            # A normalization is checked before the limits and the center are, since one centered on
            # zero supplies a center of its own.
            normalized_layers = [
                layer for layer in models if layer.get('element_normalize') is not None
            ]
            if normalized_layers:
                flags = ', '.join(
                    f"'--{_stem(layer)}-element-normalization'" for layer in normalized_layers
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
                    'individual_value_limits', 'individual-value-limits',
                    'Limits were given for', 'bound'
                ),
                (
                    'individual_value_center', 'individual-value-center',
                    'A center was given for', 'center'
                )
            ):
                affected_layers = [
                    layer for layer in models if layer.get(model_key) is not None
                ]
                if not affected_layers:
                    continue
                flags = ', '.join(
                    f"'--{_stem(layer)}-{flag_suffix}'" for layer in affected_layers
                )
                unified_suffix = flag_suffix.replace('individual-', '')
                unified_flags = (
                    f"'--reaction-{unified_suffix}'/'--compound-{unified_suffix}'"
                    if all(_role(layer) == 'fill' for layer in affected_layers) else ', '.join(
                        f"'--{_stem(layer)}-{unified_suffix}'" for layer in affected_layers
                    )
                )
                # Bars are drawn on a scale of lengths, not of colors.
                scale = 'scale' if any(_role(layer) == 'bar' for layer in affected_layers) else (
                    'color scale'
                )
                raise ConfigError(
                    f"{subject_phrase} the {scale} shared by the maps of the individual "
                    f"samples or groups ({flags}), but this run draws no such maps, so that scale "
                    f"is never drawn and this would go nowhere. Ask for those maps with "
                    f"'--draw-individual-files' and/or '--draw-grid', or {remedy_phrase} the scale "
                    f"of the 'unified' map instead, with {unified_flags}."
                )

        # The options bounding and centering the 'unified' map's own scale are the mirror case. A
        # run leaving out that map never draws that scale. The exception is a layer without samples
        # of its own. Its individual maps take that very scale, so the options still apply.
        if not draw_unified_maps:
            for model_key, flag_suffix, subject_phrase, remedy_phrase in (
                ('value_limits', 'value-limits', 'Limits were given for', 'Bound'),
                ('value_center', 'value-center', 'A center was given for', 'Center')
            ):
                affected_layers = [
                    layer for layer in models
                    if layer.get(model_key) is not None and not _unified_scale_drawn(layer)
                ]
                if not affected_layers:
                    continue
                flags = ', '.join(
                    f"'--{_stem(layer)}-{flag_suffix}'" for layer in affected_layers
                )
                individual_flags = (
                    f"'--reaction-individual-{flag_suffix}'/'--compound-individual-{flag_suffix}'"
                    if all(_role(layer) == 'fill' for layer in affected_layers) else ', '.join(
                        f"'--{_stem(layer)}-individual-{flag_suffix}'" for layer in affected_layers
                    )
                )
                scale = 'scale' if any(_role(layer) == 'bar' for layer in affected_layers) else (
                    'color scale'
                )
                raise ConfigError(
                    f"{subject_phrase} the {scale} of the 'unified' map ({flags}). That map is "
                    f"not drawn here, so its scale is never drawn either. {remedy_phrase} the "
                    f"scale shared by the maps of the individual {individual_noun}s instead, with "
                    f"{individual_flags}, or let the 'unified' map be drawn."
                )

        # Before anything is drawn, the names in the membership colors file are checked against the
        # samples, contigs databases, genomes, or groups of this run. The colors are also checked
        # against what this run would actually do with them. A layer whose contexts are all colored
        # some other way would leave them unused, which is worth an error rather than a map that
        # quietly ignores the file it was handed.
        colored_layers = [layer for layer in layers if layer.get('membership_colors') is not None]
        if colored_layers:
            if not has_individuals:
                flags = ', '.join(
                    f"'{layer['membership_colors_flag']}'" for layer in colored_layers
                )
                raise ConfigError(
                    f"Membership colors were given ({flags}). A kegg-membership-colors-txt file "
                    f"gives the color of each sample, contigs database, genome, or group. This run "
                    f"has none of these to color. It has a single source of data, so it draws one "
                    f"map rather than one per sample, contigs database, genome, or group. Set the "
                    f"color of the single map with '--reaction-color'/'--compound-color' instead."
                )
            self._resolve_membership_colors(
                colored_layers, individuals, individual_noun, draw_unified_maps=draw_unified_maps
            )
            group_ramps = grouped and (
                grouped_membership['group_colormap'] == GROUP_COLORMAP_TINT
            )
            for layer in colored_layers:
                # These colors are used in two places. The first is the 'unified' map's membership
                # scale. That scale gives a band to each sample, contigs database, genome, or group.
                # It also gives a band to each combination of them. The second is the individual
                # maps. Without groups, the map of each sample, contigs database, or genome is drawn
                # in its own color. With groups, each group's own map is drawn in a ramp running to
                # the group's color, but only when '--group-colormap tint' asked for that ramp. The
                # two checks below skip a layer whose colors reach either place. A layer reaching
                # neither would leave its colors file unused.
                if draw_unified_maps and layer['unified_mode'] == 'membership':
                    continue
                if layer['individual_mode'] == 'membership' and draw_individual_maps and (
                    group_ramps if grouped else True
                ):
                    continue
                # The file colors the individual maps only for a layer without values, and only
                # without groups. A grouped run that reaches this point did not ask for tints.
                own = not grouped and layer['individual_mode'] == 'membership'
                own_maps = (
                    f" The file also colors the map of each {individual_noun}." if own else ""
                )
                grouped_clause = (
                    f" A group's own maps are colored by the number of the group's "
                    f"{membership_singular}s containing an element, in the colors of "
                    f"'--group-colormap'. Give '--group-colormap {GROUP_COLORMAP_TINT}' to build "
                    f"each group's scale from the group's color in the file."
                ) if grouped and layer['individual_mode'] == 'membership' else ""
                raise ConfigError(
                    f"The colors file given to '{layer['membership_colors_flag']}' gives a color "
                    f"to each {individual_noun}. On the 'unified' map, these colors show exactly "
                    f"which {individual_noun}s contain an element.{own_maps} Nothing in this run "
                    f"is colored {'either way' if own else 'that way'} for the "
                    f"{layer['element_type']} layer, so the file would go unused.{grouped_clause} "
                    f"Please remove the file, or color this layer by presence."
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
            # The static color is the fill's. A border keeps its own summaries.
            border_clause = (
                " The borders of such a layer are still colored by their own summaries."
                if any(layer.get('border') is not None for layer in static_layers) else ""
            )
            self.run.warning(
                f"Groups were provided, but these will be ignored for the 'unified' map, since a "
                f"static color (or the reference map's own colors) was set: dynamic coloring based "
                f"on membership in groups is overridden by static coloring based on "
                f"presence/absence in any {membership_singular}.{group_map_clause}{border_clause}"
            )

        # The overwrite check looks at every directory this run writes map files into. Checking
        # the 'unified' directory alone is not enough. A run skipping the 'unified' map writes no
        # such directory, so that check would pass while individual maps and grids were
        # overwritten.
        unified_dir = os.path.join(output_dir, UNIFIED_SUBDIR)
        check_dirs = [unified_dir] if draw_unified_maps else []
        check_dirs += [
            os.path.join(output_dir, INDIVIDUAL_SUBDIR, individual)
            for individual in draw_individuals
        ]
        if draw_grid is not False:
            check_dirs.append(os.path.join(output_dir, GRID_SUBDIR))
        # The players are the only files written under 'by_map' that are not links. The links of an
        # earlier run are replaced as a matter of course.
        extensions = None
        if self.collate_files_by_map and self.animate_files_by_map and draw_files_individuals:
            collated_dir = os.path.join(output_dir, COLLATED_SUBDIR)
            check_dirs.append(collated_dir)
            extensions = {collated_dir: PLAYER_EXTENSION}
        pathway_numbers = self._find_maps(
            check_dirs, patterns=pathway_numbers, extensions=extensions
        )
        # A border or bar that no map of the run can show is refused before any output is written.
        self._check_channel_maps(
            borders + bars, pathway_numbers, draw_unified_maps, draw_individual_maps,
            individual_noun,
            unfilled={layer['element_type'] for layer in layers if layer['unified_mode'] == 'none'}
        )
        filesnpaths.gen_output_directory(output_dir, progress=self.progress, run=self.run)
        # The key table of a sliced layer is written once the maps are drawn. Whether it can be
        # written is checked before anything is drawn. A border with sliced rings has a key table of
        # its own.
        for layer in models:
            if layer.get('slice_by') is not None:
                filesnpaths.is_output_file_writable(
                    os.path.join(output_dir, f"slice_key_{layer['name']}.txt"),
                    ok_if_exists=self.overwrite_output
                )

        drawn: Dict[Literal['unified', 'individual', 'grid'], Dict] = {
            'unified': {}, 'individual': {}, 'grid': {}
        }

        # Finalize the coloring model of each layer. A context colored quantitatively needs a value
        # range. A single parse of each map collects the values of each layer for whichever of its
        # contexts are quantitative. These are the unified values and the values in each sample or
        # group. The values of all samples or groups are pooled for the one scale that their maps
        # share. The parse uses the same extraction and aggregation as drawing.
        norm_layers = [
            layer for layer in models
            if 'quantitative' in (layer['unified_mode'], layer['individual_mode'])
        ]
        # A border scale spans only the rings that the maps show ('_shows_ring'). A categorical
        # border is keyed only where some map shows one of its rings, so the pass records those maps
        # for it. Every border records the reaction elements with border values that a standard map
        # draws as lines, which show no ring. A border with sliced rings finds what its rings show
        # from the fill of its layer ('_ring_result'). Bars are drawn in the boxes that take rings,
        # so a bar scale spans only those, and bars record their lines too. A map with bars records
        # the boxes it leaves uncolored for want of a bar value, and the boxes it draws without
        # their labels ('_stage_entry_bar').
        ring_layers = [
            layer for layer in borders
            if 'categorical' in (layer['unified_mode'], layer['individual_mode'])
        ]
        for layer in borders:
            layer['_ring_maps'] = {'unified': set(), 'individual': set()}
            layer['_line_entries'] = {}
        for layer in bars:
            layer['_line_entries'] = {}
            layer['_missing'] = {}
            layer['_unlabeled'] = {}
        fills = {id(layer['border']): layer for layer in layers if layer.get('border') is not None}
        fills.update({id(bar): layer for layer in layers for bar in layer.get('bars', [])})

        def _unfilled(border: dict) -> str:
            # The rings of a layer in the fill mode 'none' are all that the layer draws, so the
            # warnings about rings that no map shows end by saying so.
            fill = fills[id(border)]
            if fill['unified_mode'] != 'none':
                return ""
            element_type = fill['element_type']
            return (
                f" '--{element_type}-color none' leaves the fills of the {element_type} layer "
                f"uncolored, so the layer shows nothing where it shows no border."
            )

        if norm_layers or ring_layers:
            self.progress.new("Computing the range of values across maps")
            for layer in norm_layers:
                layer['_unified_vals'] = []
                layer['_individual_vals'] = []
                # The drawing passes below share this cache. They then have nothing left to work
                # out. This pass asks for the normalized values of every element of every drawn map.
                # A normalized value covers all of the samples or groups at once.
                layer['_normalized_cache'] = {}
                # The same holds for each element's values, and for the elements whose summary is
                # undefined ('_summarize_entry_values').
                layer['_element_cache'] = {}
                layer['_undefined_unified'] = set()
                layer['_undefined_individual'] = set()
                # The 'unified' map and the maps of the groups show summaries. A summary can be
                # undefined for every element of a KEGG map. An example is a standard deviation
                # where each element has a value in a single sample. The layer then has nothing to
                # color there, so the maps where it has a summarized value are recorded.
                layer['_unified_maps'] = set()
                layer['_individual_maps'] = {
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
                    per_individual = (
                        has_individuals and layer['individual_mode'] == 'quantitative'
                        and layer['sample_values'] is not None
                    )
                    normalize = layer.get('element_normalize')
                    # A sliced layer finds the range of each scale from the slices it draws, and a
                    # border with sliced rings from the slices of its rings.
                    if layer.get('slice_by') is not None:
                        find_range = (
                            self._find_ring_range if _is_border(layer) else self._find_slice_range
                        )
                        find_range(
                            pathway_number, pathway, fills.get(id(layer), layer), unified_valued,
                            individuals if per_individual and draw_individual_maps else []
                        )
                        continue
                    for entry in self._find_element_entries(
                        pathway, use_reaction, layer['accessions']
                    ):
                        if _role(layer) != 'fill' and not self._shows_ring(
                            pathway_number, pathway, layer, entry
                        ):
                            continue
                        individual_values, unified_value = self._summarize_entry_values(
                            entry, layer
                        )
                        if unified_value is not None:
                            layer['_unified_maps'].add(pathway_number)
                            if unified_valued:
                                layer['_unified_vals'].append(unified_value)
                        for group, group_maps in layer['_individual_maps'].items():
                            if group in individual_values:
                                group_maps.add(pathway_number)
                        if per_individual and normalize is not None:
                            # A normalization sets the value in each sample or group from the values
                            # in all of them at once, so the scale must span the calculated
                            # normalized values.
                            layer['_individual_vals'].extend(
                                self._normalize_entry_value(
                                    entry, layer, individuals, normalize,
                                    cache=layer['_normalized_cache']
                                ).values()
                            )
                        elif per_individual:
                            layer['_individual_vals'].extend(
                                individual_values[individual] for individual in individuals
                                if individual in individual_values
                            )
                for layer in ring_layers:
                    unified_shown = layer['unified_mode'] == 'categorical' and (
                        draw_unified_maps
                        or (layer['sample_values'] is None and draw_individual_maps)
                    )
                    individual_shown = (
                        draw_individual_maps and layer['individual_mode'] == 'categorical'
                    )
                    if layer.get('slice_by') is not None:
                        self._find_ring_maps(
                            pathway_number, pathway, fills[id(layer)], unified_shown,
                            draw_individuals if individual_shown else []
                        )
                        continue
                    for entry in self._find_element_entries(
                        pathway, layer['use_reaction_attribute'], layer['accessions']
                    ):
                        if not self._shows_ring(pathway_number, pathway, layer, entry):
                            continue
                        individual_values, unified_value = self._summarize_entry_individuals(
                            entry, layer, count_choices=False
                        )
                        if unified_shown and unified_value is not None:
                            layer['_ring_maps']['unified'].add(pathway_number)
                        # A border without samples shows its one value on every map.
                        if individual_shown and (
                            unified_value is not None if layer['sample_values'] is None
                            else any(name in individual_values for name in draw_individuals)
                        ):
                            layer['_ring_maps']['individual'].add(pathway_number)
            self.progress.end()
            unaffected_clause = (
                " The 'unified' map is unaffected, being drawn from the unnormalized values."
            ) if draw_unified_maps else ""
            for layer in norm_layers:
                if layer['sample_values'] is not None:
                    self._warn_undefined_summaries(layer, draw_unified_maps, draw_individual_maps)
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
                if layer.get('element_normalize') is not None and not layer['_individual_vals']:
                    cause = (
                        "An offset from a circular mean is undefined where the values cancel out, "
                        "such as 6 and 18 with a period of 24." if offsets else
                        f"A ratio is undefined wherever the value it is measured against is not a "
                        f"positive number, and a z-score wherever an element is found in a single "
                        f"{individual_noun}."
                    )
                    if _is_border(layer):
                        self.run.warning(
                            f"'--{_stem(layer)}-element-normalization' gave no value to any border "
                            f"on the maps of the individual {individual_noun}s, so no scale was "
                            f"drawn for them. These maps show no border from the "
                            f"'{layer['value_column']}' column of the {layer['element_type']} "
                            f"layer. Either none of the drawn maps shows a border for its "
                            f"accessions, or the normalization is undefined for every map element. "
                            f"{cause}{unaffected_clause}{_unfilled(layer)}"
                        )
                    else:
                        self.run.warning(
                            f"Nothing on the maps of the individual {individual_noun}s could be "
                            f"colored by '--{layer['element_type']}-element-normalization', so "
                            f"those maps carry no colors from the '{layer['value_column']}' column "
                            f"of the {layer['element_type']} layer and no scale was drawn for "
                            f"them. Either none of the drawn maps contains its accessions, or the "
                            f"normalization is undefined for every map element. "
                            f"{cause}{unaffected_clause}"
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
                        consequence = (
                            'show no border' if _is_border(layer) else 'are left uncolored'
                        )
                        self.run.warning(
                            f"'--{_stem(layer)}-element-normalization' was undefined for "
                            f"{len(undefined)} map element(s). Elements with the same accessions "
                            f"count once. Their accessions in the file include these: {examples}. "
                            f"{self._undefined_cause('circular_mean')} These elements have no "
                            f"circular mean to be offset from. They {consequence} on the maps of "
                            f"the individual {individual_noun}s."
                        )
                elif not layer['_unified_vals'] and not layer['_individual_vals'] and (
                    _role(layer) == 'bar'
                ):
                    self.run.warning(
                        f"No map shows a bar from the '{layer['value_column']}' bar column of the "
                        f"{layer['element_type']} layer, so no bar key was drawn for it. Either "
                        f"none of the drawn maps shows a bar for its accessions, or the "
                        f"aggregation reducing them is undefined for every map element, as the "
                        f"standard deviation of a single value is."
                    )
                elif not layer['_unified_vals'] and not layer['_individual_vals'] and (
                    _is_border(layer)
                ):
                    self.run.warning(
                        f"No map shows a border from the '{layer['value_column']}' border column "
                        f"of the {layer['element_type']} layer, so no color scale was drawn for "
                        f"it. Either none of the drawn maps shows a border for its accessions, or "
                        f"the aggregation reducing them is undefined for every map element, as the "
                        f"standard deviation of a single value is.{_unfilled(layer)}"
                    )
                elif not layer['_unified_vals'] and not layer['_individual_vals']:
                    self.run.warning(
                        f"Nothing on the maps could be colored by the values of the "
                        f"'{layer['value_column']}' column of the {layer['element_type']} layer, "
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
                individual_value_center = layer.get('individual_value_center')
                if period is None:
                    individual_value_limits = layer.get('individual_value_limits')
                elif layer.get('element_normalize') is not None:
                    individual_value_limits = (-period / 2, period / 2)
                    individual_value_center = 0.0
                else:
                    individual_value_limits = (0.0, period)
                # A bar scale always includes 0 ('TXT_CHANNEL_ROLES').
                baseline = layer['channel']['baseline'] if 'channel' in layer else None
                if layer['unified_mode'] == 'quantitative' and _unified_scale_drawn(layer):
                    norm, vmin, vmax, clamped_low, clamped_high = self._make_quantitative_norm(
                        layer['_unified_vals'], value_limits, f"--{_stem(layer)}-value-limits",
                        center=layer.get('value_center'),
                        center_flag=f"--{_stem(layer)}-value-center",
                        subject=(
                            f"the {_noun(layer)} layer of the 'unified' map"
                            if has_individuals else f"the {_noun(layer)} layer"
                        ),
                        baseline=baseline
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
                if layer['individual_mode'] != 'quantitative':
                    continue
                if has_individuals and layer['sample_values'] is not None:
                    norm, vmin, vmax, clamped_low, clamped_high = self._make_quantitative_norm(
                        layer['_individual_vals'], individual_value_limits,
                        f"--{_stem(layer)}-individual-value-limits",
                        center=individual_value_center,
                        center_flag=layer.get(
                            'individual_value_center_flag',
                            f"--{_stem(layer)}-individual-value-center"
                        ),
                        subject=(
                            f"the {_noun(layer)} layer of the maps of the individual "
                            f"{individual_noun}s"
                        ),
                        baseline=baseline
                    )
                    layer['_individual_norm'] = norm
                    layer['_individual_range'] = (vmin, vmax)
                    layer['_individual_limited'] = tuple(
                        limit is not None for limit in (individual_value_limits or (None, None))
                    )
                    layer['_individual_clamped'] = (clamped_low, clamped_high)
                    layer['_individual_center'] = individual_value_center
                else:
                    # A layer without samples of its own is the same on every individual map.
                    layer['_individual_norm'] = layer['_unified_norm']
                    layer['_individual_range'] = layer['_unified_range']
                    layer['_individual_limited'] = layer['_unified_limited']
                    layer['_individual_clamped'] = layer['_unified_clamped']
                    layer['_individual_center'] = layer['_unified_center']
            # A categorical border whose rings no map shows gets no legend. The data left every ring
            # empty, since the options alone were checked before ('_check_channel_maps').
            for layer in ring_layers:
                if not (layer['_ring_maps']['unified'] or layer['_ring_maps']['individual']):
                    self.run.warning(
                        f"No map shows a border from the '{layer['value_column']}' border column "
                        f"of the {layer['element_type']} layer, so no legend was drawn for it. "
                        f"Either none of the drawn maps shows a border for its accessions, or the "
                        f"rules that choose among its values gave no value.{_unfilled(layer)}"
                    )
            # Lines show no ring and no bar, so a reaction element that a standard map draws as a
            # line shows only its fill. Users may expect a ring or a bar there. A layer in the fill
            # mode 'none' shows nothing there.
            for layer in borders + bars:
                lines = layer['_line_entries']
                if not lines:
                    continue
                count = sum(len(entries) for entries in lines.values())
                shown = (
                    "no border. '--reaction-color none' leaves their fills uncolored, so these "
                    "elements show nothing" if fills[id(layer)]['unified_mode'] == 'none'
                    else f"no {'border' if _is_border(layer) else 'bars'}, so these elements show "
                    f"only their fills"
                )
                self.run.warning(
                    f"On {len(lines)} standard map(s), {count} reaction element(s) with "
                    f"{layer['channel']['values']} in the column '{layer['value_column']}' are "
                    f"drawn as lines. Lines show {shown}. The maps are these: "
                    f"{', '.join(sorted(lines)[:10])}{', and more' if len(lines) > 10 else ''}."
                )

        # A bar scale of the 'fill' style with values on both sides of 0 drops the labels of the
        # boxes on its maps ('_stage_entry_bar'). The key table of unlabeled boxes is written for
        # such bars ('_write_label_keys'), and whether it can be written is checked before anything
        # is drawn.
        for bar in bars:
            bar_ranges = [bar.get('_unified_range'), bar.get('_individual_range')]
            bar['_signed'] = bar['style'] == 'fill' and any(
                bar_range[0] < 0 < bar_range[1] for bar_range in bar_ranges
                if bar_range is not None and bar_range[0] is not None
            )
            if bar['_signed']:
                filesnpaths.is_output_file_writable(
                    os.path.join(output_dir, f"label_key_{bar['name']}.txt"),
                    ok_if_exists=self.overwrite_output
                )

        # A context colored by membership needs its by-count/by-membership color scheme over the
        # samples, contigs databases, genomes, or groups. Only the 'unified' map draws such a scale,
        # and its colorbar, from the samples, contigs databases, genomes, or groups themselves. The
        # presence coloring of grouped individual maps counts each group's own sources instead. It
        # is precomputed below. A run leaving out the 'unified' map works out no such scale. The
        # limits that keep such a scale legible do not apply to it.
        membership_layers = [
            layer for layer in layers
            if draw_unified_maps and layer['unified_mode'] == 'membership'
        ]
        group_membership_layers = [
            layer for layer in layers if layer['individual_mode'] == 'membership'
        ]
        # Where a count scale is asked to stop at the highest count in the data, that count has to
        # be found first, which takes one parse of the maps about to be drawn.
        group_observed_counts: Dict[str, int] = {}
        if count_scale_max == 'observed':
            scaled_group_layers = group_membership_layers if (
                draw_individual_maps and grouped_presence
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
                    count_scale_max, layer.get('_observed_count', 0), len(individuals)
                )
                layer['_colors'] = self._membership_layer_colors(
                    layer, individuals, layer['_count_scale_top']
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
                    grouped_membership, layers, draw_individuals, group_observed_counts,
                    count_scale_max, unified_plural, individual_noun
                )
            )

        def _reaction_derived(layer, mode, cmap_key='cmap', spread=False):
            # How a reaction layer derives compound colors on a global or overview map where no
            # compound layer stages fills ('_derives_compounds'). The derived colors come from the
            # same scale as the reactions they are derived from, so the caller names the context's
            # colormap: the two can differ ('individual_cmap'). The caller also says whether the
            # context shows a spread, which is a plain number and so is averaged as one.
            if layer['element_type'] != 'reaction':
                return None
            # Categorical values are chosen, not averaged, by the layer's accession rule. The
            # chooser is built once per layer.
            if mode == 'categorical':
                if '_compound_chooser' not in layer:
                    layer['_compound_chooser'] = self._categorical_compound_chooser(layer)
                return (layer['_compound_chooser'], None)
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

        def _unfilled_spec(layer):
            # A layer in the fill mode 'none' colors no fill, so its spec has no colorer. Its
            # entries do not count as data, and it stages only the rings of its border
            # ('_with_border', '_stage_map_elements').
            return {
                'element_type': layer['element_type'],
                'use_reaction_attribute': layer['use_reaction_attribute'],
                'entry_keys': layer['accessions'],
                'colorer': None,
                'derived_compound': _reaction_derived(layer, 'none')
            }

        def _unified_spec(layer):
            mode = layer['unified_mode']
            if mode == 'none':
                return _unfilled_spec(layer)
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
            if mode == 'categorical':
                def entry_value(entry):
                    return self._summarize_entry_individuals(entry, layer)[1]
                entry_keys = layer['unified_values'] if layer['sample_values'] is None else (
                    _sample_accessions(layer, layer['sample_values'])
                )
                return {
                    'element_type': layer['element_type'],
                    'use_reaction_attribute': layer['use_reaction_attribute'],
                    'entry_keys': entry_keys,
                    'colorer': self._categorical_colorer(
                        entry_value, layer['value_color_priorities']
                    ),
                    # A rule can leave every element of a map without a value.
                    'require_color': True,
                    'derived_compound': _reaction_derived(layer, mode)
                }
            if mode == 'membership':
                _, color_priorities, membership_combos, _ = layer['_colors']
                return {
                    'element_type': layer['element_type'],
                    'use_reaction_attribute': layer['use_reaction_attribute'],
                    'entry_keys': layer['membership'],
                    'colorer': self._membership_colorer(
                        layer['membership'], color_priorities, membership_combos, group_sources,
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

        def _individual_spec(layer, individual):
            mode = layer['individual_mode']
            if mode == 'none':
                return _unfilled_spec(layer)
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
                    # normalization can still give a sample or group no value for an element. This
                    # happens where it had nothing to rescale, or where the rescaling itself was
                    # undefined.
                    if layer['group_samples'] is None:
                        entry_keys = layer['sample_values'][individual]
                    else:
                        entry_keys = _sample_accessions(layer, layer['group_samples'][individual])
                        maps = layer['_individual_maps'][individual]
                    if normalize is None:
                        def entry_value(entry):
                            return self._summarize_entry_values(entry, layer)[0].get(individual)
                    else:
                        def entry_value(entry):
                            return self._normalize_entry_value(
                                entry, layer, individuals, normalize,
                                cache=layer['_normalized_cache']
                            ).get(individual)
                return {
                    'element_type': layer['element_type'],
                    'use_reaction_attribute': layer['use_reaction_attribute'],
                    'entry_keys': entry_keys,
                    'pathway_numbers': maps,
                    'colorer': self._quantitative_colorer(
                        entry_value, layer['_individual_norm'], layer['individual_cmap'],
                        layer['reverse_overlay'], center=layer['_individual_center']
                    ),
                    'derived_compound': _reaction_derived(layer, mode, 'individual_cmap')
                }
            if mode == 'categorical':
                # A layer without samples draws the same values on every map. A sample's map finds
                # the elements of the sample's accessions. A group's map finds those of its samples'
                # accessions.
                if layer['sample_values'] is None:
                    entry_keys = layer['unified_values']

                    def entry_value(entry):
                        return self._summarize_entry_individuals(entry, layer)[1]
                else:
                    entry_keys = layer['sample_values'][individual] if (
                        layer['group_samples'] is None
                    ) else _sample_accessions(layer, layer['group_samples'][individual])

                    def entry_value(entry):
                        values = self._summarize_entry_individuals(entry, layer, individual)[0]
                        return values.get(individual)
                return {
                    'element_type': layer['element_type'],
                    'use_reaction_attribute': layer['use_reaction_attribute'],
                    'entry_keys': entry_keys,
                    'colorer': self._categorical_colorer(
                        entry_value, layer['value_color_priorities']
                    ),
                    # A rule can leave every element of a map without a value.
                    'require_color': True,
                    'derived_compound': _reaction_derived(layer, mode)
                }
            if mode == 'single':
                return {
                    'element_type': layer['element_type'],
                    'use_reaction_attribute': layer['use_reaction_attribute'],
                    'entry_keys': layer['accessions'],
                    'colorer': self._single_color_colorer(layer['color_hexcode']),
                    'derived_compound': _reaction_derived(layer, mode, 'individual_cmap')
                }
            # Ungrouped membership/static: an individual source's map colors that source's elements
            # in one fixed color (the grouped case is precomputed below). With a membership colors
            # file, that color is the source's own color from the file. A panel of a grid then shows
            # by its color which source it is. The color also matches the source's band on the
            # 'unified' map's membership scale.
            accessions = layer['source_accessions'].get(individual, set())
            membership_colors = layer.get('membership_colors')
            color_hexcode = layer['color_hexcode'] if membership_colors is None else (
                membership_colors[individual]
            )
            return {
                'element_type': layer['element_type'],
                'use_reaction_attribute': layer['use_reaction_attribute'],
                'entry_keys': accessions,
                'colorer': self._single_color_colorer(color_hexcode),
                'derived_compound': _reaction_derived(layer, mode, 'individual_cmap')
            }

        # A sliced layer slices every drawn map that shows its values. Its spec gets a slicer, which
        # stages its entries in place of the colorer ('_stage_entry_slices'). The key table lists
        # the 'unified' map first, then the maps of the samples or groups in the order they are
        # drawn.
        slice_ranks = {individual: rank for rank, individual in enumerate(draw_individuals, 1)}

        def _gene_entry_keys(layer, individual):
            # A map sliced by gene finds the entries of the accessions that have rows in its
            # samples. The values of the accessions are not used for this. They leave out an
            # accession whose genes reduce to an undefined value, such as values that cancel out on
            # a circle. The samples, and the accessions of each sample, are taken in turn, as in
            # '_sample_accessions'.
            found = layer['slice_accessions']
            if layer['sample_values'] is None:
                return found
            group_samples = layer['group_samples']
            if individual is not None:
                samples = [individual] if group_samples is None else group_samples[individual]
            elif group_samples is None:
                samples = list(found)
            else:
                samples = [sample for members in group_samples.values() for sample in members]
            return dict.fromkeys(accession for sample in samples for accession in found[sample])

        def _sliced(spec, layer, individual=None):
            mode = layer['unified_mode'] if individual is None else layer['individual_mode']
            if layer.get('slice_by') is None or mode not in ('quantitative', 'categorical'):
                return spec
            if layer['slice_by'] == 'gene':
                spec = {**spec, 'entry_keys': _gene_entry_keys(layer, individual)}
            return {
                **spec,
                'slicer': self._slicer(layer, individual),
                'slice_map': 'unified' if individual is None else individual,
                'slice_rank': 0 if individual is None else slice_ranks[individual],
                'slice_key': layer['_slice_key'],
                'slice_uncolored': layer['_slice_uncolored'],
                # In a layer with a border column, the columns of a box fit inside its ring, so the
                # cap measures the width inside the ring ('_slice_cap').
                'slice_ringed': layer.get('border') is not None
            }

        # A layer with a border adds the border's colorer to the spec of each map that shows the
        # border. The border is staged after the fill, on the shapes that show rings
        # ('_stage_entry_borders'). Its colorer comes from the spec the border model would have as a
        # layer of its own. A border with sliced rings adds a ringer, which finds the slices of each
        # ring ('_ringer'), and the keys of its own key table. Rings sliced by gene find their
        # entries from the accessions with border rows, as the fill does ('_gene_entry_keys').
        def _with_border(spec, layer, individual=None):
            border = layer.get('border')
            if border is None:
                return spec
            mode = border['unified_mode'] if individual is None else border['individual_mode']
            if mode not in ('quantitative', 'categorical'):
                return spec
            border_spec = (
                _unified_spec(border) if individual is None
                else _individual_spec(border, individual)
            )
            spec = {
                **spec,
                'border_colorer': border_spec['colorer'],
                'border_entry_keys': border_spec['entry_keys'],
                'border_pathway_numbers': border_spec.get('pathway_numbers')
            }
            if border.get('slice_by') is None:
                return spec
            if border['slice_by'] == 'gene':
                spec['border_entry_keys'] = _gene_entry_keys(border, individual)
            return {
                **spec,
                'border_ringer': self._ringer(layer, individual),
                'border_slice_map': 'unified' if individual is None else individual,
                'border_slice_rank': 0 if individual is None else slice_ranks[individual],
                'border_slice_key': border['_slice_key'],
                'border_slice_uncolored': border['_slice_uncolored'],
                'border_slice_unshown': border['_slice_unshown']
            }

        # A layer with bars adds a function to the spec of each map that shows them. It gives the
        # bar value of an Entry on the map, or None where the Entry has none there, and the spec
        # gets the bar scale of the map too. The boxes of such a map are drawn as bars in the colors
        # of their fills, or with strips under whole fills, by the style of the bars
        # ('_stage_entry_bar'). A map that shows no bar values, such as the 'unified' map without a
        # bar summary, draws its fills whole. The bars of a layer are a list of one model.
        def _with_bars(spec, layer, individual=None):
            if not layer.get('bars'):
                return spec
            bar = layer['bars'][0]
            mode = bar['unified_mode'] if individual is None else bar['individual_mode']
            if mode != 'quantitative':
                return spec
            norm = bar['_unified_norm'] if individual is None else bar['_individual_norm']
            # A layer without samples draws the same bars on every map.
            if individual is None or bar['sample_values'] is None:
                def entry_value(entry):
                    return self._summarize_entry_values(entry, bar)[1]
            else:
                def entry_value(entry):
                    return self._summarize_entry_values(entry, bar)[0].get(individual)
            return {
                **spec,
                'bar_value': entry_value,
                'bar_norm': norm,
                'bar_style': bar['style'],
                'bar_map': 'unified' if individual is None else individual,
                'bar_rank': 0 if individual is None else slice_ranks[individual],
                'bar_missing': bar['_missing'],
                'bar_unlabeled': bar['_unlabeled']
            }

        # The spec of a layer on the 'unified' map, or on the map of an individual, colors its fills,
        # slices them, and then adds its rings and bars.
        def _layer_spec(layer, individual=None):
            if individual is None:
                spec = _unified_spec(layer)
            else:
                spec = _individual_spec(layer, individual)
            spec = _sliced(spec, layer, individual)
            spec = _with_border(spec, layer, individual)
            return _with_bars(spec, layer, individual)

        # For grouped membership individual maps, precompute each group's within-group membership,
        # narrowing the layer's membership to the group's own sources so that a group's map counts
        # only them. The colors these counts take were resolved above, with every other color of the
        # run.
        group_layer_membership: Dict[str, Tuple] = {}
        if grouped_presence:
            for group in draw_individuals:
                specs = []
                for layer in group_membership_layers:
                    inner_membership = {}
                    for accession, sources in layer['membership'].items():
                        in_group = [s for s in sources if source_group.get(s) == group]
                        if in_group:
                            inner_membership[accession] = in_group
                    spec = {
                        'element_type': layer['element_type'],
                        'use_reaction_attribute': layer['use_reaction_attribute'],
                        'entry_keys': inner_membership,
                        'colorer': self._membership_colorer(
                            inner_membership, group_color_priorities[group], None, None, None,
                            layer['use_reaction_attribute']
                        ),
                        'keep_highest_priority': True,
                        'derived_compound': _reaction_derived(layer, 'membership')
                    }
                    spec = _with_border(spec, layer, group)
                    specs.append(_with_bars(spec, layer, group))
                group_layer_membership[group] = (
                    specs, group_color_priorities[group], group_scale_tops[group]
                )

        def _individual_specs(individual):
            # Grouped membership specs are precomputed; a layer colored by value or a single color
            # in its per-group context colors its own map and rides along on the same map.
            if grouped_presence:
                return group_layer_membership[individual][0] + [
                    _layer_spec(layer, individual)
                    for layer in layers if layer['individual_mode'] != 'membership'
                ]
            return [_layer_spec(layer, individual) for layer in layers]

        def _element_map_specs():
            # The specs of each kind of map that '_draw_map_elements' draws, with a phrase naming
            # it. A run in the reference map's own colors draws its other maps with
            # '_draw_map_kos_original_color'.
            if draw_unified_maps and not original_run:
                yield (
                    "the 'unified' map" if has_individuals else None, True,
                    [_layer_spec(layer) for layer in layers]
                )
            for individual in draw_individuals:
                if grouped_presence or not original_run:
                    yield f"{individual_noun} '{individual}'", False, _individual_specs(individual)

        self._check_reserved_colors(
            layers, pathway_numbers, group_color_priorities=group_color_priorities,
            group_element_types={layer['element_type'] for layer in group_membership_layers},
            drawn_individuals=draw_individuals, draw_unified_maps=draw_unified_maps
        )
        self._check_derived_compound_colors(
            layers, pathway_numbers, _element_map_specs,
            group_color_priorities=group_color_priorities
        )

        # The values of a layer can color the maps of the individual samples or groups alone, while
        # the 'unified' map shows presence across them. A run that draws no individual map then
        # colors nothing by the values. This is not an error. Presence across samples is worth a
        # map, and a value column is read either way. But the values go unused, so the user is told
        # how to show them. The warning comes after every check that can refuse the run, so it is
        # given only where maps are drawn.
        value_modes = ('quantitative', 'categorical')
        for layer in layers:
            element_type = layer['element_type']
            # A layer in the fill mode 'none' draws only rings. Where its border shows no ring on
            # the 'unified' map, that map shows nothing of the layer. A run in which no map can show
            # the rings was refused ('_check_channel_maps').
            if layer['unified_mode'] == 'none':
                border = layer['border']
                if draw_unified_maps and border['unified_mode'] is None:
                    grouped_layer = border.get('group_samples') is not None
                    summary_flag = (
                        f"--{element_type}-border-{'group' if grouped_layer else 'sample'}-summary"
                    )
                    self.run.warning(
                        f"The 'unified' map shows nothing of the {element_type} layer. "
                        f"'--{element_type}-color none' leaves its fills uncolored. Its border "
                        f"shows no ring there, since '{summary_flag}' was not given. Give that "
                        f"option an aggregation or a rule to show rings on the 'unified' map, or "
                        f"leave that map out with '--skip-unified-maps'."
                    )
                continue
            # A file of categorical values with a header and no rows has no values to show. The
            # other warnings here tell the user how to show values, so they do not apply.
            if layer.get('value_type') == 'categorical' and not layer['value_order']:
                self.run.warning(
                    f"The text file at '{layer['path']}' has a header but no rows. The "
                    f"{element_type} layer has no categorical values to color the maps by, so no "
                    f"legend of values was drawn."
                )
                continue
            # Grouped categorical values reach a map only through a rule among each group's samples
            # ('_build_categorical_txt_model'). Without one, every map counts samples or groups.
            if (
                layer.get('value_type') == 'categorical'
                and layer['individual_mode'] == 'membership'
            ):
                draw_clause = "" if draw_individual_maps else (
                    " Then draw the maps of the groups with '--draw-individual-files' or "
                    "'--draw-grid'."
                )
                unified_clause = (
                    " The 'unified' map shows presence across the groups." if draw_unified_maps
                    else ""
                )
                self.run.warning(
                    f"The {element_type} layer from the text file at '{layer['path']}' has "
                    f"categorical values in its '{layer['value_column']}' column, but nothing on "
                    f"the maps is colored by them. With groups of samples, each group's map counts "
                    f"the group's samples.{unified_clause} Give '--{element_type}-sample-summary' "
                    f"a rule, such as 'plurality_order', to color each group's map by a value "
                    f"chosen among its samples.{draw_clause}"
                )
                continue
            if draw_individual_maps or layer['individual_mode'] not in value_modes:
                continue
            if layer['unified_mode'] in value_modes:
                continue
            grouped_layer = layer.get('group_samples') is not None
            noun = 'groups' if grouped_layer else 'samples'
            summary_flag = f"--{element_type}-{'group' if grouped_layer else 'sample'}-summary"
            remedy = (
                f"set '{summary_flag}' to a categorical rule, such as 'plurality_order', to color "
                f"the 'unified' map by value. A rule assumes that a value means the same thing in "
                f"every {noun[:-1]}."
            ) if layer.get('value_type') == 'categorical' else (
                f"set '{summary_flag}' to an aggregation, such as "
                f"'{'mean' if layer.get('value_period') is None else 'circular_mean'}', to color "
                f"the 'unified' map by value."
            )
            self.run.warning(
                f"The {element_type} layer from the text file at '{layer['path']}' has a value "
                f"column, '{layer['value_column']}', but nothing on the maps is colored by it. The "
                f"'unified' map shows presence across the {noun}. The values color only the maps "
                f"of the individual {noun}, and none are drawn here. Draw them with "
                f"'--draw-individual-files' or '--draw-grid'. Or {remedy}"
            )
        self._warn_gene_ids(models, draw_unified_maps, draw_individuals)

        # Draw per-layer colorbars for the unified map. Their file names include the layer's name,
        # so the colorbars of two layers do not collide. Each layer colored quantitatively on the
        # individual maps also gets one colorbar that those maps share. Each context is keyed by its
        # own mode, so a layer summarized by presence in the unified map gets a presence colorbar
        # there. It also gets a continuous one for its individual maps. A run leaving out the
        # unified map draws neither of its bars. The exception is a layer without samples of its
        # own. Its individual maps take the unified scale, so that colorbar is drawn as the only key
        # to their colors. A border is keyed in files of its own, such as
        # 'colorbar_reactions_border.pdf'. Its labels end in ' (border)', except a label that was
        # given with its normalization. The colorbars that key the individual maps are listed, so a
        # player can show them beside those maps. A border's colorbar is listed with the element
        # type of its border, since a player shows it only beside the maps that can show the
        # border's rings ('_map_shows_borders'). Bars are keyed by length in files of their own,
        # such as 'bar_key_reactions.pdf' ('ColorbarDrawer.draw_length'). A bar key is labeled with
        # the bar column and '(bar)'. Bar keys are listed for players as the keys of borders are,
        # since bars take the shapes that rings take.
        individual_keys: List[Tuple[str, Union[str, None]]] = []
        for layer in models:
            if _role(layer) == 'bar':
                contexts = []
                if layer['unified_mode'] == 'quantitative' and _unified_scale_drawn(layer):
                    # A layer without samples of its own keys its individual maps with this key.
                    contexts.append((
                        'unified', f"bar_key_{layer['name']}.pdf",
                        draw_individual_maps and layer['individual_mode'] == 'quantitative'
                        and layer['sample_values'] is None
                    ))
                if (
                    draw_individual_maps and layer['individual_mode'] == 'quantitative'
                    and layer['sample_values'] is not None
                ):
                    contexts.append((
                        'individual', f"bar_key_{layer['name']}_{colorbar_individual_suffix}.pdf",
                        True
                    ))
                for context, basename, keys_individual_maps in contexts:
                    vmin, vmax = layer[f'_{context}_range']
                    if vmin is None:
                        continue
                    clamped_low, clamped_high = layer[f'_{context}_clamped']
                    key_path = os.path.join(output_dir, basename)
                    self.colorbar_drawer.draw_length(
                        vmin, vmax, key_path, label=f"{layer['colorbar_label']} (bar)",
                        clamped_low=clamped_low, clamped_high=clamped_high,
                        track_color=(
                            BAR_STRIP_TRACK_COLOR if layer['style'] == 'strip' else BAR_TRACK_COLOR
                        )
                    )
                    if keys_individual_maps:
                        individual_keys.append((key_path, layer['element_type']))
                continue
            suffix = ' (border)' if _is_border(layer) else ''
            border_element_type = layer['element_type'] if _is_border(layer) else None
            if layer['unified_mode'] == 'quantitative' and _unified_scale_drawn(layer):
                vmin, vmax = layer['_unified_range']
                if vmin is not None:
                    limited_low, limited_high = layer['_unified_limited']
                    clamped_low, clamped_high = layer['_unified_clamped']
                    colorbar_path = os.path.join(output_dir, f"colorbar_{layer['name']}.pdf")
                    self._draw_quantitative_colorbar(
                        layer['cmap'], vmin, vmax, colorbar_path,
                        layer.get('unified_colorbar_label', layer['colorbar_label']) + suffix,
                        limited_low=limited_low, limited_high=limited_high,
                        clamped_low=clamped_low, clamped_high=clamped_high,
                        center=layer['_unified_center']
                    )
                    if (
                        draw_individual_maps and layer['individual_mode'] == 'quantitative'
                        and layer['sample_values'] is None
                    ):
                        individual_keys.append((colorbar_path, border_element_type))
            elif layer['unified_mode'] == 'membership' and draw_unified_maps:
                scheme, color_priorities, membership_combos, presence_cmap = layer['_colors']
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
                        labels = [', '.join(combo) for combo in membership_combos]
                        label = 'groups' if grouped else membership_members_label
                    self.colorbar_drawer.draw_discrete(
                        [color for color, _ in color_priorities],
                        colorbar_path,
                        color_labels=labels,
                        label=label
                    )
            # A categorical layer has no legend where no map shows its values, or where its file has
            # a header and no rows. It has no colors in the first case and no values in the second.
            # A border has no legend where no map shows its rings.
            if layer.get('value_type') == 'categorical' and layer['value_color_priorities']:
                ring_maps = layer.get('_ring_maps', {'unified': True, 'individual': True})
                # One band per value, in the order of the values, with the first value at the top.
                # This order does not depend on 'reverse_overlay'. A colorbar draws its first band
                # at the bottom, so it is given the values in reverse. A layer without samples draws
                # its values on the individual maps too, so its legend is drawn where only those
                # maps are.
                legend_colors = [
                    color for color, _ in reversed(layer['value_color_priorities'].values())
                ]
                legend_values = list(reversed(layer['value_color_priorities']))
                legend_paths = []
                if layer['unified_mode'] == 'categorical' and (
                    draw_unified_maps or (layer['sample_values'] is None and draw_individual_maps)
                ) and ring_maps['unified']:
                    legend_paths.append(os.path.join(output_dir, f"colorbar_{layer['name']}.pdf"))
                    if layer['sample_values'] is None and draw_individual_maps:
                        individual_keys.append((legend_paths[-1], border_element_type))
                if (
                    draw_individual_maps and layer['sample_values'] is not None
                    and layer['individual_mode'] == 'categorical' and ring_maps['individual']
                ):
                    legend_paths.append(os.path.join(
                        output_dir, f"colorbar_{layer['name']}_{colorbar_individual_suffix}.pdf"
                    ))
                    individual_keys.append((legend_paths[-1], border_element_type))
                for legend_path in legend_paths:
                    self.colorbar_drawer.draw_discrete(
                        legend_colors, legend_path, color_labels=legend_values,
                        label=layer['colorbar_label'] + suffix
                    )
            if (
                draw_individual_maps and layer['individual_mode'] == 'quantitative'
                and layer['sample_values'] is not None
            ):
                vmin, vmax = layer['_individual_range']
                if vmin is not None:
                    limited_low, limited_high = layer['_individual_limited']
                    clamped_low, clamped_high = layer['_individual_clamped']
                    individual_label = layer.get(
                        'individual_colorbar_label', layer['colorbar_label']
                    )
                    if not layer.get('individual_colorbar_label_given'):
                        individual_label += suffix
                    colorbar_path = os.path.join(
                        output_dir, f"colorbar_{layer['name']}_{colorbar_individual_suffix}.pdf"
                    )
                    self._draw_quantitative_colorbar(
                        layer['individual_cmap'], vmin, vmax, colorbar_path, individual_label,
                        limited_low=limited_low, limited_high=limited_high,
                        clamped_low=clamped_low, clamped_high=clamped_high,
                        center=layer['_individual_center']
                    )
                    individual_keys.append((colorbar_path, border_element_type))

        # Draw the unified map. In a run with a single source of data, it is the only map.
        if draw_unified_maps:
            self.progress.new(
                f"Drawing 'unified' map incorporating data from all {unified_plural}"
                if has_individuals else "Drawing map"
            )
            unified_specs = None if original_run else [_layer_spec(layer) for layer in layers]
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

        if not draw_individuals:
            self._warn_missing_bars(bars)
            self._write_label_keys(bars, output_dir)
            self._write_slice_keys(models, output_dir)
            self._report_categorical_choices(
                models, draw_unified_maps, draw_individuals,
                draw_maps_lacking_data=draw_maps_lacking_data
            )
            count = sum(drawn['unified'].values()) if drawn['unified'] else 0
            self.run.info("Number of maps drawn", count)
            return drawn

        # The colorbar of each group's own count scale, where a group's maps have one. A player
        # shows it beside that group's maps.
        individual_colorbars: Dict[str, str] = {}
        for individual in draw_individuals:
            drawn_individual: Dict[str, bool] = {}
            self.progress.new(f"Drawing maps for {individual_noun} '{individual}'")
            self.progress.update("...")
            progress = self.progress
            self.progress = terminal.Progress(verbose=False)
            run = self.run
            self.run = terminal.Run(verbose=False)
            individual_output_dir = os.path.join(output_dir, INDIVIDUAL_SUBDIR, individual)
            filesnpaths.gen_output_directory(
                individual_output_dir, progress=self.progress, run=self.run
            )

            if grouped_presence:
                _, individual_color_priorities, individual_scale_top = (
                    group_layer_membership[individual]
                )
                colorbar_path = os.path.join(individual_output_dir, INDIVIDUAL_COLORBAR_BASENAME)
                if group_colormap_scheme == 'by_count_continuous':
                    # The same color per count that the discrete bands assign, shown as a gradient
                    # across the ramp from the lowest count to the highest rather than as one labeled
                    # band per count, which is what frees it from needing a distinct color for each.
                    self._draw_quantitative_colorbar(
                        group_cmaps[individual], 1, individual_scale_top, colorbar_path,
                        membership_count_label, integer_ticks=True
                    )
                else:
                    self.colorbar_drawer.draw_discrete(
                        [color for color, _ in individual_color_priorities],
                        colorbar_path,
                        color_labels=range(1, individual_scale_top + 1),
                        label=membership_count_label
                    )
                individual_colorbars[individual] = colorbar_path
                specs = _individual_specs(individual)
                for pathway_number in pathway_numbers:
                    progress.update(pathway_number)
                    drawn_individual[pathway_number] = self._draw_map_elements(
                        pathway_number, specs, individual_output_dir,
                        draw_map_lacking_data=draw_maps_lacking_data
                    )
            elif original_run:
                for pathway_number in pathway_numbers:
                    progress.update(pathway_number)
                    drawn_individual[pathway_number] = self._draw_map_kos_original_color(
                        pathway_number, layers[0]['source_accessions'].get(individual, set()),
                        individual_output_dir, draw_map_lacking_data=draw_maps_lacking_data
                    )
            else:
                specs = _individual_specs(individual)
                for pathway_number in pathway_numbers:
                    progress.update(pathway_number)
                    drawn_individual[pathway_number] = self._draw_map_elements(
                        pathway_number, specs, individual_output_dir,
                        draw_map_lacking_data=draw_maps_lacking_data
                    )

            self.progress = progress
            self.run = run
            self.progress.end()
            drawn['individual'][individual] = drawn_individual

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
                draw_individuals,
                draw_grid_individuals,
                draw_files_individuals,
                output_dir,
                drawn,
                group_scale_tops=grid_group_scale_tops,
                group_color_priorities=grid_group_color_priorities,
                group_cmaps=group_cmaps,
                group_colormap_scheme=group_colormap_scheme,
                check_maps_lacking_kos=not draw_maps_lacking_data,
                source_type=grid_source_type,
                include_unified=draw_unified_maps
            )

        # The individual maps are gathered only once nothing else will change them. Drawing the
        # grids deletes the maps that were drawn only as grid panels. It also deletes the
        # directories of the samples, contigs databases, genomes, or groups that were never asked
        # for individually.
        collated_count = 0
        animated_count = 0
        if self.collate_files_by_map and draw_files_individuals:
            self.progress.new(f"Gathering the maps of each {individual_noun}")
            self.progress.update("...")
            collated_count = self._collate_maps_by_map(output_dir, draw_files_individuals)
            self.progress.end()
            if self.animate_files_by_map:
                self.progress.new(f"Writing a player of each map across the {individual_noun}s")
                self.progress.update("...")
                animated_count = self._draw_map_players(
                    output_dir, pathway_numbers, draw_files_individuals, drawn['individual'],
                    individual_keys, individual_colorbars
                )
                self.progress.end()

        self._warn_missing_bars(bars)
        self._write_label_keys(bars, output_dir)
        self._write_slice_keys(models, output_dir)
        self._report_categorical_choices(
            models, draw_unified_maps, draw_individuals,
            draw_maps_lacking_data=draw_maps_lacking_data, draw_grid=draw_grid is not False
        )
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
            self.run.info(f"Number of maps drawn for individual {individual_noun}s", count)
            if self.collate_files_by_map:
                self.run.info(
                    f"Number of maps gathered across individual {individual_noun}s", collated_count
                )
                if self.animate_files_by_map:
                    self.run.info(
                        f"Number of maps with a player across individual {individual_noun}s",
                        animated_count
                    )
        count = sum(drawn['grid'].values()) if drawn['grid'] else 0
        self.run.info("Number of map grids drawn", count)

        return drawn

    @staticmethod
    def _shows_border(pathway: kgml.Pathway, entry: kgml.Entry, graphics: kgml.Graphics) -> bool:
        """
        Say whether a shape of a map element shows a border ring.

        A ring is drawn where a map draws an outline: around reaction boxes on standard maps, and
        around compound circles on standard and overview maps. Global maps draw circles without an
        outline, so they show no rings. No map draws a ring along a line. Both the range of a border
        scale and the staging of rings ask this, so a scale never spans a ring that is not drawn.

        Parameters
        ==========
        pathway : kgml.Pathway
            The map.

        entry : kgml.Entry
            An ortholog or compound Entry of the map.

        graphics : kgml.Graphics
            A Graphics element of the Entry.

        Returns
        =======
        bool
            True if the shape shows a ring.
        """
        if pathway.is_global_map:
            return False
        if entry.type == 'compound':
            return graphics.type == 'circle'
        if entry.type == 'ortholog':
            return graphics.type == 'rectangle' and not pathway.is_overview_map
        return False

    @staticmethod
    def _map_shows_borders(element_type: str, pathway_number: str) -> bool:
        """
        Say whether a map is of a kind that shows the rings of a border.

        A reaction border shows on standard maps, where reactions are boxes. A compound border shows
        on standard and overview maps, whose circles have an outline. Global maps show neither
        ('_shows_border').

        Parameters
        ==========
        element_type : str
            The element type of the border, 'reaction' or 'compound'.

        pathway_number : str
            Numeric ID of the map.

        Returns
        =======
        bool
            True if the map can show the rings of the border.
        """
        if re.match(GLOBAL_MAP_ID_PATTERN, pathway_number):
            return False
        return element_type == 'compound' or not re.match(OVERVIEW_MAP_ID_PATTERN, pathway_number)

    def _shows_ring(
        self,
        pathway_number: str,
        pathway: kgml.Pathway,
        border: dict,
        entry: kgml.Entry
    ) -> bool:
        """
        Say whether an Entry shows a ring on a map, and record the lines that show none.

        Users may expect a ring on a reaction element that a standard map draws as a line, so such
        an Entry is recorded in the border's '_line_entries', which a warning reads
        ('_map_elements'). Bars are drawn in the shapes that take rings, so the same holds for the
        model of bars.

        Parameters
        ==========
        pathway_number : str
            Numeric ID of the map.

        pathway : kgml.Pathway
            The map.

        border : dict
            A border or bar model with '_line_entries'.

        entry : kgml.Entry
            An ortholog or compound Entry of the map.

        Returns
        =======
        bool
            True if some shape of the Entry shows a ring ('_shows_border').
        """
        graphics = [pathway.uuid_element_lookup[uuid] for uuid in entry.children['graphics']]
        if any(self._shows_border(pathway, entry, shape) for shape in graphics):
            return True
        if not (pathway.is_global_map or pathway.is_overview_map) and any(
            shape.type == 'line' for shape in graphics
        ):
            border['_line_entries'].setdefault(pathway_number, set()).add(entry.id)
        return False

    def _check_channel_maps(
        self,
        channels: List[dict],
        pathway_numbers: List[str],
        draw_unified_maps: bool,
        draw_individual_maps: bool,
        individual_noun: Union[str, None],
        unfilled: Union[Set[str], None] = None
    ) -> None:
        """
        Refuse a border or bar that no map of a run can show, and count the maps that show none.

        A reaction border needs a standard map, where reactions are boxes. A compound border needs a
        standard or an overview map ('_map_shows_borders'). Bars are drawn in the boxes that take
        rings. A border or a bar also needs a drawn map that shows its values. One of samples shows
        none on the 'unified' map without a sample or group summary, and none on the maps of groups
        without a sample summary. These refusals come from the options alone. Where the data leave
        every ring or bar empty, the run warns instead.

        Parameters
        ==========
        channels : List[dict]
            The border and bar models of the run.

        pathway_numbers : List[str]
            The maps that the run draws.

        draw_unified_maps : bool
            True if the 'unified' map is drawn.

        draw_individual_maps : bool
            True if the maps of the individual samples or groups are drawn.

        individual_noun : Union[str, None]
            'sample' or 'group', or None for a run without samples.

        unfilled : Union[Set[str], None], None
            The element types of the layers in the fill mode 'none'. Their rings are all that they
            draw, so a refusal does not suggest dropping their border columns.
        """
        value_modes = ('quantitative', 'categorical')
        for model in channels:
            element_type = model['element_type']
            channel = model['channel']
            noun = channel['noun']
            border = channel['role'] == 'border'
            column_flag = f"--{channel['flag_stem']}-column"
            # What the channel draws, and what a map without its values shows of it.
            marks = 'borders' if border else 'bars'
            showing = [
                number for number in pathway_numbers
                if self._map_shows_borders(element_type, number)
            ]
            if not showing:
                listed = ', '.join(pathway_numbers[:10])
                more = (
                    f", and {len(pathway_numbers) - 10} more" if len(pathway_numbers) > 10 else ''
                )
                if element_type == 'reaction':
                    shapes = "reaction boxes, which only standard maps draw"
                    others = (
                        f"Global maps draw no {marks}, and overview maps draw reactions as lines, "
                        f"which have no {marks}."
                    )
                else:
                    shapes = "compound circles on standard and overview maps"
                    others = f"Global maps draw their circles without an outline, and no {marks}."
                remedy = (
                    f"Draw a map that shows borders. '--{element_type}-color none' leaves the "
                    f"fills of the {element_type} layer uncolored, so its rings are all that it "
                    f"draws." if border and element_type in (unfilled or set())
                    else f"Draw a map that shows {marks}, or drop '{column_flag}'."
                )
                draws = f"colors the borders of {shapes}" if border else f"draws bars in {shapes}"
                raise ConfigError(
                    f"'{column_flag}' {draws}. This run draws no such map. Its maps are these: "
                    f"{listed}{more}. {others} {remedy}"
                )
            if len(showing) < len(pathway_numbers):
                self.run.info(f"Maps that show no {noun}s", len(pathway_numbers) - len(showing))
            unified = draw_unified_maps and model['unified_mode'] in value_modes
            individual = draw_individual_maps and model['individual_mode'] in value_modes
            if unified or individual:
                continue
            stem = channel['flag_stem']
            grouped = model.get('group_samples') is not None
            unified_flag = f"--{stem}-{'group' if grouped else 'sample'}-summary"
            sample_flag = f"--{stem}-sample-summary"
            shows, absent = ('shows a border', 'no border') if border else ('shows bars', 'no bars')
            # A border takes a rule for categorical values, and a bar takes numbers only.
            summary = 'an aggregation or a rule' if channel['categorical'] else 'an aggregation'
            # The values of a group come from the sample summary. Without it, no map shows a border
            # or bars. Otherwise a map is missing or has no summary.
            if grouped and model['individual_mode'] is None:
                reason = (
                    f"'{sample_flag}' was not given, so neither the maps of the groups nor the "
                    f"'unified' map {shows}."
                )
                remedy = f"Give '{sample_flag}' {summary}."
            else:
                unified_reason = (
                    "No 'unified' map is drawn." if not draw_unified_maps else
                    f"'{unified_flag}' was not given, so the 'unified' map shows {absent}."
                )
                reason = f"{unified_reason} No map of an individual {individual_noun} is drawn."
                remedy = (
                    f"Draw the maps of the individual {individual_noun}s with "
                    f"'--draw-individual-files' or '--draw-grid'"
                )
                remedy += (
                    f", or give '{unified_flag}' {summary}." if draw_unified_maps else "."
                )
            raise ConfigError(
                f"'{column_flag}' names the {channel['column']} '{model['value_column']}' of the "
                f"{element_type} layer, but no map of this run would show its {marks}. {reason} "
                f"{remedy}"
            )

    def _find_slice_range(
        self,
        pathway_number: str,
        pathway: kgml.Pathway,
        layer: dict,
        unified_valued: bool,
        individuals: List[str]
    ) -> None:
        """
        Collect the values of the slices that a sliced layer of numbers draws on one KEGG map.

        This is the range-finding pass of '_map_elements' for a sliced layer. The scales span the
        values of the slices drawn, not the values of whole elements. A slice that is not drawn is
        left out, whether it has no value, has a share of 0, or belongs to an element drawn whole.
        The maps where a slice is drawn are recorded in '_unified_maps' and '_individual_maps', as
        the values of elements are for a layer that is not sliced.

        An element with more slices than fit is drawn whole in the color of its stand-in value
        ('_slice_cap', '_slice_standin_value'). Only that value is collected for the element, since
        none of its slices is drawn. The check of reserved colors then sees its color
        ('_check_reserved_colors').

        Parameters
        ==========
        pathway_number : str
            Numeric ID of the map.

        pathway : kgml.Pathway
            The map.

        layer : dict
            A sliced layer model of numeric values.

        unified_valued : bool
            True if the 'unified' scale is worked out. The values of the 'unified' map are then
            collected.

        individuals : List[str]
            The samples or groups whose values are collected for the scale their maps share. It is
            empty where those maps show no values of this layer, or none is drawn.
        """
        use_reaction = layer['use_reaction_attribute']

        def _collect(entry, slices, values, period):
            # An element drawn whole shows only its stand-in, so only that value joins the scale.
            if len(slices) > 1 and self._slice_cap(
                pathway, entry, len(slices), ringed=layer.get('border') is not None
            ):
                standin = self._slice_standin_value(
                    [value for _, value, _ in slices], [share for _, _, share in slices], period
                )
                if standin is not None:
                    values.append(standin)
                return
            values.extend(value for _, value, _ in slices)

        for entry in self._find_element_entries(pathway, use_reaction, layer['accessions']):
            accessions = tuple(self._get_entry_kegg_ids(entry, use_reaction))
            if unified_valued:
                slices = self._slice_result(accessions, layer, None)
                if slices:
                    layer['_unified_maps'].add(pathway_number)
                    _collect(
                        entry, slices, layer['_unified_vals'], self._slice_period(layer, True)
                    )
            # A layer without samples draws the values of the 'unified' scale on every map.
            if layer['sample_values'] is None:
                continue
            for individual in individuals:
                slices = self._slice_result(accessions, layer, individual)
                if not slices:
                    continue
                _collect(
                    entry, slices, layer['_individual_vals'], self._slice_period(layer, False)
                )
                if individual in layer['_individual_maps']:
                    layer['_individual_maps'][individual].add(pathway_number)

    def _find_ring_range(
        self,
        pathway_number: str,
        pathway: kgml.Pathway,
        layer: dict,
        unified_valued: bool,
        individuals: List[str]
    ) -> None:
        """
        Collect the numeric values that the sliced rings of a layer's border show on one KEGG map.

        This is the range-finding pass of '_map_elements' for a border with sliced rings, as
        '_find_slice_range' is for a sliced fill. The scales span the values of the bands and
        segments drawn, and the values of the rings drawn whole ('_ring_result'). A unit whose
        border value no ring shows is left out. The maps where a ring is drawn are recorded in the
        border's '_unified_maps' and '_individual_maps'.

        Parameters
        ==========
        pathway_number : str
            Numeric ID of the map.

        pathway : kgml.Pathway
            The map.

        layer : dict
            The model of the fill of a reaction layer whose border holds numbers and is sliced.

        unified_valued : bool
            True if the 'unified' border scale is worked out. The values of the 'unified' map are
            then collected.

        individuals : List[str]
            The samples or groups whose values are collected for the border scale that their maps
            share. It is empty where those maps show no border values, or none is drawn.
        """
        border = layer['border']
        for entry in self._find_element_entries(
            pathway, border['use_reaction_attribute'], border['accessions']
        ):
            if not self._shows_ring(pathway_number, pathway, border, entry):
                continue
            if unified_valued:
                values = self._ring_values(pathway, entry, layer, None)
                if values:
                    border['_unified_maps'].add(pathway_number)
                    border['_unified_vals'].extend(values)
            # A border without samples draws the values of the 'unified' scale on every map.
            if border['sample_values'] is None:
                continue
            for individual in individuals:
                values = self._ring_values(pathway, entry, layer, individual)
                if not values:
                    continue
                border['_individual_vals'].extend(values)
                if individual in border['_individual_maps']:
                    border['_individual_maps'][individual].add(pathway_number)

    def _find_ring_maps(
        self,
        pathway_number: str,
        pathway: kgml.Pathway,
        layer: dict,
        unified_shown: bool,
        individuals: List[str]
    ) -> None:
        """
        Record whether the sliced rings of a layer's categorical border show on one KEGG map.

        A categorical border is keyed only where some map shows one of its rings ('_map_elements').
        The choices among values are not counted here, since no map is drawn yet.

        Parameters
        ==========
        pathway_number : str
            Numeric ID of the map.

        pathway : kgml.Pathway
            The map.

        layer : dict
            The model of the fill of a reaction layer whose border holds categorical values and is
            sliced.

        unified_shown : bool
            True if a map that shows the 'unified' border values is drawn.

        individuals : List[str]
            The samples or groups whose drawn maps show border values. It is empty where none does.
        """
        border = layer['border']
        maps = border['_ring_maps']
        for entry in self._find_element_entries(
            pathway, border['use_reaction_attribute'], border['accessions']
        ):
            # Every Entry is asked, so that the lines of the map are all recorded.
            if not self._shows_ring(pathway_number, pathway, border, entry):
                continue
            if (
                unified_shown and pathway_number not in maps['unified']
                and self._ring_values(pathway, entry, layer, None)
            ):
                maps['unified'].add(pathway_number)
            if pathway_number not in maps['individual'] and any(
                self._ring_values(pathway, entry, layer, individual) for individual in individuals
            ):
                maps['individual'].add(pathway_number)

    def _warn_gene_ids(
        self,
        layers: List[dict],
        draw_unified_maps: bool,
        drawn_individuals: List[str]
    ) -> None:
        """
        Warn about the gene IDs of layers sliced by gene on the maps that pool samples.

        One gene ID names the same gene in every sample. The 'unified' map of a layer with samples,
        and the map of each group, pool a gene's values across samples on that basis. Only the drawn
        maps of this kind that show values and pool two or more samples with rows are checked. Two
        warnings can follow, each given once per layer.

        One warning names the gene IDs that have no accession in common between their rows in two
        samples of such a map. A gene of one contigs database is annotated once, so it has the same
        accessions in every sample. A sample can lack rows for some of them, so sets of accessions
        that merely differ are not reported. Sets with no accession in common suggest that one ID
        names two genes. Gene caller IDs start at 0 in every contigs database, so the IDs of
        separate assemblies collide unless they are made unique. Colliding IDs pool unrelated genes
        as one.

        The other warning names the maps where no gene ID occurs in two or more samples. Each slice
        of a gene there shows the value of a single sample, and nothing is pooled. This is expected
        where the IDs of separate assemblies were made unique.

        A border whose rings are sliced by gene pools the border values of its genes in the same
        way, so its gene IDs are checked too. The gene IDs and maps that the fill of the same file
        already reported are not reported again.

        Parameters
        ==========
        layers : List[dict]
            The layer models and border models. Only models sliced by gene are read.

        draw_unified_maps : bool
            True if the 'unified' map is drawn.

        drawn_individuals : List[str]
            The samples or groups drawn on their own maps.
        """
        def join_names(names: List[str]) -> str:
            # Two names read "'a' and 'b'". Three or more take commas and 'and' before the last.
            if len(names) > 1:
                names = names[:-1] + [f"and {names[-1]}"]
            return (', ' if len(names) > 2 else ' ').join(names)

        value_modes = ('quantitative', 'categorical')
        # The gene IDs and the maps already reported for each file.
        reported: Dict[str, Tuple[Set[str], Set[Union[str, None]]]] = {}
        for layer in layers:
            if layer.get('slice_by') != 'gene' or layer['gene_sample_accessions'] is None:
                continue
            gene_sample_accessions = layer['gene_sample_accessions']
            with_rows = {
                sample for sample_accessions in gene_sample_accessions.values()
                for sample in sample_accessions
            }
            pooled = []
            if draw_unified_maps and layer['unified_mode'] in value_modes:
                pooled.append((None, set(layer['sample_values'])))
            group_samples = layer['group_samples']
            if group_samples is not None and layer['individual_mode'] in value_modes:
                pooled.extend((group, set(group_samples[group])) for group in drawn_individuals)
            # A map that pools fewer than two samples with rows pools no genes.
            unshared = []
            colliding = set()
            for name, samples in pooled:
                if len(samples & with_rows) < 2:
                    continue
                shared = False
                for gene, sample_accessions in gene_sample_accessions.items():
                    found = [
                        accessions for sample, accessions in sample_accessions.items()
                        if sample in samples
                    ]
                    if len(found) < 2:
                        continue
                    shared = True
                    if gene in colliding:
                        continue
                    # A sample can lack rows for some accessions of a gene, so only sets with no
                    # accession in common are reported. Each distinct set is compared once.
                    distinct = {frozenset(accessions) for accessions in found}
                    if any(a.isdisjoint(b) for a, b in combinations(distinct, 2)):
                        colliding.add(gene)
                if not shared:
                    unshared.append(name)

            element_type = layer['element_type']
            path = layer['path']
            reported_genes, reported_maps = reported.setdefault(path, (set(), set()))
            colliding -= reported_genes
            reported_genes.update(colliding)
            unshared = [name for name in unshared if name not in reported_maps]
            reported_maps.update(unshared)
            # The rows of a border are those with a border value ('_channel_view').
            rows = 'rows' if layer['channel']['role'] == 'fill' else 'rows with border values'
            if colliding:
                genes = self._sort_values(colliding)
                several = len(genes) > 1
                shown = join_names([repr(gene) for gene in genes[:5]])
                if len(genes) > 5:
                    listed = f"The first 5 are {shown}."
                else:
                    listed = f"The IDs are {shown}." if several else f"The ID is {shown}."
                self.run.warning(
                    f"{len(genes)} gene {'IDs' if several else 'ID'} of the {element_type} file at "
                    f"'{path}' {'have' if several else 'has'} no accession in common between "
                    f"{'their' if several else 'its'} {rows} in two samples of a map that pools "
                    f"samples. {listed} A gene of one contigs database is annotated once, so it "
                    f"has the same accessions in every sample. A sample can lack rows for some of "
                    f"them, but rows with no accession in common suggest that one ID names two "
                    f"genes. Gene caller IDs start at 0 in every contigs database, so the IDs of "
                    f"separate assemblies must be made unique, for example with a sample prefix. "
                    f"Otherwise, unrelated genes are pooled as one."
                )
            if not unshared:
                continue
            groups = [name for name in unshared if name is not None]
            parts = ["the 'unified' map"] if None in unshared else []
            if len(groups) > 5:
                parts.append(f"the maps of {len(groups)} groups")
            elif groups:
                plural = 's' if len(groups) > 1 else ''
                joined = join_names([f"'{group}'" for group in groups])
                parts.append(f"the map{plural} of group{plural} {joined}")
            if layer['channel']['role'] == 'fill':
                ids = f"No gene ID of the {element_type} file at '{path}'"
                shows = "Each slice of a gene there shows the value of a single sample"
            else:
                ids = f"No gene ID with border values in the {element_type} file at '{path}'"
                shows = "Each slice of a ring there shows the border value of a single sample"
            self.run.warning(
                f"{ids} occurs in two or more samples of {' and '.join(parts)}. {shows}, and "
                f"nothing is pooled. This is expected where the gene IDs of separate assemblies "
                f"were made unique."
            )

    def _warn_missing_bars(self, bars: List[dict]) -> None:
        """
        Warn about the boxes that the maps with bars left uncolored for want of a bar value.

        On a map that shows bars, a box with a fill color and no bar value is left uncolored
        ('_stage_entry_bar'). Each box is counted once on each drawn map, such as the 'unified' map
        of 00010 or the map of 00010 of one sample. Bars in the 'strip' style leave no box
        uncolored, since such a box keeps its whole fill without a strip.

        Parameters
        ==========
        bars : List[dict]
            The bar models of the run. Each keeps the boxes of the drawn maps in '_missing'.
        """
        for bar in bars:
            missing = bar['_missing']
            if not missing:
                continue
            maps = {(pathway_number, drawn_map) for pathway_number, drawn_map, _ in missing}
            examples = ', '.join(sorted({
                accession for accessions in missing.values() for accession in accessions
            })[:5])
            self.run.warning(
                f"On {len(maps)} drawn map(s), {len(missing)} {bar['element_type']} box(es) with a "
                f"fill color have no bar value from the '{bar['value_column']}' bar column, so "
                f"they were left uncolored. A whole fill would read as a full bar, and a bare "
                f"track as a bar of 0. Their accessions include these: {examples}. A box has no "
                f"bar value where none of its accessions has one on the map, or where reducing "
                f"their bar values is undefined."
            )

    def _write_label_keys(self, bars: List[dict], output_dir: str) -> None:
        """
        Write the key table of the boxes that bars drew without their labels.

        A box of the 'fill' style on a bar scale with values on both sides of 0 drops its label
        ('_stage_entry_bar'). The key table lists each such box on each drawn map
        ('LABEL_KEY_COLUMNS'), so that a box can still be told by its accessions. Its name is
        'label_key_<name>.txt', such as 'label_key_reactions.txt', in the output directory. It is
        written for the bars of the 'fill' style whose scale on the 'unified' map or on the
        individual maps has values on both sides of 0, as '_map_elements' found before drawing. The
        rows are sorted as those of the key table of slices ('_write_slice_keys'): by KEGG map, by
        drawn map in the order the maps are drawn, and by entry ID.

        Parameters
        ==========
        bars : List[dict]
            The bar models of the run. Each says in '_signed' whether its table is written, and
            keeps the rows of the drawn maps in '_unlabeled'.

        output_dir : str
            The output directory.
        """
        for bar in bars:
            if not bar['_signed']:
                continue
            unlabeled = bar['_unlabeled']
            path = os.path.join(output_dir, f"label_key_{bar['name']}.txt")
            with open(path, 'w') as key_file:
                key_file.write('\t'.join(LABEL_KEY_COLUMNS) + '\n')
                for key in sorted(unlabeled):
                    key_file.write('\t'.join(unlabeled[key]) + '\n')
            self.run.info(f"Key table of the {bar['element_type']} boxes without labels", path)

    def _write_slice_keys(self, layers: List[dict], output_dir: str) -> None:
        """
        Write the key table of each sliced layer, and warn about the elements that were not sliced.

        The key table lists every slice drawn on every drawn map ('SLICE_KEY_COLUMNS'). A slice of a
        gene names the gene in 'gene_id', and its accessions in the element in 'accession', in KGML
        order and joined by commas. Its name is 'slice_key_<layer name>.txt', such as
        'slice_key_reactions.txt', in the output directory. Each row names a KEGG map and a drawn
        map ('unified', or a sample or group), not a file. An element drawn whole over the cap gets
        one row for each slice it would have had. Such a row has no position, the color of the whole
        element, and the reason in 'drawn_whole'. The rows are sorted by KEGG map, by drawn map in
        the order the maps are drawn, by entry ID, and by slice.

        Two warnings follow, each counting an element once on each drawn map. One counts the
        elements drawn whole by KEGG map. The other counts the elements left uncolored, since every
        slice with a value had a share of 0.

        A border with sliced rings has a key table of the same form,
        'slice_key_reactions_border.txt' ('_stage_entry_ring'). The positions of bands are those of
        the columns of the fill. The positions of segments count clockwise from the top left corner.
        A ring drawn whole for a reason has a row for each slice it stands for, and the reason is
        'over the cap' or 'fill not sliced'. Its warnings count the rings drawn whole over the cap,
        the rings left out by shares of 0, and the units with border values that no column of the
        fill shows.

        Parameters
        ==========
        layers : List[dict]
            The layer models and border models. Only sliced models are read.

        output_dir : str
            The output directory.
        """
        for layer in layers:
            if layer.get('slice_by') is None:
                continue
            element_type = layer['element_type']
            border = layer['channel']['role'] == 'border'
            path = os.path.join(output_dir, f"slice_key_{layer['name']}.txt")
            with open(path, 'w') as key_file:
                key_file.write('\t'.join(SLICE_KEY_COLUMNS) + '\n')
                for key in sorted(layer['_slice_key']):
                    for row in layer['_slice_key'][key]:
                        key_file.write('\t'.join(row) + '\n')
            self.run.info(f"Key table of the {layer['channel']['noun']} slices", path)

            units = 'genes' if layer['slice_by'] == 'gene' else 'accessions'
            if layer.get('value_type') == 'categorical':
                standin = f"the value chosen among the {units} of its slices"
            elif layer.get('value_period') is not None:
                standin = "the circular mean of its slice values, weighted by their shares"
            else:
                standin = "the mean of its slice values, weighted by their shares"
            reason_column = SLICE_KEY_COLUMNS.index('drawn_whole')
            whole = Counter(
                key[0] for key, rows in layer['_slice_key'].items()
                if rows[0][reason_column] and rows[0][reason_column] != 'fill not sliced'
            )
            if whole and border:
                by_map = ', '.join(f"{number} ({count})" for number, count in sorted(whole.items()))
                if layer['slice_layout'] == 'aligned':
                    why = (
                        f"since the fill of each of their boxes had more slices than fit. Such a "
                        f"ring takes the color of {standin}. Only the slices of the fill that have "
                        f"border values count"
                    )
                else:
                    why = (
                        f"since each would have had more than {MAX_RING_SLICES} slices. Such a "
                        f"ring takes the color of {standin}"
                    )
                self.run.warning(
                    f"{sum(whole.values())} ring(s) of {element_type} elements were drawn whole in "
                    f"one color, {why}. A ring is counted once on each drawn map. These are the "
                    f"counts by KEGG map: {by_map}. The key table at '{path}' lists these rings "
                    f"with the reason."
                )
            elif whole:
                by_map = ', '.join(f"{number} ({count})" for number, count in sorted(whole.items()))
                # In a layer with a border column, the slices of a box fit inside its ring.
                inside = (
                    " The slices of a box fit inside its ring, since the layer has a border column."
                    if layer.get('border') is not None else ""
                )
                self.run.warning(
                    f"{sum(whole.values())} {element_type} element(s) had more slices than fit, so "
                    f"each was drawn whole in one color. A box is drawn whole where a slice of "
                    f"equal width would be narrower than {MIN_BOX_SLICE_WIDTH} points.{inside} A "
                    f"line is drawn whole where a slice of equal length would be shorter than the "
                    f"line is wide. Only the part of the line that compound circles leave visible "
                    f"counts. Such an element takes the color of {standin}. An element is counted "
                    f"once on each drawn map. These are the counts by KEGG map: {by_map}. The key "
                    f"table at '{path}' lists these elements with the reason."
                )
            uncolored = layer['_slice_uncolored']
            if uncolored:
                examples = ', '.join(sorted({
                    accession for accessions in uncolored.values() for accession in accessions
                    if accession in layer['accessions']
                })[:5])
                if border:
                    every = 'slice of the fill' if layer['slice_layout'] == 'aligned' else (
                        f"{units[:-1]} with a border value"
                    )
                    self.run.warning(
                        f"{len(uncolored)} {element_type} element(s) show no border, since every "
                        f"{every} had a weight of 0 on that map. An element is counted once on "
                        f"each drawn map. Their accessions in the file include these: {examples}."
                    )
                else:
                    self.run.warning(
                        f"{len(uncolored)} {element_type} element(s) were left uncolored, since "
                        f"every slice with a value had a weight of 0 on that map. An element is "
                        f"counted once on each drawn map. Their accessions in the file include "
                        f"these: {examples}."
                    )
            unshown = layer.get('_slice_unshown')
            if unshown:
                found = {unit for found_units in unshown.values() for unit in found_units}
                examples = ', '.join(sorted({
                    accession for unit in found for accession in (
                        unit[1] if isinstance(unit, tuple) else (unit,)
                    )
                })[:5])
                named = (
                    "The accessions of these genes include these" if layer['slice_by'] == 'gene'
                    else "They include these"
                )
                # A legend of categorical values lists every value, shown or not.
                scale = " The border scale leaves these values out." if (
                    layer.get('value_type') != 'categorical'
                ) else ""
                self.run.warning(
                    f"On {len(unshown)} drawn map element(s), some {units} with border values had "
                    f"no slice of the fill, so no band shows their border values.{scale} An "
                    f"element is counted once on each drawn map. {named}: {examples}."
                )

    def _check_reserved_colors(
        self,
        layers: List[dict],
        pathway_numbers: List[str],
        drawn_individuals: List[str],
        group_color_priorities: Dict[str, List[Tuple[str, float]]] = None,
        group_element_types: Set[str] = None,
        draw_unified_maps: bool = True
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

        drawn_individuals : List[str]
            The samples, contigs databases, genomes, or groups drawn on their own maps. Only the
            colors that the membership colors file gives these are checked here. Without groups,
            each such color is the color of a map. With groups, a group's own maps take the group's
            color only with '--group-colormap tint'. Where the 'unified' map is drawn, its own scale
            is checked as well. That scale includes the color of every sample, contigs database,
            genome, or group.

        group_color_priorities : Dict[str, List[Tuple[str, float]]], None
            The colors of each group's individual maps ('_group_map_colors'), which belong to no one
            layer: one scale colors the within-group source counts of every layer whose per-group
            context is presence. A ramp built from a group's own color runs towards white by
            construction, so it can reach a reserved color even where a layer's own colors do not.

        group_element_types : Set[str], None
            The element types the group scale colors, deciding which reserved colors apply to it.

        draw_unified_maps : bool, True
            True if the 'unified' map is drawn. The colors of categorical values are checked only
            where a map colored by them is drawn.
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
                f"{GROUP_COLORMAP_TINT}'), start that ramp further from white with the two limits "
                f"that option also takes."
            )

        # The value scale shared by the individual maps colors only those maps.
        value_scales = [('_unified_norm', '_unified_center', '_unified_vals', 'cmap')]
        if drawn_individuals:
            value_scales.append(
                ('_individual_norm', '_individual_center', '_individual_vals', 'individual_cmap')
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
            # Every categorical value can be drawn on a map colored by the values, and the legend
            # shows them all. A layer with samples colors its 'unified' map by presence by default,
            # so a run drawing that map alone draws no value.
            value_colors: Set[str] = set()
            if layer.get('value_color_priorities') is not None and (
                (draw_unified_maps and layer['unified_mode'] == 'categorical')
                or (drawn_individuals and layer['individual_mode'] == 'categorical')
            ):
                value_colors = {
                    kgml.canonical_color(color)
                    for color, _ in layer['value_color_priorities'].values()
                }
                staged.update(value_colors)
            membership_colors = layer.get('membership_colors')
            if membership_colors is not None:
                staged.update(membership_colors[individual] for individual in drawn_individuals)
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
            # Anvi'o's default colors of categorical values are never reserved, so a clash there
            # comes from the colors file.
            clashing_values = [color for color in clashing if color in value_colors]
            if clashing_values:
                raise ConfigError(
                    f"The colors file given to '--{layer['element_type']}-value-colors' gives "
                    f"categorical values {'colors' if len(clashing_values) > 1 else 'a color'} "
                    f"that the pathway maps keep for their own unidentified elements: "
                    f"{', '.join(clashing_values)}. Anvi'o recolors the reactions and compounds a "
                    f"map does not highlight. They are gray on global maps. Elsewhere, reactions "
                    f"are black and compounds are white. A highlighted element in one of those "
                    f"colors would be invisible. Please give these values other colors in the file."
                )
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

        A global or overview map colors each compound from the reactions it touches
        ('kgml.Pathway.color_associated_compounds'), where a reaction layer stages fills and no
        compound layer does. Where the reactions are colored by value, a compound takes the color at
        the average position of its reactions on the color scale. Otherwise it takes the color of
        its reaction drawn on top. That color can be the one the map keeps for the compounds it does
        not highlight. This can happen even where no reaction has that color. 'kgml' then refuses to
        draw the map. This check finds such a clash before any map or colorbar is drawn.
        '_check_reserved_colors' does the same for the colors of the layers themselves.

        The derived colors depend on which reactions touch each compound. Each map is therefore
        colored here without being drawn. That costs about as much as drawing the map, less the
        rendering, so it is done only on the maps whose reserved color a derived color can reach at
        all. A derived color is always a color of a reaction colormap or a fixed color of the
        reactions.

        Parameters
        ==========
        layers : List[dict]
            The layer models. Compound colors are derived only where a reaction layer stages fills
            and no compound layer does. A layer in the fill mode 'none' stages no fills.

        pathway_numbers : List[str]
            The maps about to be drawn. Only global and overview maps derive compound colors.

        element_map_specs : Callable[[], Iterable[Tuple[Union[str, None], bool, List[dict]]]]
            Gives the layer specs of each kind of map that '_draw_map_elements' is about to draw,
            such as the 'unified' map or the map of one sample. Each comes with a phrase naming that
            kind of map, and with True if it is the 'unified' map. The phrase is None where the run
            has a single source of data. Such a run draws a single map. The specs are built only if
            the maps have to be colored.

        group_color_priorities : Dict[str, List[Tuple[str, float]]], None
            The colors of each group's individual maps ('_group_map_colors'). A compound on a
            group's map can take one of them.
        """
        # Compound colors are derived from a reaction layer that stages fills. A compound layer that
        # stages fills colors the compounds itself, so no compound color is derived then. A layer in
        # the fill mode 'none' stages none. A compound layer in this mode draws only rings, so the
        # compounds keep their derived fills ('_derives_compounds').
        fills = {layer['element_type'] for layer in layers if layer['unified_mode'] != 'none'}
        if 'reaction' not in fills or 'compound' in fills:
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
                ('unified_mode', 'cmap'), ('individual_mode', 'individual_cmap')
            ):
                if layer[mode_key] != 'quantitative':
                    continue
                for cmap in (layer[cmap_key], layer[cmap_key].reversed()):
                    derivable.update(mcolors.rgb2hex(color) for color in cmap(np.arange(cmap.N)))
            if '_colors' in layer:
                derivable.update(color for color, _ in layer['_colors'][1])
            if layer.get('value_color_priorities') is not None:
                derivable.update(color for color, _ in layer['value_color_priorities'].values())
            if layer.get('membership_colors') is not None:
                derivable.update(layer['membership_colors'].values())
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
                color_priority, _, _ = self._stage_map_elements(pathway_number, pathway, specs)
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
                    f"compound layer colors the compounds. The map therefore colors each compound "
                    f"from the reactions it touches. Where the reactions are colored by value, a "
                    f"compound takes the color at the average position of its reactions on the "
                    f"color scale. Where they are colored by categorical values, a rule chooses "
                    f"one of their values for it. Otherwise it takes the color of its reaction "
                    f"drawn on top. {remedy}"
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
            The layer models. No compound layer of the run stages fills, and its one reaction layer
            is here.

        unified : bool
            True for the 'unified' map. False for the map of one sample, contigs database, genome,
            or group.

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
        mode = layer['unified_mode'] if unified else layer['individual_mode']
        if reserved_color == '#ffffff':
            colormaps = (
                "Grayscale colormaps and those running to pure white, such as 'Greys', 'hot' and "
                "'bone', reach it."
            )
        else:
            colormaps = "Colormaps with grays in them, such as 'Greys' and 'RdGy', reach it."
        if not unified and group_color_priorities and mode == 'membership':
            return (
                f"These maps color the number of the group's own sources containing each reaction, "
                f"styled by '--group-colormap'. Choose a colormap that does not reach "
                f"{reserved_color}. {colormaps} A ramp from a pale tint to the group's own color "
                f"('--group-colormap {GROUP_COLORMAP_TINT}') can start further from white instead. "
                f"Raise the first of the two limits that option also takes."
            )
        if mode == 'quantitative' and not unified:
            return (
                f"Choose a colormap for '--reaction-individual-colormap' that does not reach "
                f"{reserved_color}. Without that option, these maps use the colormap given to "
                f"'--reaction-colormap'. {colormaps}"
            )
        # The colors of categorical values come from a colors file. Anvi'o's default colors are
        # never reserved. A compound takes the color of one of the values. The file gives each value
        # a color of its own, so exactly one value has the reserved color.
        if mode == 'categorical':
            value = next(
                value for value, (color, _) in layer['value_color_priorities'].items()
                if kgml.canonical_color(color) == reserved_color
            )
            return (
                f"Give the value '{value}' a color other than {reserved_color} in the file given "
                f"to '--reaction-value-colors'."
            )
        # The membership colors file gives each sample, contigs database, genome, or group a color.
        # Presence is then colored by those colors in place of a colormap.
        if mode == 'membership' and layer.get('membership_colors') is not None:
            return (
                f"Remove {reserved_color} from the file given to "
                f"'{layer['membership_colors_flag']}'."
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
        it to the total count. That is the number of samples, contigs databases, genomes, or groups
        being counted. This keeps the scale the same however few of them the data reaches, so the
        scale is comparable between runs. A number pins the scale, which is how separate figures are
        given one scale when their data differs.

        Parameters
        ==========
        count_scale_max : Union[str, int]
            'observed', 'total', or the count to stop at.

        observed : int
            The highest count anything on the drawn maps has, 0 if nothing was counted.

        total : int
            How many samples, contigs databases, genomes, or groups are counted in all.

        Returns
        =======
        int
            The count the scale runs up to.
        """
        if count_scale_max == 'observed':
            # An observed maximum of 0 means nothing on the drawn maps was counted at all. A scale
            # of no counts cannot be drawn, so the scale runs to the total count instead. The layer
            # colors nothing either way, so which of the two it is never shows on a map.
            top = observed if observed > 0 else total
        elif count_scale_max == 'total':
            top = total
        else:
            top = int(count_scale_max)

        # Whatever the choice, a scale reaches at least 1. An element colored by its count is in at
        # least one sample, contigs database, genome, or group.
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
        individual_finders = [
            self._membership_finder(
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
            for layer, find_individuals in zip(unified_layers, individual_finders):
                for entry in self._find_element_entries(
                    pathway, layer['use_reaction_attribute'], layer['membership']
                ):
                    individuals = find_individuals(entry)
                    if individuals is not None and len(individuals) > layer['_observed_count']:
                        layer['_observed_count'] = len(individuals)
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
        draw_individuals: List[str],
        group_observed_counts: Dict[str, int],
        count_scale_max: Union[str, int],
        unified_plural: str,
        individual_noun: str
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
        'GROUP_COLORMAP_TINT' is asked for, from a ramp per group running from a pale tint to that
        group's own color ('tint_hexcode'). The ramp binds the identity channel to the magnitude
        channel: every panel of a grid is built the same way and so stays comparable, while its hue
        says which group it is. A named colormap is the default because one engineered for magnitude
        reads better as magnitude than a color chosen to identify a group does.

        Parameters
        ==========
        grouped_membership : dict
            The grouping record, carrying 'group_sources' and the group style ('group_colormap'/
            'group_colormap_limits'/'group_reverse_overlay'/'group_colormap_scheme').

        layers : List[dict]
            Every layer of the run. The membership colors file of a layer, if given, sets the color
            of each group. Each group's ramp runs to that color. Since one style covers every
            layer's group maps, two layers coloring one group differently is refused.

        draw_individuals : List[str]
            The groups whose own maps are drawn.

        group_observed_counts : Dict[str, int]
            The highest within-group source count found on the drawn maps, per group, or empty when
            the scale does not stop at what was observed.

        count_scale_max : Union[str, int]
            Where each group's count scale stops ('_resolve_count_scale_top').

        unified_plural : str
            The plural of what a group's sources are, for messages, e.g. 'samples'.

        individual_noun : str
            What one group is called, for messages, e.g. 'sample group'.

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

        tinted = group_colormap == GROUP_COLORMAP_TINT
        group_colors: Dict[str, str] = {}
        group_cmap = None
        if tinted:
            # One ramp per group needs one color per group, and one style covers every layer's group
            # maps, so the layers that have colors must agree about them.
            colored_layers = [layer for layer in layers if layer.get('membership_colors')]
            if not colored_layers:
                self.progress.end()
                raise ConfigError(
                    f"The group colormap was given as '{GROUP_COLORMAP_TINT}', which colors each "
                    f"group's own maps by a ramp running from a pale tint to that group's own "
                    f"color, but no color was given for any group. Give each group a color in a "
                    f"kegg-membership-colors-txt file. Pass that file to the membership colors "
                    f"option of a layer colored by presence, such as "
                    f"'--reaction-membership-colors'. Or name a Matplotlib colormap for the group "
                    f"maps instead."
                )
            # A ramp is built only for the groups whose own maps are drawn, so they are the only
            # ones the layers must agree about: a group the run has but draws no map for has no ramp
            # to disagree over. Every one of them is guaranteed a color by
            # '_resolve_membership_colors'.
            for layer in colored_layers:
                for group in draw_individuals:
                    color = layer['membership_colors'][group]
                    if group_colors.setdefault(group, color) != color:
                        self.progress.end()
                        raise ConfigError(
                            f"The group colormap was given as '{GROUP_COLORMAP_TINT}', so each "
                            f"group's own maps are colored by a ramp running to that group's "
                            f"color. One ramp per group styles every layer's group maps at once, "
                            f"but the reaction and compound layers were given different colors for "
                            f"the {individual_noun} '{group}': '{group_colors[group]}' and "
                            f"'{color}'. Please give that {individual_noun} one color in both "
                            f"files — pointing both options at a single file is the simplest way — "
                            f"or name a Matplotlib colormap for the group maps instead."
                        )
            if colormap_limits is None:
                colormap_limits = DEFAULT_GROUP_TINT_SPAN
            if not 0.0 <= colormap_limits[0] <= colormap_limits[1] <= 1.0:
                self.progress.end()
                raise ConfigError(
                    f"The two limits of a group colormap of '{GROUP_COLORMAP_TINT}' are how far "
                    f"from white towards a group's own color its ramp starts and stops, so they "
                    f"must lie between 0.0 and 1.0 with the smaller one first. These do not: "
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
            for group in draw_individuals
        }

        group_color_priorities: Dict[str, List[Tuple[str, float]]] = {}
        group_distinct_colors: Dict[str, int] = {}
        for group in draw_individuals:
            group_scale_top = group_scale_tops[group]
            if group_scale_top == 1:
                sample_points = np.linspace(1, 1, 1)
            else:
                sample_points = np.linspace(0, 1, group_scale_top)
            if tinted:
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
            group for group in draw_individuals
            if count_scale_max != 'observed' or group_observed_counts.get(group, 0) > 0
        ] or list(draw_individuals)
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
                f"The ramp running to the color of {individual_noun} '{group}' could supply only "
                f"{group_distinct_colors[group]} colors that can be told apart, fewer than the "
                f"{group_scale_tops[group]} counts its scale runs over. Widening the ramp's span "
                f"with the two limits of the group colormap option gives it more room."
            ) if tinted else (
                f"The group colormap could supply only {group_distinct_colors[group]} colors that "
                f"can be told apart, fewer than the {group_scale_tops[group]} counts the scale of "
                f"{individual_noun} '{group}' runs over."
            )

        def _crowded_group_clause(group: str) -> str:
            return (
                f"The scale of {individual_noun} '{group}' runs over {group_scale_tops[group]} "
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
                    f"{_worst_group_clause()} The count on individual {individual_noun} maps is "
                    f"therefore drawn on a CONTINUOUS color scale rather than in discrete bands of "
                    f"one color per count. Each colorbar is a gradient running from a count of 1 "
                    f"to the top of that {individual_noun}'s scale, on which a color reads as a "
                    f"position along that range rather than as an exact count. Ask for this scale "
                    f"explicitly with '{GROUP_SCHEME_OPTIONS['by_count_continuous']}', or ask for "
                    f"'{GROUP_SCHEME_OPTIONS['by_count']}' to insist on the discrete bands and be "
                    f"told when they cannot be drawn. One scheme covers every {individual_noun}'s "
                    f"maps, so the {individual_noun} that runs furthest past a limit is the one "
                    f"reported here."
                )
        else:
            # The bands were asked for outright, so each group that runs past a limit is reported
            # once, by the limit that says the most about it.
            for group in short_groups:
                self.run.warning(
                    f"{_short_group_clause(group)} Neighboring counts therefore share a color on "
                    f"that {individual_noun}'s individual maps, and its colorbar labels more bands "
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
                    f"but that {individual_noun}'s colorbar labels will be very small. Drawing the "
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
            for group in draw_individuals:
                if tinted:
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

    def _membership_layer_given_colors(
        self,
        layer: dict,
        individuals: List[str],
        count_scale_top: int
    ) -> Tuple[str, List[Tuple[str, float]], List[Tuple[str]], None]:
        """
        Resolve a membership layer's colors from a kegg-membership-colors-txt file.

        This is the 'membership_colors' branch of '_membership_layer_colors', which see: it returns
        the same record, and the drawing and colorbar code cannot tell where the colors came from.
        Each sample, contigs database, genome, or group alone is drawn in its own color. A
        combination of them is drawn in the blend of their colors ('blend_hexcodes'). A row of the
        colors file can override the color of a combination. Drawing priority follows the order of
        the combinations, exactly as it does when the colors are sampled from a colormap, so an
        element in more samples, contigs databases, genomes, or groups is still drawn over one in
        fewer.

        Parameters
        ==========
        layer : dict
            The layer model, with 'membership_colors', 'membership_combo_colors', the flag they came
            from, and 'reverse_overlay'.

        individuals : List[str]
            The samples, contigs databases, genomes, or groups whose membership colors the layer, in
            color-assignment order.

        count_scale_top : int
            The count a color scale would run up to, named in the message that points at coloring by
            count instead.

        Returns
        =======
        Tuple[str, List[Tuple[str, float]], List[Tuple[str]], None]
            'by_membership', the (color_hexcode, priority) pairs in combination order, the
            combinations themselves, and None in place of the colormap the colors did not come from.
        """
        membership_colors: Dict[str, str] = layer['membership_colors']
        combo_colors: Dict[Tuple[str, ...], str] = layer['membership_combo_colors']
        flag = layer['membership_colors_flag']
        scheme_options = layer.get('scheme_options', PRESENCE_SCHEME_OPTIONS)
        reverse_overlay = layer.get('reverse_overlay', False)

        # A count says how many samples, contigs databases, genomes, or groups contain an element.
        # It does not say which ones, so a count cannot take the color of any one of them. Refused
        # rather than ignored, since the two were asked for together and only one of them can be
        # honored.
        colormap_scheme = layer.get('colormap_scheme')
        if colormap_scheme is not None and colormap_scheme != 'by_membership':
            self.progress.end()
            raise ConfigError(
                f"The {layer['element_type']} layer was given a color for each sample, contigs "
                f"database, genome, or group by '{flag}'. It was also asked to be colored by count "
                f"with '{scheme_options[colormap_scheme]}'. A count is how many of them contain an "
                f"element. It does not say which ones, so a count cannot take the color of any one "
                f"of them. The two requests cannot both be honored. Drop '{flag}' to color the "
                f"counts from a colormap. Or ask for '{scheme_options['by_membership']}' to color "
                f"each element by the colors that '{flag}' gives."
            )

        # Coloring by membership needs a color for each combination of the samples, contigs
        # databases, genomes, or groups. The number of combinations doubles with each one added. The
        # ceiling is checked before the combinations are enumerated, so a large number of them is
        # refused before gigabytes of memory are spent.
        combo_count = 2 ** len(individuals) - 1
        if combo_count > MAX_MEMBERSHIP_COLOR_COMBOS:
            self.progress.end()
            combos_shown = f'{combo_count}' if combo_count <= 10 ** 6 else f'{combo_count:.1e}'
            raise ConfigError(
                f"Coloring the {layer['element_type']} layer by membership needs a distinct color "
                f"for every combination of the samples, contigs databases, genomes, or groups. "
                f"There are {len(individuals)} of them. That makes {combos_shown} combinations. "
                f"Anvi'o blends the colors given by '{flag}' to derive them, and refuses past "
                f"{MAX_MEMBERSHIP_COLOR_COMBOS} combinations, well beyond the handful that any "
                f"color scale could tell apart. Color by count instead with "
                f"'{scheme_options['by_count']}' and drop '{flag}', since a count takes its colors "
                f"from a colormap and the two cannot both be honored; it needs just "
                f"{count_scale_top} colors. Alternatively, reduce the number of samples, contigs "
                f"databases, genomes, or groups. For example, group them."
            )

        membership_combos: List[Tuple[str]] = []
        for combo_size in range(1, len(individuals) + 1):
            membership_combos += list(combinations(individuals, combo_size))

        color_priorities: List[Tuple[str, float]] = []
        blended: List[Tuple[str]] = []
        for combo_index, combo in enumerate(membership_combos):
            if len(combo) == 1:
                color = membership_colors[combo[0]]
            else:
                override = combo_colors.get(tuple(sorted(combo)))
                if override is None:
                    color = blend_hexcodes(membership_colors[individual] for individual in combo)
                    blended.append(combo)
                else:
                    color = override
            # The position in the combination order is both the color's place on the scale and its
            # drawing priority, which 'reverse_overlay' inverts. Counting from the last combination
            # rather than subtracting from 1 keeps every priority non-negative, as a Pathway
            # requires.
            color_priorities.append((
                color,
                1 - combo_index / len(membership_combos) if reverse_overlay
                else (combo_index + 1) / len(membership_combos)
            ))

        # A discrete colorbar labels one band per combination, so two combinations sharing a color
        # would leave the bar with more labels than colors a reader can tell apart. Blending usually
        # causes this. Two blends can come out as one color. A blend can also come out as the color
        # of one sample, contigs database, genome, or group, so the message points at the rows that
        # would fix it.
        distinct = len({color for color, _ in color_priorities})
        if distinct != len(membership_combos):
            self.progress.end()
            blend_clause = (
                f" {len(blended)} of the {len(membership_combos)} combinations took a color "
                f"blended from their members' colors; a row of the file naming a combination, with "
                f"its names separated by '{MEMBERSHIP_COMBO_SEPARATOR}', sets that combination's "
                f"color directly instead."
            ) if blended else ""
            raise ConfigError(
                f"The colors of the {layer['element_type']} layer of this map could not be "
                f"assigned. Coloring by membership needs a distinct color for every combination of "
                f"the samples, contigs databases, genomes, or groups. There are {len(individuals)} "
                f"of them. That makes {len(membership_combos)} combinations. The colors from "
                f"'{flag}' came to only {distinct}.{blend_clause} Alternatively, color by count "
                f"with '{scheme_options['by_count']}' and drop '{flag}', since a count takes its "
                f"colors from a colormap and the two cannot both be honored; it needs just "
                f"{count_scale_top} colors."
            )

        return 'by_membership', color_priorities, membership_combos, None

    def _membership_layer_colors(
        self,
        layer: dict,
        individuals: List[str],
        count_scale_top: int = None
    ) -> Tuple[str, List[Tuple[str, float]], Union[List[Tuple[str]], None], mcolors.Colormap]:
        """
        Resolve a membership layer's coloring scheme, colors and priorities, and combinations.

        Resolves the by-count/by-membership colormap logic for one layer. 'individuals' are the
        samples, contigs databases, genomes, or groups whose count or membership colors the layer.
        'count_scale_top' is the count the color scale runs up to. It defaults to all of them
        ('_resolve_count_scale_top'). Stopping the scale where the data does spreads the colors over
        the counts that occur instead of over counts nothing reaches, which for sparse data is the
        difference between a readable map and one colored in a single shade; the cost is that the
        scale then depends on what was drawn.

        The colors come back in the order the scheme assigns them. That is by ascending count, or by
        the order of the combinations. They are not keyed by color, because a continuous count scale
        ('by_count_continuous') deliberately gives the same color to neighboring counts wherever the
        colormap has no distinguishable color left for each one. The two count schemes sample the
        same colormap and differ only in how the scale is drawn, and therefore in whether the colors
        they assign must be distinguishable: a discrete colorbar labels one band per count, so a
        count without its own color would leave the bar mislabeled, whereas a gradient from the
        first count to the last stays honest however many counts share a color. For the sequential
        colormap a count scale calls for, the colors the two assign are identical; only a
        qualitative colormap makes them differ, since 'by_count' samples it at whole positions while
        a gradient has to span a range.

        A layer given a color for each sample, contigs database, genome, or group
        ('membership_colors') is colored by membership from those colors instead of from a colormap.
        Each one alone takes its own color. A combination of them takes the blend of their colors. A
        row of the colors file can override that color. A count does not say which samples, contigs
        databases, genomes, or groups contain an element, so the count schemes have nothing to take
        from such a file, and they are refused for it.

        Returns
        =======
        Tuple[str, List[Tuple[str, float]], Union[List[Tuple[str]], None], matplotlib.colors.Colormap]
            The scheme ('by_count'/'by_count_continuous'/'by_membership'), the (color_hexcode,
            priority) pairs in assignment order, the list of combinations (for by-membership) or
            None (for the count schemes), and the trimmed colormap the colors were sampled from,
            which a continuous count scale's colorbar spans (None when the colors came from a colors
            file rather than a colormap, which only by-membership allows).
        """
        colormap = layer.get('colormap', True)
        colormap_scheme = layer.get('colormap_scheme')
        membership_colors = layer.get('membership_colors')
        # A message that tells the reader how to ask for a scheme has to name the option this
        # layer's input actually takes, since the option the other input takes is refused here.
        scheme_options = layer.get('scheme_options', PRESENCE_SCHEME_OPTIONS)
        reverse_overlay = layer.get('reverse_overlay', False)
        # Only the count schemes have a scale that can stop early. Coloring by membership needs a
        # color for every combination of the samples, contigs databases, genomes, or groups. This
        # holds however few of them the data reaches.
        if count_scale_top is None:
            count_scale_top = len(individuals)

        if membership_colors is not None:
            return self._membership_layer_given_colors(layer, individuals, count_scale_top)

        # A layer can set the scheme used where none was asked for ('default_colormap_scheme'). It
        # is still a default, so a count that cannot be drawn in bands falls back to a gradient, as
        # it does below.
        if colormap_scheme is not None:
            scheme = colormap_scheme
        elif layer.get('default_colormap_scheme') is not None:
            scheme = layer['default_colormap_scheme']
        else:
            scheme = 'by_membership' if len(individuals) < 4 else 'by_count'

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

        # Coloring by membership needs a color for each combination of the samples, contigs
        # databases, genomes, or groups. The number of combinations doubles with each one added. A
        # colormap holds at most 'cmap.N' colors, so the check below would fail anyway. Making the
        # check before the combinations are enumerated means a large number of them is refused with
        # the same error before gigabytes of memory are spent.
        if scheme == 'by_membership' and 2 ** len(individuals) - 1 > cmap.N:
            self.progress.end()
            # The number is written out in full up to a million and in scientific notation above it.
            # Past a few dozen samples, contigs databases, genomes, or groups, the exact number runs
            # to hundreds of digits that say nothing.
            combos = 2 ** len(individuals) - 1
            combos_shown = f'{combos}' if combos <= 10 ** 6 else f'{combos:.1e}'
            raise ConfigError(
                f"Coloring the {layer['element_type']} layer by membership needs a distinct color "
                f"for every combination of the samples, contigs databases, genomes, or groups. "
                f"There are {len(individuals)} of them. That makes {combos_shown} combinations. "
                f"The layer's colormap holds only {cmap.N}. Color by count instead with "
                f"'{scheme_options['by_count']}', which needs just {count_scale_top} colors, or "
                f"give a colormap with more colors. Note that no color scale can distinguish "
                f"combinations of more than a handful of samples, contigs databases, genomes, or "
                f"groups."
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

        membership_combos = None
        if scheme == 'by_membership':
            membership_combos = []
            for combo_size in range(1, len(individuals) + 1):
                membership_combos += list(combinations(individuals, combo_size))
            if qualitative:
                sample_points = range(len(membership_combos))
            else:
                sample_points = np.linspace(0, 1, len(membership_combos))
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
        needed = len(membership_combos) if scheme == 'by_membership' else count_scale_top
        distinct = len({color for color, _ in color_priorities})
        # Two things stop a band per count from being drawn: a colormap with too few distinct colors
        # to tell the bands apart, and more bands than a colorbar can label in type large enough to
        # read. Either is reason for the gradient, and the message says which it was.
        short_colors = scheme == 'by_count' and distinct != needed
        too_many_bands = scheme == 'by_count' and count_scale_top > MAX_DISCRETE_COUNT_BANDS
        if (short_colors or too_many_bands) and colormap_scheme is None:
            # Nobody asked for the discrete bands. The scheme was chosen from the number of samples,
            # contigs databases, genomes, or groups. There are more of them than the bands can show,
            # so the count is drawn as a gradient instead. The colors are resampled from fractions
            # of the colormap's range, which is what that gradient spans.
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
                    f"samples, contigs databases, genomes, or groups. There are {len(individuals)} "
                    f"of them. That makes {needed} combinations. The colormap supplied only "
                    f"{distinct}. Color by count instead with '{scheme_options['by_count']}', "
                    f"which needs just {count_scale_top} colors, or give a colormap with more "
                    f"distinct colors."
                )
            else:
                advice = (
                    f"Coloring by count in discrete bands needs a distinct color for each count up "
                    f"to {needed}, and the colormap supplied only {distinct}. Color by count on a "
                    f"continuous scale instead. It needs no distinct color per count. Ask for it "
                    f"with '{scheme_options['by_count_continuous']}'. Alternatively, give a "
                    f"colormap with more distinct colors. Or reduce the number of samples, contigs "
                    f"databases, genomes, or groups. For example, group them."
                )
            raise ConfigError(
                f"The colors of the {layer['element_type']} layer of this map could not be "
                f"assigned: {advice} Note that a color scale can hold at most a few hundred "
                f"distinguishable colors in any case, so a very large number of samples, contigs "
                f"databases, genomes, or groups cannot be told apart by color, even where it can "
                f"be drawn."
            )

        return scheme, color_priorities, membership_combos, cmap

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
            'reverse_overlay', and the colors from a kegg-membership-colors-txt file, if any, as
            'membership_colors'/'membership_combo_colors' with the 'membership_colors_flag' they
            came from.

        all_sources : List[str]
            Names of all sources (samples, contigs databases, or genomes), in color-assignment
            order.

        source_type : Literal['contigs database', 'pangenome', 'sample']
            The kind of source, selecting terminology for messages, colorbar labels, and the label
            of the grid's 'unified' panel.

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
            individuals=list(group_sources) if grouped else list(all_sources),
            # The individuals are groups when the sources are grouped, so the noun that names one
            # must then name a group. It labels the messages about each one. It also labels the
            # subdirectory each one is drawn into.
            individual_noun=group_phrase if grouped else singular,
            # The grid labels its per-group colorbars and its 'unified' panel by the kind of source,
            # which stays the same when the sources are grouped.
            grid_source_type=source_type,
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

    def _slice_cap(
        self,
        pathway: kgml.Pathway,
        entry: kgml.Entry,
        count: int,
        ringed: bool = False
    ) -> Union[str, None]:
        """
        Say whether the slices of a map element are too many to draw, and why.

        A box is drawn whole where a slice of equal width would be narrower than
        'MIN_BOX_SLICE_WIDTH'. In a layer with a border column, the columns of a box fit inside its
        ring, so the cap measures the width inside the ring ('kgml.get_box_fill_width'). This holds
        for every box of the layer, with or without a ring on the map, so a box is sliced the same
        way on every map. On global maps, compound circles cover the ends of almost every line, and
        the slices of a line fill the part that the circles leave visible ('kgml.SlicedKGMLCanvas').
        A line is drawn whole where circles hide all of one of its pieces, or where the visible part
        of a piece would give a slice of equal length shorter than the line is wide. The cap counts
        slices only, so slices with small weights can still be thin.

        Parameters
        ==========
        pathway : kgml.Pathway
            The map.

        entry : kgml.Entry
            An ortholog Entry of the map.

        count : int
            The number of slices.

        ringed : bool, False
            True for a layer with a border column.

        Returns
        =======
        Union[str, None]
            'over the cap' for a box, 'line too short' for a line, or None where the slices fit.
        """
        for index, uuid in enumerate(entry.children['graphics']):
            graphics: kgml.Graphics = pathway.uuid_element_lookup[uuid]
            if graphics.type == 'rectangle':
                width = kgml.get_box_fill_width(graphics.width) if ringed else graphics.width
                if width / count < MIN_BOX_SLICE_WIDTH:
                    return 'over the cap'
            elif graphics.type == 'line':
                # Lines are drawn 5 points wide on overview and standard maps
                # ('_stage_element_color'). A line on a global map keeps the width of the KGML file.
                # Biopython draws a line without a width 1 point wide.
                if pathway.is_global_map:
                    width = 1.0 if graphics.width is None else graphics.width
                else:
                    width = 5.0
                visible = self._get_visible_line_length(pathway, entry, index, width)
                if visible <= 0 or visible / count < width:
                    return 'line too short'
        return None

    def _get_visible_line_length(
        self,
        pathway: kgml.Pathway,
        entry: kgml.Entry,
        index: int,
        width: float
    ) -> float:
        """
        Get the length of the part of a line that compound circles leave visible.

        The cap on slices ('_slice_cap') checks a line for each drawn map, in the range-finding pass
        and again in staging. The length is therefore kept in 'visible_line_lengths' the first time
        it is found. The width of a line is set by its map, so the key leaves it out. The compound
        circles of each map are sorted into a grid once, and the grid is kept in
        'compound_circle_grids'.

        Parameters
        ==========
        pathway : kgml.Pathway
            The map.

        entry : kgml.Entry
            An ortholog Entry of the map.

        index : int
            The position of the line among the Graphics elements of the Entry.

        width : float
            The width at which the line is drawn, in points.

        Returns
        =======
        float
            The visible length in points ('kgml.get_visible_span'). It is 0 or less where circles
            hide the whole line, and 0 for a line without points.
        """
        key = (pathway.number, entry.id, index)
        try:
            return self.visible_line_lengths[key]
        except KeyError:
            pass
        graphics: kgml.Graphics = pathway.uuid_element_lookup[entry.children['graphics'][index]]
        coords = graphics.coords or ()
        if coords:
            try:
                grid = self.compound_circle_grids[pathway.number]
            except KeyError:
                grid = kgml.get_compound_circle_grid(pathway)
                self.compound_circle_grids[pathway.number] = grid
            points = list(zip(coords[0::2], coords[1::2]))
            _, visible = kgml.get_visible_span(points, grid, width)
        else:
            visible = 0.0
        self.visible_line_lengths[key] = visible
        return visible

    @staticmethod
    def _slice_label_point(pathway: kgml.Pathway, entry: kgml.Entry) -> Tuple[str, str, str]:
        """
        Find where a map element is, and its label, for the key table of slices.

        The point of a box is its center. The point of a line is where Biopython puts its label, the
        midpoint of the two coordinates at the middle of the list. With an odd number of
        coordinates, these are the middle one and the one after it. The first Graphics element of
        the Entry is read.

        Parameters
        ==========
        pathway : kgml.Pathway
            The map.

        entry : kgml.Entry
            An ortholog Entry of the map.

        Returns
        =======
        Tuple[str, str, str]
            The x and y coordinates and the label, as text. Each is blank where the Entry lacks it.
        """
        if not entry.children['graphics']:
            return '', '', ''
        graphics: kgml.Graphics = pathway.uuid_element_lookup[entry.children['graphics'][0]]
        label = graphics.name or ''
        if graphics.type == 'line':
            if not graphics.coords:
                return '', '', label
            points = list(zip(graphics.coords[0::2], graphics.coords[1::2]))
            middle = len(points) * 0.5
            if int(middle) != middle:
                first, second = int(middle - 0.5), int(middle + 0.5)
            else:
                first, second = int(middle - 1), int(middle)
            x = 0.5 * (points[first][0] + points[second][0])
            y = 0.5 * (points[first][1] + points[second][1])
        else:
            x, y = graphics.x, graphics.y
        return f'{x:g}', f'{y:g}', label

    @staticmethod
    def _slice_unit_labels(unit: Union[str, Tuple[str, Tuple[str, ...]]]) -> Tuple[str, str]:
        """
        Name what a slice stands for, in the 'accession' and 'gene_id' columns of the key table.

        Parameters
        ==========
        unit : Union[str, Tuple[str, Tuple[str, ...]]]
            An accession, or a gene and its accessions in the element ('_slice_values').

        Returns
        =======
        Tuple[str, str]
            The accession and an empty gene ID for a slice of an accession. The accessions of the
            gene in KGML order, joined by commas, and the gene ID for a slice of a gene.
        """
        if isinstance(unit, tuple):
            gene, accessions = unit
            return ','.join(accessions), gene
        return unit, ''

    def _stage_entry_slices(
        self,
        pathway_number: str,
        pathway: kgml.Pathway,
        entry: kgml.Entry,
        spec: dict,
        color_priority: dict,
        derives_compounds: bool,
        slice_records: Union[list, None] = None
    ) -> bool:
        """
        Stage an Entry of a sliced layer in slices, whole, or not at all.

        The spec's slicer finds the slices of the Entry on the map ('_slicer'). An Entry with one
        slice is drawn whole in the color of that slice. An Entry with more slices than fit
        ('_slice_cap') is drawn whole in the color of its stand-in value. An Entry whose stand-in is
        undefined is left uncolored, and so is an Entry whose slices all have a share of 0.

        Otherwise the Entry is drawn in slices. Its Graphics keep one pair of colors, those of the
        slice with the highest priority, which sets the order in which the Entry is drawn. The
        colors of every slice are registered in 'color_priority', in the form of a whole Entry
        ('_stage_element_color'). The top slice is registered last, so its color keeps its priority
        where another slice of the Entry has the same color. Each box then gets the slices as the
        cells of its glyph ('Graphics.glyph'), and each line in its 'slices' attribute, as (color,
        share) pairs in the order of the slices. On a map that derives the colors of compounds from
        reactions, the stand-in casts the Entry's one vote ('compound_vote'). Where the stand-in is
        undefined, the vote is an empty tuple, and the Entry casts no vote.

        '_find_element_entries' can find an Entry more than once. Each staging therefore overwrites
        the glyph or the slices and the vote of the last, and an Entry staged whole has none of
        them. A vote of None means that none was set.

        Parameters
        ==========
        pathway_number : str
            Numeric ID of the map.

        pathway : kgml.Pathway
            The map, whose matched entries are given their colors.

        entry : kgml.Entry
            An ortholog Entry that the layer matches.

        spec : dict
            The layer spec, with the slicing keys that '_draw_map_elements' describes.

        color_priority : dict
            The accumulating color priority dictionary of the map ('_stage_element_color').

        derives_compounds : bool
            True if the map colors compounds from the reactions that touch them.

        slice_records : Union[list, None], None
            Where to record the rows of the key table, and an Entry left uncolored by shares of 0.
            Each record is a tuple of the dictionary of the layer that keeps it ('slice_key' or
            'slice_uncolored'), the key of the Entry on this map, and the rows of the Entry or the
            accessions of an uncolored Entry. Nothing is recorded with None.

        Returns
        =======
        bool
            True if the Entry was colored.
        """
        graphics_elements: List[kgml.Graphics] = [
            pathway.uuid_element_lookup[uuid] for uuid in entry.children['graphics']
        ]
        for graphics in graphics_elements:
            graphics.slices = None
            graphics.glyph = None
            graphics.border = None
        entry.compound_vote = None
        sliced = spec['slicer'](entry)
        if sliced is None:
            return False
        slices = sliced['slices']
        key = (
            pathway_number, spec['slice_rank'], int(entry.id) if entry.id.isdigit() else -1,
            spec['slice_map'], entry.id
        )
        if not slices:
            if slice_records is not None:
                slice_records.append((spec['slice_uncolored'], key, sliced['accessions']))
            return False

        element_type = spec['element_type']
        keep_highest_priority = spec.get('keep_highest_priority', False)
        reason = self._slice_cap(
            pathway, entry, len(slices), ringed=spec.get('slice_ringed', False)
        ) if len(slices) > 1 else None
        if len(slices) > 1 and reason is None:
            top = max(range(len(slices)), key=lambda i: (slices[i][4], -i))
            for i in [i for i in range(len(slices)) if i != top] + [top]:
                self._stage_element_color(
                    pathway, entry, element_type, slices[i][3], slices[i][4], color_priority,
                    keep_highest_priority=keep_highest_priority
                )
            for graphics in graphics_elements:
                pairs = [(color, share) for _, _, share, color, _ in slices]
                if graphics.type == 'rectangle':
                    graphics.glyph = {'cells': pairs}
                else:
                    graphics.slices = pairs
            if derives_compounds:
                # An Entry whose stand-in is undefined has no value to vote with. The empty tuple
                # keeps 'kgml' from taking the color of the top slice as its vote.
                standin = sliced['standin']()
                entry.compound_vote = () if standin is None else (standin[1], standin[2])
            colors = [color for _, _, _, color, _ in slices]
        else:
            if reason is None:
                color, priority = slices[0][3], slices[0][4]
            else:
                standin = sliced['standin']()
                if standin is None:
                    return False
                _, color, priority = standin
            self._stage_element_color(
                pathway, entry, element_type, color, priority, color_priority,
                keep_highest_priority=keep_highest_priority
            )
            colors = [color] * len(slices)

        if slice_records is not None:
            x, y, label = self._slice_label_point(pathway, entry)
            rows = [
                (
                    pathway_number, spec['slice_map'], entry.id, x, y, label,
                    '' if reason else str(position + 1), *self._slice_unit_labels(unit),
                    value if isinstance(value, str) else str(float(value)), str(float(share)),
                    color, reason or ''
                )
                for position, ((unit, value, share, _, _), color)
                in enumerate(zip(slices, colors))
            ]
            slice_records.append((spec['slice_key'], key, rows))
        return True

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
        None leaves the Entry uncolored. '_quantitative_color_priority' turns the value into a color
        and a priority.
        """
        def colorer(entry: kgml.Entry) -> Union[Tuple[str, float], None]:
            value = entry_value(entry)
            if value is None:
                return None
            return Mapper._quantitative_color_priority(value, norm, cmap, reverse_overlay, center)
        return colorer

    @staticmethod
    def _quantitative_color_priority(
        value: float,
        norm: Union[mcolors.Normalize, None],
        cmap: mcolors.Colormap,
        reverse_overlay: bool,
        center: Union[float, None] = None
    ) -> Tuple[str, float]:
        """
        Find the color and priority of a continuous value on a color scale.

        The norm rescales the value to a fraction of the scale, and 'clip=True' keeps it in [0, 1].
        The color is 'cmap' sampled at that fraction. The priority is the fraction, or one minus it
        under 'reverse_overlay'. A degenerate range has no norm, so every value takes the top of the
        colormap. On a centered scale, the one value of such a range is the center itself, so it
        takes the middle color.

        Parameters
        ==========
        value : float
            The value to color.

        norm : Union[mcolors.Normalize, None]
            The scale's norm, or None for a degenerate range.

        cmap : mcolors.Colormap
            The scale's colormap.

        reverse_overlay : bool
            If True, a lower value gets a higher priority, so it is drawn on top.

        center : Union[float, None], None
            The value at the middle of a centered scale, or None for a scale that is not centered.

        Returns
        =======
        Tuple[str, float]
            The color as a hex code, and the priority.
        """
        if norm is None:
            fraction = 0.5 if center is not None else 1.0
        else:
            fraction = float(norm(value))
        priority = (1.0 - fraction) if reverse_overlay else fraction
        return mcolors.rgb2hex(cmap(fraction)), priority

    @staticmethod
    def _bar_extent(
        value: float,
        norm: Union[mcolors.Normalize, None]
    ) -> Tuple[float, float]:
        """
        Find where a bar starts and ends on a bar scale, as fractions of the length of the scale.

        A bar runs from 0 to its value. The norm rescales both to fractions of the scale, and
        'clip=True' keeps them in [0, 1], so a value past a limit reaches that end of the scale. On
        a scale with values on both sides of 0, the fraction at 0 lies inside the scale, and a bar
        of a negative value ends before it. A degenerate scale has no norm. Its one value is 0,
        which is always on a bar scale, so every bar on it has no length. A degenerate color scale
        takes its top color instead ('_quantitative_color_priority').

        Parameters
        ==========
        value : float
            The value of the bar.

        norm : Union[mcolors.Normalize, None]
            The norm of the bar scale ('_make_quantitative_norm'), or None for a degenerate scale.

        Returns
        =======
        Tuple[float, float]
            The fractions of the scale at 0 and at the value.
        """
        if norm is None:
            return 0.0, 0.0
        return float(norm(0.0)), float(norm(value))

    def _categorical_colorer(
        self,
        entry_value: Callable[[kgml.Entry], Union[str, None]],
        value_color_priorities: Dict[str, Tuple[str, float]]
    ):
        """
        Build a colorer that colors an Entry by its categorical value.

        See '_draw_map_elements' for the colorer contract. 'entry_value' gives an Entry's value on
        this map ('_summarize_entry_individuals'). None leaves the Entry uncolored. Each value has
        one color and one priority ('_build_categorical_txt_model').
        """
        def colorer(entry: kgml.Entry) -> Union[Tuple[str, float], None]:
            value = entry_value(entry)
            if value is None:
                return None
            return value_color_priorities[value]
        return colorer

    def _categorical_compound_chooser(
        self,
        layer: dict
    ) -> Callable[[str, List[Tuple[str, float]]], Union[Tuple[str, float], None]]:
        """
        Build the function that colors a compound by the categorical values of its reactions.

        On a global or overview map where no compound layer colors the fills, each compound takes a
        color from the reactions that touch it ('kgml.Pathway.color_associated_compounds'). Each
        touching reaction gives one vote for its value. The layer's accession rule chooses among the
        votes, as it does among the accessions of one map element ('_choose_categorical_value'). The
        key of a random choice is the compound's name, so a compound with the same votes gets the
        same value wherever it is drawn, at another circle or on another map. The choices are
        counted in the layer's '_choices', once per compound and set of votes.

        Parameters
        ==========
        layer : dict
            A reaction layer model of categorical values ('_build_categorical_txt_model').

        Returns
        =======
        Callable[[str, List[Tuple[str, float]]], Union[Tuple[str, float], None]]
            Takes a compound's name and the color and priority of each reaction touching it. Returns
            the compound's color and priority, or None where the rule gives no value.
        """
        value_color_priorities = layer['value_color_priorities']
        color_values = {
            kgml.canonical_color(color): value
            for value, (color, _) in value_color_priorities.items()
        }
        # kgml chooses once for each circle that reactions touch. A compound comes up again with the
        # same votes at another of its circles, on another map, and where the derived colors of a
        # map are checked before it is drawn ('_check_derived_compound_colors'). Each such choice is
        # counted once ('_report_categorical_choices').
        counted: Set[Tuple[str, Tuple[str, ...]]] = set()

        def choose(compound_name: str, color_priorities: List[Tuple[str, float]]):
            # Every colored reaction of the map has the color of a value. Only this layer colors
            # reactions where compound colors are derived ('_draw_map_elements').
            votes = [color_values[kgml.canonical_color(color)] for color, _ in color_priorities]
            value, how = self._choose_categorical_value(
                votes, layer['accession_rule'], layer['value_order'], layer['seed'],
                ('compound', compound_name)
            )
            key = (compound_name, tuple(sorted(votes)))
            if how is not None and key not in counted:
                counted.add(key)
                layer['_choices']['compound'][how] += 1
            # A rule can give no value. The compound then keeps the color given to unprioritized
            # compounds.
            return None if value is None else value_color_priorities[value]
        return choose

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

    def _membership_finder(
        self,
        membership: Dict[str, List[str]],
        group_sources: Union[Dict[str, List[str]], None],
        group_threshold: Union[float, None],
        use_reaction_attribute: bool
    ):
        """
        Build a function giving the sources or groups a map element counts as being in.

        The function returns None if the element is in none of them. Ungrouped, these are the
        sources containing the element ('_entry_sources'). With 'group_sources' they are the
        qualifying groups instead: those where the proportion of the group's own sources containing
        the element meets 'group_threshold'. Two methods call this one. '_membership_colorer' builds
        the colorer that colors an element by how many sources or groups it is in.
        '_find_observed_counts' finds the highest such count across the drawn maps, so the count a
        color stands for and the count the scale was built for cannot disagree.
        """
        grouped = group_sources is not None
        group_source_count: Dict[str, int] = {}
        source_group: Dict[str, str] = {}
        if grouped:
            for group, sources in group_sources.items():
                group_source_count[group] = len(sources)
                for source in sources:
                    source_group[source] = group

        def find_individuals(entry: kgml.Entry) -> Union[Set[str], None]:
            sources = self._entry_sources(entry, membership, use_reaction_attribute)
            if not sources:
                return None
            if not grouped:
                return sources
            group_counts = {group: 0 for group in group_source_count}
            for source in sources:
                if source in source_group:
                    group_counts[source_group[source]] += 1
            individuals = set()
            for group, count in group_counts.items():
                proportion = count / group_source_count[group]
                if (proportion > 0) if group_threshold == 0 else (proportion >= group_threshold):
                    individuals.add(group)
            return individuals or None
        return find_individuals

    def _membership_colorer(
        self,
        membership: Dict[str, List[str]],
        color_priorities: List[Tuple[str, float]],
        membership_combos: Union[List[Tuple[str]], None],
        group_sources: Union[Dict[str, List[str]], None],
        group_threshold: Union[float, None],
        use_reaction_attribute: bool
    ):
        """
        Build a colorer that colors an Entry by the sources (or groups) containing its accessions.

        See '_draw_map_elements' for the colorer contract. An Entry's containing sources are pooled
        across its accessions via 'membership'; with 'group_sources', the qualifying groups (those
        meeting 'group_threshold' among their sources) are used instead. The color and its priority
        are looked up in 'color_priorities'. With 'membership_combos' None, the lookup is by the
        number of sources or groups containing the Entry. Otherwise, when coloring by membership,
        the lookup is by the position of their exact combination in 'membership_combos'. Two counts
        can share a color, on a continuous count scale or with a cyclic colormap. A map keeps one
        priority per color, so the spec of a presence layer asks '_draw_map_elements' to keep the
        highest priority given to a color on the map ('keep_highest_priority'). An Entry in no
        source, or, when grouped, in no qualifying group, is left uncolored.
        """
        combo_lookup: Dict[Tuple[str], Tuple[str]] = {}
        if membership_combos is not None:
            for combo in membership_combos:
                combo_lookup[tuple(sorted(combo))] = combo
        find_individuals = self._membership_finder(
            membership, group_sources, group_threshold, use_reaction_attribute
        )
        count_scale_top = len(color_priorities)

        def colorer(entry: kgml.Entry) -> Union[Tuple[str, float], None]:
            individuals = find_individuals(entry)
            if individuals is None:
                return None
            if membership_combos is None:
                # A count scale can stop below the highest count there is, since '--count-scale-max'
                # takes a ceiling, so a count above the top of the scale takes the top color, as a
                # value above the top of a value scale does.
                return color_priorities[min(len(individuals), count_scale_top) - 1]
            return color_priorities[
                membership_combos.index(combo_lookup[tuple(sorted(individuals))])
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
        membership, a categorical value, or a single fixed color). All layers are staged onto one
        'Pathway' and applied in a single 'set_color_priority' call (list order sets draw order:
        earlier layers render beneath later ones), so one map can mix a quantitative layer with a
        presence layer. The map is drawn once.

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
              accessions. A color keeps the priority of the last entry given it, so if the colorer
              can give a color more than one priority, this order must not change between runs.
            - 'colorer': a function that takes an Entry. It returns a '(color_hexcode, priority)'
              pair, or None to leave the Entry uncolored. The colorer of a layer in the fill mode
              'none' is None. Such a layer stages no fill, and its entries do not count as matching
              the map.
            - 'derived_compound': how a reaction layer colors compounds by the reactions that touch
              them. It is used only on a global or overview map without a compound layer that stages
              fills ('_derives_compounds'). It is a '(mode, colormap)' pair, such as ('average',
              cmap), ('circular_average', cmap) or ('high', None). The mode is passed to
              'kgml.Pathway.set_color_priority' as 'color_associated_compounds'. A compound layer
              has None here; it is not read.

            A spec can also have these keys:
            - 'pathway_numbers': the maps on which the layer has a value to color. On any other map
              the layer colors nothing. It then does not count as matching the map's accessions. A
              summary of samples or groups can be undefined for every element of a map. This key
              lets such a map be skipped.
            - 'require_color': True to count the layer as matching the map only where its colorer
              colors an entry. A rule for categorical values can leave every entry of a map without
              a value. Such a map would otherwise be drawn with nothing colored.
            - 'keep_highest_priority': True to give each color the highest priority of the entries
              given that color on the map. Otherwise a color keeps the priority of the last entry
              given it. A presence layer sets this. Two of its counts can share a color. Its entries
              are colored in the order of its input, such as the rows of a text file. Without this
              key, that order could change the map.

            A sliced layer has these keys too ('_stage_entry_slices'):
            - 'slicer': a function used in place of 'colorer'. It takes an Entry. It returns None
              where no accession or gene of the Entry has a value on the map. Otherwise it returns a
              dictionary. Its 'accessions' are those of the Entry. Its 'slices' are the slices to
              draw, as (unit, value, share, color_hexcode, priority) tuples in the order of
              '_slice_result', or an empty list where every share is 0. The unit is an accession, or
              a gene and its accessions in the Entry. Its 'standin' is a function that takes nothing
              and returns the (value, color_hexcode, priority) of the stand-in value, or None where
              that is undefined ('_slicer').
            - 'slice_map': what the key table calls the map, 'unified' or the name of a sample or
              group.
            - 'slice_rank': where the map comes in the key table, 0 for the 'unified' map.
            - 'slice_key' and 'slice_uncolored': where the rows of the key table and the entries
              left uncolored by shares of 0 are kept, once the map is drawn. They belong to the
              layer.
            - 'slice_ringed': True for a layer with a border column. The cap on slices then
              measures the width inside the ring of a box ('_slice_cap').

            A layer whose border shows on the map has these keys too ('_with_border',
            '_stage_entry_borders'):
            - 'border_colorer': a function that takes an Entry. It returns the '(color_hexcode,
              priority)' pair of the Entry's border value, or None where the Entry has none. The
              color goes to the ring. The priority is not used.
            - 'border_entry_keys': the accessions with border values. They find the entries that
              take rings.
            - 'border_pathway_numbers': the maps on which the border has a value to color, or None
              for every map. On any other map the layer draws no ring.

            A layer whose border has sliced rings has these keys too ('_stage_entry_ring'):
            - 'border_ringer': a function used in place of 'border_colorer'. It takes the map and an
              Entry, and it returns the ring of the Entry with its colors ('_ringer').
            - 'border_slice_map' and 'border_slice_rank': what the key table of the border calls the
              map and where the map comes in it, as for the slices of the fill.
            - 'border_slice_key', 'border_slice_uncolored', and 'border_slice_unshown': where the
              rows of the key table, the entries left out by shares of 0, and the units that no
              column of the fill shows are kept. They belong to the border.

            A layer whose bars show on the map has these keys too ('_with_bars',
            '_stage_entry_bar'):
            - 'bar_value': a function that takes an Entry. It returns the bar value of the Entry on
              the map, or None where the Entry has none there.
            - 'bar_norm': the norm of the bar scale of the map, or None for a scale of a single
              value ('_bar_extent').
            - 'bar_style': how the bars are drawn in their boxes ('BAR_STYLES').
            - 'bar_map': what the warning about boxes without a bar value, and the key table of
              boxes without labels, call the map, 'unified' or the name of a sample or group.
            - 'bar_rank': where the map comes in the key table of boxes without labels, 0 for the
              'unified' map.
            - 'bar_missing' and 'bar_unlabeled': where the boxes left uncolored for want of a bar
              value, and the rows of the key table of boxes drawn without their labels, are kept,
              once the map is drawn. They belong to the model of the bars.

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

        slice_records: list = []
        color_priority, found_entries, _ = self._stage_map_elements(
            pathway_number, pathway, layer_specs, slice_records=slice_records
        )

        if not found_entries and not draw_map_lacking_data:
            return False

        # The slices of a map, and the slices of its rings, are kept for the key tables only once
        # the map is to be drawn.
        for kept, key, value in slice_records:
            kept[key] = value

        compound_layer = next(
            (spec for spec in layer_specs if spec['element_type'] == 'compound'), None
        )
        if compound_layer is not None:
            self._warn_unrenderable_compounds(pathway, compound_layer['entry_keys'])

        # On a global or overview map where a reaction layer stages fills and no compound layer
        # does, compounds are derived from the reactions that touch them, per the reaction layer's
        # 'derived_compound' ('_derives_compounds'). Otherwise, colors are applied directly.
        recolor = 'g' if pathway.is_global_map else 'w'
        if not self._derives_compounds(pathway, layer_specs):
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
        layer_specs: List[dict],
        slice_records: Union[list, None] = None
    ) -> Tuple[dict, bool, dict]:
        """
        Stage the colors of every layer on one map, without applying them.

        Each layer's colorer colors the entries the layer matches. A sliced layer stages its entries
        with its slicer instead ('_stage_entry_slices'). Drawing a map ('_draw_map_elements') and
        checking the compound colors a map derives ('_check_derived_compound_colors') both stage
        colors here. The check therefore sees the colors the drawing will use.

        Parameters
        ==========
        pathway_number : str
            Numeric ID of the map.

        pathway : kgml.Pathway
            The map, freshly loaded. Its matched entries are given their colors.

        layer_specs : List[dict]
            One spec per layer, with the keys '_draw_map_elements' describes. Earlier specs are
            drawn beneath later ones.

        slice_records : Union[list, None], None
            Where the slices of the map, and the slices of its rings, are recorded for the key
            tables ('_stage_entry_slices', '_stage_entry_borders'), and the boxes that it leaves
            uncolored for want of a bar value ('_stage_entry_bar'). Only the drawing of a map
            records them. The check of compound colors does not.

        Returns
        =======
        Tuple[dict, bool, dict]
            Three values:
            - The '{entry_type: {graphics_type: {(fg, bg): priority}}}' dictionary, to be passed to
              'kgml.Pathway.set_color_priority'.
            - Whether any layer matched an entry of the map. A layer with 'require_color' matches
              only where it colors an entry. A ring counts as a match.
            - The rings staged on the map, as the ring of each Entry keyed by its ID, in the form of
              'Graphics.border' ('_stage_entry_borders').
        """
        color_priority: dict = {}
        found_entries = False
        borders: dict = {}
        # A global or overview map can color compounds from the reactions that touch them
        # ('_derives_compounds'). A sliced entry then casts the vote of its stand-in value, or no
        # vote where that is undefined.
        derives_compounds = self._derives_compounds(pathway, layer_specs)
        for spec in layer_specs:
            # A layer with a border stages its rings after its fill. The border keeps its own maps,
            # so a map that the fill skips can still show rings.
            if spec.get('pathway_numbers') is not None and (
                pathway_number not in spec['pathway_numbers']
            ):
                if spec.get('border_colorer') is not None and self._stage_entry_borders(
                    pathway_number, pathway, spec, borders, slice_records
                ):
                    found_entries = True
                continue
            # A layer in the fill mode 'none' has no colorer. Its entries are not staged and do not
            # count as data, so it stages only its rings.
            entries = [] if spec['colorer'] is None else self._find_element_entries(
                pathway, spec['use_reaction_attribute'], spec['entry_keys']
            )
            if entries and not spec.get('require_color', False):
                found_entries = True
            for entry in entries:
                if spec.get('slicer') is not None:
                    if self._stage_entry_slices(
                        pathway_number, pathway, entry, spec, color_priority, derives_compounds,
                        slice_records
                    ):
                        found_entries = True
                    continue
                colored = spec['colorer'](entry)
                if colored is None:
                    continue
                found_entries = True
                color_hexcode, priority = colored
                # On a map that shows bars, a box without a bar value is left uncolored.
                if spec.get('bar_value') is not None and not self._stage_entry_bar(
                    pathway_number, pathway, entry, spec, color_hexcode, slice_records
                ):
                    continue
                self._stage_element_color(
                    pathway, entry, spec['element_type'], color_hexcode, priority, color_priority,
                    keep_highest_priority=spec.get('keep_highest_priority', False)
                )
            if spec.get('border_colorer') is not None and self._stage_entry_borders(
                pathway_number, pathway, spec, borders, slice_records
            ):
                found_entries = True
        return color_priority, found_entries, borders

    def _stage_entry_bar(
        self,
        pathway_number: str,
        pathway: kgml.Pathway,
        entry: kgml.Entry,
        spec: dict,
        color_hexcode: str,
        slice_records: Union[list, None] = None
    ) -> bool:
        """
        Stage the bar of an Entry whose fill is colored, on a map that shows bars.

        A bar is drawn in the boxes that take rings, which are reaction boxes on standard maps
        ('_shows_border'). Each such box gets a glyph ('kgml.Graphics.glyph') in the style of the
        bars ('BAR_STYLES'). In the 'fill' style, its one cell is the gray track
        ('BAR_TRACK_COLOR'). A mark over it, in the color of the fill, runs from 0 to the value of
        the bar along the length of the box. In the 'strip' style, its one cell is the whole fill.
        Two marks lie in a band along its bottom ('BAR_STRIP_HEIGHT'): the white track of the strip,
        and the dark bar over it from 0 to the value. An end line is centered on the free end of the
        bar, at 0 and at the full length too ('BAR_END_LINE_COLOR'). The fill color stays in the fg
        and bg colors of the Entry, so color priority, recoloring, and the colors that compounds
        derive from reactions see only the fill.

        Where 0 lies inside the bar scale of the map, the bars grow both ways from it, and a zero
        line is drawn under the end line ('BAR_ZERO_LINE_COLOR'). A box of the 'fill' style then
        drops its label, which is recorded for the key table of unlabeled boxes
        ('_write_label_keys'). A strip lies below the label, so it keeps it.

        In the 'fill' style, a box without a bar value on the map is left uncolored. A whole fill
        would read as a full bar, and a bare track as a bar of 0. The box is recorded for a warning
        ('_warn_missing_bars'). In the 'strip' style, such a box keeps its whole fill without a
        strip, which cannot be read as a value, since a strip of 0 is white. An Entry drawn as a
        line takes no bar and keeps its whole fill.

        Parameters
        ==========
        pathway_number : str
            Numeric ID of the map.

        pathway : kgml.Pathway
            The map.

        entry : kgml.Entry
            An Entry whose fill the spec's colorer has colored.

        spec : dict
            The layer spec, with the keys of bars that '_draw_map_elements' describes.

        color_hexcode : str
            The color of the fill of the Entry.

        slice_records : Union[list, None], None
            Where a box left uncolored or drawn without its label is recorded, as in
            '_stage_entry_slices'. Nothing is recorded with None.

        Returns
        =======
        bool
            False if the Entry is left uncolored. True if its bar was staged, or if it has no shape
            that takes a bar.
        """
        shapes = [pathway.uuid_element_lookup[uuid] for uuid in entry.children['graphics']]
        shapes = [shape for shape in shapes if self._shows_border(pathway, entry, shape)]
        if not shapes:
            return True
        value = spec['bar_value'](entry)
        strip = spec['bar_style'] == 'strip'
        if value is None and strip:
            return True
        # The accessions of the box in the layer's file name it in the warning about missing bar
        # values and in the key table of unlabeled boxes.
        accessions = [
            accession for accession in self._get_entry_kegg_ids(
                entry, spec['use_reaction_attribute']
            ) if accession in spec['entry_keys']
        ]
        if value is None:
            if slice_records is not None:
                slice_records.append((
                    spec['bar_missing'], (pathway_number, spec['bar_map'], entry.id), accessions
                ))
            return False
        zero, end = self._bar_extent(value, spec['bar_norm'])
        start, stop = min(zero, end), max(zero, end)
        if strip:
            band = BAR_STRIP_HEIGHT
            cells = [(color_hexcode, 1.0)]
            marks = [
                (0, 'x', 0.0, 1.0, BAR_STRIP_TRACK_COLOR, band),
                (0, 'x', start, stop, BAR_STRIP_COLOR, band)
            ]
        else:
            band = None
            cells = [(BAR_TRACK_COLOR, 1.0)]
            marks = [(0, 'x', start, stop, color_hexcode, None)]
        signed = 0 < zero < 1
        rules = []
        if signed:
            rules.append((0, 'x', zero, BAR_ZERO_LINE_COLOR, BAR_ZERO_LINE_WIDTH, band))
        rules.append((0, 'x', end, BAR_END_LINE_COLOR, BAR_END_LINE_WIDTH, band))
        glyph = {'cells': cells, 'marks': marks, 'rules': rules}
        if signed and not strip:
            glyph['label'] = False
            if slice_records is not None:
                x, y, label = self._slice_label_point(pathway, entry)
                key = (
                    pathway_number, spec['bar_rank'], int(entry.id) if entry.id.isdigit() else -1,
                    spec['bar_map'], entry.id
                )
                row = (
                    pathway_number, spec['bar_map'], entry.id, x, y, label, ','.join(accessions),
                    str(float(value))
                )
                slice_records.append((spec['bar_unlabeled'], key, row))
        for shape in shapes:
            shape.glyph = glyph
        return True

    @staticmethod
    def _derives_compounds(pathway: kgml.Pathway, layer_specs: List[dict]) -> bool:
        """
        Say whether a map colors its compounds from the reactions that touch them.

        A global or overview map does this where a reaction layer stages fills and no compound layer
        does. A layer in the fill mode 'none' stages none. A compound layer in this mode draws only
        rings, so the compounds keep the fills derived from the reactions
        ('kgml.Pathway.color_associated_compounds'). '_check_derived_compound_colors' makes the same
        test of the layers.

        Parameters
        ==========
        pathway : kgml.Pathway
            The map.

        layer_specs : List[dict]
            One spec per layer, with the keys '_draw_map_elements' describes.

        Returns
        =======
        bool
            True if the compounds of the map take colors from its reactions.
        """
        if not (pathway.is_global_map or pathway.is_overview_map):
            return False
        fills = {spec['element_type'] for spec in layer_specs if spec['colorer'] is not None}
        return 'reaction' in fills and 'compound' not in fills

    def _stage_entry_borders(
        self,
        pathway_number: str,
        pathway: kgml.Pathway,
        spec: dict,
        borders: dict,
        slice_records: Union[list, None] = None
    ) -> bool:
        """
        Stage the rings of a layer's border on one map.

        The border's colorer colors the entries of its accessions. A ring goes on each shape of an
        Entry that shows rings ('_shows_border'). Its color goes to the Python-only
        'Graphics.border', never to the fg and bg colors of the Entry. Color priority, recoloring,
        and the colors that compounds derive from reactions therefore never see a ring. An Entry
        whose fill has no color keeps an unprioritized fill. 'kgml' recolors that fill, and the ring
        is drawn on it. A border with sliced rings stages each ring with its ringer instead
        ('_stage_entry_ring').

        Parameters
        ==========
        pathway_number : str
            Numeric ID of the map.

        pathway : kgml.Pathway
            The map.

        spec : dict
            The layer spec, with the border keys that '_draw_map_elements' describes.

        borders : dict
            The ring of each staged Entry, keyed by its ID, in the form of 'Graphics.border'. It is
            updated here.

        slice_records : Union[list, None], None
            Where the slices of rings are recorded for the key table of the border, as in
            '_stage_entry_slices'. Nothing is recorded with None.

        Returns
        =======
        bool
            True if a ring was staged.
        """
        # Sliced rings are asked on every map, so that the units and rings that are not shown are
        # recorded too. A map without a value to show gets no ring from the ringer.
        maps = spec.get('border_pathway_numbers')
        if maps is not None and pathway_number not in maps and spec.get('border_ringer') is None:
            return False
        staged = False
        for entry in self._find_element_entries(
            pathway, spec['use_reaction_attribute'], spec['border_entry_keys']
        ):
            shapes = [
                pathway.uuid_element_lookup[uuid] for uuid in entry.children['graphics']
            ]
            shapes = [shape for shape in shapes if self._shows_border(pathway, entry, shape)]
            if not shapes:
                continue
            if spec.get('border_ringer') is not None:
                if self._stage_entry_ring(
                    pathway_number, pathway, entry, shapes, spec, borders, slice_records
                ):
                    staged = True
                continue
            colored = spec['border_colorer'](entry)
            if colored is None:
                continue
            color, _ = colored
            for shape in shapes:
                shape.border = color
            borders[entry.id] = color
            staged = True
        return staged

    def _stage_entry_ring(
        self,
        pathway_number: str,
        pathway: kgml.Pathway,
        entry: kgml.Entry,
        shapes: List[kgml.Graphics],
        spec: dict,
        borders: dict,
        slice_records: Union[list, None] = None
    ) -> bool:
        """
        Stage the sliced ring of an Entry on one map.

        The spec's ringer finds the ring ('_ringer'). The bands of a ring lie above and below the
        columns of the fill, which the fill of the layer staged before ('_stage_entry_slices'). The
        segments of a ring, and a whole ring, do not depend on these columns. The rows of the key
        table of the border are recorded, as the rows of slices are ('_stage_entry_slices'). An
        Entry whose ring is left out by shares of 0 is recorded too, and so are the units whose
        border values no column of the fill shows.

        Parameters
        ==========
        pathway_number : str
            Numeric ID of the map.

        pathway : kgml.Pathway
            The map.

        entry : kgml.Entry
            An ortholog Entry with border values.

        shapes : List[kgml.Graphics]
            The shapes of the Entry that show rings ('_shows_border').

        spec : dict
            The layer spec, with the keys of sliced rings that '_draw_map_elements' describes.

        borders : dict
            The ring of each staged Entry, as in '_stage_entry_borders'. It is updated here.

        slice_records : Union[list, None], None
            Where the rows are recorded, as in '_stage_entry_slices'. Nothing is recorded with None.

        Returns
        =======
        bool
            True if a ring was staged.
        """
        ring = spec['border_ringer'](pathway, entry)
        if ring is None:
            return False
        key = (
            pathway_number, spec['border_slice_rank'], int(entry.id) if entry.id.isdigit() else -1,
            spec['border_slice_map'], entry.id
        )
        if slice_records is not None and ring['unshown']:
            slice_records.append((spec['border_slice_unshown'], key, ring['unshown']))
        if ring['kind'] == 'uncolored' and slice_records is not None:
            slice_records.append((spec['border_slice_uncolored'], key, ring['accessions']))
        if ring['kind'] in ('uncolored', 'none'):
            return False
        for shape in shapes:
            shape.border = ring['border']
        borders[entry.id] = ring['border']
        if slice_records is not None:
            x, y, label = self._slice_label_point(pathway, entry)
            rows = [
                (
                    pathway_number, spec['border_slice_map'], entry.id, x, y, label,
                    '' if position is None else str(position), *self._slice_unit_labels(unit),
                    value if isinstance(value, str) else str(float(value)), str(float(share)),
                    color, ring['reason'] or ''
                )
                for (position, unit, value, share), color in zip(ring['slices'], ring['colors'])
            ]
            slice_records.append((spec['border_slice_key'], key, rows))
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

    def _find_maps(
        self,
        check_dirs: List[str],
        patterns: List[str] = None,
        extensions: Dict[str, str] = None
    ) -> List[str]:
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

        extensions : Dict[str, str], None
            The extension of the files written into a check directory, where it is not '.pdf'. Keys
            are check directories, and values are extensions, such as '.html'.
        """
        if extensions is None:
            extensions = {}

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
                    out_path = self._map_path(
                        check_dir, pathway_number, extension=extensions.get(check_dir, '.pdf')
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
        # recolor the circles that prioritized reactions touch, to reflect their colors.
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

        Reaction elements and compound circles are labeled if the map is of a kind in
        'reaction_label_map_kinds' and 'compound_label_map_kinds', respectively.

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

        if pathway.is_global_map:
            map_kind = 'global'
        elif pathway.is_overview_map:
            map_kind = 'overview'
        else:
            map_kind = 'standard'
        # The reaction elements of a KO map are its ortholog entries.
        self.drawer.draw_map(
            pathway,
            out_path,
            label_orthologs=map_kind in self.reaction_label_map_kinds,
            label_compounds=map_kind in self.compound_label_map_kinds
        )

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

    def _collate_maps_by_map(self, output_dir: str, individuals: List[str]) -> int:
        """
        Gather the individual maps into a directory per map.

        The individual maps are those drawn for each sample, contigs database, genome, or group.
        Each map gets a directory. It holds one file for each sample, contigs database, genome, or
        group, named after it, so a single map can be compared across them by stepping through one
        directory with a file browser's preview. There is no need to open one file in the directory
        of each sample, contigs database, genome, or group.

        Every file gathered is a link to a map already drawn under 'individual', so this second
        arrangement costs no disk space of its own.

        Which maps to gather is read from the drawn files themselves. It is not rebuilt from pathway
        numbers, so it cannot disagree with what was drawn. Nor can it lose track of where
        '--categorize-files' put a map. A map's place within the directory of a sample, contigs
        database, genome, or group becomes its place here. This includes any BRITE subdirectories.

        Parameters
        ==========
        output_dir : str
            Path to the output directory holding the 'individual' subdirectory.

        individuals : List[str]
            Names of the samples, contigs databases, genomes, or groups whose maps are gathered.
            Each is a subdirectory of 'individual'.

        Returns
        =======
        int
            The number of maps gathered, which is the number of directories written.
        """
        collated_dir = os.path.join(output_dir, COLLATED_SUBDIR)
        map_dirs: Set[str] = set()
        for individual in individuals:
            self.progress.update(individual)
            individual_output_dir = os.path.join(output_dir, INDIVIDUAL_SUBDIR, individual)
            for dir_path, subdir_names, basenames in os.walk(individual_output_dir):
                # Everything in the directory of a sample, contigs database, genome, or group is one
                # of its maps. The exceptions are the colorbar of its maps and the flat directory of
                # links to those same maps.
                if FLAT_SUBDIR in subdir_names:
                    subdir_names.remove(FLAT_SUBDIR)
                relative_dir = os.path.relpath(dir_path, individual_output_dir)
                for basename in basenames:
                    map_name, extension = os.path.splitext(basename)
                    if extension != '.pdf' or basename == INDIVIDUAL_COLORBAR_BASENAME:
                        continue
                    map_dir = os.path.normpath(os.path.join(collated_dir, relative_dir, map_name))
                    map_dirs.add(map_dir)
                    self._link_map(
                        os.path.join(dir_path, basename),
                        os.path.join(map_dir, f'{individual}{extension}')
                    )
        return len(map_dirs)

    def _draw_map_players(
        self,
        output_dir: str,
        pathway_numbers: List[str],
        individuals: List[str],
        drawn_individual: Dict[str, Dict[str, bool]],
        keys: List[Tuple[str, Union[str, None]]],
        individual_colorbars: Dict[str, str]
    ) -> int:
        """
        Write a player of each map gathered by map.

        The individual maps are those drawn for each sample, contigs database, genome, or group. A
        player is an HTML file that shows a map's individual maps one after another, with controls
        to pause, step, and move to any frame. It is written beside the directory of the map under
        'by_map', so that directory still holds only the maps.

        A map gets a player where this run drew it for at least one sample, contigs database,
        genome, or group. Where one of them lacks the map, its frame shows the bare map, so every
        player of a run has the same frames in the same order.

        The frames are found from the record of the maps drawn in this run, not from the files on
        disk. The record leaves out the blank maps of grids and any maps of an earlier run.
        Collation reads the files on disk, so it also gathers maps left by an earlier run, and those
        get no player.

        Parameters
        ==========
        output_dir : str
            Path to the output directory holding the 'individual' and 'by_map' subdirectories.

        pathway_numbers : List[str]
            Numeric IDs of the maps of the run.

        individuals : List[str]
            Names of the samples, contigs databases, genomes, or groups whose maps are shown, in the
            order of the frames.

        drawn_individual : Dict[str, Dict[str, bool]]
            Whether each map was drawn. Outer keys are the names of the samples, contigs databases,
            genomes, or groups. Inner keys are pathway numbers.

        keys : List[Tuple[str, Union[str, None]]]
            Paths to the colorbars and bar keys that key the individual maps, each with the element
            type of the border or bars it keys, or None for a colorbar of fills. A colorbar of fills
            is shown beside every frame. The key of a border or of bars is shown only in the players
            of maps that can show the border's rings, which are the shapes that bars take too
            ('_map_shows_borders'). Global maps show none, and overview maps show only the rings of
            compounds.

        individual_colorbars : Dict[str, str]
            Paths to colorbars that key the maps of one group only, where groups have their own
            count scales. Keys are group names. Each colorbar is shown beside its group's frames.

        Returns
        =======
        int
            The number of players written.
        """
        collated_dir = os.path.join(output_dir, COLLATED_SUBDIR)
        frame_keys = [individual_colorbars.get(individual) for individual in individuals]
        if not any(frame_keys):
            frame_keys = None
        count = 0
        for pathway_number in pathway_numbers:
            has_map = [
                drawn_individual[individual].get(pathway_number, False)
                for individual in individuals
            ]
            if not any(has_map):
                continue
            self.progress.update(pathway_number)
            map_keys = [
                key_path for key_path, border_element_type in keys
                if border_element_type is None
                or self._map_shows_borders(border_element_type, pathway_number)
            ]
            blank_dir = None
            try:
                if not all(has_map):
                    self.progress.update(f"{pathway_number}: blank map")
                    blank_dir = self._draw_blank_map(pathway_number)
                in_paths: List[str] = []
                notes: List[str] = []
                for individual, has in zip(individuals, has_map):
                    if has:
                        in_paths.append(self._map_path(
                            os.path.join(output_dir, INDIVIDUAL_SUBDIR, individual), pathway_number
                        ))
                        notes.append('')
                    else:
                        in_paths.append(self._map_path(blank_dir, pathway_number))
                        notes.append(PLAYER_BARE_NOTE)
                self.progress.update(f"{pathway_number}: player")
                title = f"ko{pathway_number} {self.pathway_names.get(pathway_number, '')}".strip()
                self.player_drawer.draw(
                    in_paths,
                    self._map_path(collated_dir, pathway_number, extension=PLAYER_EXTENSION),
                    individuals,
                    notes=notes,
                    title=title,
                    keys=map_keys,
                    frame_keys=frame_keys
                )
            finally:
                if blank_dir is not None:
                    shutil.rmtree(blank_dir, ignore_errors=True)
            count += 1
        return count

    def _draw_blank_map(self, pathway_number: str) -> str:
        """
        Draw the bare map of a pathway in a temporary directory.

        The bare map is drawn as the blank maps of grids are drawn ('_draw_map_grids'). Nothing on
        it is highlighted, and the reactions and compounds are recolored as they are wherever a map
        highlights nothing.

        Parameters
        ==========
        pathway_number : str
            ID number of the pathway, which is five digits.

        Returns
        =======
        str
            The temporary directory, which holds the map at '_map_path(directory, pathway_number)'.
            The caller removes it.
        """
        blank_dir = filesnpaths.get_temp_directory_path()
        progress = self.progress
        run = self.run
        self.progress = terminal.Progress(verbose=False)
        self.run = terminal.Run(verbose=False)
        try:
            self._map_kos_fixed_colors(
                [], blank_dir, [pathway_number], draw_maps_lacking_data=True
            )
        except BaseException:
            shutil.rmtree(blank_dir, ignore_errors=True)
            raise
        finally:
            self.progress = progress
            self.run = run
        return blank_dir

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
        draw_individuals: List[str],
        draw_grid_individuals: List[str],
        draw_files_individuals: List[str],
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
        Make map grids from the individual maps of samples, contigs databases, genomes, or groups.
        Each individual map is a cell in the grid. For example, the cells are the maps of genomes
        for the ungrouped 'pangenome' data source type. They are the maps of contigs databases for
        the ungrouped 'contigs database' type. They are the maps of groups for the grouped
        'pangenome' or 'contigs database' types.

        This method picks up in map methods for multiple data sources after unified and individual
        map files are drawn.

        Parameters
        ==========
        pathway_numbers : List[str]
            IDs of pathway maps to draw.

        draw_individuals : List[str]
            All samples, contigs databases, genomes, or groups for which map files are attempted to
            be drawn.

        draw_grid_individuals : List[str]
            All samples, contigs databases, genomes, or groups for which map grids are attempted to
            be drawn.

        draw_files_individuals : List[str]
            All samples, contigs databases, genomes, or groups for which individual map files were
            drawn.

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
            If True, give a blank map to each panel of a grid whose map was not drawn. This covers
            the individual maps and the 'unified' map. The blank map of a pathway is drawn once and
            copied to each panel that needs it. The blank maps are deleted after the grids are
            drawn.

        source_type : Literal['contigs database', 'pangenome', 'sample'], 'unknown'
            The kind of source in the run. It labels the per-group colorbars of source counts. The
            'unified' panel is labeled 'pangenome' for a pangenome and 'all' otherwise. Each engine
            passes its kind of source even for a grouped run, since a group is made of sources of
            that kind.

        include_unified : bool, True
            If True, each grid leads with the 'unified' map. It summarizes every sample, contigs
            database, genome, or group. It is labeled 'all' or 'pangenome'. If False, a grid holds
            the individual maps alone. The run has not drawn the 'unified' map.
        """
        self.progress.new("Drawing map grid")
        self.progress.update("...")

        # Draw empty maps needed to fill in grids.
        paths_to_remove: List[str] = []
        unified_fillers = False
        if check_maps_lacking_kos:
            # Make a new dictionary. Its outer keys are pathway numbers. Its inner dictionaries say
            # whether the map was drawn for each sample, contigs database, genome, or group.
            drawn_pathway_number: Dict[str, Dict[str, bool]] = {}
            for individual, drawn_individual in drawn['individual'].items():
                for pathway_number, drawn_map in drawn_individual.items():
                    try:
                        drawn_pathway_number[pathway_number][individual] = drawn_map
                    except KeyError:
                        drawn_pathway_number[pathway_number] = {individual: drawn_map}

            # Find the panels that need a blank map. A sample, contigs database, genome, or group in
            # the grid needs one where it has no map of a pathway that another one has.
            unified_dir = os.path.join(output_dir, UNIFIED_SUBDIR)
            blank_dirs: Dict[str, List[str]] = {}
            for pathway_number in pathway_numbers:
                drawn_individual = drawn_pathway_number.get(pathway_number, {})
                map_dirs: List[str] = []
                if set(drawn_individual.values()) == set([True, False]):
                    map_dirs = [
                        os.path.join(output_dir, INDIVIDUAL_SUBDIR, individual)
                        for individual in draw_grid_individuals
                        if not drawn_individual.get(individual, False)
                    ]
                # A 'unified' map is not drawn where it has nothing to color. For example, a rule
                # for categorical values can leave every element of the map without a value. A grid
                # with an individual map of the pathway then gets a blank 'unified' panel, as it
                # does for a missing individual panel. Otherwise the grid would not be drawn.
                if include_unified and any(
                    drawn_individual.get(individual, False) for individual in draw_grid_individuals
                ) and not os.path.exists(self._map_path(unified_dir, pathway_number)):
                    map_dirs.append(unified_dir)
                    unified_fillers = True
                if map_dirs:
                    blank_dirs[pathway_number] = map_dirs

            # The blank map of a pathway is the same in every panel, so it is drawn once and copied
            # to the other panels.
            progress = self.progress
            self.progress = terminal.Progress(verbose=False)
            run = self.run
            self.run = terminal.Run(verbose=False)

            for pathway_number, map_dirs in blank_dirs.items():
                progress.update(f"{pathway_number}: blank map")
                self._map_kos_fixed_colors(
                    [], map_dirs[0], [pathway_number], draw_maps_lacking_data=True
                )
                blank_path = self._map_path(map_dirs[0], pathway_number)
                paths_to_remove.append(blank_path)
                if self.categorize_files:
                    paths_to_remove.append(
                        os.path.join(map_dirs[0], FLAT_SUBDIR, os.path.basename(blank_path))
                    )
                for map_dir in map_dirs[1:]:
                    out_path = self._map_path(map_dir, pathway_number)
                    os.makedirs(os.path.dirname(out_path), exist_ok=True)
                    shutil.copyfile(blank_path, out_path)
                    paths_to_remove.append(out_path)

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
            for group in draw_individuals:
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
            # nothing to show. The pathway's data is in a sample, contigs database, genome, or group
            # that this grid does not hold. Fillers even out a grid that has content. They do not
            # make a grid out of nothing.
            if not include_unified and not any(
                drawn['individual'].get(individual, {}).get(pathway_number, False)
                for individual in draw_grid_individuals
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

            for individual in draw_grid_individuals:
                if self.pathway_categorization is None:
                    individual_map_path = os.path.join(
                        output_dir, INDIVIDUAL_SUBDIR, individual, pathway_basename
                    )
                else:
                    out_dir = os.path.join(
                        output_dir, INDIVIDUAL_SUBDIR, individual,
                        *self.pathway_categorization[pathway_number]
                    )
                    individual_map_path = os.path.join(out_dir, pathway_basename)
                if not os.path.exists(individual_map_path):
                    break
                in_paths.append(individual_map_path)
                labels.append(individual)
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

        # Remove the maps and their links that were only needed for map grids.
        for path in paths_to_remove:
            os.remove(path)
        # Remove the directories left empty by the blank 'unified' maps. Where no 'unified' map was
        # drawn, this removes the 'unified' directory itself.
        if unified_fillers:
            for dir_path, _, _ in os.walk(unified_dir, topdown=False):
                if not os.listdir(dir_path):
                    os.rmdir(dir_path)
        for individual in set(draw_individuals).difference(set(draw_files_individuals)):
            shutil.rmtree(os.path.join(output_dir, INDIVIDUAL_SUBDIR, individual))
            drawn['individual'].pop(individual)

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
        for individual in draw_files_individuals:
            individual_output_dir = os.path.join(output_dir, INDIVIDUAL_SUBDIR, individual)
            for dir_path, subdir_names, filenames in os.walk(individual_output_dir, topdown=False):
                if dir_path != individual_output_dir and not os.listdir(dir_path):
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

    def _map_basename(self, pathway_number: str, extension: str = '.pdf') -> str:
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

        extension : str, '.pdf'
            The extension of the file. A player of the map is an '.html' file.

        Returns
        =======
        str
            The file name, e.g. 'ko00010.pdf', or 'ko00010_Glycolysis_Gluconeogenesis.pdf' where
            file names carry pathway names.
        """
        pathway_name = f'_{self._name_pathway(pathway_number)}' if self.name_files else ''
        return f'ko{pathway_number}{pathway_name}{extension}'

    def _map_path(self, map_dir: str, pathway_number: str, extension: str = '.pdf') -> str:
        """
        Find the path of a pathway's file in a directory of maps.

        With 'categorize_files', the file is in the BRITE subdirectories of the map under that
        directory. Otherwise it is in the directory itself.

        Parameters
        ==========
        map_dir : str
            A directory of maps, such as the directory of one sample under 'individual'.

        pathway_number : str
            ID number of the pathway, which is five digits.

        extension : str, '.pdf'
            The extension of the file. A player of the map is an '.html' file.

        Returns
        =======
        str
            The path of the file, which need not exist.
        """
        basename = self._map_basename(pathway_number, extension=extension)
        if self.pathway_categorization is None:
            return os.path.join(map_dir, basename)
        return os.path.join(map_dir, *self.pathway_categorization[pathway_number], basename)

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

    length_figsize : Tuple[float, float]
        Dimensions in inches of the figure of a key of bar lengths ('draw_length'). The box of the
        key has about the shape of a reaction box.
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
        self.length_figsize: Tuple[float, float] = (6, 2.2)

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
            tick_labels = self._limit_tick_labels(
                cb.formatter, ticks, span, clamped_low, clamped_high
            )
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

    @staticmethod
    def _limit_tick_labels(
        formatter: mticker.ScalarFormatter,
        ticks: List[float],
        span: float,
        clamped_low: bool,
        clamped_high: bool
    ) -> List[str]:
        """
        Typeset the ticks of a scale whose limits, center, or start are labeled.

        The ticks are run through the scale's own formatter, so that such a scale is typeset exactly
        like a plain one ('draw_continuous', 'draw_length').

        Parameters
        ==========
        formatter : matplotlib.ticker.ScalarFormatter
            The formatter of the axis of the scale.

        ticks : List[float]
            The ticks, sorted.

        span : float
            The length of the scale.

        clamped_low, clamped_high : bool
            Whether values lie past the limit at the bottom or the top of the scale. The first or
            last label is then marked with 'CLAMPED_MIN_PREFIX' or 'CLAMPED_MAX_PREFIX'.

        Returns
        =======
        List[str]
            The label of each tick.
        """
        formatter.set_locs(ticks)
        if formatter.get_offset():
            # The formatter wants to set a shared offset or multiplier beside the bar, which the
            # fixed labels below would drop, taking the meaning of every tick with it. Spelling the
            # numbers out in full is longer but says what it means without it.
            formatter.set_useOffset(False)
            formatter.set_scientific(False)
            formatter.set_locs(ticks)
        # The formatter's own call prints any value below 1e-8 as zero. Every label of a bar of tiny
        # values would then read zero. Its format string is applied directly instead. The offset and
        # multiplier are both zero here, so nothing else differs. Only a value within rounding error
        # of zero, as tick arithmetic leaves it, is printed as zero.
        tick_labels = [
            formatter.fix_minus(formatter.format % (0 if abs(tick) < 1e-9 * span else tick))
            for tick in ticks
        ]
        if clamped_low:
            tick_labels[0] = f'{CLAMPED_MIN_PREFIX}{tick_labels[0]}'
        if clamped_high:
            tick_labels[-1] = f'{CLAMPED_MAX_PREFIX}{tick_labels[-1]}'
        return tick_labels

    def draw_length(
        self,
        vmin: float,
        vmax: float,
        out_path: str,
        label: str = None,
        clamped_low: bool = False,
        clamped_high: bool = False,
        track_color: str = BAR_TRACK_COLOR
    ) -> None:
        """
        Save a standalone key of a bar scale to a file.

        A bar scale keys lengths, not colors. The key is a box in the color of the track of bars,
        with a ruler under it. The width of the box stands for the scale. A bar starts at 0, which
        is always on the scale ('Mapper._bar_extent'), so 0 is always labeled. Both ends of the
        scale are labeled too, since a full bar stands for its end, whether the data or a limit sets
        it. An end is marked '≤' or '≥' where values lie past a limit, as on a continuous colorbar
        ('draw_continuous'). A scale with values on both sides of 0 marks 0 with a gray line across
        the box, as its bars do ('BAR_ZERO_LINE_COLOR'). A scale whose one value is 0 is drawn as a
        box labeled 0 alone, since its bars have no length.

        Parameters
        ==========
        vmin : float
            The bottom of the scale.

        vmax : float
            The top of the scale.

        out_path : str
            Path to PDF output file.

        label : str, None
            The label of the key.

        clamped_low, clamped_high : bool, False
            As in 'draw_continuous'.

        track_color : str, BAR_TRACK_COLOR
            The color of the box: the gray track of bars, or the white track of strips
            ('BAR_STRIP_TRACK_COLOR').
        """
        fig = Figure(figsize=self.length_figsize)
        ax = fig.subplots()
        # The axes themselves are the box: the color of the track inside, with a black outline as a
        # reaction box has.
        ax.set_facecolor(track_color)
        for spine in ax.spines.values():
            spine.set_color('#000000')
            spine.set_linewidth(1)
        ax.yaxis.set_visible(False)
        # The labels take the size of the labels of colorbars. They are sized before the automatic
        # ticks are found, which are then spaced to fit them.
        ax.tick_params(
            axis='x', labelsize=(
                self.max_tick_fontsize if self.tick_fontsize is None else self.tick_fontsize
            )
        )
        if vmin == vmax:
            ax.set_xlim(vmin, vmin + 1)
            ticks = [vmin]
            tick_labels = self._limit_tick_labels(
                ax.xaxis.get_major_formatter(), ticks, 1, False, False
            )
        else:
            ax.set_xlim(vmin, vmax)
            # The line at 0 is as wide as the outline of the box, so that it shows at the size of
            # the key.
            if vmin < 0 < vmax:
                ax.axvline(0, color=BAR_ZERO_LINE_COLOR, linewidth=1)
            span = vmax - vmin
            # The start of the bars and the ends of the scale are the labels a reader needs, so an
            # automatic tick all but on top of one of them is dropped.
            marked = sorted({0.0, vmin, vmax})
            kept = [
                tick for tick in ax.get_xticks()
                if vmin <= tick <= vmax and all(
                    abs(tick - value) >= MIN_LENGTH_TICK_SEPARATION_FRACTION * span
                    for value in marked
                )
            ]
            ticks = sorted(kept + marked)
            tick_labels = self._limit_tick_labels(
                ax.xaxis.get_major_formatter(), ticks, span, clamped_low, clamped_high
            )
        ax.set_xticks(ticks)
        ax.set_xticklabels(tick_labels)
        if label:
            ax.set_xlabel(label, labelpad=self.labelpad, fontsize=self.label_fontsize)

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


class HTMLPlayerDrawer:
    """
    Writes HTML files that show input PDF files one after another.

    The page is a player. It can be paused, stepped frame by frame, and moved to any frame with a
    slider. It needs no other file and no network connection, so it can be opened from the disk and
    passed around as a single file.

    Each frame is the first page of a PDF file. It is rendered at the resolution of the largest
    image on the page and stored in the page as a lossless WebP image. A KEGG map draws its base
    image at one pixel per point, so the frame keeps the pixels of the base image exactly. Keys,
    such as colorbars, are shown beside the frames.

    Attributes
    ==========
    overwrite_output : bool
        If True, methods in this class overwrite existing output files.

    seconds_per_frame : float
        The number of seconds each frame is first shown during playback. The viewer can change it on
        the page.

    loop : bool
        If True, playback starts again after the last frame. The viewer can turn this on or off on
        the page.

    key_zoom : float
        The zoom at which a key is rendered. A key is shown at one CSS pixel per point, so a zoom of
        2 keeps it sharp on a display with two device pixels per CSS pixel.

    webp_options : dict
        Arguments of Pillow's WebP encoder. 'quality' sets the effort of lossless compression.
    """
    def __init__(self, overwrite_output: bool = FORCE_OVERWRITE) -> None:
        """
        Parameters
        ==========
        overwrite_output : bool, FORCE_OVERWRITE
            If True, methods in this class overwrite existing output files.
        """
        self.overwrite_output = overwrite_output

        self.seconds_per_frame: float = PLAYER_SECONDS_PER_FRAME
        self.loop: bool = True
        self.key_zoom: float = 2.0
        self.webp_options: dict = {'lossless': True, 'quality': 80, 'method': 4}

    def draw(
        self,
        in_paths: List[str],
        out_path: str,
        labels: List[str],
        notes: List[str] = None,
        title: str = '',
        keys: List[str] = None,
        frame_keys: List[Union[str, None]] = None
    ) -> None:
        """
        Write a player.

        The page is written as its frames are rendered, so no more than one frame is held at once.

        Parameters
        ==========
        in_paths : List[str]
            Paths to the PDF file of each frame, in the order of play. The same path may be given
            for several frames, and its image is then stored once. Every page must have the same
            size.

        out_path : str
            Path to the output HTML file.

        labels : List[str]
            The label of each frame. It is shown above the frame with the frame's position.

        notes : List[str], None
            A note shown beside the label of each frame, or an empty string for no note.

        title : str, ''
            The title of the page.

        keys : List[str], None
            Paths to PDF files of keys, such as colorbars, shown beside every frame.

        frame_keys : List[Union[str, None]], None
            Paths to PDF files of keys shown beside one frame each, or None for a frame without a
            key of its own.
        """
        assert len(in_paths) == len(labels) > 0
        assert notes is None or len(notes) == len(in_paths)
        assert frame_keys is None or len(frame_keys) == len(in_paths)
        filesnpaths.is_output_file_writable(out_path, ok_if_exists=self.overwrite_output)

        with open(PLAYER_TEMPLATE_PATH, encoding='utf-8') as template_file:
            template = template_file.read()
        head, rest = template.split('<!--FRAMES-->')
        middle, tail = rest.split('<!--KEYS-->')
        head = head.replace('__TITLE__', html.escape(title, quote=False))
        head = head.replace('__LAST__', str(len(in_paths) - 1))
        head = head.replace('__SECONDS__', f'{self.seconds_per_frame:g}')
        head = head.replace('__LOOP__', 'true' if self.loop else 'false')
        # The first frame is shown before the script runs, as it is by the script itself.
        first_note = notes[0] if notes is not None else ''
        head = head.replace('__LABEL__', html.escape(f"{labels[0]} (1 of {len(labels)})"))
        head = head.replace('__NOTE__', html.escape(first_note))

        frames = [
            {'label': label, 'image': None, 'key': None, 'note': note}
            for label, note in zip(labels, notes if notes is not None else [''] * len(labels))
        ]
        with open(out_path, 'w', encoding='utf-8', newline='\n') as out_file:
            out_file.write(head)
            # An image shown by several frames is stored once. Images are told apart by their bytes,
            # so two PDF files that render the same are stored once too. The label of a frame is
            # text on the page, so an image can be shown under any label. The first image belongs to
            # the first frame, so it is the one image that is not hidden.
            image_indices: Dict[str, int] = {}
            path_indices: Dict[str, int] = {}
            size = None
            for frame, in_path in zip(frames, in_paths):
                if in_path not in path_indices:
                    image = self._render(in_path)
                    if size is None:
                        size = image.size
                    assert image.size == size, "Every frame of a player must have the same size."
                    data = self._encode(image)
                    digest = hashlib.sha256(data).hexdigest()
                    if digest not in image_indices:
                        image_indices[digest] = len(image_indices)
                        self._write_image(
                            out_file, data, image.size, 'frame', 'map',
                            hidden=image_indices[digest] > 0
                        )
                    path_indices[in_path] = image_indices[digest]
                frame['image'] = path_indices[in_path]
            out_file.write(middle)

            for key_path in keys or []:
                image = self._render(key_path, zoom=self.key_zoom)
                self._write_image(
                    out_file, self._encode(image), image.size, 'key', 'key', scale=self.key_zoom,
                    hidden=False
                )
            if frame_keys is not None:
                out_file.write('<div id="frame-keys">\n')
                key_indices: Dict[str, int] = {}
                for frame, key_path in zip(frames, frame_keys):
                    if key_path is None:
                        continue
                    if key_path not in key_indices:
                        key_indices[key_path] = len(key_indices)
                        image = self._render(key_path, zoom=self.key_zoom)
                        self._write_image(
                            out_file, self._encode(image), image.size, 'key', 'key',
                            scale=self.key_zoom, hidden=frame is not frames[0]
                        )
                    frame['key'] = key_indices[key_path]
                out_file.write('</div>')

            # The frames are described in JSON inside a script element. Escaping '<', '>', and '&'
            # keeps a label such as '</script>' from ending the element early.
            frames_json = json.dumps(frames, ensure_ascii=False, separators=(',', ':'))
            for character, escape in (('<', '\\u003c'), ('>', '\\u003e'), ('&', '\\u0026')):
                frames_json = frames_json.replace(character, escape)
            out_file.write(tail.replace('__FRAMES__', frames_json))

    @staticmethod
    def _render(path: str, zoom: float = None) -> Image.Image:
        """
        Render the first page of a PDF file as RGB pixels on white.

        Parameters
        ==========
        path : str
            Path to the PDF file.

        zoom : float, None
            The zoom at which the page is rendered. If None, the page is rendered at the resolution
            of its largest image, so that the image keeps its pixels. MuPDF resamples an image that
            does not start on the edge of a pixel. The page is therefore also shifted by less than a
            pixel, so that the image's corner falls on one.

        Returns
        =======
        PIL.Image.Image
            The rendered page.
        """
        with fitz.open(path) as pdf_doc:
            page = pdf_doc.load_page(0)
            shift_x = shift_y = 0.0
            if zoom is None:
                zoom = 1.0
                infos = [info for info in page.get_images(full=True) if info[2] * info[3] > 0]
                if infos:
                    info = max(infos, key=lambda info: info[2] * info[3])
                    bbox = page.get_image_bbox(info)
                    if bbox.is_valid and not bbox.is_empty and not bbox.is_infinite:
                        zoom = info[2] / bbox.width
                        shift_x = math.ceil(bbox.x0 * zoom) - bbox.x0 * zoom
                        shift_y = math.ceil(bbox.y0 * zoom) - bbox.y0 * zoom
            pixmap = page.get_pixmap(
                matrix=fitz.Matrix(zoom, 0, 0, zoom, shift_x, shift_y), alpha=False
            )
            image = Image.frombytes('RGB', (pixmap.width, pixmap.height), pixmap.samples)
        # MuPDF keeps the decoded image of each closed PDF file in its store, though no later frame
        # can use it. Emptying the store keeps memory from growing with the number of frames.
        fitz.TOOLS.store_shrink(100)
        return image

    def _encode(self, image: Image.Image) -> bytes:
        """
        Encode an image as a lossless WebP image.

        Parameters
        ==========
        image : PIL.Image.Image
            The image to encode.

        Returns
        =======
        bytes
            The bytes of the WebP image.
        """
        buffer = BytesIO()
        image.save(buffer, format='WEBP', **self.webp_options)
        return buffer.getvalue()

    @staticmethod
    def _write_image(
        out_file,
        data: bytes,
        size: Tuple[int, int],
        role: str,
        alt: str,
        scale: float = 1.0,
        hidden: bool = True
    ) -> None:
        """
        Write an image element whose source is the data of a WebP image.

        Parameters
        ==========
        out_file : file object
            The HTML file being written.

        data : bytes
            The bytes of the WebP image.

        size : Tuple[int, int]
            The width and height of the image in pixels.

        role : str
            The class of the image element, 'frame' or 'key'.

        alt : str
            The text that describes the image.

        scale : float, 1.0
            The number of pixels of the image per CSS pixel, so that an image rendered at a zoom of
            2 is shown at the size of the page it was rendered from. Other than 1, it is declared in
            the 'srcset' of the image, which sets the size at which the browser shows the image.

        hidden : bool, True
            If True, the image is hidden until the page shows it.
        """
        width, height = size
        source = 'src' if scale == 1 else 'srcset'
        out_file.write(
            f'<img class="{role}" width="{round(width / scale)}" height="{round(height / scale)}" '
            f'alt="{html.escape(alt)}" decoding="sync"{" hidden" if hidden else ""} '
            f'{source}="data:image/webp;base64,'
        )
        out_file.write(base64.b64encode(data).decode('ascii'))
        out_file.write('">\n' if scale == 1 else f' {scale:g}x">\n')
