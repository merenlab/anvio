#!/usr/bin/env python
"""
Manipulate KEGG KGML files, which store certain KEGG pathway map data and can be used to create
customized map images.

The XMLOps class loads KGML (XML) files into memory in an object-oriented framework, with a Pathway
element object containing all data from a file via subelements. The XMLOps class also converts a
Pathway object back to a properly formatted string that can be written to an XML file.

The KGML framework implemented in the class is based on the schema:
https://www.kegg.jp/kegg/xml/docs/
"""

from __future__ import annotations

import os
import re
import Bio
import math
import uuid
import numpy as np
import xml.etree.ElementTree as ET

from io import StringIO
from argparse import Namespace
from reportlab.lib import colors as reportlab_colors
from matplotlib.colors import Colormap, rgb2hex, to_hex
from typing import Callable, Dict, Iterable, List, Literal, Tuple, Union
from Bio.KEGG.KGML import KGML_parser, KGML_pathway
from Bio.Graphics.KGML_vis import KGMLCanvas, color_to_reportlab, darken

import anvio.terminal as terminal

from anvio.errors import ConfigError
from anvio import FORCE_OVERWRITE, __version__ as VERSION
from anvio.filesnpaths import is_file_exists, is_output_file_writable

from anvio.metabolism.context import KeggContext
from anvio.metabolism.downloads import download_org_pathway_image_files
from anvio.metabolism.constants import GLOBAL_MAP_ID_PATTERN, OVERVIEW_MAP_ID_PATTERN

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2026, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = VERSION
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"
__status__ = "Development"

# The color of the thin lines drawn between the slices of a box. Gray keeps them distinct from the
# black borders of KEGG boxes.
SLICE_DIVIDER_COLOR = '#808080'

# Half the width of the outline of a compound circle, in points. Biopython draws the outline 1 point
# wide, centered on the edge of the fill, so a circle reaches this far past half its width.
CIRCLE_OUTLINE_HALF_WIDTH = 0.5

# The ring that colors the border of a box or a compound circle ('SlicedKGMLCanvas'). A ring lies
# just inside the black outline of its shape, which Biopython draws 1 point wide, centered on the
# edge. The shape therefore keeps its size, and a circle hides as much of a line as it did before. A
# thin black line lies between the ring and the fill. It keeps a ring visible where its color is
# close to the color of the fill. Widths are in points. A circle wider than 'SMALL_CIRCLE_MAX_WIDTH'
# takes the wider ring. A box ring of 2.5 points stays clear of the label at the top of a box.
BOX_RING_WIDTH = 2.5
SMALL_CIRCLE_RING_WIDTH = 2.0
LARGE_CIRCLE_RING_WIDTH = 3.0
SMALL_CIRCLE_MAX_WIDTH = 24
RING_LINE_WIDTH = 0.5
RING_LINE_COLOR = '#000000'
BOX_OUTLINE_HALF_WIDTH = 0.5

# The largest width and height of a box or circle that Drawer leaves unlabeled, in pixels of the
# base image ('Drawer._unlabel_tiny_graphics'). A shape this small cannot be seen, but Biopython
# would still draw its label. KEGG puts boxes of 2 x 2 pixels in the top left corner of a few 2x
# maps, such as 04723, which also has compound circles of 0 x 0 pixels there. anvi'o shrinks
# compound rectangles to 0 x 0 ('keggmapping.Mapper._zero_out_compound_rectangles').
TINY_GRAPHICS_MAX_SIZE = 2

# A distance in points small enough to count as 0 in the walk along a line into compound circles
# ('get_hidden_length'). A point this close to a circle counts as inside it, so rounding errors do
# not stop the walk where two circles meet.
CIRCLE_TOLERANCE = 1e-6

# How far the end of an ortholog line can lie outside a compound circle and still touch it, in
# points ('Pathway.color_associated_compounds'). The lines of global maps run into the circles of
# their compounds, and a few stop up to a point short of the edge. The arrows of overview maps stop
# short of their circles in the 2x KGML files from which anvi'o draws these maps, most by 2 to 12
# points and some by up to 40 points where many arrows meet at one circle. No arrow ends within 56
# points of a circle that it does not touch. In the KGML files of October 2025, any reach from 40 to
# 56 points therefore finds the same circles on overview maps, and any reach from 1 to 4.5 points
# finds the same circles on global maps.
GLOBAL_MAP_LINE_END_REACH = 2.0
OVERVIEW_MAP_LINE_END_REACH = 48.0


def canonical_color(color: Union[str, None]) -> Union[str, None]:
    """
    Return a color in the single form used for storing and comparing colors.

    Hex codes are case-insensitive, and one color can be written in more than one way besides: with
    an alpha channel, in three-digit shorthand, or by name. Colors are canonicalized wherever they
    are stored in the 'color_priority' attribute of a Pathway or compared against it, because two
    spellings of one color must not be taken for two colors. Without this, an entry colored
    '#FFFFFF' would not match a priority registered as '#ffffff': it would be treated as
    unprioritized, and then recolored into the background, rendering it invisible.

    Parameters
    ==========
    color : Union[str, None]
        A color, or None for a Graphics attribute that is not set.

    Returns
    =======
    Union[str, None]
        The color as a lowercase six-digit hex code. None is returned unchanged, as is a string
        Matplotlib does not recognize as a color, so that an unusable value is reported by whatever
        consumes it rather than here.
    """
    # Only hex codes are canonicalized. Matplotlib would also read a bare number as a shade of gray
    # and a word as a named color, but KGML files use neither for colors and do use such values for
    # other purposes -- '0' marks a reaction that is not highlighted -- so anything that is not a
    # hex code is left exactly as it is rather than being reinterpreted as a color.
    if not isinstance(color, str) or not color.startswith('#'):
        return color
    try:
        return to_hex(color)
    except ValueError:
        return color.lower()

def reserved_recolor_colors(
    recolor: Literal['w', 'g'],
    is_overview_map: bool
) -> Dict[str, str]:
    """
    Return the colors that automatic recoloring gives unprioritized entries, by Entry type.

    These colors are reserved: an entry that is colored one of them cannot be told apart from the
    unprioritized entries of the base map, so 'recolor_unprioritized_entries' refuses to recolor
    when a prioritized entry already holds one. Exposing them here lets a caller check its own
    colors against them BEFORE any map is drawn, rather than discovering the clash part way through.

    Parameters
    ==========
    recolor : Literal['w', 'g']
        The automatic recoloring scheme: 'w' erases unprioritized entries into the base map's own
        white or black, and 'g' colors them light gray.

    is_overview_map : bool
        Whether the map is an overview map, whose unidentified reactions are black rather than
        white.

    Returns
    =======
    Dict[str, str]
        Keys are Entry types ('ortholog', 'compound'), values are the reserved color hex codes.
    """
    assert recolor in ('w', 'g')
    if recolor == 'g':
        return {'ortholog': '#E0E0E0', 'compound': '#E0E0E0'}
    return {'ortholog': '#000000' if is_overview_map else '#FFFFFF', 'compound': '#FFFFFF'}

class Element:
    """
    Representation of an XML element from a KGML file.

    Attributes
    ==========
    uuid : str
        Unique ID, which can be used to look up child elements in a Pathway object.
    """
    # Subclass names are the same as the capitalized tag attribute.
    tag: str
    # Element attributes are required or not, according to the KGML schema.
    attribute_required: Dict[str, bool]

    def __init__(self) -> None:
        self.uuid = str(uuid.uuid4())

class Pathway(Element):
    """
    A pathway element is the parent of all other elements in a KGML file.

    Attributes
    ==========
    subelement_tags : Tuple[str]
        Possible child element tags.

    name : str, None
        KEGG ID of the pathway map.

    org : str, None,
        ko/ec/rn/[org prefix] in ID.

    number : str, None
        Map number in ID.

    title : str, None
        Map title.

    image : str, None
        URL of map image file.

    link : str, None
        URL of map information.

    xml_declaration : str, None
        XML declaration line from file metadata. This is the first line of a reference KGML file.

    xml_doctype : str, None
        Doctype line from file metadata. This is the second line of a reference KGML file.

    xml_comment : str, None
        Comment line from file metadata. This is the third line of a reference KGML file.

    children : Dict[str, List[str]]
        Keys are subelement tags, values are lists of subelement UUIDs.

    uuid_element_lookup : Dict[str, Element], {}
        Keys are unique IDs of elements, values are element objects.

    kegg_id_element_lookup : Dict[str, List[Element]] = {}
        Keys are KEGG IDs, values are lists of element objects with the ID in the name attribute.
        An ID is not necessarily unique to an element, and an element can have multiple KEGG IDs.

    is_global_map : bool, None
        True if the pathway map is a global map, as indicated by the map number.

    is_overview_map : bool, None
        True if the pathway map is an overview map, as indicated by the map number.

    color_priority : Dict[str, Dict[str, Dict[Tuple[str, str], float]]]
        This defines the order of entry graphics by foreground and background color. Set this
        attribute with the method, set_color_priority, which is what keeps its color keys canonical
        ('canonical_color'): the lookups against it canonicalize the colors they are given, so keys
        written in another form -- upper case, an alpha channel, a color name — would not be found,
        and the entries holding those colors would lose their priority and be recolored away.

        Outermost dictionary keys are Entry types, any of the possible values of the type attribute
        of the Entry class, e.g., 'ortholog' and 'compound'. Middle dict keys are Graphics types,
        any of the possible values of the type attribute of the Graphics class, e.g., 'rectangle'
        and 'line'. Inner dict keys are length-2 tuples of fgcolor and bgcolor hex codes,
        respectively. Inner dict values are non-negative numbers indicating the priority of Entry
        Graphics with the given foreground and background colors: higher numbers indicate higher
        priority colors.

        The following is an example of a valid color priority dictionary: ortholog entries with a
        white background (#FFFFFF) take precedence over those with a gray background (#EDEDED).
        Orthologs may be drawn as rectangles or lines.
        {
            'ortholog': {
                'rectangle': {
                    ('#000000', '#FFFFFF'): 1.0,
                    ('#000000', '#EDEDED'): 0.0
                },
                'line': {
                    ('#000000', '#FFFFFF'): 1.0,
                    ('#000000', '#EDEDED'): 0.0
                }
            }
        }
    """
    tag = 'pathway'
    attribute_required = {
        'name': True,
        'org': True,
        'number': True,
        'title': False,
        'image': False,
        'link': False
    }
    subelement_tags: Tuple[str] = (
        'entry',
        'relation',
        'reaction'
    )

    def __init__(self) -> None:
        self.name: str = None
        self.org: str = None
        self.number: str = None
        self.title: str = None
        self.image: str = None
        self.link: str = None

        # Store the XML metadata of the KGML file from which the pathway was loaded.
        self.xml_declaration: str = None
        self.xml_doctype: str = None
        self.xml_comment: str = None

        self.children: Dict[str, List[str]] = {tag: [] for tag in self.subelement_tags}

        self.uuid_element_lookup: Dict[str, Element] = {}
        self.kegg_id_element_lookup: Dict[str, List[Element]] = {}

        self._is_global_map: bool = None
        self._is_overview_map: bool = None

        self.color_priority: Dict[str, Dict[str, Dict[Tuple[str, str], float]]] = {}

        super().__init__()

    @property
    def is_global_map(self):
        if self.number is None:
            return None
        return True if re.match(GLOBAL_MAP_ID_PATTERN, self.number) else False

    @property
    def is_overview_map(self):
        if self.number is None:
            return None
        return True if re.match(OVERVIEW_MAP_ID_PATTERN, self.number) else False

    def set_color_priority(
        self,
        new_color_priority: Dict[str, Dict[str, Dict[Tuple[str, str], float]]],
        recolor_unprioritized_entries: Union[str, Dict[str, Tuple[str, str]]] = False,
        color_associated_compounds: Union[
            Literal['high', 'low', 'average', 'circular_average'],
            Callable[[str, List[Tuple[str, float]]], Union[Tuple[str, float], None]]
        ] = None,
        colormap: Colormap = None
    ) -> None:
        """
        Set the color_priority attribute. Entry elements in the children attribute are automatically
        reordered.

        A single Entry (e.g., representing KOs or compounds) can occur multiple times on a map
        (e.g., as different rectangles or circles), and thus have multiple Graphics elements. It is
        required here that Graphics elements of the same type (e.g., rectangle type Graphics or line
        type Graphics) for an Entry must all have the same foreground and background colors if they
        are to be ordered.

        Entries with higher priority fg/bg colors are placed last in the children attribute and in
        KGML files, and they are rendered in the foreground of the map. The lowest priority entries
        are always those without fg/bg colors defined in the color_priority attribute; these entries
        are placed first in the children attribute and KGML files, and they are rendered in the
        background of the map and thus can be overlaid by higher priority entries.

        Parameters
        ==========
        new_color_priority : Dict[str, Dict[str, Dict[Tuple[str, str], float]]]
            This dictionary is the basis of the color_priority attribute.

            It has the same structure as the color_priority attribute. Outermost dict keys are Entry
            types, any of the possible values of the type attribute of the Entry class, e.g.,
            'ortholog' and 'compound'. Middle dict keys are Graphics types, any of the possible
            values of the type attribute of the Graphics class, e.g., 'rectangle' and 'line'. Inner
            dict keys are length-2 tuples of fgcolor and bgcolor hex codes, respectively. Inner dict
            values are non-negative numbers indicating the priority of Entry Graphics with the given
            foreground and background colors: higher numbers indicate higher priority colors.

            What is actually used to set color_priority is a deep copy of the argument in which
            fg/bg color combinations (entries in each inner dict) are reordered by priority value
            ascending, so that the lowest priority colors appear first in each inner dict. The order
            of Entry and Graphics types (outer and middle dict entries) do not not change: so if,
            for example, 'ortholog' appears before 'compound' in the outermost dict keys, then
            ortholog entries will occur before compound entries in the KGML file, and compounds can
            be drawn over orthologs.

        recolor_unprioritized_entries : Union[str, Dict[str, Dict[str, Tuple[str, str]]]], False
            Recolor unprioritized entries, either automatically with a string argument, or with a
            custom dictionary argument for fine-tuning foreground and background colors by Entry
            type. The valid string arguments for automatic recoloring are 'w' and 'g'.

            It is assumed that global maps contain reaction lines and compound circles, so automatic
            recoloring is tailored to ortholog Entry line Graphics and compound Entry circle
            Graphics. 'w' erases unprioritized lines and circles by coloring them entirely white.
            'g' colors them a light gray (#E0E0E0), consistent with other "unidentified" reactions
            in the base map.

            It is assumed that overview maps contain reaction lines (drawn as arrows) and compound
            circles. Unlike global maps, 'w' colors unprioritized arrows black, consistent with
            other "unidentified" reactions in the base map. 'w' colors unprioritized circles white
            (with a black border). 'g' colors arrows and circles light gray.

            It is assumed that standard maps contain ortholog boxes or lines and compound circles.
            'w' colors unprioritized boxes white, with black text; lines black; and circles white.
            'g' colors unprioritized boxes light gray, with black text; lines light gray; and
            circles light gray. Compound rectangles are also found in a small number of KGML files
            (see 00121, 00621, 01052, 01054), and it is best if the size of these is set to zero so
            they don't obscure the chemical structure drawings to which they correspond on maps (see
            `anvio.keggmapping.Mapper._zero_out_compound_rectangles`).

            A custom dictionary argument can be used to set fg/bg colors in detail. Outer dict keys
            are Entry types, e.g., 'ortholog', 'compound'. Inner dict keys are Graphics types, e.g.,
            'rectangle', 'line'. Inner dict values are length-2 tuples of color hex codes for fg and
            bg colors, respectively. This is shown in the following example, which sets the fg
            (text) color of unprioritized ortholog Entry rectangle Graphics to dark gray and the bg
            to light gray; the fg color of unprioritized ortholog Entry line Graphics to black and
            the bg to white; and the fg (border) of unprioritized compounds to black and the bg to
            white.
            {
                'ortholog': {
                    'rectangle': ('#A9A9A9', '#E0E0E0'),
                    'line': ('#000000', '#FFFFFF')
                },
                'compound': {
                    'circle': ('#000000', '#FFFFFF')
                }
            }

        color_associated_compounds : Union[
            Literal['high', 'low', 'average', 'circular_average'],
            Callable[[str, List[Tuple[str, float]]], Union[Tuple[str, float], None]]
        ], None
            Automatically set the background color of compound entries based on the color priority
            of the ortholog entries that touch them ('color_associated_compounds'). By default,
            compounds participating in reactions are circles, and orthologs are lines on
            global/overview maps and boxes or lines on standard maps.

            An argument of 'high' or 'low' sets the compound background color to the bg color of the
            ortholog with the highest or lowest priority fg/bg color combination. 'average' sets the
            bg color to the average bg color of orthologs with prioritized colors -- unprioritized
            orthologs are not taken into account. 'average' should only be used if priority values
            are normalized to the interval [0, 1] and can be converted to a color given by the
            colormap argument. The average priority value of the orthologs with prioritized colors
            is mapped to a bg color for compound Entry circle Graphics.

            'circular_average' is for a colormap with the same color at both ends, such as one for
            clock times. Priorities of 0 and 1 are then the same point on a circle. The compound
            takes the mean direction of the priorities on that circle. So priorities of 0.99 and
            0.01 give 0, not 0.5. Some priorities have no mean direction, such as 0.25 and 0.75. The
            compound then takes the 'average' color.

            A function can be given instead of a name. It chooses the color of each compound
            ('color_associated_compounds').

            Automatically colored compound entries are added to the color_priority attribute, and
            Entry elements in the children attribute are reordered accordingly. Compound entries
            that are already in the color_priority attribute are exempt from recoloring and given
            higher priority than automatically recolored compound entries.

        colormap : matplotlib.colors.Colormap, None
            If 'average' or 'circular_average' is used as the color_associated_compounds argument, a
            colormap must be provided to map averaged priority values on the interval [0, 1] to a
            background color for compound Entry circle Graphics.
        """
        # Check that new_color_priority only contains positive priority values.
        for new_entry_color_priority in new_color_priority.values():
            for new_graphics_color_priority in new_entry_color_priority.values():
                for priority in new_graphics_color_priority.values():
                    assert priority >= 0

        # Make the color_priority attribute dict, reordering colors from lowest to highest priority.
        color_priority = {}
        for entry_type, new_entry_color_priority in new_color_priority.items():
            color_priority[entry_type] = entry_type_color_priority = {}
            for graphics_type, new_graphics_color_priority in new_entry_color_priority.items():
                entry_type_color_priority[graphics_type] = graphics_type_color_priority = {}
                for colors, priority in sorted(
                    new_graphics_color_priority.items(), key=lambda item: item[1]
                ):
                    graphics_type_color_priority[
                        tuple(canonical_color(color) for color in colors)
                    ] = priority
        self.color_priority = color_priority

        # Reorder Entry elements in the children attribute from lowest to highest priority.
        unprioritized_entry_uuids = self.order_entries_by_color_priority()

        if recolor_unprioritized_entries:
            if isinstance(recolor_unprioritized_entries, str):
                assert recolor_unprioritized_entries in ('w', 'g')

                reserved = reserved_recolor_colors(
                    recolor_unprioritized_entries, self.is_overview_map
                )
                self.recolor_unprioritized_ortholog_entries(
                    unprioritized_entry_uuids, reserved['ortholog']
                )
                self.recolor_unprioritized_compound_entries(
                    unprioritized_entry_uuids, reserved['compound']
                )
            else:
                self.recolor_unprioritized_entries(
                    unprioritized_entry_uuids, recolor_unprioritized_entries
                )

        if color_associated_compounds is None:
            return
        self.color_associated_compounds(color_associated_compounds, colormap=colormap)

        # Reorder compound entries in the children attribute according to color priority.
        unprioritized_entry_uuids = self.order_entries_by_color_priority()

        # Recolor compounds that have not been assigned a color priority.
        if recolor_unprioritized_entries:
            if isinstance(recolor_unprioritized_entries, str):
                self.recolor_unprioritized_compound_entries(
                    unprioritized_entry_uuids,
                    reserved_recolor_colors(
                        recolor_unprioritized_entries, self.is_overview_map
                    )['compound']
                )
            else:
                self.recolor_unprioritized_entries(
                    unprioritized_entry_uuids, recolor_unprioritized_entries
                )

    def order_entries_by_color_priority(self) -> List[str]:
        """
        Reorder Entry (e.g., 'ortholog' and 'compound') UUIDs by color priority in the children
        attribute of the Pathway object. This determines how entries are ordered in KGML files and
        rendered in maps.

        Among the entries without a color priority, those with a ring ('Graphics.border') are put
        after the others, which otherwise keep their order. A ring on an element without a fill
        color is then drawn over an element at the same position that shows nothing, such as a box
        of map 00230 that lies on top of another box.

        Returns
        =======
        List[str]
            UUIDs of Entry elements without a color priority.
        """
        # Entries have different types ('ortholog', 'compound', etc.). Group entries into two
        # classes. "Qualifying" entries have types in the color priority outermost dict. Other
        # entries do not have types in the dict and are assigned the lowest nominal priority of
        # -1.0. These entries keep their order, except that those with rings come last.
        def _rings_last(entry_uuids: List[str]) -> List[str]:
            ringed = {
                entry_uuid for entry_uuid in entry_uuids if any(
                    self.uuid_element_lookup[graphics_uuid].border
                    for graphics_uuid in self.uuid_element_lookup[entry_uuid].children['graphics']
                )
            }
            return (
                [entry_uuid for entry_uuid in entry_uuids if entry_uuid not in ringed]
                + [entry_uuid for entry_uuid in entry_uuids if entry_uuid in ringed]
            )

        reordered_entry_uuids: List[str] = []
        unprioritized_entry_uuids: List[str] = []
        qualifying_entry_uuids: Dict[str, List[str]] = {
            entry_type: [] for entry_type in self.color_priority
        }
        for entry_uuid in self.children['entry']:
            entry: Entry = self.uuid_element_lookup[entry_uuid]
            if entry.type in self.color_priority:
                qualifying_entry_uuids[entry.type].append(entry_uuid)
            else:
                reordered_entry_uuids.append(entry_uuid)
                unprioritized_entry_uuids.append(entry_uuid)
        reordered_entry_uuids = _rings_last(reordered_entry_uuids)

        # Sort "qualifying" entries. Loop through each Entry type in the color priority dict.
        for entry_type, entry_type_color_priority in self.color_priority.items():
            # Retrieve each Entry object of the type. Its priority is determined from fg and bg
            # colors.
            type_qualifying_entry_uuids = qualifying_entry_uuids[entry_type]
            type_priority_entry_uuids: Dict[float, List[str]] = {}
            for entry_uuid in type_qualifying_entry_uuids:
                entry: Entry = self.uuid_element_lookup[entry_uuid]

                # Ensure that all of the Entry Graphics elements are of the same type and have the
                # same fg and bg colors.
                graphics_types: List[str] = []
                fgcolors: List[str] = []
                bgcolors: List[str] = []
                for graphics_uuid in entry.children['graphics']:
                    graphics_element: Graphics = self.uuid_element_lookup[graphics_uuid]
                    graphics_types.append(graphics_element.type)
                    fgcolors.append(graphics_element.fgcolor)
                    bgcolors.append(graphics_element.bgcolor)
                if len(set(graphics_types)) != 1:
                    graphics_type_message = ', '.join([f"'{gt}'" for gt in graphics_types])
                    raise AssertionError(
                        f"The Graphics elements for the Entry with UUID '{entry_uuid}' do not "
                        "have the same type, which is required for ordering entries based on "
                        f"color. Graphics have types: {graphics_type_message}"
                    )
                if len(set(fgcolors)) != 1 or len(set(bgcolors)) != 1:
                    raise AssertionError(
                        f"The Graphics elements in the Entry with UUID '{entry_uuid}' do not "
                        "have consistent foreground and background colors, which is required "
                        "for ordering entries based on color."
                    )

                graphics_type = graphics_types[0]
                try:
                    priority = entry_type_color_priority[graphics_type][
                        (canonical_color(fgcolors[0]), canonical_color(bgcolors[0]))
                    ]
                except KeyError:
                    # The Entry does not have colors in the priority dictionary.
                    priority = -1.0

                try:
                    type_priority_entry_uuids[priority].append(entry_uuid)
                except KeyError:
                    type_priority_entry_uuids[priority] = [entry_uuid]

            # Add the reordered UUIDs of the Entry type to the new list of Entry UUIDs and to the
            # dict mapping priority values to UUIDs of entries of all types. The entries of the type
            # without a color priority come first, with those with rings last among them.
            for priority, entry_uuids in sorted(type_priority_entry_uuids.items()):
                reordered_entry_uuids += (
                    _rings_last(entry_uuids) if priority == -1.0 else entry_uuids
                )

            try:
                unprioritized_entry_uuids += type_priority_entry_uuids[-1.0]
            except KeyError:
                pass

        self.children['entry'] = reordered_entry_uuids

        return unprioritized_entry_uuids

    def recolor_unprioritized_ortholog_entries(
        self,
        unprioritized_entry_uuids: List[str],
        color_hex_code: str
    ) -> None:
        """
        Recolor orthologs without a color priority.

        Ortholog entries are expected to have Graphics elements of type 'line' in global and
        overview maps and type 'rectangle' or 'line' in standard maps. The color is applied to the
        foreground of a line, and the background is made white. The color is applied to the
        background of a rectangle, and the foreground (text) is made black.

        Parameters
        ==========
        unprioritized_entry_uuids : List[str]
            List of UUIDs of all entries without a color priority.

        color_hex_code : str
            Hex code of the color for ortholog graphics.
        """
        if self.is_global_map or self.is_overview_map:
            self.recolor_unprioritized_entries(
                unprioritized_entry_uuids, {'ortholog': {'line': (color_hex_code, '#FFFFFF')}}
            )
        else:
            self.recolor_unprioritized_entries(
                unprioritized_entry_uuids,
                {'ortholog': {
                    'rectangle': ('#000000', color_hex_code), 'line': (color_hex_code, '#000000')
                }}
            )

    def recolor_unprioritized_compound_entries(
        self,
        unprioritized_entry_uuids: List[str],
        color_hex_code: str
    ) -> None:
        """
        Recolor compounds without a color priority.

        Compound entries are assumed to have Graphics elements of type 'circle'. In global maps, the
        color is applied to both the background (fill) and foreground (border) of the circle. In
        overview and standard maps, the background is colored, and the foreground is made black.
        Compound rectangles are also found in a small number of KGML files (see 00121, 00621, 01052,
        01054), and it is best if the size of these is set to zero so they don't obscure the
        chemical structure drawings to which they correspond on maps (see
        `anvio.keggmapping.Mapper._zero_out_compound_rectangles`).

        Parameters
        ==========
        unprioritized_entry_uuids : List[str]
            List of UUIDs of all entries without a color priority.

        color_hex_code : str
            Hex code of the color for compound graphics.
        """
        if self.is_global_map:
            self.recolor_unprioritized_entries(
                unprioritized_entry_uuids,
                {'compound': {'circle': (color_hex_code, color_hex_code)}}
            )
        else:
            self.recolor_unprioritized_entries(
                unprioritized_entry_uuids, {'compound': {'circle': ('#000000', color_hex_code)}}
            )

    def recolor_unprioritized_entries(
        self,
        unprioritized_entry_uuids: List[str],
        type_colors: Dict[str, Dict[str, Tuple[str, str]]]
    ) -> None:
        """
        Entries without a color priority are recolored by Entry type.

        Parameters
        ==========
        unprioritized_entry_uuids : List[str]
            List of UUIDs of all entries without a color priority.

        type_colors : Dict[str, Dict[str, Tuple[str, str]]]
            Outer dictionary keys are Entry types, e.g., 'ortholog', 'compound'. Inner dict keys are
            Graphics types, e.g., 'rectangle', 'line', 'circle'. Inner dict values are length-2
            tuples of foreground and background color hex codes, respectively. This is shown in the
            following example, which sets the fg (text) color of unprioritized ortholog Entry
            rectangle Graphics to dark gray and the bg to light gray; the fg color of unprioritized
            ortholog Entry line Graphics to black and the bg to white; and the fg (border) of
            unprioritized compounds to black and the bg to white.
            {
                'ortholog': {
                    'rectangle': ('#A9A9A9', '#E0E0E0'),
                    'line': ('#000000', '#FFFFFF')
                },
                'compound': {
                    'circle': ('#000000', '#FFFFFF')
                }
            }
        """
        # Prevent unprioritized entries from being assigned prioritized colors.
        for entry_type, entry_type_colors in type_colors.items():
            try:
                entry_type_color_priority = self.color_priority[entry_type]
            except KeyError:
                continue
            for graphics_type, graphics_type_colors in entry_type_colors.items():
                try:
                    graphics_type_color_priority = entry_type_color_priority[graphics_type]
                except KeyError:
                    continue
                # The keys of 'color_priority' are canonical ('set_color_priority'), so the recolor
                # pair is canonicalized to match: the reserved colors here are written in upper case
                # while colormap-sampled colors arrive in lower case.
                if tuple(
                    canonical_color(color) for color in graphics_type_colors
                ) in graphics_type_color_priority:
                    raise ConfigError(
                        f"Unprioritized entry graphics cannot be assigned the same combination of "
                        f"foreground and background colors as prioritized entries of the same entry "
                        f"and graphics types. This occurred for '{entry_type}' entries with "
                        f"'{graphics_type}' graphics."
                    )

        for entry_uuid in unprioritized_entry_uuids:
            entry: Entry = self.uuid_element_lookup[entry_uuid]
            try:
                entry_type_colors = type_colors[entry.type]
            except KeyError:
                continue
            for graphics_uuid in entry.children['graphics']:
                graphics: Graphics = self.uuid_element_lookup[graphics_uuid]
                try:
                    fgcolor_hex_code, bgcolor_hex_code = entry_type_colors[graphics.type]
                except KeyError:
                    continue
                graphics.fgcolor = fgcolor_hex_code
                graphics.bgcolor = bgcolor_hex_code

    def color_associated_compounds(
        self,
        transfer: Union[
            Literal['high', 'low', 'average', 'circular_average'],
            Callable[[str, List[Tuple[str, float]]], Union[Tuple[str, float], None]]
        ],
        colormap: Colormap = None
    ) -> None:
        """
        Set the color of compound entries based on the color priority of the ortholog entries that
        touch them.

        Compound entries are assumed to have Graphics elements of type 'circle'. If the map is
        global, color both the background (interior) and foreground (border) of the circle. In
        overview and standard maps, the background is colored, and the foreground is made black.
        Compound rectangles are also found in a small number of KGML files (see 00121, 00621, 01052,
        01054), and it is best if the size of these is set to zero so they don't obscure the
        chemical structure drawings to which they correspond on maps (see
        `anvio.keggmapping.Mapper._zero_out_compound_rectangles`).

        Each ortholog entry casts one vote, a color and priority, for each compound Entry it
        touches. An entry whose 'compound_vote' attribute is set casts that vote, or no vote if it
        is an empty tuple. Any other entry votes with the color of its Graphics elements and the
        priority of that color. It does not vote if the color has no priority.

        On global and overview maps, an ortholog entry touches the circles at the ends of its lines.
        Each end touches the circle whose edge lies nearest it, if the edge lies within reach
        ('GLOBAL_MAP_LINE_END_REACH', 'OVERVIEW_MAP_LINE_END_REACH'). A compound drawn in several
        places has an Entry for each circle, so each circle takes the votes of only the lines that
        end at it. On other maps, where orthologs are mostly boxes, an ortholog entry touches the
        compound entries that its Reaction element names. That is the Reaction element with the ID
        of the ortholog Entry, whose Substrate and Product elements have the IDs of compound
        entries.

        Parameters
        ==========
        transfer : Union[
            Literal['high', 'low', 'average', 'circular_average'],
            Callable[[str, List[Tuple[str, float]]], Union[Tuple[str, float], None]]
        ]
            An argument of 'high' or 'low' sets the compound color to the bg color of the ortholog
            with the highest or lowest priority fg/bg color combination. 'average' sets the compound
            color to the average bg color of orthologs with prioritized colors -- unprioritized
            orthologs are not taken into account. 'average' should only be used if priority values
            are normalized to the interval [0, 1] and can be converted to a color given by the
            colormap argument. The average priority value of the orthologs with prioritized colors
            is mapped to a color for the compound Entry.

            'circular_average' is for a colormap with the same color at both ends, such as one for
            clock times. Priorities of 0 and 1 are then the same point on a circle. The compound
            takes the mean direction of the priorities on that circle. So priorities of 0.99 and
            0.01 give 0, not 0.5. Some priorities have no mean direction, such as 0.25 and 0.75. The
            compound then takes the 'average' color.

            A function chooses the color itself. It is called once for each compound Entry that has
            votes and no color yet. It takes the Entry's name, such as 'cpd:C00031', and the color
            and priority of each ortholog. It returns the color and priority of the compound. It can
            return None instead. The compound then keeps the color given to unprioritized entries.
            The name lets the function make the same choice for a compound wherever it appears with
            the same votes.

            Compound entries that are already in the color_priority attribute are exempt from
            recoloring and given higher priority than recolored compound entries.

        colormap : matplotlib.colors.Colormap, None
            If 'average' or 'circular_average' is used as an argument to transfer, a colormap must
            be provided to map averaged priority values on the interval [0, 1] to a background color
            for compound entries.
        """
        # Make compound entries searchable by ID, which should be a unique pathway element ID.
        id_compound_entry: Dict[str, Entry] = {}
        for entry_uuid in self.children['entry']:
            entry: Entry = self.uuid_element_lookup[entry_uuid]
            if entry.type != 'compound':
                continue
            id_compound_entry[entry.id] = entry

        # On global and overview maps, the ends of ortholog lines touch compound circles. The
        # circles are sorted into a grid, so that the one nearest the end of a line is found fast.
        if self.is_global_map:
            line_end_reach = GLOBAL_MAP_LINE_END_REACH
        elif self.is_overview_map:
            line_end_reach = OVERVIEW_MAP_LINE_END_REACH
        else:
            line_end_reach = None
        if line_end_reach is not None:
            # A circle belongs to more than one compound Entry only where a KGML file draws them in
            # the same place.
            circle_compound_entries: Dict[Tuple[float, float, float], List[Entry]] = {}
            for compound_entry in id_compound_entry.values():
                for graphics_uuid in compound_entry.children['graphics']:
                    graphics: Graphics = self.uuid_element_lookup[graphics_uuid]
                    if graphics.type != 'circle' or graphics.width is None:
                        continue
                    circle = (graphics.x, graphics.y, get_drawn_circle_radius(graphics.width))
                    circle_compound_entries.setdefault(circle, []).append(compound_entry)
            circle_grid = CircleGrid(circle_compound_entries.keys(), reach=line_end_reach)

        # Make Reaction elements searchable by ID, which is the ID of the ortholog Entry they belong
        # to. Several ortholog entries can have Reaction elements of the same name, each naming the
        # compound entries of its own place on the map.
        id_reaction: Dict[str, Reaction] = {}
        for reaction_uuid in self.children['reaction']:
            reaction: Reaction = self.uuid_element_lookup[reaction_uuid]
            id_reaction[reaction.id] = reaction

        # For each compound Entry touched by color-prioritized ortholog entries, record the colors
        # and priorities of these entries.
        compound_uuid_color_priorities: Dict[str, List[Tuple[str, float]]] = {}
        # Loop through each ortholog Entry.
        for entry_uuid in self.children['entry']:
            entry: Entry = self.uuid_element_lookup[entry_uuid]
            if entry.type != 'ortholog':
                continue

            # Ensure that all of the ortholog Entry Graphics elements have the same fg/bg colors.
            graphics_types: List[str] = []
            fgcolors: List[str] = []
            bgcolors: List[str] = []
            line_ends: List[Tuple[float, float]] = []
            for graphics_uuid in entry.children['graphics']:
                graphics: Graphics = self.uuid_element_lookup[graphics_uuid]
                graphics_types.append(graphics.type)
                fgcolors.append(graphics.fgcolor)
                bgcolors.append(graphics.bgcolor)
                if graphics.type == 'line' and graphics.coords:
                    line_ends.append((graphics.coords[0], graphics.coords[1]))
                    line_ends.append((graphics.coords[-2], graphics.coords[-1]))
            if len(set(graphics_types)) != 1:
                graphics_type_message = ', '.join([f"'{gt}'" for gt in graphics_types])
                raise AssertionError(
                    f"The Graphics elements for the Entry with UUID '{entry_uuid}' do not "
                    "have the same type, which is required for ordering entries based on "
                    f"color. Graphics have types: {graphics_type_message}"
                )
            if len(set(fgcolors)) != 1 or len(set(bgcolors)) != 1:
                raise AssertionError(
                    "The Graphics elements in the ortholog Entry with the following UUID do not "
                    "have consistent foreground and background colors, which is required for "
                    f"ordering entries based on color: {entry_uuid}"
                )

            graphics_type = graphics_types[0]
            fgcolor = fgcolors[0]
            bgcolor = bgcolors[0]
            if entry.compound_vote is None:
                try:
                    priority = self.color_priority['ortholog'][graphics_type][
                        (canonical_color(fgcolor), canonical_color(bgcolor))
                    ]
                except KeyError:
                    # Unprioritized ortholog entries do not affect the color of the compounds they
                    # touch.
                    continue
                if graphics_type == 'line':
                    ortholog_color = fgcolor
                else:
                    ortholog_color = bgcolor
            elif entry.compound_vote == ():
                # A sliced entry with no value to stand for its slices casts no vote.
                continue
            else:
                # A sliced entry casts the one vote that was set for it.
                ortholog_color, priority = entry.compound_vote

            # Find the compound entries that the ortholog touches, each once.
            compound_entries: List[Entry] = []
            if line_end_reach is not None and line_ends:
                for x, y in line_ends:
                    circle = circle_grid.get_nearest_circle(x, y)
                    if circle is None:
                        continue
                    for compound_entry in circle_compound_entries[circle]:
                        if compound_entry not in compound_entries:
                            compound_entries.append(compound_entry)
            elif entry.id in id_reaction:
                reaction = id_reaction[entry.id]
                for reactant_uuid in reaction.children['substrate'] + reaction.children['product']:
                    reactant: Union[Substrate, Product] = self.uuid_element_lookup[reactant_uuid]
                    compound_entry = id_compound_entry[reactant.id]
                    if compound_entry not in compound_entries:
                        compound_entries.append(compound_entry)

            for compound_entry in compound_entries:
                try:
                    compound_uuid_color_priorities[compound_entry.uuid].append(
                        (ortholog_color, priority)
                    )
                except KeyError:
                    compound_uuid_color_priorities[compound_entry.uuid] = [
                        (ortholog_color, priority)
                    ]

        # Define functions for finding compound Entry color.
        def _get_high_color(color_priorities: List[Tuple[str, float]]) -> Tuple[str, float]:
            return sorted(color_priorities, key=lambda t: -t[1])[0]

        def _get_low_color(color_priorities: List[Tuple[str, float]]) -> Tuple[str, float]:
            return sorted(color_priorities, key=lambda t: t[1])[0]

        def _get_average_color(color_priorities: List[Tuple[str, float]]) -> Tuple[str, float]:
            priority = np.mean([t[1] for t in color_priorities])
            color = rgb2hex(colormap(priority))
            return color, priority

        def _get_circular_average_color(
            color_priorities: List[Tuple[str, float]]
        ) -> Tuple[str, float]:
            # Each priority is a point on a circle, with 0 and 1 the same point. The means of the
            # sines and cosines of the points give their mean direction.
            priorities = np.array([t[1] for t in color_priorities])
            angles = 2 * np.pi * priorities
            sine = np.mean(np.sin(angles))
            cosine = np.mean(np.cos(angles))
            # Equal priorities, as from a single ortholog, take the plain average. This gives them
            # the same color and priority as 'average'. Priorities that cancel out, such as 0.25 and
            # 0.75, have no mean direction. They take the plain average too.
            if np.all(priorities == priorities[0]) or np.hypot(sine, cosine) < 1e-9:
                return _get_average_color(color_priorities)
            # Rounding removes floating-point error. Without it, a mean on the boundary between two
            # colors of the colormap could take the lower color.
            priority = round(float(np.arctan2(sine, cosine) / (2 * np.pi)), 12) % 1.0
            return rgb2hex(colormap(priority)), priority

        if callable(transfer):
            get_color_priority = None
        elif transfer == 'high':
            get_color_priority = _get_high_color
        elif transfer == 'low':
            get_color_priority = _get_low_color
        elif transfer == 'average':
            get_color_priority = _get_average_color
        elif transfer == 'circular_average':
            get_color_priority = _get_circular_average_color
        else:
            raise AssertionError

        # Set compound Entry color.
        for compound_uuid, color_priorities in compound_uuid_color_priorities.items():
            compound_entry: Entry = self.uuid_element_lookup[compound_uuid]

            # Get all of the Graphics elements for the Entry.
            graphics_elements: List[Graphics] = []
            for graphics_uuid in compound_entry.children['graphics']:
                graphics_elements.append(self.uuid_element_lookup[graphics_uuid])

            set_color = True
            for graphics in graphics_elements:
                try:
                    # The compound Entry has already been assigned a color priority, so don't
                    # recolor it automatically.
                    self.color_priority['compound']['circle'][
                        (canonical_color(graphics.fgcolor), canonical_color(graphics.bgcolor))
                    ]
                    set_color = False
                except KeyError:
                    continue
            if not set_color:
                continue

            if get_color_priority is None:
                chosen = transfer(compound_entry.name, color_priorities)
                if chosen is None:
                    continue
                compound_color, compound_priority = chosen
            else:
                compound_color, compound_priority = get_color_priority(color_priorities)
            # Set the color of each Graphics element.
            for graphics in graphics_elements:
                graphics.bgcolor = compound_color
                if self.is_global_map:
                    graphics.fgcolor = compound_color

            # Record the compound Element color priority.
            try:
                entry_type_color_priority = self.color_priority['compound']
            except KeyError:
                self.color_priority['compound'] = entry_type_color_priority = {}
            try:
                graphics_type_color_priority = entry_type_color_priority['circle']
            except KeyError:
                entry_type_color_priority['circle'] = graphics_type_color_priority = {}
            # Canonical keys, as 'set_color_priority' writes: this dictionary is written here rather
            # than through that method, and the lookups against it canonicalize the colors they are
            # given, so a key in any other form would never be found.
            if self.is_global_map:
                graphics_type_color_priority[
                    (canonical_color(compound_color), canonical_color(compound_color))
                ] = compound_priority
            else:
                graphics_type_color_priority[
                    (canonical_color(graphics.fgcolor), canonical_color(compound_color))
                ] = compound_priority

    def get_entries(
        self,
        entry_type: str = None,
        kegg_ids: Iterable[str] = None,
        expect_kegg_ids: bool = False
    ) -> List[Entry]:
        """
        Get Entry elements from the pathway.

        Parameters
        ==========
        entry_type : str, None
            The type of Entry to return. By default entries of all types are returned. Permitted
            Entry types are given by the types attribute of the Entry class.

            The box Entry types (line Entry types in global and overview maps) are as follows given
            the map name prefix: maps starting with 'ko' have box/line entries of type 'ortholog',
            'ec' have type 'enzyme', 'rn' have type 'reaction', and organism-specific maps with
            <org prefix> have type 'ko'.

        kegg_ids : Iterable[str], None
            If KEGG IDs are provided, then only entries with these IDs in their 'name' attribute are
            sought. With the default argument of None, all entries of the type are returned. KEGG
            IDs should not be the full ID found in the Entry name attribute, but the part after the
            colon. For example, instead of 'ko:K01080', 'cpd:C12144', and 'path:map00604', which is
            how they appear in the KGML file, use 'K01080', 'C12144', and 'map00604'.

        expect_kegg_ids : bool, False
            If KEGG IDs are provided and this argument is True, then an exception is raised if they
            are not found among the entries in the pathway.

        Returns
        =======
        List[Entry]
            A list of entry element objects contained in the pathway.
        """
        if entry_type is not None:
            assert entry_type in Entry.types
        if kegg_ids is not None:
            assert entry_type is None

        entries: List[Entry] = []

        if kegg_ids is None:
            for uuid in self.children['entry']:
                entry: Entry = self.uuid_element_lookup[uuid]
                if entry_type is not None and entry.type != entry_type:
                    continue
                entries.append(entry)
            return entries

        missing_kegg_ids: List[str] = []
        for kegg_id in kegg_ids:
            try:
                elements = self.kegg_id_element_lookup[kegg_id]
            except KeyError:
                missing_kegg_ids.append(kegg_id)
                continue
            for element in elements:
                if element.tag == 'entry':
                    entries.append(element)
        if missing_kegg_ids and expect_kegg_ids:
            raise ValueError(
                "The following 'kegg_ids' that were provided are not found among entries in the "
                f"pathway: {', '.join(missing_kegg_ids)}"
            )
        return entries

    def scale_graphics(self, factor: float, entry_type: str = None):
        """
        Change the scale of entry graphics.

        Rescaling all of the graphics is useful in fitting the KGML file to map images with
        different resolutions. For example, 1x and 2x resolution image files can be downloaded from
        KEGG, but only KGML files fitting the 1x images.

        Parameters
        ==========
        factor : float
            Factor by which to rescale all graphical elements in the pathway.

        entry_type : str, None
            Only rescale graphics for a certain type of entry in the KGML file, such as "ortholog"
            or "compound". The argument must be from the types attribute of the Entry class.
        """
        for entry in self.get_entries(entry_type=entry_type):
            for graphics_uuid in entry.children['graphics']:
                graphics: Graphics = self.uuid_element_lookup[graphics_uuid]
                for attrib in ('x', 'y', 'width', 'height'):
                    value = getattr(graphics, attrib, None)
                    if value is None:
                        continue
                    setattr(graphics, attrib, value * factor)
                value = getattr(graphics, 'coords', None)
                if value is None:
                    continue
                setattr(graphics, 'coords', tuple([coord * factor for coord in value]))

class Entry(Element):
    """
    An entry element contains information about a node of the pathway.

    Attributes
    ==========
    types : Tuple[str]
        Possible entry types.

    subelement_tags : Tuple[str]
        Possible subelement tags.

    id : str, None
        ID unique to map.

    name : str, None
        KEGG ID(s) represented by the element.

    type : str, None
        Entry type.

    reaction : str, None
        KEGG ID(s) of reaction(s) represented by the element.

    link : str, None
        URL of entry information.

    children : Dict[str, List[str]]
        Keys are subelement tags, values are lists of subelement UUIDs.

    compound_vote : Tuple[str, float], Tuple[()], None
        A color and priority that replace those of the Graphics elements of an ortholog entry when
        the entry colors its compounds ('Pathway.color_associated_compounds'). A sliced entry has
        one, so it casts one vote rather than one per slice. An empty tuple means that the entry
        casts no vote. None means that no vote was set, so the entry votes with the colors of its
        Graphics elements. This attribute is not part of KGML, and XMLOps does not write it.
    """
    tag = 'entry'
    attribute_required = {
        'id': True,
        'name': True,
        'type': True,
        'reaction': False,
        'link': False
    }
    types: Tuple[str] = (
        'ortholog',
        'enzyme',
        'reaction',
        'gene',
        'group',
        'compound',
        'map',
        'brite',
        'other'
    )
    subelement_tags: Tuple[str] = (
        'graphics',
        'component'
    )

    def __init__(self) -> None:
        self.id: str = None
        self.name: str = None
        self.type: str = None
        self.reaction: str = None
        self.link: str = None

        self.children: Dict[str, List[str]] = {n: [] for n in self.subelement_tags}

        self.compound_vote: Union[Tuple[str, float], Tuple[()]] = None

        super().__init__()

class Graphics(Element):
    """
    A graphics element contains drawing information on the entry parent element.

    Attributes
    ==========
    types : Tuple[str]
        Possible shapes of graphical objects.

    name : str, None
        Label of graphical object on map.

    fgcolor : str, None
        Foreground color of graphical object on map.

    bgcolor : str, None
        Background color of graphical object on map.

    type : str, None
        Shape of graphical object on map.

    x : float, None
        X axis position of graphical object on map.

    y : float, None
        Y axis position of graphical object on map.

    coords : Tuple[float], None
        Polyline coordinates of "line"-type graphical object on map.

    width : float, None
        Width of graphical object on map.

    height : float, None
        Height of graphical object on map.

    slices : List[Tuple[str, float]], None
        The color and share of each slice of a line, in order from its first coordinate. There are
        at least 2 slices, and the shares sum to 1. Drawer draws them as dashes with
        'SlicedKGMLCanvas'. A rectangle shows its slices in its 'glyph' instead. This attribute is
        not part of KGML, and XMLOps does not write it.

    glyph : Dict[str, Union[str, list]], None
        What a rectangle or a circle shows inside its outline, as a display list. The shapes of the
        glyph are painted first, kind by kind in the order listed here. The ring ('border') comes
        next, then the thin gray lines between columns, the black line inside the ring, the outline,
        and the label. The glyph of a rectangle has these kinds:
        - 'cells': the color and share of each column of the fill, left to right, as (color, share)
          pairs. Each column is the full height of the rectangle. The columns share its width, or
          the width inside its ring ('get_box_fill_width'). The shares sum to 1.
        - 'marks': rectangles over the columns, such as bars, as (cell, axis, start, end, color,
          band) tuples. Each lies in the part of a column that the box shows, and its extent along
          the axis is given in fractions of that part ('SlicedKGMLCanvas._draw_box_marks').
        - 'rules': lines over the marks, such as the end lines of bars, as (cell, axis, position,
          color, line_width, band) tuples in the same part of a column.
        A glyph can also hold 'label': False, which draws the rectangle without its label, where the
        label would hide what the glyph shows. Without this key, the rectangle is labeled as any
        other. A rectangle without a glyph is one column in its bg color. The glyph of a circle,
        such as a pie, has these kinds instead, with positions in fractions of a turn clockwise from
        12 o'clock:
        - 'track': the color of the whole circle, under everything else.
        - 'wedges': sectors over the track, as (start, end, color) tuples. Each runs from the center
          to the outline.
        - 'rules': lines over the wedges along radii, such as the end line of a pie, as (position,
          color, line_width) tuples. Each runs from the center to the edge of the part of the circle
          that it shows ('SlicedKGMLCanvas._draw_circle_glyph').
        Drawer draws glyphs with 'SlicedKGMLCanvas'. This attribute is not part of KGML, and XMLOps
        does not write it.

    border : Union[str, Tuple[str, list]], None
        The ring drawn just inside the outline of a rectangle or a circle ('BOX_RING_WIDTH'). A
        color draws a whole ring. A rectangle can instead take one of these:
        - ('bands', colors): one color per column of the 'cells' of 'glyph', in the same order, or
          None for a column without a band. Each column takes a band above and below it. The first
          and last columns also take the left and right ends of the ring.
        - ('segments', [(color, share), ...]): the ring is divided among at least 2 segments,
          clockwise from the top left corner. The shares sum to 1.
        Drawer draws rings with 'SlicedKGMLCanvas'. This attribute is not part of KGML, and XMLOps
        does not write it.
    """
    tag = 'graphics'
    attribute_required = {
        'name': False,
        'fgcolor': False,
        'bgcolor': False,
        'type': False,
        'x': False,
        'y': False,
        'coords': False,
        'width': False,
        'height': False
    }
    types: Tuple[str] = (
        'rectangle',
        'circle',
        'roundrectangle',
        'line'
    )

    def __init__(self) -> None:
        self.name: str = None
        self.fgcolor: str = None
        self.bgcolor: str = None
        self.type: str = None
        self.x: float = None
        self.y: float = None
        self.coords: Tuple[float] = None
        self.width: float = None
        self.height: float = None

        self.slices: List[Tuple[str, float]] = None
        self.glyph: Dict[str, Union[str, list]] = None
        self.border: Union[str, Tuple[str, list]] = None

        super().__init__()

class Component(Element):
    """
    A component element is only applicable to "group"-type entry elements representing a complex.

    Attributes
    ==========
    id : str, None
        ID of component element unique to map.
    """
    tag = 'component'
    attribute_required = {
        'id': True
    }
    def __init__(self) -> None:
        self.id: str = None

        super().__init__()

class Relation(Element):
    """
    A relation element is an edge between proteins, gene products, compounds, and pathways.

    Attributes
    ==========
    types : Tuple[str]
        Possible types of relations.

    subelement_tags : Tuple[str]
        Possible subelement tags.

    entry1 : str, None
        ID unique to map representing a node in the relationship.

    entry2 : str, None
        ID unique to map representing the other node in the relationship.

    children : Dict[str, List[str]]
        Keys are subelement tags, values are lists of subelement UUIDs.
    """

    tag = 'relation'
    attribute_required = {
        'entry1': True,
        'entry2': True,
        'type': True
    }
    types: Tuple[str] = (
        'ECrel',
        'PPrel',
        'GErel',
        'PCrel',
        'maplink'
    )
    subelement_tags: Tuple[str] = (
        'subtype',
    )
    def __init__(self) -> None:
        self.entry1: str = None
        self.entry2: str = None
        self.type: str = None

        self.children: Dict[str, List[str]] = {n: [] for n in self.subelement_tags}

        super().__init__()

class Subtype(Element):
    """
    A subtype element specifies more detailed information about the relation.

    Attributes
    ==========
    names : Tuple[str]
        Possible names of subcategories of relation.

    name : str, None
        Name of the subcategory of relation.

    value : str, None
        The value represents information on the subcategory relation.
    """
    tag = 'subtype'
    attribute_required = {
        'name': True,
        'value': True
    }
    names: Tuple[str] = (
        'compound',
        'hidden compound',
        'activation',
        'inhibition',
        'expression',
        'repression',
        'indirect effect',
        'state change',
        'binding/association',
        'dissociation',
        'missing information',
        'phosphorylation',
        'dephosphorylation',
        'glycosylation',
        'ubiquitination',
        'methylation'
    )

    def __init__(self) -> None:
        self.name: str = None
        self.value: str = None

        super().__init__()

class Reaction(Element):
    """
    A chemical reaction element.

    Attributes
    ==========
    types : Tuple[str]
        Possible types of reactions.

    subelement_tags : Tuple[str]
        Possible subelement tags.

    id : str, None
        ID of reaction unique to map.

    name : str, None
        KEGG ID(s) represented by the reaction.

    type : str, None
        Reversible vs. irreversible reaction, as drawn on the map.

    children : Dict[str, List[str]]
        Keys are subelement tags, values are lists of subelement UUIDs.
    """
    tag = 'reaction'
    attribute_required = {
        'id': True,
        'name': True,
        'type': True
    }
    types: Tuple[str] = (
        'reversible',
        'irreversible'
    )
    subelement_tags: Tuple[str] = (
        'substrate',
        'product'
    )

    def __init__(self) -> None:
        self.id: str = None
        self.name: str = None
        self.type: str = None

        self.children: Dict[str, List[str]] = {n: [] for n in self.subelement_tags}

        super().__init__()

class Substrate(Element):
    """
    A substrate element represents a substrate node in a parent reaction element.

    Attributes
    ==========
    subelement_tags : Tuple[str]
        Possible subelement tags.

    id : str, None
        ID of substrate unique to map corresponding to a compound entry.

    name : str, None
        KEGG ID of the compound.

    children : Dict[str, List[str]]
        Keys are subelement tags, values are lists of subelement UUIDs.
    """
    tag = 'substrate'
    attribute_required = {
        'id': True,
        'name': True
    }
    subelement_tags: Tuple[str] = (
        'alt',
    )

    def __init__(self) -> None:
        self.id: str = None
        self.name: str = None

        self.children: Dict[str, List[str]] = {n: [] for n in self.subelement_tags}

        super().__init__()

class Product(Element):
    """
    A product element represents a product node in a parent reaction element.

    Attributes
    ==========
    subelement_tags : Tuple[str]
        Possible subelement tags.

    id : str, None
        ID of product unique to map corresponding to a compound entry.

    name : str, None
        KEGG ID of the compound.

    children : Dict[str, List[str]]
        Keys are subelement tags, values are lists of subelement UUIDs.
    """
    tag = 'product'
    attribute_required = {
        'id': True,
        'name': True
    }
    subelement_tags: Tuple[str] = (
        'alt',
    )

    def __init__(self) -> None:
        self.id: str = None
        self.name: str = None

        self.children: Dict[str, List[str]] = {n: [] for n in self.subelement_tags}

        super().__init__()

class Alt(Element):
    """
    An alt element specifies an alternative name of a parent substrate or product element.

    Attributes
    ==========
    name : str, None
        Alternative KEGG ID of the compound.
    """
    tag = 'alt'
    attribute_required = {
        'name': True
    }

    def __init__(self) -> None:
        self.name: str = None

        super().__init__()

class XMLOps:
    """
    This class loads KGML (XML) files into memory in an object-oriented framework, and converts KGML
    objects into a properly formatted string that can be written to an XML file.

    Attributes
    ==========
    subelement_indentation_increment : int
        Class variable setting the indentation increment of subelements relative to parents in an
        output KGML XML file. The value of 4 spaces is that used in KGML reference files.

    attribute_indentations : Dict[Tuple[str, str, str], int]
        Class variable setting the absolute indentation of element attributes placed on new lines in
        an output KGML XML file. Keys are tuples of element tag, name of the attribute before the
        line break, and name of the attribute after the line break; values are the number of spaces.
        The attributes placed on new lines and numbers of spaces are those used in KGML reference
        files.
    """
    subelement_indentation_increment: int = 4
    attribute_indentations: Dict[Tuple[str, str, str], int] = {
        ('pathway', 'number', 'title'): 9,
        ('pathway', 'title', 'image'): 9,
        ('pathway', 'image', 'link'): 9,
        ('entry', 'reaction', 'link'): 8,
        ('entry', 'type', 'link'): 8,
        ('graphics', 'bgcolor', 'type'): 13
    }

    def __init__(self) -> None:
        pass

    def load(self, kgml_filepath: str) -> Pathway:
        """
        Load a KGML file as element objects.

        Parameters
        ==========
        kgml_filepath : str
            Path to a KGML file.

        Returns
        =======
        Pathway
            KGML pathway element object containing all data from the file via subelements.
        """
        assert os.path.exists(kgml_filepath)

        with open(kgml_filepath, 'rb') as file:
            kgml_bytes = file.read()
        root = ET.fromstring(kgml_bytes)
        assert root.tag == Pathway.tag

        pathway: Pathway = self.load_element(root)
        pathway.xml_declaration, pathway.xml_doctype, pathway.xml_comment = [
            line.decode('utf-8') for line in kgml_bytes.split(b'\n')[: 3]
        ]

        return pathway

    def load_element(self, xml_element: ET.Element, pathway: Pathway = None) -> Element:
        """
        Load a KGML element object representing an XML element from a KGML file.

        Parameters
        ==========
        xml_element : xml.etree.ElementTree.Element
            XML element loaded from KGML file.

        pathway : Pathway
            The pathway object containing the loaded KGML element, None if the element being loaded
            is the pathway itself.

        Returns
        =======
        Element
            Object representing KGML element.
        """
        kgml_element_class = globals()[xml_element.tag.capitalize()]
        kgml_element: Element = kgml_element_class()

        # Consider each possible attribute of the KGML element.
        for attribute, is_required in kgml_element.attribute_required.items():
            try:
                value = xml_element.attrib[attribute]
            except KeyError:
                if is_required:
                    # The required attribute was not present.
                    error_message = ""
                    for a, v in xml_element.attrib.items():
                        error_message += f" '{a}': '{v}'"
                    raise AssertionError(
                        "An XML element was encountered that should but does not contain an "
                        f"attribute, '{attribute}'. Here is a list of the element's attributes "
                        f"read from the KGML file:{error_message}")
                else:
                    # The optional attribute was not present.
                    continue

            # Convert certain non-string values stored in KGML element objects.
            if kgml_element.tag == 'graphics':
                if attribute in ('x', 'y', 'width', 'height'):
                    value = float(value)
                elif attribute == 'coords':
                    value = tuple([int(coord) for coord in value.split(',')])

            setattr(kgml_element, attribute, value)

        # Recursively load subelements.
        for xml_subelement in xml_element:
            assert xml_subelement.tag in kgml_element.subelement_tags
            if pathway is None:
                kgml_subelement = self.load_element(xml_subelement, pathway=kgml_element)
            else:
                kgml_subelement = self.load_element(xml_subelement, pathway=pathway)
            kgml_element.children[kgml_subelement.tag].append(kgml_subelement.uuid)

        if pathway is None:
            return kgml_element

        # Map unique IDs to element objects.
        pathway.uuid_element_lookup[kgml_element.uuid] = kgml_element

        # Map KEGG IDs to element objects; IDs are not necessarily unique to objects. Note that the
        # first part of the KEGG ID as it appears in the KGML file is stripped. For example,
        # 'ko:K01080', 'cpd:C12144', and 'path:map00604' become 'K01080', 'C12144', and 'map00604'
        # in the dictionary keys.
        if kgml_element.tag in (
            'pathway',
            'entry',
            'reaction',
            'substrate',
            'product',
            'alt'
        ):
            for kegg_name in kgml_element.name.split():
                # Some Entry names are "undefined".
                split_kegg_name = kegg_name.split(':')
                if len(split_kegg_name) != 2:
                    continue
                kegg_id = split_kegg_name[1]
                try:
                    pathway.kegg_id_element_lookup[kegg_id].append(kgml_element)
                except KeyError:
                    pathway.kegg_id_element_lookup[kegg_id] = [kgml_element]

        return kgml_element

    def write(self, pathway: Pathway, output_filepath: str) -> None:
        """
        Write a KGML object representation as a formatted KGML (XML) file.

        Parameters
        ==========
        pathway : Pathway
            KGML pathway element object.

        output_filepath : str
            Path to KGML (XML) output file.

        Returns
        =======
        None
        """
        assert is_output_file_writable(output_filepath)
        with open(output_filepath, 'w') as file:
            file.write(self.get_str(pathway))

    def get_str(self, pathway: Pathway) -> str:
        """
        Convert a KGML object representation to a formatted KGML XML string.

        Parameters
        ==========
        pathway : Pathway
            KGML pathway element object.

        Returns
        =======
        str
            Formatted XML string ready to write as a KGML file.
        """
        assert pathway.tag == Pathway.tag
        tree = self.get_tree(pathway, str_values=True)

        output_str = ''

        # Record KGML XML metadata, stored in the pathway object. Note that special HTML characters
        # (&, >, <) are not escaped, as these were not observed in the metadata of KGML reference
        # files.
        xml_declaration = pathway.xml_declaration
        if xml_declaration is not None:
            assert xml_declaration[: 6] == '<?xml ' and xml_declaration[-2: ] == '?>'
            output_str += f'<?xml {xml_declaration[6: -2]}?>\n'

        xml_doctype = pathway.xml_doctype
        if xml_doctype is not None:
            assert xml_doctype[: 10] == '<!DOCTYPE ' and xml_doctype[-1: ] == '>'
            output_str += f'<!DOCTYPE {xml_doctype[10: -1]}>\n'

        xml_comment = pathway.xml_comment
        if xml_comment is not None:
            assert xml_comment[: 5] == '<!-- ' and xml_comment[-4: ] == ' -->'
            output_str += f'<!-- {xml_comment[5: -4]} -->\n'

        output_str += self.get_indented_str(tree.getroot())

        return output_str

    def get_tree(self, pathway: Pathway, str_values: bool = False) -> ET.ElementTree:
        """
        Convert a KGML object representation to an XML tree.

        Parameters
        ==========
        pathway : Pathway
            KGML pathway element object.

        str_values : bool, False
            If True, convert non-string attribute values stored in KGML element objects to strings
            like those in KGML reference files.

        Returns
        =======
        xml.etree.ElementTree.ElementTree
            XML representation of KGML elements, with the root pathway element.
        """
        root = self.get_element(pathway, str_values=str_values)
        tree = ET.ElementTree(root)
        return tree

    def get_element(
        self,
        kgml_element: Element,
        pathway: Pathway = None,
        str_values: bool = False
    ) -> ET.Element:
        """
        Convert a KGML element object representation to an XML element.

        Parameters
        ==========
        kgml_element : Element
            KGML element object to convert into an XML element.

        pathway : Pathway
            Pathway object containing the KGML element, needed if the element is not a pathway
            element.

        str_values : bool, False
            If True, convert non-string attribute values stored in the KGML element object to
            strings like those in KGML reference files.

        Returns
        =======
        xml.etree.ElementTree.Element
            XML representation of the KGML element.
        """
        if pathway is None:
            if kgml_element.tag != 'pathway':
                raise ValueError(
                    "The 'pathway' element containing 'kgml_element' must be given as an argument."
                )
            pathway = kgml_element

        xml_element = ET.Element(kgml_element.tag)

        for attribute in kgml_element.attribute_required:
            # The KGML element object should have a value of each possible attribute, even those not
            # required, for which the default value is None.
            value = getattr(kgml_element, attribute)
            if value is None:
                continue
            if str_values:
                if kgml_element.tag == 'graphics':
                    if attribute in ('x', 'y', 'width', 'height'):
                        value = str(round(value))
                    elif attribute == 'coords':
                        value = ','.join([str(round(coord)) for coord in value])
            xml_element.attrib[attribute] = value

        if not hasattr(kgml_element, 'children'):
            return xml_element

        # Recursively add XML subelements.
        for subelement_uuids in kgml_element.children.values():
            for uuid in subelement_uuids:
                kgml_subelement = pathway.uuid_element_lookup[uuid]
                xml_subelement = self.get_element(kgml_subelement, pathway=pathway)
                xml_element.append(xml_subelement)

        return xml_element

    def get_indented_str(self, xml_element: ET.Element, tag_indentation: int = 0) -> str:
        """
        Convert a KGML XML element to a string with formatting, including indentation, that is
        consistent with reference KGML files.

        Parameters
        ==========
        xml_element : xml.etree.ElementTree
            XML element loaded from KGML file.

        tag_indentation : int, 0
            Indentation of the element tag, which, by default, is 0 for the root element, and
            increments by the class variable, `subelement_indentation_increment`, in the recursive
            calls to this class for subelements.

        Returns
        =======
        str
            The formatted string for the element, ready to be written to a KGML file.
        """
        tag = xml_element.tag
        attributes = xml_element.attrib
        children = list(xml_element)

        if attributes:
            indented_output = f'{" " * tag_indentation}<{tag}'
            for i, (attr, value) in enumerate(attributes.items()):
                if isinstance(value, str):
                    value: str
                    # ">" and "<" are encountered in attribute values and are represented by HTML
                    # escape characters in KGML files.
                    assert '&' not in value
                    v = value.replace(">", "&gt;").replace("<", "&lt;")
                elif isinstance(value, float):
                    value: float
                    v = str(round(value))
                elif isinstance(value, tuple):
                    value: Tuple[float]
                    v = ','.join([str(round(coord)) for coord in value])
                else:
                    raise AssertionError(
                        f"The attribute, '{attr}', had a value, '{value}', of unrecognized type, "
                        f"'{type(value)}'."
                    )
                if i == 0:
                    indented_output += f' {attr}="{v}"'
                else:
                    attr_indentation = self.attribute_indentations.get((tag, prev_attr, attr))
                    if attr_indentation is None:
                        indented_output += f' {attr}="{v}"'
                    else:
                        indented_output += f'\n{" " * attr_indentation}{attr}="{v}"'
                prev_attr = attr
        else:
            indented_output = f'{" " * tag_indentation}<{tag}'
            raise AssertionError(
                "It is assumed that all KGML XML elements have attributes, but an element with the "
                f"tag, '{tag}', did not have any."
            )

        if children:
            indented_output += '>\n'
            for child in children:
                indented_output += self.get_indented_str(
                    child, tag_indentation=tag_indentation + self.subelement_indentation_increment
                )
            # End tag
            indented_output += f'{" " * tag_indentation}</{tag}>\n'
        else:
            # Tags without possible subelements, plus substrate and product, which have possible alt
            # subelement, which is never present in reference files.
            if tag in ('graphics', 'component', 'subtype', 'substrate', 'product', 'alt'):
                # Self-closing tag
                indented_output += '/>\n'
            else:
                indented_output += '>\n'
                # End tag
                indented_output += f'{" " * tag_indentation}</{tag}>\n'

        if xml_element.text:
            raise AssertionError(
                "It is assumed that KGML XML elements do not contain 'text', but an element with "
                f"the tag, '{tag}', contained the text, '{xml_element.text}'."
            )

        return indented_output

class CircleGrid:
    """
    Circles sorted into the square cells of a grid, so that the circles near a point are found fast.

    A cell is wider than the largest circle and the reach of the grid. A circle that holds a point,
    or whose edge lies within reach of it, therefore has its center in the cell of the point or in
    one of the eight cells around it.

    Attributes
    ==========
    reach : float
        How far outside a circle a point can lie and still be near it, in points.

    cell_width : float
        The width of a cell, in points.

    cells : Dict[Tuple[int, int], List[Tuple[float, float, float]]]
        Keys are the column and row of a cell. Each value holds the circles centered in the cell,
        each as its x, y, and radius.
    """
    def __init__(self, circles: Iterable[Tuple[float, float, float]], reach: float = 0.0) -> None:
        """
        Parameters
        ==========
        circles : Iterable[Tuple[float, float, float]]
            The x, y, and radius of each circle, in points.

        reach : float, 0.0
            How far outside a circle a point can lie and still be near it, in points.
        """
        circles = list(circles)
        self.reach = reach
        # A cell is 1 point wider than the largest circle and the reach, so it has a width even with
        # no circles.
        self.cell_width = 2 * max([radius for _, _, radius in circles], default=0.0) + reach + 1
        self.cells: Dict[Tuple[int, int], List[Tuple[float, float, float]]] = {}
        for x, y, radius in circles:
            self.cells.setdefault(self._get_cell(x, y), []).append((x, y, radius))

    def _get_cell(self, x: float, y: float) -> Tuple[int, int]:
        """Get the column and row of the cell that holds a point."""
        return int(x // self.cell_width), int(y // self.cell_width)

    def get_nearby_circles(self, x: float, y: float) -> List[Tuple[float, float, float]]:
        """Get the circles centered in the cell of a point and in the eight cells around it."""
        column, row = self._get_cell(x, y)
        return [
            circle
            for nearby_column in (column - 1, column, column + 1)
            for nearby_row in (row - 1, row, row + 1)
            for circle in self.cells.get((nearby_column, nearby_row), ())
        ]

    def get_nearest_circle(self, x: float, y: float) -> Union[Tuple[float, float, float], None]:
        """
        Get the circle whose edge lies nearest a point, if the edge lies within reach of the point.

        A point inside a circle lies a negative distance outside its edge.

        Parameters
        ==========
        x : float
            The x coordinate of the point.

        y : float
            The y coordinate of the point.

        Returns
        =======
        Union[Tuple[float, float, float], None]
            The x, y, and radius of the circle, or None if no circle lies within reach.
        """
        nearest_circle = None
        nearest_distance = math.inf
        for circle in self.get_nearby_circles(x, y):
            circle_x, circle_y, radius = circle
            distance = math.hypot(circle_x - x, circle_y - y) - radius
            if distance < nearest_distance:
                nearest_circle = circle
                nearest_distance = distance
        if nearest_distance > self.reach:
            return None
        return nearest_circle

def get_drawn_circle_radius(width: float) -> float:
    """
    Get the radius of a circle as Biopython draws it, which includes half of its outline.

    Parameters
    ==========
    width : float
        The width of the circle in KGML.

    Returns
    =======
    float
        The radius of the circle with its outline, in points.
    """
    return width * 0.5 + CIRCLE_OUTLINE_HALF_WIDTH

def get_box_fill_width(width: float) -> float:
    """
    Get the width of the fill of a box inside a ring, where the columns of its slices are drawn.

    Parameters
    ==========
    width : float
        The width of the box in KGML.

    Returns
    =======
    float
        The width between the black lines inside the ring, in points ('BOX_RING_WIDTH').
    """
    return width - 2 * (BOX_OUTLINE_HALF_WIDTH + BOX_RING_WIDTH + RING_LINE_WIDTH)

def get_circle_ring_width(width: float) -> float:
    """
    Get the width of the ring of a compound circle.

    Parameters
    ==========
    width : float
        The width of the circle in KGML.

    Returns
    =======
    float
        The width of the ring in points, wider for a circle wider than 'SMALL_CIRCLE_MAX_WIDTH'.
    """
    return LARGE_CIRCLE_RING_WIDTH if width > SMALL_CIRCLE_MAX_WIDTH else SMALL_CIRCLE_RING_WIDTH

def get_compound_circle_grid(pathway: Pathway) -> CircleGrid:
    """
    Sort the circles of the compound entries of a map into a grid.

    Every circle with a width is put in the grid, whatever its color. 'SlicedKGMLCanvas' leaves out
    the circles that it draws see-through, so its own grid holds no circle that this grid lacks.

    Parameters
    ==========
    pathway : Pathway
        Object representation of a KGML file.

    Returns
    =======
    CircleGrid
        The compound circles, each with the radius that Biopython draws.
    """
    circles: List[Tuple[float, float, float]] = []
    for entry in pathway.get_entries(entry_type='compound'):
        for graphics_uuid in entry.children['graphics']:
            graphics: Graphics = pathway.uuid_element_lookup[graphics_uuid]
            if graphics.type == 'circle' and graphics.width is not None:
                circles.append((graphics.x, graphics.y, get_drawn_circle_radius(graphics.width)))
    return CircleGrid(circles)

def get_hidden_length(
    points: List[Tuple[float, float]],
    grid: CircleGrid,
    line_width: float
) -> float:
    """
    Measure how much of a line the circles of a grid hide, starting from its first point.

    A walk starts at the first point if a circle holds it. It follows the line through bends, and
    from one circle into another that overlaps it, until the line leaves the last circle. The edges
    of a wide line leave a circle a little before the middle of the line does, so the hidden length
    ends where the edges leave. This is worked out for a line that runs along a radius of the last
    circle, as a line from the center of a compound circle does. For a line that leaves a circle
    away from its center, the hidden length is a little too long, by up to about 3 points on global
    maps.

    Parameters
    ==========
    points : List[Tuple[float, float]]
        The x and y coordinates of the points of the line, in order.

    grid : CircleGrid
        The circles drawn over the line.

    line_width : float
        The width of the line, in points.

    Returns
    =======
    float
        The length along the line from its first point to where the line shows, in points. This is 0
        where no circle holds the first point, and the length of the whole line where circles hide
        all of it.
    """
    hidden = 0.0
    exit_radius = 0.0
    point = points[0]
    i = 1
    while i < len(points):
        next_point = points[i]
        step = math.dist(point, next_point)
        if step == 0:
            i += 1
            continue
        dx = (next_point[0] - point[0]) / step
        dy = (next_point[1] - point[1]) / step
        # Find how far the line runs before it leaves every circle that holds the point. The run in
        # one circle is the distance along the line to where it crosses the edge of the circle.
        run = 0.0
        run_radius = 0.0
        for x, y, radius in grid.get_nearby_circles(*point):
            offset_x = point[0] - x
            offset_y = point[1] - y
            squared_distance = offset_x * offset_x + offset_y * offset_y
            if squared_distance > (radius + CIRCLE_TOLERANCE) ** 2:
                continue
            projection = offset_x * dx + offset_y * dy
            circle_run = -projection + math.sqrt(
                max(projection * projection - squared_distance + radius * radius, 0.0)
            )
            if circle_run > run:
                run = circle_run
                run_radius = radius
        if run <= CIRCLE_TOLERANCE:
            break
        exit_radius = run_radius
        if run < step:
            hidden += run
            point = (point[0] + dx * run, point[1] + dy * run)
        else:
            hidden += step
            point = next_point
            i += 1
    if hidden == 0 or i == len(points):
        # No circle holds the first point, or circles hide the whole line.
        return hidden
    half_width = line_width * 0.5
    back_off = exit_radius - math.sqrt(max(exit_radius ** 2 - half_width ** 2, 0.0))
    return max(hidden - back_off, 0.0)

def get_visible_span(
    points: List[Tuple[float, float]],
    grid: CircleGrid,
    line_width: float
) -> Tuple[float, float]:
    """
    Find the part of a line that the circles of a grid leave visible.

    Parameters
    ==========
    points : List[Tuple[float, float]]
        The x and y coordinates of the points of the line, in order.

    grid : CircleGrid
        The circles drawn over the line.

    line_width : float
        The width of the line, in points.

    Returns
    =======
    Tuple[float, float]
        The length along the line from its first point to the start of the visible part, and the
        length of the visible part, in points ('get_hidden_length'). The second is 0 or less where
        circles hide the whole line.
    """
    length = sum([math.dist(points[i], points[i + 1]) for i in range(len(points) - 1)])
    start = get_hidden_length(points, grid, line_width)
    end = get_hidden_length(points[::-1], grid, line_width)
    return start, length - start - end

class SlicedKGMLCanvas(KGMLCanvas):
    """
    Biopython's KGML canvas, extended to draw ortholog entries in slices and to draw rings on boxes
    and compound circles.

    A rectangle with a glyph is drawn from that display list ('Graphics.glyph'). Its cells are
    columns, left to right, with a thin gray line between columns. It has one outline, and one label
    where ortholog entries are labeled ('label_orthologs'), unless its glyph drops the label. A
    compound circle with a glyph is drawn from its track, wedges, and rules, inside its ring and its
    outline. A sliced line is drawn as repeating dashes with no gaps, one dash per slice in each
    cycle. The dashes fill the part of the line that compound circles leave visible.

    A ring colors the border of a rectangle or a compound circle ('Graphics.border'). It lies just
    inside the black outline, with a thin black line between the ring and the fill
    ('BOX_RING_WIDTH'). A rectangle with a ring is drawn as one column or as the columns of its
    glyph, which then fit inside the ring. Its ring can be whole, or divided into bands that pair
    with the columns, or into segments of its own. A circle with a ring is drawn whole. Lines have
    no rings.

    Everything else is drawn by Biopython's own methods, as in KGMLCanvas. The canvas overrides
    three private methods of KGMLCanvas and calls two others. Biopython can rename private methods
    in any release, so the canvas checks that all five exist.

    Attributes
    ==========
    entry_glyphs : Dict[int, List[Union[Dict[str, Union[str, list]], None]]]
        Keys are the IDs of entries with glyphs. Each value holds the glyph of each Graphics element
        of the entry, in order, as in 'Graphics.glyph'. An element without a glyph has None.

    entry_slices : Dict[int, List[Union[List[Tuple[str, float]], None]]]
        Keys are the IDs of entries with sliced lines. Each value holds the slices of each Graphics
        element of the entry, in order, as in 'Graphics.slices'. An element drawn whole has None.

    entry_borders : Dict[int, List[Union[str, Tuple[str, list], None]]]
        Keys are the IDs of entries with rings. Each value holds the ring of each Graphics element
        of the entry, in order, as in 'Graphics.border'. An element without a ring has None.

    circle_grid : Union[CircleGrid, None]
        The compound circles drawn opaque over lines ('_get_circle_grid'). It is None until the
        first sliced line is drawn.
    """
    private_methods: Tuple[str] = (
        '_KGMLCanvas__add_orthologs',
        '_KGMLCanvas__add_reaction_entries',
        '_KGMLCanvas__add_compounds',
        '_KGMLCanvas__add_graphics',
        '_KGMLCanvas__add_labels'
    )

    def __init__(
        self,
        pathway: KGML_pathway.Pathway,
        entry_slices: Dict[int, List[Union[List[Tuple[str, float]], None]]],
        entry_borders: Union[Dict[int, List[Union[str, Tuple[str, list], None]]], None] = None,
        entry_glyphs: Union[Dict[int, List[Union[Dict[str, Union[str, list]], None]]], None] = None,
        **kwargs
    ) -> None:
        """
        Parameters
        ==========
        pathway : Bio.KEGG.KGML.KGML_pathway.Pathway
            Biopython representation of a KGML file.

        entry_slices : Dict[int, List[Union[List[Tuple[str, float]], None]]]
            The slices of each entry with sliced lines, as in the attribute.

        entry_borders : Union[Dict[int, List[Union[str, Tuple[str, list], None]]], None], None
            The rings of each entry with rings, as in the attribute, or None for no rings.

        entry_glyphs : Union[Dict[int, List[Union[Dict[str, Union[str, list]], None]]], None], None
            The glyphs of each entry with glyphs, as in the attribute, or None for no glyphs.

        **kwargs
            Arguments to KGMLCanvas.
        """
        self.check_biopython()
        super().__init__(pathway, **kwargs)
        self.entry_slices = entry_slices
        self.entry_borders = {} if entry_borders is None else entry_borders
        self.entry_glyphs = {} if entry_glyphs is None else entry_glyphs
        # The grid is built when the first sliced line is drawn, since Drawer sets
        # 'non_reactant_transparency' after it makes the canvas.
        self.circle_grid: Union[CircleGrid, None] = None

    @classmethod
    def check_biopython(cls) -> None:
        """
        Raise a ConfigError if KGMLCanvas lacks a private method that this canvas uses.

        The canvas runs the check when it is made. Callers can also run it before any map is drawn,
        so that no output is written before the error.
        """
        missing = [name for name in cls.private_methods if not hasattr(KGMLCanvas, name)]
        if missing:
            missing_message = ', '.join([f"'{name}'" for name in missing])
            raise ConfigError(
                f"Biopython {Bio.__version__} cannot draw map elements in slices, with borders, or "
                f"with bars. Slices, borders, and bars are drawn with private methods of "
                f"Biopython's KGMLCanvas, and this version lacks these: {missing_message}. Drop "
                f"'--reaction-slice-by', '--reaction-border-column', '--compound-border-column', "
                f"'--reaction-bar-column', and '--compound-bar-column' to draw the maps without "
                f"slices, borders, or bars, or install a version of Biopython that has these "
                f"methods, such as 1.87."
            )

    def _KGMLCanvas__add_reaction_entries(self) -> None:
        """
        Draw reaction entries as Biopython does, except for entries with glyphs, slices, or rings.

        Biopython draws an ortholog entry with a reaction twice, here and in the ortholog pass. An
        entry with glyphs, slices, or rings is drawn once, in the ortholog pass.
        """
        for reaction in self.pathway.reaction_entries:
            if (
                reaction.id in self.entry_glyphs or reaction.id in self.entry_slices
                or reaction.id in self.entry_borders
            ):
                continue
            for g in reaction.graphics:
                self.drawing.setStrokeColor(color_to_reportlab(g.fgcolor))
                self.drawing.setFillColor(color_to_reportlab(g.bgcolor))
                self._KGMLCanvas__add_graphics(g)
                if self.label_reaction_entries:
                    self.drawing.setFillColor(darken(g.fgcolor))
                    self._KGMLCanvas__add_labels(g)

    def _KGMLCanvas__add_orthologs(self) -> None:
        """
        Draw ortholog entries as Biopython does, except for Graphics with glyphs, slices, or rings.

        Only rectangles take glyphs and rings, and only lines take slices. Anything else is drawn as
        Biopython draws it. A rectangle whose glyph drops its label is drawn without one.
        """
        for ortholog in self.pathway.orthologs:
            graphics_glyphs = self.entry_glyphs.get(ortholog.id)
            graphics_slices = self.entry_slices.get(ortholog.id)
            graphics_borders = self.entry_borders.get(ortholog.id)
            for i, g in enumerate(ortholog.graphics):
                glyph = None if graphics_glyphs is None else graphics_glyphs[i]
                slices = None if graphics_slices is None else graphics_slices[i]
                border = None if graphics_borders is None else graphics_borders[i]
                self.drawing.setStrokeColor(color_to_reportlab(g.fgcolor))
                self.drawing.setFillColor(color_to_reportlab(g.bgcolor))
                if (glyph or border) and g.type == 'rectangle':
                    self._draw_box_glyph(g, glyph, border)
                elif slices and g.type == 'line':
                    self._draw_sliced_line(g, slices)
                else:
                    self._KGMLCanvas__add_graphics(g)
                # A glyph can drop the label of its rectangle.
                if self.label_orthologs and (glyph is None or glyph.get('label', True)):
                    # The label is drawn in a darker shade of the fg color, as Biopython does.
                    self.drawing.setFillColor(darken(g.fgcolor))
                    self._KGMLCanvas__add_labels(g)

    def _KGMLCanvas__add_compounds(self) -> None:
        """
        Draw compound entries as Biopython does, and the glyphs and rings of circles that have them.

        This follows Biopython's own method step by step, so a map without glyphs or rings is drawn
        as Biopython draws it. A ring is drawn after its circle and before the label. A circle with
        a glyph is drawn from it ('_draw_circle_glyph').
        """
        for compound in self.pathway.compounds:
            graphics_glyphs = self.entry_glyphs.get(compound.id)
            graphics_borders = self.entry_borders.get(compound.id)
            # Biopython finds whether a compound takes part in a reaction by searching every
            # reaction of the map, so it is asked once.
            is_reactant = compound.is_reactant
            for i, g in enumerate(compound.graphics):
                glyph = None if graphics_glyphs is None else graphics_glyphs[i]
                border = None if graphics_borders is None else graphics_borders[i]
                self.drawing.setStrokeColor(color_to_reportlab(g.fgcolor))
                # Biopython draws a compound that takes part in no reaction see-through, by
                # 'non_reactant_transparency'.
                if glyph and g.type == 'circle':
                    alpha = 1.0 if is_reactant else self.non_reactant_transparency
                    self._draw_circle_glyph(g, glyph, border, alpha)
                else:
                    fillcolor = color_to_reportlab(g.bgcolor)
                    if not is_reactant:
                        fillcolor.alpha *= self.non_reactant_transparency
                    self.drawing.setFillColor(fillcolor)
                    self._KGMLCanvas__add_graphics(g)
                    if border and g.type == 'circle':
                        self._draw_circle_ring(g, border)
                if self.label_compounds:
                    if not is_reactant:
                        t = 0.3
                    else:
                        t = 1
                    self.drawing.setFillColor(reportlab_colors.Color(0.2, 0.2, 0.2, t))
                    self._KGMLCanvas__add_labels(g)

    def _draw_circle_ring(self, graphics: KGML_pathway.Graphics, color: str) -> None:
        """
        Draw a whole ring just inside the outline of a circle, and the black line inside the ring.

        A circle wider than 'SMALL_CIRCLE_MAX_WIDTH' takes the wider ring. The outline of the circle
        stays where it is ('get_drawn_circle_radius').

        Parameters
        ==========
        graphics : Bio.KEGG.KGML.KGML_pathway.Graphics
            A circle.

        color : str
            The color of the ring.
        """
        drawing = self.drawing
        ring = get_circle_ring_width(graphics.width)
        radius = graphics.width * 0.5 - CIRCLE_OUTLINE_HALF_WIDTH
        drawing.saveState()
        drawing.setLineWidth(ring)
        drawing.setStrokeColor(color_to_reportlab(color))
        drawing.circle(graphics.x, graphics.y, radius - ring * 0.5, stroke=1, fill=0)
        drawing.setLineWidth(RING_LINE_WIDTH)
        drawing.setStrokeColor(color_to_reportlab(RING_LINE_COLOR))
        drawing.circle(
            graphics.x, graphics.y, radius - ring - RING_LINE_WIDTH * 0.5, stroke=1, fill=0
        )
        drawing.restoreState()

    def _draw_circle_glyph(
        self,
        graphics: KGML_pathway.Graphics,
        glyph: Dict[str, Union[str, list]],
        border: Union[str, None],
        alpha: float
    ) -> None:
        """
        Draw a circle from its glyph ('Graphics.glyph'), and its ring.

        The track fills the circle. Each wedge runs clockwise from 12 o'clock over its part of a
        turn, from the center to the middle of the outline, so that no seam shows between the wedge
        and the outline drawn over it. A wedge of a whole turn is a disc. The track and the wedges
        reach under the ring, which is drawn over them. The rules are drawn over the wedges with
        square ends that stop at their own ends. Each runs along a radius from the center to the
        edge of the part of the circle that it shows: the inside of the outline, or the inside of
        the black line within the ring. The ring is drawn next ('_draw_circle_ring'), and the
        outline last, in the stroke color of the canvas, as Biopython draws it. The track and the
        wedges take the alpha of the circle, as Biopython's fill does, and the rules, the ring, and
        the outline are drawn opaque. A circle too small to show anything inside its outline, such
        as one of width 0, draws its rules with no length.

        Parameters
        ==========
        graphics : Bio.KEGG.KGML.KGML_pathway.Graphics
            A circle. The stroke color of the canvas is already set to its fg color.

        glyph : Dict[str, Union[str, list]]
            The display list of the circle, as in 'Graphics.glyph'.

        border : Union[str, None]
            The color of the ring, or None for no ring.

        alpha : float
            The opacity of the fill of the circle, which is below 1 for a compound that takes part
            in no reaction ('non_reactant_transparency').
        """
        drawing = self.drawing
        x, y = graphics.x, graphics.y
        radius = graphics.width * 0.5
        # The canvas draws downward, so 12 o'clock is at -90 degrees, and angles grow clockwise.
        drawing.saveState()
        track = color_to_reportlab(glyph['track'])
        track.alpha *= alpha
        drawing.setFillColor(track)
        drawing.circle(x, y, radius, stroke=0, fill=1)
        for start, end, color in glyph.get('wedges', []):
            if end <= start:
                continue
            fill = color_to_reportlab(color)
            fill.alpha *= alpha
            drawing.setFillColor(fill)
            if end - start >= 1:
                drawing.circle(x, y, radius, stroke=0, fill=1)
                continue
            path = drawing.beginPath()
            path.moveTo(x, y)
            path.arcTo(
                x - radius, y - radius, x + radius, y + radius, startAng=-90 + 360 * start,
                extent=360 * (end - start)
            )
            path.close()
            drawing.drawPath(path, stroke=0, fill=1)
        shown_radius = radius - CIRCLE_OUTLINE_HALF_WIDTH
        if border:
            shown_radius -= get_circle_ring_width(graphics.width) + RING_LINE_WIDTH
        shown_radius = max(shown_radius, 0.0)
        drawing.setLineCap(0)
        for position, color, line_width in glyph.get('rules', []):
            angle = 2 * math.pi * position
            drawing.setLineWidth(line_width)
            drawing.setStrokeColor(color_to_reportlab(color))
            drawing.line(
                x, y, x + shown_radius * math.sin(angle), y - shown_radius * math.cos(angle)
            )
        drawing.restoreState()
        if border:
            self._draw_circle_ring(graphics, border)
        drawing.circle(x, y, radius, stroke=1, fill=0)

    def _draw_box_glyph(
        self,
        graphics: KGML_pathway.Graphics,
        glyph: Union[Dict[str, list], None],
        border: Union[str, Tuple[str, list], None] = None
    ) -> None:
        """
        Draw a rectangle from its glyph ('Graphics.glyph'), and its ring.

        The cells of the glyph are columns, left to right. Without a ring, the columns span the
        whole rectangle. With a ring, each column takes its share of the width of the fill inside
        the ring ('get_box_fill_width'), so each column shows its exact share. Every column is still
        drawn the full height of the rectangle, and the first and last columns reach its left and
        right sides. Where a ring has no band, it shows the color of the column under it, not the
        base map.

        The marks and rules of the glyph are drawn over the columns, in the part of each column that
        the box shows: inside the outline, or inside the black line within the ring. An edge of a
        mark that lies on the edge of that part reaches 0.5 pt past it, under the outline or the
        black line, so that no seam shows there. Its other edges are exact, so a bar keeps its exact
        length. Rules are drawn with square ends that stop at their own ends.

        The ring is drawn over the columns, marks, and rules. Then come the gray lines between
        columns, the black line inside the ring, and the outline.

        Parameters
        ==========
        graphics : Bio.KEGG.KGML.KGML_pathway.Graphics
            A rectangle. The stroke color of the canvas is already set to its fg color.

        glyph : Union[Dict[str, list], None]
            The display list of the rectangle, as in 'Graphics.glyph', or None for one column in the
            bg color of the rectangle.

        border : Union[str, Tuple[str, list], None], None
            The ring, as in 'Graphics.border', or None for no ring.
        """
        drawing = self.drawing
        # KGML gives the center of a rectangle, and reportlab takes a corner.
        x0 = graphics.x - graphics.width * 0.5
        y0 = graphics.y - graphics.height * 0.5
        x1 = x0 + graphics.width
        cells = [(graphics.bgcolor, 1.0)] if glyph is None else glyph['cells']
        if border is None:
            fill_x0, fill_width, fill_x1 = x0, graphics.width, x1
        else:
            fill_width = get_box_fill_width(graphics.width)
            fill_x0 = x0 + (graphics.width - fill_width) * 0.5
            fill_x1 = fill_x0 + fill_width

        drawing.saveState()
        # Each column reaches 0.5 pt into the next. Without the overlap, anti-aliasing leaves thin
        # seams that show the base map. No column reaches past the right side of the box.
        start = 0.0
        boundaries: List[float] = []
        columns: List[Tuple[float, float]] = []
        for index, (color, share) in enumerate(cells):
            left = fill_x0 + start * fill_width
            right = min(left + share * fill_width + 0.5, fill_x1)
            drawing.setFillColor(color_to_reportlab(color))
            if border is None:
                drawing.rect(left, y0, right - left, graphics.height, stroke=0, fill=1)
            else:
                draw_left = x0 if index == 0 else left
                draw_right = x1 if index == len(cells) - 1 else right
                drawing.rect(
                    draw_left, y0, draw_right - draw_left, graphics.height, stroke=0, fill=1
                )
            if start > 0:
                boundaries.append(left)
            columns.append((left, fill_x0 + (start + share) * fill_width))
            start += share

        if glyph is not None and (glyph.get('marks') or glyph.get('rules')):
            self._draw_box_marks(
                x0, y0, graphics.width, graphics.height, columns, glyph, border is not None
            )

        if isinstance(border, str):
            self._draw_box_ring(x0, y0, graphics.width, graphics.height, border)
        elif border is not None and border[0] == 'bands':
            self._draw_box_ring_bands(x0, y0, graphics.width, graphics.height, columns, border[1])
        elif border is not None:
            self._draw_box_ring_segments(x0, y0, graphics.width, graphics.height, border[1])

        # A thin line between columns shows each column when neighbors have the same color. It
        # crosses bands, which share the boundaries of the columns, and stops at any other ring.
        line_y0, line_y1 = y0, y0 + graphics.height
        if border is not None and (isinstance(border, str) or border[0] != 'bands'):
            inset = (graphics.width - fill_width) * 0.5
            line_y0, line_y1 = y0 + inset, y0 + graphics.height - inset
        drawing.setLineWidth(0.5)
        drawing.setStrokeColor(color_to_reportlab(SLICE_DIVIDER_COLOR))
        for x in boundaries:
            drawing.line(x, line_y0, x, line_y1)

        # A thin black line encloses the fill inside the ring.
        if border is not None:
            inset = BOX_OUTLINE_HALF_WIDTH + BOX_RING_WIDTH + RING_LINE_WIDTH * 0.5
            drawing.setLineWidth(RING_LINE_WIDTH)
            drawing.setStrokeColor(color_to_reportlab(RING_LINE_COLOR))
            drawing.rect(
                x0 + inset, y0 + inset, graphics.width - 2 * inset, graphics.height - 2 * inset,
                stroke=1, fill=0
            )
        drawing.restoreState()

        drawing.rect(x0, y0, graphics.width, graphics.height, stroke=1, fill=0)

    def _draw_box_marks(
        self,
        x0: float,
        y0: float,
        width: float,
        height: float,
        columns: List[Tuple[float, float]],
        glyph: Dict[str, list],
        ringed: bool
    ) -> None:
        """
        Draw the marks and rules of the glyph of a rectangle, over its columns.

        Each mark and rule lies in the area of one column that the box shows. Across the column, the
        area runs inside the outline, or inside the black line within the ring. Along the column, it
        runs between the sides of the column, which are cut back to the same lines at the sides of
        the box. A band of a mark or rule is the bottom of the area, so many points tall.

        Parameters
        ==========
        x0, y0 : float
            The corner of the rectangle with the smallest coordinates.

        width, height : float
            The size of the rectangle.

        columns : List[Tuple[float, float]]
            The left and right sides of each column of the cells of the glyph.

        glyph : Dict[str, list]
            The display list of the rectangle ('Graphics.glyph'). Its marks and rules are drawn:
            - 'marks': [(cell, axis, start, end, color, band), ...], rectangles in the color. 'cell'
              is the position of the column. 'axis' is 'x', along the width from the left, or 'y',
              along the height from the bottom. 'start' and 'end' are fractions of the area along
              the axis. The mark spans the area across the axis. 'band' is None for the whole area,
              or the height in points of a band along its bottom.
            - 'rules': [(cell, axis, position, color, line_width, band), ...], lines across the area
              at a fraction along the axis, such as the end line of a bar.

        ringed : bool
            True if the rectangle has a ring.
        """
        drawing = self.drawing
        # The box shows its fill inside the outline, which Biopython draws 1 point wide, centered on
        # the edge. With a ring, it shows it inside the black line within the ring. The width of
        # either line inside that edge is how far an edge of a mark reaches under it.
        if ringed:
            inset = BOX_OUTLINE_HALF_WIDTH + BOX_RING_WIDTH + RING_LINE_WIDTH
            reach = RING_LINE_WIDTH
        else:
            inset = BOX_OUTLINE_HALF_WIDTH
            reach = BOX_OUTLINE_HALF_WIDTH
        shown_x0, shown_x1 = x0 + inset, x0 + width - inset
        shown_y0, shown_y1 = y0 + inset, y0 + height - inset

        def area(cell: int, band: Union[float, None]) -> Tuple[float, float, float, float]:
            left, right = columns[cell]
            top = shown_y0 if band is None else shown_y1 - band
            return max(left, shown_x0), top, min(right, shown_x1), shown_y1

        for cell, axis, start, end, color, band in glyph.get('marks', []):
            if end <= start:
                continue
            area_x0, area_y0, area_x1, area_y1 = area(cell, band)
            if axis == 'x':
                mark_x0 = area_x0 + start * (area_x1 - area_x0)
                mark_x1 = area_x0 + end * (area_x1 - area_x0)
                mark_y0, mark_y1 = area_y0, area_y1
            else:
                mark_x0, mark_x1 = area_x0, area_x1
                mark_y0 = area_y1 - end * (area_y1 - area_y0)
                mark_y1 = area_y1 - start * (area_y1 - area_y0)
            # An edge on the edge of the shown fill reaches under the line that bounds it.
            if mark_x0 <= shown_x0:
                mark_x0 -= reach
            if mark_x1 >= shown_x1:
                mark_x1 += reach
            if mark_y0 <= shown_y0:
                mark_y0 -= reach
            if mark_y1 >= shown_y1:
                mark_y1 += reach
            drawing.setFillColor(color_to_reportlab(color))
            drawing.rect(mark_x0, mark_y0, mark_x1 - mark_x0, mark_y1 - mark_y0, stroke=0, fill=1)

        drawing.setLineCap(0)
        for cell, axis, position, color, line_width, band in glyph.get('rules', []):
            area_x0, area_y0, area_x1, area_y1 = area(cell, band)
            drawing.setLineWidth(line_width)
            drawing.setStrokeColor(color_to_reportlab(color))
            if axis == 'x':
                x = area_x0 + position * (area_x1 - area_x0)
                drawing.line(x, area_y0, x, area_y1)
            else:
                y = area_y1 - position * (area_y1 - area_y0)
                drawing.line(area_x0, y, area_x1, y)

    def _draw_box_ring(self, x0: float, y0: float, width: float, height: float, color: str) -> None:
        """
        Draw a whole ring just inside the outline of a rectangle.

        Parameters
        ==========
        x0, y0 : float
            The corner of the rectangle with the smallest coordinates.

        width, height : float
            The size of the rectangle.

        color : str
            The color of the ring.
        """
        inset = BOX_OUTLINE_HALF_WIDTH + BOX_RING_WIDTH * 0.5
        self.drawing.setLineWidth(BOX_RING_WIDTH)
        self.drawing.setStrokeColor(color_to_reportlab(color))
        self.drawing.rect(
            x0 + inset, y0 + inset, width - 2 * inset, height - 2 * inset, stroke=1, fill=0
        )

    def _draw_box_ring_bands(
        self,
        x0: float,
        y0: float,
        width: float,
        height: float,
        columns: List[Tuple[float, float]],
        band_colors: List[Union[str, None]]
    ) -> None:
        """
        Draw a ring of a rectangle as bands that pair with the columns of its slices.

        Each column takes a band above it and a band below it, in the color that 'band_colors' gives
        the column. The first column also takes the left end of the ring, and the last column the
        right end. A column with None in 'band_colors' takes no band.

        Parameters
        ==========
        x0, y0 : float
            The corner of the rectangle with the smallest coordinates.

        width, height : float
            The size of the rectangle.

        columns : List[Tuple[float, float]]
            The left and right sides of each column in the fill.

        band_colors : List[Union[str, None]]
            The color of the bands of each column, or None for a column without bands.
        """
        drawing = self.drawing
        ring = BOX_RING_WIDTH
        outer = BOX_OUTLINE_HALF_WIDTH
        x1 = x0 + width
        top = y0 + outer
        bottom = y0 + height - outer - ring
        last = len(columns) - 1
        for index, ((left, right), color) in enumerate(zip(columns, band_colors)):
            if color is None:
                continue
            band_left = x0 + outer if index == 0 else left
            # A band reaches 0.5 pt under the next band, so that no seam shows between them.
            if index == last:
                band_right = x1 - outer
            elif band_colors[index + 1] is not None:
                band_right = right + 0.5
            else:
                band_right = right
            drawing.setFillColor(color_to_reportlab(color))
            drawing.rect(band_left, top, band_right - band_left, ring, stroke=0, fill=1)
            drawing.rect(band_left, bottom, band_right - band_left, ring, stroke=0, fill=1)
            if index == 0:
                drawing.rect(x0 + outer, top, ring, height - 2 * outer, stroke=0, fill=1)
            if index == last:
                drawing.rect(x1 - outer - ring, top, ring, height - 2 * outer, stroke=0, fill=1)

    def _draw_box_ring_segments(
        self,
        x0: float,
        y0: float,
        width: float,
        height: float,
        segments: List[Tuple[str, float]]
    ) -> None:
        """
        Draw a ring of a rectangle as segments, clockwise from the top left corner.

        Each segment takes its share of the length of the middle line of the ring. The ring is laid
        out as four strips, one for each side, clockwise from the top. Each strip starts with the
        corner where its side starts, and it is as long as its side of the middle line. A length
        along the middle line is then the same length along the strips. Each segment is drawn as one
        rectangle on each strip that it covers. A thin gray line crosses the ring between
        neighboring segments, as between columns, so that neighbors of one color can be told apart.

        Parameters
        ==========
        x0, y0 : float
            The corner of the rectangle with the smallest coordinates.

        width, height : float
            The size of the rectangle.

        segments : List[Tuple[str, float]]
            The color and share of each segment.
        """
        drawing = self.drawing
        ring = BOX_RING_WIDTH
        half = ring * 0.5
        # The corners of the middle line of the ring.
        a = x0 + BOX_OUTLINE_HALF_WIDTH + half
        b = y0 + BOX_OUTLINE_HALF_WIDTH + half
        c = x0 + width - BOX_OUTLINE_HALF_WIDTH - half
        d = y0 + height - BOX_OUTLINE_HALF_WIDTH - half

        # Each strip gives the rectangle that covers the lengths from 't0' to 't1' along it, and the
        # line across it at length 't'. The top strip runs right, the right strip down, the bottom
        # strip left, and the left strip up.
        def top(t0: float, t1: float) -> Tuple[float, float, float, float]:
            return a - half + t0, b - half, t1 - t0, ring

        def right(t0: float, t1: float) -> Tuple[float, float, float, float]:
            return c - half, b - half + t0, ring, t1 - t0

        def bottom(t0: float, t1: float) -> Tuple[float, float, float, float]:
            return c + half - t1, d - half, t1 - t0, ring

        def left(t0: float, t1: float) -> Tuple[float, float, float, float]:
            return a - half, d + half - t1, ring, t1 - t0

        def top_line(t: float) -> Tuple[float, float, float, float]:
            return a - half + t, b - half, a - half + t, b + half

        def right_line(t: float) -> Tuple[float, float, float, float]:
            return c - half, b - half + t, c + half, b - half + t

        def bottom_line(t: float) -> Tuple[float, float, float, float]:
            return c + half - t, d - half, c + half - t, d + half

        def left_line(t: float) -> Tuple[float, float, float, float]:
            return a - half, d + half - t, a + half, d + half - t

        strips = []
        start = 0.0
        for length, rect, line in (
            (c - a, top, top_line), (d - b, right, right_line), (c - a, bottom, bottom_line),
            (d - b, left, left_line)
        ):
            strips.append((start, start + length, rect, line))
            start += length
        perimeter = start

        # A segment with a share of 0 is not drawn. The last segment ends at the top left corner,
        # where the first segment starts.
        segments = [(color, share) for color, share in segments if share > 0]
        boundaries: List[float] = []
        position = 0.0
        for index, (color, share) in enumerate(segments):
            segment_start = position
            segment_end = perimeter if index == len(segments) - 1 else position + share * perimeter
            position = segment_end
            if index > 0:
                boundaries.append(segment_start)
            drawing.setFillColor(color_to_reportlab(color))
            for strip_start, strip_end, rect, _ in strips:
                piece_start = max(segment_start, strip_start)
                piece_end = min(segment_end, strip_end)
                if piece_end <= piece_start:
                    continue
                # Each rectangle reaches 0.5 pt into the next one, which is drawn later, so that
                # no seam shows. The first segment is drawn first, so no rectangle reaches into it.
                reach = min(0.5, perimeter - piece_end)
                x, y, w, h = rect(piece_start - strip_start, piece_end - strip_start + reach)
                drawing.rect(x, y, w, h, stroke=0, fill=1)
        # The last segment also meets the first, at the top left corner.
        if boundaries:
            boundaries.append(perimeter)

        drawing.setLineWidth(0.5)
        drawing.setStrokeColor(color_to_reportlab(SLICE_DIVIDER_COLOR))
        for boundary in boundaries:
            # A boundary at a corner lies at the end of the strip before the corner. Rounding can
            # put a boundary meant for a corner a tiny distance past it, so a small tolerance keeps
            # it on that strip. The strips are searched in order.
            for strip_start, strip_end, _, line in strips:
                if boundary <= strip_end + 1e-9:
                    drawing.line(*line(min(boundary, strip_end) - strip_start))
                    break

    def _draw_sliced_line(
        self,
        graphics: KGML_pathway.Graphics,
        slices: List[Tuple[str, float]]
    ) -> None:
        """
        Draw a line as repeating dashes. Each cycle holds one dash per slice, in order.

        On global maps, compound circles cover the ends of almost every line. A whole number of
        cycles fills the part of the line that the circles leave visible ('get_visible_span'), so
        each slice shows its exact share of that part. The number of cycles is the largest that
        keeps the shortest dash at least twice as long as the line is wide, and there is always at
        least one cycle. Where circles hide the whole line, the cycles fill it from end to end.

        Parameters
        ==========
        graphics : Bio.KEGG.KGML.KGML_pathway.Graphics
            A line. The stroke color of the canvas is already set to its fg color.

        slices : List[Tuple[str, float]]
            The color and share of the length of each dash in a cycle.
        """
        drawing = self.drawing
        line_width = 1 if graphics.width is None else graphics.width
        coords = graphics.coords
        # The path is built as Biopython builds it, starting with a segment of zero length.
        path = drawing.beginPath()
        path.moveTo(*coords[0])
        for x, y in coords:
            path.lineTo(x, y)
        length = sum([math.dist(coords[i], coords[i + 1]) for i in range(len(coords) - 1)])
        if self.circle_grid is None:
            self.circle_grid = self._get_circle_grid()
        start, visible = get_visible_span(coords, self.circle_grid, line_width)
        if visible <= 0:
            start, visible = 0.0, length

        drawing.saveState()
        drawing.setLineWidth(line_width)
        # The line is first drawn solid in the fg color. This fills notches at corners where a dash
        # ends near a vertex.
        drawing.drawPath(path)
        if visible > 0:
            cycle_count = max(
                1, math.floor(visible * min([share for _, share in slices]) / (2 * line_width))
            )
            cycle = visible / cycle_count
            drawing.setLineCap(0)
            # Dashes overlap by 0.5 pt, so anti-aliasing leaves no seam. An overlap always lies
            # under the dash drawn later, so each dash shows its exact length. Each dash but the
            # last reaches into the next dash. The first dash also reaches back under the last dash
            # of the cycle before. The dashes repeat along the whole path, under the circles too,
            # and the first cycle starts where the line shows.
            for i, (color, share) in enumerate(slices):
                head = 0.5 if i == 0 else 0.0
                tail = 0.5 if i < len(slices) - 1 else 0.0
                dash = head + share * cycle + tail
                # reportlab writes the phase with Python's str(), which puts a number below 0.0001
                # in exponent form, such as '2e-05'. PDF readers reject that form, so the phase is
                # rounded.
                phase = round((cycle - start + head) % cycle, 4)
                drawing.setStrokeColor(color_to_reportlab(color))
                drawing.setDash([dash, max(cycle - dash, 0.001)], phase)
                drawing.drawPath(path)
                start += share * cycle
        # Restoring the state ends the dashes. Otherwise later lines and circles would be dashed.
        drawing.restoreState()

    def _get_circle_grid(self) -> CircleGrid:
        """
        Sort the compound circles that the canvas draws opaque into a grid.

        Biopython draws compounds after ortholog entries, so their circles cover the ends of lines.
        A circle with a see-through fill is left out. Where 'non_reactant_transparency' is below 1,
        Biopython draws a compound that takes part in no reaction see-through. No circle is drawn
        where 'show_compounds' is False.

        Returns
        =======
        CircleGrid
            The circles, each with the radius that Biopython draws ('get_drawn_circle_radius').
        """
        circles: List[Tuple[float, float, float]] = []
        if not self.show_compounds:
            return CircleGrid(circles)
        for compound in self.pathway.compounds:
            for graphics in compound.graphics:
                if graphics.type != 'circle' or color_to_reportlab(graphics.bgcolor).alpha < 1:
                    continue
                # Biopython searches every reaction of the map to find whether a compound takes part
                # in one. This is slow, so it is done only where the answer matters.
                if self.non_reactant_transparency < 1 and not compound.is_reactant:
                    continue
                circles.append((graphics.x, graphics.y, get_drawn_circle_radius(graphics.width)))
        return CircleGrid(circles)

class Drawer:
    """
    Write pathway map image files incorporating KGML data.

    A map with Graphics elements that have glyphs, slices, or rings is drawn with
    'SlicedKGMLCanvas'. Any other map is drawn with Biopython's 'KGMLCanvas'. Neither canvas labels
    a box or circle too small to see ('TINY_GRAPHICS_MAX_SIZE').

    Attributes
    ==========
    kegg_context : anvio.metabolism.context.KeggContext
        This contains anvi'o KEGG database attributes, such as filepaths.

    xml_ops : XMLOps
        Loads KGML files.

    overwrite_output : bool
        If True, methods in this class overwrite existing output files.

    run : anvio.terminal.Run
        This object prints run information to the terminal.

    progress : anvio.terminal.Progress
        This object prints transient progress information to the terminal.

    non_reactant_transparency : float, 1.0
        This controls the transparency, or alpha, of the background color of compound graphics
        rendered from KGML for non-reactants, or compounds that don't participate in reactions. This
        value is used to set the attribute of 'Bio.Graphics.KGML_vis.KGMLCanvas'; the default value
        of 0.3 allows the color of compound circles in underlying global base maps to bleed through,
        which is not desirable.
    """
    def __init__(
        self,
        kegg_dir: str = None,
        overwrite_output: bool = FORCE_OVERWRITE,
        run: terminal.Run = terminal.Run(),
        progress: terminal.Progress = terminal.Progress()
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

        run : anvio.terminal.Run, anvio.terminal.Run()
            This object prints run information to the terminal.

        progress : anvio.terminal.Progress, anvio.terminal.Progress()
            This object prints transient progress information to the terminal.
        """
        args = Namespace()
        args.kegg_data_dir = kegg_dir
        self.kegg_context = KeggContext(args)

        self.xml_ops = XMLOps()

        self.overwrite_output = overwrite_output
        self.run = run
        self.progress = progress

        self.non_reactant_transparency = 1.0

    def draw_map(
        self,
        pathway: Pathway,
        output_filepath: str,
        map_filepath: str = None,
        use_org_map: bool = False,
        **kwargs
    ) -> None:
        """
        Draw a pathway map with KGML data as a PDF file.

        Parameters
        ==========
        pathway : Pathway
            Object representation of a KGML file.

        output_filepath : str
            Path to PDF output file containing the pathway map.

        map_filepath : str, None
            Path to pathway map image file to use as the base image of the output file. If None,
            then a PNG image file is automatically sought in the KEGG data directory, and it is
            assumed that the KGML data is scaled to fit the image size; here is more information on
            the files that are sought.

            For a standard or overview (not global) map, a 2x resolution 'map' file is sought. If
            the org attribute of the pathway object is an organism code and the argument,
            use_org_map, is True, then a 1x resolution organism-specific file is sought.

            For a global map, a 1x resolution file is sought; it is assumed the KGML data is scaled
            to fit this image size. The org attribute of the pathway object is used to seek the
            corresponding file, i.e., a 'ko' pathway containing reactions with KO IDs results in the
            reference 'ko' file being sought, whereas an 'ec' pathway containing reactions with EC
            number IDs corresponds to the reference 'ec' file. If the org attribute is an organism
            code and the argument, use_org_map, is True, then an organism-specific image file is
            sought, and if the argument is False, then a 'ko' file is sought.

        use_org_map : bool, False
            If True and the org attribute of the pathway object is an organism code, such as 'eco'
            for E. coli, then an organism-specific 1x resolution file is used if available locally
            in the KEGG directory or online for download to that directory. If False and the org
            attribute is an organism code, then the 1x 'ko' file is used.

        **kwargs
            Valid kwargs are arguments to a biopython.Bio.Graphics.KGML_vis.KGMLCanvas object. These
            control what is displayed on the map from the KGML file.

            Arguments include the following. See the KGMLCanvas class definition in the source code
            for a full list.
            https://github.com/biopython/biopython/blob/master/Bio/Graphics/KGML_vis.py

            import_imagemap : bool
                By default True. Setting to False prevents the base map image from being rendered
                beneath KGML graphics, which is especially useful for decluttering global maps.

            label_compounds : bool
                By default, Drawer sets to False to reduce clutter. Setting to True labels each
                compound with the name of its Graphics element, a KEGG ID such as a COMPOUND or
                GLYCAN ID. Biopython also labels gene entries when this is True. Organism-specific
                KGML files have gene entries, but KO KGML files do not.

            label_orthologs : bool
                By default, Drawer sets to False for global and overview maps to reduce clutter next
                to reaction arrows and to True for standard maps. Setting to True labels each
                ortholog entry with the name of its Graphics element. In KO KGML, that is the first
                KO of the entry, followed by '...' when the entry has several. A rectangle whose
                glyph drops its label is not labeled ('Graphics.glyph').

            label_reaction_entries : bool
                By default, Drawer sets to False for KO maps and for global and overview maps. Other
                standard maps keep Biopython's default of True. Setting to True labels each entry
                that has a reaction with the name of its Graphics element, not with a KEGG REACTION
                ID. In KO KGML, these entries are ortholog entries, so the label repeats the one
                that 'label_orthologs' draws. Biopython draws it in a pass before the ortholog pass,
                which draws boxes over it. With 'label_orthologs' also True, the PDF holds the label
                twice. In EC and RN KGML, these entries are enzyme and reaction entries, which no
                other argument labels. In organism-specific KGML, they are gene entries, which
                'label_compounds' also labels.

            fontname : str
                KGML label font name, with the default being Helvetica.

            fontsize : float
                KGML label font size. Drawer sets the default to 9 for 1x resolution maps and 18
                for 2x, if the map base image is not provided explicitly by map_filepath. If it is
                provided explicitly, then the default is 9, erring on the side of fitting the text
                in an ortholog box on a standard 1x map.
        """
        is_output_file_writable(output_filepath, ok_if_exists=self.overwrite_output)

        # These canvas parameters apply to both standard and global/overview maps.
        if kwargs.get('import_imagemap') is None:
            kwargs['import_imagemap'] = True
        if kwargs.get('label_compounds') is None:
            kwargs['label_compounds'] = False
        # In KO KGML, the ortholog pass already labels every entry that has a reaction.
        if pathway.org == 'ko' and kwargs.get('label_reaction_entries') is None:
            kwargs['label_reaction_entries'] = False
        if map_filepath is not None and kwargs.get('fontsize') is None:
            kwargs['fontsize'] = 9

        if pathway.is_global_map:
            self._draw_global_map(
                pathway,
                output_filepath,
                map_filepath=map_filepath,
                use_org_map=use_org_map,
                **kwargs
            )
        elif pathway.is_overview_map:
            self._draw_overview_map(
                pathway,
                output_filepath,
                map_filepath=map_filepath,
                use_org_map=use_org_map,
                **kwargs
            )
        else:
            self._draw_standard_map(
                pathway,
                output_filepath,
                map_filepath=map_filepath,
                use_org_map=use_org_map,
                **kwargs
            )

    @staticmethod
    def _get_entry_glyphs(
        pathway: Pathway
    ) -> Dict[int, List[Union[Dict[str, Union[str, list]], None]]]:
        """
        Get the glyphs of the entries that have any.

        Parameters
        ==========
        pathway : Pathway
            Object representation of a KGML file.

        Returns
        =======
        Dict[int, List[Union[Dict[str, Union[str, list]], None]]]
            Keys are the IDs of entries with any Graphics element that has a glyph. Each value holds
            the 'glyph' attribute of each Graphics element of the entry, in order. The dictionary is
            empty if no entry has a glyph.
        """
        entry_glyphs: Dict[int, List[Union[Dict[str, Union[str, list]], None]]] = {}
        for entry_uuid in pathway.children['entry']:
            entry: Entry = pathway.uuid_element_lookup[entry_uuid]
            graphics_glyphs = [
                pathway.uuid_element_lookup[graphics_uuid].glyph
                for graphics_uuid in entry.children['graphics']
            ]
            if any(graphics_glyphs):
                entry_glyphs[int(entry.id)] = graphics_glyphs
        return entry_glyphs

    @staticmethod
    def _get_entry_slices(
        pathway: Pathway
    ) -> Dict[int, List[Union[List[Tuple[str, float]], None]]]:
        """
        Get the slices of the entries that have any.

        Parameters
        ==========
        pathway : Pathway
            Object representation of a KGML file.

        Returns
        =======
        Dict[int, List[Union[List[Tuple[str, float]], None]]]
            Keys are the IDs of entries with any sliced Graphics element. Each value holds the
            'slices' attribute of each Graphics element of the entry, in order. The dictionary is
            empty if no entry has slices.
        """
        entry_slices: Dict[int, List[Union[List[Tuple[str, float]], None]]] = {}
        for entry_uuid in pathway.children['entry']:
            entry: Entry = pathway.uuid_element_lookup[entry_uuid]
            graphics_slices = [
                pathway.uuid_element_lookup[graphics_uuid].slices
                for graphics_uuid in entry.children['graphics']
            ]
            if any(graphics_slices):
                entry_slices[int(entry.id)] = graphics_slices
        return entry_slices

    @staticmethod
    def _get_entry_borders(
        pathway: Pathway
    ) -> Dict[int, List[Union[str, Tuple[str, list], None]]]:
        """
        Get the rings of the entries that have any.

        Parameters
        ==========
        pathway : Pathway
            Object representation of a KGML file.

        Returns
        =======
        Dict[int, List[Union[str, Tuple[str, list], None]]]
            Keys are the IDs of entries with any Graphics element that has a ring. Each value holds
            the 'border' attribute of each Graphics element of the entry, in order. The dictionary
            is empty if no entry has a ring.
        """
        entry_borders: Dict[int, List[Union[str, Tuple[str, list], None]]] = {}
        for entry_uuid in pathway.children['entry']:
            entry: Entry = pathway.uuid_element_lookup[entry_uuid]
            graphics_borders = [
                pathway.uuid_element_lookup[graphics_uuid].border
                for graphics_uuid in entry.children['graphics']
            ]
            if any(graphics_borders):
                entry_borders[int(entry.id)] = graphics_borders
        return entry_borders

    @staticmethod
    def _unlabel_tiny_graphics(bio_pathway: KGML_pathway.Pathway) -> None:
        """
        Blank the names of boxes and circles too small to see, so that they are not labeled.

        A shape no wider and no taller than 'TINY_GRAPHICS_MAX_SIZE' gets an empty name. Lines are
        left alone, since their width is the width of the stroke. Only Biopython's copy of the
        pathway changes, so the anvi'o Pathway keeps every name.

        Parameters
        ==========
        bio_pathway : Bio.KEGG.KGML.KGML_pathway.Pathway
            Biopython representation of a KGML file.
        """
        for entry in bio_pathway.entries.values():
            for graphics in entry.graphics:
                if graphics.type == 'line':
                    continue
                if (
                    graphics.width <= TINY_GRAPHICS_MAX_SIZE
                    and graphics.height <= TINY_GRAPHICS_MAX_SIZE
                ):
                    graphics.name = ''

    def _draw_global_map(
        self,
        pathway: Pathway,
        output_filepath: str,
        map_filepath: str = None,
        use_org_map: bool = False,
        **kwargs
    ) -> None:
        """
        Draw a global pathway map with KGML data as a PDF file.

        Parameters
        ==========
        pathway : Pathway
            Object representation of a KGML file.

        output_filepath : str
            Path to PDF output file containing the pathway map.

        map_filepath : str, None
            Path to pathway map image file to use as the base image of the output file. If None,
            then a 1x resolution PNG file stored in the KEGG directory is used; it is assumed the
            KGML data is scaled to fit this image size. The org attribute of the pathway object is
            used to seek the corresponding PNG file, i.e., a 'ko' pathway containing reactions with
            KO IDs results in the reference 'ko' map being sought, whereas an 'ec' pathway
            containing reactions with EC number IDs corresponds to the reference 'ec' map.

        use_org_map : bool, False
            If True and the org attribute of the pathway object is an organism code, such as 'eco'
            for E. coli, then an organism-specific 1x resolution file is used if available locally
            in the KEGG directory or available online for download to that directory. If False and
            the org attribute is an organism code, then the 1x 'ko' file is used.

        **kwargs
            Valid kwargs are arguments to a biopython.Bio.Graphics.KGML_vis.KGMLCanvas object.
            These control what is displayed on the map from the KGML file.
        """
        if kwargs.get('label_orthologs') is None:
            kwargs['label_orthologs'] = False
        if kwargs.get('label_reaction_entries') is None:
            kwargs['label_reaction_entries'] = False
        if kwargs.get('fontsize') is None:
            kwargs['fontsize'] = 9

        # Glyphs, slices, and rings are not KGML, so they are gathered before the KGML round trip.
        entry_glyphs = self._get_entry_glyphs(pathway)
        entry_slices = self._get_entry_slices(pathway)
        entry_borders = self._get_entry_borders(pathway)
        bio_pathway = KGML_parser.read(StringIO(self.xml_ops.get_str(pathway)))
        self._unlabel_tiny_graphics(bio_pathway)

        if map_filepath is None:
            if pathway.org == 'ko':
                map_filepath = os.path.join(
                    self.kegg_context.png_1x_ko_dir, f'{pathway.org}{pathway.number}.png'
                )
            elif pathway.org == 'ec':
                map_filepath = os.path.join(
                    self.kegg_context.png_1x_ec_dir, f'ec{pathway.number}.png'
                )
            elif pathway.org == 'rn':
                map_filepath = os.path.join(
                    self.kegg_context.png_1x_rn_dir, f'rn{pathway.number}.png'
                )
            elif use_org_map:
                map_filepath = os.path.join(
                    self.kegg_context.png_1x_org_dir, f'{pathway.org}{pathway.number}.png'
                )
            else:
                map_filepath = os.path.join(
                    self.kegg_context.png_1x_ko_dir, f'ko{pathway.number}.png'
                )
        else:
            assert not use_org_map
            is_file_exists(map_filepath)

        if use_org_map and not is_file_exists(map_filepath, dont_raise=True):
            download_org_pathway_image_files(f'{pathway.org}{pathway.number}', self.kegg_dir)

        bio_pathway.image = map_filepath

        if entry_glyphs or entry_slices or entry_borders:
            canvas = SlicedKGMLCanvas(
                bio_pathway, entry_slices, entry_borders, entry_glyphs=entry_glyphs, **kwargs
            )
        else:
            canvas = KGMLCanvas(bio_pathway, **kwargs)
        canvas.non_reactant_transparency = self.non_reactant_transparency
        canvas.draw(output_filepath)

    def _draw_overview_map(
        self,
        pathway: Pathway,
        output_filepath: str,
        map_filepath: str = None,
        use_org_map: bool = False,
        **kwargs
    ) -> None:
        """
        Draw an overview pathway map with KGML data as a PDF file.

        Parameters
        ==========
        pathway : Pathway
            Object representation of a KGML file.

        output_filepath : str
            Path to PDF output file containing the pathway map.

        map_filepath : str, None
            Path to pathway map image file to use as the base image of the output file. If None,
            then a 2x resolution 'map' PNG file stored in the KEGG directory is used; it is assumed
            the KGML data is scaled to fit this image size.

        use_org_map : bool, False
            If True and the org attribute of the pathway object is an organism code, such as 'eco'
            for E. coli, then an organism-specific 1x resolution file is used if available locally
            in the KEGG directory or available online for download to that directory. If False and
            the org attribute is an organism code, then the 2x 'ko' file is used.

        **kwargs
            Valid kwargs are arguments to a biopython.Bio.Graphics.KGML_vis.KGMLCanvas object.
            These control what is displayed on the map from the KGML file.
        """
        if kwargs.get('label_orthologs') is None:
            kwargs['label_orthologs'] = False
        if kwargs.get('label_reaction_entries') is None:
            kwargs['label_reaction_entries'] = False

        # Glyphs, slices, and rings are not KGML, so they are gathered before the KGML round trip.
        entry_glyphs = self._get_entry_glyphs(pathway)
        entry_slices = self._get_entry_slices(pathway)
        entry_borders = self._get_entry_borders(pathway)
        bio_pathway = KGML_parser.read(StringIO(self.xml_ops.get_str(pathway)))
        self._unlabel_tiny_graphics(bio_pathway)

        if map_filepath is None:
            if use_org_map:
                map_filepath = os.path.join(
                    self.kegg_context.png_1x_org_dir, f'{pathway.org}{pathway.number}.png'
                )
                if kwargs.get('fontsize') is None:
                    kwargs['fontsize'] = 9
            else:
                map_dir = self.kegg_context.png_2x_map_dir
                map_filepath = os.path.join(
                    self.kegg_context.png_2x_map_dir, f'map{pathway.number}.png'
                )
                if kwargs.get('fontsize') is None:
                    kwargs['fontsize'] = 18
        else:
            assert not use_org_map
            is_file_exists(map_filepath)
            if kwargs.get('fontsize') is None:
                kwargs['fontsize'] = 9

        if use_org_map and not is_file_exists(map_filepath, dont_raise=True):
            download_org_pathway_image_files(f'{pathway.org}{pathway.number}', self.kegg_dir)

        bio_pathway.image = map_filepath

        if entry_glyphs or entry_slices or entry_borders:
            canvas = SlicedKGMLCanvas(
                bio_pathway, entry_slices, entry_borders, entry_glyphs=entry_glyphs, **kwargs
            )
        else:
            canvas = KGMLCanvas(bio_pathway, **kwargs)
        canvas.non_reactant_transparency = self.non_reactant_transparency
        canvas.draw(output_filepath)

    def _draw_standard_map(
        self,
        pathway: Pathway,
        output_filepath: str,
        map_filepath: str = None,
        use_org_map: bool = False,
        **kwargs
    ) -> None:
        """
        Draw a standard (not global/overview) pathway map with KGML data as a PDF file.

        Parameters
        ==========
        pathway : Pathway
            Object representation of a KGML file.

        output_filepath : str
            Path to PDF output file containing the pathway map.

        map_filepath : str, None
            Path to pathway map image file to use as the base image of the output file. If None,
            then the 2x resolution 'map' PNG file stored in the KEGG directory is used; it is
            assumed the KGML data is scaled to fit this image size.

        use_org_map : bool, False
            If True and the org attribute of the pathway object is an organism code, such as 'eco'
            for E. coli, then an organism-specific 1x resolution file is used if available locally
            in the KEGG directory or available online for download to that directory. If False and
            the org attribute is an organism code, then the 2x 'ko' file is used.

        **kwargs
            Valid kwargs are arguments to a biopython.Bio.Graphics.KGML_vis.KGMLCanvas object.
            These control what is displayed on the map from the KGML file.
        """
        if kwargs.get('label_orthologs') is None:
            kwargs['label_orthologs'] = True

        # Glyphs, slices, and rings are not KGML, so they are gathered before the KGML round trip.
        entry_glyphs = self._get_entry_glyphs(pathway)
        entry_slices = self._get_entry_slices(pathway)
        entry_borders = self._get_entry_borders(pathway)
        bio_pathway = KGML_parser.read(StringIO(self.xml_ops.get_str(pathway)))
        self._unlabel_tiny_graphics(bio_pathway)

        if map_filepath is None:
            if use_org_map:
                map_filepath = os.path.join(
                    self.kegg_context.png_1x_org_dir, f'{pathway.org}{pathway.number}.png'
                )
                if kwargs.get('fontsize') is None:
                    kwargs['fontsize'] = 9
            else:
                map_filepath = os.path.join(
                    self.kegg_context.png_2x_map_dir, f'map{pathway.number}.png'
                )
                if kwargs.get('fontsize') is None:
                    kwargs['fontsize'] = 18
        else:
            assert not use_org_map
            is_file_exists(map_filepath)
            if kwargs.get('fontsize') is None:
                kwargs['fontsize'] = 9

        if use_org_map and not is_file_exists(map_filepath, dont_raise=True):
            download_org_pathway_image_files(f'{pathway.org}{pathway.number}', self.kegg_dir)

        bio_pathway.image = map_filepath

        if entry_glyphs or entry_slices or entry_borders:
            canvas = SlicedKGMLCanvas(
                bio_pathway, entry_slices, entry_borders, entry_glyphs=entry_glyphs, **kwargs
            )
        else:
            canvas = KGMLCanvas(bio_pathway, **kwargs)
        canvas.non_reactant_transparency = self.non_reactant_transparency
        canvas.draw(output_filepath)

class Tester:
    """
    Tests KGML operations.

    Attributes
    ==========
    xml_ops : XMLOps()
        Loads KMGL files.

    run : anvio.terminal.Run
        This object prints run information to the terminal.

    progress : anvio.terminal.Progress
        This object prints transient progress information to the terminal.
    """
    def __init__(
        self,
        run: terminal.Run = terminal.Run(),
        progress: terminal.Progress = terminal.Progress(),
    ) -> None:
        """
        Parameters
        ==========
        run : anvio.terminal.Run, anvio.terminal.Run()
            This object prints run information to the terminal.

        progress : anvio.terminal.Progress, anvio.terminal.Progress()
            This object prints transient progress information to the terminal.
        """
        self.xml_ops = XMLOps()

        self.run = run
        self.progress = progress

    def load_all_anvio_kgml_files(self, kegg_dirpath: str = None) -> None:
        """
        Load each KGML file within a superdirectory formatted like an anvi'o KEGG directory into
        memory as a KGML pathway object, and test that the object can be converted back to a string
        equivalent to the contents of the file.

        Parameters
        ==========
        kegg_dirpath : str, None
            A directory of KEGG files like that installed by `anvi-setup-kegg-data`. By default, the
            default anvi'o KEGG installation location is sought.

        Returns
        =======
        None
        """
        args = Namespace()
        if kegg_dirpath is not None:
            args.kegg_data_dir = kegg_dirpath
        kegg_context = KeggContext(args)
        for dirname in os.listdir(kegg_context.map_image_kgml_dir):
            dirpath = os.path.join(kegg_context.map_image_kgml_dir, dirname)
            if not os.path.isdir(dirpath):
                continue
            self.load_kgml_files_in_dir(dirpath)

    def load_kgml_files_in_dir(self, dirpath: str) -> None:
        """
        Load each KGML file in a single directory into memory as a KGML pathway object, and test
        that the object can be converted back to a string equivalent to the contents of the file.

        Parameters
        ==========
        dirpath : str
            Path to a directory in which to look for KGML files, assumed to have a '.xml' extension.

        Returns
        =======
        None
        """
        filepaths: List[str] = []
        for filename in os.listdir(dirpath):
            filepath = os.path.join(dirpath, filename)
            if not os.path.isfile(filepath) or not os.path.splitext(filepath)[1] == '.xml':
                continue
            filepaths.append(filepath)
        self.progress.new(
            f"Testing KGML files in {os.path.dirname(dirpath)}",
            progress_total_items=len(filepaths)
        )
        for filepath in filepaths:
            self.progress.update("...", increment=True)
            self.load_kgml_file(filepath)
        self.progress.end()
        self.run.info_single(f"Tested {len(filepaths)} KGML files in '{dirpath}'")

    def load_kgml_file(self, filepath: str, buffer: int = 100) -> None:
        """
        Load a KGML file into memory as a KGML pathway object, and test that the object can be
        converted back to a string equivalent to the contents of the file.

        Parameters
        ==========
        filepath : str
            Path to a KGML file.

        buffer : int, 100
            If an inconsistency is found between the string representing the contents of the KGML
            file and the string representing a KGML file reconstructed from the loaded pathway
            object, then display the text where these strings diverge, including the number of
            characters given by `buffer` around this point.

        Returns
        =======
        None
        """
        pathway = self.xml_ops.load(filepath)
        reconstructed_xml_str = self.xml_ops.get_str(pathway)
        with open(filepath) as file:
            xml_str = file.read()
        if xml_str != reconstructed_xml_str:
            for i in range(1, min(len(xml_str), len(reconstructed_xml_str))):
                if xml_str[: i] == reconstructed_xml_str[: i]:
                    continue
                error_message = (
                    "Anvi'o loaded a KGML file as a kgml.Pathway object. It then tried to convert "
                    "the object back into a string equivalent to the text of the KGML file, and "
                    "failed. Here is the area of the KGML text where an inconsistency was "
                    "detected. First, the file text is displayed, and then the text reconstructed "
                    "from the object is displayed.\n"
                )
                error_message += xml_str[max(0, i - buffer): min(i + buffer, len(xml_str))] + "\n"
                error_message += reconstructed_xml_str[
                    max(0, i - buffer): min(i + buffer, len(reconstructed_xml_str))
                ]
                raise AssertionError(error_message)
