#!/usr/bin/env python
# -*- coding: utf-8
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
import uuid
import numpy as np
import xml.etree.ElementTree as ET

from io import StringIO
from copy import deepcopy
from argparse import Namespace
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas
from matplotlib.colors import Colormap, rgb2hex

from typing import Iterable, Literal, NewType, Union

import anvio.kegg as kegg
import anvio.terminal as terminal

from anvio.errors import ConfigError
from anvio import FORCE_OVERWRITE, __version__ as VERSION
from anvio.filesnpaths import is_file_exists, is_output_file_writable

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2024, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = VERSION
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"
__status__ = "Development"

class Element:
    """
    Representation of an XML element from a KGML file.

    Attributes
    ==========
    uuid : str
        Unique ID assigned by anvi'o, not KEGG. Since Element objects do not directly reference each
        other, UUIDs are used to record element-subelement relationships.
    """
    # Subclass names are the same as the capitalized tag attribute.
    tag: str
    # Element attributes are required or not, according to the KGML schema.
    attribute_required: dict[str, bool]

    def __init__(self) -> None:
        self.uuid = str(uuid.uuid4())

Priority = NewType('Priority', float)
Priority.__doc__ = "Priorities are used to reorder graphics along the z-axis in pathway maps."

GraphicsColor = NewType('GraphicsColor', tuple[str, str])
GraphicsColor.__doc__ = "Graphics element foreground and background color hex codes, respectively."

class PathwayColorPriority:
    """
    Used to reorder graphics in a pathway map by their colors.

    Pathway Entry Graphics colors are assigned Priority values.

    Attributes
    ==========
    entry_type : dict[str, EntryColorPriority], {}
        Keys are Entry type strings, any of the possible values of the type attribute of the Entry
        class, e.g., 'ortholog', 'reaction', 'compound'. Values are EntryColorPriority objects.
    """
    def __init__(self) -> None:
        self.entry_type: dict[str, EntryColorPriority] = {}

class EntryColorPriority:
    """
    Used to reorder graphics for an entry type by their colors.

    Entry Graphics colors are assigned Priority values.

    Attributes
    ==========
    graphics_type : dict[str, GraphicsColorPriority], {}
        Keys are Graphics type strings, any of the possible values of the type attribute of the Graphics
        class, e.g., 'rectangle', 'line'. Values are GraphicsColorPriority objects.
    """
    def __init__(self) -> None:
        self.graphics_type: dict[str, GraphicsColorPriority] = {}

class GraphicsColorPriority:
    """
    Used to reorder graphics by their colors.

    Graphics colors (unique combinations of foreground and background colors) are assigned Priority
    values.

    Attributes
    ==========
    color : dict[GraphicsColor, Priority], {}
        Keys are GraphicsColor tuples. Values are Priority floats: larger numbers indicate higher
        priority graphics given the foreground and background colors.
    """
    def __init__(self) -> None:
        self.color: dict[GraphicsColor, Priority] = {}

class PathwayUnprioritizedColor:
    """
    Used to set colors of unprioritized graphics that are drawn behind prioritized graphics.

    Unprioritized colors are set for each type of Entry Graphics element.

    Attributes
    ==========
    entry_type : dict[str, EntryUnprioritizedColor], {}
        Keys are Entry type strings, any of the possible values of the type attribute of the Entry
        class, e.g., 'ortholog', 'reaction', 'compound'. Values are EntryUnprioritizedColor objects.
    """
    def __init__(self) -> None:
        self.entry_type: dict[str, EntryUnprioritizedColor] = {}

class EntryUnprioritizedColor:
    """
    Used to set colors of unprioritized graphics of an entry type that are drawn behind prioritized
    graphics.

    Unprioritized colors are set for each type of Graphics element.

    Attributes
    ==========
    graphics_type : dict[str, GraphicsColor], {}
        Keys are Graphics type strings, any of the possible values of the type attribute of the Graphics
        class, e.g., 'rectangle', 'line'. Values are GraphicsColor tuples.
    """
    def __init__(self) -> None:
        self.graphics_type: dict[str, GraphicsColor] = {}

class PathwayThicknessPriority:
    """
    Used to reorder line graphics in a global or overview pathway map by line thickness.

    Line thicknesses are assigned Priority values.

    Attributes
    ==========
    width : dict[float, Priority], {}
        Keys are width floats, which are compared to values of the width attribute of Graphics
        elements of type 'line'. Values are Priority floats: larger numbers indicate higher priority
        graphics given the width.
    """
    def __init__(self) -> None:
        self.width: dict[float, Priority] = {}

class Pathway(Element):
    """
    The single top-level element in a KGML file.

    Entry subelements representing reactions and compounds, among other items, are shown in maps.

    The order of entries affects the order in which they are drawn, with entries occurring first in
    the subelements attribute being drawn under entries occurring last. The z-axis position is
    especially important in global maps, which have numerous crossing reaction lines. Entries can be
    assigned priorities to control how they are overlaid, with higher-priority entries overlaying
    lower-priority entries. Priorities can be set from other properties, such as the color of Entry
    Graphics, which can represent data such as presence/absence, abundances, and fluxes. This class
    contains methods related to the colors and priorities of entries.

    Attributes
    ==========
    subelement_tags : tuple[str]
        Tags of possible Pathway subelements.

    name : str, None
        KEGG ID of the pathway map, e.g., 'ko00010', 'ec01100', 'eco00010'.

    org : str, None,
        ko/ec/rn/[org prefix] in ID, e.g., 'ko', 'ec', 'eco'.

    number : str, None
        Map number in ID, e.g., '00010', '001100'.

    title : str, None
        Map title, e.g., 'Glycolysis / Gluconeogenesis', 'Metabolic pathways'.

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

    subelements : dict[str, list[str]]
        Keys are subelement tags, values are lists of subelement UUIDs.

    uuid_element_lookup : dict[str, Element], {}
        Keys are UUIDs of all elements in the pathway, values are Element objects.

    kegg_id_element_lookup : dict[str, list[Element]] = {}
        Keys are KEGG IDs (e.g., 'K00001', 'R00010'), values are lists of Element objects with the
        KEGG ID in the name attribute. A KEGG ID is not necessarily unique to an element, and an
        element can have multiple KEGG IDs.

    is_global_map : bool, None
        True if the pathway map is a global map, as indicated by the map number, rather than a
        standard or overview map.

    is_overview_map : bool, None
        True if the pathway map is an overview map, as indicated by the map number, rather than a
        standard or global map.

    color_priority : PathwayColorPriority, PathwayColorPriority()
        Defines the order of Entry Graphics by combination of foreground and background colors. The
        order is set along with the attribute value using the method, set_priority.

    unprioritized_color : PathwayUnprioritizedColor, PathwayUnprioritizedColor()
        Combinations of foreground and background colors of unprioritized Entry Graphics drawn below
        prioritized entries. The attribute value is set using the method, set_priority, when a
        recolor_unprioritized_entries argument is provided.
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
    subelement_tags: tuple[str] = (
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

        self.subelements: dict[str, list[str]] = {tag: [] for tag in self.subelement_tags}

        self.uuid_element_lookup: dict[str, Element] = {}
        self.kegg_id_element_lookup: dict[str, list[Element]] = {}

        self._is_global_map: bool = None
        self._is_overview_map: bool = None

        self.color_priority: PathwayColorPriority = PathwayColorPriority()
        self.unprioritized_color: PathwayUnprioritizedColor = PathwayUnprioritizedColor()
        self.thickness_priority: PathwayThicknessPriority = PathwayThicknessPriority()
        self.unprioritized_width: float = None

        super().__init__()

    @property
    def is_global_map(self):
        if self.number is None:
            return None
        return True if re.match(kegg.GLOBAL_MAP_ID_PATTERN, self.number) else False

    @property
    def is_overview_map(self):
        if self.number is None:
            return None
        return True if re.match(kegg.OVERVIEW_MAP_ID_PATTERN, self.number) else False

    def set_priority(
        self,
        new_priority: Union[
            PathwayColorPriority,
            PathwayThicknessPriority,
            tuple[PathwayColorPriority, PathwayThicknessPriority],
            tuple[PathwayThicknessPriority, PathwayColorPriority]
        ],
        recolor_unprioritized_entries: Union[str, PathwayUnprioritizedColor] = None,
        color_associated_compounds: Literal['high', 'low', 'average'] = None,
        colormap: Colormap = None,
        reset_unprioritized_thickness: Union[bool, float] = False
    ) -> None:
        """
        Set the color_priority and/or thickness_priority attributes. Entry elements in the
        subelements attribute are automatically reordered.

        Regarding color priorities, a single Entry (e.g., representing KOs, ECs, RNs, compounds) can
        occur multiple times on a map (e.g., as different rectangles or circles), and thus have
        multiple Graphics elements. It is required here that Graphics elements of the same type
        (e.g., rectangle type Graphics or line type Graphics) for an Entry must all have the same
        foreground and background colors if they are to be ordered.

        Thickness priorities apply to line Graphics of reaction-like entries (KO, EC, RN, gene,
        group) in global and overview pathways; thickness priorities cannot be set for standard
        pathways.

        Entries with higher priorities are placed last in the subelements attribute and in KGML
        files, and they are rendered in the foreground of the map. The lowest Priority entries are
        always those without fg/bg colors defined in the color_priority attribute or line widths
        defined in the thickness_priority attribute. These entries are placed first in the
        subelements attribute and KGML files, and they are rendered in the background of the map and
        thus can be overlaid by higher Priority entries.

        If both color and thickness priorities are set, one takes precedence over the other given
        the new_priority argument. For example, with color taking precedence over thickness, if an
        Entry has a lower Priority line width but higher Priority fg/bg colors than another Entry,
        then the former will be placed before the latter in the subelements attribute and in KGML
        files, and the former will be rendered above the latter on the map.

        Attributes
        ==========
        new_priority : Union[
            PathwayColorPriority,
            PathwayThicknessPriority,
            tuple[PathwayColorPriority, PathwayThicknessPriority],
            tuple[PathwayThicknessPriority, PathwayColorPriority]
        ]
            Priority objects used to set the color_priority and/or thickness_priority attributes.
            Entry elements in the subelements attribute are reordered accordingly.

            Only global and overview maps can have line entries reordered by thickness. If entries
            are prioritized by both color and thickness, a tuple is provided as an argument, and the
            order of the priority objects in the tuple determines which type of priority takes
            precedence, with the first priority object in the tuple taking precedence over the
            second.

            Here is what is actually set as color_priority. A deep copy is made of the
            PathwayColorPriority object. GraphicsColorPriority objects (the innermost entries in the
            PathwayColorPriority dictionary structure) are reordered by Priority value ascending, so
            that the lowest Priority colors occur first for each Entry Graphics type. Otherwise, the
            dictionary structure is not altered: if, for example, 'ortholog' appears before
            'compound' as a key in PathwayColorPriority, then ortholog entries will occur before
            compound entries in a KGML file output from the Pathway object, and compound graphics
            can overlay ortholog graphics in the map drawn from the Pathway object.

            Here is what is actually set as thickness_priority. A deep copy is made of the
            PathwayThicknessPriority object. The width dictionary attribute of the object is
            reordered by Priority value ascending, so that the lowest Priority widths occur first.

        recolor_unprioritized_entries : Union[str, PathwayUnprioritizedColor], None
            Recolor unprioritized entries, either automatically with a string argument, or with a
            custom dictionary argument for fine-tuning foreground and background colors by Entry
            type. The valid string arguments for automatic recoloring are 'w' for black and white
            and 'g' for gray. The unprioritized_color attribute is set given the argument value,
            either automatically with the string argument or directly from the custom
            PathwayUnprioritizedColor.

        color_associated_compounds : Literal['high', 'low', 'average'], None
            Automatically set the background color of compound entries based on the color Priority
            of reaction-like entries involving the compounds. By default, compounds participating in
            reactions are circles, and reactions are lines on global/overview maps and boxes or
            lines on standard maps.

            An argument of 'high' or 'low' sets the compound background color to the bg color of the
            reaction with the highest or lowest Priority fg/bg color combination.

            'average' sets the bg color to the average bg color of reactions with prioritized
            colors: unprioritized reactions are not taken into account. 'average' should only be
            used if Priority values are normalized to the interval [0, 1] and can be converted to a
            color given by the colormap argument. The average Priority value of the reactions with
            prioritized colors is mapped to a bg color for compound Entry circle Graphics.

            Automatically colored compound entries are added to the color_priority attribute, and
            Entry elements in the subelements attribute are reordered accordingly. Compound entries
            that are already in the color_priority attribute are exempt from recoloring and given
            higher Priority than automatically recolored compound entries.

        colormap : Colormap, None
            If 'average' is used as the color_associated_compounds argument, a colormap must be
            provided to map averaged Priority values on the interval [0, 1] to a background color
            for compound Entry circle Graphics.

        reset_unprioritized_thickness : Union[bool, float], False
            If not False, set unprioritized reaction-like Entry line Graphics to a uniform width.
            With an argument value of True, a width of 6.0 is used for global maps and 1.0 for
            overview maps. Alternatively, a custom width can be specified with a float value.
        """
        # Determine which graphical property takes precedence.
        if isinstance(new_priority, PathwayColorPriority):
            precedence = {'color': new_priority}
        elif isinstance(new_priority, PathwayThicknessPriority):
            precedence = {'thickness': new_priority}
        elif isinstance(new_priority, tuple):
            if (
                isinstance(new_priority[0], PathwayColorPriority) and
                isinstance(new_priority[1], PathwayThicknessPriority)
            ):
                precedence = {'color': new_priority[0], 'thickness': new_priority[1]}
            elif (
                isinstance(new_priority[0], PathwayThicknessPriority) and
                isinstance(new_priority[1], PathwayColorPriority)
            ):
                precedence = {'thickness': new_priority[0], 'color': new_priority[1]}
            else:
                raise ConfigError(
                    "A tuple value of the 'new_priority' argument must be a "
                    "(<PathwayColorPriority>, <PathwayThicknessPriority>) or a "
                    "(<PathwayThicknessPriority>, <PathwayColorPriority>)."
                )
        else:
            raise ConfigError(
                "Valid 'new_priority' argument values are a <PathwayColorPriority>, "
                "<PathwayThicknessPriority>, (<PathwayColorPriority>, <PathwayThicknessPriority>), "
                "or (<PathwayThicknessPriority>, <PathwayColorPriority>)."
            )

        if 'thickness' in precedence:
            if not (self.is_global_map or self.is_overview_map):
                raise ConfigError(
                    "Only global and overview maps can have line entries reordered by thickness."
                )

        for prioritized_property in reversed(precedence):
            if prioritized_property == 'color':
                self.set_color_priority(
                    precedence['color'],
                    recolor_unprioritized_entries=recolor_unprioritized_entries,
                    color_associated_compounds=color_associated_compounds,
                    colormap=colormap
                )
            elif prioritized_property == 'thickness':
                self.set_thickness_priority(
                    precedence['thickness'],
                    reset_unprioritized_thickness=reset_unprioritized_thickness
                )
            else:
                raise AssertionError

    def set_color_priority(
        self,
        new_color_priority: PathwayColorPriority,
        recolor_unprioritized_entries: Union[str, PathwayUnprioritizedColor] = False,
        color_associated_compounds: Literal['high', 'low', 'average'] = None,
        colormap: Colormap = None
    ) -> None:
        """
        Set the color_priority attribute. Entry elements in the subelements attribute are
        automatically reordered.

        A single Entry (e.g., representing KOs, ECs, KEGG reactions, or compounds) can occur
        multiple times on a map (e.g., as different rectangles or circles), and thus have multiple
        Graphics elements. It is required here that Graphics elements of the same type (e.g.,
        rectangle type Graphics or line type Graphics) for an Entry must all have the same
        foreground and background colors if they are to be ordered.

        Entries with higher priority fg/bg colors are placed last in the subelements attribute and
        in KGML files, and they are rendered in the foreground of the map. The lowest priority
        entries are always those without fg/bg colors defined in the color_priority attribute; these
        entries are placed first in the subelements attribute and KGML files, and they are rendered
        in the background of the map and thus can be overlaid by higher priority entries.

        Parameters
        ==========
        new_color_priority : PathwayColorPriority
            This dictionary, a PathwayColorPriority object, is used to set the color_priority
            attribute. Entry elements in the subelements attribute are reordered accordingly.

            Here is what is actually set as color_priority. A deep copy is made of the argument.
            GraphicsColorPriority objects (the innermost entries in the PathwayColorPriority
            dictionary structure) are reordered by Priority value ascending, so that the lowest
            Priority colors occur first for each Entry Graphics type. Otherwise, the dictionary
            structure is not altered: if, for example, 'ortholog' appears before 'compound' as a key
            in PathwayColorPriority, then ortholog entries will occur before compound entries in a
            KGML file output from the Pathway object, and compound graphics can overlay ortholog
            graphics in the map drawn from the Pathway object.

        recolor_unprioritized_entries : Union[str, PathwayUnprioritizedColor], False
            Recolor unprioritized entries, either automatically with a string argument, or with a
            custom dictionary argument for fine-tuning foreground and background colors by Entry
            type. The valid string arguments for automatic recoloring are 'w' and 'g'.

            It is assumed that global maps contain lines for reaction-like entries and circles for
            compounds, so automatic recoloring affects reaction-like Entry line Graphics and
            compound Entry circle Graphics. 'w' erases unprioritized lines and circles by coloring
            them entirely white. 'g' colors them a light gray (#E0E0E0), consistent with other
            "unidentified" reactions in the base map.

            It is assumed that overview maps contain reaction lines (drawn as arrows) and compound
            circles. Unlike global maps, 'w' colors unprioritized arrows black, consistent with
            other "unidentified" reactions in the base map. 'w' colors unprioritized circles white
            (with a black border). 'g' colors arrows and circles light gray.

            It is assumed that standard maps contain reaction boxes or lines and compound circles.
            'w' colors unprioritized boxes white, with black text; lines black; and circles white.
            'g' colors unprioritized boxes light gray, with black text; lines light gray; and
            circles light gray. Compound rectangles are also found in a small number of KGML files
            (see 00121, 00621, 01052, 01054), and it is best if the size of these is set to zero so
            that they don't obscure the chemical structure drawings to which they correspond on maps
            (see `anvio.keggmapping.Mapper._zero_out_compound_rectangles`).

            A custom dictionary argument, a PathwayUnprioritizedColor object, can be used to set
            fg/bg colors in detail.

            The unprioritized_color attribute is set given the argument value, either automatically
            with the string ('w', 'g') argument or directly from the custom
            PathwayUnprioritizedColor.

        color_associated_compounds : Literal['high', 'low', 'average'], None
            Automatically set the background color of compound entries based on the color Priority
            of reaction-like entries involving the compounds. By default, compounds participating in
            reactions are circles, and reactions are lines on global/overview maps and boxes or
            lines on standard maps.

            An argument of 'high' or 'low' sets the compound background color to the bg color of the
            reaction with the highest or lowest Priority fg/bg color combination.

            'average' sets the bg color to the average bg color of reactions with prioritized
            colors: unprioritized reactions are not taken into account. 'average' should only be
            used if Priority values are normalized to the interval [0, 1] and can be converted to a
            color given by the colormap argument. The average Priority value of the reactions with
            prioritized colors is mapped to a bg color for compound Entry circle Graphics.

            Automatically colored compound entries are added to the color_priority attribute, and
            Entry elements in the subelements attribute are reordered accordingly. Compound entries
            that are already in the color_priority attribute are exempt from recoloring and given
            higher Priority than automatically recolored compound entries.

        colormap : matplotlib.colors.Colormap, None
            If 'average' is used as the color_associated_compounds argument, a colormap must be
            provided to map averaged Priority values on the interval [0, 1] to a background color
            for compound Entry circle Graphics.
        """
        # Check that new_color_priority only contains nonnegative priority values.
        for entry_color_priority in new_color_priority.entry_type.values():
            for graphics_color_priority in entry_color_priority.graphics_type.values():
                for priority in graphics_color_priority.color.values():
                    assert priority > 0

        # Make the object assigned to the color_priority attribute, reordering colors from lowest to
        # highest priority.
        color_priority = PathwayColorPriority()
        for entry_type, ecp in new_color_priority.entry_type.items():
            entry_color_priority = EntryColorPriority()
            color_priority.entry_type[entry_type] = entry_color_priority
            for graphics_type, gcp in ecp.graphics_type.items():
                graphics_color_priority = GraphicsColorPriority()
                entry_color_priority.graphics_type[graphics_type] = graphics_color_priority
                for graphics_color, priority in sorted(
                    graphics_color_priority.color.items(), key=lambda item: item[1]
                ):
                    graphics_color_priority[graphics_color] = priority
        self.color_priority = color_priority

        # Reorder Entry elements in the subelements attribute from lowest to highest priority.
        unprioritized_entry_uuids = self.order_entries_by_color_priority()

        if recolor_unprioritized_entries:
            if isinstance(recolor_unprioritized_entries, str):
                assert recolor_unprioritized_entries in ('w', 'g')

                # Recolor reactions that have not been assigned a color Priority.
                if recolor_unprioritized_entries == 'w':
                    if self.is_overview_map:
                        color_hex_code = '#000000'
                    else:
                        color_hex_code = '#FFFFFF'
                elif recolor_unprioritized_entries == 'g':
                    color_hex_code = '#E0E0E0'
                reaction_unprioritized_color = self.recolor_unprioritized_reaction_entries(
                    unprioritized_entry_uuids, color_hex_code
                )

                if self.color_associated_compounds is None:
                    # Recolor compounds that have not been assigned a color Priority. Avoid
                    # redundancy if this is to happen again after coloring associated compounds.
                    if recolor_unprioritized_entries == 'w':
                        color_hex_code = '#FFFFFF'
                    elif recolor_unprioritized_entries == 'g':
                        color_hex_code = '#E0E0E0'
                    compound_unprioritized_color = self.recolor_unprioritized_compound_entries(
                        unprioritized_entry_uuids, color_hex_code
                    )

                    # Set the unprioritized_color attribute before returning.
                    unprioritized_color = deepcopy(reaction_unprioritized_color)
                    unprioritized_color.entry_type['compound'] = deepcopy(
                        compound_unprioritized_color.entry_type['compound']
                    )
                    self.unprioritized_color = unprioritized_color

                    return
            elif isinstance(recolor_unprioritized_entries, PathwayUnprioritizedColor):
                if self.color_associated_compounds is None:
                    self.recolor_unprioritized_entries(
                        unprioritized_entry_uuids, recolor_unprioritized_entries
                    )

                    # Set the unprioritized_color attribute before returning.
                    self.unprioritized_color = deepcopy(recolor_unprioritized_entries)

                    return
                else:
                    # Avoid redundancy in coloring compounds that have not been assigned a color
                    # Priority if this is to happen again after coloring associated compounds.
                    noncompound_unprioritized_color = PathwayColorPriority()
                    for (
                        entry_type, entry_type_unprioritized_color
                    ) in recolor_unprioritized_entries.entry_type.items():
                        if entry_type != 'compound':
                            noncompound_unprioritized_color.entry_type[
                                entry_type
                            ] = entry_type_unprioritized_color
                    self.recolor_unprioritized_entries(
                        unprioritized_entry_uuids, noncompound_unprioritized_color
                    )
            else:
                raise ConfigError(
                    "Valid 'recolor_unprioritized_entries' argument values are the strings, 'w' or "
                    "'g', or a <PathwayUnprioritizedColor>."
                )

        self.color_associated_compounds(color_associated_compounds, colormap=colormap)

        # Reorder compound entries in the subelements attribute according to color priority.
        unprioritized_entry_uuids = self.order_entries_by_color_priority()

        if not recolor_unprioritized_entries:
            return

        # Recolor compounds that have not been assigned a color priority.
        if isinstance(recolor_unprioritized_entries, str):
            if recolor_unprioritized_entries == 'w':
                color_hex_code = '#FFFFFF'
            elif recolor_unprioritized_entries == 'g':
                color_hex_code = '#E0E0E0'
            compound_unprioritized_color = self.recolor_unprioritized_compound_entries(
                unprioritized_entry_uuids, color_hex_code
            )

            # Set the unprioritized_color attribute before returning.
            unprioritized_color = deepcopy(reaction_unprioritized_color)
            unprioritized_color.entry_type['compound'] = deepcopy(
                compound_unprioritized_color.entry_type['compound']
            )
            self.unprioritized_color = unprioritized_color
        else:
            self.recolor_unprioritized_entries(
                unprioritized_entry_uuids, recolor_unprioritized_entries
            )

            # Set the unprioritized_color attribute before returning.
            self.unprioritized_color = deepcopy(recolor_unprioritized_entries)

    # def set_color_priority1(
    #     self,
    #     new_color_priority: dict[str, dict[str, dict[tuple[str, str], float]]],
    #     recolor_unprioritized_entries: Union[str, dict[str, tuple[str, str]]] = False,
    #     color_associated_compounds: Literal['high', 'low', 'average'] = None,
    #     colormap: Colormap = None
    # ) -> None:
    #     """
    #     Set the color_priority attribute. Entry elements in the subelements attribute are
    #     automatically reordered.

    #     A single Entry (e.g., representing KOs or compounds) can occur multiple times on a map
    #     (e.g., as different rectangles or circles), and thus have multiple Graphics elements. It is
    #     required here that Graphics elements of the same type (e.g., rectangle type Graphics or line
    #     type Graphics) for an Entry must all have the same foreground and background colors if they
    #     are to be ordered.

    #     Entries with higher priority fg/bg colors are placed last in the subelements attribute and
    #     in KGML files, and they are rendered in the foreground of the map. The lowest priority
    #     entries are always those without fg/bg colors defined in the color_priority attribute; these
    #     entries are placed first in the subelements attribute and KGML files, and they are rendered
    #     in the background of the map and thus can be overlaid by higher priority entries.

    #     Parameters
    #     ==========
    #     new_color_priority : dict[str, dict[str, dict[tuple[str, str], float]]]
    #         This dictionary is the basis of the color_priority attribute.

    #         It has the same structure as the color_priority attribute. Outermost dict keys are Entry
    #         types, any of the possible values of the type attribute of the Entry class, e.g.,
    #         'ortholog' and 'compound'. Middle dict keys are Graphics types, any of the possible
    #         values of the type attribute of the Graphics class, e.g., 'rectangle' and 'line'. Inner
    #         dict keys are length-2 tuples of fgcolor and bgcolor hex codes, respectively. Inner dict
    #         values are non-negative numbers indicating the priority of Entry Graphics with the given
    #         foreground and background colors: higher numbers indicate higher priority colors.

    #         What is actually used to set color_priority is a deep copy of the argument in which
    #         fg/bg color combinations (entries in each inner dict) are reordered by priority value
    #         ascending, so that the lowest priority colors appear first in each inner dict. The order
    #         of Entry and Graphics types (outer and middle dict entries) do not not change: so if,
    #         for example, 'ortholog' appears before 'compound' in the outermost dict keys, then
    #         ortholog entries will occur before compound entries in the KGML file, and compounds can
    #         be drawn over orthologs.

    #     recolor_unprioritized_entries : Union[str, dict[str, dict[str, tuple[str, str]]]], False
    #         Recolor unprioritized entries, either automatically with a string argument, or with a
    #         custom dictionary argument for fine-tuning foreground and background colors by Entry
    #         type. The valid string arguments for automatic recoloring are 'w' and 'g'.

    #         It is assumed that global maps contain reaction lines and compound circles, so automatic
    #         recoloring is tailored to ortholog Entry line Graphics and compound Entry circle
    #         Graphics. 'w' erases unprioritized lines and circles by coloring them entirely white.
    #         'g' colors them a light gray (#E0E0E0), consistent with other "unidentified" reactions
    #         in the base map.

    #         It is assumed that overview maps contain reaction lines (drawn as arrows) and compound
    #         circles. Unlike global maps, 'w' colors unprioritized arrows black, consistent with
    #         other "unidentified" reactions in the base map. 'w' colors unprioritized circles white
    #         (with a black border). 'g' colors arrows and circles light gray.

    #         It is assumed that standard maps contain ortholog boxes or lines and compound circles.
    #         'w' colors unprioritized boxes white, with black text; lines black; and circles white.
    #         'g' colors unprioritized boxes light gray, with black text; lines light gray; and
    #         circles light gray. Compound rectangles are also found in a small number of KGML files
    #         (see 00121, 00621, 01052, 01054), and it is best if the size of these is set to zero so
    #         they don't obscure the chemical structure drawings to which they correspond on maps (see
    #         `anvio.keggmapping.Mapper._zero_out_compound_rectangles`).

    #         A custom dictionary argument can be used to set fg/bg colors in detail. Outer dict keys
    #         are Entry types, e.g., 'ortholog', 'compound'. Inner dict keys are Graphics types, e.g.,
    #         'rectangle', 'line'. Inner dict values are length-2 tuples of color hex codes for fg and
    #         bg colors, respectively. This is shown in the following example, which sets the fg
    #         (text) color of unprioritized ortholog Entry rectangle Graphics to dark gray and the bg
    #         to light gray; the fg color of unprioritized ortholog Entry line Graphics to black and
    #         the bg to white; and the fg (border) of unprioritized compounds to black and the bg to
    #         white.
    #         {
    #             'ortholog': {
    #                 'rectangle': ('#A9A9A9', '#E0E0E0'),
    #                 'line': ('#000000', '#FFFFFF')
    #             },
    #             'compound': {
    #                 'circle': ('#000000', '#FFFFFF')
    #             }
    #         }

    #     color_associated_compounds : Literal['high', 'low', 'average'], None
    #         Automatically set the background color of compound entries based on the color priority
    #         of ortholog entries involving the compounds. By default, compounds participating in
    #         reactions are circles, and orthologs are lines on global/overview maps and boxes or
    #         lines on standard maps.

    #         An argument of 'high' or 'low' sets the compound background color to the bg color of the
    #         ortholog with the highest or lowest priority fg/bg color combination. 'average' sets the
    #         bg color to the average bg color of orthologs with prioritized colors -- unprioritized
    #         orthologs are not taken into account. 'average' should only be used if priority values
    #         are normalized to the interval [0, 1] and can be converted to a color given by the
    #         colormap argument. The average priority value of the orthologs with prioritized colors
    #         is mapped to a bg color for compound Entry circle Graphics.

    #         Automatically colored compound entries are added to the color_priority attribute, and
    #         Entry elements in the subelements attribute are reordered accordingly. Compound entries
    #         that are already in the color_priority attribute are exempt from recoloring and given
    #         higher priority than automatically recolored compound entries.

    #     colormap : matplotlib.colors.Colormap, None
    #         If 'average' is used as the color_associated_compounds argument, a colormap must be
    #         provided to map averaged priority values on the interval [0, 1] to a background color
    #         for compound Entry circle Graphics.
    #     """
    #     # Check that new_color_priority only contains positive priority values.
    #     for new_entry_color_priority in new_color_priority.values():
    #         for new_graphics_color_priority in new_entry_color_priority.values():
    #             for priority in new_graphics_color_priority.values():
    #                 assert priority >= 0

    #     # Make the color_priority attribute dict, reordering colors from lowest to highest priority.
    #     color_priority = {}
    #     for entry_type, new_entry_color_priority in new_color_priority.items():
    #         color_priority[entry_type] = entry_type_color_priority = {}
    #         for graphics_type, new_graphics_color_priority in new_entry_color_priority.items():
    #             entry_type_color_priority[graphics_type] = graphics_type_color_priority = {}
    #             for colors, priority in sorted(
    #                 new_graphics_color_priority.items(), key=lambda item: item[1]
    #             ):
    #                 graphics_type_color_priority[colors] = priority
    #     self.color_priority = color_priority

    #     # Reorder Entry elements in the subelements attribute from lowest to highest priority.
    #     unprioritized_entry_uuids = self.order_entries_by_color_priority()

    #     if recolor_unprioritized_entries:
    #         if isinstance(recolor_unprioritized_entries, str):
    #             assert recolor_unprioritized_entries in ('w', 'g')

    #             # Recolor orthologs.
    #             if recolor_unprioritized_entries == 'w':
    #                 if self.is_overview_map:
    #                     color_hex_code = '#000000'
    #                 else:
    #                     color_hex_code = '#FFFFFF'
    #             elif recolor_unprioritized_entries == 'g':
    #                 color_hex_code = '#E0E0E0'
    #             self.recolor_unprioritized_ortholog_entries(
    #                 unprioritized_entry_uuids, color_hex_code
    #             )

    #             # Recolor compounds.
    #             if recolor_unprioritized_entries == 'w':
    #                 color_hex_code = '#FFFFFF'
    #             elif recolor_unprioritized_entries == 'g':
    #                 color_hex_code = '#E0E0E0'
    #             self.recolor_unprioritized_compound_entries(
    #                 unprioritized_entry_uuids, color_hex_code
    #             )
    #         else:
    #             self.recolor_unprioritized_entries(
    #                 unprioritized_entry_uuids, recolor_unprioritized_entries
    #             )

    #     if color_associated_compounds is None:
    #         return
    #     self.color_associated_compounds(color_associated_compounds, colormap=colormap)

    #     # Reorder compound entries in the subelements attribute according to color priority.
    #     unprioritized_entry_uuids = self.order_entries_by_color_priority()

    #     # Recolor compounds that have not been assigned a color priority.
    #     if recolor_unprioritized_entries:
    #         if isinstance(recolor_unprioritized_entries, str):
    #             if recolor_unprioritized_entries == 'w':
    #                 color_hex_code = '#FFFFFF'
    #             elif recolor_unprioritized_entries == 'g':
    #                 color_hex_code = '#E0E0E0'
    #             self.recolor_unprioritized_compound_entries(
    #                 unprioritized_entry_uuids, color_hex_code
    #             )
    #         else:
    #             self.recolor_unprioritized_entries(
    #                 unprioritized_entry_uuids, recolor_unprioritized_entries
    #            )

    def order_entries_by_color_priority(self) -> list[str]:
        """
        Reorder Entry (e.g., 'ortholog', 'compound') UUIDs by color Priority in the subelements
        attribute of the Pathway. This determines how entries are ordered in KGML files and rendered
        in maps.

        Returns
        =======
        list[str]
            UUIDs of Entry elements without a color Priority.
        """
        # Entries have different types ('ortholog', 'enzyme', 'reaction', 'compound', etc.). Group
        # entries into two classes. "Qualifying" entries have types represented in the
        # color_priority attribute. Other entries do not have types represented in the
        # color_priority attribute, and are assigned the lowest nominal Priority of -1.0. No effort
        # is made to sort these other entries in any way.
        reordered_entry_uuids: list[str] = []
        unprioritized_entry_uuids: list[str] = []
        qualifying_entry_uuids: dict[str, list[str]] = {
            entry_type: [] for entry_type in self.color_priority.entry_type
        }
        for entry_uuid in self.subelements['entry']:
            entry: Entry = self.uuid_element_lookup[entry_uuid]
            if entry.type in self.color_priority.entry_type:
                qualifying_entry_uuids[entry.type].append(entry_uuid)
            else:
                reordered_entry_uuids.append(entry_uuid)
                unprioritized_entry_uuids.append(entry_uuid)

        # Sort "qualifying" entries. Loop through each Entry type in the color_priority attribute.
        for entry_type, entry_type_color_priority in self.color_priority.entry_type.items():
            # Retrieve each Entry object of the type. Its Priority is determined from Graphics
            # foreground and background colors.
            type_qualifying_entry_uuids = qualifying_entry_uuids[entry_type]
            type_priority_entry_uuids: dict[float, list[str]] = {}
            for entry_uuid in type_qualifying_entry_uuids:
                entry: Entry = self.uuid_element_lookup[entry_uuid]

                # Ensure that all of the Entry Graphics elements are of the same type and have the
                # same fg and bg colors.
                graphics_types: list[str] = []
                fgcolors: list[str] = []
                bgcolors: list[str] = []
                for graphics_uuid in entry.subelements['graphics']:
                    graphics_element: Graphics = self.uuid_element_lookup[graphics_uuid]
                    graphics_types.append(graphics_element.type)
                    fgcolors.append(graphics_element.fgcolor)
                    bgcolors.append(graphics_element.bgcolor)
                if len(set(graphics_types)) != 1:
                    graphics_type_message = ', '.join([f"'{gt}'" for gt in graphics_types])
                    raise ConfigError(
                        f"The Graphics elements for the Entry with UUID '{entry_uuid}' do not "
                        "have the same type, which is required for ordering entries based on "
                        f"color. Graphics have the following types: {graphics_type_message}"
                    )
                if len(set(fgcolors)) != 1 or len(set(bgcolors)) != 1:
                    raise ConfigError(
                        f"The Graphics elements in the Entry with UUID '{entry_uuid}' do not "
                        "have consistent foreground and background colors, which is required "
                        "for ordering entries based on color."
                    )

                graphics_type = graphics_types[0]
                try:
                    graphics_type_color_priority = entry_type_color_priority.graphics_type[graphics_type]
                except KeyError:
                    raise ConfigError(
                        f"The Graphics type, '{graphics_type}', does not have an entry in the "
                        f"color priority for the Entry type, '{entry_type}'."
                    )
                try:
                    priority = graphics_type_color_priority.color[(fgcolors[0], bgcolors[0])]
                except KeyError:
                    # The Entry does not have prioritized colors.
                    priority = -1.0

                try:
                    type_priority_entry_uuids[priority].append(entry_uuid)
                except KeyError:
                    type_priority_entry_uuids[priority] = [entry_uuid]

            # Add the reordered UUIDs of the Entry type to the new list of Entry UUIDs and to the
            # dict mapping Priority values to UUIDs of entries of all types.
            for priority, entry_uuids in sorted(type_priority_entry_uuids.items()):
                reordered_entry_uuids += entry_uuids

            try:
                unprioritized_entry_uuids += type_priority_entry_uuids[-1.0]
            except KeyError:
                pass

        self.subelements['entry'] = reordered_entry_uuids

        return unprioritized_entry_uuids

    def recolor_unprioritized_reaction_entries(
        self,
        unprioritized_entry_uuids: list[str],
        color_hex_code: str
    ) -> PathwayUnprioritizedColor:
        """
        Recolor reaction-like entries without a color Priority.

        Reaction entries are expected to have Graphics elements of type 'line' in global and
        overview maps and type 'rectangle' or 'line' in standard maps. The color_hex_code argument
        color is applied to the foreground of a line, and the background is made white. The color is
        applied to the background of a rectangle, and the foreground (text) is made black.

        Parameters
        ==========
        unprioritized_entry_uuids : list[str]
            List of UUIDs of all entries without a color Priority.

        color_hex_code : str
            Hex code of the color for reaction Graphics.

        Returns
        =======
        reaction_unprioritized_color : PathwayUnprioritizedColor
            Combinations of foreground and background colors of unprioritized reaction Entry
            Graphics drawn below prioritized entries.
        """
        reaction_unprioritized_color = PathwayUnprioritizedColor()
        reaction_like_types = Entry.reaction_like_types
        if self.is_global_map or self.is_overview_map:
            for entry_type in reaction_like_types:
                entry_unprioritized_color = EntryUnprioritizedColor()
                reaction_unprioritized_color.entry_type[entry_type] = entry_unprioritized_color
                entry_unprioritized_color.graphics_type['line'] = (color_hex_code, '#FFFFFF')
            self.recolor_unprioritized_entries(
                unprioritized_entry_uuids, reaction_unprioritized_color
            )
        else:
            for entry_type in reaction_like_types:
                entry_unprioritized_color = EntryUnprioritizedColor()
                reaction_unprioritized_color.entry_type[entry_type] = entry_unprioritized_color
                entry_unprioritized_color.graphics_type['rectangle'] = ('#000000', color_hex_code)
                entry_unprioritized_color.graphics_type['line'] = (color_hex_code, '#FFFFFF')
            self.recolor_unprioritized_entries(
                unprioritized_entry_uuids, reaction_unprioritized_color
            )

        return reaction_unprioritized_color

    # def recolor_unprioritized_ortholog_entries(
    #     self,
    #     unprioritized_entry_uuids: list[str],
    #     color_hex_code: str
    # ) -> None:
    #     """
    #     Recolor orthologs without a color priority.

    #     Ortholog entries are expected to have Graphics elements of type 'line' in global and
    #     overview maps and type 'rectangle' or 'line' in standard maps. The color is applied to the
    #     foreground of a line, and the background is made white. The color is applied to the
    #     background of a rectangle, and the foreground (text) is made black.

    #     Parameters
    #     ==========
    #     unprioritized_entry_uuids : list[str]
    #         List of UUIDs of all entries without a color priority.

    #     color_hex_code : str
    #         Hex code of the color for ortholog graphics.
    #     """
    #     if self.is_global_map or self.is_overview_map:
    #         self.recolor_unprioritized_entries(
    #             unprioritized_entry_uuids, {'ortholog': {'line': (color_hex_code, '#FFFFFF')}}
    #         )
    #     else:
    #         self.recolor_unprioritized_entries(
    #             unprioritized_entry_uuids,
    #             {'ortholog': {
    #                 'rectangle': ('#000000', color_hex_code), 'line': (color_hex_code, '#000000')
    #             }}
    #         )

    def recolor_unprioritized_compound_entries(
        self,
        unprioritized_entry_uuids: list[str],
        color_hex_code: str
    ) -> PathwayUnprioritizedColor:
        """
        Recolor entries of the 'compound' type without a color Priority.

        Compound entries are expected to have Graphics elements of type 'circle'. In global maps, the
        color is applied to both the background (fill) and foreground (border) of the circle. In
        overview and standard maps, the background is colored, and the foreground is made black.

        Compound rectangles are also found in a small number of KGML files (see 00121, 00621, 01052,
        01054), and it is best if the size of these is set to zero so they don't obscure the
        chemical structure drawings to which they correspond on maps (see
        `anvio.keggmapping.Mapper._zero_out_compound_rectangles`).

        Parameters
        ==========
        unprioritized_entry_uuids : list[str]
            List of UUIDs of all entries without a color Priority.

        color_hex_code : str
            Hex code of the color for compound Graphics.

        Returns
        =======
        compound_unprioritized_color : PathwayUnprioritizedColor
            Combinations of foreground and background colors of unprioritized compound Entry
            Graphics drawn below prioritized entries.
        """
        compound_unprioritized_color = PathwayUnprioritizedColor()
        entry_unprioritized_color = EntryUnprioritizedColor()
        compound_unprioritized_color.entry_type['compound'] = entry_unprioritized_color
        if self.is_global_map:
            entry_unprioritized_color.graphics_type['circle'] = (color_hex_code, color_hex_code)
            self.recolor_unprioritized_entries(
                unprioritized_entry_uuids, compound_unprioritized_color
            )
        else:
            entry_unprioritized_color.graphics_type['circle'] = ('#000000', color_hex_code)
            self.recolor_unprioritized_entries(
                unprioritized_entry_uuids, compound_unprioritized_color
            )

        return compound_unprioritized_color

    def recolor_unprioritized_entries(
        self,
        unprioritized_entry_uuids: list[str],
        pathway_unprioritized_color: PathwayUnprioritizedColor
    ) -> None:
        """
        Entries without a color Priority are recolored by Entry type.

        Parameters
        ==========
        unprioritized_entry_uuids : list[str]
            List of UUIDs of all entries without a color Priority.

        pathway_unprioritized_color : PathwayUnprioritizedColor
            Combinations of foreground and background colors of unprioritized Entry Graphics drawn
            below prioritized entries.
        """
        # Prevent unprioritized entries from being assigned prioritized colors.
        for entry_type, entry_unprioritized_color in pathway_unprioritized_color.entry_type.items():
            try:
                entry_type_color_priority = self.color_priority.entry_type[entry_type]
            except KeyError:
                # The Entry type has no prioritized colors.
                continue

            for (
                graphics_type, graphics_unprioritized_color
            ) in entry_unprioritized_color.graphics_type.items():
                try:
                    graphics_type_color_priority = entry_type_color_priority.graphics_type[
                        graphics_type
                    ]
                except KeyError:
                    # The Entry Graphics type has no prioritized colors.
                    continue

                if graphics_unprioritized_color in graphics_type_color_priority.color:
                    raise ConfigError(
                        "Unprioritized Entry Graphics cannot be assigned the same combination of "
                        "foreground and background colors as prioritized entries of the same Entry "
                        f"and Graphics types. The Entry type is '{entry_type}'. The Graphics type "
                        f"is '{graphics_type}'. The unprioritized foreground color is "
                        f"'{graphics_unprioritized_color[0]}' and background color is "
                        f"'{graphics_unprioritized_color[1]}'."
                    )

        # Set colors of unprioritized Entry Graphics.
        for entry_uuid in unprioritized_entry_uuids:
            entry: Entry = self.uuid_element_lookup[entry_uuid]
            try:
                entry_unprioritized_color = pathway_unprioritized_color.entry_type[entry.type]
            except KeyError:
                # The type of the Entry object is not under consideration for recoloring.
                continue

            for graphics_uuid in entry.subelements['graphics']:
                graphics: Graphics = self.uuid_element_lookup[graphics_uuid]
                try:
                    fgcolor_hex_code, bgcolor_hex_code = entry_unprioritized_color.graphics_type[
                        graphics.type
                    ]
                except KeyError:
                    # The type of the Graphics object is not under consideration for recoloring.
                    continue
                graphics.fgcolor = fgcolor_hex_code
                graphics.bgcolor = bgcolor_hex_code

    def color_associated_compounds(
        self,
        transfer: Literal['high', 'low', 'average'],
        colormap: Colormap = None
    ) -> None:
        """
        Set the color of compound entries based on the color Priority of reaction-like ('ortholog',
        'enzyme', 'reaction', 'gene', 'group' type) entries involving the compounds.

        Compound entries are assumed to have Graphics elements of type 'circle'. If the map is
        global, color both the background (interior) and foreground (border) of the circle. In
        overview and standard maps, the background is colored, and the foreground is made black.

        Compound rectangles are also found in a small number of KGML files (see 00121, 00621, 01052,
        01054), and it is best if the size of these is set to zero so they don't obscure the
        chemical structure drawings to which they correspond on maps (see
        `anvio.keggmapping.Mapper._zero_out_compound_rectangles`).

        Parameters
        ==========
        transfer : Literal['high', 'low', 'average']
            An argument of 'high' or 'low' sets the compound color to the bg color of the reaction
            Entry with the highest or lowest Priority fg/bg color combination.

            'average' sets the compound color to the average bg color of reactions with prioritized
            colors: unprioritized reactions are not taken into account. 'average' should only be
            used if Priority values are normalized to the interval [0, 1] and can be converted to a
            color given by the colormap argument. The average Priority value of the reactions with
            prioritized colors is mapped to a color for the compound Entry.

            Compound entries that are already in the color_priority attribute are exempt from
            recoloring and given higher priority than recolored compound entries.

        colormap : matplotlib.colors.Colormap, None
            If 'average' is used as the transfer argument value, a colormap must be provided to map
            averaged Priority values on the interval [0, 1] to a background color for compound
            entries.
        """
        # Make Reaction elements searchable by name. Reaction elements link compound entries to
        # reaction entries. (Note that Reaction elements are not the same as Entry elements with
        # types that we call reaction types.)
        name_reaction: dict[str, Reaction] = {}
        for entry_uuid in self.subelements['reaction']:
            reaction: Reaction = self.uuid_element_lookup[entry_uuid]
            name_reaction[reaction.name] = reaction

        # For each compound Entry with associated color-prioritized reaction entries, record the
        # colors and priorities of these entries.
        compound_uuid_color_priorities: dict[str, list[tuple[str, float]]] = {}
        # Loop through each reaction Entry.
        reaction_like_types = Entry.reaction_like_types
        for entry_uuid in self.subelements['entry']:
            entry: Entry = self.uuid_element_lookup[entry_uuid]
            entry_type = entry.type
            if entry_type not in reaction_like_types:
                continue

            # Ensure that all of the reaction Entry Graphics elements have the same combination of
            # fg/bg colors.
            graphics_types: list[str] = []
            fgcolors: list[str] = []
            bgcolors: list[str] = []
            for graphics_uuid in entry.subelements['graphics']:
                graphics: Graphics = self.uuid_element_lookup[graphics_uuid]
                graphics_types.append(graphics.type)
                fgcolors.append(graphics.fgcolor)
                bgcolors.append(graphics.bgcolor)
            if len(set(graphics_types)) != 1:
                graphics_type_message = ', '.join([f"'{gt}'" for gt in graphics_types])
                raise ConfigError(
                    f"The Graphics elements for the {entry_type} Entry with UUID '{entry_uuid}' do "
                    "not have the same type, which is required for ordering entries based on "
                    f"color. Graphics have types: {graphics_type_message}"
                )
            if len(set(fgcolors)) != 1 or len(set(bgcolors)) != 1:
                raise ConfigError(
                    f"The Graphics elements for the {entry_type} Entry with UUID '{entry_uuid}' do "
                    "not have consistent foreground and background colors, which is required for "
                    f"ordering entries based on color."
                )

            graphics_type = graphics_types[0]
            fgcolor = fgcolors[0]
            bgcolor = bgcolors[0]

            try:
                entry_type_color_priority = self.color_priority.entry_type[entry_type]
            except KeyError:
                # No color priority has been defined for the type of Entry.
                continue

            try:
                graphics_type_color_priority = entry_type_color_priority.graphics_type[
                    graphics_type
                ]
            except KeyError:
                # No color priority has been defined for the type of Graphics.
                continue

            try:
                priority = graphics_type_color_priority.color[(fgcolor, bgcolor)]
            except KeyError:
                # Unprioritized reaction entries do not affect the color of associated compounds.
                continue

            reaction_name = entry.reaction
            if reaction_name is None:
                # The Entry does not have a reaction attribute.
                continue

            try:
                reaction = name_reaction[reaction_name]
            except KeyError:
                # No Reaction element has the reaction name attached to the Entry. This can be
                # explained by an error in KGML file contents.
                continue

            if graphics.type == 'line':
                reaction_entry_color = fgcolor
            else:
                reaction_entry_color = bgcolor

            # Record the color and priority of the reaction for substrate compound entries.
            for substrate_uuid in reaction.subelements['substrate']:
                substrate: Substrate = self.uuid_element_lookup[substrate_uuid]
                split_substrate_names = [
                    split_name.split(':') for split_name in substrate.name.split()
                ]
                for split_name in split_substrate_names:
                    for compound_entry in self.kegg_id_element_lookup[split_name[1]]:
                        if not isinstance(compound_entry, Entry):
                            continue
                        compound_entry: Entry
                        try:
                            compound_uuid_color_priorities[compound_entry.uuid].append(
                                (reaction_entry_color, priority)
                            )
                        except KeyError:
                            compound_uuid_color_priorities[compound_entry.uuid] = [
                                (reaction_entry_color, priority)
                            ]

            # Record the color and priority of the reaction for product compound entries.
            for product_uuid in reaction.subelements['product']:
                product: Product = self.uuid_element_lookup[product_uuid]
                split_product_names = [split_name.split(':') for split_name in product.name.split()]
                for split_name in split_product_names:
                    for compound_entry in self.kegg_id_element_lookup[split_name[1]]:
                        if not isinstance(compound_entry, Entry):
                            continue
                        compound_entry: Entry
                        try:
                            compound_uuid_color_priorities[compound_entry.uuid].append(
                                (reaction_entry_color, priority)
                            )
                        except KeyError:
                            compound_uuid_color_priorities[compound_entry._uuid] = [
                                (reaction_entry_color, priority)
                            ]

        # Make compound entries searchable by ID (not the UUID assigned by anvi'o), which should be
        # a unique element ID in the Pathway.
        id_compound_entry: dict[str, Entry] = {}
        for entry_uuid in self.subelements['entry']:
            entry: Entry = self.uuid_element_lookup[entry_uuid]
            if entry.type != 'compound':
                continue
            id_compound_entry[entry.id] = entry

        # Define functions for finding compound Entry color.
        def _get_high_color(color_priorities: list[tuple[str, float]]) -> tuple[str, float]:
            return sorted(color_priorities, key=lambda t: -t[1])[0]

        def _get_low_color(color_priorities: list[tuple[str, float]]) -> tuple[str, float]:
            return sorted(color_priorities, key=lambda t: t[1])[0]

        def _get_average_color(color_priorities: list[tuple[str, float]]) -> tuple[str, float]:
            priority = np.mean([t[1] for t in color_priorities])
            color = rgb2hex(colormap(priority))
            return color, priority

        if transfer == 'high':
            get_color_priority = _get_high_color
        elif transfer == 'low':
            get_color_priority = _get_low_color
        elif transfer == 'average':
            get_color_priority = _get_average_color
        else:
            raise ConfigError(
                "The 'transfer' argument value must be one of the strings, 'high', 'low', or "
                "'average'."
            )

        # Set compound Entry color.
        for compound_uuid, color_priorities in compound_uuid_color_priorities.items():
            compound_entry: Entry = self.uuid_element_lookup[compound_uuid]

            # Get all of the Graphics elements for the Entry.
            graphics_elements: list[Graphics] = []
            for graphics_uuid in compound_entry.subelements['graphics']:
                graphics_elements.append(self.uuid_element_lookup[graphics_uuid])

            set_color = True
            for graphics in graphics_elements:
                try:
                    compound_color_priority = self.color_priority.entry_type['compound']
                except KeyError:
                    # No color priority has been defined for compound entries.
                    continue

                try:
                    circle_color_priority = compound_color_priority.graphics_type['circle']
                except KeyError:
                    # No color priority has been defined for circle Graphics in compound entries.
                    continue

                try:
                    # The compound Entry has already been assigned a color Priority, so don't
                    # recolor it automatically.
                    circle_color_priority.color[(graphics.fgcolor, graphics.bgcolor)]
                except KeyError:
                    continue
                set_color = False
            if not set_color:
                continue

            compound_color, compound_priority = get_color_priority(color_priorities)
            # Set the color of each Graphics element.
            for graphics in graphics_elements:
                graphics.bgcolor = compound_color
                if self.is_global_map:
                    graphics.fgcolor = compound_color

            # Record the compound Element color priority.
            try:
                compound_color_priority = self.color_priority.entry_type['compound']
            except KeyError:
                compound_color_priority = EntryColorPriority()
                self.color_priority.entry_type['compound'] = compound_color_priority

            try:
                circle_color_priority = compound_color_priority.graphics_type['circle']
            except KeyError:
                circle_color_priority = GraphicsColorPriority()
                compound_color_priority.graphics_type['circle'] = circle_color_priority

            if self.is_global_map:
                circle_color_priority.color[(compound_color, compound_color)] = compound_priority
            else:
                circle_color_priority[(graphics.fgcolor, compound_color)] = compound_priority

    # def color_associated_compounds1(
    #     self,
    #     transfer: Literal['high', 'low', 'average'],
    #     colormap: Colormap = None
    # ) -> None:
    #     """
    #     Set the color of compound entries based on the color priority of ortholog entries involving
    #     the compounds.

    #     Compound entries are assumed to have Graphics elements of type 'circle'. If the map is
    #     global, color both the background (interior) and foreground (border) of the circle. In
    #     overview and standard maps, the background is colored, and the foreground is made black.
    #     Compound rectangles are also found in a small number of KGML files (see 00121, 00621, 01052,
    #     01054), and it is best if the size of these is set to zero so they don't obscure the
    #     chemical structure drawings to which they correspond on maps (see
    #     `anvio.keggmapping.Mapper._zero_out_compound_rectangles`).

    #     Parameters
    #     ==========
    #     transfer : Literal['high', 'low', 'average']
    #         An argument of 'high' or 'low' sets the compound color to the bg color of the ortholog
    #         with the highest or lowest priority fg/bg color combination. 'average' sets the compound
    #         color to the average bg color of orthologs with prioritized colors -- unprioritized
    #         orthologs are not taken into account. 'average' should only be used if priority values
    #         are normalized to the interval [0, 1] and can be converted to a color given by the
    #         colormap argument. The average priority value of the orthologs with prioritized colors
    #         is mapped to a color for the compound Entry.

    #         Compound entries that are already in the color_priority attribute are exempt from
    #         recoloring and given higher priority than recolored compound entries.

    #     colormap : matplotlib.colors.Colormap, None
    #         If 'average' is used as an argument to transfer, a colormap must be provided to map
    #         averaged priority values on the interval [0, 1] to a background color for compound
    #         entries.
    #     """
    #     # Make Reaction elements searchable by name (KEGG IDs). Reaction elements link Compound
    #     # elements to ortholog Entry elements.
    #     name_reaction: dict[str, Reaction] = {}
    #     for entry_uuid in self.subelements['reaction']:
    #         reaction: Reaction = self.uuid_element_lookup[entry_uuid]
    #         name_reaction[reaction.name] = reaction

    #     # For each compound Entry with associated color-prioritized ortholog entries, record the
    #     # colors and priorities of these entries.
    #     compound_uuid_color_priorities: dict[str, list[tuple[str, float]]] = {}
    #     # Loop through each ortholog Entry.
    #     for entry_uuid in self.subelements['entry']:
    #         entry: Entry = self.uuid_element_lookup[entry_uuid]
    #         if entry.type != 'ortholog':
    #             continue

    #         # Ensure that all of the ortholog Entry Graphics elements have the same fg/bg colors.
    #         graphics_types: list[str] = []
    #         fgcolors: list[str] = []
    #         bgcolors: list[str] = []
    #         for graphics_uuid in entry.subelements['graphics']:
    #             graphics: Graphics = self.uuid_element_lookup[graphics_uuid]
    #             graphics_types.append(graphics.type)
    #             fgcolors.append(graphics.fgcolor)
    #             bgcolors.append(graphics.bgcolor)
    #         if len(set(graphics_types)) != 1:
    #             graphics_type_message = ', '.join([f"'{gt}'" for gt in graphics_types])
    #             raise AssertionError(
    #                 f"The Graphics elements for the Entry with UUID '{entry_uuid}' do not "
    #                 "have the same type, which is required for ordering entries based on "
    #                 f"color. Graphics have types: {graphics_type_message}"
    #             )
    #         if len(set(fgcolors)) != 1 or len(set(bgcolors)) != 1:
    #             raise AssertionError(
    #                 "The Graphics elements in the ortholog Entry with the following UUID do not "
    #                 "have consistent foreground and background colors, which is required for "
    #                 f"ordering entries based on color: {entry_uuid}"
    #             )

    #         graphics_type = graphics_types[0]
    #         fgcolor = fgcolors[0]
    #         bgcolor = bgcolors[0]
    #         try:
    #             priority = self.color_priority['ortholog'][graphics_type][(fgcolor, bgcolor)]
    #         except KeyError:
    #             # Unprioritized ortholog entries do not affect the color of associated compounds.
    #             continue

    #         reaction_name = entry.reaction
    #         if reaction_name is None:
    #             # The ortholog is not associated with a reaction.
    #             continue

    #         try:
    #             reaction = name_reaction[reaction_name]
    #         except KeyError:
    #             # No Reaction element is present with the name of the ortholog reaction.
    #             continue

    #         if graphics.type == 'line':
    #             ortholog_color = fgcolor
    #         else:
    #             ortholog_color = bgcolor

    #         for substrate_uuid in reaction.subelements['substrate']:
    #             substrate: Substrate = self.uuid_element_lookup[substrate_uuid]
    #             split_substrate_names = [
    #                 split_name.split(':') for split_name in substrate.name.split()
    #             ]
    #             for split_name in split_substrate_names:
    #                 for compound_entry in self.kegg_id_element_lookup[split_name[1]]:
    #                     if not isinstance(compound_entry, Entry):
    #                         continue
    #                     compound_entry: Entry
    #                     compound_uuid = compound_entry.uuid
    #                     try:
    #                         compound_uuid_color_priorities[compound_uuid].append(
    #                             (ortholog_color, priority)
    #                         )
    #                     except KeyError:
    #                         compound_uuid_color_priorities[compound_uuid] = [
    #                             (ortholog_color, priority)
    #                         ]

    #         for product_uuid in reaction.subelements['product']:
    #             product: Product = self.uuid_element_lookup[product_uuid]
    #             split_product_names = [split_name.split(':') for split_name in product.name.split()]
    #             for split_name in split_product_names:
    #                 for compound_entry in self.kegg_id_element_lookup[split_name[1]]:
    #                     if not isinstance(compound_entry, Entry):
    #                         continue
    #                     compound_entry: Entry
    #                     compound_uuid = compound_entry.uuid
    #                     try:
    #                         compound_uuid_color_priorities[compound_uuid].append(
    #                             (ortholog_color, priority)
    #                         )
    #                     except KeyError:
    #                         compound_uuid_color_priorities[compound_uuid] = [
    #                             (ortholog_color, priority)
    #                         ]

    #     # Make compound entries searchable by ID, which should be a unique pathway element ID.
    #     id_compound_entry: dict[str, Entry] = {}
    #     for entry_uuid in self.subelements['entry']:
    #         entry: Entry = self.uuid_element_lookup[entry_uuid]
    #         if entry.type != 'compound':
    #             continue
    #         id_compound_entry[entry.id] = entry

    #     # Define functions for finding compound Entry color.
    #     def _get_high_color(color_priorities: list[tuple[str, float]]) -> tuple[str, float]:
    #         return sorted(color_priorities, key=lambda t: -t[1])[0]

    #     def _get_low_color(color_priorities: list[tuple[str, float]]) -> tuple[str, float]:
    #         return sorted(color_priorities, key=lambda t: t[1])[0]

    #     def _get_average_color(color_priorities: list[tuple[str, float]]) -> tuple[str, float]:
    #         priority = np.mean([t[1] for t in color_priorities])
    #         color = rgb2hex(colormap(priority))
    #         return color, priority

    #     if transfer == 'high':
    #         get_color_priority = _get_high_color
    #     elif transfer == 'low':
    #         get_color_priority = _get_low_color
    #     elif transfer == 'average':
    #         get_color_priority = _get_average_color
    #     else:
    #         raise AssertionError

    #     # Set compound Entry color.
    #     for compound_uuid, color_priorities in compound_uuid_color_priorities.items():
    #         compound: Union[Substrate, Product] = self.uuid_element_lookup[compound_uuid]
    #         compound_entry: Entry = id_compound_entry[compound.id]

    #         # Get all of the Graphics elements for the Entry.
    #         graphics_elements: list[Graphics] = []
    #         for graphics_uuid in compound_entry.subelements['graphics']:
    #             graphics_elements.append(self.uuid_element_lookup[graphics_uuid])

    #         set_color = True
    #         for graphics in graphics_elements:
    #             try:
    #                 # The compound Entry has already been assigned a color priority, so don't
    #                 # recolor it automatically.
    #                 self.color_priority['compound']['circle'][(graphics.fgcolor, graphics.bgcolor)]
    #                 set_color = False
    #             except KeyError:
    #                 continue
    #         if not set_color:
    #             continue

    #         compound_color, compound_priority = get_color_priority(color_priorities)
    #         # Set the color of each Graphics element.
    #         for graphics in graphics_elements:
    #             graphics.bgcolor = compound_color
    #             if self.is_global_map:
    #                 graphics.fgcolor = compound_color

    #         # Record the compound Element color priority.
    #         try:
    #             entry_type_color_priority = self.color_priority['compound']
    #         except KeyError:
    #             self.color_priority['compound'] = entry_type_color_priority = {}
    #         try:
    #             graphics_type_color_priority = entry_type_color_priority['circle']
    #         except KeyError:
    #             entry_type_color_priority['circle'] = graphics_type_color_priority = {}
    #         if self.is_global_map:
    #             graphics_type_color_priority[(compound_color, compound_color)] = compound_priority
    #         else:
    #             graphics_type_color_priority[(graphics.fgcolor, compound_color)] = compound_priority

    def set_thickness_priority(
        self,
        new_thickness_priority: PathwayThicknessPriority,
        reset_unprioritized_thickness: Union[bool, float] = False
    ) -> None:
        """
        Set the thickness_priority attribute. Entry elements in the subelements attribute are
        automatically reordered.

        Thickness priorities apply to line Graphics of reaction-like entries (KO, EC, RN, gene,
        group) in global and overview pathways; thickness priorities cannot be set for standard
        pathways.

        Attributes
        ==========
        new_thickness_priority: PathwayThicknessPriority
            Here is what is actually set as thickness_priority. A deep copy is made of the
            PathwayThicknessPriority object. The width dictionary attribute of the object is
            reordered by Priority value ascending, so that the lowest Priority widths occur first.

        reset_unprioritized_thickness : Union[bool, float], False
            If not False, set unprioritized reaction-like Entry line Graphics to a uniform width.
            With an argument value of True, a width of 6.0 is used for global maps and 1.0 for
            overview maps. Alternatively, a custom width can be specified with a float value.
        """
        # Check that new_thickness_priority only contains nonnegative priority values.
        for priority in new_thickness_priority.width.values():
            assert priority >= 0

        # Make the object assigned to the thickness_priority attribute, reordering widths from
        # lowest to highest priority.
        thickness_priority = PathwayThicknessPriority()
        for width, priority in sorted(
            new_thickness_priority.width.items(), key=lambda item: item[1]
        ):
            thickness_priority.width[width] = priority
        self.thickness_priority = thickness_priority

        # Reorder Entry elements in the subelements attribute from lowest to highest priority.
        unprioritized_line_entry_uuids = self.order_entries_by_thickness_priority()

        if not reset_unprioritized_thickness:
            return

        # Standardize the widths of reaction-like entries with line Graphics.
        self.unprioritized_width = self.reset_thickness(
            unprioritized_line_entry_uuids,
            width=reset_unprioritized_thickness if isinstance(
                reset_unprioritized_thickness, float
            ) else None
        )

    def order_entries_by_thickness_priority(self) -> list[str]:
        """
        Reorder Entry UUIDs by thickness Priority in the subelements attribute of the Pathway. This
        determines how entries are ordered in KGML files and rendered in maps. All reaction-like
        entries with line Graphics are inspected for a prioritized width.

        Returns
        =======
        list[str]
            UUIDs of line Entry elements without a thickness Priority.
        """
        # Group entries into two classes. Prioritized entries are reaction-like entries with line
        # Graphics and widths with Priority values in the thickness_priority attribute. All other
        # entries are called unprioritized entries. Other reaction-like entries with line Graphics
        # are assigned the lowest nominal Priority of -1.0. No effort is made to sort these other
        # entries in any way.
        unprioritized_entry_uuids: list[str] = []
        unprioritized_line_entry_widths: list[float] = []
        prioritized_entry_uuids: dict[str, float] = {}
        reaction_like_types = Entry.reaction_like_types
        for entry_uuid in self.subelements['entry']:
            entry: Entry = self.uuid_element_lookup[entry_uuid]
            if entry.type not in reaction_like_types:
                unprioritized_entry_uuids.append(entry_uuid)
                break

            widths: list[float] = []
            for graphics_uuid in entry.subelements['graphics']:
                graphics: Graphics = self.uuid_element_lookup[graphics_uuid]
                if graphics.type == 'line':
                    widths.append(graphics.width)

            if not widths:
                # The Entry does not have any line Graphics.
                unprioritized_entry_uuids.append(entry_uuid)
                continue

            if len(set(widths)) > 1:
                raise ConfigError(
                    f"The Entry with UUID '{entry_uuid}' has line Graphics with different widths. "
                    "Line widths for reaction-like entries should be the same."
                )

            width = widths[0]
            try:
                priority = self.thickness_priority.width[width]
            except KeyError:
                # The line width does not correspond to a Priority value.
                unprioritized_entry_uuids.append(entry_uuid)
                unprioritized_line_entry_widths.append(width)
                continue

            prioritized_entry_uuids[entry_uuid] = priority

        for width in set(unprioritized_line_entry_widths):
            self.thickness_priority.width[width] = -1.0

        reordered_entry_uuids: list[str] = unprioritized_entry_uuids
        reordered_entry_uuids += [
            item[0] for item in sorted(prioritized_entry_uuids.items(), key=lambda item: item[1])
        ]
        self.subelements['entry'] = reordered_entry_uuids

        return unprioritized_entry_uuids

    def reset_thickness(self, unprioritized_entry_uuids: list[str], width: float = None) -> float:
        """
        In global and overview maps, set line Graphics of provided unprioritized entries to a
        uniform width.

        Parameters
        ==========
        unprioritized_entry_uuids : list[str]
            UUIDs of Entry elements without a thickness Priority.

        width : float, None
            Width value set for unprioritized Entry line Graphics. With an argument value of None, a
            width of 6.0 is used for global maps and 1.0 for overview maps.

        Returns
        =======
        float
            Width value used for unprioritized Entry line Graphics.
        """
        if width is None:
            if self.is_global_map:
                width = 6.0
            elif self.is_overview_map:
                width = 1.0
            else:
                raise ConfigError(
                    "The 'reset_thickness' method should not be called with a standard Pathway, "
                    "only with a global or overview Pathway."
                )

        for entry_uuid in unprioritized_entry_uuids:
            entry: Entry = self.uuid_element_lookup[entry_uuid]

            for graphics_uuid in entry.subelements['graphics']:
                graphics: Graphics = self.uuid_element_lookup[graphics_uuid]
                if graphics.type != 'line':
                    continue

                graphics.width = width

        return width

    def get_entries(
        self,
        entry_type: str = None,
        kegg_ids: Iterable[str] = None,
        expect_kegg_ids: bool = False
    ) -> list[Entry]:
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
        list[Entry]
            A list of Entry element objects contained in the pathway.
        """
        if entry_type is not None:
            assert entry_type in Entry.types
        if kegg_ids is not None:
            assert entry_type is None

        entries: list[Entry] = []

        if kegg_ids is None:
            for uuid in self.subelements['entry']:
                entry: Entry = self.uuid_element_lookup[uuid]
                if entry_type is not None and entry.type != entry_type:
                    continue
                entries.append(entry)
            return entries

        missing_kegg_ids: list[str] = []
        for kegg_id in kegg_ids:
            try:
                elements = self.kegg_id_element_lookup[kegg_id]
            except KeyError:
                missing_kegg_ids.append(kegg_id)
                continue
            for element in elements:
                if isinstance(element, Entry):
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
            for graphics_uuid in entry.subelements['graphics']:
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
    types : tuple[str]
        Possible entry types.

    subelement_tags : tuple[str]
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

    subelements : dict[str, list[str]]
        Keys are subelement tags, values are lists of subelement UUIDs.

    reaction_like_types : tuple[str]
        Reaction-like entry types designated by anvi'o.
    """
    tag = 'entry'
    attribute_required = {
        'id': True,
        'name': True,
        'type': True,
        'reaction': False,
        'link': False
    }
    types: tuple[str] = (
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
    subelement_tags: tuple[str] = (
        'graphics',
        'component'
    )
    reaction_like_types: tuple[str] = (
        'ortholog',
        'enzyme',
        'reaction',
        'gene',
        'group'
    )

    def __init__(self) -> None:
        self.id: str = None
        self.name: str = None
        self.type: str = None
        self.reaction: str = None
        self.link: str = None

        self.subelements: dict[str, list[str]] = {n: [] for n in self.subelement_tags}

        super().__init__()

class Graphics(Element):
    """
    A graphics element contains drawing information on the entry parent element.

    Attributes
    ==========
    types : tuple[str]
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

    coords : tuple[float], None
        Polyline coordinates of "line"-type graphical object on map.

    width : float, None
        Width of graphical object on map.

    height : float, None
        Height of graphical object on map.
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
    types: tuple[str] = (
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
        self.coords: tuple[float] = None
        self.width: float = None
        self.height: float = None

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
    types : tuple[str]
        Possible types of relations.

    subelement_tags : tuple[str]
        Possible subelement tags.

    entry1 : str, None
        ID unique to map representing a node in the relationship.

    entry2 : str, None
        ID unique to map representing the other node in the relationship.

    subelements : dict[str, list[str]]
        Keys are subelement tags, values are lists of subelement UUIDs.
    """

    tag = 'relation'
    attribute_required = {
        'entry1': True,
        'entry2': True,
        'type': True
    }
    types: tuple[str] = (
        'ECrel',
        'PPrel',
        'GErel',
        'PCrel',
        'maplink'
    )
    subelement_tags: tuple[str] = (
        'subtype',
    )
    def __init__(self) -> None:
        self.entry1: str = None
        self.entry2: str = None
        self.type: str = None

        self.subelements: dict[str, list[str]] = {n: [] for n in self.subelement_tags}

        super().__init__()

class Subtype(Element):
    """
    A subtype element specifies more detailed information about the relation.

    Attributes
    ==========
    names : tuple[str]
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
    names: tuple[str] = (
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
    types : tuple[str]
        Possible types of reactions.

    subelement_tags : tuple[str]
        Possible subelement tags.

    id : str, None
        ID of reaction unique to map.

    name : str, None
        KEGG ID(s) represented by the reaction.

    type : str, None
        Reversible vs. irreversible reaction, as drawn on the map.

    subelements : dict[str, list[str]]
        Keys are subelement tags, values are lists of subelement UUIDs.
    """
    tag = 'reaction'
    attribute_required = {
        'id': True,
        'name': True,
        'type': True
    }
    types: tuple[str] = (
        'reversible',
        'irreversible'
    )
    subelement_tags: tuple[str] = (
        'substrate',
        'product'
    )

    def __init__(self) -> None:
        self.id: str = None
        self.name: str = None
        self.type: str = None

        self.subelements: dict[str, list[str]] = {n: [] for n in self.subelement_tags}

        super().__init__()

class Substrate(Element):
    """
    A substrate element represents a substrate node in a parent reaction element.

    Attributes
    ==========
    subelement_tags : tuple[str]
        Possible subelement tags.

    id : str, None
        ID of substrate unique to map corresponding to a compound entry.

    name : str, None
        KEGG ID of the compound.

    subelements : dict[str, list[str]]
        Keys are subelement tags, values are lists of subelement UUIDs.
    """
    tag = 'substrate'
    attribute_required = {
        'id': True,
        'name': True
    }
    subelement_tags: tuple[str] = (
        'alt',
    )

    def __init__(self) -> None:
        self.id: str = None
        self.name: str = None

        self.subelements: dict[str, list[str]] = {n: [] for n in self.subelement_tags}

        super().__init__()

class Product(Element):
    """
    A product element represents a product node in a parent reaction element.

    Attributes
    ==========
    subelement_tags : tuple[str]
        Possible subelement tags.

    id : str, None
        ID of product unique to map corresponding to a compound entry.

    name : str, None
        KEGG ID of the compound.

    subelements : dict[str, list[str]]
        Keys are subelement tags, values are lists of subelement UUIDs.
    """
    tag = 'product'
    attribute_required = {
        'id': True,
        'name': True
    }
    subelement_tags: tuple[str] = (
        'alt',
    )

    def __init__(self) -> None:
        self.id: str = None
        self.name: str = None

        self.subelements: dict[str, list[str]] = {n: [] for n in self.subelement_tags}

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

    attribute_indentations : dict[tuple[str, str, str], int]
        Class variable setting the absolute indentation of element attributes placed on new lines in
        an output KGML XML file. Keys are tuples of element tag, name of the attribute before the
        line break, and name of the attribute after the line break; values are the number of spaces.
        The attributes placed on new lines and numbers of spaces are those used in KGML reference
        files.
    """
    subelement_indentation_increment: int = 4
    attribute_indentations: dict[tuple[str, str, str], int] = {
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
        try:
            assert os.path.exists(kgml_filepath)
        except:
            print(kgml_filepath)
            raise Exception

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
            kgml_element.subelements[kgml_subelement.tag].append(kgml_subelement.uuid)

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

        if not hasattr(kgml_element, 'subelements'):
            return xml_element

        # Recursively add XML subelements.
        for subelement_uuids in kgml_element.subelements.values():
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
        subelements = list(xml_element)

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
                    value: tuple[float]
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

        if subelements:
            indented_output += '>\n'
            for subelement in subelements:
                indented_output += self.get_indented_str(
                    subelement,
                    tag_indentation=tag_indentation + self.subelement_indentation_increment
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

class Drawer:
    """
    Write pathway map image files incorporating KGML data.

    Attributes
    ==========
    kegg_context : anvio.kegg.KeggContext
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
        self.kegg_context = kegg.KeggContext(args)

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
                By default, Drawer sets to False to reduce clutter. Setting to True displays KEGG
                COMPOUND IDs.

            label_orthologs : bool
                By default, Drawer sets to False for global and overview maps to reduce clutter next
                to reaction arrows and to True for standard maps. Setting to True displays KO IDs.

            label_reaction_entries : bool
                By default, Drawer sets to False to reduce clutter. Setting to True displays KEGG
                REACTION IDs.

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

        bio_pathway = KGML_parser.read(StringIO(self.xml_ops.get_str(pathway)))

        if map_filepath is None:
            if pathway.org == 'ko':
                map_filepath = os.path.join(
                    self.kegg_context.png_1x_ko_dir, f'ko{pathway.number}.png'
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
            kegg.download_org_pathway_image_files(f'{pathway.org}{pathway.number}', self.kegg_dir)

        bio_pathway.image = map_filepath

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

        bio_pathway = KGML_parser.read(StringIO(self.xml_ops.get_str(pathway)))

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
            kegg.download_org_pathway_image_files(f'{pathway.org}{pathway.number}', self.kegg_dir)

        bio_pathway.image = map_filepath

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
        bio_pathway = KGML_parser.read(StringIO(self.xml_ops.get_str(pathway)))

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
            kegg.download_org_pathway_image_files(f'{pathway.org}{pathway.number}', self.kegg_dir)

        bio_pathway.image = map_filepath

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
        kegg_context = kegg.KeggContext(args)
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
        filepaths: list[str] = []
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
