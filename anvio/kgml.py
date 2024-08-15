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
from argparse import Namespace
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas
from matplotlib.colors import Colormap, rgb2hex

from typing import Dict, Iterable, List, Literal, Tuple, Union

import anvio.kegg as kegg
import anvio.terminal as terminal

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
    
    color_priority : Dict[str, Dict[Tuple[str, str], float]]
        This defines the order of entries by the foreground and background colors of their graphics.
        Set this attribute with the method, set_color_priority.
        
        Outer dictionary keys are Entry types, corresponding to the possible values of the type
        attribute of the Entry class, e.g., 'ortholog' and 'compound'. Keys of inner dictionaries
        are length-2 tuples of fgcolor and bgcolor hex codes, respectively. Inner dictionary values
        are non-negative numbers indicating the priority of entries with the given foreground and
        background colors. Higher numbers indicate higher priority colors.
        
        The following is an example of a valid color priority dictionary: orthologs with a
        white background take precedence over those with a gray background.
        {
            'ortholog': {
                ('#000000', '#FFFFFF'): 1.0,
                ('#000000', '#EDEDED'): 0.0
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
        
        self.color_priority: Dict[str, Dict[Tuple[str, str], float]] = {}
        
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

    def set_color_priority(
        self,
        new_color_priority: Dict[str, Dict[Tuple[str, str], float]],
        recolor_unprioritized_entries: Union[str, Dict[str, Tuple[str, str]]] = False,
        color_associated_compounds: Literal['high', 'low', 'average'] = None,
        colormap: Colormap = None
    ) -> None:
        """
        Set the color_priority attribute. Entry elements in the children attribute are automatically
        reordered.
        
        Though entries (e.g., KOs and compounds) can have multiple Graphics elements due to multiple
        occurrences on a map, Entry Graphics must all have the same foreground and background colors
        if they are to be ordered. Entries with higher priority fg/bg colors are placed first in the
        children attribute and in KGML files, and they are rendered in the foreground of the map.
        The lowest priority entries are always those without fg/bg colors defined in the
        color_priority attribute; these entries are placed first in the Pathway and KGML files, and
        are rendered in the background of the map and thus can be overlaid by higher priority
        entries.
        
        Parameters
        ==========
        new_color_priority : Dict[str, Dict[Tuple[str, str], float]]
            This dictionary is the basis of the color_priority attribute.
            
            It has the same structure as the color_priority attribute. Outer dict keys are Entry
            types, corresponding to the possible values of the type attribute of the Entry class,
            e.g., 'ortholog', 'compound'. Inner dict keys are length-2 tuples of fgcolor and bgcolor
            hex codes, respectively. Inner dict values are non-negative numbers indicating the
            priority of entries with the given fg/bg colors. Higher numbers indicate higher priority
            colors.
            
            What is actually used to set color_priority is a deep copy of the argument in which
            fg/bg color combinations (entries in each inner dict) are reordered in ascending order
            of priority value so that the lowest priority colors appear first in each inner dict.
            The order of Entry types (outer dict entries) does not not change, so if, for example,
            'ortholog' appears before 'compound' in the dict keys, then ortholog entries will occur
            before compound entries in the KGML file, and compounds can be drawn over orthologs.
        
        recolor_unprioritized_entries : Union[str, Dict[str, Tuple[str, str]]], False
            Recolor unprioritized entries, either automatically with a string argument, or with a
            custom dictionary argument for fine-tuning foreground and background colors by Entry
            type.
            
            The valid string arguments for automatic recoloring are 'w' and 'g'. In global maps,
            automatic recoloring is tailored to reaction lines and compound circles. 'w' erases
            unprioritized lines and circles by coloring them white. 'g' colors them a light gray
            (#E0E0E0) consistent with other "unidentified" reactions in the base map. In overview
            maps, automatic recoloring is likewise tailored to arrows and compound circles. However,
            here 'w' colors unprioritized arrows black, consistent with other "unidentified"
            reactions in the base map, and circles' backgrounds white and their (foreground) borders
            black. 'g' colors arrows and circles' backgrounds light gray. In standard maps,
            automatic recoloring is tailored to ortholog boxes, not arrows, and compound circles.
            'w' color unprioritized boxes' and circles' backgrounds white and their (foreground)
            text black. 'g' colors boxes' and circles' backgrounds light gray.
            
            To fine-tune fg/bg colors of unprioritized entries, the arg is a custom dict. Keys are
            Entry type (e.g., 'ortholog', 'compound') and values are length-2 tuples of color hex
            codes for fg and bg colors, respectively. This is shown in the following example, which
            sets the fg color of unprioritized orthologs to dark gray and the bg to light gray, and
            the fg of unprioritized compounds to black and the bg to white.
            {
                'ortholog': ('#A9A9A9', '#E0E0E0'),
                'compound': ('#000000', '#FFFFFF')
            }
        
        color_associated_compounds : Literal['high', 'low', 'average'], None
            Automatically set the background color of compound entries based on the color priority
            of ortholog entries involving the compounds. By default, compounds are circles, and
            orthologs are arrows on global/overview maps and boxes on standard maps.
            
            An argument of 'high' or 'low' sets the compound background color to the bg color of the
            ortholog with the highest or lowest priority fg/bg color combination. 'average' sets the
            bg color to the average bg color of orthologs with prioritized colors -- unprioritized
            orthologs are not taken into account. 'average' should only be used if priority values
            are normalized to the interval [0, 1] and can be converted to a color given by the
            colormap argument. The average priority value of the orthologs with prioritized colors
            is mapped to a bg color for the compound Entry.
            
            Automatically colored compound entries are added to the color_priority attribute, and
            Entry elements in the children attribute are reordered accordingly. Compound entries
            that are already in the color_priority attribute are exempt from recoloring and given
            higher priority than automatically recolored compound entries.
        
        colormap : matplotlib.colors.Colormap, None
            If 'average' is used as the color_associated_compounds argument, a colormap must be
            provided to map averaged priority values on the interval [0, 1] to a background color
            for compound entries.
        """
        # Check that new_color_priority only contains positive priority values.
        for new_type_color_priority in new_color_priority.values():
            for priority in new_type_color_priority.values():
                assert priority >= 0

        # Make the color_priority dict, reordering colors from lowest to highest priority.
        color_priority = {}
        for entry_type, new_type_color_priority in new_color_priority.items():
            color_priority[entry_type] = type_color_priority = {}
            for colors, priority in sorted(
                new_type_color_priority.items(), key=lambda item: item[1]
            ):
                type_color_priority[colors] = priority
        self.color_priority = color_priority
        
        # Reorder Entry elements in the children attribute from lowest to highest priority.
        unprioritized_entry_uuids = self.order_entries_by_color_priority()
        
        if recolor_unprioritized_entries:
            if isinstance(recolor_unprioritized_entries, str):
                assert recolor_unprioritized_entries in ('w', 'g')
                
                # Recolor orthologs.
                if recolor_unprioritized_entries == 'w':
                    if self.is_overview_map:
                        color_hex_code = '#000000'
                    else:
                        color_hex_code = '#FFFFFF'
                elif recolor_unprioritized_entries == 'g':
                    color_hex_code = '#E0E0E0'
                self.recolor_unprioritized_ortholog_entries(
                    unprioritized_entry_uuids, color_hex_code
                )
                
                # Recolor compounds.
                if recolor_unprioritized_entries == 'w':
                    color_hex_code = '#FFFFFF'
                elif recolor_unprioritized_entries == 'g':
                    color_hex_code = '#E0E0E0'
                self.recolor_unprioritized_compound_entries(
                    unprioritized_entry_uuids, color_hex_code
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
                if recolor_unprioritized_entries == 'w':
                    color_hex_code = '#FFFFFF'
                elif recolor_unprioritized_entries == 'g':
                    color_hex_code = '#E0E0E0'
                self.recolor_unprioritized_compound_entries(
                    unprioritized_entry_uuids, color_hex_code
                )
            else:
                self.recolor_unprioritized_entries(
                    unprioritized_entry_uuids, recolor_unprioritized_entries
                )
    
    def recolor_unprioritized_ortholog_entries(
        self,
        unprioritized_entry_uuids: List[str],
        color_hex_code: str
    ) -> None:
        """
        Recolor orthologs without a color priority. Orthologs are expected to be lines in global
        and overview maps and boxes in standard maps.
        
        In global and overview maps, the foreground of lines is colored and the background is made
        white. In standard maps, the background (fill) of boxes is colored and the foreground (text)
        is made black.
        
        Parameters
        ==========
        unprioritized_entry_uuids : List[str]
            List of UUIDs of all entries without a color priority.
        
        color_hex_code : str
            Hex code of the color to give ortholog graphics.
        """
        if self.is_global_map or self.is_overview_map:
            self.recolor_unprioritized_entries(
                unprioritized_entry_uuids, {'ortholog': (color_hex_code, '#FFFFFF')}
            )
        else:
            self.recolor_unprioritized_entries(
                unprioritized_entry_uuids, {'ortholog': ('#000000', color_hex_code)}
            )

    def recolor_unprioritized_compound_entries(
        self,
        unprioritized_entry_uuids: List[str],
        color_hex_code: str
    ) -> None:
        """
        Recolor compounds, expected to be circles, without a color priority.
        
        In global maps, both the fill and border are colored. In overview and standard maps, the
        fill is colored, and the border is made black.
        
        Parameters
        ==========
        unprioritized_entry_uuids : List[str]
            List of UUIDs of all entries without a color priority.
            
        color_hex_code : str
            Hex code of the color to give compound graphics.
        """
        if self.is_global_map:
            self.recolor_unprioritized_entries(
                unprioritized_entry_uuids, {'compound': (color_hex_code, color_hex_code)}
            )
        else:
            self.recolor_unprioritized_entries(
                unprioritized_entry_uuids, {'compound': ('#000000', color_hex_code)}
            )
    
    def recolor_unprioritized_entries(
        self,
        unprioritized_entry_uuids: List[str],
        type_colors: Dict[str, Tuple[str, str]]
    ) -> None:
        """
        Entries without a color priority are recolored by Entry type.
        
        Parameters
        ==========
        unprioritized_entry_uuids : List[str]
            List of UUIDs of all entries without a color priority.
        
        type_colors : Dict[str, Tuple[str, str]]
            Keys are Entry type (e.g., 'ortholog', 'compound') and values are length-2 tuples of
            color hex codes for foreground and background colors, respectively. This is shown in the
            following example, which sets the fg color of unprioritized orthologs to dark gray and
            the bg to light gray, and the fg of unprioritized compounds to black and the bg to
            white.
            {
                'ortholog': ('#A9A9A9', '#E0E0E0'),
                'compound': ('#000000', '#FFFFFF')
            }
        """
        for entry_uuid in unprioritized_entry_uuids:
            entry: Entry = self.uuid_element_lookup[entry_uuid]
            try:
                fgcolor_hex_code, bgcolor_hex_code = type_colors[entry.type]
            except KeyError:
                continue
            for graphics_uuid in entry.children['graphics']:
                graphics: Graphics = self.uuid_element_lookup[graphics_uuid]
                graphics.fgcolor = fgcolor_hex_code
                graphics.bgcolor = bgcolor_hex_code

    def order_entries_by_color_priority(self) -> List[str]:
        """
        Reorder Entry (e.g., 'ortholog' and 'compound') UUIDs by color priority in the children
        attribute of the Pathway object. This determines how entries are ordered in KGML files and
        rendered in maps.
        
        Returns
        =======
        List[str]
            UUIDs of Entry elements without a color priority.
        """
        # Entries have different types ('ortholog', 'compound', etc.). Group entries into two
        # classes. "Qualifying" entries have Entry types in the color priority dict. Other entries
        # do not have Entry types in the dict and are assigned the lowest nominal priority of -1.0.
        # No effort is made to sort these entries in any way.
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
        
        # Sort "qualifying" entries. Loop through each Entry type in the color priority dict.
        for entry_type, type_color_priority in self.color_priority.items():
            # Retrieve each Entry object of the type. Its fg and bg colors determine its priority.
            type_qualifying_entry_uuids = qualifying_entry_uuids[entry_type]
            type_priority_entry_uuids: Dict[float, List[str]] = {}
            for entry_uuid in type_qualifying_entry_uuids:
                entry: Entry = self.uuid_element_lookup[entry_uuid]
                # Ensure that all of the Entry Graphics elements have the same fg and bg colors.
                fgcolors: List[str] = []
                bgcolors: List[str] = []
                for graphics_uuid in entry.children['graphics']:
                    graphics_element: Graphics = self.uuid_element_lookup[graphics_uuid]
                    fgcolors.append(graphics_element.fgcolor)
                    bgcolors.append(graphics_element.bgcolor)
                if len(set(fgcolors)) != 1 or len(set(bgcolors)) != 1:
                    raise AssertionError(
                        "The Graphics elements in the Entry with the following UUID do not have "
                        "consistent foreground and background colors, which is required for "
                        f"ordering entries based on color: {entry_uuid}"
                    )
                try:
                    priority = type_color_priority[(fgcolors[0], bgcolors[0])]
                except KeyError:
                    # The Entry does not have colors in the priority dictionary.
                    priority = -1.0
                
                try:
                    type_priority_entry_uuids[priority].append(entry_uuid)
                except KeyError:
                    type_priority_entry_uuids[priority] = [entry_uuid]
                    
            # Add the reordered UUIDs of the Entry type to the new list of Entry UUIDs and to the
            # dict mapping priority values to UUIDs of entries of all types.
            for priority, entry_uuids in sorted(type_priority_entry_uuids.items()):
                reordered_entry_uuids += entry_uuids
            
            try:
                unprioritized_entry_uuids += type_priority_entry_uuids[-1.0]
            except KeyError:
                pass
            
        self.children['entry'] = reordered_entry_uuids

        return unprioritized_entry_uuids

    def color_associated_compounds(
        self,
        transfer: Literal['high', 'low', 'average'],
        colormap: Colormap = None
    ) -> None:
        """
        Set the color of compound entries based on the color priority of ortholog entries involving
        the compounds. By default, compounds are circles, and orthologs are arrows on
        global/overview maps and boxes on standard maps.
        
        If the map is global, color both the background (interior) and foreground (border) of the
        circle. Otherwise, just color the interior and leave the border black.
        
        Parameters
        ==========
        transfer : Literal['high', 'low', 'average']
            An argument of 'high' or 'low' sets the compound color to the bg color of the ortholog
            with the highest or lowest priority fg/bg color combination. 'average' sets the compound
            color to the average bg color of orthologs with prioritized colors -- unprioritized
            orthologs are not taken into account. 'average' should only be used if priority values
            are normalized to the interval [0, 1] and can be converted to a color given by the
            colormap argument. The average priority value of the orthologs with prioritized colors
            is mapped to a color for the compound Entry.
            
            Compound entries that are already in the color_priority attribute are exempt from
            recoloring and given higher priority than recolored compound entries.
        
        colormap : matplotlib.colors.Colormap, None
            If 'average' is used as an argument to transfer, a colormap must be provided to map
            averaged priority values on the interval [0, 1] to a background color for compound
            entries.
        """
        # Make Reaction elements searchable by name (KEGG IDs). Reaction elements link Compound
        # elements to ortholog Entry elements.
        name_reaction: Dict[str, Reaction] = {}
        for entry_uuid in self.children['reaction']:
            reaction: Reaction = self.uuid_element_lookup[entry_uuid]
            name_reaction[reaction.name] = reaction
        
        # For each compound Entry with associated color-prioritized ortholog entries, record the
        # colors and priorities of these entries.
        compound_uuid_color_priorities: Dict[str, List[Tuple[str, float]]] = {}
        # Loop through each ortholog Entry.
        for entry_uuid in self.children['entry']:
            entry: Entry = self.uuid_element_lookup[entry_uuid]
            if entry.type != 'ortholog':
                continue
            
            # Ensure that all of the ortholog Entry Graphics elements have the same fg/bg colors.
            fgcolors: List[str] = []
            bgcolors: List[str] = []
            for graphics_uuid in entry.children['graphics']:
                graphics: Graphics = self.uuid_element_lookup[graphics_uuid]
                fgcolors.append(graphics.fgcolor)
                bgcolors.append(graphics.bgcolor)
            if len(set(fgcolors)) != 1 or len(set(bgcolors)) != 1:
                raise AssertionError(
                    "The Graphics elements in the ortholog Entry with the following UUID do not "
                    "have consistent foreground and background colors, which is required for "
                    f"ordering entries based on color: {entry_uuid}"
                )
            
            try:
                priority = self.color_priority['ortholog'][(fgcolors[0], bgcolors[0])]
            except KeyError:
                # Unprioritized ortholog entries do not affect the color of associated compounds.
                continue
            
            reaction_name = entry.reaction
            if reaction_name is None:
                # The ortholog is not associated with a reaction.
                continue
            
            try:
                reaction = name_reaction[reaction_name]
            except KeyError:
                # No Reaction element is present with the name of the ortholog reaction.
                continue
            
            if graphics.type == 'line':
                ortholog_color = fgcolors[0]
            else:
                ortholog_color = bgcolors[0]
            for substrate_uuid in reaction.children['substrate']:
                substrate: Substrate = self.uuid_element_lookup[substrate_uuid]
                split_substrate_name = substrate.name.split(':')
                assert len(split_substrate_name) == 2
                for compound_entry in self.kegg_id_element_lookup[split_substrate_name[1]]:
                    compound_entry: Entry
                    compound_uuid = compound_entry.uuid
                    try:
                        compound_uuid_color_priorities[compound_uuid].append(
                            (ortholog_color, priority)
                        )
                    except KeyError:
                        compound_uuid_color_priorities[compound_uuid] = [(ortholog_color, priority)]
            for product_uuid in reaction.children['product']:
                product: Product = self.uuid_element_lookup[product_uuid]
                split_product_name = product.name.split(':')
                assert len(split_product_name) == 2
                for compound_entry in self.kegg_id_element_lookup[split_product_name[1]]:
                    compound_entry: Entry
                    compound_uuid = compound_entry.uuid
                    try:
                        compound_uuid_color_priorities[compound_uuid].append(
                            (ortholog_color, priority)
                        )
                    except KeyError:
                        compound_uuid_color_priorities[compound_uuid] = [(ortholog_color, priority)]
                            
        # Make compound entries searchable by ID, which should be a unique pathway element ID.
        id_compound_entry: Dict[str, Entry] = {}
        for entry_uuid in self.children['entry']:
            entry: Entry = self.uuid_element_lookup[entry_uuid]
            if entry.type != 'compound':
                continue
            id_compound_entry[entry.id] = entry
        
        # Define functions for finding compound Entry color.
        def _get_high_color(color_priorities: List[Tuple[str, float]]) -> Tuple[str, float]:
            return sorted(color_priorities, key=lambda t: -t[1])[0]
        
        def _get_low_color(color_priorities: List[Tuple[str, float]]) -> Tuple[str, float]:
            return sorted(color_priorities, key=lambda t: t[1])[0]
        
        def _get_average_color(color_priorities: List[Tuple[str, float]]) -> Tuple[str, float]:
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
            raise AssertionError

        # Set compound Entry color.
        for compound_uuid, color_priorities in compound_uuid_color_priorities.items():
            compound: Union[Substrate, Product] = self.uuid_element_lookup[compound_uuid]
            compound_entry: Entry = id_compound_entry[compound.id]
            # Get all of the Graphics elements for the Entry.
            graphics_elements: List[Graphics] = []
            for graphics_uuid in compound_entry.children['graphics']:
                graphics_elements.append(self.uuid_element_lookup[graphics_uuid])
            set_color = True
            for graphics in graphics_elements:
                try:
                    # The compound Entry has already been assigned a color priority, so don't
                    # recolor it automatically.
                    self.color_priority['compound'][(graphics.fgcolor, graphics.bgcolor)]
                    set_color = False
                except KeyError:
                    continue
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
                type_color_priority = self.color_priority['compound']
            except KeyError:
                self.color_priority['compound'] = type_color_priority = {}
            if self.is_global_map:
                type_color_priority[(compound_color, compound_color)] = compound_priority
            else:
                type_color_priority[(graphics.fgcolor, compound_color)] = compound_priority

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
                entries += self.kegg_id_element_lookup[kegg_id]
            except KeyError:
                missing_kegg_ids.append(kegg_id)
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
                        f"The attribute, '{attr}', had a value, '{value}', of unrecognized type "
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
        # with open('/Users/sam/Software/anvio/anvio/data/misc/KEGG/map_images/kgml/1x/ko/ko01100.xml') as f:
        #     bio_pathway = KGML_parser.read(f)
        
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
            kegg.download_org_pathway_image_files(f'{pathway.org}{pathway.number}', self.kegg_dir)
            
        bio_pathway.image = map_filepath
        
        canvas = KGMLCanvas(bio_pathway, **kwargs)
        
        # Non-reactants are compounds that don't participate in reactions. By default, this
        # attribute, the alpha of the bg color of the compound circle rendered from the KGML, is set
        # to 0.3. This allows the color of the circle in the underlying base map to bleed through,
        # which is not desirable.
        canvas.non_reactant_transparency = 1
        
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
        
        # Non-reactants are compounds that don't participate in reactions. By default, this
        # attribute, the alpha of the bg color of the compound circle rendered from the KGML, is set
        # to 0.3. This allows the color of the circle in the underlying base map to bleed through,
        # which is not desirable.
        canvas.non_reactant_transparency = 1
        
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
        
        # Non-reactants are compounds that don't participate in reactions. By default, this
        # attribute, the alpha of the bg color of the compound circle rendered from the KGML, is set
        # to 0.3. This allows the color of the circle in the underlying base map to bleed through,
        # which is not desirable.
        canvas.non_reactant_transparency = 1
        
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
