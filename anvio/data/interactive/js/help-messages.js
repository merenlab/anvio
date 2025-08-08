/**
 * Constants for tooltips
 *
 *  Authors: A. Murat Eren <a.murat.eren@gmail.com>
 *           Ozcan Esen
 *
 * Copyright 2015-2021, The anvi'o project (http://anvio.org)
 *
 * Anvi'o is a free software. You can redistribute this program
 * and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * You should have received a copy of the GNU General Public License
 * along with anvi'o. If not, see <http://opensource.org/licenses/GPL-3.0>.
 *
 * @license GPL-3.0+ <http://opensource.org/licenses/GPL-3.0>
 */

var help_contents = {
	'load-state-button': 'Load a previously stored visual state from the profile database',
    'save-state-button': 'Save a snapshot of all user settings for layers into the profile database',
    'layers-tab': 'Access to all visualization options',
    'bins-tab': 'Create, load, and save collection of splits. Access completion and redundancy estimates.',
    'tooltips-tab': 'Information on splits under the mouse pointer while browsing the display',
    'search-tab': 'Search splits, highlight them on the display, manipulate matching splits',
    'drawing-type': "Choose one of the two displays anvi'o supports: circular, and perpendicular tree",
    'order-by': "Order splits based on an anvi'o clustering",
    'view': "Select a data display",
    'draw-angle': "Start and the end degrees for the circular tree",
    'tree-radius': "Size of the circular tree from the center. Size will be automatically determined if this is 0",
    'tree-height': "Height of the perpendicular tree",
    'tree-width': "Width of the perpendicular tree",
    'layer-margins': "The distance between each layer on the tree",
    'bins-layer-height': "Size of the color bar that appears for each bin",
    'custom-layer-margins': "Set layer margin values for each layer separately",
    'edge-length-norm': "Log-normalize edge lengths of the tree",
    'show-grids': "Use grids instead of panels to identify bins",
    'layer-name': "Name of the layer.",
    'layer-color': "Color of the layer (two colors will be required for intensity layers)",
    'layer-type': "Layer type: numeric layers can be represented as bars, or heatmaps",
    'layer-norm': "Normalization function for numeric layers",
    'layer-height': "Layer height in pixels",
    'layer-margin': "Margin with the previous layer",
    'layer-min': "A value below which will be assumed zero during vizualization",
    'layer-max': "A value above which will be treated as maximum during vizualization",
    'edit-multiple-layers': "Select layers by their names using the input box belpw, or by clicking checkboxes next to them, and edit attributes of all selected layers by changing corresponding input items.",
    'draw-button': "This draws the tree on the right side of the screen. It may take a while for large datasets",
    "zoom-in": "Zoom in",
    "zoom-out": "Zoom out",
    "return-to-natural": "Retun to the natural (1:1) view", 
    "save-svg": "Download the display as an SVG file",
    'search-expression': "Build your search expression by selecting proper items from combo boxes and typing in the search term",
    'search-layer': "In which layer your search term should be searched for?",
    'search-value': "Search term. It is case sensitive.",
    'search-operator': "The search operator to specify matching criterion",
    'search-button': "Perform search",
    'search-show-results': "Click to show matchin splits and values down below in text format",
    'search-highlight': "Highlight search results on the tree",
    'search-color': "Highlight color for search results",
    'search-clear': "Clear search highlights from the tree",
    'search-append': "Append matching splits the to the active bin (see which bin is active from the Bins tab)",
    'search-remove': "Remove matching splits from the acrive bin",
    'max-font-size': "Maximum font size for the text layer",
    '': "",
    '': "",
}

// default tooltip placement is 'top', add new entry below if you want to manually set the placement.

var tooltip_placements = {
    'layers-tab': 'bottom',
    'bins-tab': 'bottom',
    'tooltips-tab': 'bottom',
    'search-tab': 'bottom',
    'search-expression': 'right',
}
