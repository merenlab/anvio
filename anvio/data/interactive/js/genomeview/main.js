/**
 * Javascript library for anvi'o genome view
 *
 *  Authors: Isaac Fink <iafink@uchicago.edu>
 *           Matthew Klein <mtt.l.kln@gmail.com>
 *           A. Murat Eren <a.murat.eren@gmail.com>
 *
 * Copyright 2021, The anvi'o project (http://anvio.org)
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

/**
 * File Overview : This file is the entrypoint for genomeview. Here, genomic + state data are retrieved from the backend, processed, and passed to the various other
 * genomeview modules to build out UI and render visualizations. Generally speaking, this file should stay pretty lean and purpose-driven.  
 */

 var VIEWER_WIDTH = window.innerWidth || document.documentElement.clientWidth || document.getElementsByTagName('body')[0].clientWidth;
 var VIEWER_HEIGHT = window.innerHeight || document.documentElement.clientHeight || document.getElementsByTagName('body')[0].clientHeight;

 var canvas;
 var genomeLabelsCanvas;
 var genomeMax = 0;

 var stateData = {};
 var settings = {} // packaged obj sent off to GenomeDrawer
 var mainCanvasHeight;
 var spacing = 50; // vertical spacing between genomes
 var yOffset = 0 // vertical space between additional data layers
 var marginTop = 20; // vertical margin at the top of the genome display
 var xDisplacement = 0; // x-offset of genome start, activated if genome labels are shown
 var showLabels = true; // show genome labels?
 var genomeLabelSize = 15; // font size of genome labels
 var showGeneLabels = true; // show gene labels?
 var geneLabelSize = 40; // gene label font size
 var geneLabelPos = "above"; // gene label position; one of "above", "slanted", "inside"
 var labelSpacing = 30;  // spacing default for genomeLabel canvas
 var scaleInterval = 100; // nt scale intervals
 var dynamicScaleInterval = true; // if true, scale interval automatically adjusts to zoom level
 var adlPtsPerLayer = 10000; // number of data points to be subsampled per ADL. TODO: more meaningful default?
 var scaleFactor = 1; // widths of all objects are scaled by this value to zoom in/out
 var maxGroupSize = 2 // used to calculate group height. base of 1 as each group will contain at minimum a genome layer + group ruler.
 let percentScale = false; // if true, scale measured in proportions (0-1) of total sequence breadth rather than NT ranges.
 var renderWindow = [];
 var brush;
 var drawer 

 var alignToGC = null;

 var arrowStyle = 1; // gene arrow cosmetics. 1 (default) = 'inspect-page', 2 = thicker arrows, 3 = pentagon, 4 = rect

 var color_db;

var genomeData;
var xDisps = {};

$(document).ready(function () {
  initData();
  loadState();
  loadAll();
});

function initData() {
  // initialize the bulk of the data.
  $.ajax({
    type: 'POST',
    cache: false,
    url: '/data/get_genome_view_data',
    async: false,
    success: function (data) {
      genomeData = data;
      console.log("Saved the following data:");
      console.log(data);

      let genomes = Object.entries(genomeData.genomes) // an array of 2d arrays, where each genome[0] is the object key, and genome[1] is the value
      genomeData.genomes = genomes
    }
  });
}

function loadState() {
  // $.ajax({
  //   type: 'GET',
  //   cache: false,
  //   url: '/data/genome_view/state/get/' + state_name,
  //   success: function (response) {
  //     try {
  //       // processState(state_name, response[0]); process actual response from backend
  //       processState(state_name, mockStateData); // process mock state data
  //     } catch (e) {
  //       console.error("Exception thrown", e.stack);
  //       toastr.error('Failed to parse state data, ' + e);
  //       defer.reject();
  //       return;
  //     }
  //     // waitingDialog.hide();
  //   }
  // })
  processState('default', stateData) // moved here until state route is hooked in from backend

}

function serializeSettings() {
  // TODO same process as the serializeSettings() function for anvi-interactive
  // first we run through all of the UI element default values and store them as state
  // then we update them as necessary below in processState
}

function processState(stateName, stateData) {
  calculateMaxGenomeLength()
  if (stateData.hasOwnProperty('group-layer-order')) {
    settings['group-layer-order'] = stateData['group-layer-order']
  } else {
    settings['group-layer-order'] = ['Genome', 'Ruler']
  }

  if (stateData.hasOwnProperty('additional-data-layers')) {
    settings['additional-data-layers'] = stateData['additional-data-layers']
  } else {
    stateData['additional-data-layers'] = []
    generateMockADL()
    settings['additional-data-layers'] = stateData['additional-data-layers']
  }

  // working under the assumption that all genome groups with contain the same additional data layers,
  // we can query the first genome group for specific ADL and go from there
  buildGroupLayersTable('Genome')

  if (stateData['additional-data-layers'][0]['ruler']) {
    buildGroupLayersTable('Ruler')
  }
  if (stateData['additional-data-layers'][0]['coverage']) {
    buildGroupLayersTable('Coverage')
    settings['group-layer-order'].push('Coverage')
    maxGroupSize += 1 // increase group size if coverage layer exists
  }

  if (stateData['additional-data-layers'][0]['gcContent']) {
    buildGroupLayersTable('GC_Content')
    settings['group-layer-order'].push('GC_Content')
    maxGroupSize += 1 // increase group size if GC layer exists
  }

  if (stateData.hasOwnProperty('genome-order-method')) {
    settings['genome-order-method'] = stateData['genome-order-method']
    settings['genome-order-method'].forEach(orderMethod => {
      $('#genome_order_select').append((new Option(orderMethod["name"], orderMethod["name"]))) // set display + value of new select option.
    })
  } else {
    generateMockGenomeOrder()
    settings['genome-order-method'] = stateData['genome-order-method']
    settings['genome-order-method'].forEach(orderMethod => {
      $('#genome_order_select').append((new Option(orderMethod["name"], orderMethod["name"]))) // set display + value of new select option.
    })
  }

  if (stateData.hasOwnProperty('display')) {
    // TODO process
    settings['display'] = stateData['display']
  } else {
    stateData['display'] = {}
    stateData['display']['additionalDataLayers'] = {}
    settings['display'] = stateData['display']
  }
  if (stateData['display'].hasOwnProperty('bookmarks')) {
    settings['display']['bookmarks'] = stateData['bookmarks']
    settings['display']['bookmarks'].map(bookmark => {
      $('#bookmarks-select').append((new Option(bookmark['name'], [bookmark["start"], bookmark['stop']])))
    })
    respondToBookmarkSelect()
  } else {
    calculateMaxGenomeLength() // remove later, only here to set max length for bookmark
    stateData['display']['bookmarks'] = generateMockBookmarks() // for testing purposes
    settings['display']['bookmarks'] = stateData['display']['bookmarks']
    settings['display']['bookmarks'].map(bookmark => {
      $('#bookmarks-select').append((new Option(bookmark['name'], [bookmark["start"], bookmark['stop']])))
    })
    respondToBookmarkSelect() // set listener for user bookmark selection
  }
}

function loadAll() {
  buildGenomesTable(genomeData.genomes, 'alphabetical') // hardcode order method until backend order data is hooked in
  canvas = new fabric.Canvas('myCanvas');
  canvas.setWidth(VIEWER_WIDTH * 0.85);

  $('.container').css({ 'height': VIEWER_HEIGHT + 'px', 'overflow-y': 'auto' })
  xDisplacement = showLabels ? 120 : 0;
  for (genome of genomeData.genomes) {
    xDisps[genome[0]] = xDisplacement;
  }

  // Find max length genome
  calculateMaxGenomeLength()

  drawScale();

  if (showGeneLabels && arrowStyle != 3) {
    marginTop = 60;
    spacing = 200; // TODO maybe we refactor this out into a setSpacing() method for clarity?
    $("#genome_spacing").val(spacing);
  }

  $('#gene_color_order').append($('<option>', {
    value: 'Source',
    text: 'Source'
  }));
  for (fn of getFunctionalAnnotations()) {
    if (!['COG_CATEGORY', 'KEGG_CATEGORY'].includes(fn)) continue; // TODO: support any option
    $('#gene_color_order').append($('<option>', {
      value: fn,
      text: fn
    }));
  }
  color_db = $('#gene_color_order').val();
  generateColorTable(fn_colors = null, fn_type = color_db);

  brush.extent([parseInt($('#brush_start').val()), parseInt($('#brush_end').val())]);
  brush(d3.select(".brush"));
  updateRenderWindow();
  // draw();

  setEventListeners()
  settings = Object.assign(settings, genomeData, stateData)
  drawer = new GenomeDrawer(settings)
  drawer.draw() 
}