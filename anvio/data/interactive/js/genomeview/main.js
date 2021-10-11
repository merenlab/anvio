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

var current_state_name = null
var alignToGC = null;
var stateData = {};
var settings = {} // packaged obj sent off to GenomeDrawer
var xDisps = {};
var renderWindow = [];
var genomeMax = 0;
var yOffset = 0 // vertical space between additional data layers
var xDisplacement = 0; // x-offset of genome start, activated if genome labels are shown
var scaleFactor = 1; // widths of all objects are scaled by this value to zoom in/out
var arrowStyle = 1; // gene arrow cosmetics. 1 (default) = 'inspect-page', 2 = thicker arrows, 3 = pentagon, 4 = rect
var maxGroupSize = 2 // used to calculate group height. base of 1 as each group will contain at minimum a genome layer + group ruler.
var genomeLabelSize = 15; // font size of genome labels
var marginTop = 20; // vertical margin at the top of the genome display
var labelSpacing = 30;  // spacing default for genomeLabel canvas
var geneLabelSize = 40; // gene label font size
var spacing = 50; // vertical spacing between genomes
var scaleInterval = 100; // nt scale intervals
var adlPtsPerLayer = 10000; // number of data points to be subsampled per ADL. TODO: more meaningful default?
var showLabels = true; // show genome labels?
var showGeneLabels = true; // show gene labels?
var dynamicScaleInterval = true; // if true, scale interval automatically adjusts to zoom level
var percentScale = false; // if true, scale measured in proportions (0-1) of total sequence breadth rather than NT ranges.
var geneLabelPos = "above"; // gene label position; one of "above", "slanted", "inside"
var thresh_count_gene_colors = 4; // min # occurences of annotation for filtering gene color table
var order_gene_colors_by_count = true; // if true, order annotations on gene color table by descending order of count, otherwise order alphabetically
var filter_gene_colors_to_window = false; // if true, only display gene colors in current render window, otherwise show all gene colors in split
var mainCanvasHeight;
var canvas;
var genomeLabelsCanvas;
var brush;
var drawer
var color_db;
var counts; // stores # occurences for each category in the current function type
var genomeData;

$(document).ready(function () {
  toastr.options = {
    "closeButton": true,
    "debug": false,
    "newestOnTop": true,
    "progressBar": false,
    "positionClass": "toast-top-right",
    "preventDuplicates": false,
    "onclick": null,
    "showDuration": "500",
    "hideDuration": "2000",
    "timeOut": "12000",
    "extendedTimeOut": "1000",
    "showEasing": "swing",
    "hideEasing": "linear",
    "showMethod": "fadeIn",
    "hideMethod": "fadeOut",
  }
  initData();
  loadAdditionalDataLayers()
  processState('default', stateData)
  loadAll('init');
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
      settings['genomeData'] = genomeData
      settings['genomeData']['genomes'] = genomes
    }
  });
}

function loadAdditionalDataLayers(){
  $.ajax({
    type: 'POST',
    cache: false,
    url: '/data/get_genome_view_adl',
    async: false,
    success: function (resp) {
      settings['additional-data-layers'] = resp
      settings['additional-data-layers']['layers'].push('ruler') // add ruler by default
      settings['group-layer-order'] = ['Genome', 'Ruler']

      if (settings['additional-data-layers']['layers'].includes('Coverage')) {
        settings['group-layer-order'].unshift('Coverage')
        maxGroupSize += 1 // increase group size if coverage layer exists
      }

      if (settings['additional-data-layers']['layers'].includes('GC_content')) {
        settings['group-layer-order'].unshift('GC_Content')
        maxGroupSize += 1 // increase group size if GC layer exists
      }
    }
  });
}

function loadState() {

  var defer = $.Deferred();
    $('#modLoadState').modal('hide');
    if ($('#loadState_list').val() == null) {
        defer.reject();
        return;
    }

    var state_name = $('#loadState_list').val();

  $.ajax({
    type: 'GET',
    cache: false,
    url: '/state/get/' + state_name,
    success: function (response) {
      try {
        processState(state_name, response['content']);
        loadAll('reload')
      } catch (e) {
        console.error("Exception thrown", e.stack);
        toastr.error('Failed to parse state data, ' + e);
        defer.reject();
        return;
      }
    },
    error: function(resp){
      console.log(resp)
    }
  })
}

function serializeSettings() {
  // TODO same process as the serializeSettings() function for anvi-interactive
  // first we run through all of the UI element default values and store them as state
  // then we update them as necessary below in processState
  let state = {}

  state['genome-spacing'] = $('#genome_spacing').val()
  state['order-method'] = $('#genome_order_select').val()
  state['dynamic-scale-interval'] = $('#show_dynamic_scale_box').is(':checked')
  state['genome-scale-interval'] = $('#genome_scale_interval').val()
  state['arrow-stype'] = $('#arrow_style').val()
  state['group-layer-order'] = settings['group-layer-order']
  state['bookmarks'] = settings['display']['bookmarks']
  state['genome-order'] = settings['genomeData']['genomes']
  state['gene-link-style'] = $('#link_style').val()
  state['gene-shade-style'] = $('#shade_by').val()
  state['show-genome-labels'] = $('#show_genome_labels_box').is(':checked')
  state['genome-label-size'] = $('#genome_label').val()
  state['genome-label-color'] = $('#genome_label_color').attr(':color')
  state['show-gene-labels'] = $('#show_gene_labels_box').is(':checked')
  state['gene-label-size'] = $('#gene_label').val()
  state['gene-label-color'] = $('#gene_label_color').attr(':color')
  state['gene-text-position'] = $('#gene_text_pos').val()
  state['gc-window-size'] = $('#gc_window_size').val()
  state['gc-step-size'] = $('#gc_step_size').val()
  state['gc-overlay-color'] = $('#gc_overlay_color').attr(':color')
  state['gene-color-order'] = $('#gene_color_order').val()
  state['annotation-color-dict'] = []

  $('.annotation_color').each((idx, row) => {
    let color = $(row).attr('color')
    let id = ($(row).attr('id').split('_')[1])
    state['annotation-color-dict'].push({
      id : id,
      color : color
    })
  })

  return state
}

function processState(stateName, stateData) {
  settings['state-name'] = stateName

  calculateMaxGenomeLength()
  if (stateData.hasOwnProperty('group-layer-order')) {
    settings['group-layer-order'] = stateData['group-layer-order']
  }

  if (stateData.hasOwnProperty('additional-data-layers')) {
    settings['additional-data-layers'] = stateData['additional-data-layers']
  }

  $('#tbody_additionalDataLayers').html('') // clear div before reprocess
  settings['group-layer-order'].map(layer => {
    buildGroupLayersTable(layer)
  })

  if (stateData.hasOwnProperty('genome-order-method')) {
    settings['genome-order-method'] = stateData['genome-order-method']
    settings['genome-order-method'].forEach(orderMethod => {
      $('#genome_order_select').append((new Option(orderMethod["name"], orderMethod["name"]))) // set display + value of new select option.
    })
  } else {
    generateMockGenomeOrder()
    settings['genome-order-method'].forEach(orderMethod => {
      $('#genome_order_select').append((new Option(orderMethod["name"], orderMethod["name"]))) // set display + value of new select option.
    })
  }

  if (stateData.hasOwnProperty('genome-order')){
    settings['genomeData']['genomes'] = stateData['genome-order']
  }

  if (stateData.hasOwnProperty('display')) {
    settings['display'] = stateData['display']
    if(stateData['display']['additional-data-layers']['coverage']){
      settings['display']['additional-data-layers']['coverage'] = stateData['display']['additional-data-layers']['coverage']
    }
    if(stateData['display']['additional-data-layers']['gc-content']){
      settings['display']['additional-data-layers']['gc-content'] = stateData['display']['additional-data-layers']['gc-content']
    }
  } else {
    settings['display'] = {}
    settings['display']['additional-data-layers'] = {}
  }

  if (stateData.hasOwnProperty('genome-spacing')){
    settings['display']['genome-spacing'] = stateData['genome-spacing']
  }

  if (stateData['display'] && stateData['display'].hasOwnProperty('bookmarks')) {
    settings['display']['bookmarks'] = stateData['bookmarks']
    settings['display']['bookmarks'].map(bookmark => {
      $('#bookmarks-select').append((new Option(bookmark['name'], [bookmark["start"], bookmark['stop']])))
    })
    respondToBookmarkSelect()
  } else {
    calculateMaxGenomeLength() // remove later, only here to set max length for bookmark
    settings['display']['bookmarks'] = generateMockBookmarks()
    settings['display']['bookmarks'].map(bookmark => {
      $('#bookmarks-select').append((new Option(bookmark['name'], [bookmark["start"], bookmark['stop']])))
    })
    respondToBookmarkSelect() // set listener for user bookmark selection
  }
}

function loadAll(loadType) {
  canvas = new fabric.Canvas('myCanvas');
  canvas.clear() // clear existing canvas if loadAll is being called from loadState
  canvas.setWidth(VIEWER_WIDTH * 0.85);

  $('.container').css({ 'height': VIEWER_HEIGHT + 'px', 'overflow-y': 'auto' })
  xDisplacement = showLabels ? 120 : 0;
  for (genome of settings['genomeData']['genomes']) {
    xDisps[genome[0]] = xDisplacement;
  }

  calculateMaxGenomeLength()

  if (showGeneLabels && arrowStyle != 3) {
    marginTop = 60;
    spacing = settings['display']['genome-spacing'] ? settings['display']['genome-spacing'] : 200; // TODO maybe we refactor this out into a setSpacing() method for clarity?
    $("#genome_spacing").val(spacing);
  }

  $('#gene_color_order').append($('<option>', {
    value: 'Source',
    text: 'Source'
  }));
  for (fn of getFunctionalAnnotations()) {
    $('#gene_color_order').append($('<option>', {
      value: fn,
      text: fn
    }));
  }

  color_db = $('#gene_color_order').val();

  if(loadType == 'init'){
    buildGenomesTable(settings['genomeData']['genomes'], 'alphabetical') // hardcode order method until backend order data is hooked in
    drawScale();
    setEventListeners()
    generateColorTable(fn_colors = null, fn_type = color_db);
  }

  brush.extent([parseInt($('#brush_start').val()), parseInt($('#brush_end').val())]);
  brush(d3.select(".brush"));
  updateRenderWindow();

  console.log('Sending this data obj to GenomeDrawer', settings)
  drawer = new GenomeDrawer(settings)
  drawer.draw()
}
