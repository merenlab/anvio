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
var yOffset = 0 // y-location for the current additional data layer
var scaleFactor = 1; // widths of all objects are scaled by this value to zoom in/out
var maxGroupSize = 2 // used to calculate group height. base of 2 as each group will contain at minimum a genome layer + group ruler.
var marginTop = 20; // vertical margin at the top of the genome display
var groupLayerPadding = 10 // padding between each layer in a given genome group
var groupMargin = 100 // space between each genome group
var spacing = 50; // multiplied by maxGroupSize to determine group height allocation
var labelCanvasWidth = 70 // default width of genomeLabelCanvas
var percentScale = false; // if true, scale measured in proportions (0-1) of total sequence breadth rather than NT ranges.
var order_gene_colors_by_count = true; // if true, order annotations on gene color table by descending order of count, otherwise order alphabetically
var filter_gene_colors_to_window = false; // if true, only display gene colors in current render window, otherwise show all gene colors in split
var firstDraw = true // flag to determine whether to set zoom to initial zoom level
var canvas;
var labelCanvas;
var brush;
var drawer
var color_db; // functional annotation type currently being used to color genes
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
      settings['display'] = {}
      settings['display']['hidden'] = {}
      settings['display']['colors'] = {}
      settings['display']['colors']['genes'] = {}
      settings['display']['colors']['genes']['annotations'] = {}
      settings['display']['colors']['Batch'] = '#FFFFFF';
      settings['display']['layers'] = {}
      settings['display']['labels'] = {}
      settings['display']['labels']['set-labels'] = {}
      settings['display']['labels']['gene-sets'] = {}
      settings['display']['layers']['Ruler'] = true
      settings['display']['layers']['Genome'] = true

      if (settings['additional-data-layers']['layers'].includes('Coverage')) {
        settings['group-layer-order'].unshift('Coverage')
        settings['display']['layers']['Coverage'] = true
        settings['display']['colors']['Coverage'] = '#000000'
        maxGroupSize += 1
      }

      if (settings['additional-data-layers']['layers'].includes('GC_content')) {
        settings['group-layer-order'].unshift('GC_Content')
        settings['display']['layers']['GC_Content'] = true
        settings['display']['colors']['GC_Content'] = '#000000'
        maxGroupSize += 1
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

function showSaveStateWindow(){
  $.ajax({
      type: 'GET',
      cache: false,
      url: '/state/all',
      success: function(state_list) {
          $('#saveState_list').empty();

          for (let state_name in state_list) {
              var _select = "";
              if (state_name == current_state_name)
              {
                  _select = ' selected="selected"';
              }
              $('#saveState_list').append('<option ' + _select + '>' + state_name + '</option>');
          }

          $('#modSaveState').modal('show');
          if ($('#saveState_list').val() === null) {
              $('#saveState_name').val('default');
          } else {
              $('#saveState_list').trigger('change');
          }
      },
      error: function(error){
        console.log('got an error', error)
      }
  });
}

function showLoadStateWindow(){
  $.ajax({
      type: 'GET',
      cache: false,
      url: '/state/all',
      success: function(state_list) {
          $('#loadState_list').empty();

          for (let state_name in state_list) {
              $('#loadState_list').append('<option lastmodified="' + state_list[state_name]['last_modified'] + '">' + state_name + '</option>');
          }

          $('#modLoadState').modal('show');
      }
  });
}

function saveState()
{
  var name = $('#saveState_name').val();

  if (name.length==0) {
      $('#saveState_name').focus();
      return;
  }

  var state_exists = false;

  $.ajax({
      type: 'GET',
      cache: false,
      async: false,
      url: '/state/all',
      success: function(state_list) {
          for (let state_name in state_list) {
              if (state_name == name)
              {
                  state_exists = true;
              }
          }

      }
  });

  if (state_exists && !confirm('"' + name + '" already exist, do you want to overwrite it?')) {
      return;
  }

  $.ajax({
      type: 'POST',
      cache: false,
      url: '/state/save/' + name,
      data: {
          'content': JSON.stringify(serializeSettings())
      },
      success: function(response) {
          if (typeof response != 'object') {
              response = JSON.parse(response);
          }

          if (response['status_code']==0)
          {
              toastr.error("Failed, Interface running in read only mode.");
          }
          else if (response['status_code']==1)
          {
              // successfull
              $('#modSaveState').modal('hide');

              current_state_name = name;
              toastr.success("State '" + current_state_name + "' successfully saved.");
          }
      }
  });
}

function serializeSettings() {
  // TODO same process as the serializeSettings() function for anvi-interactive
  // first we run through all of the UI element default values and store them as state
  // then we update them as necessary below in processState
  let state = {}

  state['group-layer-order'] = settings['group-layer-order']
  state['genome-order'] = settings['genomeData']['genomes']
  state['display'] = settings['display']
  state['display']['bookmarks'] = settings['display']['bookmarks']
  state['display']['metadata'] = settings['display']['metadata']

  state['display']['order-method'] = $('#genome_order_select').val()
  state['display']['dynamic-scale-interval'] = $('#show_dynamic_scale_box').is(':checked') // if true, scale interval automatically adjusts to zoom level
  state['display']['genome-scale-interval'] = parseInt($('#genome_scale_interval').val())
  state['display']['genome-spacing'] = $('#genome_spacing').val()
  state['display']['genome-margin'] = $('#genome_margin').val()
  state['display']['arrow-style'] = $('#arrow_style').val()
  state['display']['gene-link-style'] = $('#link_style').val()
  state['display']['gene-shade-style'] = $('#shade_by').val()
  state['display']['show-genome-labels'] = $('#show_genome_labels_box').is(':checked')
  state['display']['genome-label-size'] = $('#genome_label').val() 
  state['display']['colors']['genome-label'] = $('#genome_label_color').attr('color')
  state['display']['show-gene-labels'] = $('#show_gene_labels_box').is(':checked')
  state['display']['gene-label-size'] = $('#gene_label').val() 
  state['display']['colors']['gene-label'] = $('#gene_label_color').attr('color')
  state['display']['gene-text-position'] = $('#gene_text_pos').val() // gene label position; one of "above", "inside"
  state['display']['gene-text-angle'] = $('#gene_text_angle').val()
  state['display']['gc-window-size'] = $('#gc_window_size').val()
  state['display']['gc-step-size'] = $('#gc_step_size').val()
  state['display']['gc-overlay-color'] = $('#gc_overlay_color').attr('color')
  state['display']['color_db'] = $('#gene_color_order').val()
  state['display']['gene-label-source'] = $('#gene_label_source').val()
  state['display']['link-gene-label-color-source'] = $('#link_gene_label_color_source').is(':checked') // by default, allow users to display different gene arrow color / gene label source (false)
  state['display']['annotation-color-dict'] = []
  state['display']['viewportTransform'] = canvas.viewportTransform[4];
  state['display']['nt_window'] = percentScale ? [parseFloat($('#brush_start').val()), parseFloat($('#brush_end').val())] : [parseInt($('#brush_start').val()), parseInt($('#brush_end').val())];
  state['display']['scaleFactor'] = scaleFactor;
  state['display']['percentScale'] = percentScale;
  state['display']['xDisps'] = xDisps;
  state['display']['thresh-count-gene-colors'] = $('#thresh_count').val() // min # occurences of annotation for filtering gene color table
  state['display']['adlPtsPerLayer'] = $('#adl_pts_per_layer').val() // number of data points to be subsampled per ADL. TODO: more meaningful default?
  state['display']['user-defined-colors'] = $('#user_defined_colors').is(':checked')

  $('.annotation_color').each((idx, row) => {
    let color = $(row).attr('color')
    let id = ($(row).attr('id').split('_')[1])
    state['display']['annotation-color-dict'].push({
      id : id,
      color : color
    })
  })

  return state
}

function processState(stateName, stateData) {
  settings['state-name'] = stateName
  console.log('processing this state obj', stateData)

  if (stateData.hasOwnProperty('group-layer-order')) {
    settings['group-layer-order'] = stateData['group-layer-order']
  }

  if (stateData.hasOwnProperty('additional-data-layers')) {
    settings['additional-data-layers'] = stateData['additional-data-layers']
  }

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

  if (stateData?.['display']) {
    settings['display'] = stateData['display']
  }

  if (stateData?.['display']?.['arrow-style']){
    $('#arrow_style').val(stateData['display']['arrow-style'])
  } else {
    $('#arrow_style').val(1)
  }

  if (stateData?.['display']?.['bookmarks']) {
    settings['display']['bookmarks'].map(bookmark => {
      $('#bookmarks-select').append((new Option(bookmark['name'], [bookmark["start"], bookmark['stop']])))
    })
    respondToBookmarkSelect()
  } else {
    settings['display']['bookmarks'] = []
    respondToBookmarkSelect() // set listener for user bookmark selection
  }

  if (stateData?.['display']?.['link-gene-label-color-source']){
    $('#link_gene_label_color_source').prop('checked', settings['display']['link-gene-label-color-source'])
  }

  $('#tbody_additionalDataLayers').html('') // clear div before reprocess
  settings['group-layer-order'].map(layer => {
    buildGroupLayersTable(layer)
  })

  if(stateData?.['display']?.['percentScale']) {
    percentScale = stateData['display']['percentScale']
  } else {
    percentScale = false;
  }

  if(stateData?.['display']?.['scaleFactor']) {
    scaleFactor = stateData['display']['scaleFactor']
  } else {
    scaleFactor = 1
  }

  if(stateData?.['display']?.['nt_window']) {
    settings['display']['nt_window'] = stateData['display']['nt_window']
  } else {
    settings['display']['nt_window'] = percentScale ? [parseFloat($('#brush_start').val()), parseFloat($('#brush_end').val())] : [parseInt($('#brush_start').val()), parseInt($('#brush_end').val())];
  }

  if(stateData?.['display']?.['viewportTransform']) {
    settings['display']['viewportTransform'] = stateData['display']['viewportTransform']
  }

  if(stateData?.['display']?.['xDisps']) {
    settings['display']['xDisps'] = stateData['display']['xDisps']
  }

  if(stateData?.['display']?.['genome-label-size']) {
    settings['display']['genome-label-size'] = stateData['display']['genome-label-size'];
  } else {
    settings['display']['genome-label-size'] = 15
  }
  $('#genome_label').val(settings['display']['genome-label-size'])

  if(stateData?.['display']?.['gene-label-size']) {
    settings['display']['gene-label-size'] = stateData['display']['gene-label-size']
  } else {
    settings['display']['gene-label-size'] = 40
  }
  $('#gene_label').val(settings['display']['gene-label-size'])

  if(stateData?.['display']?.['gene-text-position']) {
    settings['display']['gene-text-position'] = stateData['display']['gene-text-position']
  } else {
    settings['display']['gene-text-position'] = 'above'
  }
  $('#gene_text_pos').val(settings['display']['gene-text-position'])

  if(stateData?.['display']?.['gene-text-angle']) {
    settings['display']['gene-text-angle'] = stateData['display']['gene-text-angle']
  } else {
    settings['display']['gene-text-angle'] = 0
  }
  $('#gene_text_angle').val(settings['display']['gene-text-angle'])

  if(stateData?.['display']?.hasOwnProperty('dynamic-scale-interval')) {
    settings['display']['dynamic-scale-interval'] = stateData['display']['dynamic-scale-interval']
  } else {
    settings['display']['dynamic-scale-interval'] = true;
  }

  if(stateData?.['display']?.hasOwnProperty('genome-scale-interval')) {
    settings['display']['genome-scale-interval'] = stateData['display']['genome-scale-interval']
  } else {
    settings['display']['genome-scale-interval'] = 100 
  }
  $('#genome_scale_interval').val(settings['display']['genome-scale-interval'])

  if(stateData?.['display']?.hasOwnProperty('show-genome-labels')) {
    settings['display']['show-genome-labels'] = stateData['display']['show-genome-labels']
  } else {
    settings['display']['show-genome-labels'] = true 
  }

  if(stateData?.['display']?.hasOwnProperty('show-gene-labels')) {
    settings['display']['show-gene-labels'] = stateData['display']['show-gene-labels']
  } else {
    settings['display']['show-gene-labels'] = true 
  }

  if(stateData?.['display']?.hasOwnProperty('thresh-count-gene-colors')) {
    settings['display']['thresh-count-gene-colors'] = stateData['display']['thresh-count-gene-colors']
  } else {
    settings['display']['thresh-count-gene-colors'] = 4
  }
  $('#thresh_count').val(settings['display']['thresh-count-gene-colors'])

  if(stateData?.['display']?.hasOwnProperty('adlPtsPerLayer')) {
    settings['display']['adlPtsPerLayer'] = stateData['display']['adlPtsPerLayer']
  } else {
    settings['display']['adlPtsPerLayer'] = 10000 // number of data points to be subsampled per ADL. TODO: more meaningful default?
  }
  $('#adl_pts_per_layer').val(settings['display']['adlPtsPerLayer'])
  
  if(stateData?.['display']?.hasOwnProperty('link-gene-label-color-source')) {
    settings['display']['link-gene-label-color-source'] = stateData['display']['link-gene-label-color-source']
  } else {
    settings['display']['link-gene-label-color-source'] = false
  }

  if(stateData?.['display']?.['colors']?.hasOwnProperty('genome-label')) {
    settings['display']['colors']['genome-label'] = stateData['display']['colors']['genome-label'];
  } else {
    settings['display']['colors']['genome-label'] = '#000000'
  }

  if(stateData?.['display']?.['colors']?.hasOwnProperty('gene-label')) {
    settings['display']['colors']['gene-label'] = stateData['display']['colors']['gene-label'];
  } else {
    settings['display']['colors']['gene-label'] = '#000000'
  }

  if(stateData?.['display']?.hasOwnProperty('genome-spacing')) {
    spacing = parseInt(stateData['display']['genome-spacing']);
  } else {
    spacing = 200
  }
  $("#genome_spacing").val(spacing);

  if(stateData?.['display']?.hasOwnProperty('genome-margin')) {
    groupMargin = parseInt(stateData['display']['genome-margin']);
  } else {
    groupMargin = 100
  }
  $("#genome_margin").val(groupMargin);

  $('#show_genome_labels_box').prop("checked", settings['display']['show-genome-labels']);
  $('#show_gene_labels_box').prop("checked", settings['display']['show-gene-labels']);
  $('#show_dynamic_scale_box').prop("checked", settings['display']['dynamic-scale-interval']);
  $('#link_gene_label_color_source_box').prop("checked", settings['display']['link-gene-label-color-source']);
  $('#user_defined_colors').prop("checked", settings['display']['user-defined-colors']);
}

function loadAll(loadType) {
  if(canvas) canvas.dispose(); // clear old canvas in fabric if we are calling loadAll from loadState (bugs here cause mismatch b/w 'this' and 'canvas')
  canvas = new fabric.Canvas('myCanvas');

  if(labelCanvas) labelCanvas.dispose();
  labelCanvas = new fabric.Canvas('genomeLabelCanvas');
  setLabelCanvas(); // also sets canvas width

  $('.container').css({ 'height': VIEWER_HEIGHT + 'px', 'overflow-y': 'auto' })

  $('#genome_label_color').css('background-color', settings['display']['colors']['genome-label']);
  $('#genome_label_color').attr('color', settings['display']['colors']['genome-label']);
  drawGenomeLabels();

  if(settings?.['display']?.['xDisps']) {
    xDisps = settings?.['display']['xDisps'];
  } else {
    for (genome of settings['genomeData']['genomes']) {
      xDisps[genome[0]] = 0;
    }
  }

  calculateMaxGenomeLength()

  if(loadType == 'init') {
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
    buildGeneLabelsSelect()
  }

  $('#gene_label_color').css('background-color', settings['display']['colors']['gene-label']);
  $('#gene_label_color').attr('color', settings['display']['colors']['gene-label']);

  color_db = settings?.['display']?.['color_db'] ? settings['display']['color_db'] : $('#gene_color_order').val();
  $('#gene_color_order').val(color_db);
  generateColorTable(fn_colors = settings?.display?.colors?.genes?.annotations[color_db], fn_type = color_db);

  buildGenomesTable(settings['genomeData']['genomes'], 'alphabetical') // hardcode order method until backend order data is hooked in
  drawScale();
  if(firstDraw) setEventListeners();
  setCanvasListeners();

  let [start, stop] = settings['display']['nt_window'];
  $('#brush_start').val(start);
  $('#brush_end').val(stop);
  if(loadType == 'reload') {
    canvas.viewportTransform[4] = settings['display']['viewportTransform'];
    canvas.setViewportTransform(canvas.viewportTransform);
  }
  brush.extent([start, stop]);
  brush(d3.select(".brush"));
  updateRenderWindow();
  setLabelCanvas(); // set a second time to adjust to new brush extent

  console.log('Sending this data obj to GenomeDrawer', settings)
  drawer = new GenomeDrawer(settings)
  drawer.draw()
  console.log('draw from loadAll')
  toastr.success(`Successfully loaded from ${settings['state-name']} state`)
}
