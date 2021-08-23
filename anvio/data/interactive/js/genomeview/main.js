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
 * File Overview : This file is the entrypoint for genomeview. Here, genomic + state data are are retrieved from the backend, processed, and passed to the various other
 * genomeview modules to build out UI and render visualizations. Generally speaking, this file should stay pretty lean and purpose-driven.  
 */

$(document).ready(function () {
  initData();
  // loadState();
  processState('default', stateData) // lifted here until loadState is hooked in from backend
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
  $.ajax({
    type: 'GET',
    cache: false,
    url: '/data/genome_view/state/get/' + state_name,
    success: function (response) {
      try {
        // processState(state_name, response[0]); process actual response from backend
        processState(state_name, mockStateData); // process mock state data
      } catch (e) {
        console.error("Exception thrown", e.stack);
        toastr.error('Failed to parse state data, ' + e);
        defer.reject();
        return;
      }
      // waitingDialog.hide();
    }
  })
}

function serializeSettings() {
  // TODO same process as the serializeSettings() function for anvi-interactive
  // first we run through all of the UI element default values and store them as state
  // then we update them as necessary below in processState
}

function processState(stateName, stateData) {
  calculateMaxGenomeLength()
  if (stateData.hasOwnProperty('group-layer-order')) {
    // TODO process
  } else {
    stateData['group-layer-order'] = ['Genome', 'Ruler']
  }

  if (stateData.hasOwnProperty('additional-data-layers')) {
    // TODO process
  } else {
    stateData['additional-data-layers'] = []
    generateMockADL()
  }

  // working under the assumption that all genome groups with contain the same additional data layers,
  // we can query the first genome group for specific ADL and go from there
  buildGroupLayersTable('Genome')
  if (stateData['additional-data-layers'][0]['ruler']) {
    buildGroupLayersTable('Ruler')
    // don't increase group size for ruler since it requires less space
  }
  if (stateData['additional-data-layers'][0]['coverage']) {
    buildGroupLayersTable('Coverage')
    stateData['group-layer-order'].push('Coverage')
    maxGroupSize += 1 // increase group size if coverage layer exists
  }

  if (stateData['additional-data-layers'][0]['gcContent']) {
    buildGroupLayersTable('GC_Content')
    stateData['group-layer-order'].push('GC_Content')
    maxGroupSize += 1 // increase group size if GC layer exists
  }

  if (stateData.hasOwnProperty('genome-order-method')) {
    stateData['genome-order-method'].forEach(orderMethod => {
      $('#genome_order_select').append((new Option(orderMethod["name"], orderMethod["name"]))) // set display + value of new select option.
    })
  } else {
    generateMockGenomeOrder()
    stateData['genome-order-method'].forEach(orderMethod => {
      $('#genome_order_select').append((new Option(orderMethod["name"], orderMethod["name"]))) // set display + value of new select option.
    })
  }

  if (stateData.hasOwnProperty('display')) {
    // TODO process
  } else {
    stateData['display'] = {}
    stateData['display']['additionalDataLayers'] = {}
  }
  if (stateData['display'].hasOwnProperty('bookmarks')) {
    stateData['display']['bookmarks'].map(bookmark => {
      $('#bookmarks-select').append((new Option(bookmark['name'], [bookmark["start"], bookmark['stop']])))
    })
    respondToBookmarkSelect()
  } else {
    calculateMaxGenomeLength() // remove later, only here to set max length for bookmark
    stateData['display']['bookmarks'] = [ // gen mock data
      {
        name: 'entire seq',
        start: '0',
        stop: genomeMax,
        description: 'a mighty fine placeholder'
      },
      {
        name: 'shindig',
        start: '5000',
        stop: '9000',
        description: 'a beautiful placeholder'
      },
      {
        name: 'fiesta',
        start: '15000',
        stop: '19000',
        description: 'an adequate placeholder'
      },
      {
        name: 'party',
        start: '25000',
        stop: '29000',
        description: 'the very best placeholder'
      },
    ]
    stateData['display']['bookmarks'].map(bookmark => {
      $('#bookmarks-select').append((new Option(bookmark['name'], [bookmark["start"], bookmark['stop']])))
    })
    // set listener for user bookmark selection
    respondToBookmarkSelect()
  }

  function generateMockADL() {
    for (let i = 0; i < genomeData.genomes.length; i++) { // generate mock additional data layer content
      let gcContent = []
      let coverage = []
      for (let j = 0; j < genomeMax; j++) {
        gcContent.push(Math.floor(Math.random() * 45))
        coverage.push(Math.floor(Math.random() * 45))
      }
      let genomeLabel = Object.keys(genomeData.genomes[i][1]['contigs']['info'])[0];
      let additionalDataObject = {
        'genome': genomeLabel,
        'coverage': coverage,
        'coverage-color': 'blue',
        'gcContent': gcContent,
        'gcContent-color': 'purple',
        'ruler': true // TODO: store any genome-specific scale data here
      }
      stateData['additional-data-layers'].push(additionalDataObject)
    }
  }
  function generateMockGenomeOrder() {
    stateData['genome-order-method'] = [{
      'name': 'cats',
      'ordering': 'some order'
    }, {
      'name': 'dogs',
      'ordering': 'some other order'
    }, {
      'name': 'birds',
      'ordering': 'beaks to tails'
    }]
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

  $('#brush_start').val(0);
  $('#brush_end').val(Math.floor(canvas.getWidth()));

  $('#tooltip-body').hide() // set initual tooltip hide value
  $('#show_genome_labels_box').attr("checked", showLabels);
  $('#show_gene_labels_box').attr("checked", showGeneLabels);
  $('#show_dynamic_scale_box').attr("checked", dynamicScaleInterval);

  // can either set arrow click listener on the canvas to check for all arrows, or when arrow is created.

  // disable group selection
  canvas.on('selection:created', (e) => {
    if (e.target.type === 'activeSelection') {
      canvas.discardActiveObject();
    }
  })

  canvas.on('mouse:over', (event) => {
    if (event.target && event.target.id === 'arrow') {
      showToolTip(event)
    }
  })

  canvas.on('mouse:out', (event) => {
    $('#tooltip-body').html('').hide()
  })

  if (showGeneLabels && arrowStyle != 3) {
    marginTop = 60;
    spacing = 200; // TODO maybe we refactor this out into a setSpacing() method for clarity?
    $("#genome_spacing").val(spacing);
  }

  function showToolTip(event) {
    $('#tooltip-body').show().append(`
            <p></p>
            <style type="text/css">
              .tftable {font-size:12px;color:#333333;width:100%;border-width: 1px;border-color: #729ea5;border-collapse: collapse;}
              .tftable th {font-size:12px;background-color:#acc8cc;border-width: 1px;padding: 8px;border-style: solid;border-color: #729ea5;text-align:left;}
              .tftable tr {background-color:#d4e3e5;}
              .tftable td {font-size:12px;border-width: 1px;padding: 8px;border-style: solid;border-color: #729ea5;}
              .tftable tr:hover {background-color:#ffffff;}
            </style>
      
            <table class="tftable" border="1">
              <tr><th>Data</th><th>Value</th></tr>
              <tr><td>Split</td><td>${event.target.gene.contig}</td></tr>
              <tr><td>Start in Contig</td><td>${event.target.gene.start}</td></tr>
              <tr><td>Length</td><td>${event.target.gene.stop - event.target.gene.start}</td></tr>
              <tr><td>Gene Callers ID</td><td>${event.target.geneID}</td></tr>
              <tr><td>Gene Cluster</td><td>${genomeData.gene_associations["anvio-pangenome"] ? genomeData.gene_associations["anvio-pangenome"]["genome-and-gene-names-to-gene-clusters"][event.target.genomeID][event.target.geneID] : "None"}</td></tr>
            </table>
            <button>some action</button>
            <button>some other action</button>
            `).css({ 'position': 'absolute', 'left': event.e.clientX, 'top': event.e.clientY })
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
  draw();

  // panning
  canvas.on('mouse:down', function (opt) {
    var evt = opt.e;
    if (evt.shiftKey === true) {
      this.isDragging = true;
      this.selection = false;
      this.lastPosX = evt.clientX;
    }
    this.shades = true;
    if (opt.target && opt.target.groupID) this.prev = opt.target.left;
  });
  canvas.on('mouse:move', function (opt) {
    if (this.isDragging) {
      var e = opt.e;
      var vpt = this.viewportTransform;
      vpt[4] += e.clientX - this.lastPosX;
      bindViewportToWindow();
      this.requestRenderAll();
      this.lastPosX = e.clientX;

      let [l, r] = getFracForVPT();
      if (l < renderWindow[0] || r > renderWindow[1]) {
        updateRenderWindow();
        draw();
      }
    }
  });
  canvas.on('mouse:up', function (opt) {
    this.setViewportTransform(this.viewportTransform);
    if (this.isDragging) updateScalePos();
    this.isDragging = false;
    this.selection = true;
    if (!this.shades) {
      // slid a genome
      this.shades = true;
      drawTestShades();
      bindViewportToWindow();
      updateScalePos(); // adjust scale box to new sequence breadth
      updateRenderWindow();

      redrawGenome(opt.target.groupID);
    }
  });
  canvas.on('object:moving', function (opt) {
    var gid = opt.target ? opt.target.groupID : null;
    if (gid == null) return;

    if (opt.target.id == 'genomeLine' || (opt.target.id == 'arrow' && arrowStyle == 3)) canvas.sendBackwards(opt.target);
    if (this.shades) {
      clearShades();
      this.shades = false;
    }

    let objs = canvas.getObjects().filter(obj => obj.groupID == gid);

    var delta = opt.target.left - this.prev;
    canvas.getObjects().filter(obj => obj.groupID == gid).forEach(o => {
      if (o !== opt.target) o.left += delta;
    });
    xDisps[gid] += delta;

    this.setViewportTransform(this.viewportTransform);
    setPercentScale();
    this.prev = opt.target.left;
  });
  canvas.on('mouse:wheel', function (opt) {
    opt.e.preventDefault();
    opt.e.stopPropagation();

    var delta = opt.e.deltaY;
    let tmp = scaleFactor * (0.999 ** delta);
    let diff = tmp - scaleFactor;
    let [start, end] = [parseInt($('#brush_start').val()), parseInt($('#brush_end').val())];
    let [newStart, newEnd] = [Math.floor(start - diff * genomeMax), Math.floor(end + diff * genomeMax)];
    if (newStart < 0) newStart = 0;
    if (newEnd > genomeMax) newEnd = genomeMax;
    if (newEnd - newStart < 50) return;

    brush.extent([newStart, newEnd]);
    brush(d3.select(".brush").transition()); // if zoom is slow or choppy, try removing .transition()
    brush.event(d3.select(".brush"));
    $('#brush_start').val(newStart);
    $('#brush_end').val(newEnd);
  });

  $('#alignClusterInput').on('keydown', function (e) {
    if (e.keyCode == 13) { // 13 = enter key
      alignToCluster($(this).val());
      $(this).blur();
    }
  });
  $('#panClusterInput').on('keydown', function (e) {
    if (e.keyCode == 13) { // 13 = enter key
      viewCluster($(this).val());
      $(this).blur();
    }
  });
  document.body.addEventListener("keydown", function (ev) {
    if (ev.which == 83 && ev.target.nodeName !== 'TEXTAREA' && ev.target.nodeName !== 'INPUT') { // S = 83
      toggleSettingsPanel();
    }
  });
  $('#genome_spacing').on('keydown', function (e) {
    if (e.keyCode == 13) { // 13 = enter key
      setGenomeSpacing($(this).val());
      $(this).blur();
    }
  });
  $('#gene_label').on('keydown', function (e) {
    if (e.keyCode == 13) { // 13 = enter key
      setGeneLabelSize($(this).val());
      $(this).blur();
    }
  });
  $('#genome_label').on('keydown', function (e) {
    if (e.keyCode == 13) { // 13 = enter key
      setGenomeLabelSize($(this).val());
      $(this).blur();
    }
  });
  $('#genome_scale_interval').on('keydown', function (e) {
    if (e.keyCode == 13) { // 13 = enter key
      setScaleInterval($(this).val());
      $(this).blur();
    }
  });
  $('#gene_color_order').on('change', function () {
    color_db = $(this).val();
    generateColorTable(null, color_db); // TODO: include highlight_genes, fn_colors etc from state
    draw();
    $(this).blur();
  });
  $('#arrow_style').on('change', function () {
    arrowStyle = parseInt($(this).val());
    draw();
    $(this).blur();
  });
  $('#gene_text_pos').on('change', function () {
    geneLabelPos = $(this).val();
    if (!(geneLabelPos == "inside" && arrowStyle != 3)) draw();
    $(this).blur();
  });
  $('#show_genome_labels_box').on('change', function () {
    showLabels = !showLabels;
    xDisplacement = showLabels ? 120 : 0;
    alignToGC = null;
    draw();
  });
  $('#show_gene_labels_box').on('change', function () {
    showGeneLabels = !showGeneLabels;
    draw();
  });
  $('#show_dynamic_scale_box').on('change', function () {
    dynamicScaleInterval = !dynamicScaleInterval;
  });
  $('#adl_pts_per_layer').on('change', function () {
    setPtsPerADL($(this).val());
    $(this).blur();
  });
  $('#brush_start, #brush_end').keydown(function (ev) {
    if (ev.which == 13) { // enter key
      let [start, end] = percentScale ? [parseFloat($('#brush_start').val()), parseFloat($('#brush_end').val())]
        : [parseInt($('#brush_start').val()), parseInt($('#brush_end').val())];
      let endBound = percentScale ? 1 : genomeMax;

      if (isNaN(start) || isNaN(end) || start < 0 || start > endBound || end < 0 || end > endBound) {
        alert(`Invalid value, value needs to be in range 0-${endBound}.`);
        return;
      }

      if (start >= end) {
        alert('Starting value cannot be greater or equal to the ending value.');
        return;
      }

      brush.extent([start, end]);
      brush(d3.select(".brush").transition());
      brush.event(d3.select(".brush").transition());
    }
  });
}