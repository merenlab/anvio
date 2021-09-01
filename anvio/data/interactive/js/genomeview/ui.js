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
 * File Overview : This file contains functions related to building + updating UI elements and responding to user interaction with those elements. 
 * As a general rule, processes that invoke jQuery should probably live here. 
 */

/*
 *  set event listeners for DOM elements, user input, default jquery values
 */
function setEventListeners(){
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
        // draw();
        drawer.draw()
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
      drawer.redrawSingleGenome(opt.target.groupID);
    }
  });
  canvas.on('object:moving', function (opt) {
    var gid = opt.target ? opt.target.groupID : null;
    if (gid == null) return;

    if (opt.target.id == 'genomeLine' || (opt.target.id == 'arrow' && arrowStyle == 3)) canvas.sendBackwards(opt.target);
    if (this.shades) {
      drawer.clearShades();
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
      drawer.alignToCluster($(this).val());
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
      drawer.setGenomeSpacing($(this).val());
      $(this).blur();
    }
  });
  $('#gene_label').on('keydown', function (e) {
    if (e.keyCode == 13) { // 13 = enter key
      drawer.setGeneLabelSize($(this).val());
      $(this).blur();
    }
  });
  $('#genome_label').on('keydown', function (e) {
    if (e.keyCode == 13) { // 13 = enter key
      drawer.setGenomeLabelSize($(this).val());
      $(this).blur();
    }
  });
  $('#genome_scale_interval').on('keydown', function (e) {
    if (e.keyCode == 13) { // 13 = enter key
      drawer.setScaleInterval($(this).val());
      $(this).blur();
    }
  });
  $('#gene_color_order').on('change', function () {
    color_db = $(this).val();
    generateColorTable(null, color_db); // TODO: include highlight_genes, fn_colors etc from state
    // draw();
    drawer.draw()
    $(this).blur();
  });
  $('#arrow_style').on('change', function () {
    arrowStyle = parseInt($(this).val());
    // draw();
    drawer.draw()
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
    // draw();
    drawer.draw()
  });
  $('#show_gene_labels_box').on('change', function () {
    showGeneLabels = !showGeneLabels;
    // draw();
    drawer.draw()
  });
  $('#show_dynamic_scale_box').on('change', function () {
    dynamicScaleInterval = !dynamicScaleInterval;
  });
  $('#adl_pts_per_layer').on('change', function () {
    drawer.setPtsPerADL($(this).val());
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
}

function showToolTip(event){
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
  `).css({'position' : 'absolute', 'left' : event.e.clientX, 'top' : event.e.clientY })
}

function toggleSettingsPanel() {
  $('#settings-panel').toggle();

  if ($('#settings-panel').is(':visible')) {
    $('#toggle-panel-settings').addClass('toggle-panel-settings-pos');
    $('#toggle-panel-settings-inner').html('&#9658;');
  } else {
    $('#toggle-panel-settings').removeClass('toggle-panel-settings-pos');
    $('#toggle-panel-settings-inner').html('&#9664;');
  }
}

function buildGenomesTable(genomes, order){
  genomes.map(genome => {
    var height = '50';
    var margin = '15';
    var template = '<tr id={genomeLabel}>' +
                  '<td><img src="images/drag.gif" class="drag-icon" id={genomeLabel} /></td>' +
                  '<td> {genomeLabel} </td>' +
                  '<td>n/a</td>' +
                  '<td>n/a</td>' +
                  '<td>n/a</td>' +
                  '<td><input class="input-height" type="text" size="3" id="height{id}" value="{height}"></input></td>' +
                  '<td class="column-margin"><input class="input-margin" type="text" size="3" id="margin{id}" value="{margin}"></input></td>' +
                  '<td>n/a</td>' +
                  '<td>n/a</td>' +
                  '<td><input type="checkbox" class="layer_selectors"></input></td>' +
                  '</tr>';
    let genomeLabel= Object.keys(genome[1]['contigs']['info']);

    template = template.replace(new RegExp('{height}', 'g'), height)
                        .replace(new RegExp('{margin}', 'g'), margin)
                        .replace(new RegExp('{genomeLabel}', 'g'), genomeLabel);

    $('#tbody_genomes').append(template);
  })
  
  $("#tbody_genomes").sortable({helper: fixHelperModified, handle: '.drag-icon', items: "> tr"}).disableSelection();

  $("#tbody_genomes").on("sortupdate", (event, ui) => {
    changeGenomeOrder($("#tbody_genomes").sortable('toArray'))
  })
}
  
function buildGroupLayersTable(layerLabel){
  var height = '50';
  var margin = '25';
  var template = '<tr id={layerLabel}>' +
                  '<td><img src="images/drag.gif" class="drag-icon" id={layerLabel} /></td>' +
                  '<td> {layerLabel} </td>' +
                  '<td><div id="{layerLabel}_color" style="margin-left: 5px;" class="colorpicker" style="background-color: #FFFFFF" color="#FFFFFF"></div></td>' +
                  '<td>n/a</td>' +
                  '<td>n/a</td>' +
                  '<td><input type="checkbox" class="additional_selectors" id={layerLabel}-show onclick="toggleAdditionalDataLayer(event)" checked=true></input></td>' +
                  '</tr>';
  template = template.replace(new RegExp('{height}', 'g'), height)
                      .replace(new RegExp('{margin}', 'g'), margin)
                      .replace(new RegExp('{layerLabel}', 'g'), layerLabel);
  $('#tbody_additionalDataLayers').append(template);
  $("#tbody_additionalDataLayers").sortable({helper: fixHelperModified, handle: '.drag-icon', items: "> tr"}).disableSelection();

  $("#tbody_additionalDataLayers").on("sortupdate", (event, ui) => {
    changeGroupLayersOrder($("#tbody_additionalDataLayers").sortable('toArray'))
  })
}
  
function toggleAdditionalDataLayer(e){
  let layer = e.target.id.split('-')[0]

  if(e.target.checked){
    stateData['display']['additionalDataLayers'][layer] = true
    maxGroupSize += 1
  } else {
    stateData['display']['additionalDataLayers'][layer] = false
    maxGroupSize -= 1 // decrease group height if hiding the layer
  }
  // draw()
  drawer.draw()
}

/*
 *  capture user bookmark values and store obj in state
 */
function createBookmark(){
  if(!$('#create_bookmark_input').val()){
    alert('please provide a name for your bookmark :)')
    return
  }
  try {
    stateData['display']['bookmarks'].push(
      {
        name : $('#create_bookmark_input').val(),
        start : $('#brush_start').val(),
        stop : $('#brush_end').val(),
        description : $('#create_bookmark_description').val(),
      }
    )
    alert('bookmark successfully created :)')
  } catch (error) {
    alert(`anvi'o was unable to save your bookmark because of an error ${error} :/`)
    // throw to error landing page?
  }
}
/*
 *  update sequence position, bookmark description upon user select from dropdown
 */
function respondToBookmarkSelect(){
  $('#bookmarks-select').change(function(e){
    let [start, stop] = [$(this).val().split(',')[0], $(this).val().split(',')[1] ]
    $('#brush_start').val(start);
    $('#brush_end').val(stop);
    brush.extent([start, stop]);
        brush(d3.select(".brush").transition());
        brush.event(d3.select(".brush").transition());
    let selectedBookmark = stateData['display']['bookmarks'].find(bookmark => bookmark.start == start && bookmark.stop == stop)
    $('#bookmark-description').text(selectedBookmark['description'])
  })
}

/*
 *  respond to ui, redraw with updated group layer order
 */
function changeGroupLayersOrder(updatedOrder){
  stateData['group-layer-order'] = updatedOrder
  // draw()
  drawer.draw()
}

/*
 *  respond to ui, redraw with updated genome group order
 */
function changeGenomeOrder(updatedOrder){
  let newGenomeOrder = []
  updatedOrder.map(label => {
    genomeData.genomes.map(genome => {
        if(label == Object.keys(genome[1]['contigs']['info'])[0]){ // matching label text to first contig name of each genome
          newGenomeOrder.push(genome)
        }
    })
  })
  genomeData.genomes = newGenomeOrder
  // draw()
  drawer.draw()
}


/*
 *  [TO BE ADDED TO 'regular' utils.js]
 *  Generates functional annotation color table for a given color palette.
 *
 *  @param fn_colors :       dict matching each category to a hex color code to override defaults
 *  @param fn_type :         string indicating function category type: currently one of "COG_CATEGORY", "KEGG_CATEGORY", "Source"
 *  @param highlight_genes : array of format [{genomeID: 'g01', geneID: 3, color: '#FF0000'}, ...] to override other coloring for specific genes
 *  @param filter_to_split : if true, filters categories to only those shown in the split
 */
function generateColorTable(fn_colors, fn_type, highlight_genes=null, filter_to_split=true) {
  // TODO: consider call_type? see inspectionalutils.js for how this was dealt with earlier

  let db = getColorDefaults(fn_type ? fn_type : 'Source');
  if(db == null) return;
  // Override default values with any values supplied to fn_colors
  if(fn_colors) db = Object.keys(db).map(cag => Object.keys(fn_colors).includes(cag) ? fn_colors[cag] : db[cag]);

  if(filter_to_split && fn_type != 'Source') {
    let save = [];
    for(genome of genomeData.genomes) {
      let geneFuns = Object.values(genome[1].genes.functions);
      for(funs of geneFuns) {
        let cag = getCagForType(funs, fn_type);
        if(cag && !save.includes(cag)) save.push(cag);
      }
    }
    Object.keys(db).forEach((cag, color) => {
      if(!save.includes(cag)) delete db[cag];
    });
  }

  $('#tbody_function_colors').empty();
  Object.keys(db).forEach(category => appendColorRow(getCagName(category, fn_type), category, db[category]) );

  $('.colorpicker').colpick({
      layout: 'hex',
      submit: 0,
      colorScheme: 'light',
      onChange: function(hsb, hex, rgb, el, bySetColor) {
          $(el).css('background-color', '#' + hex);
          $(el).attr('color', '#' + hex);
          // TODO: save new color once state is implemented
          //state[$('#gene_color_order').val().toLowerCase() + '-colors'][el.id.substring(7)] = '#' + hex;
          if (!bySetColor) $(el).val(hex);
      }
  }).keyup(function() {
      $(this).colpickSetColor(this.value);
  });

  if(highlight_genes) {
    let genomes = Object.entries(genomeData.genomes).map(g => g[1][0]);
    for(entry of highlight_genes) {
      let genomeID = entry['genomeID'];
      let geneID = entry['geneID'];
      let color = entry['color'];

      if(!genomes.includes(genomeID)) continue;

      let ind = genomeData.genomes.findIndex(g => g[0] == genomeID);
      let genes = Object.keys(genomeData.genomes[ind][1].genes.gene_calls);
      if(!(geneID in genes)) continue;

      let label = 'Genome: ' + genomeID + ', Gene: ' + geneID;
      appendColorRow(label, genomeID + '-' + geneID, color, prepend=true);
    }
    $('colorpicker').colpick({
        layout: 'hex',
        submit: 0,
        colorScheme: 'light',
        onChange: function(hsb, hex, rgb, el, bySetColor) {
            $(el).css('background-color', '#' + hex);
            $(el).attr('color', '#' + hex);
            //state['highlight-genes'][el.id.substring(7)] = '#' + hex;
            if (!bySetColor) $(el).val(hex);
        }
    }).keyup(function() {
        $(this).colpickSetColor(this.value);
    });
  }
}

/*
 *  [TO BE ADDED TO 'regular' utils.js]
 */
function appendColorRow(label, cag, color, prepend=false) {
  let code = cag.split(' ').join('_').split('(').join('_').split(')').join('_').split(':').join('_');
  var tbody_content =
   '<tr id="picker_row_' + code + '"> \
      <td></td> \
      <td> \
        <div id="picker_' + code + '" class="colorpicker" color="' + color + '" background-color="' + color + '" style="background-color: ' + color + '; margin-right:16px; margin-left:16px"></div> \
      </td> \
      <td>' + label + '</td> \
    </tr>';

  if(prepend) {
    $('#tbody_function_colors').prepend(tbody_content);
  } else {
    $('#tbody_function_colors').append(tbody_content);
  }
}

/*
 * [TO BE ADDED to genomeview/UI.js OR 'regular' utils.js]
 * Resets function color table to the default set
 *
 * @param fn_colors: If set, resets state to this dictionary instead of the defaults.
 */
function resetFunctionColors(fn_colors=null) {
  // TODO: this should reset color dictionaries in state and then redraw table, once state is implemented
  if($('#gene_color_order') == null) return;
  generateColorTable(fn_colors, color_db);
}
