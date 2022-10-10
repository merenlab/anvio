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
  canvas.on('mouse:over', (event) => {
    if (event.target && event.target.id === 'arrow') {
      showToolTip(event)
    }
  })

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

      let [l, r] = percentScale ? getFracForVPT() : getNTRangeForVPT();
      if (l < renderWindow[0] || r > renderWindow[1]) {
        updateRenderWindow();
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
      // slide a genome
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

    if (opt.target.id == 'genomeLine' || (opt.target.id == 'arrow' && settings['display']?.['arrowStyle'] == 3)) canvas.sendBackwards(opt.target);
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
  canvas.on('mouse:over', function(event) {
    // if(event.target.class == 'ruler'){
    //   console.log(event.target)
    // }
  })

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
  document.body.addEventListener("keydown", function (ev) {
    if (ev.which == 77 && ev.target.nodeName !== 'TEXTAREA' && ev.target.nodeName !== 'INPUT') { // M = 77
      toggleRightPanel('#mouseover-panel')
    }
  });
  document.body.addEventListener("keydown", function (ev) {
    if (ev.which == 81 && ev.target.nodeName !== 'TEXTAREA' && ev.target.nodeName !== 'INPUT') { // Q = 81
      toggleRightPanel('#query-panel')
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
    generateColorTable(settings.display.colors.genes.annotations[color_db], color_db); // TODO: include highlight_genes, fn_colors etc from state

    if(link_gene_label_color_source){
      $('#gene_label_source').val($(this).val())
    }

    drawer.draw()
    $(this).blur();
  });
  $('#arrow_style').on('change', function () {
    // arrowStyle = parseInt($(this).val());
    settings['display']['arrow-style'] = parseInt($(this).val());
    drawer.draw()
    $(this).blur();
  });
  $('#gene_text_pos').on('change', function () {
    geneLabelPos = $(this).val();
    if (!(geneLabelPos == "inside" && settings['display']['arrow-style'] != 3)) drawer.draw();
    $(this).blur();
  });
  $('#gene_text_angle').on('change', function () {
    geneLabelAngle = $(this).val();
    if (geneLabelPos != "inside") drawer.draw();
    $(this).blur();
  });
  $('#show_genome_labels_box').on('change', function () {
    showLabels = !showLabels;
    xDisplacement = showLabels ? 120 : 0;
    alignToGC = null;
    drawer.draw()
  });
  $('#gene_label_source').on('change', function(){
    if(link_gene_label_color_source){
      color_db = $(this).val();
      $('#gene_color_order').val($(this).val())
      generateColorTable(settings.display.colors.genes.annotations[color_db], color_db)
    }
    drawer.draw()
  })
  $('#show_gene_labels_box').on('change', function () {
    showGeneLabels = !showGeneLabels;
    drawer.draw()
  });
  $('#link_gene_label_color_source').on('change', function(){
    link_gene_label_color_source = !link_gene_label_color_source;
    drawer.draw()
  })
  $('#show_only_cags_in_window').on('change', function () {
    filter_gene_colors_to_window = !filter_gene_colors_to_window;
    generateColorTable(settings.display.colors.genes.annotations[color_db], color_db);
  });
  $('#thresh_count').on('keydown', function (e) {
    if (e.keyCode == 13) { // 13 = enter key
      filterColorTable($(this).val());
      $(this).blur();
    }
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
  $('#batch_colorpicker').colpick({
      layout: 'hex',
      submit: 0,
      colorScheme: 'light',
      onChange: function(hsb, hex, rgb, el, bySetColor) {
          $(el).css('background-color', '#' + hex);
          $(el).attr('color', '#' + hex);
          settings['display']['colors']['Batch'] = '#' + hex;
          if (!bySetColor) $(el).val(hex);
      }
  }).keyup(function() {
      $(this).colpickSetColor(this.value);
  });
  $('#genome_label_color').colpick({
      layout: 'hex',
      submit: 0,
      colorScheme: 'light',
      onChange: function(hsb, hex, rgb, el, bySetColor) {
          $(el).css('background-color', '#' + hex);
          $(el).attr('color', '#' + hex);
          if (!bySetColor) $(el).val(hex);
      }
  }).keyup(function() {
      $(this).colpickSetColor(this.value);
  });
  $('#gene_label_color').colpick({
      layout: 'hex',
      submit: 0,
      colorScheme: 'light',
      onChange: function(hsb, hex, rgb, el, bySetColor) {
          $(el).css('background-color', '#' + hex);
          $(el).attr('color', '#' + hex);
          if (!bySetColor) $(el).val(hex);
      }
  }).keyup(function() {
      $(this).colpickSetColor(this.value);
  });

  $('#brush_start').val(0);
  $('#brush_end').val(Math.floor(canvas.getWidth()));

  $('#deepdive-tooltip-body').hide() // set initual tooltip hide value
  $('#tooltip-body').hide()
  $('#show_genome_labels_box').attr("checked", showLabels);
  $('#show_gene_labels_box').attr("checked", showGeneLabels);
  $('#show_dynamic_scale_box').attr("checked", dynamicScaleInterval);

  $("#tabular-modal-body").on('shown.bs.modal', function(){
    showTabularModal()
  });
  $("#tabular-modal-body").on('hide.bs.modal', function(){
    $('#tabular-modal-nav-tabs').html('')
    $('#modal-tab-content').html('')
  });

  // can either set arrow click listener on the canvas to check for all arrows, or when arrow is created.

  // drag box selection
  canvas.on('selection:created', (e) => {
    let selected_genes = e.selected.filter(obj => obj.id == 'arrow');

    if(selected_genes.length > 1) {
      showLassoMenu(selected_genes, e.e.clientX, e.e.clientY);
    }

    // disable group selection
    if (e.target.type === 'activeSelection') {
      canvas.discardActiveObject();
    }
  })

  canvas.on('mouse:down', (event) => {
    if(event.target && event.target.id === 'arrow'){
      showDeepDiveToolTip(event)
    }
    $('#lasso-modal-body').modal('hide')
  })
}
function showDeepDiveToolTip(event){
  $('.canvas-container').dblclick(function(){
    if($("#deepdive-tooltip-body").is(":visible")){
      $("#deepdive-tooltip-body").hide()
    }
  })

  $('#deepdive-modal-tab-content').html('')
  let totalMetadataString = String()
  let totalAnnotationsString = String()
  let metadataLabel = String()
  let button = `<button type='button' id='metadata-query' class='btn btn-default btn-sm'>Query Sequence for matches</button>`

  if(settings['display']?.['metadata']){
    let geneMetadata = settings['display']['metadata'].filter(metadata => metadata.genome == event.target.genomeID && metadata.gene == event.target.geneID )
    const createMetadataContent = () => {
      geneMetadata.map(metadata => {
        metadataLabel = metadata.label

        totalMetadataString += `
        <tr>
        <td>${metadata.label}</td>
        <td>${metadata.type == 'tag'? button : 'n/a'}</td>
        </tr>
        `
      })
    }
    createMetadataContent()
  } else { // create metadata array on first tooltip load if none exists
    settings['display']['metadata'] = []
  }

  if(event.target.functions){
    Object.entries(event.target.functions).map(func => {
      totalAnnotationsString += `
      <tr>
      <td>${func[0]}</td>
      <td>${func[1] ? func[1][0] : 'n/a'}</td>
      <td>${func?.[1]?.[1]}</td>
      </tr>
      `
    })
  }

  $('#deepdive-modal-body').modal('show')
  $('#deepdive-modal-tab-content').append(`
  <span class="popover-close-button" onclick="$(this).closest(\'.popover\').popover(\'hide\');"></span>
  <h2>Gene Call</h2>
  <table class="table table-striped" style="width: 100%; text-align: center;">
    <thead>
      <th>ID</th>
      <th>Source</th>
      <th>Length</th>
      <th>Direction</th>
      <th>Start</th>
      <th>Stop</th>
      <th>Call type</th>
    </thead>
    <tbody>
      <tr>
        <td>${event.target.geneID}</td>
        <td>${event.target.gene?.source}</td>
        <td>${event.target.gene.stop - event.target.gene.start}</td>
        <td>${event.target.gene.direction}</td>
        <td>${event.target.gene.start}</td>
        <td>${event.target.gene.stop}</td>
        <td>${event.target.gene?.call_type}</td>
      </tr>
    </tbody>
  </table>;

  <button type="button" id="gene-dna-sequence-button"class="btn btn-default btn-sm" >DNA</button>
  <button type="button" id="gene-aa-sequence-button" class="btn btn-default btn-sm" >AA</button>
  <button type="button" id="gene-blastn-at-nr-button"class="btn btn-default btn-sm" >blastn @ nr</button>
  <button type="button" id="gene-blastn-at-refseq-button" class="btn btn-default btn-sm">blastn @ refseq_genomic</button>
  <button type="button" id="gene-blastx-at-nr-button"class="btn btn-default btn-sm" >blastx @ nr</button>
  <button type="button" id="gene-blastx-at-refseq-button" class="btn btn-default btn-sm">blastx @ refseq_genomic</button>
  <br>

  <h2>Visible Geen Range</h2>
  <label for='hide-range-low'>Downstream Visible Gene Range</label>
  <input type='number' id='hide-range-low'  min='1' step='1'/>
  <label for='hide-range-high'>Upstream Visible Gene Range</label>
  <input type='number' id='hide-range-high' min='1' step='1'/>
  <br>
  <button id='gene-visibility-range-set' class="btn btn-default btn-sm" >Set Visible Range</button>
  <br>

  <h2>color</h2>
  <div id="picker_tooltip" class="colorpicker" color="#808080" background-color="#808080" style="background-color: #808080; margin-right:16px; margin-left:16px"></div>
  <p>set gene arrow color</p>
  <br>
  <h2>metadata</h2>
  <div id='metadata-container'>
      <input    id='metadata-gene-label' type='text' placeholder='metadata tag'>
      <button   id='metadata-gene-label-add' type='button' class="btn btn-default btn-sm">add tag</button>
      <br>
      <textarea id='metadata-gene-description' placeholder='metadata description'></textarea>
      <button   id='metadata-gene-description-add' type='button' class="btn btn-default btn-sm">add description</button>
      <br>
      <table class="table table-striped">
        <thead><th>metadata</th><th>action</th></thead>
        <tbody id="metadata-body">
         ${totalMetadataString}
        </tbody>
      </table>
  </div>

  <h2>Annotations</h2>;
  <table class="table table-striped">;
    <thead>
      <th>Source</th>;
      <th>Accession</th>;
      <th>Annotation</th>
    </thead>;
    <tbody>;
      ${totalAnnotationsString}
    </tbody>
  </table>;
  `)

  $('#metadata-query').on('click', function(){
    drawer.queryMetadata(metadataLabel)
  })

  $('#gene-visibility-range-set').click(function(){
    setGeneVisibilityRange(event.target.geneID, event.target.genomeID)
  })

  $('#gene-dna-sequence-button').on('click', function(){
    show_sequence_modal('DNA Sequence', settings['genomeData']['genomes'].filter(genome => genome[0] == event.target.genomeID)[0][1]['genes']['dna'][event.target.geneID], event.target.geneID, event.target.genomeID)
  })

  $('#gene-aa-sequence-button').on('click', function(){
    show_sequence_modal('AA Sequence', settings['genomeData']['genomes'].filter(genome => genome[0] == event.target.genomeID)[0][1]['genes']['aa'][event.target.geneID], event.target.geneID, event.target.genomeID)
  })

  $('#gene-blastx-at-nr-button').on('click', function(){
    let sequence = settings['genomeData']['genomes'].filter(genome => genome[0] == event.target.genomeID)[0][1]['genes']['dna'][event.target.geneID]['sequence']
    let sequenceConcat = '>' + 'DNA_SEQUENCE' + '\n' + sequence
    fire_up_ncbi_blast(sequenceConcat, 'blastx', 'nr', 'gene')
  })
  $('#gene-blastn-at-nr-button').on('click', function(){
    let sequence = settings['genomeData']['genomes'].filter(genome => genome[0] == event.target.genomeID)[0][1]['genes']['dna'][event.target.geneID]['sequence']
    let sequenceConcat = '>' + 'DNA_SEQUENCE' + '\n' + sequence
    fire_up_ncbi_blast(sequenceConcat, 'blastn', 'nr', 'gene')
  })

  $('#gene-blastx-at-refseq-button').on('click', function(){
    let sequence = settings['genomeData']['genomes'].filter(genome => genome[0] == event.target.genomeID)[0][1]['genes']['dna'][event.target.geneID]['sequence']
    let sequenceConcat = '>' + 'DNA_SEQUENCE' + '\n' + sequence
    fire_up_ncbi_blast(sequenceConcat, 'blastx', 'refseq_genomic', 'gene')
  })
  $('#gene-blastn-at-refseq-button').on('click', function(){
    let sequence = settings['genomeData']['genomes'].filter(genome => genome[0] == event.target.genomeID)[0][1]['genes']['dna'][event.target.geneID]['sequence']
    let sequenceConcat = '>' + 'DNA_SEQUENCE' + '\n' + sequence
    fire_up_ncbi_blast(sequenceConcat, 'blastn', 'refseq_genomic', 'gene')
  })
  // TODO consider metadata option to include 'author' field
  $('#metadata-gene-label-add').on('click', function(){
    let metadataObj = {
      label  : $('#metadata-gene-label').val(),
      genome : event.target.genomeID,
      gene   : event.target.geneID,
      type   : 'tag'
    }
    settings['display']['metadata'].push(metadataObj)
    $('#metadata-gene-label').val('')
    $('#metadata-body').append(`
      <tr>
        <td>${metadataObj.label}</td>
        <td>${button}</td>
      </tr>
    `)
    $('#metadata-query').on('click', function(){ // re-trigger listener for new DOM buttons
      drawer.queryMetadata(metadataLabel)
    })
  })
  $('#metadata-gene-description-add').on('click', function(){
    let metadataObj = {
      label  : $('#metadata-gene-description').val(),
      genome : event.target.genomeID,
      gene   : event.target.geneID,
      type   : 'description'
    }
    settings['display']['metadata'].push(metadataObj)
    $('#metadata-gene-description').val('')
    $('#metadata-body').append(`
      <tr>
        <td>${metadataObj.label}</td>
        <td>n/a</td>
      </tr>
    `)
    $('#metadata-query').on('click', function(){
      drawer.queryMetadata(metadataLabel)
    })
  })
  $('#picker_tooltip').colpick({
    layout: 'hex',
    submit: 0,
    colorScheme: 'light',
    onChange: function(hsb, hex, rgb, el, bySetColor) {
        $(el).css('background-color', '#' + hex);
        $(el).attr('color', '#' + hex);
        if (!bySetColor) $(el).val(hex);
        let gene = canvas.getObjects().find(obj => obj.id == 'arrow' &&
                       event.target.genomeID == obj.genomeID &&
                       event.target.geneID == obj.geneID);
        gene.fill = `#${hex}`
        gene.dirty = true
        if(!settings['display']['colors']['genes']?.[event.target.genomeID]?.[event.target.geneID]){
          if(!settings['display']['colors']['genes'][event.target.genomeID]) {
            settings['display']['colors']['genes'][event.target.genomeID] = [];
          }
          settings['display']['colors']['genes'][event.target.genomeID][event.target.geneID] = `#${hex}`
        } else {
          settings['display']['colors']['genes'][event.target.genomeID][event.target.geneID] = `#${hex}`
        }
        canvas.renderAll()
    }
  }).keyup(function() {
      $(this).colpickSetColor(this.value);
  });
}

function showToolTip(event){
  $('#mouseover-panel-table-body').html('')
  $('#mouseover-panel-annotations-table-body').html('')
  let totalAnnotationsString = ''
  if(event.target.functions){
    Object.entries(event.target.functions).map(func => {
      totalAnnotationsString += `
      <tr>
        <td>${func[0]}</td>
        <td>${func[1] ? func[1][0] : 'n/a'}</td>
        <td>${func?.[1]?.[1]}</td>
      </tr>
      `
    })
  }
  $('#mouseover-panel-table-body').append(`
    <tr>
      <td> ${event.target.geneID}</td>
      <td> ${event.target.gene?.source}</td>
      <td> ${event.target.gene.stop - event.target.gene.start}</td>
      <td> ${event.target.gene.direction}</td>
      <td> ${event.target.gene.start}</td>
      <td> ${event.target.gene.stop}</td>
      <td> ${event.target.gene?.call_type}</td>
    </tr>
  `)
  $('#mouseover-panel-annotations-table-body').append(totalAnnotationsString)
}

function show_sequence_modal(title, gene, geneID, genomeID) {
 // remove previous modal window
  $('.modal-sequence').modal('hide');
  $('.modal-sequence').remove();

  let header, content
  if(title == 'DNA Sequence'){
    header = `>${geneID}|contig:${gene['contig']}|start:${gene['start']}|stop:${gene['stop']}|direction:${gene['direction']}|rev_compd:${gene['rev_compd']}|length:${gene['length']}`
    content = header + '\n' + gene['sequence']
  } else if (title == 'AA Sequence' ){
    // get header data from equivalent gene dna object, as gene aa object only contains sequence
    let dna_gene_obj = settings['genomeData']['genomes'].filter(g => g[0] == genomeID)[0][1]['genes']['dna'][geneID]
    header = `>${geneID}|contig:${dna_gene_obj['contig']}|start:${dna_gene_obj['start']}|stop:${dna_gene_obj['stop']}|direction:${dna_gene_obj['direction']}|rev_compd:${dna_gene_obj['rev_compd']}|length:${dna_gene_obj['length']}`
    content = header + '\n' + gene['sequence']
  }

  $('body').append('<div class="modal modal-sequence" style="z-index: 10000;"> \
      <div class="modal-dialog"> \
          <div class="modal-content"> \
              <div class="modal-header"> \
                  <button class="close" data-dismiss="modal" type="button"><span>&times;</span></button> \
                  <h4 class="modal-title">' + title + '</h4> \
              </div> \
              <div class="modal-body"> \
                      <textarea class="form-control" style="width: 100%; height: 100%; font-family: monospace;" rows="16" onclick="$(this).select();" readonly>' + content + '</textarea> \
              </div> \
              <div class="modal-footer"> \
                  <button class="btn btn-default" data-dismiss="modal" type="button">Close</button> \
              </div> \
          </div> \
      </div> \
  </div>');
  $('.modal-sequence').modal('show');
  $('.modal-sequence textarea').trigger('click');
}

function showTabularModal(){
  var arrows = canvas.getObjects().filter(obj => obj.id == 'arrow')
  var genomesObj = Object()
  var functionSourcesArr = Array()

  arrows.map(arrow => {
    if(!genomesObj[arrow['genomeID']]){
      genomesObj[arrow['genomeID']] = [arrow]
    } else {
      genomesObj[arrow['genomeID']].push(arrow)
    }
    if(arrow['functions']){
      Object.keys(arrow['functions']).forEach(source => {
        if(!functionSourcesArr.includes(source)){
          functionSourcesArr.push(source)
        }
      })
    }
  })

  handleSequenceSourceSelect = (target, type) => {
    if($(target).is(':checked')){
      Object.entries(genomesObj).map(genome => {
        $(`#${genome[0]}-tabular-modal-table-header-tr`).append(`<th id='th-source-${genome[0]}-${type}'>${type}</th>`)
        genome[1].forEach(gene => {
          if(type == 'dna'){
            $(`#${genome[0]}-table-row-${gene['geneID']}`).append(`<td id='${genome[0]}-${gene['geneID']}-sequence'>${gene['dnaSequence']}</td>`)
          } else if(type == 'aa'){
            $(`#${genome[0]}-table-row-${gene['geneID']}`).append(`<td id='${genome[0]}-${gene['geneID']}-sequence'>${gene['aaSequence']}</td>`)
          }
        })
      })
    } else {
      Object.entries(genomesObj).map(genome => {
        $(`#th-source-${genome[0]}-${type}`).remove()
        genome[1].forEach(gene => {
          $(`#${genome[0]}-${gene['geneID']}-sequence`).remove()
        })
      })

    }
  }

  handleAnnotationSourceSelect = (target) => {
    if($(target).is(":checked")){
      Object.entries(genomesObj).map(genome => {
        $(`#${genome[0]}-tabular-modal-table-header-tr`).append(`<th id='th-source-${genome[0]}-${target.value}'>${target.value}</th>`)
        genome[1].forEach(gene => {
          if(gene?.['functions']?.[target.value]){
            $(`#${genome[0]}-table-row-${gene['geneID']}`).append(`<td id='${genome[0]}-${gene['geneID']}-${target.value}'>${gene['functions'][target.value][1]}</td>`)
          } else {
            $(`#${genome[0]}-table-row-${gene['geneID']}`).append(`<td id='${genome[0]}-${gene['geneID']}-${target.value}'>n/a</td>`)
          }
        })
      })
    } else {
      Object.entries(genomesObj).map(genome => {
        $(`#th-source-${genome[0]}-${target.value}`).remove()
        genome[1].forEach(gene => {
          $(`#${genome[0]}-${gene['geneID']}-${target.value}`).remove()
        })
      })
    }
  }

  handleGeneShowHideToggle = (target) => {
    let [genomeID, geneID] = target.value.split('-')
    let arrow = canvas.getObjects().filter(obj => obj.id == 'arrow').find(arrow => arrow.geneID == geneID && arrow.genomeID == genomeID)
    if($(target).is(':checked')){
      if(genomeID in settings['display']['hidden']){
        settings['display']['hidden'][genomeID][geneID] = true
      } else {
        settings['display']['hidden'][genomeID] = {}
        settings['display']['hidden'][genomeID][geneID] = true
      }
      arrow.opacity = 0.1
    } else {
      delete settings['display']['hidden'][genomeID][geneID]
      arrow.opacity = 1
    }
    // re-render arrow objects with updated opacity values
    arrow.dirty = true
    canvas.renderAll()
  }

  $('#tabular-modal-annotation-checkboxes').empty()
  functionSourcesArr.map(s => {
    $('#tabular-modal-annotation-checkboxes').append(
      `
      <input type='checkbox' name=${s} value=${s} onclick='handleAnnotationSourceSelect(this)'/>
      <label for=${s}>${s}</label><br>
      `
    )
  })

  $('#tabular-modal-sequence-checkboxes').empty()
  $('#tabular-modal-sequence-checkboxes').append(
	`
	<label for="'tabular-dna-seq-checkbox">DNA Sequence</label>
	<input type='checkbox' id='tabular-dna-seq-checkbox' onclick="handleSequenceSourceSelect(this, 'dna')"/>
	<br>
	<label for="'tabular-aa-seq-checkbox" onclick="handleSequenceSourceSelect(this, 'aa')">AA Sequence</label>
	<input type='checkbox' id='tabular-aa-seq-checkbox' onclick="handleSequenceSourceSelect(this, 'aa')"/>
	<br>
	`
  )

  Object.keys(genomesObj).map((genome, idx) => {
    $('#tabular-modal-nav-tabs').append(`
      <li class="${idx == 0 ? 'active' : ''}"><a data-toggle="tab" href="#${genome}">${genome}</a></li>
    `)
    $('#modal-tab-content').append(`
      <div id="${genome}" class="tab-pane fade in ${idx == 0 ? 'active' : ''}">
        <h3>${genome}</h3>
        <table id="${genome}-table"class="table table-striped" style="width: 100%; text-align: left; font-size: 12px; background: white">
          <thead id='tabular-modal-table-head'>
            <tr id='${genome}-tabular-modal-table-header-tr'>
              <th id="${genome}-select"><input class="form-check-input" type='checkbox'></input>Select</th>
              <th>Gene Caller ID</th>
              <th>Start</th>
              <th>Stop</th>
              <th>Direction</th>
              <th>Contig</th>
              <th id="${genome}-deepdive">Deepdive</th>
              <th id="${genome}-color">Color</th>
              <th id="${genome}-hidden">Hidden?</th>
            </tr>
          </thead>
          <tbody id="${genome}-table-body">
          </tbody>
        </table>
      </div>
    `)
  })
  Object.entries(genomesObj).map((genome, idx) => {
    let totalTableString = String()
    let totalAnnotationsString = String()
    genome[1].forEach(gene => {
      if(gene['functions']){
        Object.entries(gene['functions']).forEach((func, idx) => {
          totalAnnotationsString += `
            <td>
              ${func[1] ? func[1][0].slice(0,20) : 'n/a'} || ${func?.[1]?.[1].slice(0,20)} <button class="btn" type="button" data-toggle="collapse" data-target="#collapse-${idx}" aria-expanded="false" aria-controls="collapse-${idx}">▼</button>
              <div class="collapse" id="collapse-${idx}">${func?.[1]?.[0]}|| ${func?.[1]?.[1]}</div>
            </td>
          `
        })
      }
      let geneHex
      if(settings['display']['colors']['genes']?.[genome[0]]?.[gene['geneID']]){
        geneHex = settings['display']['colors']['genes']?.[genome[0]]?.[gene['geneID']]
      }
      totalTableString += `
      <tr id='${genome[0]}-table-row-${gene['geneID']}'>
        <td class='select'><input class="form-check-input" id="${genome[0]}-${gene['geneID']}" value="${genome[0]}-${gene['geneID']}" type='checkbox'></input></td>
        <td>${gene['geneID']}</td>
        <td>${gene['gene']['start']}</td>
        <td>${gene['gene']['stop']}</td>
        <td>${gene['gene']['direction']}</td>
        <td>${gene['gene']['contig']}</td>
        <td class='deepdive'><button class="btn btn-default btn-sm" id="${genome[0]}-${gene['geneID']}" onclick=transitionTabularModalToDeepdive(event)>Deep Dive</button></td>
        <td class='color'><div id="${gene['geneID']}-${genome[0]}-picker-tabular-modal" class="colorpicker" color="${geneHex ? geneHex : '#808080'}" background-color="${geneHex ? geneHex : '#808080'}" style="background-color: ${geneHex ? geneHex : '#808080'}; margin-right:16px; margin-left:16px"></div></td>
        <td class='hidden?'><input class="form-hidden-check-input"
                   id="${genome[0]}-${gene['geneID']}-hidden" value="${genome[0]}-${gene['geneID']}"
                   onclick='handleGeneShowHideToggle(this)'
                   type='checkbox' ${settings['display']['hidden']?.[genome[0]]?.[gene['geneID']] ? 'checked' : null}>
        </input></td>
      </tr>`
    })
    $(`#${genome[0]}-table-body`).append(totalTableString)
  })

  $('.colorpicker').colpick({
    layout: 'hex',
    submit: 0,
    colorScheme: 'light',
    onChange: function(hsb, hex, rgb, el, bySetColor) {
        $(el).css('background-color', '#' + hex);
        $(el).attr('color', '#' + hex);
        if (!bySetColor) $(el).val(hex);

        let [geneID, genomeID] = [el.id.split('-')[0], el.id.split('-')[1]]

        if(settings['display']['colors']['genes']?.[genomeID]?.[geneID]){
          settings['display']['colors']['genes'][genomeID][geneID] = '#' + hex
        } else {
          settings['display']['colors']['genes'][genomeID] = {}
          settings['display']['colors']['genes'][genomeID][geneID] = '#' + hex
        }

        let arrow = canvas.getObjects().filter(obj => obj.id == 'arrow').find(arrow => arrow.geneID == geneID && arrow.genomeID == genomeID)
        arrow.fill = '#' + hex
        arrow.dirty = true
        canvas.renderAll();
    }
  }).keyup(function() {
      $(this).colpickSetColor(this.value);
  });
}

function gatherTabularModalSelectedItems(action){
  let targetedGenes = []
  $('#modal-tab-content :checked').each(function(){
    let [genome, gene] = $(this).val().split('-')
    targetedGenes.push({genomeID: genome, geneID: gene})
  })

  switch (action) {
    case 'glow':
      drawer.glowGenes(targetedGenes)
      $('#tabular-modal-body').fadeOut(1000).delay(3000).fadeIn(1000)
      break;
    case 'color':
      break;
    case 'metadata':
      break;
    default:
      break;
  }
}

function transitionTabularModalToDeepdive(event){
  let [genomeID, geneID] = event.target.id.split('-')
  let genomeOfInterest = this.settings['genomeData']['genomes'].filter(genome => genome[0] == genomeID)

  // this object generation is to mimick the expected behavior of a click event that the deepdive tooltip expects
  let generatedEventObj = {}
  generatedEventObj.target = {}
  generatedEventObj.target.gene = genomeOfInterest[0][1]['genes']['dna'][geneID]
  generatedEventObj.target.functions = genomeOfInterest[0][1]['genes']['functions'][geneID]
  $('#tabular-modal-body').modal('hide')
  showDeepDiveToolTip(generatedEventObj)
}

function exportTabularModalToTSV(){

  // TODO function to strip out UI columns before processing
  let titles = new Array;
  let data = new Array;
  let active = $("div#modal-tab-content div.active")[0].id

  $(`#${active} th`).each(function() {
    if(['Select', 'Color', 'Hidden?', 'Deepdive'].includes($(this).text())){
      return
    }
    titles.push($(this).text());
  });

  $(`#${active} td`).each(function() {
    if(['select', 'color', 'hidden?', 'deepdive'].includes($(this).attr('class'))){
      return
    }
    data.push($(this).text());
  });

  let CSVString = prepTSVRow(titles, titles.length, '');
  CSVString = prepTSVRow(data, titles.length, CSVString);
  let downloadLink = document.createElement("a");
  let blob = new Blob(["\ufeff", CSVString]);
  let url = URL.createObjectURL(blob);
  downloadLink.href = url;
  downloadLink.download = `${active}-tabular-modal-table.tsv`;
  downloadLink.click()
}

function prepTSVRow(arr, columnCount, initial) { // https://stackoverflow.com/questions/40428850/how-to-export-data-from-table-to-csv-file-using-jquery
  var row = ''; // this will hold data
  var delimeter = '\t'; // data slice separator, in excel it's `;`, in usual CSv it's `,`
  var newLine = '\r\n'; // newline separator for CSV row

  /*
   * Convert [1,2,3,4] into [[1,2], [3,4]] while count is 2
   * @param _arr {Array} - the actual array to split
   * @param _count {Number} - the amount to split
   * return {Array} - splitted array
   */
  function splitArray(_arr, _count) {
    var splitted = [];
    var result = [];
    _arr.forEach(function(item, idx) {
      if ((idx + 1) % _count === 0) {
        splitted.push(item);
        result.push(splitted);
        splitted = [];
      } else {
        splitted.push(item);
      }
    });
    return result;
  }
  var plainArr = splitArray(arr, columnCount);
  plainArr.forEach(function(arrItem) {
    arrItem.forEach(function(item, idx) {
      row += item + ((idx + 1) === arrItem.length ? '' : delimeter);
    });
    row += newLine;
  });
  return initial + row;
}

function showLassoMenu(selected_genes, x, y) {
  //x = 600;
  //y = 200;
  let start, stop, length;
  let showSetLabel = false;
  if(selected_genes.every(obj => obj.genomeID == selected_genes[0].genomeID)) {
    start = selected_genes[0].gene.start;
    stop = selected_genes[selected_genes.length-1].gene.stop;
    length = stop - start;
    showSetLabel = true;
  } else {
    start = stop = length = "N/A";
  }
  $('#lasso-modal-body').modal('show')
  $('#lasso-modal-content').empty().append(
    `
    <table class="table table-striped" style="width: 100%; text-align: center; font-size: 12px; background: white"> \
      <tr><td>Start</td><td>${start}</td></tr> \
      <tr><td>Stop</td><td>${stop}</td></tr> \
      <tr><td>Length</td><td>${length}</td></tr> \
      <tr><td>Gene count</td><td>${selected_genes.length}</td></tr> \
    </table>';

    <div id="picker_lasso" class="colorpicker" color="#808080" background-color="#808080" style="background-color: #808080; margin-right:16px; margin-left:16px"></div>
    <br><br>
    <input id="create_gene_set_label" class="form-control input-sm" type="text" placeholder="New gene set label" style="margin-bottom: 15px; display:${showSetLabel ? "block" : "none"}" size="4">

    <button type="button" id="lasso-done" class="btn btn-default btn-sm" style="float:right" onClick="applyLasso()">Apply</button>
    `
  )

    $('#picker_lasso').colpick({
      layout: 'hex',
      submit: 0,
      colorScheme: 'light',
      onChange: function(hsb, hex, rgb, el, bySetColor) {
          $(el).css('background-color', '#' + hex);
          $(el).attr('color', '#' + hex);
          if (!bySetColor) $(el).val(hex);

          if(!settings['display']['colors']['genes'][selected_genes[0].genomeID])
            settings['display']['colors']['genes'][selected_genes[0].genomeID] = {};
          selected_genes.forEach(gene => {
            gene.fill = '#' + hex;
            gene.dirty = true;
            settings['display']['colors']['genes'][selected_genes[0].genomeID][gene.geneID] = '#' + hex;
          });
          canvas.renderAll();
      }
    }).keyup(function() {
        $(this).colpickSetColor(this.value);
        selected_genes.forEach(gene => {
          gene.fill = this.value;
        });
    });

    if(showSetLabel) {
      $('#create_gene_set_label').on('change', function() {
        createGeneSetLabel(selected_genes, $(this).val());
      });
    }
}

function createGeneSetLabel(selected_genes, title) {
  // assume selected_genes are from the same genomeID
  let genomeID = selected_genes[0].genomeID;
  let geneIDs = selected_genes.map(gene => gene.geneID);

  // if no set labels yet for this genome, initialize empty array
  if(!settings['display']['labels']['set-labels'][genomeID]) {
    let numGenes = Object.keys(settings['genomeData']['genomes'].find(genome => genome[0] == genomeID)[1].genes.gene_calls).length;
    settings['display']['labels']['set-labels'][genomeID] = new Array(numGenes).fill(false);
    settings['display']['labels']['gene-sets'][genomeID] = [];
  } else {
    // if a gene in this selection is already part of a set, do not create a new label
    if(geneIDs.some(
      geneID => settings['display']['labels']['set-labels'][genomeID][geneID]
    )) return;
  }

  // save label to settings for redrawing
  geneIDs.forEach(id => {
    settings['display']['labels']['set-labels'][genomeID][id] = true
  });
  settings['display']['labels']['gene-sets'][genomeID].push([title, geneIDs]);

  drawer.draw();
}

function applyLasso() {
  $('#lasso-menu-body').empty().hide();

  // TODO: apply any relevant changes here.
  // Alternatively, we could apply changes immediately and remove this button.
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
  $("#tbody_genomes").empty() // clear table before redraw
  genomes.map(genome => {
    var height = '50';
    var margin = '15';
    let genomeLabel= genome[0];
    var template = `<tr id=${genomeLabel}>
                  <td><img src="images/drag.gif" class="drag-icon" id=${genomeLabel} /></td>
                  <td> ${genomeLabel} </td>
                  <td>n/a</td>
                  <td>n/a</td>
                  <td>n/a</td>
                  <td><input class="input-height" type="text" size="3" id="height_${genomeLabel}" value="${height}"></input></td>
                  <td class="column-margin"><input class="input-margin" type="text" size="3" id="margin_${genomeLabel}" value="${margin}"></input></td>
                  <td>n/a</td>
                  <td>n/a</td>
                  <td><input type="checkbox" class="layer_selectors"></input></td>
                  </tr>;`

    $('#tbody_genomes').append(template);
  })

  $("#tbody_genomes").sortable({helper: fixHelperModified, handle: '.drag-icon', items: "> tr"}).disableSelection();

  $("#tbody_genomes").on("sortupdate", (event, ui) => {
    changeGenomeOrder($("#tbody_genomes").sortable('toArray'))
  })
}

function buildGroupLayersTable(layerLabel){
  let color;
  let show = settings['display']['layers'][layerLabel]

  if(layerLabel == 'GC_Content') color = settings['display']['colors']['GC_Content'] // TODO: replace hardcoded values w state data
  if(layerLabel == 'Coverage') color = settings['display']['colors']['coverage']

  if(layerLabel === 'Ruler' || layerLabel === 'Genome'){
    var template =  `<tr id=${layerLabel}>
                    <td><img src="images/drag.gif" class="drag-icon" id=${layerLabel} /></td>
                    <td> ${layerLabel} </td>
                    <td style="margin-left: 5px;"> n/a </td>
                    <td><input type="checkbox" class="additional_selectors" id="${layerLabel}-show" onclick="toggleAdditionalDataLayer(event)" checked=${show}></input></td>
                    </tr>;`
  } else {
    var template =  `<tr id=${layerLabel}>
                     <td><img src="images/drag.gif" class="drag-icon" id=${layerLabel} /></td>
                     <td> ${layerLabel} </td>' +
                     <td><div id="${layerLabel}_color" style="margin-left: 5px;" class="colorpicker" style="background-color: ${color}" color=${color} background-color=${color} ></div></td>
                     <td><input type="checkbox" class="additional_selectors" id="${layerLabel}-show" onclick="toggleAdditionalDataLayer(event)" checked=''></input></td>
                     </tr>;`
  }
  $('#tbody_additionalDataLayers').append(template);
  $("#tbody_additionalDataLayers").sortable({helper: fixHelperModified, handle: '.drag-icon', items: "> tr"}).disableSelection();

  $(`#${layerLabel}-show`).prop('checked', show) // needs to trigger after the template layer is appended to the DOM

  $("#tbody_additionalDataLayers").on("sortupdate", (event, ui) => {
    changeGroupLayersOrder($("#tbody_additionalDataLayers").sortable('toArray'))
  })

  $(`#${layerLabel}_color`).colpick({
    layout: 'hex',
    submit: 0,
    colorScheme: 'light',
    onChange: function(hsb, hex, rgb, el, bySetColor) {
        $(el).css('background-color', '#' + hex);
        $(el).attr('color', '#' + hex);
        if(layerLabel == 'Coverage' && settings['display']['colors']['coverage'] != `#${hex}`){
          settings['display']['colors']['coverage'] = `#${hex}`
          drawer.draw()
        }
        else if(layerLabel == 'GC_Content' && settings['display']['colors']['GC_Content'] != `#${hex}`){
          settings['display']['colors']['GC_Content'] = `#${hex}`
          drawer.draw()
        } else {
          $(el).val(hex);
        }
    }
  }).keyup(function() {
      $(this).colpickSetColor(this.value);
  });
}

function toggleAdditionalDataLayer(e){
  let layer = e.target.id.split('-')[0]

  if(e.target.checked){
    settings['display']['layers'][layer] = true
    maxGroupSize += 1
  } else {
    settings['display']['layers'][layer] = false
    maxGroupSize -= 1 // decrease group height if hiding the layer
  }
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
    settings['display']['bookmarks'].push(
      {
        name : $('#create_bookmark_input').val(),
        start : $('#brush_start').val(),
        stop : $('#brush_end').val(),
        description : $('#create_bookmark_description').val(),
      }
    )
    toastr.success('bookmark successfully created :)')
    $('#create_bookmark_input').val('')
    $('#create_bookmark_description').val('')
    $('#bookmarks-select').empty()
    $('#bookmarks-select').prepend('<option>Bookmarks</option>')

    settings['display']['bookmarks'].map(bookmark => {
      $('#bookmarks-select').append((new Option(bookmark['name'], [bookmark["start"], bookmark['stop']])))
    })
  } catch (error) {
    toastr.error(`anvi'o was unable to save your bookmark because of an error ${error} :/`)
    // throw to error landing page?
  }
}
/*
 *  update sequence position, bookmark description upon user select from dropdown
 */
function respondToBookmarkSelect(){
  $('#bookmarks-select').change(function(e){
    let [start, stop] = [$(this).val().split(',')[0], $(this).val().split(',')[1] ]
    if(!stop)return // edge case for empty 'Bookmarks' placeholder select value
    $('#brush_start').val(start);
    $('#brush_end').val(stop);
    try {
      brush.extent([start, stop]);
          brush(d3.select(".brush").transition());
          brush.event(d3.select(".brush").transition());
      let selectedBookmark = settings['display']['bookmarks'].find(bookmark => bookmark.start == start && bookmark.stop == stop)
      $('#bookmark-description').text(selectedBookmark['description'])
      toastr.success("Bookmark successfully loaded")
    } catch (error) {
      toastr.warn(`Unable to load bookmark because of an error ${error}`)
    }
  })
}

/*
 *  respond to ui, redraw with updated group layer order
 */
function changeGroupLayersOrder(updatedOrder){
  settings['group-layer-order'] = updatedOrder
  drawer.draw()
}

/*
 *  respond to ui, redraw with updated genome group order
 */
function changeGenomeOrder(updatedOrder){
  let newGenomeOrder = []
  updatedOrder.map(label => {
    newGenomeOrder.push(settings['genomeData']['genomes'].find(genome => genome[0] == label))
  })
  settings['genomeData']['genomes'] = newGenomeOrder
  drawer.draw()
}


/*
 *  Generates functional annotation color table for a given color palette.
 *
 *  @param fn_colors :       dict matching each category to a hex color code to override defaults
 *  @param fn_type :         string indicating function category type
 *  @param highlight_genes : array of format [{genomeID: 'g01', geneID: 3, color: '#FF0000'}, ...] to override other coloring for specific genes
 *  @param filter_to_window : if true, filters categories to only those shown in the current render window
 *  @param sort_by_count   : if true, sort annotations by # occurrences, otherwise sort alphabetically
 *  @param thresh_count    : int indicating min # occurences required for a given category to be included in the table
 */
function generateColorTable(fn_colors, fn_type, highlight_genes=null, filter_to_window=filter_gene_colors_to_window, sort_by_count=order_gene_colors_by_count, thresh_count = thresh_count_gene_colors) {
  let db;
  if(fn_type == 'Source') {
    db = default_source_colors;
  } else {
    counts = [];

    // Traverse categories
    for(genome of settings['genomeData']['genomes']) {
      let gene_calls = genome[1].genes.gene_calls;
      let gene_funs = genome[1].genes.functions;

      for(let i = 0; i < Object.keys(gene_calls)[Object.keys(gene_calls).length-1]; i++) {
        if(filter_to_window && (gene_calls[i].start < renderWindow[0] || gene_calls[i].stop > renderWindow[1])) continue;
        let cag = getCagForType(gene_funs[i], fn_type);
        counts.push(cag ? cag : "None");
      }
    }

    // Get counts for each category
    counts = counts.reduce((counts, val) => {
      counts[val] = counts[val] ? counts[val]+1 : 1;
      return counts;
    }, {});

    // Filter by count
    let count_removed = 0;
    counts = Object.fromEntries(
      Object.entries(counts).filter(([cag,count]) => {
        if(count < thresh_count) count_removed += count;
        return count >= thresh_count || cag == "None";
      })
    );

    // Save pre-sort order
    let order = {};
    for(let i = 0; i < Object.keys(counts).length; i++) {
      order[Object.keys(counts)[i]] = i;
    }

    // Sort categories
    counts = Object.fromEntries(
      Object.entries(counts).sort(function(first, second) {
        return sort_by_count ? second[1] - first[1] : first[0].localeCompare(second[0]);
      })
    );
    if(count_removed > 0) counts["Other"] = count_removed;

    // Create custom color dict from categories
    db = getCustomColorDict(fn_type, cags=Object.keys(counts), order=order);
    if(Object.keys(db).includes("Other") && !db["Other"]) db["Other"] = "#FFFFFF"; //  bug fix
  }

  // Override default values with any values supplied to fn_colors
  if(fn_colors) {
    Object.keys(db).forEach(cag => { if(Object.keys(fn_colors).includes(cag)) db[cag] = fn_colors[cag] });
  }

  $('#tbody_function_colors').empty();
  Object.keys(db).forEach(category => {
    appendColorRow(fn_type == 'Source' ? category : category + " (" + counts[category] + ")", category, db[category]);
    $('#picker_' + getCleanCagCode(category)).colpick({
        layout: 'hex',
        submit: 0,
        colorScheme: 'light',
        onChange: function(hsb, hex, rgb, el, bySetColor) {
            $(el).css('background-color', '#' + hex);
            $(el).attr('color', '#' + hex);
            // TODO: save new color once state is implemented
            //state[$('#gene_color_order').val().toLowerCase() + '-colors'][el.id.substring(7)] = '#' + hex;
            if (!bySetColor) $(el).val(hex);
            if(!settings.display.colors.genes.annotations[fn_type]) settings.display.colors.genes.annotations[fn_type] = {};
            settings.display.colors.genes.annotations[fn_type][category] = '#' + hex;
        }
    }).keyup(function() {
        $(this).colpickSetColor(this.value);
    });
  });

  /*$('.colorpicker').colpick({
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
      settings['display']['colors']['genes']['annotations'][fn_type][category] = '#' + hex;
  });*/

  if(highlight_genes) {
    let genomes = Object.entries(settings['genomeData']['genomes']).map(g => g[1][0]);
    for(entry of highlight_genes) {
      let genomeID = entry['genomeID'];
      let geneID = entry['geneID'];
      let color = entry['color'];

      if(!genomes.includes(genomeID)) continue;

      let ind = settings['genomeData']['genomes'].findIndex(g => g[0] == genomeID);
      let genes = Object.keys(settings['genomeData']['genomes'][ind][1].genes.gene_calls);
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
 * [TO BE ADDED to genomeview/UI.js OR 'regular' utils.js]
 * Resets function color table to the default set
 *
 * @param fn_colors: If set, resets state to this dictionary instead of the defaults.
 */
function resetFunctionColors(fn_colors=null) {
  // TODO: this should reset color dictionaries in state and then redraw table, once state is implemented
  if($('#gene_color_order') == null) return;
  delete settings['display']['colors']['genes']['annotations'][color_db];
  generateColorTable(fn_colors, color_db);
}

/*
 *  Responds to 'Apply' button in Settings panel under batch coloring
 */
 function batchColor() {
     var rule = $('[name=batch_rule]:checked').val()
     var color = settings['display']['colors']['Batch'];
     var randomize_color = $('#batch_randomcolor').is(':checked');

     let fn_type = $('#gene_color_order').val();
     let dict = getCustomColorDict(fn_type);
     if(counts && Object.keys(counts).includes("Other")) dict["Other"] = $("#picker_Other").attr('color') ? $("#picker_Other").attr('color') : "#FFFFFF";

     Object.keys(dict).forEach(category => {
       if(randomize_color) {
         // color = randomColor();
         // TODO: instead of using default color dict for random colors, dynamically create random color here to avoid similarity to chosen group color
       }
       let code = getCleanCagCode(category);
       if(rule == 'all') {
         $("#picker_" + code).colpickSetColor(color);
       } else if(rule == 'name') {
         if(category.toLowerCase().indexOf($('#name_rule').val().toLowerCase()) > -1) {
           $("#picker_" + code).colpickSetColor(color);
         }
       } else if(rule == 'count') {
         let count = counts[category];
         if (eval(count + unescape($('#count_rule').val()) + " " + parseFloat($('#count_rule_value').val()))) {
           $("#picker_" + code).colpickSetColor(color);
         }
       }
     });

     drawer.draw();
 }

function buildGeneLabelsSelect(){
  let sourcesObj = {}
  settings['genomeData']['genomes'].map(genome => {
    Object.values(genome[1]['genes']['functions']).map(func => {
      Object.keys(func).forEach(source => {
        if(sourcesObj[source]){
          return
        } else {
          sourcesObj[source] = true
        }
      })
    })
  })
  Object.keys(sourcesObj).forEach(source => {
    $("#gene_label_source").append(new Option(source, source));
    // we can also build the dropdown UI element for source-selection in functional querying
    $("#function_search_category").append(new Option(source, source))
  })
  // while we're here, we add 'metadata' to the function_search_category dropdown select
  // TODO refactor naming convention to sequence_search_category, bc we're not just querying functions!
  $('#function_search_category').append(new Option('metadata', 'metadata'))
}

function setGeneVisibilityRange(targetGene, targetGenome){
  let genomeOfInterest = settings['genomeData']['genomes'].filter(genome => genome[0] == targetGenome)
  let maxGeneID = Object.keys(genomeOfInterest[0][1]['genes']['dna']).length
  settings['display']['hidden'][targetGenome] = {}
  let startingPoint = parseInt(targetGene)
  let lowRangeStart = startingPoint - parseInt($('#hide-range-low').val())
  let highRangeStart = startingPoint + parseInt($('#hide-range-high').val())

  lowRangeStart < 0 ? lowRangeStart = 0 : null
  highRangeStart > parseInt(maxGeneID) ? highRangeStart = maxGeneID : null

  for(let i = lowRangeStart; i >= 0; i--){
    settings['display']['hidden'][targetGenome][i] = true
  }
  for(let i = highRangeStart; i <= parseInt(maxGeneID); i++){
    settings['display']['hidden'][targetGenome][i] = true
  }
  $('#deepdive-modal-body').modal('toggle');
  drawer.draw()
}

function showAllHiddenGenes(){
  settings['display']['hidden'] = {}
  drawer.draw()
}