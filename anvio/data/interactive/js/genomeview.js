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

 var VIEWER_WIDTH = window.innerWidth || document.documentElement.clientWidth || document.getElementsByTagName('body')[0].clientWidth;

 var canvas;
 var genomeLabelsCanvas;
 var brush;
 var genomeMax = 0;

 // Settings vars
 // TODO migrate below variables to kvp in state
 var stateData = {};
 var calculatedSpacing; // like spacing, but calculated ;)
 var spacing = 30; // vertical spacing between genomes
 var yOffset = 0 // vertical space between additional data layers
 var showLabels = true; // show genome labels?
 var genomeLabelSize = 15; // font size of genome labels
 var showGeneLabels = true; // show gene labels?
 var geneLabelSize = 40; // font size of gene labels
 var labelSpacing = 30;  // spacing default for genomeLabel canvas
 var draggableGridUnits = 35; // 'snap' to grid for better user feedback
 var showScale = true; // show nt scale?
 var scaleInterval = 100; // nt scale intervals
 var dynamicScaleInterval = true; // if true, scale interval automatically adjusts to zoom level
 var scaleFactor = 1; // widths of all objects are scaled by this value to zoom in/out

 var alignToGC = null;

 var arrowStyle = 1; // gene arrow cosmetics. 1 (default) = 'inspect-page', 2 = thicker arrows, 3 = pentagon

 var color_db;
 var cog_annotated = true, kegg_annotated = false;
 // doesn't make sense to iterate through every gene with dozens of genomes...
 // will need to find an efficient way to find these automatically

var genomeData;

$(document).ready(function() {
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
        async:false,
        success: function(data) {
            genomeData = data;
            console.log("Saved the following data:");
            console.log(data);

            let genomes = Object.entries(genomeData.genomes) // an array of 2d arrays, where each genome[0] is the object key, and genome[1] is the value
            genomeData.genomes = genomes
        }
    });
}

function loadState(){
  let mockStateData = {
    "state-name" : 'default',
    "zoom" : "1.0",
    "arrow-style" : "1",
    "color-genes-by" : "Source",
    "genome-order-method" : [
      {
        "name" : "some value",
        "order" : ['genome1', 'genome2', 'genome3']
      },
      {
        "name" : "some other value",
        "order" : ['genome2', 'genome1', 'genome3']
      }
    ]
  };

  var defer = $.Deferred();
  $('#modLoadState').modal('hide');
  if ($('#loadState_list').val() == null) {
      defer.reject();
      return;
  }

  var state_name = $('#loadState_list').val();
  waitingDialog.show('Requesting state data from the server ...',
      {
          dialogSize: 'sm',
          onHide: function() {
              defer.resolve();
          },
          onShow: function() {
              $.ajax({
                      type: 'GET',
                      cache: false,
                      url: '/data/genome_view/state/get/' + state_name,
                      success: function(response) {
                          try{
                              // processState(state_name, response[0]); process actual response from backend
                              processState(state_name, mockStateData); // process mock state data
                          }catch(e){
                              console.error("Exception thrown", e.stack);
                              toastr.error('Failed to parse state data, ' + e);
                              defer.reject();
                              return;
                          }
                          waitingDialog.hide();
                      }
                  });
          },
      }
  );

  return defer.promise();
}

function processState(stateName, stateData){
  // set genome order options from state

  // mock data
  stateData['genome-order-method'] = [{
      'name' : 'cats',
      'ordering' : 'some order'
    }, {
      'name' : 'dogs',
      'ordering' : 'some other order'
    }, {
      'name' : 'birds',
      'ordering' : 'beaks to tails'
    }
  ]
  stateData['additional-data-layers'] = []

  calculateMaxGenomeLength()

  for(let i = 0; i < genomeData.genomes.length; i++){ // generate mock additional data layer content
    let gcContent = []
    let coverage = []

    for(let j = 0; j < genomeMax; j++){
      gcContent.push(Math.floor(Math.random() * 45))
      coverage.push(Math.floor(Math.random() * 45))
    }
    let genomeLabel = Object.keys(genomeData.genomes[i][1]['contigs']['info'])[0];
    let additionalDataObject = {
      'genome' : genomeLabel,
      'coverage' : coverage,
      'gcContent' : gcContent
    }
    stateData['additional-data-layers'].push(additionalDataObject)
  }

  if(stateData['genome-order-method']){
    stateData['genome-order-method'].forEach(orderMethod => {
      $('#genome_order_select').append((new Option(orderMethod["name"], orderMethod["name"]))) // set display + value of new select option.
    })
  }

  if(stateData['some-data']){
    state['some-data'] = stateData['some-data']
  }
  if(stateData['some-other-data']){
    state['some-other-data'] = stateData['some-other-data']
  }
}

function loadAll() {
  buildGenomesTable(genomeData.genomes, 'alphabetical') // hardcode order method until backend order data is hooked in
  canvas = new fabric.Canvas('myCanvas');
  canvas.setWidth(VIEWER_WIDTH * 0.85);

  // Find max length genome
  calculateMaxGenomeLength()
  calculatedSpacing = calculateSpacingForGroups()

  var scaleWidth = canvas.getWidth();
  var scaleHeight = 200;

  var xScale = d3.scale.linear().range([0, scaleWidth]).domain([0,genomeMax]);

  var scaleAxis = d3.svg.axis()
              .scale(xScale)
              .tickSize(scaleHeight);

  var scaleArea = d3.svg.area()
              .interpolate("monotone")
              .x(function(d) { return xScale(d); })
              .y0(scaleHeight)
              .y1(0);

  brush = d3.svg.brush()
              .x(xScale)
              .on("brushend", onBrush);

  $("#scaleSvg").attr("width", scaleWidth + 10);

  var scaleBox = d3.select("#scaleSvg").append("g")
              .attr("id", "scaleBox")
              .attr("class","scale")
              .attr("y", 230) // rather than 80 from 50?
              .attr("transform", "translate(5,0)");

  scaleBox.append("g")
              .attr("class", "x axis top noselect")
              .attr("transform", "translate(0,0)")
              .call(scaleAxis);

  scaleBox.append("g")
              .attr("class", "x brush")
              .call(brush)
              .selectAll("rect")
              .attr("y", 0)
              .attr("height", scaleHeight);

  $('#brush_start').val(0);
  $('#brush_end').val(Math.floor(scaleWidth));

  function onBrush(){
      var b = brush.empty() ? xScale.domain() : brush.extent();

      if (brush.empty()) {
          $('.btn-selection-sequence').addClass('disabled').prop('disabled', true);
      } else {
          $('.btn-selection-sequence').removeClass('disabled').prop('disabled', false);
      }

      b = [Math.floor(b[0]), Math.floor(b[1])];

      $('#brush_start').val(b[0]);
      $('#brush_end').val(b[1]);

      let ntsToShow = b[1] - b[0];
      scaleFactor = canvas.getWidth()/ntsToShow;

      if(dynamicScaleInterval) adjustScaleInterval();

      draw();
      let moveToX = (showLabels?120:0) + b[0];
      canvas.absolutePan({x: scaleFactor*moveToX, y: 0});

      // TODO: restrict min view to 300 NTs? (or e.g. scaleFactor <= 4)
  }

  $('#tooltip-body').hide() // set initual tooltip hide value
  $('#show_genome_labels_box').attr("checked", showLabels);
  $('#show_gene_labels_box').attr("checked", showGeneLabels);
  $('#show_scale_box').attr("checked", showScale);
  $('#show_dynamic_scale_box').attr("checked", dynamicScaleInterval);

  // can either set it on the canvas to check for all arrows, or when arrow is created.
  canvas.on('mouse:down', function(options) {
    if(options.target && options.target.id == 'arrow') {
      options.target.set('fill', options.target.fill=="red" ? "blue" : 'red');
      var testid = options.target.gene.gene_callers_id;
      console.log(options.target.gene);
    }
  });

  canvas.on('mouse:over', (event) => {
    if(event.target && event.target.id === 'arrow'){
      showToolTip(event)
    }
  })

  canvas.on('mouse:out', (event) => {
    $('#tooltip-body').html('').hide()
  })

  if(showGeneLabels && arrowStyle != 3) {
    spacing = 60;
    $("#genome_spacing").val(spacing);
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
        <tr><td>Gene Cluster</td><td>${genomeData.gene_associations["anvio-pangenome"]["genome-and-gene-names-to-gene-clusters"][event.target.genomeID][event.target.geneID]}</td></tr>
      </table>
      <button>some action</button>
      <button>some other action</button>
      `).css({'position' : 'absolute', 'left' : event.e.clientX, 'top' : event.e.clientY })
  }

  $('#gene_color_order').append($('<option>', {
    value: 'Source',
    text: 'Source'
  }));
  if(cog_annotated) {
    $('#gene_color_order').append($('<option>', {
      value: 'COG',
      text: 'COG'
    }));
  }
  if(kegg_annotated) {
    $('#gene_color_order').append($('<option>', {
      value: 'KEGG',
      text: 'KEGG'
    }));
  }

  draw();

  // zooming and panning
  // http://fabricjs.com/fabric-intro-part-5#pan_zoom
  canvas.on('mouse:down', function(opt) {
    var evt = opt.e;
    if (evt.shiftKey === true) {
      this.isDragging = true;
      this.selection = false;
      this.lastPosX = evt.clientX;
    }
  });
  canvas.on('mouse:move', function(opt) {
    if (this.isDragging) {
      var e = opt.e;
      var vpt = this.viewportTransform;
      vpt[4] += e.clientX - this.lastPosX;
      this.requestRenderAll();
      this.lastPosX = e.clientX;
    }
  });
  canvas.on('mouse:up', function(opt) {
    // on mouse up we want to recalculate new interaction
    // for all objects, so we call setViewportTransform
    this.setViewportTransform(this.viewportTransform);
    this.isDragging = false;
    this.selection = true;
  });
  canvas.on('mouse:wheel', function(opt) {
    opt.e.preventDefault();
    opt.e.stopPropagation();

    var delta = opt.e.deltaY;
    scaleFactor *= 0.999 ** delta;
    if (scaleFactor > 4) scaleFactor = 4;
    if (scaleFactor < 0.01) scaleFactor = 0.01;
    if(dynamicScaleInterval) adjustScaleInterval();
    draw();
  });

  $('#geneClusterInput').on('keydown', function(e) {
    if(e.keyCode == 13) { // 13 = enter key
      alignToCluster($(this).val());
      $(this).blur();
    }
  });
  document.body.addEventListener("keydown", function(ev) {
    if(ev.which == 83 && ev.target.nodeName !== 'TEXTAREA' && ev.target.nodeName !== 'INPUT') { // S = 83
      toggleSettingsPanel();
    }
  });
  $('#genome_spacing').on('keydown', function(e) {
    if(e.keyCode == 13) { // 13 = enter key
      setGenomeSpacing($(this).val());
      $(this).blur();
    }
  });
  $('#gene_label').on('keydown', function(e) {
    if(e.keyCode == 13) { // 13 = enter key
      setGeneLabelSize($(this).val());
      $(this).blur();
    }
  });
  $('#genome_label').on('keydown', function(e) {
    if(e.keyCode == 13) { // 13 = enter key
      setGenomeLabelSize($(this).val());
      $(this).blur();
    }
  });
  $('#genome_scale_interval').on('keydown', function(e) {
    if(e.keyCode == 13) { // 13 = enter key
      setScale($(this).val());
      $(this).blur();
    }
  });
  $('#gene_color_order').on('change', function() {
      color_db = $(this).val();
      draw();
      $(this).blur();
  });
  $('#arrow_style').on('change', function() {
      arrowStyle = parseInt($(this).val());
      draw();
      $(this).blur();
  });
  $('#show_genome_labels_box').on('change', function() {
    showLabels = !showLabels;
    alignToGC = null;
    draw();
  });
  $('#show_gene_labels_box').on('change', function() {
    showGeneLabels = !showGeneLabels;
    draw();
  });
  $('#show_scale_box').on('change', function() {
    showScale = !showScale;
    draw();
  });
  $('#show_dynamic_scale_box').on('change', function() {
    dynamicScaleInterval = !dynamicScaleInterval;
  });
  $('#brush_start, #brush_end').keydown(function(ev) {
      if (ev.which == 13) { // enter key
          let start = parseInt($('#brush_start').val());
          let end = parseInt($('#brush_end').val());

          if (isNaN(start) || isNaN(end) || start < 0 || start > genomeMax || end < 0 || end > genomeMax) {
              alert(`Invalid value, value needs to be in range 0-${genomeMax}.`);
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

function draw(scaleX=scaleFactor) {
  canvas.clear()
  labelSpacing = 30 // reset to default value upon each draw() call
  yOffset = 0 // reset 
  var y = 1;
  for(genome of genomeData.genomes) {
    let label = genome[1].genes.gene_calls[0].contig;
    addGenome(label, genome[1].genes.gene_calls, genome[0], y, scaleX=scaleX)
    addLayers(label, genome[1], genome[0])
    labelSpacing += 30
    y++;
  }
  drawScale(y);
  shadeGeneClusters(["GC_00000034","GC_00000097","GC_00000002"],{"GC_00000034":"green","GC_00000097":"red","GC_00000002":"purple"},spacing);
  checkGeneLabels();
}

function drawScale(y, scaleX=scaleFactor) {
  if(!showScale) return;

  for(var w = 0; w < genomeMax; w+=scaleInterval) {
    canvas.add(new fabric.Line([0,0,0,20], {left: (w+(showLabels?120:0))*scaleX,
          top: y*(spacing)-24,
          stroke: 'black',
          strokeWidth: 1,
          fontSize: 10,
          fontFamily: 'sans-serif',
          selectable: false}));

    canvas.add(new fabric.Text(w/1000 + " kB", {left: (w+5+(showLabels?120:0))*scaleX,
          top: y*(spacing)-24,
          stroke: 'black',
          strokeWidth: .25,
          fontSize: 15,
          fontFamily: 'sans-serif',
          selectable: false}));
  }

  canvas.add(new fabric.Line([0,0,100,0], {left: (showLabels?120:0)*scaleX,
        top: y*(1.25*spacing)-4,
        stroke: 'black',
        strokeWidth: 2,
        selectable: false}));
  canvas.add(new fabric.Text("100 nts", {left: (15+(showLabels?120:0))*scaleX,
        top: y*(1.25*spacing)-4,
        stroke: 'black',
        strokeWidth: 1,
        fontSize: 20,
        fontFamily: 'sans-serif',
        selectable: false}));
}

function zoomIn() {
  scaleFactor += (scaleFactor < 0.2 ? .01 : .1);
  if(scaleFactor > 4) scaleFactor = 4;
  if(dynamicScaleInterval) adjustScaleInterval();

  draw();
}

function zoomOut() {
  scaleFactor -= (scaleFactor < 0.2 ? .01 : .1);
  if(scaleFactor < 0.01) scaleFactor = 0.01;
  if(dynamicScaleInterval) adjustScaleInterval();

  draw();
}

function shadeGeneClusters(geneClusters, colors, y) {
  for(var i = 0; i < genomeData.genomes.length-1; i++) {
    let genomeA = genomeData.genomes[i][1].genes.gene_calls;
    let genomeB = genomeData.genomes[i+1][1].genes.gene_calls;
    let genomeID_A = genomeData.genomes[i][0];
    let genomeID_B = genomeData.genomes[i+1][0];

    for(gc of geneClusters) {
      let g1 = [], g2 = [];
      for(geneID of genomeData.gene_associations["anvio-pangenome"]["gene-cluster-name-to-genomes-and-genes"][gc][genomeID_A]) {
        g1.push(genomeA[geneID].start, genomeA[geneID].stop);
      }
      for(geneID of genomeData.gene_associations["anvio-pangenome"]["gene-cluster-name-to-genomes-and-genes"][gc][genomeID_B]) {
        g2.push(genomeB[geneID].start, genomeB[geneID].stop);
      }

      g1 = g1.map(val => val*scaleFactor);
      g2 = g2.map(val => val*scaleFactor);

      /* TODO: implementation for multiple genes of the same genome in the same gene cluster */
      var path = new fabric.Path("M " + g1[0] + " " + y + " L " + g1[1] + " " + y + " L " + g2[1] + " " + (y+spacing) + " L " + g2[0] + " " + (y+spacing) + " z", {
        fill: colors[gc],
        opacity: 0.25,
        selectable: false
      });
      path.left += (showLabels?120:0)*scaleFactor;
      path.sendBackwards();
      canvas.add(path);
    }
    y += spacing
  }
}

function checkGeneLabels() {
  var labels = canvas.getObjects().filter(obj => obj.id == 'geneLabel');
  for(var i = 0; i < labels.length-1; i++) {
    if(arrowStyle == 3) {
      if(labels[i].width/2 > canvas.getObjects().filter(obj => obj.id == 'arrow')[i].width) {
        labels[i].visible = false;
        continue;
      }
      labels[i].visible = true;
    }
    var p = i+1;
    while(labels[i].intersectsWithObject(labels[p]) && p < labels.length-1) {
      labels[p].visible = false;
      p++;
    }
    labels[p].visible = true;
    i = p - 1;
  }
}

function alignToCluster(gc) {
  if(!gc || gc in genomeData.gene_associations["anvio-pangenome"]["gene-cluster-name-to-genomes-and-genes"]) {
    alignToGC = gc;

    for(genome of genomeData.genomes) {
      var genomeGCs = genomeData.gene_associations["anvio-pangenome"]["gene-cluster-name-to-genomes-and-genes"][alignToGC][genome[0]];
      if(genomeGCs.length == 0) continue;
      var targetGeneID = genomeGCs[0]; /* TODO: implementation for multiple matching gene IDs */
      var targetGene = genome[1].genes.gene_calls[targetGeneID];
      var genePos = targetGene.start + (targetGene.stop - targetGene.start) / 2;
      var windowCenter = canvas.getWidth()/2 - canvas.viewportTransform[4]; // canvas.getWidth()/2 is clientX of center of screen
      var shift = windowCenter - genePos;
      canvas.viewportTransform[4] += (shift*scaleFactor);
      break;
    }

  } else {
    console.log('Warning: ' + gc + ' is not a gene cluster in data structure');
  }
  draw();
}

function setGenomeSpacing(newSpacing) {
  if(isNaN(newSpacing)) return;
  newSpacing = parseInt(newSpacing);
  if(newSpacing < 0 || newSpacing > 1000) {
    alert(`Invalid value, genome spacing must be in range 0-1000.`);
    return;
  }
  spacing = newSpacing;
  draw();
}

function setScale(newScale) {
  if(isNaN(newScale)) return;
  newScale = parseInt(newScale);
  if(newScale < 50) {
    alert(`Invalid value, scale interval must be >=50.`);
    return;
  }
  scaleInterval = newScale;
  draw();
}

function setGeneLabelSize(newSize) {
  if(isNaN(newSize)) return;
  newSize = parseInt(newSize);
  if(newSize < 0 || newSize > 1000) {
    alert(`Invalid value, gene label size must be in range 0-1000.`);
    return;
  }
  geneLabelSize = newSize;
  if(showGeneLabels) draw();
}

function setGenomeLabelSize(newSize) {
  if(isNaN(newSize)) return;
  newSize = parseInt(newSize);
  if(newSize < 0 || newSize > 1000) {
    alert(`Invalid value, genome label size must be in range 0-1000.`);
    return;
  }
  genomeLabelSize = newSize;
  if(showLabels) draw();
}

function addGenome(label, gene_list, genomeID, y, scaleX=1) {
  if(showLabels) {
    canvas.add(new fabric.Text(label, {top: spacing*y-5, selectable: false, fontSize: genomeLabelSize, fontFamily: 'sans-serif', fontWeight: 'bold'}));
  }

  // line
  canvas.add(new fabric.Line([0,0,genomeMax*scaleX,0], {left: (showLabels?120:0),
        top: spacing*y + 4,
        stroke: 'black',
        strokeWidth: 2,
        selectable: false}));

  // here we can either draw genes individually, or collectively as a group
  // grouping allows you to transform them all at once, but makes selecting individual arrows difficult
  // so for now they are drawn individually

  //var geneGroup = new fabric.Group();
  for(let geneID in gene_list) {
    let gene = gene_list[geneID];
    //geneGroup.addWithUpdate(geneArrow(gene,y));   // IMPORTANT: only way to select is to select the group or use indices. maybe don't group them but some alternative which lets me scale them all at once?
    var geneIndex = genomeData.genomes.findIndex(g => genomeID == g[0]);
    var geneObj = geneArrow(gene,geneID,genomeData.genomes[geneIndex][1].genes.functions[geneID],y,genomeID,arrowStyle,scaleX=scaleX);
    if(showLabels) {
      geneObj.left += 120*scaleX;
    }
    canvas.add(geneObj);

    if(showGeneLabels) {
      var label = new fabric.IText("geneID: "+geneID, {
        id: 'geneLabel',
        fontSize: geneLabelSize,
        hasControls: false,
        lockMovementX: true,
        lockMovementY: true,
        lockScaling: true,
        hoverCursor: 'text'
      });

      if(arrowStyle == 3) {
        label.set({
          top: -5+spacing*y,
          left: (150+gene.start)*scaleX,
          scaleX: 0.5,
          scaleY: 0.5,
          selectionColor:'rgba(128,128,128,.5)'
        });
      } else {
        label.set({
          scaleX: 0.5,
          scaleY: 0.5,
          top: -30+spacing*y,
          left: (200+gene.start)*scaleX,
          angle: -10,
          selectionColor:'rgba(128,128,128,.2)'
        });
      }
      canvas.add(label);
    }
  }
  //canvas.add(geneGroup.set('scaleX',canvas.getWidth()/genomeMax/3));
  //geneGroup.destroy();
}

function addLayers(label, genome, genomeID){ // this will work alongside addGenome to render out any additional data layers associated with each group (genome)
  let additionalDataLayers = stateData['additional-data-layers'].find(group =>  group.genome = label)
  if(additionalDataLayers['coverage']){
    let maxCoverageValue = 0
    for(let i = 0; i < additionalDataLayers['coverage'].length; i++){ // TODO this seems like an inefficient way to find the max range for coverage
      additionalDataLayers['coverage'][i] > maxCoverageValue ? maxCoverageValue = additionalDataLayers['coverage'][i] : null 
    }
    for(let i = 0; i < 1000; i++){ // TODO hardcoded to 1k for now because canvas ~really~ doesn't like rendering out 130k * 3 objects (maybe)
      canvas.add(new fabric.Rect({
        height : 8, 
        width : 1*scaleFactor, 
        fill : 'blue', 
        opacity : 1 * [additionalDataLayers['coverage'][i] / maxCoverageValue],
        selectable : false, 
        top : 70 + yOffset, 
        left : 120 + i
      }) )
    }
  } 
  if(additionalDataLayers['gcContent']){
    let maxGCValue = 0
    for(let i = 0; i < additionalDataLayers['gcContent'].length; i++){ 
      additionalDataLayers['gcContent'][i] > maxGCValue ? maxGCValue = additionalDataLayers['gcContent'][i] : null 
    }
    for(let i = 0; i < 1000; i++){ // 
      canvas.add(new fabric.Rect({
        height : 8, 
        width : 1*scaleFactor, 
        fill : 'pink', 
        opacity : 1 * [additionalDataLayers['gcContent'][i] / maxGCValue],
        selectable : false, 
        top : 70 + yOffset + 12, 
        left : 120 + i
      }) )
    }
  } 
  yOffset += 60
}

function geneArrow(gene, geneID, functions, y, genomeID, style, scaleX=1) {
  var cag = null;
  var color = 'gray';
  if(functions) {
    switch(color_db) {
      case 'COG':
        if(functions["COG_CATEGORY"]) cag = functions["COG_CATEGORY"][1];
        color = cag in default_COG_colors ? default_COG_colors[cag] : 'gray';
        break;
      case 'KEGG':
        if(functions.hasOwnProperty("KEGG_Class") && functions.KEGG_Class != null) {
          cag = getCategoryForKEGGClass(functions["KEGG_Class"][1]);
        }
        color = cag in default_KEGG_colors ? default_KEGG_colors[cag] : 'gray';
        break;
      default:
        if (gene.source.startsWith('Ribosomal_RNA')) {
          cag = 'rRNA';
        } else if (gene.source == 'Transfer_RNAs') {
          cag = 'tRNA';
        } else if (gene.functions !== null) {
          cag = 'Function';
        }
        color = cag in default_source_colors ? default_source_colors[cag] : 'gray';
    }
  }
  /* Issue here: each genome might be differentially annotated... how to make sure all have COG annotations for example? */

  let length = (gene.stop-gene.start)*scaleX;

  var arrowPathStr;
  switch(style) {
    case 2: // thicker arrows
      arrowPathStr = 'M ' + (length-25) + ' -5 L 0 -5 L 0 15 L ' + (length-25) + ' 15 L ' + (length-25) + ' 15 L ' + (length-25) + ' 20 L ' + length + ' 5 L ' + (length-25) + ' -10 z';
      break;
    case 3: // pentagon arrows
      arrowPathStr = 'M 0 0 L ' + (length-25) + ' 0 L ' + length + ' 20 L ' + (length-25) + ' 40 L 0 40 L 0 0 z';
      break;
    default: // 'inspect page' arrows
      arrowPathStr = 'M ' + (length-25) + ' 0 L 0 0 L 0 10 L ' + (length-25) + ' 10 L ' + (length-25) + ' 10 L ' + (length-25) + ' 20 L ' + length + ' 5 L ' + (length-25) + ' -10 z';
      break;
  }

  var arrow = new fabric.Path(arrowPathStr);
  arrow.set({
    id: 'arrow',
    selectable: false,
    gene: gene,
    geneID: geneID,
    genomeID: genomeID,
    top: style == 3 ? -17+spacing*y : -11+spacing*y,
    left: (1.5+gene.start)*scaleX,
    fill: color,
    stroke: 'gray',
    strokeWidth: style == 3 ? 3 : 1.5
  });
  if(gene.direction == 'r') arrow.rotate(180);

  return arrow;
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

function getCategoryForKEGGClass(class_str) {
  if(class_str == null) return null;

  var category_name = getClassFromKEGGAnnotation(class_str);
  return getKeyByValue(KEGG_categories, category_name);
}

function getClassFromKEGGAnnotation(class_str) {
  return class_str.substring(17, class_str.indexOf(';', 17));
}

// https://stackoverflow.com/questions/9907419/how-to-get-a-key-in-a-javascript-object-by-its-value/36705765
function getKeyByValue(object, value) {
  return Object.keys(object).find(key => object[key] === value);
}

function clamp(num, min, max) {
  return Math.min(Math.max(num, min), max);
}

function resetScale(){
  canvas.setZoom(1)
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
  draw()
}

function calculateMaxGenomeLength(){
  for(genome of genomeData.genomes) {
    genome = genome[1].genes.gene_calls;
    let genomeEnd = genome[Object.keys(genome).length-1].stop;
    if(genomeEnd > genomeMax) genomeMax = genomeEnd;
  }
}

function calculateSpacingForGroups(){ // to be used for setting vertical spacing
  let maxGroupSize = 1; // default, as each group will always have at minimum a 'genome' layer
  stateData['additional-data-layers'].map(group => {
    Object.keys(group).length > maxGroupSize ? maxGroupSize = Object.keys(group).length : null
  })
  let spacing = 500 / [maxGroupSize * genomeData.genomes.length] // 500 is hardcoded main canvas height
  return spacing
}

function adjustScaleInterval() { // dynamically set scale interval based on scaleFactor
  let val = Math.floor(100/scaleFactor);
  let roundToDigits = Math.floor(Math.log10(val)) - 1;
  let newInterval = Math.floor(val/(10**roundToDigits)) * (10**roundToDigits);
  scaleInterval = newInterval;
  $('#genome_scale_interval').val(scaleInterval);
}

var fixHelperModified = function(e, tr) { // ripped from utils.js instead of importing the whole file
  var $originals = tr.children();
  var $helper = tr.clone();
  $helper.children().each(function(index) {
      $(this).width($originals.eq(index).width());
  });
  return $helper;
};
