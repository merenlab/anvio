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
 var scaleCanvas;
 var genomeLabelsCanvas;
 var genomeMax = 0;

 // Settings vars
 // TODO migrate below variables to kvp in state
 var state = {}; 
 var spacing = 30; // vertical spacing between genomes
 var showLabels = true; // show genome labels?
 var showGeneLabels = true; // show gene labels?
 var labelSpacing = 30;  // spacing default for genomeLabel canvas
 var draggableGridUnits = 35; // 'snap' to grid for better user feedback
 var showScale = true; // show nt scale?
 var scale = 100; // nt scale intervals

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
            
            // implementing the following will require a refactor of associated methods
            // genomeData.genomes = genomes
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

  if(stateData['genome-order-method']){
    stateData['genome-order-method'].forEach(orderMethod => {
      $('#genomeOrderSelect').append((new Option(orderMethod["name"], orderMethod["name"]))) // set display + value of new select option. 
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
  canvas = new fabric.Canvas('myCanvas');
  genomeLabelsCanvas = new fabric.Canvas('genomeLabels');
  scaleCanvas = new fabric.Canvas('scale')
  scaleCanvas.add(new fabric.Rect({
    width: 1200,
    height : 200,
    fill : 'gray',
    opacity : .6,
    selectable : false,

  }))

  scaleCanvas.add(new fabric.Line([0,200,0,0],{ // 'cursor' for current zoom location
    strokeWidth : 10,
    stroke : 'green',
    top : 0,
    left : 10,
    opacity : .6, 
    selectable : false 
  }))
  let scaleDragStartingX;
  let totalScaleX

  scaleCanvas.on('mouse:down', (event) => {
    scaleDragStartingX = event.pointer.x
  })

  scaleCanvas.on('mouse:up', (event) => {
    let scaleDragEndingX = event.pointer.x // click + drag ending x position
    totalScaleX = event.target.aCoords.tr.x // total x axis length

    if(scaleDragStartingX === scaleDragEndingX){ //account for accidental drag, click and release
      return
    }
    let mainCanvasWidth = canvas.getWidth()
    let [percentileStart, percentileEnd] = [scaleDragStartingX/totalScaleX, scaleDragEndingX/totalScaleX]
    let percentileToShow = (percentileEnd - percentileStart).toFixed(2)
 
    let moveTo = new fabric.Point(mainCanvasWidth * percentileStart, 0)
    // process for contextual zooming based on drag event from scale canvas.
    // 1) zoom all the way out so the main canvas theoretically displays the entire sequence length
    // 2) pan to new x coordinate, a percentile of total sequence length
    // 3) set zoom back to original value, centering on new x coordinate

    scaleCanvas.remove(scaleCanvas.getObjects()[1]) // this should be changed
    scaleCanvas.add(new fabric.Line([0,200,0,0],{ // 'cursor' for current zoom location
      strokeWidth : percentileToShow > .1 ? 1200*percentileToShow : 10,
      stroke : 'green',
      top : 0,
      left : scaleCanvas.width * percentileStart,
      opacity : .6,
      selectable : false 
    }))
    canvas.setZoom(.01) 
    canvas.absolutePan({x: moveTo.x, y: moveTo.y})
    canvas.setZoom(percentileToShow > .1 ? 1 - percentileToShow : 1) //TODO ask someone who knows how to math to make this better :)
  })

  $('#tooltip-body').hide() // set initual tooltip hide value
  $('#toggle_label_box').attr("checked", showLabels);

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

  genomeLabelsCanvas.on('object:moving', function(options) {
    // console.log(genomeLabelsCanvas.getBoundingRect())
    // console.log(options.target.getBoundingRect())
    options.target.set({
      left: Math.round(options.target.left / draggableGridUnits) * draggableGridUnits,
      top: Math.round(options.target.top / draggableGridUnits) * draggableGridUnits
    });
    let labelsArr = genomeLabelsCanvas.getObjects()
    console.log('pre sort ==>', labelsArr)
    labelsArr.sort((a,b) => (a.top > b.top) ? 1 : -1)
    console.log('after sort ==>', labelsArr)
  });

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
      this.lastPosY = evt.clientY;
    }
  });
  canvas.on('mouse:move', function(opt) {
    if (this.isDragging) {
      var e = opt.e;
      var vpt = this.viewportTransform;
      vpt[4] += e.clientX - this.lastPosX;
      vpt[5] += e.clientY - this.lastPosY;
      this.requestRenderAll();
      this.lastPosX = e.clientX;
      this.lastPosY = e.clientY;
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
    var delta = opt.e.deltaY;
    var zoom = canvas.getZoom();
    zoom *= 0.999 ** delta;
    if (zoom > 20) zoom = 20;
    if (zoom < 0.01) zoom = 0.01;

    canvas.zoomToPoint({ x: opt.e.offsetX, y: opt.e.offsetY }, zoom);
    opt.e.preventDefault();
    opt.e.stopPropagation();
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
  $('#gene_color_order').on('change', function() {
      color_db = $(this).val();
      canvas.clear();
      draw();
      $(this).blur();
  });
  $('#arrow_style').on('change', function() {
      arrowStyle = parseInt($(this).val());
      canvas.clear();
      draw();
      $(this).blur();
  });
  $('#toggle_label_box').on('change', function() {
    showLabels = !showLabels;
    alignToGC = null;
    canvas.clear();
    draw();
  });
}

function draw() {
  // Find max length genome
  for(let genomeID in genomeData.genomes) {
    let genome = genomeData.genomes[genomeID].genes.gene_calls;
    let g = genome[Object.keys(genome).length-1].stop;
    if(g > genomeMax) genomeMax = g;
  }

  if(showGeneLabels && arrowStyle != 3) {
    spacing = 60;
  }

  var y = 1;
  for(let genomeID in genomeData.genomes) {
    let genome = genomeData.genomes[genomeID].genes.gene_calls;
    let label = genome[0].contig;
    addGenome(label, genome, genomeID, y)
    drawGenomeLabels(label)
    labelSpacing += 30
    y++;
  }

  if(showScale) {
    for(var w = 0; w < genomeMax; w+=scale) {
      canvas.add(new fabric.Line([0,0,0,20], {left: w+(showLabels?120:0),
            top: y*(spacing)-24,
            stroke: 'black',
            strokeWidth: 1,
            fontSize: 10,
            fontFamily: 'sans-serif',
            selectable: false}));

      canvas.add(new fabric.Text(w/1000 + " kB", {left: w+5+(showLabels?120:0),
            top: y*(spacing)-24,
            stroke: 'black',
            strokeWidth: .25,
            fontSize: 15,
            fontFamily: 'sans-serif',
            selectable: false}));
    }

    canvas.add(new fabric.Line([0,0,100,0], {left: (showLabels?120:0),
          top: y*(1.25*spacing)-4,
          stroke: 'black',
          strokeWidth: 2,
          selectable: false}));
    canvas.add(new fabric.Text("100 nts", {left: 15+(showLabels?120:0),
          top: y*(1.25*spacing)-4,
          stroke: 'black',
          strokeWidth: 1,
          fontSize: 20,
          fontFamily: 'sans-serif',
          selectable: false}));
  }
}

function shadeGeneClusters(geneClusters, colors, y) {
  for(var i = 0; i < Object.keys(genomeData.genomes).length-1; i++) {
    let genomeID_A = Object.keys(genomeData.genomes)[i];
    let genomeID_B = Object.keys(genomeData.genomes)[i+1];
    //let genomeA = genomeData.genomes[genomeID_A].genes.gene_calls;
    //let genomeB = genomeData.genomes[genomeID_B].genes.gene_calls;

    for(gc of geneClusters) {
      let g1 = [], g2 = [];
      for(geneID of genomeData.gene_associations["anvio-pangenome"]["gene-cluster-name-to-genomes-and-genes"][gc][genomeID_A]) {
        let g = genomeData.genomes[genomeID_A].genes.gene_calls[geneID];
        g1.push(g.start, g.stop);
      }
      for(geneID of genomeData.gene_associations["anvio-pangenome"]["gene-cluster-name-to-genomes-and-genes"][gc][genomeID_B]) {
        let g = genomeData.genomes[genomeID_B].genes.gene_calls[geneID];
        g2.push(g.start, g.stop);
      }

      /* TODO: implementation for multiple genes of the same genome in the same gene cluster */
      var path = new fabric.Path("M " + g1[0] + " " + y + " L " + g1[1] + " " + y + " L " + g2[1] + " " + (y+spacing) + " L " + g2[0] + " " + (y+spacing) + " z", {
        fill: colors[gc],
        opacity: 0.25,
        scaleX: 0.5,
        selectable: false
      });
      path.left += (showLabels?120:0);
      path.sendBackwards();
      canvas.add(path);
    }
    y += spacing
  }
}

function alignToCluster(gc) {
  if(!gc || gc in genomeData.gene_associations["anvio-pangenome"]["gene-cluster-name-to-genomes-and-genes"]) {
    alignToGC = gc;
    showLabels = !gc; // only show labels if changing to default view
    $('#toggle_label_box').attr("checked", showLabels);
  } else {
    console.log('Warning: ' + gc + ' is not a gene cluster in data structure');
  }
  canvas.clear();
  draw();
}

function addGenome(label, gene_list, genomeID, y) {
  var offsetX = 0;
  if(alignToGC) { /* TEMPORARILY DISABLED until proper data structure (geneID -> GC) added */
    var genomeGCs = genomeData.gene_associations["anvio-pangenome"]["gene-cluster-name-to-genomes-and-genes"][alignToGC][genomeID];
    var targetGeneID = genomeGCs[0]; /* TODO: implementation for multiple matching gene IDs */
    var targetGene = gene_list[targetGeneID];
    var genePos = targetGene.start + (targetGene.stop - targetGene.start) / 2;
    //var windowCenter = fabric.util.transformPoint({x:canvas.getWidth()/2,y:0}, canvas.viewportTransform)['x'];
    var windowCenter = canvas.getWidth()/2 - canvas.viewportTransform[4]; // canvas.getWidth()/2 is clientX of center of screen
    offsetX = windowCenter - genePos;
    console.log('offsetX: ' + offsetX + ', genePos: ' + genePos + ', windowCenter: ' + windowCenter);
  }

  // label
  if(showLabels) {
    canvas.add(new fabric.Text(label, {top: spacing*y-10, selectable: false, fontSize: 15, fontFamily: 'sans-serif', fontWeight: 'bold'}));
  }

  // line
  canvas.add(new fabric.Line([0,0,genomeMax,0], {left: offsetX+(showLabels?120:0),
        top: spacing*y - 4,
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
    var geneObj = geneArrow(gene,geneID,genomeData.genomes[genomeID].genes.functions[geneID],y,genomeID,arrowStyle);
    if(showLabels) {
      geneObj.left += 120;
    }
    if(alignToGC) {
      geneObj.left += offsetX;
    }
    canvas.add(geneObj);

    if(showGeneLabels) {
      var label = new fabric.IText("geneID: "+geneID, {
        hasControls:false,
        lockMovementX: true,
        lockMovementY: true,
        lockScaling: true,
        hoverCursor:'text'
      });

      if(arrowStyle == 3) {
        label.set({
          top: -10+spacing*y,
          left: 150+gene.start,
          scaleX: 0.25,
          scaleY: 0.25,
          selectionColor:'rgba(128,128,128,.5)'
        });
      } else {
        label.set({
          scaleX: 0.5,
          scaleY: 0.5,
          top: -30+spacing*y,
          left: 200+gene.start,
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

function geneArrow(gene, geneID, functions, y, genomeID, style) {
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

  var length = gene.stop-gene.start;

  var arrowPathStr;
  switch(style) {
    case 2: // thicker arrows
      arrowPathStr = 'M 0 -5 L ' + length + ' -5 L ' + length + ' 15 L 0 15 M ' + length + ' -5 L ' + length + ' 20 L ' + (25+length) + ' 5 L ' + length + ' -10 z';
      break;
    case 3: // pentagon arrows
      arrowPathStr = 'M 0 0 L ' + (length-25) + ' 0 L ' + length + ' 20 L ' + (length-25) + ' 40 L 0 40 L 0 0 z';
      break;
    default: // 'inspect page' arrows
      arrowPathStr = 'M 0 0 L ' + length + ' 0 L ' + length + ' 10 L 0 10 M ' + length + ' 0 L ' + length + ' 20 L ' + (25+length) + ' 5 L ' + length + ' -10 z';
      break;
  }

  var arrow = new fabric.Path(arrowPathStr);
  arrow.set({
    id: 'arrow',
    selectable: false,
    gene: gene,
    geneID: geneID,
    genomeID: genomeID,
    top: style == 3 ? -14+spacing*y : -11+spacing*y,
    left: 1.5+gene.start,
    scaleX: 0.5,
    scaleY: 0.5,
    fill: color,
    stroke: 'gray',
    strokeWidth: style == 3 ? 3 : 1.5,
    zoomX: 0.2,
    zoomY: 0.2
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

function resetScale(){
  canvas.setZoom(1)
}

function drawGenomeLabels(label){
  genomeLabelsCanvas.add(new fabric.Text(label, {
    top : labelSpacing,
    fontSize : 15,
    fontFamily: 'sans-serif',
    fontWeight: 'bold',
    selectable : true,
    hasControls : false,
    lockMovementX : true
  }))
}

function changeGenomeOrder(event){

  if(event.target.value === "Alphabetical"){
    console.log(typeof genomeData.genomes)

  }
}