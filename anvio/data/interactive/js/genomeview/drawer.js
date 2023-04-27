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
 * File Overview : The Drawer class defined here is responsible for rendering genomic + associated data passed from main.js to an interactive
 * browser canvas. This is where most of the heavy lifting should happen, and where most of our Fabric.js library interactions should occur.
 */
var GenomeDrawer = function (settings) {
  this.settings = settings;
};

GenomeDrawer.prototype.draw = function () {
  canvas.clear()
  canvas.setHeight(calculateMainCanvasHeight()) // set canvas height dynamically
  labelCanvas.setHeight(calculateMainCanvasHeight());
console.log("canvas height set")
  let orderIndex = 0
  this.settings['genomeData']['genomes'].map((genome, idx) => {
    console.log("inside addLayers loop iteration: " + idx)
    if($('#' + genome[0] + '-show').is(':checked')) {
      this.addLayers(idx, orderIndex)
      orderIndex++
    }
  })

  checkGeneLabels();
  drawTestShades();
  $('.loading-screen').hide();

  if(firstDraw){
    this.setInitialZoom()
    firstDraw = false
  }
}

/*
 *  For each genome group, iterate additional all layers and render where appropriate
 */
GenomeDrawer.prototype.addLayers = function (genomeIndex, orderIndex) {
  console.log("inside addLayers orderIndex" + orderIndex)
  let [dataLayerHeight, rulerHeight] = [this.calculateLayerSizes()[0], this.calculateLayerSizes()[1]]

  yOffset = orderIndex * spacing + (orderIndex * maxGroupSize * groupLayerPadding);
  let layerPos = 0
  let genomeID = this.settings['genomeData']['genomes'][genomeIndex][0];
  // let genome = this.settings['genomeData']['genomes'][genomeIndex][1];
  // let label = genome.genes.gene_calls[0].contig;

  let additionalDataLayers = this.settings['additional-data-layers']['data'][genomeID]

  let ptInterval = Math.floor(globalGenomeMax / settings['display']['adlPtsPerLayer']);

  this.settings['group-layer-order'].map((layer, idx) => {  // render out layers, ordered via group-layer-order array
    console.log("groupLayerOrder index "+ idx)
    if (layer == 'Genome' && $('#Genome-show').is(':checked')) {
      this.addGenome(orderIndex, dataLayerHeight, layerPos, genomeIndex)
      layerPos += dataLayerHeight + groupLayerPadding
      console.log("genome added successfully")
    }
    if (layer == 'Coverage' && this.settings['additional-data-layers']['layers'].includes('Coverage') && $('#Coverage-show').is(':checked')) {
      this.buildNumericalDataLayer('Coverage', layerPos, genomeID, additionalDataLayers, ptInterval, 'blue', dataLayerHeight, orderIndex, genomeIndex)
      layerPos += dataLayerHeight + groupLayerPadding
      console.log("coverage added successfully")
    }
    if (layer == 'GC_Content' && this.settings['additional-data-layers']['layers'].includes('GC_content') && $('#GC_Content-show').is(':checked')) {
      this.buildNumericalDataLayer('GC_content', layerPos, genomeID, additionalDataLayers, ptInterval, 'purple', dataLayerHeight, orderIndex, genomeIndex)
      layerPos += dataLayerHeight + groupLayerPadding
      console.log("GC Content added successfully")
    }
    if (layer == 'Ruler' && this.settings['additional-data-layers']['layers'].includes('ruler') && $('#Ruler-show').is(':checked')) {
      this.buildGroupRulerLayer(genomeID, layerPos, rulerHeight, orderIndex)
      layerPos += rulerHeight + groupLayerPadding
      console.log("Ruler added successfully")
    }
  })
console.log("groupLayer loop done")
  canvas.remove(canvas.getObjects().find(obj => obj.id == 'groupBorder' + genomeIndex));
  this.addGroupBorder(yOffset, orderIndex)
}

/*
 *  add a stylish and visually significant border around each group
 */
GenomeDrawer.prototype.addGroupBorder = function (yOffset, orderIndex) {

  if(maxGroupSize == 1) return;

  let top = yOffset + marginTop - 10 + (orderIndex * groupMargin)
  let left = -20
  let width = genomeMax[genomeID]*scaleFactor + 50
  let height = spacing + 50

  let rect = new fabric.Rect({
    id: 'groupBorder' + orderIndex,
    top: top,
    left: left,
    width: width,
    height: height,
    stroke: 'black',
    strokeWidth: 0.5,
    fill: '',
    selectable: false,
    hoverCursor: 'default'
  })

  canvas.sendToBack(rect)
}

/*
 *  programmatically calculate layer height values, given that the ruler layer should be allocated comparatively less space
 */
GenomeDrawer.prototype.calculateLayerSizes = function () {
  let parityHeight = spacing / maxGroupSize
  let rulerHeight = settings['display']['layers']['Ruler'] ? Math.floor(parityHeight * .5) : 0 // some arbitrary percentage of parity since ruler should get less y-axis space

  // with the extra space carved out by a smaller rulerHeight, distribute the excess evenly amongst all layers that are NOT rulers
  let dataLayerHeight = Math.floor(parityHeight + rulerHeight / (maxGroupSize - 1))
  if(maxGroupSize == 1) dataLayerHeight = parityHeight;
  return [dataLayerHeight, rulerHeight]
}

GenomeDrawer.prototype.addGenome = function (orderIndex, layerHeight, layerPos, genomeIndex) {
  console.log("addGenome function")
  let genome = this.settings['genomeData']['genomes'][genomeIndex];
  let gene_list = genome[1].genes.gene_calls;
  let genomeID = genome[0];
  let y = marginTop + yOffset + layerPos + (layerHeight / 2) + (orderIndex * groupMargin) // render arrows in the center of genome layer's allotted vertical space

  let [start, stop] = percentScale ? getRenderXRangeForFrac() : renderWindow.map(x => x * scaleFactor + xDisps[genomeID]);
  start = clamp(start > xDisps[genomeID] ? start : xDisps[genomeID], calcXBoundsForGenome(genomeID)[0], calcXBoundsForGenome(genomeID)[1]);
  stop = clamp(stop, calcXBoundsForGenome(genomeID)[0], calcXBoundsForGenome(genomeID)[1]);

  // line
  let lineObj = new fabric.Line([start, 0, stop, 0], {
    id: 'genomeLine',
    groupID: genomeID,
    top: y + 4,
    stroke: 'black',
    strokeWidth: 2,
    selectable: false,
    lockMovementY: true,
    hasControls: false,
    hasBorders: false,
    lockScaling: true,
    hoverCursor: 'default'
  });
  canvas.add(lineObj);
  //this.addBackgroundShade((marginTop + yOffset + layerPos + (orderIndex * groupMargin)), start, genomeMax[genomeID], layerHeight, orderIndex, genomeIndex)

  // draw set labels
  if (settings['display']['show-gene-labels'] && settings['display']['labels']['gene-sets'][genomeID]) {
    settings['display']['labels']['gene-sets'][genomeID].forEach(obj => {
      drawSetLabel(obj[0], genomeID, obj[1]);
    });
  }

  for (let geneID in gene_list) {
    console.log("inside for geneID in gene_list geneID: " + geneID)
    let gene = gene_list[geneID];
    let [ntStart, ntStop] = getRenderNTRange(genomeID);
    if (gene.start < ntStart) continue;
    if (gene.stop > ntStop) return;
    var geneObj = this.geneArrow(gene, geneID, y, genomeID, this.settings['display']['arrow-style']);
    canvas.add(geneObj)
    canvas.bringToFront(geneObj);

    if (settings['display']['show-gene-labels'] && ($('#brush_end').val() - $('#brush_start').val()) < 10000) {
      var label = new fabric.Text(setGeneLabelFromSource(geneID, genomeID), {
        id: 'geneLabel',
        groupID: genomeID,
        fontSize: settings['display']['gene-label-size'],
        fontFamily: 'sans-serif',
        angle: settings['display']['gene-text-position'] == "above" ? -1 * settings['display']['gene-text-angle'] : 0,
        //left: xDisps[genomeID] + (((gene.start + gene.stop) / 2) - $('#gene_label').val()*2) * scaleFactor,
        fill: $('#gene_label_color').attr('color'),
        scaleX: 0.5,
        scaleY: 0.5,
        editable: true,
        hasControls: false,
        opacity: settings['display']['hidden']?.[genomeID]?.[geneID] ? .2 : 1.0,
        lockMovementX: true,
        lockMovementY: true,
        lockScaling: true,
        hoverCursor: 'text'
      });
      
      label.left = xDisps[genomeID] + (((gene.start + gene.stop) / 2) - label.width/2) * scaleFactor;

      if (this.settings['display']['arrow-style'] == 3) {
        label.set({
          top: settings['display']['gene-text-position'] == "inside" ? y + 13 - settings['display']['gene-label-size'] / 2 : y - 20 - settings['display']['gene-label-size'] / 2,
          selectionColor: 'rgba(128,128,128,.5)'
        });
      } else {
        label.set({
          top: y - 10 - settings['display']['gene-label-size'] / 2,
          selectionColor: 'rgba(128,128,128,.2)'
        });
      }
      label.on("editing:exited", function (e) {
        console.log(label.text)
      });

      canvas.add(label);
    }

    function setGeneLabelFromSource(geneID, genomeID) {
      let genomeOfInterest = this.settings['genomeData']['genomes'].filter(genome => genome[0] == genomeID)
      let source = $('#gene_label_source').val()

      if (source == 'default') {
        return `${geneID}`
      }
      if (source == 'user') {
        if (this.settings['display']?.hasOwnProperty('gene-labels')) {
          return this.settings['display']['gene-labels'][genomeID][geneID]
        } else {
          return ''
        }
      } else {
        // operating under the assumption that
        // A) the relevant source value lives at index 1 of the source array
        // B) the selected value from the #gene_label_source dropdown === the source object itself (and it should, since we build the dropdown programmatically!)
        // this below logic should retrieve the correct annotation value, regardless of source
        if (genomeOfInterest[0][1]['genes']['functions'][geneID]?.hasOwnProperty(source) && genomeOfInterest[0][1]['genes']['functions'][geneID][source]) {
          return ellipsisMachine(genomeOfInterest[0][1]['genes']['functions'][geneID][source][1])
        } else {
          return 'None'
        }
      }
    }

    function ellipsisMachine(string) { // add ellipsis only to truncated gene label values
      if (string.substring(0, 20).length == 20) {
        return `${string.substring(0, 20)}...`
      } else {
        return string
      }
    }
  }

  function drawSetLabel(title, genomeID, geneIDs) {
    // assume gene IDs form contiguous list
    let geneObjs = geneIDs.map(geneID => settings['genomeData']['genomes'].find(obj => obj[0] == genomeID)[1].genes.gene_calls[geneID]);
    let x_set_label = geneObjs[0].start + (geneObjs[geneObjs.length - 1].stop - geneObjs[0].start) / 2;
    let y_set_label = y - 10 - settings['display']['gene-label-size'];
    var set_label = new fabric.IText(title, {
      id: 'setLabel',
      groupID: genomeID,
      fontSize: settings['display']['gene-label-size'],
      left: xDisps[genomeID] + x_set_label * scaleFactor,
      top: y_set_label,
      scaleX: 0.5,
      scaleY: 0.5,
      editable: true,
      hasControls: false,
      lockMovementX: true,
      lockMovementY: true,
      lockScaling: true,
      hoverCursor: 'text',
    });
    canvas.add(set_label);
  }
}

/*
 *  Process to generate numerical ADL for genome groups (ie Coverage, GC Content )
 */
GenomeDrawer.prototype.buildNumericalDataLayer = function (layer, layerPos, genomeID, additionalDataLayers, ptInterval, defaultColor, layerHeight, orderIndex, genomeIndex) {
  // TODO this will need to be refactored once we begin testing genomes comprised of multiple contigs
  let contigObj = Object.values(additionalDataLayers)[0]
  let contigArr = Object.values(contigObj)[0]
  let stroke = 'black'

  // if(layer == 'Coverage'){
  //   this.settings['display']['additional-data-layers']['coverage'] ? stroke = this.settings['display']['additional-data-layers']['coverage'] : stroke = 'black'
  // }
  if (layer == 'GC_content') { // we will need to refactor and get our variable case/formatting nonsense sorted.
    this.settings['display']['colors']['GC_Content'] ? stroke = this.settings['display']['colors']['GC_Content'] : stroke = 'red'
  }

  let maxDataLayerValue = 0
  let startingTop = marginTop + yOffset + layerPos + (orderIndex * groupMargin)
  let startingLeft = xDisps[genomeID]

  let globalPathDirective = [`L ${startingLeft} ${layerHeight}`]
  let layer_end_final_coordinates

  let pathDirective = [`M ${startingLeft} 0 L ${startingLeft} ${layerHeight}`]

  for (let i = 0; i < contigArr.length; i++) {
    contigArr[i] > maxDataLayerValue ? maxDataLayerValue = contigArr[i] : null
  }

  let nGroups = 20
  let j = 0
  let final_l = 0 //used to create final line segments to 'close out' path obj for shading purposes.
  let [l, r] = getRenderNTRange(genomeID);
  for (let i = 0; i <= nGroups; i++) {
    for (; j < i * contigArr.length / nGroups; j += ptInterval) {
      if (j < l) continue;
      if (j > r) break;

      let left = j * scaleFactor + startingLeft
      let top = layerHeight - ((contigArr[j] / maxDataLayerValue) * layerHeight)
      let segment = `L ${left} ${top}`
      final_l = left // final_l is always last-seen x coordinate
      pathDirective.push(segment)
      globalPathDirective.push(segment)
    }
    // TODO resolve performance-related aspects of the chunking done below

    // let graphObj = new fabric.Path(pathDirective.join(' '))
    // graphObj.set({
    //   top : startingTop,
    //   stroke : stroke,
    //   fill : '',
    //   selectable: false,
    //   objectCaching: false,
    //   id : `${layer} graph`,
    //   groupID : genomeID,
    //   genome : genomeID
    // })
    // canvas.bringToFront(graphObj)
    pathDirective = []
  }
  layer_end_final_coordinates = `L ${final_l} ${layerHeight} L ${startingLeft} ${layerHeight} z`
  globalPathDirective.push(layer_end_final_coordinates)

  let shadedObj = new fabric.Path(globalPathDirective.join(' '))
  shadedObj.set({
    top: startingTop,
    stroke: stroke,
    fill: stroke,
    selectable: false,
    objectCaching: false,
    hoverCursor: 'default',
    id: `${layer}-graph-shaded`,
    groupID: genomeID,
    genome: genomeID
  })
  canvas.bringToFront(shadedObj)
  this.addBackgroundShade(startingTop, startingLeft, genomeMax[genomeID]*scaleFactor, layerHeight, orderIndex, genomeIndex)
}

/*
 *  Generate individual genome group rulers
 */
GenomeDrawer.prototype.buildGroupRulerLayer = function (genomeID, layerPos, layerHeight, orderIndex) {
  let startingTop = marginTop + yOffset + layerPos + (orderIndex * groupMargin)
  let startingLeft = xDisps[genomeID]
  // let layerHeight = (spacing / maxGroupSize)

  // split ruler into several objects to avoid performance cost of large object pixel size
  let nRulers = 20;
  let w = 0;
  let [l, r] = getRenderNTRange(genomeID);
  for (let i = 0; i < nRulers; i++) {
    let ruler = new fabric.Group();
    for (; w < (i + 1) * genomeMax[genomeID] / nRulers; w += settings['display']['genome-scale-interval']) {
      if (w < l) continue;
      if (w > r) break;
      let tick = new fabric.Line([0, 0, 0, 20], {
        left: (w * scaleFactor),
        stroke: 'black',
        strokeWidth: 1,
        fontSize: 10,
        fontFamily: 'sans-serif'
      });
      let lbl = new fabric.Text(w / 1000 + " kB", {
        left: (w * scaleFactor + 5),
        stroke: 'black',
        strokeWidth: .25,
        fontSize: 15,
        fontFamily: 'sans-serif'
      });
      ruler.add(tick);
      ruler.add(lbl);
    }
    ruler.set({
      left: startingLeft,
      top: startingTop + (layerHeight / 2),
      lockMovementY: true,
      hasControls: false,
      hasBorders: false,
      lockScaling: true,
      objectCaching: false,
      selectable: false,
      hoverCursor: 'default',
      groupID: genomeID,
      class: 'ruler'
    });
    ruler.addWithUpdate();
    canvas.add(ruler);
  }
  //this.addBackgroundShade(startingTop, startingLeft, genomeMax[genomeID], layerHeight+5, orderIndex, genomeIndex)
}

/*
 *  adds an alternating shade to each genome group for easier visual distinction amongst adjacent groups
 */
GenomeDrawer.prototype.addBackgroundShade = function (top, left, width, height, orderIndex, genomeIndex) {
  let backgroundShade;
  orderIndex % 2 == 0 ? backgroundShade = '#b8b8b8' : backgroundShade = '#f5f5f5'

  let border = new fabric.Rect({
    groupID: this.settings['genomeData']['genomes'][genomeIndex][0],
    top: top,
    left: left,
    width: width,
    height: height,
    stroke: 'black',
    strokeWidth: 2,
    fill: "rgba(0,0,0,0.0)",
    selectable: false,
    // opacity : .5
  });
  let background = new fabric.Rect({
    groupID: this.settings['genomeData']['genomes'][genomeIndex][0],
    top: top,
    left: left,
    width: width,
    height: height,
    fill: "#b0b0b0",
    stroke: 'black',
    selectable: false,
    hoverCursor: 'default',
    opacity: .25
  });
  // canvas.sendToBack(border)
  canvas.sendToBack(background)
}

GenomeDrawer.prototype.geneArrow = function (gene, geneID, y, genomeID, style) {
  let ind = this.settings['genomeData']['genomes'].findIndex(g => g[0] == genomeID);
  let functions = this.settings['genomeData']['genomes'][ind][1].genes.functions[geneID];
  let color;

  if($('#user_defined_colors').is(':checked')) {
    color = '#808080' // default = grey

    // check if gene is highlighted
    // possibly deprecated - requires new table for user-defined colors
    let pickerCode = genomeID + '-' + geneID;
    if ($('#picker_' + pickerCode).length > 0) {
      console.log(pickerCode)
      color = $('#picker_' + pickerCode).attr('color');
    }

    // check for set colors
    if (settings['display']['colors']['genes'][genomeID] && settings['display']['colors']['genes'][genomeID][geneID]) {
      color = settings['display']['colors']['genes'][genomeID][geneID];
    }
  } else {
    color = $('#picker_None').length > 0 ? $('#picker_None').attr('color') : 'gray';
    let cag = getCagForType(functions, color_db);
    if (cag) {
      cag = getCleanCagCode(cag);
      let color_other = $('#picker_Other').length > 0 ? $('#picker_Other').attr('color') : 'white';
      color = $('#picker_' + cag).length > 0 ? $('#picker_' + cag).attr('color') : color_other;
    } else {
      if (gene.source.startsWith('Ribosomal_RNA')) {
        cag = 'rRNA';
      } else if (gene.source == 'Transfer_RNAs') {
        cag = 'tRNA';
      } else if (functions != null) {
        cag = 'Function';
      }
      if ($('#picker_' + cag).length > 0) color = $('#picker_' + cag).attr('color');
    }
  }

  let length = (gene.stop - gene.start) * scaleFactor;
  let stemLength = length - 25 > 0 ? length - 25 : 0;

  var arrowPathStr;
  switch (parseInt(style)) {
    case 2: // thicker arrows
      arrowPathStr = `M ${stemLength} -10
                      L 0 -10
                      L 0 20
                      L ${stemLength} 20
                      L ${stemLength} 20
                      L ${stemLength} 20
                      L ${length} 5
                      L ${stemLength} -10 z`;
      break;
    case 3: // pentagon arrows
      arrowPathStr = `M 0 0
                      L ${stemLength} 0
                      L ${length} 20
                      L ${stemLength} 40
                      L 0 40
                      L 0 0 z`;
      break;
    case 4: // rect arrows
      arrowPathStr = `M ${length} -5
                      L 0 -5
                      L 0 15
                      L ${length} 15 z`;
      break;
    default: // 'inspect page' arrows
      arrowPathStr = `M ${stemLength} -2.5
                      L 0 -2.5
                      L 0 12.5
                      L ${stemLength} 12.5
                      L ${stemLength} 20
                      L ${length} 5
                      L ${stemLength} -10 z`;
      break;
  }

  var arrow = new fabric.Path(arrowPathStr);
  let genomeOfInterest = settings['genomeData']['genomes'].filter(genome => genome[0] == genomeID)
  arrow.set({
    id: 'arrow',
    groupID: genomeID,
    lockMovementX: true,
    lockMovementY: true,
    selectable: true,
    hasControls: false,
    hasBorders: false,
    lockScaling: true,
    hoverCursor: 'pointer',
    opacity: settings['display']['hidden']?.[genomeID]?.[geneID] ? .2 : 1.0,
    aaSequence: genomeOfInterest[0][1]['genes']['aa'][geneID]['sequence'],
    dnaSequence: genomeOfInterest[0][1]['genes']['dna'][geneID]['sequence'],
    gene: gene,
    functions: functions,
    geneID: geneID,
    genomeID: genomeID,
    top: style == 3 ? y - 17 : y - 11, // TODO update this offset to reflect genome layer height (we want to render this arrow in the middle of its allocated height)
    left: xDisps[genomeID] + (1.5 + gene.start) * scaleFactor,
    fill: color,
    stroke: 'gray',
    strokeWidth: 1.5
  });
  if (gene.direction == 'r') arrow.rotate(180);

  return arrow;
}

/*
 *  Draw background shades between genes of the same cluster.
 *  TODO: generalize this function to take in [start,stop,val] NT sequence ranges to shade any arbitrary metric
 *        - add a separate function for retrieving [start,stop,val] given gene cluster IDs
 *
 *  @param geneClusters : array of GC IDs to be shaded
 *  @param colors       : dict defining color of each shade, in the form {geneClusterID : hexColor}
 */
GenomeDrawer.prototype.shadeGeneClusters = function (geneClusters, colors) {
  if (!genomeData.gene_associations["anvio-pangenome"]) return;

  let y = marginTop;
  for (var i = 0; i < genomeData.genomes.length - 1; i++) {
    let genomeA = this.settings['genomeData']['genomes'][i][1].genes.gene_calls;
    let genomeB = this.settings['genomeData']['genomes'][i + 1][1].genes.gene_calls;
    let genomeID_A = this.settings['genomeData']['genomes'][i][0];
    let genomeID_B = this.settings['genomeData']['genomes'][i + 1][0];
    let [l1, r1] = getRenderNTRange(genomeID_A);
    let [l2, r2] = getRenderNTRange(genomeID_B);

    for (gc of geneClusters) {
      let g1 = [], g2 = [];
      for (geneID of genomeData.gene_associations["anvio-pangenome"]["gene-cluster-name-to-genomes-and-genes"][gc][genomeID_A]) {
        g1.push(genomeA[geneID].start, genomeA[geneID].stop);
      }
      for (geneID of genomeData.gene_associations["anvio-pangenome"]["gene-cluster-name-to-genomes-and-genes"][gc][genomeID_B]) {
        g2.push(genomeB[geneID].start, genomeB[geneID].stop);
      }

      // if shades outside render bounds, don't draw them
      if (g1[1] < l1 && g2[1] < l2) continue;
      if (g1[0] > r1 && g2[0] > r2) break;

      g1 = g1.map(val => val * scaleFactor + xDisps[genomeID_A]);
      g2 = g2.map(val => val * scaleFactor + xDisps[genomeID_B]);

      /* TODO: implementation for multiple genes of the same genome in the same gene cluster */
      var path = new fabric.Path("M " + g1[0] + " " + y + " L " + g1[1] + " " + y + " L " + g2[1] + " " + (y + spacing) + " L " + g2[0] + " " + (y + spacing) + " z", {
        id: 'link',
        fill: colors[gc],
        opacity: 0.25,
        selectable: false
      });
      path.sendBackwards();
      canvas.sendBackwards(path);
    }
    y += spacing
  }
}

/*
 *  Add a temporary glow effect around given gene(s).
 *
 *  @param geneParams : array of dicts, in one of two formats:
 *    (1) [{genomeID: gid_1, geneID: [id_1, id_2, ...]}, ...]
 *    (2) [{genomeID: gid_1, geneID: id_1}, ...]
 */
GenomeDrawer.prototype.glowGenes = function (geneParams, indefinite=false, timeInterval=5000) {
  // convert geneParams format (1) to format (2)
  if (Array.isArray(geneParams[0].geneID)) {
    let newParams = [];
    for (genome of geneParams) {
      for (gene of genome.geneID) newParams.push({ genomeID: genome.genomeID, geneID: gene });
    }
    geneParams = newParams;
  }

  this.removeAllGeneGlows();

  var shadow = new fabric.Shadow({
    color: 'red',
    blur: 30
  });
  var arrows = canvas.getObjects().filter(obj => obj.id == 'arrow' && geneParams.some(g => g.genomeID == obj.genomeID && g.geneID == obj.geneID));

  if(indefinite) {
    for(arrow of arrows) {
      shadow.blur = 20;
      arrow.set('stroke', 'black');
      arrow.set('shadow', shadow);
    }
    canvas.renderAll();
    return;
  }

  for (arrow of arrows) {
    arrow.set('shadow', shadow);
    arrow.animate('shadow.blur', 0, {
      duration: timeInterval,
      onChange: canvas.renderAll.bind(canvas),
      onComplete: function () {},
      easing: fabric.util.ease['easeInQuad']
    });
  }
  canvas.renderAll();
}

/*
 *  Removes all active gene glow effects.
 */
GenomeDrawer.prototype.removeAllGeneGlows = function () {
  canvas.getObjects().filter(obj => obj.id == 'arrow' && obj.shadow).forEach(arrow => {
    arrow.set('stroke', 'gray');
    delete arrow.shadow;
  });
  canvas.renderAll();
}

/*
 *  Shift genomes horizontally to align genes around the target gene cluster.
 *
 *  @param gc : target gene cluster ID
 */
GenomeDrawer.prototype.alignToCluster = function (gc) {
  if (!this.settings['genomeData']['gene_associations']["anvio-pangenome"]) return;

  let targetGeneInfo = viewCluster(gc);
  if (targetGeneInfo == null) return;
  let [firstGenomeID, targetGeneMid] = targetGeneInfo;
  if (firstGenomeID != null) {
    alignToGC = gc;
    let index = this.settings['genomeData']['genomes'].findIndex(g => g[0] == firstGenomeID);
    for (var i = index + 1; i < this.settings['genomeData']['genomes'].length; i++) {
      let gid = genomeData.genomes[i][0];
      let geneMids = getGenePosForGenome(genomeData.genomes[i][0], alignToGC);
      if (geneMids == null) continue;
      let geneMid = geneMids[0]; /* TODO: implementation for multiple matching gene IDs */
      let shift = scaleFactor * (targetGeneMid - geneMid) + (xDisps[firstGenomeID] - xDisps[gid]);
      let objs = canvas.getObjects().filter(obj => obj.groupID == gid);
      for (o of objs) o.left += shift;
      xDisps[gid] += shift;
      canvas.setViewportTransform(canvas.viewportTransform);

      // clear and redraw shades
      clearShades();
      drawTestShades();
    }
  }
}

/*
 *  Clear all gene links from the canvas.
 */
GenomeDrawer.prototype.clearShades = function () {
  canvas.getObjects().filter(obj => obj.id == 'link').forEach((l) => { canvas.remove(l) });
}

GenomeDrawer.prototype.setPtsPerADL = function (newResolution) {
  if (isNaN(newResolution)) return;
  newResolution = parseInt(newResolution);
  if (newResolution < 0 || newResolution > globalGenomeMax) {
    alert(`Invalid value, points per additional data layer must be in range 0-${globalGenomeMax}.`);
    return;
  }
  settings['display']['adlPtsPerLayer'] = newResolution;
  this.draw();
}

GenomeDrawer.prototype.showAllADLPts = function () {
  this.setPtsPerADL(globalGenomeMax); // TODO: account for different genome maxes
  $('#adl_pts_per_layer').val(globalGenomeMax);
  $('#showAllADLPtsBtn').blur();
}

GenomeDrawer.prototype.alignRulers = function () {
  for (genome of this.settings['genomeData']['genomes']) {
    xDisps[genome[0]] = 0;
  }
  percentScale = false;
  drawScale();
  bindViewportToWindow();
  updateScalePos();
  updateRenderWindow();
  this.draw();
  $('#alignRulerBtn').blur();
}

GenomeDrawer.prototype.setGenomeSpacing = function (newSpacing) {
  if (isNaN(newSpacing)) return;
  newSpacing = parseInt(newSpacing);
  if (newSpacing < 0 || newSpacing > 1000) {
    alert(`Invalid value, genome spacing must be in range 0-1000.`);
    return;
  }
  spacing = newSpacing;
  this.draw();
  drawGenomeLabels(); // update label font size
  setLabelCanvas(); // update label canvas width
}

GenomeDrawer.prototype.setGenomeMargin = function (newMargin) {
  if (isNaN(newMargin)) return;
  newMargin = parseInt(newMargin);
  if (newMargin < 0 || newMargin > 1000) {
    alert(`Invalid value, genome margin must be in range 0-1000.`);
    return;
  }
  groupMargin = newMargin;
  this.draw();
  drawGenomeLabels(); // update label font size
}

GenomeDrawer.prototype.setScaleInterval = function (newScale) {
  if (isNaN(newScale)) return;
  newScale = parseInt(newScale);
  if (newScale < 50) {
    alert(`Invalid value, scale interval must be >=50.`);
    return;
  }
  settings['display']['genome-scale-interval'] = newScale;
  this.draw();
}

GenomeDrawer.prototype.setGeneLabelSize = function (newSize) {
  if (isNaN(newSize)) return;
  newSize = parseInt(newSize);
  if (newSize < 0 || newSize > 1000) {
    alert(`Invalid value, gene label size must be in range 0-1000.`);
    return;
  }
  settings['display']['gene-label-size'] = newSize;
  if (settings['display']['show-gene-labels']) this.draw();
}

GenomeDrawer.prototype.setGenomeLabelSize = function (newSize) {
  if (isNaN(newSize)) return;
  newSize = parseInt(newSize);
  if (newSize < 0 || newSize > 1000) {
    alert(`Invalid value, genome label size must be in range 0-1000.`);
    return;
  }
  settings['display']['genome-label-size'] = newSize;
  if (settings['display']['show-genome-labels']) this.draw();
}

GenomeDrawer.prototype.redrawSingleGenome = function (genomeID) {
  canvas.getObjects().filter(o => o.groupID == genomeID).forEach(obj => canvas.remove(obj));
  let idx = this.settings['genomeData']['genomes'].findIndex(obj => obj[0] == genomeID);
  this.addLayers(idx);
  checkGeneLabels();
}

/*
 *  Dynamically set scale tick interval based on scaleFactor.
 */
GenomeDrawer.prototype.adjustScaleInterval = function () {
  let val = Math.floor(100 / scaleFactor);
  let roundToDigits = Math.floor(Math.log10(val)) - 1;
  let newInterval = Math.floor(val / (10 ** roundToDigits)) * (10 ** roundToDigits);
  settings['display']['genome-scale-interval'] = newInterval;
  $('#genome_scale_interval').val(settings['display']['genome-scale-interval']);
}

GenomeDrawer.prototype.queryFunctions = async function () {
  $('#query-results-table').empty()
  $('#query-results-span').empty()
  let query = $('#function_search_query').val().toLowerCase()
  let category = $('#function_search_category').val()
  let distinctQueryMatches = Object()
  let glowPayload = Array()
  let foundInGenomes = Object()

  if (!query || !category) {
    alert('please provide values for function category and/or query')
    return
  }
  if (category == 'metadata'){
    drawer.queryMetadata(query)
    return
  }
  this.settings['genomeData']['genomes'].map(genome => {
    for (const [key, value] of Object.entries(genome[1]['genes']['functions'])) {
      if (value[category]?.[0].toLowerCase().includes(query)) {
        let glowObject = {
          genomeID: genome[0],
          geneID: key,
          matchedQuery: value[category][0]
        }
        glowPayload.push(glowObject)
        if (!(genome[0] in foundInGenomes)) {
          foundInGenomes[genome[0]] = true
        }
        if(!(value[category][0].toLowerCase() in distinctQueryMatches)){
          distinctQueryMatches[value[category][0]] = true
        }
      }  // check for accession and annotation values separately, as we want to capture the exact match for sorting results
      else if (value[category]?.[1].toLowerCase().includes(query)) {
        let glowObject = {
          genomeID: genome[0],
          geneID: key,
          matchedQuery: value[category][1]
        }
        glowPayload.push(glowObject)
        if (!(genome[0] in foundInGenomes)) {
          foundInGenomes[genome[0]] = true
        }
        if(!(value[category][1].toLowerCase() in distinctQueryMatches)){
          distinctQueryMatches[value[category][1]] = true
        }
      }
    }
  })

  if (glowPayload.length < 1) {
    alert(`No hits were found matching ${query} in ${category}`)
    return
  }
  $('#query-results-span').append(`<select class="form-control" id="query-results-select"></select>`)
  $('#query-results-select').append(new Option('Show All', 'all'))
  Object.keys(distinctQueryMatches).map(k => {
    $('#query-results-select').append(new Option(k.length > 80? k.slice(0,80) + '...' : k, k))
  })
  $('#query-results-select').on('change', function(){
    $('#query-results-table').empty()
    let matchedQuery = this.value
    if(matchedQuery == 'all'){
      renderAllGenes()
    } else {
      glowPayload.filter(gene => gene['matchedQuery'] == matchedQuery).map(gene => {
        let genomeOfInterest = settings['genomeData']['genomes'].filter(genome => genome[0] == gene['genomeID'])
        let start = genomeOfInterest[0][1]['genes']['gene_calls'][gene['geneID']]['start']
        let end = genomeOfInterest[0][1]['genes']['gene_calls'][gene['geneID']]['stop']
        $('#query-results-table').append(`
          <tr>
            <td>${gene['geneID']}</td>
            <td>${gene['genomeID']}</td>
            <td>${start}</td>
            <td>${end}</td>
            <td><button onclick="goToGene(${gene['genomeID']}[0]['id'],${gene['geneID']},${start},${end})">go to</button</td>
          </tr>
        `)
      })
    }
  })
  let lowestStart, highestEnd = null
  function renderAllGenes(){
    if(globalGenomeMax > 35000){
      glowPayload.map(gene => {
        let genomeOfInterest = this.settings['genomeData']['genomes'].filter(genome => genome[0] == gene['genomeID'])
        let start = genomeOfInterest[0][1]['genes']['gene_calls'][gene['geneID']]['start']
        let end = genomeOfInterest[0][1]['genes']['gene_calls'][gene['geneID']]['stop']
        if (start < lowestStart || lowestStart == null) lowestStart = start
        if (end > highestEnd || highestEnd == null) highestEnd = end
        $('#query-results-table').append(`
          <tr>
            <td>${gene['geneID']}</td>
            <td>${gene['genomeID']}</td>
            <td>${start}</td>
            <td>${end}</td>
            <td><button onclick="goToGene(${gene['genomeID']}[0]['id'],${gene['geneID']},${start},${end})">go to</button</td>
          </tr>
        `)
      })
    } else {
      lowestStart = 0
      highestEnd = globalGenomeMax
    }
  }
  renderAllGenes()

  $('#function-query-results-statement').text(`Retreived ${glowPayload.length} hit(s) from ${Object.keys(foundInGenomes).length} genomes`)
  await zoomOutAndWait('partial', lowestStart, highestEnd, 350)
  this.glowGenes(glowPayload, true)
}

GenomeDrawer.prototype.queryMetadata = async function(metadataLabel){
  $('#query-results-table').empty()
  let glowPayload = Array()
  let foundInGenomes = Object()
  let matches = settings['display']['metadata'].filter( m => m.label.toLowerCase().includes(metadataLabel.toLowerCase()))
                                               .filter( m => m.type == 'tag')
  matches.map(metadata => {
    glowPayload.push({
      geneID: metadata.gene,
      genomeID: metadata.genome
    })
    if (!(metadata.genome in foundInGenomes)) {
      foundInGenomes[metadata.genome] = true
    }
  })
  if (glowPayload.length < 1) {
    alert(`No hits were found matching ${metadataLabel} in metadata`)
    return
  }
  let lowestStart, highestEnd = null
  if (globalGenomeMax > 35000) { // TODO: instead of globalGenomeMax, use genomeMax[genomeID] for currently selected genomeID
    glowPayload.map(gene => {
      let genomeOfInterest = this.settings['genomeData']['genomes'].filter(genome => genome[0] == gene['genomeID'])
      let start = genomeOfInterest[0][1]['genes']['gene_calls'][gene['geneID']]['start']
      let end = genomeOfInterest[0][1]['genes']['gene_calls'][gene['geneID']]['stop']
      if (start < lowestStart || lowestStart == null) lowestStart = start
      if (end > highestEnd || highestEnd == null) highestEnd = end
      $('#query-results-table').append(`
        <tr>
          <td>${gene['geneID']}</td>
          <td>${gene['genomeID']}</td>
          <td>${start}</td>
          <td>${end}</td>
          <td><button onclick="goToGene(${gene['genomeID']}[0]['id'],${gene['geneID']},${start},${end})">go to</button</td>
        </tr>
      `)
    })
  } else {
    lowestStart = 0
    highestEnd = globalGenomeMax
  }
  await zoomOutAndWait('partial', lowestStart, highestEnd, 350)
  this.glowGenes(glowPayload, true)
}

GenomeDrawer.prototype.setInitialZoom = function(){
    let start = 0
    let stop = globalGenomeMax > 35000 ? 35000 : globalGenomeMax
    zoomOut('partial', start, stop)
}
