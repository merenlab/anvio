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

  let orderIndex = 0
  this.settings['genomeData']['genomes'].map((genome, idx) => {
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
  let [dataLayerHeight, rulerHeight] = [this.calculateLayerSizes()[0], this.calculateLayerSizes()[1]]

  yOffset = orderIndex * spacing + (orderIndex * maxGroupSize * groupLayerPadding);
  let layerPos = 0
  let genomeID = this.settings['genomeData']['genomes'][genomeIndex][0];
  // let genome = this.settings['genomeData']['genomes'][genomeIndex][1];
  // let label = genome.genes.gene_calls[0].contig;

  let additionalDataLayers = this.settings['additional-data-layers']['data'][genomeID]

  let ptInterval = globalGenomeMax > settings['display']['adlPtsPerLayer'] ? Math.floor(globalGenomeMax / settings['display']['adlPtsPerLayer']) : 1;

  this.settings['group-layer-order'].map((layer, idx) => {  // render out layers, ordered via group-layer-order array
    if (layer == 'Genome' && $('#Genome-show').is(':checked')) {
      this.addGenome(orderIndex, dataLayerHeight, layerPos, genomeIndex)
      layerPos += dataLayerHeight + groupLayerPadding
    }
    if (layer == 'Coverage' && this.settings['additional-data-layers']['layers'].includes('Coverage') && $('#Coverage-show').is(':checked')) {
      this.buildNumericalDataLayer('Coverage', layerPos, genomeID, additionalDataLayers, ptInterval, 'blue', dataLayerHeight, orderIndex, genomeIndex)
      layerPos += dataLayerHeight + groupLayerPadding
    }
    if (layer == 'GC_Content' && this.settings['additional-data-layers']['layers'].includes('GC_content') && $('#GC_Content-show').is(':checked')) {
      this.buildNumericalDataLayer('GC_content', layerPos, genomeID, additionalDataLayers, ptInterval, 'purple', dataLayerHeight, orderIndex, genomeIndex)
      layerPos += dataLayerHeight + groupLayerPadding
    }
    if (layer == 'Ruler' && this.settings['additional-data-layers']['layers'].includes('ruler') && $('#Ruler-show').is(':checked')) {
      this.buildGroupRulerLayer(genomeID, layerPos, rulerHeight, orderIndex)
      layerPos += rulerHeight + groupLayerPadding
    }
  })

  canvas.remove(canvas.getObjects().find(obj => obj.id == 'groupBorder' + genomeIndex));
  this.addGroupBorder(yOffset, orderIndex, genomeIndex)
}

/*
 *  add a stylish and visually significant border around each group
 */
GenomeDrawer.prototype.addGroupBorder = function (yOffset, orderIndex, genomeIndex) {

  if(maxGroupSize == 1) return;

  let genomeID = this.settings['genomeData']['genomes'][genomeIndex][0]

  let top = yOffset + marginTop - 10 + (orderIndex * groupMargin)
  let left = -20 + nt_disps[genomeID]*scaleFactor
  let width = genomeMax[genomeID]*scaleFactor + 50
  let height = spacing + 50

  let rect = new fabric.Rect({
    id: 'groupBorder' + orderIndex,
    groupID: genomeID,
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
  let genome = this.settings['genomeData']['genomes'][genomeIndex];
  let gene_list = genome[1].genes.gene_calls;
  let genomeID = genome[0];
  let y = marginTop + yOffset + layerPos + (layerHeight / 2) + (orderIndex * groupMargin) // render arrows in the center of genome layer's allotted vertical space

  // get x-range for this genome
  //let [startX, stopX] = percentScale ? getRenderXRangeForFrac() : getRenderXRangeForGenome(genomeID)
  let [startX, stopX] = renderWindow.map(pos => pos * scaleFactor)
  startX = clamp(startX > nt_disps[genomeID]*scaleFactor ? startX : nt_disps[genomeID]*scaleFactor, calcXBoundsForGenome(genomeID)[0], calcXBoundsForGenome(genomeID)[1]);
  stopX = clamp(stopX, calcXBoundsForGenome(genomeID)[0], calcXBoundsForGenome(genomeID)[1]);

  // genome line
  let lineObj = new fabric.Line([startX, 0, stopX, 0], {
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
  //this.addBackgroundShade((marginTop + yOffset + layerPos + (orderIndex * groupMargin)), startX, genomeMax[genomeID], layerHeight, orderIndex, genomeIndex)

  // draw set labels
  if (settings['display']['show-gene-labels'] && settings['display']['labels']['gene-sets'][genomeID]) {
    settings['display']['labels']['gene-sets'][genomeID].forEach(obj => {
      drawSetLabel(obj[0], genomeID, obj[1]);
    });
  }

  for (let geneID in gene_list) {
    // TODO: renderWindow needs to compensate for genomes slid left of 0...don't clamp to 0, end etc
    // rather than clamping to [0,max], clamp to [-PAD, max-PAD]
    let gene = gene_list[geneID];
    let [ntStart, ntStop] = getGenomeRenderWindow(genomeID);
    if (gene.stop < ntStart) continue;
    if (gene.start > ntStop) return;
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
        //left: nt_disps[genomeID] + (((gene.start + gene.stop) / 2) - $('#gene_label').val()*2) * scaleFactor,
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
      
      label.left = (scaleFactor * (nt_disps[genomeID] + (gene.start + gene.stop)/2)) - label.width/4;

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
      } else if (source == 'user') {
        if (this.settings['display']?.['metadata']?.filter(m => m.genome == genomeID && m.gene == geneID && m.type == 'annotation')?.length > 0) {
          return this.settings['display']['metadata'].filter(m => m.genome == genomeID && m.gene == geneID && m.type == 'annotation')[0].annotation;
        } else {
          return ''
        }
      } else {
        // operating under the assumption that
        // A) the relevant source value lives at index 1 of the source array
        // B) the selected value from the #gene_label_source dropdown === the source object itself (and it should, since we build the dropdown programmatically!)
        // this below logic should retrieve the correct annotation value, regardless of source
        if (genomeOfInterest[0][1]['genes']['functions'][geneID]?.hasOwnProperty(source) && genomeOfInterest[0][1]['genes']['functions'][geneID][source]) {
          let gene = genomeOfInterest[0][1]['genes']['gene_calls'][geneID]
          let geneWidth = scaleFactor * (gene.stop - gene.start)
          let charWidthEstimate = 0.3 * settings['display']['gene-label-size'];
          let charLimit = Math.floor(geneWidth / charWidthEstimate)
          if(charLimit < 5) return ""; 
          return ellipsisMachine(genomeOfInterest[0][1]['genes']['functions'][geneID][source][1], charLimit)
        } else {
          return 'None'
        }
      }
    }

    function ellipsisMachine(string, maxLen=20) { // add ellipsis only to truncated gene label values
      if (string.substring(0, maxLen).length == maxLen) {
        return `${string.substring(0, maxLen)}...`
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
      fontFamily: 'sans-serif',
      left: nt_disps[genomeID]*scaleFactor + x_set_label * scaleFactor,
      top: y_set_label+10,
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
  let startingLeft = nt_disps[genomeID]*scaleFactor

  let globalPathDirective = [`L ${startingLeft} ${layerHeight}`]
  let layer_end_final_coordinates

  let pathDirective = [`M ${startingLeft} 0 L ${startingLeft} ${layerHeight}`]

  for (let i = 0; i < contigArr.length; i++) {
    contigArr[i] > maxDataLayerValue ? maxDataLayerValue = contigArr[i] : null
  }
  
  if(ptInterval == 0) {
    console.log("error: ptInterval for numerical additional data layer was set to 0. if you are getting this error, ptInterval was likely calculated wrong due to some rounding error (perhaps you are using a short genome), and you should contact the developers.")
    return;
  }

  let nGroups = 20
  let j = 0
  let final_l = 0 //used to create final line segments to 'close out' path obj for shading purposes.
  let [l, r] = getGenomeRenderWindow(genomeID);
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
    lockMovementY: true,
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
  let startingLeft = nt_disps[genomeID]*scaleFactor
  // let layerHeight = (spacing / maxGroupSize)

  // split ruler into several objects to avoid performance cost of large object pixel size
  let nRulers = 20;
  let w = 0;
  let [l, r] = getGenomeRenderWindow(genomeID);
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
      lockMovementX: true,
      hasControls: false,
      hasBorders: false,
      lockScaling: true,
      objectCaching: false,
      selectable: true,
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
    lockMovementY: true,
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
    case 2: // pentagon arrows
      arrowPathStr = `M ${stemLength} -10
                      L 0 -10
                      L 0 20
                      L ${stemLength} 20
                      L ${stemLength} 20
                      L ${stemLength} 20
                      L ${length} 5
                      L ${stemLength} -10 z`;
      break;
    case 3: // pentagon (wide) arrows
      arrowPathStr = `M 0 0
                      L ${stemLength} 0
                      L ${length} 20
                      L ${stemLength} 40
                      L 0 40
                      L 0 0 z`;
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
    left: (1.5 + gene.start + nt_disps[genomeID]) * scaleFactor,
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
    let [l1, r1] = getGenomeRenderWindow(genomeID_A);
    let [l2, r2] = getGenomeRenderWindow(genomeID_B);

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

      g1 = g1.map(val => (val + nt_disps[genomeID_A]) * scaleFactor)
      g2 = g2.map(val => (val + nt_disps[genomeID_B]) * scaleFactor);

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
 *  Shift genomes horizontally to align around the target gene's centers. Pan viewport to first target gene (on first genome).
 *  Note: if there are multiple genes for a single genome, it will center around the first (lowest geneID) gene in the genome.
 *
 *  @param genes: list of target genes in format [{genomeID: 'ABC', geneID: 1}]
 */
GenomeDrawer.prototype.centerGenes = function (genes, centerToGeneStart=false) {
  let firstGenome = true;
  let basePad = 0;
  let centeredGenes = [];
  this.settings['genomeData']['genomes'].map(g => g[0]).forEach(genomeID => {
    let geneIDs = genes.filter(gene => gene.genomeID == genomeID).map(gene => gene.geneID).sort();
    if(geneIDs.length == 0) return;

    // select targetGeneID
    let targetGeneID = geneIDs[0];
    if(geneIDs.length > 1) {
      toastr.warning(`${geneIDs.length} gene hits were found on genome ${genomeID}, so the first one (gene ID ${targetGeneID}) was selected as the anchor gene.`);
    }
    let genePos = centerToGeneStart ? getGeneStart(genomeID, targetGeneID) : getGeneMid(genomeID, targetGeneID);
    centeredGenes.push({genomeID: genomeID, geneID: targetGeneID});
    
    if(firstGenome) {
      firstGenome = false;
      basePad = nt_disps[genomeID] + genePos;
      // pan to this gene
      let gene = this.settings['genomeData']['genomes'].find(g=>g[0]==genomeID)[1].genes.gene_calls[targetGeneID];
      let len = gene.stop - gene.start;
      let ntsToShow = parseInt($('#brush_end').val()) - parseInt($('#brush_start').val());
      let extraNts = ntsToShow - len;
      moveTo(gene.start - extraNts/2, gene.stop + extraNts/2, genomeID);
      return;
    }

    nt_disps[genomeID] += (basePad - (nt_disps[genomeID]+genePos));
  });

  // reenable scale if all genomes are aligned
  let disps = Object.values(nt_disps);
  if(disps.every(x => x==disps[0])) {
    for (genome of this.settings['genomeData']['genomes']) {
      nt_disps[genome[0]] = 0;
    }
    moveTo(parseInt($('#brush_start').val())-disps[0], parseInt($('#brush_end').val())-disps[0]) // shift viewport to same location
    drawScale();
    bindViewportToWindow();
    updateScalePos();
    updateRenderWindow();
    slidingActive = false;
    toggleScaleAttributes();
    toastr.success('Genomes aligned perfectly! Scale and bookmarks have been reenabled :)');
  } else {
    if(!slidingActive) {
      slidingActive = true;
      toggleScaleAttributes();
    }
  }

  this.draw();
  this.glowGenes(centeredGenes);
}

GenomeDrawer.prototype.centerGenesToProp = function (category, type, value, centerToGeneStart=false) {
  let targetGenes;
  switch(category) {
    case 'annotation':
      targetGenes = getGenesWithAnnotation(type, value);
      break;
    case 'metadata':
      targetGenes = getGenesWithMetadata(type, value);
      break;
    default:
      break;
  }
  if(targetGenes.length > 0) {
    this.centerGenes(targetGenes, centerToGeneStart);
  }
}

GenomeDrawer.prototype.centerGenesFromPropUI = function () {
  let category, type;
  switch($('#center_genomes_category').val()) {
    case 'metadata tag':
      category = 'metadata';
      type = 'tag';
      break;
    case 'user-defined annotation':
      category = 'metadata';
      type = 'annotation';
      break;
    default:
      category = 'annotation';
      type = $('#center_genomes_category').val();
      break;
  }

  let value = $('#center_genes_prop_value').val();

  this.centerGenesToProp(category, type, value, $('#center_to_gene_start_box').is(':checked'));
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

      let nt_shift = targetGeneMid - geneMid + nt_disps[firstGenomeID] - nt_disps[gid];
      let x_shift = nt_shift * scaleFactor;
      let objs = canvas.getObjects().filter(obj => obj.groupID == gid);
      for (o of objs) o.left += x_shift;
      nt_disps[gid] += nt_shift;
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
  if(!slidingActive) return;

  slidingActive = false;
  toggleScaleAttributes();

  for (genome of this.settings['genomeData']['genomes']) {
    nt_disps[genome[0]] = 0;
  }
  //percentScale = false;
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
  drawGenomeLabels(settings['display']['genome-label-size']);
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
  if (newSize < 0 || newSize > 100) {
    alert(`Invalid value, genome label size must be in range 0-100.`);
    return;
  }
  drawGenomeLabels(newSize);
  setLabelCanvas();
}

GenomeDrawer.prototype.redrawSingleGenome = function (genomeID) {
  if(!$('#' + genomeID + '-show').is(':checked')) return;

  canvas.getObjects().filter(o => o.groupID == genomeID).forEach(obj => canvas.remove(obj));
  let idx = this.settings['genomeData']['genomes'].findIndex(obj => obj[0] == genomeID);
  let orderIndex = idx;
  for(genome of this.settings['genomeData']['genomes']) {
    if(genome[0] == genomeID) break;
    if(!$('#' + genome[0] + '-show').is(':checked')) orderIndex--;
  }
  this.addLayers(idx, orderIndex);
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
  $('#query-results-head').empty().append(`
    <tr>
      <th>Gene ID</th>
      <th>Genome</th>
      <th>Start</th>
      <th>Stop</th>
      <th>Go To</th>
    </tr>
  `);
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
  if (category.split(' ')[0] == 'metadata'){
    type = category.split(' ')[1]
    drawer.queryMetadata(query, type)
    return;
  }

  let matchByAccession = false;
  this.settings['genomeData']['genomes'].map(genome => {
    // if we start with finding accession, stay with accession
    for (const [key, value] of Object.entries(genome[1]['genes']['functions'])) {
      if (value[category]?.[0].toLowerCase().includes(query)) { // accession
        if(!matchByAccession) {
          matchByAccession = true;
          glowPayload = [];
          foundInGenomes = {};
          distinctQueryMatches = {};
        }
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
      else if (!matchByAccession && value[category]?.[1].toLowerCase().includes(query)) { // annotation
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
        let [globalStart, globalEnd] = [start, end].map(pos => pos + nt_disps[gene['genomeID']])
        if (globalStart < lowestStart || lowestStart == null) lowestStart = globalStart
        if (globalEnd > highestEnd || highestEnd == null) highestEnd = globalEnd
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
      [lowestStart, highestEnd] = calcNTBounds();
    }
  }
  renderAllGenes()

  $('#function-query-results-statement').text(`Retreived ${glowPayload.length} hit(s) from ${Object.keys(foundInGenomes).length} genomes`)
  await zoomOutAndWait('partial', lowestStart, highestEnd, 350)
  this.glowGenes(glowPayload, true)
}

GenomeDrawer.prototype.queryMetadata = async function(metadataLabel, type){
  $('#query-results-table').empty()
  $('#query-results-head').empty().append(`
    <tr>
      <th>Gene ID</th>
      <th>Genome</th>
      <th>Start</th>
      <th>Stop</th>
      <th>Go To</th>
    </tr>
  `);

  if(!settings['display']['metadata']) {
    alert(`No hits were found matching ${metadataLabel} in ${type == 'annotation' ? 'user-defined annotations' : 'metadata'}`)
    return
  }
  let glowPayload = Array()
  let foundInGenomes = Object()
  let typeMatches = settings['display']['metadata'].filter(m => m.type == type);
  let matches;
  if(type == 'annotation') {
    matches = typeMatches.filter(m => m.accession.toLowerCase().includes(metadataLabel.toLowerCase()))
    if(matches.length == 0) {
      matches = typeMatches.filter(m => m.annotation.toLowerCase().includes(metadataLabel.toLowerCase()))
    }
  } else {
    matches = typeMatches.filter(m => m.label.toLowerCase().includes(metadataLabel.toLowerCase()))
  }
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
    alert(`No hits were found matching ${metadataLabel} in ${type == 'annotation' ? 'user-defined annotations' : 'metadata'}`)
    return
  }
  let lowestStart, highestEnd = null
  if (globalGenomeMax > 35000) { // TODO: instead of globalGenomeMax, use genomeMax[genomeID] for currently selected genomeID
    glowPayload.map(gene => {
      let genomeOfInterest = this.settings['genomeData']['genomes'].filter(genome => genome[0] == gene['genomeID'])
      let start = genomeOfInterest[0][1]['genes']['gene_calls'][gene['geneID']]['start']
      let end = genomeOfInterest[0][1]['genes']['gene_calls'][gene['geneID']]['stop']
      let [globalStart, globalEnd] = [start, end].map(pos => pos + nt_disps[gene['genomeID']])
      if (globalStart < lowestStart || lowestStart == null) lowestStart = globalStart
      if (globalEnd > highestEnd || highestEnd == null) highestEnd = globalEnd
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
    [lowestStart, highestEnd] = calcNTBounds()
  }
  await zoomOutAndWait('partial', lowestStart, highestEnd, 350)
  this.glowGenes(glowPayload, true)
}

GenomeDrawer.prototype.showAllMetadata = function(type){
  if(!settings['display']['metadata'] || settings['display']['metadata'].filter(metadata => metadata.type == type).length == 0) {
    alert(`No metadata ${type}s currently exist`);
    return;
  }

  $('#query-results-table').empty();
  $('#query-results-head').empty().append(`
    <tr>
      ${type == 'tag' ? '<th>Tag</th>' : ''}
      <th>Gene ID</th>
      <th>Genome</th>
      <th>Start</th>
      <th>Stop</th>
      <th>Go To</th>
      <th>Remove</th>
    </tr>
  `);

  settings['display']['metadata'].filter(metadata => metadata.type == type).forEach(metadata => {
    let genome = this.settings['genomeData']['genomes'].filter(genome => genome[0] == metadata.genome)
    let start = genome[0][1]['genes']['gene_calls'][metadata.gene]['start']
    let end = genome[0][1]['genes']['gene_calls'][metadata.gene]['stop']

    $('#query-results-table').append(`
      <tr>
        ${type == 'tag' ? '<td>' + metadata.label + '</td>' : ''}
        <td>${metadata.gene}</td>
        <td>${metadata.genome}</td>
        <td>${start}</td>
        <td>${end}</td>
        <td><button onclick="goToGene('${metadata.genome}',${metadata.gene},${start},${end})">Go To</button</td>
        <td><button onclick="$(this).closest('tr').remove(); removeMetadataFromQueryTable('${metadata.genome}',${metadata.gene},'${type=='tag' ? metadata.label : null}')">Remove ${type}</button</td>
      </tr>
    `)
  });

  removeMetadataFromQueryTable = (genomeID, geneID, label=null) => {
    let index = settings['display']['metadata'].findIndex(m => (type=='tag' ? m.label == label : true) && m.gene == geneID && m.genome == genomeID && m.type == type);
    settings['display']['metadata'].splice(index, 1);
    if($('#query-results-table').children().length == 0) {
      $('#query-results-head').empty();
    }
  }
}

GenomeDrawer.prototype.showAllUserDefined = function(){
  if(!settings['display']['metadata'] || settings['display']['metadata'].filter(metadata => metadata.type == 'annotation').length == 0) {
    alert('No user-defined annotations currently exist');
    return;
  }

  $('#query-results-table').empty();
  $('#query-results-head').empty().append(`
    <tr>
      <th>Gene ID</th>
      <th>Genome</th>
      <th>Start</th>
      <th>Stop</th>
      <th>Go To</th>
      <th>Remove</th>
    </tr>
  `);

  settings['display']['metadata'].filter(metadata => metadata.type == 'annotation').forEach(annotation => {
    let genome = this.settings['genomeData']['genomes'].filter(genome => genome[0] == annotation.genome)
    let start = genome[0][1]['genes']['gene_calls'][annotation.gene]['start']
    let end = genome[0][1]['genes']['gene_calls'][annotation.gene]['stop']

    $('#query-results-table').append(`
      <tr>
        <td>${annotation.gene}</td>
        <td>${annotation.genome}</td>
        <td>${start}</td>
        <td>${end}</td>
        <td><button onclick="goToGene('${annotation.genome}',${annotation.gene},${start},${end})">Go To</button</td>
        <td><button onclick="$(this).closest('tr').remove(); removeAnnotation('${annotation.genome}',${annotation.gene})">Remove Tag</button</td>
      </tr>
    `)
  });

  removeAnnotation = (genomeID, geneID) => {
    let index = settings['display']['metadata'].findIndex(m => m.gene == geneID && m.genome == genomeID && m.type == 'annotation');
    settings['display']['metadata'].splice(index, 1);
    if($('#gene_label_source').val() == 'user') this.redrawSingleGenome(genomeID);
  }
}

GenomeDrawer.prototype.setInitialZoom = function(){
    let start = 0
    let stop = globalGenomeMax > 35000 ? 35000 : globalGenomeMax
    zoomOut('partial', start, stop)
}
