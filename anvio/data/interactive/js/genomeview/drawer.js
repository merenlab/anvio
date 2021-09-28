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
var GenomeDrawer = function(settings) {
  this.settings = settings;
};

GenomeDrawer.prototype.draw = function(){
  canvas.clear()
  labelSpacing = 30 // reset to default value upon each draw() call
  canvas.setHeight(calculateMainCanvasHeight()) // set canvas height dynamically

  this.settings['genomeData']['genomes'].map((genome, idx) => {
    this.addLayers(idx)
    labelSpacing += 30
  })

  checkGeneLabels();
  drawTestShades();
}

/*
 *  For each genome group, iterate additional all layers and render where appropriate
 */
GenomeDrawer.prototype.addLayers = function(orderIndex){
  let [dataLayerHeight, rulerHeight] = [this.calculateLayerSizes()[0], this.calculateLayerSizes()[1]]

  yOffset = orderIndex * spacing;
  let layerPos = 0
  let genomeID = this.settings['genomeData']['genomes'][orderIndex][0];
  let genome = this.settings['genomeData']['genomes'][orderIndex][1];
  let label = genome.genes.gene_calls[0].contig;

  let additionalDataLayers = this.settings['additional-data-layers']['data'][genomeID]

  let ptInterval = Math.floor(genomeMax / adlPtsPerLayer);

  this.settings['group-layer-order'].map((layer, idx) => {  // render out layers, ordered via group-layer-order array
    if(layer == 'Genome' && $('#Genome-show').is(':checked')){
      this.addGenome(orderIndex, dataLayerHeight, layerPos)
      layerPos += dataLayerHeight
    }
    if(layer == 'Coverage' && this.settings['additional-data-layers']['layers'].includes('Coverage') && $('#Coverage-show').is(':checked')){
      this.buildNumericalDataLayer('Coverage', layerPos, genomeID, additionalDataLayers, ptInterval, 'blue', dataLayerHeight, orderIndex)
      layerPos += dataLayerHeight
    }
    if(layer == 'GC_Content' && this.settings['additional-data-layers']['layers'].includes('GC_content') && $('#GC_Content-show').is(':checked')){
      this.buildNumericalDataLayer('GC_content', layerPos, genomeID, additionalDataLayers, ptInterval, 'purple', dataLayerHeight, orderIndex)
      layerPos += dataLayerHeight
    }
    if(layer == 'Ruler' && this.settings['additional-data-layers']['layers'].includes('ruler') && $('#Ruler-show').is(':checked')) {
      this.buildGroupRulerLayer(genomeID, layerPos, rulerHeight, orderIndex)
      layerPos += rulerHeight
    }
  })
}

/*
 *  programmatically calculate layer height values, given that the ruler layer should be allocated comparatively less space
 */
GenomeDrawer.prototype.calculateLayerSizes = function(){
  let parityHeight = spacing / maxGroupSize
  let rulerHeight = Math.floor(parityHeight * .5) // some arbitrary percentage of parity since ruler should get less y-axis space

  // with the extra space carved out by a smaller rulerHeight, distribute the excess evenly amongst all layers that are NOT rulers
  let dataLayerHeight = Math.floor(parityHeight * (1 + (.5 / (maxGroupSize -1))) )
  return [dataLayerHeight, rulerHeight]
}

GenomeDrawer.prototype.addGenome = function(orderIndex, layerHeight, layerPos){
  let genome = this.settings['genomeData']['genomes'][orderIndex];
  let gene_list = genome[1].genes.gene_calls;
  let genomeLabel = gene_list[0].contig;
  let genomeID = genome[0];
  let y = marginTop + yOffset + layerPos + (layerHeight / 2) // render arrows in the center of genome layer's allotted vertical space

  if(showLabels) {
    canvas.add(new fabric.Text(genomeLabel, {top: y-5, selectable: false, fontSize: genomeLabelSize, fontFamily: 'sans-serif', fontWeight: 'bold'}));
  }

  let [start, stop] = percentScale ? getRenderXRangeForFrac() : renderWindow.map(x => x*scaleFactor + xDisps[genomeID]);
  start = clamp(start > xDisps[genomeID] ? start : xDisps[genomeID], calcXBounds()[0], calcXBounds()[1]);
  stop = clamp(stop, calcXBounds()[0], calcXBounds()[1]);

  // line
  let lineObj = new fabric.Line([start,0,stop,0], {
        id: 'genomeLine',
        groupID: genomeID,
        top: y + 4,
        stroke: 'black',
        strokeWidth: 2,
        lockMovementY: true,
        hasControls: false,
        hasBorders: false,
        lockScaling: true});
  canvas.add(lineObj);
  this.addBackgroundShade((marginTop + yOffset + layerPos), start, genomeMax, layerHeight, orderIndex)

  for(let geneID in gene_list) {
    let gene = gene_list[geneID];
    let [ntStart, ntStop] = getRenderNTRange(genomeID);
    if(gene.start < ntStart) continue;
    if(gene.stop > ntStop) return;
    var geneObj = this.geneArrow(gene,geneID,y,genomeID,arrowStyle);
    canvas.add(geneObj);

    if(showGeneLabels) {
      var label = new fabric.IText("geneID: "+geneID, {
        id: 'geneLabel',
        groupID: genomeID,
        fontSize: geneLabelSize,
        angle: geneLabelPos == "slanted" ? -10 : 0,
        left: xDisps[genomeID]+(gene.start+50)*scaleFactor,
        scaleX: 0.5,
        scaleY: 0.5,
        hasControls: false,
        lockMovementX: true,
        lockMovementY: true,
        lockScaling: true,
        hoverCursor: 'text'
      });
      if(arrowStyle == 3) {
        label.set({
          top: geneLabelPos == "inside" ? y-5 : y-30,
          selectionColor:'rgba(128,128,128,.5)'
        });
      } else {
        label.set({
          top: y-30,
          selectionColor:'rgba(128,128,128,.2)'
        });
      }
      canvas.add(label);
    }
  }
}
/*
 *  Process to generate numerical ADL for genome groups (ie Coverage, GC Content )
 */
GenomeDrawer.prototype.buildNumericalDataLayer = function(layer, layerPos, genomeID, additionalDataLayers, ptInterval, defaultColor, layerHeight, orderIndex){
    // TODO this will need to be refactored once we begin testing genomes comprised of multiple contigs
    let contigObj = Object.values(additionalDataLayers)[0]
    let contigArr = Object.values(contigObj)[0]

    let maxDataLayerValue = 0
    let startingTop = marginTop + yOffset + layerPos
    let startingLeft = xDisps[genomeID]
    let pathDirective = [`M 0 0`]

    for(let i = 0; i < contigArr.length; i++){
      contigArr[i] > maxDataLayerValue ? maxDataLayerValue = contigArr[i] : null
    }

    let nGroups = 20
    let j = 0
    let [l,r] = getRenderNTRange(genomeID);
    for(let i = 0; i < nGroups; i++) {
      for(; j < i*genomeMax/nGroups; j+=ptInterval){
        if(j < l) continue;
        if(j > r) break;
        let left = j * scaleFactor + startingLeft
        let top = [contigArr[j] / maxDataLayerValue] * layerHeight
        let segment = `L ${left} ${top}`
        pathDirective.push(segment)
      }
      let graphObj = new fabric.Path(pathDirective.join(' '))
      graphObj.set({
        top : startingTop,
        stroke : 'black', //additionalDataLayers[layer] ? additionalDataLayers[`${layer}-color`] : defaultColor,
        fill : '', //additionalDataLayers['gcContent-color'] ? additionalDataLayers['gcContent-color'] : 'black',
        selectable: false,
        objectCaching: false,
        id : `${layer} graph`,
        groupID : genomeID,
        genome : genomeID
      })
      canvas.bringToFront(graphObj)
      pathDirective = []
    }
    this.addBackgroundShade(startingTop, startingLeft, genomeMax, layerHeight, orderIndex)
}

/*
 *  Generate individual genome group rulers
 */
GenomeDrawer.prototype.buildGroupRulerLayer = function(genomeID, layerPos, layerHeight, orderIndex){
  let startingTop = marginTop + yOffset + layerPos
  let startingLeft = xDisps[genomeID]
  // let layerHeight = (spacing / maxGroupSize)

  // split ruler into several objects to avoid performance cost of large object pixel size
  let nRulers = 20;
  let w = 0;
  let [l,r] = getRenderNTRange(genomeID);
  for(let i = 0; i < nRulers; i++) {
    let ruler = new fabric.Group();
    for(; w < (i+1)*genomeMax/nRulers; w+=scaleInterval) {
      if(w < l) continue;
      if(w > r) break;
      let tick = new fabric.Line([0,0,0,20], {left: (w*scaleFactor),
            stroke: 'black',
            strokeWidth: 1,
            fontSize: 10,
            fontFamily: 'sans-serif'});
      let lbl = new fabric.Text(w/1000 + " kB", {left: (w*scaleFactor+5),
            stroke: 'black',
            strokeWidth: .25,
            fontSize: 15,
            fontFamily: 'sans-serif'});
      ruler.add(tick);
      ruler.add(lbl);
    }
      ruler.set({
        left: startingLeft,
        top: startingTop,
        lockMovementY: true,
        hasControls: false,
        hasBorders: false,
        lockScaling: true,
        objectCaching: false,
        groupID: genomeID
      });
      ruler.addWithUpdate();
      canvas.add(ruler);
  }
  this.addBackgroundShade(startingTop, startingLeft, genomeMax, layerHeight, orderIndex)
}

/*
 *  adds an alternating shade to each genome group for easier visual distinction amongst adjacent groups
 */
GenomeDrawer.prototype.addBackgroundShade = function(top, left, width, height, orderIndex){
  let backgroundShade;
  orderIndex % 2 == 0 ? backgroundShade = '#b8b8b8' : backgroundShade = '#f5f5f5'

  let background = new fabric.Rect({
    groupID: this.settings['genomeData']['genomes'][orderIndex][0],
    top: top,
    left: left,
    width: width,
    height: height,
    fill: backgroundShade,
    selectable: false,
    opacity : .5
  });
  canvas.sendToBack(background)
}

GenomeDrawer.prototype.geneArrow = function(gene, geneID, y, genomeID, style){
  let ind = this.settings['genomeData']['genomes'].findIndex(g => g[0] == genomeID);
  let functions = this.settings['genomeData']['genomes'][ind][1].genes.functions[geneID];

  let color = $('#picker_None').length > 0 ? $('#picker_None').attr('color') : 'gray';
  let cag = getCagForType(functions, color_db);

  // TODO: use state instead of hardcoded color pickers

  // check if gene is highlighted
  let pickerCode = genomeID + '-' + geneID;
  if($('#picker_' + pickerCode).length > 0) {
    color = $('#picker_' + pickerCode).attr('color');
  } else {
    if(cag) {
       cag = getCleanCagCode(cag);
       let color_other = $('#picker_Other').length > 0 ? $('#picker_Other').attr('color') : 'white';
       color = $('#picker_' + cag).length > 0 ? $('#picker_' + cag).attr('color') : color_other;
    } else {
      if (gene.source.startsWith('Ribosomal_RNA')) {
        cag = 'rRNA';
      } else if (gene.source == 'Transfer_RNAs') {
        cag = 'tRNA';
      } else if (gene.functions !== null) {
        cag = 'Function';
      }
      if($('#picker_' + cag).length > 0) color = $('#picker_' + cag).attr('color');
    }
  }

  let length = (gene.stop-gene.start)*scaleFactor;
  let stemLength = length-25 > 0 ? length-25 : 0;

  var arrowPathStr;
  switch(style) {
    case 2: // thicker arrows
      arrowPathStr = 'M ' + stemLength + ' -5 L 0 -5 L 0 15 L ' + stemLength + ' 15 L ' + stemLength + ' 15 L ' + stemLength + ' 20 L ' + length + ' 5 L ' + stemLength + ' -10 z';
      break;
    case 3: // pentagon arrows
      arrowPathStr = 'M 0 0 L ' + stemLength + ' 0 L ' + length + ' 20 L ' + stemLength + ' 40 L 0 40 L 0 0 z';
      break;
    case 4: // rect arrows
      arrowPathStr = 'M ' + length + ' -5 L 0 -5 L 0 15 L ' + length + ' 15 z';
      break;
    default: // 'inspect page' arrows
      arrowPathStr = 'M ' + stemLength + ' 0 L 0 0 L 0 10 L ' + stemLength + ' 10 L ' + stemLength + ' 10 L ' + stemLength + ' 20 L ' + length + ' 5 L ' + stemLength + ' -10 z';
      break;
  }

  var arrow = new fabric.Path(arrowPathStr);
  arrow.set({
    id: 'arrow',
    groupID: genomeID,
    lockMovementY: true,
    hasControls: false,
    hasBorders: false,
    lockScaling: true,
    gene: gene,
    geneID: geneID,
    genomeID: genomeID,
    top: style == 3 ? y-17 : y-11, // TODO update this offset to reflect genome layer height (we want to render this arrow in the middle of its allocated height)
    left: xDisps[genomeID] + (1.5+gene.start)*scaleFactor,
    fill: color,
    stroke: 'gray',
    strokeWidth: style == 3 ? 3 : 1.5
  });
  if(gene.direction == 'r') arrow.rotate(180);

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
GenomeDrawer.prototype.shadeGeneClusters = function(geneClusters, colors){
  if(!genomeData.gene_associations["anvio-pangenome"]) return;

  let y = marginTop;
  for(var i = 0; i < genomeData.genomes.length-1; i++) {
    let genomeA = this.settings['genomeData']['genomes'][i][1].genes.gene_calls;
    let genomeB = this.settings['genomeData']['genomes'][i+1][1].genes.gene_calls;
    let genomeID_A = this.settings['genomeData']['genomes'][i][0];
    let genomeID_B = this.settings['genomeData']['genomes'][i+1][0];
    let [l1,r1] = getRenderNTRange(genomeID_A);
    let [l2,r2] = getRenderNTRange(genomeID_B);

    for(gc of geneClusters) {
      let g1 = [], g2 = [];
      for(geneID of genomeData.gene_associations["anvio-pangenome"]["gene-cluster-name-to-genomes-and-genes"][gc][genomeID_A]) {
        g1.push(genomeA[geneID].start, genomeA[geneID].stop);
      }
      for(geneID of genomeData.gene_associations["anvio-pangenome"]["gene-cluster-name-to-genomes-and-genes"][gc][genomeID_B]) {
        g2.push(genomeB[geneID].start, genomeB[geneID].stop);
      }

      // if shades outside render bounds, don't draw them
      if(g1[1] < l1 && g2[1] < l2) continue;
      if(g1[0] > r1 && g2[0] > r2) break;

      g1 = g1.map(val => val*scaleFactor + xDisps[genomeID_A]);
      g2 = g2.map(val => val*scaleFactor + xDisps[genomeID_B]);

      /* TODO: implementation for multiple genes of the same genome in the same gene cluster */
      var path = new fabric.Path("M " + g1[0] + " " + y + " L " + g1[1] + " " + y + " L " + g2[1] + " " + (y+spacing) + " L " + g2[0] + " " + (y+spacing) + " z", {
        id: 'link',
        fill: colors[gc],
        opacity: 0.25,
        selectable: false
      });
      path.sendBackwards();
      canvas.sendToBack(path);
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
GenomeDrawer.prototype.glowGenes = function(geneParams){
   // convert geneParams format (1) to format (2)
   if(Array.isArray(geneParams[0].geneID)) {
    let newParams = [];
    for(genome of geneParams) {
      for(gene of genome.geneID) newParams.push({genomeID:genome.genomeID, geneID:gene});
    }
    geneParams = newParams;
  }

  var shadow = new fabric.Shadow({
    color: 'red',
    blur: 30
  });
  var arrows = canvas.getObjects().filter(obj => obj.id == 'arrow' && geneParams.some(g => g.genomeID == obj.genomeID && g.geneID == obj.geneID));
  for(arrow of arrows) {
    arrow.set('shadow', shadow);
    arrow.animate('shadow.blur', 0, {
      duration: 3000,
      onChange: canvas.renderAll.bind(canvas),
      onComplete: function(){ arrow.shadow = null; },
      easing: fabric.util.ease['easeInQuad']
    });
  }
  canvas.renderAll();
}

/*
 *  Shift genomes horizontally to align genes around the target gene cluster.
 *
 *  @param gc : target gene cluster ID
 */
GenomeDrawer.prototype.alignToCluster = function(gc){
  if(!this.settings['genomeData']['gene_associations']["anvio-pangenome"]) return;

  let targetGeneInfo = viewCluster(gc);
  if(targetGeneInfo == null) return;
  let [firstGenomeID, targetGeneMid] = targetGeneInfo;
  if(firstGenomeID != null) {
    alignToGC = gc;
    let index = this.settings['genomeData']['genomes'].findIndex(g => g[0] == firstGenomeID);
    for(var i = index+1; i < this.settings['genomeData']['genomes'].length; i++) {
      let gid = genomeData.genomes[i][0];
      let geneMids = getGenePosForGenome(genomeData.genomes[i][0], alignToGC);
      if(geneMids == null) continue;
      let geneMid = geneMids[0]; /* TODO: implementation for multiple matching gene IDs */
      let shift = scaleFactor * (targetGeneMid - geneMid) + (xDisps[firstGenomeID] - xDisps[gid]);
      let objs = canvas.getObjects().filter(obj => obj.groupID == gid);
      for(o of objs) o.left += shift;
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
GenomeDrawer.prototype.clearShades = function(){
  canvas.getObjects().filter(obj => obj.id == 'link').forEach((l) => { canvas.remove(l) });
}

GenomeDrawer.prototype.setPtsPerADL = function(newResolution){
  if(isNaN(newResolution)) return;
  newResolution = parseInt(newResolution);
  if(newResolution < 0 || newResolution > genomeMax) {
    alert(`Invalid value, genome spacing must be in range 0-${genomeMax}.`);
    return;
  }
  adlPtsPerLayer = newResolution;
  this.draw();
}

GenomeDrawer.prototype.showAllADLPts = function(){
  this.setPtsPerADL(genomeMax);
  $('#showAllADLPtsBtn').blur();
}

GenomeDrawer.prototype.alignRulers = function(){
  for(genome of this.settings['genomeData']['genomes']) {
    xDisps[genome[0]] = xDisplacement;
  }
  percentScale = false;
  drawScale();
  bindViewportToWindow();
  updateScalePos();
  updateRenderWindow();
  this.draw();
  $('#alignRulerBtn').blur();
}

GenomeDrawer.prototype.setGenomeSpacing = function(newSpacing){
  if(isNaN(newSpacing)) return;
  newSpacing = parseInt(newSpacing);
  if(newSpacing < 0 || newSpacing > 1000) {
    alert(`Invalid value, genome spacing must be in range 0-1000.`);
    return;
  }
  spacing = newSpacing;
  this.draw();
}

GenomeDrawer.prototype.setScaleInterval = function(newScale){
  if(isNaN(newScale)) return;
  newScale = parseInt(newScale);
  if(newScale < 50) {
    alert(`Invalid value, scale interval must be >=50.`);
    return;
  }
  scaleInterval = newScale;
  this.draw();
}

GenomeDrawer.prototype.setGeneLabelSize = function(newSize){
  if(isNaN(newSize)) return;
  newSize = parseInt(newSize);
  if(newSize < 0 || newSize > 1000) {
    alert(`Invalid value, gene label size must be in range 0-1000.`);
    return;
  }
  geneLabelSize = newSize;
  if(showGeneLabels) this.draw();
}

GenomeDrawer.prototype.setGenomeLabelSize = function(newSize){
  if(isNaN(newSize)) return;
  newSize = parseInt(newSize);
  if(newSize < 0 || newSize > 1000) {
    alert(`Invalid value, genome label size must be in range 0-1000.`);
    return;
  }
  genomeLabelSize = newSize;
  if(showLabels) this.draw();
}

GenomeDrawer.prototype.redrawSingleGenome = function(genomeID){
  canvas.getObjects().filter(o => o.groupID == genomeID).forEach(obj => canvas.remove(obj));
  let idx = this.settings['genomeData']['genomes'].findIndex(obj => obj[0] == genomeID);
  this.addLayers(idx);
  checkGeneLabels();
}

/*
 *  Dynamically set scale tick interval based on scaleFactor.
 */
GenomeDrawer.prototype.adjustScaleInterval = function(){
  let val = Math.floor(100/scaleFactor);
  let roundToDigits = Math.floor(Math.log10(val)) - 1;
  let newInterval = Math.floor(val/(10**roundToDigits)) * (10**roundToDigits);
  scaleInterval = newInterval;
  $('#genome_scale_interval').val(scaleInterval);
}
