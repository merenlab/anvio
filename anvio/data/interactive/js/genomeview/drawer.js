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
var gDrawer = function(settings) {
  this.settings = settings;
};

gDrawer.prototype.draw = function(){
  canvas.clear()
  labelSpacing = 30 // reset to default value upon each draw() call
  canvas.setHeight(calculateMainCanvasHeight()) // set canvas height dynamically

  genomeData['genomes'].map((genome, idx) => {
    addGenome(idx)
    addLayers(idx)
    labelSpacing += 30
  })

  checkGeneLabels();
  drawTestShades();
}

gDrawer.prototype.addGenome = function(orderIndex){
  let genome = genomeData.genomes[orderIndex];
  let gene_list = genome[1].genes.gene_calls;
  let genomeLabel = gene_list[0].contig;
  let genomeID = genome[0];
  let y = marginTop + orderIndex*spacing;

  let layerHeight = spacing / maxGroupSize
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
  addBackgroundShade(y, start, genomeMax, layerHeight, orderIndex)

  for(let geneID in gene_list) {
    let gene = gene_list[geneID];
    let [ntStart, ntStop] = getRenderNTRange(genomeID);
    if(gene.start < ntStart) continue;
    if(gene.stop > ntStop) return;
    var geneObj = geneArrow(gene,geneID,y,genomeID,arrowStyle);
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
 *  Draw background shades between genes of the same cluster.
 *  TODO: generalize this function to take in [start,stop,val] NT sequence ranges to shade any arbitrary metric
 *        - add a separate function for retrieving [start,stop,val] given gene cluster IDs
 *
 *  @param geneClusters : array of GC IDs to be shaded
 *  @param colors       : dict defining color of each shade, in the form {geneClusterID : hexColor}
 */
gDrawer.prototype.shadeGeneClusters = function(geneClusters, colors){
  if(!genomeData.gene_associations["anvio-pangenome"]) return;

  let y = marginTop;
  for(var i = 0; i < genomeData.genomes.length-1; i++) {
    let genomeA = genomeData.genomes[i][1].genes.gene_calls;
    let genomeB = genomeData.genomes[i+1][1].genes.gene_calls;
    let genomeID_A = genomeData.genomes[i][0];
    let genomeID_B = genomeData.genomes[i+1][0];
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

gDrawer.prototype.glowGenes = function(geneParams){
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
 *  Clear all gene links from the canvas.
 */
gDrawer.prototype.clearShades = function(){
  canvas.getObjects().filter(obj => obj.id == 'link').forEach((l) => { canvas.remove(l) });
}