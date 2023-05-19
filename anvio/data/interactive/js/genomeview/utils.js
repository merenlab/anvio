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
 * File Overview : This file contains utility functions used throughout genomeview. As a general rule,
 * functions defined here explicitly return some value or mutate an existing global variable. Functions
 * defined here should not interact directly with elements of UI, state, or canvas objects.
 */


/*
 *  return height value for main canvas element
 */
function calculateMainCanvasHeight(){
  let additionalSpacing = 100 // arbitrary additional spacing for cosmetics
  let scaleHeight = 230
  let nGenomes = $('.genome_selectors:checked').length
  let mainCanvasHeight = marginTop + nGenomes*(spacing+maxGroupSize*groupLayerPadding) + groupMargin*(nGenomes-1) + scaleHeight + additionalSpacing
  return mainCanvasHeight
}

/*
 *  Save NT length of each genome in `genomeMax` dict, and save length of the largest genome in `globalGenomeMax`.
 */
function calculateGenomeLengths(){
  globalGenomeMax = -999999;
  for(genome of genomeData.genomes) {
    let genomeID = genome[0]
    genome = genome[1].genes.gene_calls;
    let genomeEnd = Object.values(genome)[Object.keys(genome).length-1].stop;
    genomeMax[genomeID] = genomeEnd;
    if(genomeEnd > globalGenomeMax) globalGenomeMax = genomeEnd;
  }
}

/*
 *  @returns NT position of the middle of each gene in a given genome with a specified gene cluster
 */
function getGenePosForGenome(genomeID, gc) {
  var targetGenes = getGenesOfGC(genomeID, gc);
  if(targetGenes == null) return null;

  let genome = settings['genomeData']['genomes'].find(g => g[0] == genomeID);
  let mids = [];
  for(geneID of targetGenes) {
    let gene = genome[1].genes.gene_calls[geneID];
    let geneMid = gene.start + (gene.stop - gene.start) / 2;
    mids.push(geneMid);
  }
  return mids;
}

/*
 *  @returns array of geneIDs in a given genome with a specified gene cluster
 */
function getGenesOfGC(genomeID, gc) {
  var targetGenes = settings['genomeData']['gene_associations']["anvio-pangenome"]["gene-cluster-name-to-genomes-and-genes"][gc][genomeID];
  return targetGenes.length > 0 ? targetGenes : null;
}

/* 
 *  @returns NT position of the middle of a given gene
 */
function getGeneMid(genomeID, geneID) {
  let gene = settings['genomeData']['genomes'].find(g => g[0]==genomeID)[1].genes.gene_calls[geneID];
  if(!gene) return null;
  return (gene.start + gene.stop) / 2;
}

function getGeneStart(genomeID, geneID) {
  let gene = settings['genomeData']['genomes'].find(g => g[0]==genomeID)[1].genes.gene_calls[geneID];
  if(!gene) return null;
  return gene.start;
}

/*
 *  @returns array of genes in the form [{genomeID: 'ABC', geneID: 1}] with a specified functional annotation
 */
function getGenesWithAnnotation(annotation_type, annotation_value) {
  let targetGenes = [];
  settings['genomeData']['genomes'].forEach(genome => {
    let geneFunctions = genome[1]['genes']['functions'];
    if(geneFunctions) {
      for(geneID in geneFunctions) {
        let annotation_info = geneFunctions[geneID][annotation_type];
        if(annotation_info && annotation_info[0].toLowerCase().includes(annotation_value.toLowerCase()) || annotation_info[1].toLowerCase().includes(annotation_value.toLowerCase())) {
          targetGenes.push({genomeID: genome[0], geneID: geneID});
        }
      }
    }
  });
  return targetGenes;
}

/*
 *  @returns array of genes in the form [{genomeID: 'ABC', geneID: 1}] with a specified metadata value
 */
function getGenesWithMetadata(metadata_type, metadata_value) {
  let targetGenes = [];

  if(!settings.display.metadata) return [];
  settings['display']['metadata'].filter(m => m.type == metadata_type).forEach(metadata => {
    if(metadata_type == 'tag') {
      if(metadata.label.toLowerCase() == metadata_value.toLowerCase()) {
        targetGenes.push({genomeID: metadata.genome, geneID: metadata.gene});
      }
    } else if(metadata_type == 'annotation') {
      if(metadata.accession == metadata_value || metadata.annotation == metadata_value) {
        targetGenes.push({genomeID: metadata.genome, geneID: metadata.gene});
      }
    }
  });
  
  return targetGenes;
}

/*
 *  Show/hide gene labels to show the max amount possible s.t. none overlap.
 */
function checkGeneLabels() {
  var labels = canvas.getObjects().filter(obj => obj.id == 'geneLabel');
  for(var i = 0; i < labels.length-1; i++) {
    if(this.settings['display']['arrow-style'] == 3) {
      // hide labels that don't fit inside pentagon arrows
      if(labels[i].width/2 > canvas.getObjects().filter(obj => obj.id == 'arrow')[i].width) {
        labels[i].visible = false;
        continue;
      }
      labels[i].visible = true;
    }
    var p = i+1;
    while(p < labels.length && labels[i].intersectsWithObject(labels[p])) {
      labels[p].visible = false;
      p++;
    }
    if(p == labels.length) return;
    labels[p].visible = true;
    i = p - 1;
  }
}

 /*
 *  @returns [start, stop] nt range for the current viewport and scaleFactor
 */
 function getNTRangeForVPT() {
  let vpt = canvas.viewportTransform;
  let window_left = Math.floor((-1*vpt[4])/scaleFactor);
  let window_right = Math.floor(window_left + canvas.getWidth()/scaleFactor);
  let [lbound, rbound] = calcNTBounds();
  // if window is out of bounds, shift to be in bounds
  if(window_left < lbound) {
    window_right += (lbound - window_left);
    window_left = lbound;
  }
  if(window_right > rbound) {
    window_left -= (window_right - rbound);
    window_right = rbound;
  }
  return [window_left, window_right];
}

/* 
 *  @returns x-range of a specific genome based on its current nt range (since they can be dragged)
 *  - used to draw genome line, ruler
 */
function getRenderXRangeForGenome(genomeID) {
  return renderWindow.map(pos => Math.floor((pos+nt_disps[genomeID])*scaleFactor));
}

/* 
 *  @returns nt range of a specific genome that appears in the current render window (since they can be dragged)
 *  - used to determine which gene arrows to draw for a given genome while genome sliding is activated
 */
function getGenomeRenderWindow(genomeID) {
  if(!slidingActive) return renderWindow;

  let [start, end] = renderWindow.map(pos => Math.floor(pos - nt_disps[genomeID]));
  return [clamp(start,0,genomeMax[genomeID]), clamp(end,0,genomeMax[genomeID])];
}

/*
 *  @returns array [min, max] where
 *    min = x-start of the given genome, max = x-end of the given genome
 */
function calcXBoundsForGenome(genomeID) {
  return [nt_disps[genomeID]*scaleFactor, (nt_disps[genomeID]+genomeMax[genomeID])*scaleFactor];
}

/*
 *  @returns array [min, max] where
 *    min = x-start of the leftmost genome, max = x-end of the rightmost genome
 */
function calcXBounds() {
  let min = 9*(10**9), max = -9*(10**9);
  for(genome of genomeData.genomes) {
    let [start, end] = calcXBoundsForGenome(genome[0]);
    if(start < min) min = start;
    if(end > max) max = end;
  }
  return [min, max];
}

/*
 *  @returns array [lowestStart, highestEnd] of nt positions for all genomes
 *    note if slidingActive = false, this will always equal [0, globalGenomeMax]
 */
function calcNTBounds() {
  return calcXBounds().map(x => x/scaleFactor);
}

/*
 *  @returns array of functional annotation types from table in `genomeData`
 */
function getFunctionalAnnotations() {
  let fns = genomeData.genomes[0][1].genes.functions;
  if(fns && Object.values(fns).length > 0) return Object.keys(Object.values(fns)[0]);
  return []; // return empty array if no functional annotations found
}

/*
 *  @returns arbitrary category:color dict given a list of categories
 */
function getCustomColorDict(fn_type, cags=null, order=null) {
  if(!getFunctionalAnnotations().includes(fn_type)) return null;

  if(!cags) {
    cags = [];
    genomeData.genomes.forEach(genome => {
      Object.values(genome[1].genes.functions).forEach(fn => {
        let cag = getCagForType(fn, fn_type);
        if(cag && !cags.includes(cag)) cags.push(cag);
        if(!cag && !cags.includes("None")) cags.push("None");
      });
    });
  }

  // move "Other" and "None" to end of list
  if(cags.includes("Other")) cags.push(cags.splice(cags.indexOf("Other"), 1)[0]);
  if(cags.includes("None")) cags.push(cags.splice(cags.indexOf("None"), 1)[0]);

  let dict = custom_cag_colors.reduce((out, field, index) => {
    out[cags[index]] = field;
    return out;
  }, {});

  // sort using order
  if(order) {
    let colors = Object.values(dict);
    Object.keys(dict).forEach(cag => { dict[cag] = colors[order[cag]] });
  }

  if(dict["Other"]) dict["Other"] = "#FFFFFF";
  if(dict["None"]) dict["None"] = "#808080";
  delete dict["undefined"];
  return dict;
}

function orderColorTable(order) {
  order_gene_colors_by_count = order == 'count';
  generateColorTable(null, $("#gene_color_order").val());
}

function filterColorTable(thresh) {
  if(isNaN(thresh)) {
    alert("Error: filtering threshold must be numeric");
    return;
  } else if(thresh < 1) {
    alert("Error: filtering threshold must be an integer >= 1");
    return;
  }
  settings['display']['thresh-count-gene-colors'] = thresh;
  generateColorTable(null, $("#gene_color_order").val());
  drawer.draw();
}

function toggleRightPanel(name) {
  ['#mouseover-panel', '#settings-panel', '#query-panel'].forEach(function(right_panel) {
      if (right_panel == name)
          return;

      $(right_panel).hide();
  });
  $(name).toggle();

  if ($('#mouseover-panel').is(':visible')) {
      $('#toggle-panel-mouseover').addClass('toggle-panel-mouseover-pos');
      $('#toggle-panel-mouseover-inner').html('&#9658;');
  } else {
      $('#toggle-panel-mouseover').removeClass('toggle-panel-mouseover-pos');
      $('#toggle-panel-mouseover-inner').html('&#9664;');
  }

  if ($('#settings-panel').is(':visible')) {
      $('#toggle-panel-settings').addClass('toggle-panel-settings-pos');
      $('#toggle-panel-settings-inner').html('&#9658;');
  } else {
      $('#toggle-panel-settings').removeClass('toggle-panel-settings-pos');
      $('#toggle-panel-settings-inner').html('&#9664;');
  }

  if ($('#query-panel').is(':visible')) {
      $('#toggle-panel-query').addClass('toggle-panel-query-pos');
      $('#toggle-panel-query-inner').html('&#9658;');
  } else {
      $('#toggle-panel-query').removeClass('toggle-panel-query-pos');
      $('#toggle-panel-query-inner').html('&#9664;');
  }
}





// ================================================================================================================================
// DEPRECATED FUNCTIONS
// ================================================================================================================================

/*
 *  [NOTE: DEPRECATED]
 *  @returns [start, stop] proportional (0-1) range, used with scale for non-aligned genomes
 */
function getFracForVPT() {
  let resolution = 4; // number of decimals to show
  let [x1, x2] = calcXBounds();
  let window_left = Math.round(10**resolution * (-1*canvas.viewportTransform[4] - x1) / (x2 - x1)) / 10**resolution;
  let window_right = Math.round(10**resolution * (window_left + (canvas.getWidth()) / (x2 - x1))) / 10**resolution;
  // if window is out of bounds, shift to be in bounds
  if(window_left < 0) {
    window_right -= window_left;
    window_left = 0;
  }
  if(window_right > 1) {
    window_left -= (window_right - 1);
    window_right = 1;
  }
  return [window_left, window_right];
}

/*
 *  [NOTE: DEPRECATED]
 *  @returns absolute x-position viewport range for a given proportional selection range
 */
function getVPTForFrac() {
  if(!percentScale) return null;
  let [l,r] = calcXBounds();
  return [$('#brush_start').val(),$('#brush_end').val()].map(x => l+x*(r-l));
}

/*
 *  [NOTE: DEPRECATED]
 *  @returns absolute x-position viewport range for a given proportional renderWindow range
 */
function getRenderXRangeForFrac() {
  if(!percentScale) return null;
  let [l,r] = calcXBounds();
  let [x1, x2] = renderWindow.map(x => l+x*(r-l));
  return [x1, x2];
}