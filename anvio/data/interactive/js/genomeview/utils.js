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


/** 
 * ripped from primary utils.js instead of importing the whole file
*/
var fixHelperModified = function(e, tr) { 
  var $originals = tr.children();
  var $helper = tr.clone();
  $helper.children().each(function(index) {
      $(this).width($originals.eq(index).width());
  });
  return $helper;
};

/*
 *  return height value for main canvas element
 */
function calculateMainCanvasHeight(){
  let additionalSpacing = 100 // arbitrary additional spacing for cosmetics
      let mainCanvasHeight =  spacing * settings['genomeData']['genomes'].length + additionalSpacing
  return mainCanvasHeight
}

/*
 *  Save NT length of the largest genome in `genomeMax`.
 */
function calculateMaxGenomeLength(){
  for(genome of genomeData.genomes) {
    genome = genome[1].genes.gene_calls;
    let genomeEnd = genome[Object.keys(genome).length-1].stop;
    if(genomeEnd > genomeMax) genomeMax = genomeEnd;
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
 *  Show/hide gene labels to show the max amount possible s.t. none overlap.
 */
function checkGeneLabels() {
  var labels = canvas.getObjects().filter(obj => obj.id == 'geneLabel');
  for(var i = 0; i < labels.length-1; i++) {
    if(arrowStyle == 3) {
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
  let window_left = Math.floor((-1*vpt[4]-xDisplacement)/scaleFactor);
  let window_right = Math.floor(window_left + canvas.getWidth()/scaleFactor);
  // if window is out of bounds, shift to be in bounds
  if(window_left < 0) {
    window_right -= window_left;
    window_left = 0;
  }
  if(window_right > genomeMax) {
    window_left -= (window_right - genomeMax);
    window_right = genomeMax;
  }
  return [window_left, window_right];
}

/*
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
 *  @returns range of renderWindow x-positions for a given proportional range
 */
function getRenderXRangeForFrac() {
  if(!percentScale) return null;
  let [l,r] = calcXBounds();
  let [x1, x2] = renderWindow.map(x => l+x*(r-l));
  return [x1, x2];
}

function getRenderNTRange(genomeID) {
  if(!percentScale) return renderWindow;
  let [l,r] = calcXBounds();
  let [start, end] = getRenderXRangeForFrac().map(x => (x-xDisps[genomeID])/scaleFactor);
  return [clamp(start,0,genomeMax), clamp(end,0,genomeMax)];
}

/*
 *  @returns array [min, max] where
 *    min = x-start of the leftmost genome, max = x-end of the rightmost genome
 */
function calcXBounds() {
  let min = 9*(10**9), max = -9*(10**9);
  for(let g in xDisps) {
    if(xDisps[g] > max) max = xDisps[g];
    if(xDisps[g] < min) min = xDisps[g];
  }
  return [min, max + scaleFactor*genomeMax];
}

/*
 *  [TO BE ADDED TO 'regular' utils.js since it can be used in the inspect page]
 *  @returns array of functional annotation types from table in `genomeData`
 */
function getFunctionalAnnotations() {
  return Object.keys(genomeData.genomes[0][1].genes.functions[0]);
}

/*
 *  [TO BE ADDED TO 'regular' utils.js]
 *  @returns relevant category:color dict from constants.js for a given functional annotation
 */
function getColorDefaults(fn_type) {
  switch(fn_type) {
    case 'COG_CATEGORY':
      return default_COG_colors;
    case 'KEGG_CATEGORY':
      return default_KEGG_colors;
    case 'Source':
      return default_source_colors;
    default:
      // user-supplied color table
      return getCustomColorDict(fn_type);
  }
}

/*
 *  [TO BE ADDED TO 'regular' utils.js]
 *  @returns arbitrary category:color dict given a list of categories
 */
function getCustomColorDict(fn_type) {
  if(!Object.keys(genomeData.genomes[0][1].genes.functions[0]).includes(fn_type)) return null;

  let cags = [];
  genomeData.genomes.forEach(genome => {
    Object.values(genome[1].genes.functions).forEach(fn => {
      let cag = getCagForType(fn, fn_type);
      if(cag && !cags.includes(cag)) cags.push(cag);
    });
  });

  let out = custom_cag_colors.reduce((out, field, index) => {
    out[cags[index]] = field;
    return out;
  }, {});
  delete out["undefined"];

  return out;
}

/*
 *  [TO BE ADDED TO genomeview/UI.js ... OR potentially 'regular' utils.js]
 *  @returns target gene's category code for a given functional annotation type.
 */
function getCagForType(geneFunctions, fn_type) {
  switch(fn_type) {
    case 'COG_CATEGORY':
      return geneFunctions && geneFunctions[fn_type] ? geneFunctions[fn_type][1][0] : null;
    case 'KEGG_CATEGORY':
      // TODO
      return null;
    default:
      let out = geneFunctions != null && geneFunctions[fn_type] != null ? geneFunctions[fn_type][1] : null;
      if(out && out.indexOf(',') != -1) out = out.substr(0,out.indexOf(',')); // take first cag in case of a comma-separated list
      if(out && out.indexOf(';') != -1) out = out.substr(0,out.indexOf(';')); // or semicolon-separated
      return out;
  }
}

/*
 *  [TO BE ADDED TO 'regular' utils.js]
 *  @returns category name corresponding to a given single-character code, for the approporiate functional annotation type
 */
function getCagName(category, fn_type) {
  switch(fn_type) {
    case 'COG_CATEGORY':
      return COG_categories[category];
    case 'KEGG_CATEGORY':
      return KEGG_categories[category];
    default:
      return category;
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