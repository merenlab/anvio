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
 * File Overview : This file contains code related to the functionality of genomeview's scale element and other general navigation items.
 * Scale is great! Scale is dynamic! Scale is responsive.
 */

 function zoomIn() {
    let start = percentScale ? parseFloat($('#brush_start').val()) : parseInt($('#brush_start').val());
    let end = percentScale ? parseFloat($('#brush_end').val()) : parseInt($('#brush_end').val());
    let newStart, newEnd;

    let len = end - start;

    if(percentScale){
      if(len > 0.08) {
        newStart = start+0.02, newEnd = end-0.02;
      } else {
        if(len < 0.02) return;
        newStart = start+len/4;
        newEnd = end-len/4;
        if(newEnd - newStart <= 0) return;
      }
    } else {
      if(len > 4*genomeMax/50) {
        newStart = Math.floor(start + genomeMax/50), newEnd = Math.floor(end - genomeMax/50);
      } else {
        if(len < 50) return;
        newStart = Math.floor(start + len/4);
        newEnd = Math.floor(end - len/4);
        if(newEnd - newStart <= 0) return;
      }
    }

    brush.extent([newStart, newEnd]);
    brush(d3.select(".brush").transition());
    brush.event(d3.select(".brush").transition());
  }

  function zoomOut(type, start, end) {
    let newStart, newEnd

    if(type && type == 'fully'){
      newStart = 0
      newEnd = percentScale ? 1 : genomeMax
    } else if(type && type == 'partial'){
      newStart = start - (percentScale ? .02 : 10000)
      newEnd = end + (percentScale ? .02 : 10000)
    }else {
      if(percentScale) {
        let start = parseFloat($('#brush_start').val());
        let end = parseFloat($('#brush_end').val());
        newStart = start - 0.02, newEnd = end + 0.02;
      } else {
        let start = parseInt($('#brush_start').val());
        let end = parseInt($('#brush_end').val());
        newStart = start - genomeMax/50, newEnd = end + genomeMax/50;
      }
      if(newStart == 0 && newEnd == genomeMax) { // for extra-zoomed-out view
        scaleFactor = 0.01;
        if(dynamicScaleInterval) adjustScaleInterval();
      drawer.draw()
        return;
      }
    }
    if(newStart < 0) newStart = 0;
    if(newEnd > genomeMax) newEnd = percentScale ? 1 : genomeMax;

    brush.extent([newStart, newEnd]);
    brush(d3.select(".brush").transition());
    brush.event(d3.select(".brush").transition());
  }

  async function zoomOutAndWait(type, start, end, time) {
    zoomOut(type, start, end);

    return new Promise(resolve => {
      setTimeout(() => {
        resolve();
      }, time);
    });
  }

  async function goToGene(genomeID, geneID, start, end) {
    await zoomOutAndWait('partial', start, end, 500);
    drawer.glowGenes([{genomeID: genomeID, geneID: geneID}]);
  }

/*
 *  Resets viewport if outside bounds of the view window, with padding on each end
 */
function bindViewportToWindow() {
  let vpt = canvas.viewportTransform;
  let [l,r] = calcXBounds();
  if(vpt[4] > 250 - l) {
    vpt[4] = 250 - l;
  } else if(vpt[4] < canvas.getWidth() - r - 125) {
    vpt[4] = canvas.getWidth() - r - 125;
  }
}

/*
 *  Replaces nt scale with a 0-1 proportional scale
 */
function setPercentScale() {
  percentScale = true;
  drawScale();
}

function drawScale() {
  let scaleWidth = canvas.getWidth();
  let scaleHeight = 100;
  let domain = percentScale ? [0,1] : [0,genomeMax];
  let xScale = d3.scale.linear().range([0, scaleWidth]).domain(domain);
  let scaleAxis = d3.svg.axis()
              .scale(xScale)
              .tickSize(scaleHeight);
  let scaleArea = d3.svg.area()
              .interpolate("monotone")
              .x(function(d) { return xScale(d); })
              .y0(scaleHeight)
              .y1(0);
  brush = d3.svg.brush()
              .x(xScale)
              .on("brushend", onBrush);

  $('#scaleSvg').empty();
  let scaleBox = d3.select("#scaleSvg").append("g")
              .attr("id", "scaleBox")
              .attr("class","scale")
              .attr("y", 230)
              .attr("transform", percentScale ? "translate(10,0)" : "translate(5,0)");

  scaleBox.append("g")
              .attr("id", "scaleMarkers")
              .attr("class", "x axis top noselect")
              .attr("transform", "translate(0,0)")
              .call(scaleAxis);

  scaleBox.append("g")
              .attr("class", "x brush")
              .call(brush)
              .selectAll("rect")
              .attr("y", 0)
              .attr("height", scaleHeight);

  $("#scaleSvg").attr("width", percentScale ? scaleWidth + 20 : scaleWidth + 10);

  function onBrush(){
      var b = brush.empty() ? xScale.domain() : brush.extent();

      if (brush.empty()) {
          $('.btn-selection-sequence').addClass('disabled').prop('disabled', true);
      } else {
          $('.btn-selection-sequence').removeClass('disabled').prop('disabled', false);
      }

      if(!percentScale) b = [Math.floor(b[0]), Math.floor(b[1])];

      $('#brush_start').val(b[0]);
      $('#brush_end').val(b[1]);

      let ntsToShow = b[1] - b[0];
      scaleFactor = percentScale ? canvas.getWidth()/(ntsToShow*genomeMax) : canvas.getWidth()/ntsToShow;
      updateRenderWindow();

      if(dynamicScaleInterval) drawer.adjustScaleInterval();

      drawer.draw()
      let moveToX = percentScale ? getRenderXRangeForFrac()[0] : xDisplacement+scaleFactor*b[0];
      canvas.absolutePan({x: moveToX, y: 0});

      // TODO: restrict min view to 300 NTs? (or e.g. scaleFactor <= 4)
  }
}

/*
 *  Update scale info to match viewport location.
 */
function updateScalePos() {
  let [newStart, newEnd] = percentScale ? getFracForVPT() : getNTRangeForVPT();
  brush.extent([newStart, newEnd]);
  brush(d3.select(".brush").transition());
  $('#brush_start').val(newStart);
  $('#brush_end').val(newEnd);
}

function updateRenderWindow() {
  if(percentScale) {
    let resolution = 4; // # decimals to show for renderw window
    let [start, end] = [parseFloat($('#brush_start').val()), parseFloat($('#brush_end').val())];
    let diff = end - start > 0.1 ? Math.round(10**resolution*(end - start) / 2)/10**resolution : 0.05;
    renderWindow = [clamp(start - diff, 0, 1), clamp(end + diff, 0, 1)];
  } else {
    let [start, end] = [parseInt($('#brush_start').val()), parseInt($('#brush_end').val())];
    let diff = end - start > 10000 ? Math.floor((end - start)/2) : 5000;
    renderWindow = [clamp(start - diff, 0, genomeMax), clamp(end + diff, 0, genomeMax)];

    if(filter_gene_colors_to_window) {
      generateColorTable(null, color_db);
      // TODO: filter to window for percent scale, too
    }
  }
}

/*
 * Pan viewport to the first gene in the target gene cluster.
 *
 * @param gc : target gene cluster ID
 * @returns tuple [a,b] where
 *  a is genomeID of the first genome containing `gc` and
 *  b is NT position of the middle of the target gene
*/
function viewCluster(gc) {
  if(!genomeData.gene_associations["anvio-pangenome"]) return;

  let genes = [];
  let geneMid;
  let first = true;
  let firstGenomeID;

  if(!gc || gc in genomeData.gene_associations["anvio-pangenome"]["gene-cluster-name-to-genomes-and-genes"]) {
    for(genome of genomeData.genomes) {
      var targetGenes = getGenesOfGC(genome[0], gc);
      if(targetGenes == null) continue;
      var targetGeneID = targetGenes[0]; /* TODO: implementation for multiple matching gene IDs */
      var targetGene = genome[1].genes.gene_calls[targetGeneID];
      genes.push({genomeID:genome[0],geneID:targetGeneID});
      if(first) {
        geneMid = targetGene.start + (targetGene.stop - targetGene.start) / 2;
        canvas.absolutePan({x: scaleFactor*geneMid + xDisps[genome[0]] - canvas.getWidth()/2, y: 0});
        canvas.viewportTransform[4] = clamp(canvas.viewportTransform[4], canvas.getWidth() - genomeMax*scaleFactor - xDisps[genome[0]] - 125, 125);
        firstGenomeID = genome[0];
        first = false;
      }
    }
    updateScalePos();
    updateRenderWindow();
    drawer.draw();
    drawer.glowGenes(genes, true);
    return (first ? null : [firstGenomeID, geneMid]);
  } else {
    console.log('Warning: ' + gc + ' is not a gene cluster in data structure');
    return null;
  }
}
