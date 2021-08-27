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
    let start = parseInt($('#brush_start').val()), end = parseInt($('#brush_end').val());
    let newStart, newEnd;
  
    let len = end - start;
    if(len > 4*genomeMax/50) {
      newStart = Math.floor(start + genomeMax/50), newEnd = Math.floor(end - genomeMax/50);
    } else {
      if(len < 50) return;
      newStart = Math.floor(start + len/4);
      newEnd = Math.floor(end - len/4);
      if(newEnd - newStart <= 0) return;
    }
  
    brush.extent([newStart, newEnd]);
    brush(d3.select(".brush").transition());
    brush.event(d3.select(".brush").transition());
  }
  
  function zoomOut() {
    let start = parseInt($('#brush_start').val()), end = parseInt($('#brush_end').val());
  
    let newStart = start - genomeMax/50, newEnd = end + genomeMax/50;
    if(newStart == 0 && newEnd == genomeMax) { // for extra-zoomed-out view
      scaleFactor = 0.01;
      if(dynamicScaleInterval) adjustScaleInterval();
    // draw();
    drawer.draw()
      return;
    }
    if(newStart < 0) newStart = 0;
    if(newEnd > genomeMax) newEnd = genomeMax;
  
    brush.extent([newStart, newEnd]);
    brush(d3.select(".brush").transition());
    brush.event(d3.select(".brush").transition());
  }

  function setPtsPerADL(newResolution) {
    if(isNaN(newResolution)) return;
    newResolution = parseInt(newResolution);
    if(newResolution < 0 || newResolution > genomeMax) {
      alert(`Invalid value, genome spacing must be in range 0-${genomeMax}.`);
      return;
    }
    adlPtsPerLayer = newResolution;
    // draw();
    drawer.draw()
  }
  
  function showAllADLPts() {
    setPtsPerADL(genomeMax);
    $('#showAllADLPtsBtn').blur();
  }
  
  function alignRulers() {
    for(genome of genomeData.genomes) {
      xDisps[genome[0]] = xDisplacement;
    }
    percentScale = false;
    drawScale();
    bindViewportToWindow();
    updateScalePos();
    updateRenderWindow();
    // draw();
    drawer.draw()
    $('#alignRulerBtn').blur();
  }
  
  function setGenomeSpacing(newSpacing) {
    if(isNaN(newSpacing)) return;
    newSpacing = parseInt(newSpacing);
    if(newSpacing < 0 || newSpacing > 1000) {
      alert(`Invalid value, genome spacing must be in range 0-1000.`);
      return;
    }
    spacing = newSpacing;
    // draw();
    drawer.draw()
  }
  
  function setScaleInterval(newScale) {
    if(isNaN(newScale)) return;
    newScale = parseInt(newScale);
    if(newScale < 50) {
      alert(`Invalid value, scale interval must be >=50.`);
      return;
    }
    scaleInterval = newScale;
    // draw();
    drawer.draw()
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

      if(dynamicScaleInterval) adjustScaleInterval();

      // draw();
      drawer.draw()
      let moveToX = percentScale ? getRenderXRangeForFrac()[0] : xDisplacement+scaleFactor*b[0];
      canvas.absolutePan({x: moveToX, y: 0});

      // TODO: restrict min view to 300 NTs? (or e.g. scaleFactor <= 4)
  }
}

/*
 *  Dynamically set scale tick interval based on scaleFactor.
 */
function adjustScaleInterval() {
  let val = Math.floor(100/scaleFactor);
  let roundToDigits = Math.floor(Math.log10(val)) - 1;
  let newInterval = Math.floor(val/(10**roundToDigits)) * (10**roundToDigits);
  scaleInterval = newInterval;
  $('#genome_scale_interval').val(scaleInterval);
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
  }
}