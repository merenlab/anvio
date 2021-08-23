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
    draw();
      return;
    }
    if(newStart < 0) newStart = 0;
    if(newEnd > genomeMax) newEnd = genomeMax;
  
    brush.extent([newStart, newEnd]);
    brush(d3.select(".brush").transition());
    brush.event(d3.select(".brush").transition());
  }