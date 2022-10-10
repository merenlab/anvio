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
      let mainCanvasHeight =  spacing * settings['genomeData']['genomes'].length + additionalSpacing + (spacing * maxGroupSize * groupLayerPadding)
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
 *  @returns array of functional annotation types from table in `genomeData`
 */
function getFunctionalAnnotations() {
  return Object.keys(genomeData.genomes[0][1].genes.functions[0]);
}

/*
 *  @returns arbitrary category:color dict given a list of categories
 */
function getCustomColorDict(fn_type, cags=null, order=null) {
  if(!Object.keys(genomeData.genomes[0][1].genes.functions[0]).includes(fn_type)) return null;

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
  thresh_count_gene_colors = thresh;
  generateColorTable(null, $("#gene_color_order").val());
  drawer.draw();
}

function showSaveStateWindow(){
    $.ajax({
        type: 'GET',
        cache: false,
        url: '/state/all',
        success: function(state_list) {
            $('#saveState_list').empty();

            for (let state_name in state_list) {
                var _select = "";
                if (state_name == current_state_name)
                {
                    _select = ' selected="selected"';
                }
                $('#saveState_list').append('<option ' + _select + '>' + state_name + '</option>');
            }

            $('#modSaveState').modal('show');
            if ($('#saveState_list').val() === null) {
                $('#saveState_name').val('default');
            } else {
                $('#saveState_list').trigger('change');
            }
        },
        error: function(error){
          console.log('got an error', error)
        }
    });
}

function showLoadStateWindow(){
    $.ajax({
        type: 'GET',
        cache: false,
        url: '/state/all',
        success: function(state_list) {
            $('#loadState_list').empty();

            for (let state_name in state_list) {
                $('#loadState_list').append('<option lastmodified="' + state_list[state_name]['last_modified'] + '">' + state_name + '</option>');
            }

            $('#modLoadState').modal('show');
        }
    });
}

function saveState()
{
    var name = $('#saveState_name').val();

    if (name.length==0) {
        $('#saveState_name').focus();
        return;
    }

    var state_exists = false;

    $.ajax({
        type: 'GET',
        cache: false,
        async: false,
        url: '/state/all',
        success: function(state_list) {
            for (let state_name in state_list) {
                if (state_name == name)
                {
                    state_exists = true;
                }
            }

        }
    });

    if (state_exists && !confirm('"' + name + '" already exist, do you want to overwrite it?')) {
        return;
    }

    $.ajax({
        type: 'POST',
        cache: false,
        url: '/state/save/' + name,
        data: {
            'content': JSON.stringify(serializeSettings())
        },
        success: function(response) {
            if (typeof response != 'object') {
                response = JSON.parse(response);
            }

            if (response['status_code']==0)
            {
                toastr.error("Failed, Interface running in read only mode.");
            }
            else if (response['status_code']==1)
            {
                // successfull
                $('#modSaveState').modal('hide');

                current_state_name = name;
                toastr.success("State '" + current_state_name + "' successfully saved.");
            }
        }
    });
}

function toggleRightPanel(name) {
  ['#mouseover_panel', '#settings-panel', '#query-panel'].forEach(function(right_panel) {
      if (right_panel == name)
          return;

      $(right_panel).hide();
  });
  console.log(name);
  $(name).toggle();

  if ($('#mouseover_panel').is(':visible')) {
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
      $('#toggle-panel-settings').html('&#9664;');
  }

  if ($('#query-panel').is(':visible')) {
      $('#toggle-panel-query').addClass('toggle-panel-query-pos');
      $('#toggle-panel-query-inner').html('&#9658;');
  } else {
      $('#toggle-panel-query').removeClass('toggle-panel-query-pos');
      $('#toggle-panel-query-inner').html('&#9664;');
  }
}

