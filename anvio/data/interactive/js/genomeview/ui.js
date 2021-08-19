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
 * File Overview : This file contains functions related to building + updating UI elements and responding to user interaction with those elements. 
 * As a general rule, processes that invoke jQuery should probably live here. 
 */

function showToolTip(event){
  $('#tooltip-body').show().append(`
    <p></p>
    <style type="text/css">
      .tftable {font-size:12px;color:#333333;width:100%;border-width: 1px;border-color: #729ea5;border-collapse: collapse;}
      .tftable th {font-size:12px;background-color:#acc8cc;border-width: 1px;padding: 8px;border-style: solid;border-color: #729ea5;text-align:left;}
      .tftable tr {background-color:#d4e3e5;}
      .tftable td {font-size:12px;border-width: 1px;padding: 8px;border-style: solid;border-color: #729ea5;}
      .tftable tr:hover {background-color:#ffffff;}
    </style>

    <table class="tftable" border="1">
      <tr><th>Data</th><th>Value</th></tr>
      <tr><td>Split</td><td>${event.target.gene.contig}</td></tr>
      <tr><td>Start in Contig</td><td>${event.target.gene.start}</td></tr>
      <tr><td>Length</td><td>${event.target.gene.stop - event.target.gene.start}</td></tr>
      <tr><td>Gene Callers ID</td><td>${event.target.geneID}</td></tr>
      <tr><td>Gene Cluster</td><td>${genomeData.gene_associations["anvio-pangenome"] ? genomeData.gene_associations["anvio-pangenome"]["genome-and-gene-names-to-gene-clusters"][event.target.genomeID][event.target.geneID] : "None"}</td></tr>
    </table>
    <button>some action</button>
    <button>some other action</button>
  `).css({'position' : 'absolute', 'left' : event.e.clientX, 'top' : event.e.clientY })
}

function toggleSettingsPanel() {
  $('#settings-panel').toggle();

  if ($('#settings-panel').is(':visible')) {
    $('#toggle-panel-settings').addClass('toggle-panel-settings-pos');
    $('#toggle-panel-settings-inner').html('&#9658;');
  } else {
    $('#toggle-panel-settings').removeClass('toggle-panel-settings-pos');
    $('#toggle-panel-settings-inner').html('&#9664;');
  }
}

function buildGenomesTable(genomes, order){
  genomes.map(genome => {
    var height = '50';
    var margin = '15';
    var template = '<tr id={genomeLabel}>' +
                  '<td><img src="images/drag.gif" class="drag-icon" id={genomeLabel} /></td>' +
                  '<td> {genomeLabel} </td>' +
                  '<td>n/a</td>' +
                  '<td>n/a</td>' +
                  '<td>n/a</td>' +
                  '<td><input class="input-height" type="text" size="3" id="height{id}" value="{height}"></input></td>' +
                  '<td class="column-margin"><input class="input-margin" type="text" size="3" id="margin{id}" value="{margin}"></input></td>' +
                  '<td>n/a</td>' +
                  '<td>n/a</td>' +
                  '<td><input type="checkbox" class="layer_selectors"></input></td>' +
                  '</tr>';
    let genomeLabel= Object.keys(genome[1]['contigs']['info']);

    template = template.replace(new RegExp('{height}', 'g'), height)
                        .replace(new RegExp('{margin}', 'g'), margin)
                        .replace(new RegExp('{genomeLabel}', 'g'), genomeLabel);

    $('#tbody_genomes').append(template);
  })
  
  $("#tbody_genomes").sortable({helper: fixHelperModified, handle: '.drag-icon', items: "> tr"}).disableSelection();

  $("#tbody_genomes").on("sortupdate", (event, ui) => {
    changeGenomeOrder($("#tbody_genomes").sortable('toArray'))
  })
}
  
function buildGroupLayersTable(layerLabel){
  var height = '50';
  var margin = '25';
  var template = '<tr id={layerLabel}>' +
                  '<td><img src="images/drag.gif" class="drag-icon" id={layerLabel} /></td>' +
                  '<td> {layerLabel} </td>' +
                  '<td><div id="{layerLabel}_color" style="margin-left: 5px;" class="colorpicker" style="background-color: #FFFFFF" color="#FFFFFF"></div></td>' +
                  '<td>n/a</td>' +
                  '<td>n/a</td>' +
                  '<td><input type="checkbox" class="additional_selectors" id={layerLabel}-show onclick="toggleAdditionalDataLayer(event)" checked=true></input></td>' +
                  '</tr>';
  template = template.replace(new RegExp('{height}', 'g'), height)
                      .replace(new RegExp('{margin}', 'g'), margin)
                      .replace(new RegExp('{layerLabel}', 'g'), layerLabel);
  $('#tbody_additionalDataLayers').append(template);
  $("#tbody_additionalDataLayers").sortable({helper: fixHelperModified, handle: '.drag-icon', items: "> tr"}).disableSelection();

  $("#tbody_additionalDataLayers").on("sortupdate", (event, ui) => {
    changeGroupLayersOrder($("#tbody_additionalDataLayers").sortable('toArray'))
  })
}
  
function toggleAdditionalDataLayer(e){
  let layer = e.target.id.split('-')[0]

  if(e.target.checked){
    stateData['display']['additionalDataLayers'][layer] = true
    maxGroupSize += 1
  } else {
    stateData['display']['additionalDataLayers'][layer] = false
    maxGroupSize -= 1 // decrease group height if hiding the layer
  }
  draw()
}