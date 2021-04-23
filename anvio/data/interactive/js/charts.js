/**
 * Javascript library to visualize anvi'o charts
 *
 *  Author: A. Murat Eren <a.murat.eren@gmail.com>
 *  Credits: Özcan Esen, Gökmen Göksel, Tobias Paczian.
 *  Copyright 2015, The anvio Project
 *
 * This file is part of anvi'o (<https://github.com/meren/anvio>).
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

var request_prefix = getParameterByName('request_prefix');
var VIEWER_WIDTH = window.innerWidth || document.documentElement.clientWidth || document.getElementsByTagName('body')[0].clientWidth;

var current_state_name = "";

var layers;
var coverage;
var variability;
var maxVariability = 0;
var maxCountOverCoverage = 0;
var indels;
var geneParser;
var contextSvg;
var state;
var layers_ordered;
var visible_layers;
var contig_id;
var page_header;
var highlight_gene;
var gene_mode;
var sequence;
var charts;
var brush;
var inspect_mode;
var highlightBoxes;
var indels_enabled;
var show_nucleotides = true;
var maxNucleotidesInWindow = 300;
var minNucleotidesInWindow = 30;
var gene_offset_y = 0;
var select_boxes = {};
var curr_height;
var show_cags_in_split = true;

var mcags;

var cog_annotated = false, kegg_annotated = false;

function loadAll() {
    console.log("Initiated: " + Math.round(Date.now()/1000));
    $.ajaxPrefilter(function(options) {
        if (request_prefix) {
            options.url = request_prefix + options.url;
            if (options.type.toLowerCase() == 'post')
            {
                options.data += '&csrfmiddlewaretoken=' + getCookie('csrftoken');
            }
        }
        return options;
    });

    $('.colorpicker').colpick({
        layout: 'hex',
        submit: 0,
        colorScheme: 'light',
        onChange: function(hsb, hex, rgb, el, bySetColor) {
            $(el).css('background-color', '#' + hex);
            $(el).attr('color', '#' + hex);

            if (!bySetColor) $(el).val(hex);
        }
    }).keyup(function() {
        $(this).colpickSetColor(this.value);
    });

    contig_id = getParameterByName('id');
    highlight_gene = getParameterByName('highlight_gene') == 'true';
    gene_mode = getParameterByName('gene_mode') == 'true';

    if (typeof localStorage.state === 'undefined')
    {
        state = {}
    }
    else
    {
        state = JSON.parse(localStorage.state);
    }

    if(state['snvs_enabled'] == null) {
        state['snvs_enabled'] = getParameterByName('show_snvs') == 'true';
    }

    var endpoint = (gene_mode ? 'charts_for_single_gene' : 'charts');

    $.ajax({
            type: 'POST',
            cache: false,
            url: '/data/' + endpoint + '/' + state['order-by'] + '/' + contig_id,
            data: {'state': JSON.stringify(state)},
            success: function(contig_data) {
                state = contig_data['state'];
                page_header = contig_data.title;
                layers = contig_data.layers;
                coverage = contig_data.coverage;
                sequence = contig_data.sequence;
                variability = [];
                indels = [];

                for (var i=0; i<coverage.length; i++) {
                    variability[i] = [];
                    for (var l=0; l<4; l++) {
                        variability[i][l] = [];
                        for (var h=0; h<coverage[i].length; h++) {
                            if (contig_data.variability[i][l].hasOwnProperty(h)) {
                                variability[i][l].push(contig_data.variability[i][l][h]);
                                if (state['snvs_enabled'] && contig_data.variability[i][l][h] > maxVariability) {
                                    maxVariability = contig_data.variability[i][l][h];
                                }
                            } else {
                                variability[i][l].push(0);
                            }
                        }
                    }
                }

                competing_nucleotides = contig_data.competing_nucleotides;
                indels = contig_data.indels;

                for(var i=0; i<indels.length; i++) {
                  var ikeys = Object.keys(indels[i]);
                  for(var j=0; j<ikeys.length; j++) {
                    let ccVal = indels[i][ikeys[j]]["count"]/indels[i][ikeys[j]]["coverage"];
                    if(ccVal > maxCountOverCoverage) maxCountOverCoverage = ccVal;
                    if(maxCountOverCoverage >= 1) {
                      maxCountOverCoverage = 1;
                      i = indels.length;
                      break;
                    }
                  }
                }
                console.log("Checking SNV and indels tables from state data (done: " + Math.round(Date.now()/1000) + ")");

                previous_contig_name = contig_data.previous_contig_name;
                next_contig_name = contig_data.next_contig_name;
                index = contig_data.index;
                total = contig_data.total;
                genes = contig_data.genes;

                if(layers.length == 0){
                    console.log('Warning: no layers returned')
                }

                next_str = " | next &gt;&gt;&gt;";
                prev_str = "&lt;&lt;&lt; prev | ";
                position = index + " of " + total;

                // anvi-server uses iframes for prettier urls, links need to be open _top
                var target_str = '';

                if (self != top) {
                    target_str = 'target="_top"';
                }

                inspect_mode = 'inspect';

                if (gene_mode) {
                    inspect_mode = 'inspect_gene';
                }
                else if (highlight_gene) {
                    inspect_mode = 'inspect_context';
                }

                if(next_contig_name)
                    next_str = '<a onclick="localStorage.state = JSON.stringify(state);" href="' + generate_inspect_link({'type': inspect_mode, 'item_name': next_contig_name, 'show_snvs': state['snvs_enabled']}) +'" '+target_str+'> | next &gt;&gt;&gt;</a>';

                if(previous_contig_name)
                    prev_str = '<a onclick="localStorage.state = JSON.stringify(state);" href="' + generate_inspect_link({'type': inspect_mode, 'item_name': previous_contig_name, 'show_snvs': state['snvs_enabled']}) + '" '+target_str+'>&lt;&lt;&lt; prev | </a>';

                $('#window-title').html("anvi-inspect: " + page_header);
                $('#header').append("<strong>" + page_header + "</strong> detailed");
                $('#split-settings').append("<p style='text-align: center'>" + prev_str + position + next_str + "</p>");
                $('#split-settings').append("<p style='text-align: center; margin-top: -10px; margin-bottom: 15px;'><a href='#' onclick='showSearchItemsDialog();'>Select or Search Item</a></p>");

                $('#range-box').append(`<div style="display: inline-block; margin-bottom:10px;" class="form-inline"> \
                                        Selection range from
                                            <input class="form-control input-sm" id="brush_start" type="text" value="0" size="5">
                                        to
                                            <input class="form-control input-sm" id="brush_end" type="text" value="${sequence.length}" size="5">\
                                    </div>`);

                geneParser = new GeneParser(genes);
                geneParser["data"].forEach(function(gene) {
                  if(gene.functions != null) {
                    if(gene.functions.hasOwnProperty("COG20_CATEGORY")) {
                      gene.functions["COG_CATEGORY"] = gene.functions["COG20_CATEGORY"];
                      gene.functions["COG_FUNCTION"] = gene.functions["COG20_FUNCTION"];
                    } else if(gene.functions.hasOwnProperty("COG14_CATEGORY")) {
                      gene.functions["COG_CATEGORY"] = gene.functions["COG14_CATEGORY"];
                      gene.functions["COG_FUNCTION"] = gene.functions["COG14_FUNCTION"];
                    }

                    if(gene.functions.hasOwnProperty("COG_CATEGORY")) cog_annotated = true;
                    if(gene.functions.hasOwnProperty("KEGG_Class")) kegg_annotated = true;
                    if(cog_annotated && kegg_annotated) return;
                  }
                });
                console.log("Checking for gene functional annotations (done: " + Math.round(Date.now()/1000) + ")");

                if(!state['highlight-genes']) state['highlight-genes'] = {};
                state['large-indel'] = 10;
                $("#largeIndelInput").val(state['large-indel']);
                state['min-indel-coverage'] = 0;
                //$("#minIndelInput").val(state['min-indel-coverage']);

                if(state['show_snvs'] == null) {
                  state['show_snvs'] = state['snvs_enabled']!=null ? (state['snvs_enabled'] && maxVariability!=0) : false;
                }
                indels_enabled = maxCountOverCoverage != 0;
                if(!indels_enabled || state['show_indels'] == null) state['show_indels'] = indels_enabled;

                // adjust menu options
                if(!indels_enabled && (!state['snvs_enabled'] || maxVariability==0)) {
                  $('#toggleSNVIndelTable').hide();
                  $("#indels").hide();
                  $('#settings-section-info-SNV-warning').append("Note: SNVs and indels are disabled for this split.");
                  $('#settings-section-info-SNV-warning').show();
                } else {
                  if(!indels_enabled) {
                    $('#indels').hide();
                    $('#indels_picker').hide();
                    $('#settings-section-info-SNV-warning').append("Note: indels are disabled for this split.");
                    $('#settings-section-info-SNV-warning').show();
                  }
                  if(!state['snvs_enabled'] || maxVariability==0) {
                    $('#snv_picker').hide();
                    $('#settings-section-info-SNV-warning').append("Note: SNVs are disabled for this split.");
                    $('#settings-section-info-SNV-warning').show();
                  }
                }

                // create function color menu and table; set default color states
                $('#gene_color_order').append($('<option>', {
                  value: 'Source',
                  text: 'Source'
                }));
                if(!state.hasOwnProperty('source-colors')) {
                  state['source-colors'] = default_source_colors;
                }
                generateFunctionColorTable(state['source-colors'], "Source", highlight_genes=state['highlight-genes'], show_cags_in_split);
                mcags = Object.keys(default_source_colors);
                if(cog_annotated) {
                  $('#gene_color_order').append($('<option>', {
                    value: 'COG',
                    text: 'COG'
                  }));

                  if(!state.hasOwnProperty('cog-colors')) {
                    state['cog-colors'] = default_COG_colors
                  }
                }
                if(kegg_annotated) {
                  $('#gene_color_order').append($('<option>', {
                    value: 'KEGG',
                    text: 'KEGG'
                  }));

                  if(!state.hasOwnProperty('kegg-colors')) {
                    state['kegg-colors'] = default_KEGG_colors
                  }
                }

                // show SNVs and indels?
                if(state['snvs_enabled']) {
                  let numSNVs = 0;
                  for(var i = 0; i < competing_nucleotides.length; i++) {
                    for(var key in competing_nucleotides[i]) {
                      if(competing_nucleotides[i].hasOwnProperty(key)) numSNVs++;
                    }
                  }
                  if(state['show_snvs'] && numSNVs > 1000) {
                    state['show_snvs'] = false;
                    $("div.snvs-disabled").append("WARNING: A total of " + numSNVs + " SNVs were dedected on this page and are not shown to optimize perfomance. Use the settings panel to show them.");
                    $("div.snvs-disabled").fadeIn(300);
                  }
                }
                if(indels_enabled) {
                  let numIndels = 0;
                  for(var i = 0; i < indels.length; i++) {
                    for(var key in indels[i]) {
                      if(indels[i].hasOwnProperty(key)) numIndels++;
                    }
                  }
                  if(state['show_indels'] && numIndels > 1000) {
                    state['show_indels'] = false;
                    $("div.indels-disabled").append("WARNING: A total of " + numIndels + " INDELs were dedected on this page and are not shown to optimize perfomance. Use the settings panel to show them.");
                    $("div.indels-disabled").fadeIn(300);
                  }
                }

                if(state['show_snvs']) $('#toggle_snv_box').attr("checked", "checked");
                if(state['show_indels']) $('#toggle_indel_box').attr("checked", "checked");
                $('#toggle_highlight_box').attr("checked", "checked");
                $('#toggle_nucl_box').attr("checked", "checked");

                createCharts(state);
                console.log("Creating charts (done: " + Math.round(Date.now()/1000) + ")");
                $('.loading-screen').hide();

                // on initial load from main interface
                if(state['state-name'] != current_state_name) {
                  $.ajax({
                          type: 'GET',
                          cache: false,
                          url: '/state/get/' + state['state-name'],
                          success: function(response) {
                              try{
                                  clusteringData = response[1]['data'];
                                  loadOrderingAdditionalData(response[1]);
                                  processState(state['state-name'], response[0]);
                                  console.log("Processing state data from the server (done: " + Math.round(Date.now()/1000) + ")");
                              }catch(e){
                                  console.error("Exception thrown", e.stack);
                                  toastr.error('Failed to parse state data, ' + e);
                                  defer.reject();
                                  return;
                              }
                              waitingDialog.hide();
                          }
                  });
                }

                $('#brush_start, #brush_end').keydown(function(ev) {
                    if (ev.which == 13) {
                        let start = parseInt($('#brush_start').val());
                        let end = parseInt($('#brush_end').val());

                        if (!isNumber(start) || !isNumber(end) || start < 0 || start > sequence.length || end < 0 || end > sequence.length) {
                            alert(`Invalid value, value needs to be in range 0-${sequence.length}.`);
                            return;
                        }

                        if (start >= end) {
                            alert('Starting value cannot be greater or equal to the ending value.');
                            return;
                        }

                        brush.extent([start, end]);
                        brush(d3.select(".brush").transition());
                        brush.event(d3.select(".brush").transition());
                    }
                });

                document.body.addEventListener("keydown", function(ev) {
                    if(ev.which == 83) { // S = 83
                      toggleSettingsPanel();
                    } else if(!($('#brush_start', '#brush_end').is(':focus')) && (ev.which == 37 || ev.which == 39)) {
                      let start = parseInt($('#brush_start').val());
                      let end = parseInt($('#brush_end').val());

                      if(ev.which == 37 && ev.shiftKey) {
                        start--;
                        end--;
                        if(start < 0) return;
                      } else if(ev.which == 39 && ev.shiftKey) {
                        start++;
                        end++;
                        if(end > sequence.length) return;
                      } else if(ev.which == 37) { // Left arrow = 37
                        start -= 3;
                        end -= 3;
                        if(start < 0) return;
                      } else if(ev.which == 39) { // Right arrow = 39
                        start += 3;
                        end += 3;
                        if(end > sequence.length) return;
                      }

                      $('#brush_start').val(start);
                      $('#brush_end').val(end);

                      brush.extent([start, end]);
                      brush(d3.select(".brush").transition());
                      brush.event(d3.select(".brush").transition());

                      display_nucleotides();
                      drawArrows(start, end, $('#gene_color_order').val(), gene_offset_y, Object.keys(state['highlight-genes']));
                    }
                });

                $('#gene_color_order').on('focus', function() {
                    Object.keys(state['highlight-genes']).forEach(gene_id => {
                      state['highlight-genes'][gene_id] = $('#picker_' + gene_id).attr('color');
                    });

                }).change(function() {
                    state['gene-fn-db'] = $(this).val();
                    switch($(this).val()) {
                      case "COG":
                        mcags = Object.keys(COG_categories);
                        break;
                      case "KEGG":
                        mcags = Object.keys(KEGG_categories);
                        break;
                      case "Source":
                        mcags = Object.keys(default_source_colors);
                    }
                    resetFunctionColors(state[$(this).val().toLowerCase() + '-colors']);
                    redrawArrows();
                    $(this).blur();
                });

                $('#largeIndelInput').on('keydown', function(e) {
                  if(e.keyCode == 13) { // 13 = enter key
                    if($(this).val() < 1 || $(this).val() > 9999) {
                      alert("Invalid value, value needs to be in range 1-9999");
                    } else {
                      state['large-indel'] = $(this).val();
                      $(this).blur();
                    }
                  }
                });

                /*$('#minIndelInput').on('keydown', function(e) {
                  if(e.keyCode == 13) { // 13 = enter key
                    if($(this).val() < 0 || $(this).val() > 9999) {
                      alert("Invalid value, value needs to be in range 0-9999");
                    } else {
                      state['min-indel-coverage'] = $(this).val();
                      $(this).blur();
                    }
                  }
                });*/

                $('#toggle_snv_box').on('change', function() {
                  waitingDialog.show('Drawing ...',
                      {
                          dialogSize: 'sm',
                          onShow: function() {
                              toggleSNVs();
                              waitingDialog.hide();
                          },
                      });
                  if($('div.snvs-disabled').length > 0) $('div.snvs-disabled').remove();
                });
                $('#toggle_indel_box').on('change', function() {
                  waitingDialog.show('Drawing ...',
                      {
                          dialogSize: 'sm',
                          onShow: function() {
                              toggleIndels();
                              waitingDialog.hide();
                          },
                      });
                  if($('div.indels-disabled').length > 0) $('div.indels-disabled').remove();
                });
                $('#toggle_highlight_box').on('change', function() {
                  toggleHighlightBoxes();
                });
                $('#toggle_nucl_box').on('change', function() {
                  toggle_nucleotide_display();
                });
                console.log("Setting event listeners (done: " + Math.round(Date.now()/1000) + ")");
            }
        });
}

function drawHighlightBoxes() {
  var nucl_shown = $("#DNA_sequence").length > 0;

  var width = VIEWER_WIDTH * .80;

  var start = $('#brush_start').val(), end = $('#brush_end').val();

  highlightBoxes.attr("height", curr_height);
  $("#highlightBoxesSvg").empty();
  var nBoxes = nucl_shown ? end - start : 100;
  var endpts = nucl_shown ? getGeneEndpts(start, end) : [];
  for(var i = 0; i < nBoxes; i++) {
    highlightBoxes.append("rect")
                  .attr("id", "highlight_" + i)
                  .attr("class", "highlightbox")
                  .attr("x", i*(width/nBoxes))
                  .attr("width", (width/nBoxes))
                  .attr("height", curr_height)
                  .attr("fill", endpts.includes(i) ? "red" : "#989898")
                  .attr("fill-opacity", 0)
                  .attr("transform", "translate(50,20)");
  }
  console.log("Drawing vertical highlight boxes (done: " + Math.round(Date.now()/1000) + ")");
}

function drawAAHighlightBoxes() {
  var endpts = getGeneEndpts($('#brush_start').val(), $('#brush_end').val());

  $('#context-container').on('mouseover', function(e) {
    if(!e.target.id.startsWith("AA_")) return;
    var box_num = parseFloat(get_box_id_for_AA(e.target, "highlight_").substring(10));
    var marked = endpts.includes(box_num) || endpts.includes(box_num+2);
    $('#highlight_' + box_num + ', #highlight_' + (box_num+1) + ', #highlight_' + (box_num+2)).attr('fill-opacity', 0.25).attr('fill', marked ? 'red' : '#989898');
  }).mouseout(function(e) {
    if(!e.target.id.startsWith("AA_")) return;
    var box_num = parseFloat(get_box_id_for_AA(e.target, "highlight_").substring(10));
    $('#highlight_' + box_num + ', #highlight_' + (box_num+1) + ', #highlight_' + (box_num+2)).attr('fill-opacity', 0);

    if(endpts.includes(box_num)) {
      $('#highlight_' + (box_num+1) + ', #highlight_' + (box_num+2)).attr('fill', '#989898');
    } else if(endpts.includes(box_num+2)) {
      $('#highlight_' + box_num + ', #highlight_' + (box_num+1)).attr('fill', '#989898');
    }
  });
  console.log("Drawing amino acid highlight boxes (done: " + Math.round(Date.now()/1000) + ")");
}

function get_box_id_for_AA(aa, id_start) {
  var aa_x = Math.round($(aa).attr("x")*100)/100;
  var id = $('[id^=' + id_start + ']').filter(function() {
    return Math.round($(this).attr("x")*100)/100 == aa_x;
  }).attr("id");
  return id;
}

/*
 *  Generates KEGG color table html given a color palette.
 *
 *  Params:
 *   - fn_colors:         a dictionary matching each category to a hex color code
 *   - fn_type:           a string indicating function category type: one of "COG", "KEGG", "Source"
 *   - highlight_genes:   a dictionary matching specific gene caller IDs to hex colors to override other coloring
 *   - filter_to_split: if true, filters categories to only those shown in the split
 */
function generateFunctionColorTable(fn_colors, fn_type, highlight_genes={}, filter_to_split) {
  var db = (function(type){
    switch(type) {
      case "COG":
        return COG_categories;
      case "KEGG":
        return KEGG_categories;
      case "Source":
        return default_source_colors;
      default:
        console.log("Warning: Invalid type for function color table");
        return null;
    }
  })(fn_type);
  if(db == null) return;

  $('#tbody_function_colors').empty();

  if(fn_type != "Source" && filter_to_split) {
    var new_colors = {};
    for(gene of genes) {
      if(gene.functions && gene.functions[fn_type + "_CATEGORY"]) {
        new_colors[gene.functions[fn_type + "_CATEGORY"][0][0]] = fn_colors[gene.functions[fn_type + "_CATEGORY"][0][0]]
      }
    }
    if(Object.keys(new_colors).length > 0) fn_colors = new_colors;
  }

  Object.keys(db).forEach(function(category){
    if(!(category in fn_colors)) {
      return;
    }
    var category_name = (fn_type == "Source" ? category : db[category]);

    var tbody_content =
     '<tr id="picker_row_' + category + '"> \
        <td></td> \
        <td> \
          <div id="picker_' + category + '" class="colorpicker" color="' + fn_colors[category] + '" background-color="' + fn_colors[category] + '" style="background-color: ' + fn_colors[category] + '; margin-right:16px; margin-left:16px"></div> \
        </td> \
        <td>' + category_name + '</td> \
      </tr>';

    $('#tbody_function_colors').append(tbody_content);
  });

  $('.colorpicker').colpick({
      layout: 'hex',
      submit: 0,
      colorScheme: 'light',
      onChange: function(hsb, hex, rgb, el, bySetColor) {
          $(el).css('background-color', '#' + hex);
          $(el).attr('color', '#' + hex);

          state[$('#gene_color_order').val().toLowerCase() + '-colors'][el.id.substring(7)] = '#' + hex;
          if (!bySetColor) $(el).val(hex);
      }
  }).keyup(function() {
      $(this).colpickSetColor(this.value);
  });

  if(!isEmpty(highlight_genes)) {

    for(gene of genes) {
      var gene_callers_id = "" + gene.gene_callers_id;
      if(Object.keys(highlight_genes).includes(gene_callers_id)) {
        var tbody_content =
         '<tr id="picker_row_' + gene_callers_id + '"> \
            <td></td> \
            <td> \
              <div id="picker_' + gene_callers_id + '" class="colorpicker" color="' + highlight_genes[gene_callers_id] + '" background-color="' + highlight_genes[gene_callers_id] + '" style="background-color: ' + highlight_genes[gene_callers_id] + '; margin-right:16px; margin-left:16px"></div> \
            </td> \
            <td>Gene ID: ' + gene_callers_id + '</td> \
          </tr>';

        $('#tbody_function_colors').prepend(tbody_content);

        $('#picker_' + gene_callers_id).colpick({
            layout: 'hex',
            submit: 0,
            colorScheme: 'light',
            onChange: function(hsb, hex, rgb, el, bySetColor) {
                $(el).css('background-color', '#' + hex);
                $(el).attr('color', '#' + hex);

                state['highlight-genes'][el.id.substring(7)] = '#' + hex;
                if (!bySetColor) $(el).val(hex);
            }
        }).keyup(function() {
            $(this).colpickSetColor(this.value);
        });
      }
    }

  }
  console.log("Generating gene functional annotation color table (done: " + Math.round(Date.now()/1000) + ")");
}

function toggleSNVs() {
  state['show_snvs'] = !state['show_snvs'];
  createCharts(state);
  console.log("Toggling SNV markers (done: " + Math.round(Date.now()/1000) + ")");
}

function toggleIndels() {
  state['show_indels'] = !state['show_indels'];
  createCharts(state);
  console.log("Toggling indel markers (done: " + Math.round(Date.now()/1000) + ")");
}

function toggleHighlightBoxes() {
  if(state['show_highlights']) {
    $('#highlightBoxesSvg').empty();
    $('#highlight-boxes').css('pointer-events', 'none');
    $('#context-container').off('mouseover mouseout');
  } else {
    drawHighlightBoxes();
    if($('#DNA_sequence').length == 1) drawAAHighlightBoxes();
    $('#highlight-boxes').css('pointer-events', 'all');
  }
  state['show_highlights'] = !state['show_highlights'];
  console.log("Togging highlight boxes (done: " + Math.round(Date.now()/1000) + ")");
}

function toggleGeneIDColor(gene_id, color="#FF0000") {
  if(gene_id in state['highlight-genes']) {
    removeGeneIDColor(gene_id);
  } else {
    addGeneIDColor(gene_id, color);
  }
}

function addGeneIDColor(gene_id, color="#FF0000") {
  state['highlight-genes'][gene_id] = color;
  var tbody_content =
   '<tr id="picker_row_' + gene_id + '"> \
      <td></td> \
      <td> \
        <div id="picker_' + gene_id + '" class="colorpicker" color="' + color + '" background-color="' + color + '" style="background-color: ' + color + '; margin-right:16px; margin-left:16px"></div> \
      </td> \
      <td>Gene ID: ' + gene_id + '</td> \
    </tr>';

  $('#tbody_function_colors').prepend(tbody_content);

  $('#picker_' + gene_id).colpick({
      layout: 'hex',
      submit: 0,
      colorScheme: 'light',
      onChange: function(hsb, hex, rgb, el, bySetColor) {
          $(el).css('background-color', '#' + hex);
          $(el).attr('color', '#' + hex);

          state['highlight-genes'][el.id.substring(7)] = '#' + hex;
          if (!bySetColor) $(el).val(hex);
      }
  }).keyup(function() {
      $(this).colpickSetColor(this.value);
  });

  // define arrow marker
  contextSvg.select('#contextSvgDefs')
      .append('svg:marker')
      .attr('id', 'arrow_' + gene_id)
      .attr('markerHeight', 2)
      .attr('markerWidth', 2)
      .attr('orient', 'auto')
      .attr('refX', 0)
      .attr('refY', 0)
      .attr('viewBox', '-5 -5 10 10')
      .append('svg:path')
        .attr('d', 'M 0,0 m -5,-5 L 5,0 L -5,5 Z')
        .attr('fill', $('#picker_' + gene_id).attr('color'));

  redrawArrows();
}

function removeGeneIDColor(gene_id) {
  delete state['highlight-genes'][gene_id];
  $('#picker_row_' + gene_id).remove();
  contextSvg.select('#contextSvgDefs').select('#arrow_' + gene_id).remove();

  redrawArrows();
}

function redrawArrows() {
  resetArrowMarkers();
  drawArrows(parseInt($('#brush_start').val()), parseInt($('#brush_end').val()), $('#gene_color_order').val(), gene_offset_y, Object.keys(state['highlight-genes']));
  console.log("Redrawing gene arrows (done: " + Math.round(Date.now()/1000) + ")");
}

function resetArrowMarkers() {
  $('#contextSvgDefs').empty();

  ["none"].concat(mcags).concat(Object.keys(state['highlight-genes'])).forEach(function(category){
    contextSvg.select('#contextSvgDefs')
        .append('svg:marker')
        .attr('id', 'arrow_' + category)
        .attr('markerHeight', 2)
        .attr('markerWidth', 2)
        .attr('orient', 'auto')
        .attr('refX', 0)
        .attr('refY', 0)
        .attr('viewBox', '-5 -5 10 10')
        .append('svg:path')
          .attr('d', 'M 0,0 m -5,-5 L 5,0 L -5,5 Z')
          .attr('fill', category != "none" ? $('#picker_' + category).attr('color') : "gray");
  });
  console.log("Resetting arrow markers (done: " + Math.round(Date.now()/1000) + ")");
}

/*
 * Sets gene colors for the selected category type to the default set
 *
 * Params:
 * - fn_colors: If set, resets state to this dictionary instead of the defaults.
 */
function resetFunctionColors(fn_colors=null) {
  if($('#gene_color_order') == null) return;

  switch($('#gene_color_order').val()) {
    case 'Source':
      Object.assign(state['source-colors'], fn_colors ? fn_colors : default_source_colors);
      break;
    case 'COG':
      Object.assign(state['cog-colors'], fn_colors ? fn_colors : default_COG_colors);
      break;
    case 'KEGG':
      Object.assign(state['kegg-colors'], fn_colors ? fn_colors : default_KEGG_colors);
      break;
  }

  generateFunctionColorTable(state[$('#gene_color_order').val().toLowerCase() + '-colors'],
                             $('#gene_color_order').val(),
                             state['highlight-genes'],
                             show_cags_in_split);
  console.log("Resetting functional annotation colors (done: " + Math.round(Date.now()/1000) + ")");
}

function toggleShowCagsInSplit() {
  show_cags_in_split = !show_cags_in_split;
  generateFunctionColorTable(state[$('#gene_color_order').val().toLowerCase() + '-colors'],
                             $('#gene_color_order').val(),
                             state['highlight-genes'],
                             show_cags_in_split);
}

function toggle_nucleotide_display() {
  show_nucleotides = !show_nucleotides;
  if(show_nucleotides) {
    display_nucleotides();
  } else {
    contextSvg.select("#DNA_sequence").remove();
    contextSvg.select("#AA_sequence").remove();
    contextSvg.select("#solids").remove();
    contextSvg.attr("height", 150);
    $("#gene-chart").attr("transform", "translate(50, 10)");
    $("#context-chart").attr("transform", "translate(50, 80)");
    $("#gene-arrow-chart").attr("transform", "translate(50, -10)");
    gene_offset_y = 0;
    drawHighlightBoxes();
    $('#context-container').off('mouseover mouseout');
  }
  console.log("Toggling nucleotide display (done: " + Math.round(Date.now()/1000) + ")");
}

/*
 * Sequence styling inspired by the integrated genomics viewer:
 * http://software.broadinstitute.org/software/igv/
 */
function display_nucleotides() {
  if(!show_nucleotides) return;

  contextSvg.select("#DNA_sequence").remove();
  contextSvg.select("#AA_sequence").remove();
  contextSvg.select("#solids").remove();

  var margin = {top: 20, right: 50, bottom: 150, left: 50};
  var width = VIEWER_WIDTH * .80;

  let start = parseInt($('#brush_start').val());
  let end = parseInt($('#brush_end').val());

  if(end - start > maxNucleotidesInWindow || end - start < minNucleotidesInWindow) {
    contextSvg.attr("height", 150);
    $("#gene-chart").attr("transform", "translate(50, 10)");
    $("#context-chart").attr("transform", "translate(50, 80)");
    $("#gene-arrow-chart").attr("transform", "translate(50, -10)");
    gene_offset_y = 0;
    $('#context-container').off('mouseover mouseout');
    return;
  }

  let curSeq = sequence.substring(start, end);

  var nucleotideDefs = contextSvg.append('defs');
  var linearGradient = nucleotideDefs.append('linearGradient')
                                     .attr('id', 'solids')
                                     .attr('x1', "0")
                                     .attr('y1', "0")
                                     .attr('x2', "1")
                                     .attr('y2', "0");

  // define nucleotide color gradient
  for(var i = 0; i < curSeq.length; i++) {
    var txtColor;
    switch(curSeq[i]) {
      case "A":
        txtColor = "rgb(0,144,0)";
        break;
      case "G":
        txtColor = "rgb(208,104,7)";
        break;
      case "T":
        txtColor = "rgb(255,41,26)";
        break;
      case "C":
        txtColor = "rgb(53,48,220)";
        break;
      default:
        txtColor = "black";
    }

    linearGradient.append('svg:stop')
                  .attr('offset', "" + i/curSeq.length)
                  .attr('style', "stop-color:" + txtColor + ";stop-opacity:1");

    if(i < curSeq.length-1) {
      linearGradient.append('svg:stop')
                    .attr('offset', "" + (i+1)/curSeq.length)
                    .attr('style', "stop-color:" + txtColor + ";stop-opacity:1");
    } else {
      linearGradient.append('svg:stop')
                    .attr('offset', "1")
                    .attr('style', "stop-color:" + txtColor + ";stop-opacity:1");
    }
  }

  /* width of monospaced character per font size */
  var nucl_text_font = width/((end-start)*.6002738402061856);

  var nucl_sequence = contextSvg.append("text")
                                .text(sequence.substring(start, end))
                                .attr("id", "DNA_sequence")
                                .attr("fill", "url(#solids)")
                                .attr("class", "noselect")
                                .attr("font-size", nucl_text_font);
  var dna_seq_height = contextSvg.select("#DNA_sequence")[0][0].getBBox().height;
  var nucl_text_y = .75*dna_seq_height;
                   nucl_sequence.attr("y", nucl_text_y)
                                .attr("font-family", "monospace")
                                .attr("transform", "translate(" + (margin.left) + ", 10)");

  var show_AAs = false;
  geneParser["data"].forEach(function(gene){
    if(gene.start_in_split < end-2 && gene.stop_in_split > start+2) show_AAs = true;
  });

  if(show_AAs) {
    var aa_sequence = contextSvg.append("g")
                                .attr("id", "AA_sequence")
                                .attr("class", "noselect")
                                .attr('transform', 'translate(50, 10)');

    var textWidth = width/(end-start);
    var prev_gene_stop = 0;
    var offset_y = 0;
    var overlapping_genes = false;

    geneParser["data"].sort((a,b) => a.start_in_split > b.start_in_split? 1 : -1).forEach(function(gene){
      if(gene.start_in_split < end-2 && gene.stop_in_split > start+2) {
        display_AA(gene);
      }
    });

    function display_AA(gene) {
      var aa_string = "";
      var rect_x = 0;
      var buffer, aa_i, aa_f;
      if(gene.start_in_split > start) {
        buffer = gene.start_in_split - start;
        aa_i = 0;
      } else {
        buffer = (gene.stop_in_split - start) % 3;
        aa_i = Math.ceil((start-gene.start_in_split)/3);
      }
      aa_f = Math.floor((Math.min(end, gene.stop_in_split) - gene.start_in_split)/3);
      rect_x += textWidth * buffer;

      if(gene.start_in_split < prev_gene_stop) {
        offset_y = offset_y == 0 ? 5+dna_seq_height : 0;
        overlapping_genes = true;
      }
      if(gene.call_type == 1 && typeof gene.aa_sequence != "undefined") {
        var aas_in_window = gene.aa_sequence.substring(aa_i, aa_f);
        aa_string += "\xa0" + aas_in_window.split('').join('\xa0\xa0') + "\xa0";
        if(aa_f == gene.aa_sequence.length+1) {
          aa_string += "STP";
          aas_in_window += "\xa0";
        }
        for(var i = 0; i < aas_in_window.length; i++) {
          aa_sequence.append("rect")
                     .attr("id", "AA_" + i)
                     .attr("height", dna_seq_height)
                     .attr("width", 3*textWidth)
                     .attr("x", rect_x)
                     .attr("y", dna_seq_height + offset_y)
                     .attr("fill", i % 2 == 0 ? "rgb(144,137,250)" : "rgb(81,68,211)");
          rect_x += (3*textWidth);
        }

        aa_sequence.append("text")
                  .text(aa_string)
                  .attr('font-size', nucl_text_font)
                  .attr("font-family", "monospace")
                  .attr("fill", "white")
                  .attr("x", buffer*textWidth)
                  .attr("y", dna_seq_height + .67*dna_seq_height + offset_y);

        if(gene.direction == "r") {
          aa_sequence.append("rect")
                     .attr("height", 1.1*dna_seq_height)
                     .attr("width", aas_in_window.length*3*textWidth)
                     .attr("x", buffer*textWidth)
                     .attr("y", dna_seq_height + offset_y)
                     .attr("fill", "red")
                     .attr("fill-opacity", 0.1);
        }
      }
      prev_gene_stop = gene.stop_in_split;
    }
  }

  var extra_y = dna_seq_height +
                (show_AAs? dna_seq_height : 0) +
                (overlapping_genes? 5+dna_seq_height : 0);

  contextSvg.attr("height", 150 + extra_y);

  // reposition gene arrow chart
  $("#gene-chart").attr("transform", "translate(50, " + (20+extra_y) + ")");
  $("#context-chart").attr("transform", "translate(50, " + (90+extra_y) + ")");
  $("#gene-arrow-chart").attr("transform", "translate(50, " + (10+extra_y) + ")");
  gene_offset_y = 10+extra_y;

  if(state['show_highlights']) {
    drawHighlightBoxes();
    drawAAHighlightBoxes();
  }
  console.log("Drawing nucleotides (done: " + Math.round(Date.now()/1000) + ")");
}


function show_selected_sequence() {
    let range = charts[0].xScale.domain();

    // who knows?
    range[0] = Math.max(range[0], 0);
    range[1] = Math.min(range[1], sequence.length);

    show_sequence_modal(`Sequence [${range[0]}, ${range[1]}]`,
        `${page_header} range:${range[0]},${range[1]}\n` + sequence.substring(range[0], range[1]));
}

function computeGCContent(window_size, step_size) {
    let gc_array = [];
    let padding = parseInt(window_size / 2);

    if (step_size > 0) {
        for (let i=padding; i < sequence.length - padding; i = i + step_size) {
            let gc_count = 0;

            for (let j=i; j < i + window_size; j++) {
                let pos = j - padding;

                if (sequence[pos] == 'C' || sequence[pos] == 'G' || sequence[pos] == 'c' || sequence[pos] == 'g') {
                    gc_count++;
                }
            }

            gc_array.push(gc_count / window_size);
        }
    }

    return gc_array;
    console.log("Computing GC content (done: " + Math.round(Date.now()/1000) + ")");
}


function showOverlayGCContentDialog() {
    if (typeof sessionStorage.gc_overlay_settings !== 'undefined') {
        let gc_overlay_settings = JSON.parse(sessionStorage.gc_overlay_settings);

        $('#gc_window_size').val(gc_overlay_settings['gc_window_size']);
        $('#gc_step_size').val(gc_overlay_settings['gc_step_size']);

        $('#gc_overlay_color').attr('color', gc_overlay_settings['gc_overlay_color']);
        $('#gc_overlay_color').css('background-color', gc_overlay_color['gc_overlay_color']);

    }

    $('#GCContentOverlayDialog').modal('show');
}


function applyOverlayGCContent() {
    let gc_overlay_settings = {
        'gc_window_size': $('#gc_window_size').val(),
        'gc_step_size': $('#gc_step_size').val(),
        'gc_overlay_color': $('#gc_overlay_color').attr('color')
    }

    sessionStorage.gc_overlay_settings = JSON.stringify(gc_overlay_settings);
    createCharts(state);
    console.log("Applying GC content overlay (done: " + Math.round(Date.now()/1000) + ")");
}


function resetOverlayGCContent() {
    delete sessionStorage.gc_overlay_settings;
    createCharts(state);
    console.log("Resetting GC content overlay (done: " + Math.round(Date.now()/1000) + ")");
}


function showSetMaxValuesDialog() {
    var table = '<table class="table borderless"><thead class="thead-light"><tr><th>Sample</th><th>Max Coverage</th><th>Limit Max Coverage</th></tr></thead><tbody>';

    var max_coverage_values;
    var has_max_coverage_values = (typeof sessionStorage.max_coverage !== 'undefined');
    if (has_max_coverage_values) {
        max_coverage_values = JSON.parse(sessionStorage.max_coverage);
    }

    var j=0;
    for (i in layers_ordered) {
        var layer_name = layers_ordered[i];
        var layer_index = layers.indexOf(layer_name);

        if (!(state['layers'].hasOwnProperty(layer_name) && parseFloat(state['layers'][layer_name]['height']) == 0)) {
            var max_val
            var actual_max_val = GetMaxMin(coverage[layer_index])['Max'];
            if (has_max_coverage_values) {
                max_val = max_coverage_values[j];
            } else {
                max_val = 0;
            }

            table += '<tr> \
                        <td>' + layer_name + '</td> \
                        <td><a href="#" onclick="$(\'#max_multiple\').val(\'' + actual_max_val + '\')">' + actual_max_val + '</a></td> \
                        <td style="text-align: center;"><input class="form-control input-sm max-coverage-input" type="text" size="5" value="' + max_val + '"/></td> \
                      </tr>';

            j++;
        }
    }

    table += '<tr> \
                <td>Set Max for all samples:</td> \
                <td>&nbsp;</td> \
                <td style="text-align: center;"><div class="input-group">\
                    <input class="form-control input-sm" id="max_multiple" type="text" size="5" value="0"/> \
                        <span class="input-group-btn"> \
                            <button type="button" class="btn btn-default btn-sm" onclick="$(\'.max-coverage-input\').val($(\'#max_multiple\').val());">Set</button> \
                        </span> \
                    </div> \
                </td> \
            </tr>';


    $('#setMaxValuesDialog .modal-body').empty().append(table + '</tbody></table>');
    $('#setMaxValuesDialog').modal('show');
}

function applyMaxValues() {
    var max_values = []
    $('#setMaxValuesDialog .modal-body tbody tr').each(function(index, row) {
        max_values.push(parseInt($(row).find('td:last input').val()));
    });

    sessionStorage.max_coverage = JSON.stringify(max_values);
    createCharts(state);
    console.log("Applying max values (done: " + Math.round(Date.now()/1000) + ")");
}


function resetMaxValues() {
    delete sessionStorage.max_coverage;
    createCharts(state);
    console.log("Resetting max values (done: " + Math.round(Date.now()/1000) + ")");
}


function showSearchItemsDialog() {
    $('#searchItemDialog').modal('show');
    $('#searchItemDialog .search-results').empty();
}

function search_items(search_query, page) {
    if (typeof page === 'undefined') {
        page = 0;
    }

    // anvi-server uses iframes for prettier urls, links need to be open _top
    var target_str = '';

    if (self != top) {
        target_str = 'target="_top"';
    }

    $('#searchItemDialog .search-results').empty();

    $.ajax({
            type: 'POST',
            cache: false,
            url: '/data/search_items',
            data: {'search-query': search_query, 'page': page},
            success: function(data) {
                page = parseInt(data['page']);
                search_query = data['search-query']
                let total_page = parseInt(data['total_page']);

                let results = data['results'];
                let results_html = '';

                for (let i=0; i < results.length; i++) {
                    let item_name = results[i];
                    let item_name_pretty = item_name;

                    if (search_query.length > 0) {
                        let begin = item_name.toLowerCase().indexOf(search_query.toLowerCase());
                        item_name_pretty = [item_name.slice(0, begin),
                                           '<mark>',
                                           item_name.slice(begin, begin + search_query.length),
                                           '</mark>',
                                           item_name.slice(begin + search_query.length, item_name.length)
                                           ].join("");

                    }

                    let link = '<a onclick="localStorage.state = JSON.stringify(state);" href="' + generate_inspect_link({'type': inspect_mode, 'item_name': item_name, 'show_snvs': state['snvs_enabled']}) +'" '+target_str+'>' + item_name_pretty + '</a>';
                    results_html += link + '<br />';

                }

                results_html += '<br /><br /><center>';

                if (results.length > 0) {
                    if (page > 0) {
                        results_html += `<a href="#" onclick="search_items('${search_query}', ${page-1});">&lt;&lt;&lt; prev</a> | `;
                    }

                    results_html += " page " + (page + 1) + " of " + total_page;

                    if (page + 1 < total_page) {
                        results_html += ` | <a href="#" onclick="search_items('${search_query}', ${page+1});"> next &gt;&gt;&gt;</a>`;
                    }
                }
                else
                {
                    results_html+="<b>No results found.</b>";
                }

                results_html += '</center>';

                $('#searchItemDialog .search-results').append(results_html);
            }
        });
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

function showLoadStateWindow()
{
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

function showSaveStateWindow()
{
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
        }
    });
}

function loadState()
{
    var defer = $.Deferred();
    $('#modLoadState').modal('hide');
    if ($('#loadState_list').val() == null) {
        defer.reject();
        return;
    }

    var state_name = $('#loadState_list').val();
    waitingDialog.show('Requesting state data from the server ...',
        {
            dialogSize: 'sm',
            onHide: function() {
                defer.resolve();
            },
            onShow: function() {
                $.ajax({
                        type: 'GET',
                        cache: false,
                        url: '/state/get/' + state_name,
                        success: function(response) {
                            try{
                                clusteringData = response[1]['data'];
                                loadOrderingAdditionalData(response[1]);
                                processState(state_name, response[0]);
                            }catch(e){
                                console.error("Exception thrown", e.stack);
                                toastr.error('Failed to parse state data, ' + e);
                                defer.reject();
                                return;
                            }
                            waitingDialog.hide();
                        }
                    });
            },
        }
    );

    return defer.promise();
}

function loadOrderingAdditionalData(order) {
    collapsedNodes = [];

    if (order.hasOwnProperty('additional')) {
        let orders_additional = order['additional'];

        if (typeof orders_additional === 'string') {
            orders_additional = JSON.parse(orders_additional);
        }

        if (orders_additional.hasOwnProperty('collapsedNodes')) {
            collapsedNodes = orders_additional['collapsedNodes'];
        }
    }
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
            'content': JSON.stringify(serializeSettings(), null, 4)
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

 /*
  *  This function takes the state passed from the main interface and
  *  updates only the function colors and indel variables.
  */
function processState(state_name, state) {
    // set color defaults
    if(!state['cog-colors']) state['source-colors'] = default_source_colors;
    if(!state['cog-colors']) state['cog-colors'] = default_COG_colors;
    if(!state['cog-colors']) state['kegg-colors'] = default_KEGG_colors;

    if(JSON.parse(localStorage.state) && JSON.parse(localStorage.state)['gene-fn-db']) {
      state['gene-fn-db'] = JSON.parse(localStorage.state)['gene-fn-db'];
    } else {
      state['gene-fn-db'] = 'Source';
    }
    $('#gene_color_order').val(state['gene-fn-db']);

    generateFunctionColorTable(state[state['gene-fn-db'].toLowerCase() + '-colors'], state['gene-fn-db'], highlight_genes=state['highlight-genes'], show_cags_in_split);
    mcags = Object.keys(state[state['gene-fn-db'].toLowerCase() + '-colors']);
    this.state = state;

    if(!state['highlight-genes']) {
      state['highlight-genes'] = {};
    }
    if(!state['show_highlights']) state['show_highlights'] = true;

    // define arrow markers for highlighted gene ids
    Object.keys(state['highlight-genes']).forEach(function(gene_id){
      contextSvg.select('#contextSvgDefs')
          .append('svg:marker')
          .attr('id', 'arrow_' + gene_id)
          .attr('markerHeight', 2)
          .attr('markerWidth', 2)
          .attr('orient', 'auto')
          .attr('refX', 0)
          .attr('refY', 0)
          .attr('viewBox', '-5 -5 10 10')
          .append('svg:path')
            .attr('d', 'M 0,0 m -5,-5 L 5,0 L -5,5 Z')
            .attr('fill', $('#picker_' + gene_id).attr('color'));
    });
    redrawArrows();

    if(state['show_indels']) {
      if(!state.hasOwnProperty('large-indel')) state['large-indel'] = 10;
      if(!state.hasOwnProperty('min-indel-coverage')) state['large-indel'] = 0;
      $('#largeIndelInput').val(20);
      //$('#minIndelInput').val(0);
    }

    state['state-name'] = current_state_name = state_name;

    toastr.success("State '" + current_state_name + "' successfully loaded.");
}

 /*
  *  This function saves only function colors to the state passed
  *  from the main interface.
  */
function serializeSettings() {
    var state = this.state;

    if(!isEmpty(state['highlight-genes'])) {
      Object.keys(state['highlight-genes']).forEach(function(gene_id){
        state['highlight-genes'][gene_id] = $('#picker_' + gene_id).attr('color');
      });
    }

    state['large-indel'] = $("#largeIndelInput").val();
    //state['min-indel-coverage'] = $("#minIndelInput").val();

    return state;
}

function createCharts(state){
    /* Adapted from Tyler Craft's Multiple area charts with D3.js article:
    http://tympanus.net/codrops/2012/08/29/multiple-area-charts-with-d3-js/  */
    $('#chart-container, #context-container, #highlight-boxes').empty();

    if (state['current-view'] == "single"){
        // if we are working with a non-merged single profile, we need to do some ugly hacks here,
        // simply because the sample name does not appear among 'layers' found in the state variable.
        layers_ordered = layers;
        state['layers'][layers[0]] = state['layers']['mean_coverage'];
    } else if (state['current-view'] == "blank_view") {
        layers_ordered = layers;
        state['layers'][layers[0]] = state['layers']['length'];
    } else {
        // this is the usual path for merged profiles:
        layers_ordered = state['layer-order'].filter(function (value) { if (layers.indexOf(value)>-1) return true; return false; });
    }

    visible_layers = 0;
    for (i in layers_ordered)
    {
      var layer_id = layers_ordered[i];

      if (parseFloat(state['layers'][layer_id]['height']) > 0)
        visible_layers++;
    }

    var margin = {top: 20, right: 50, bottom: 150, left: 50};
    var width = VIEWER_WIDTH * .80;
    var chartHeight = 200;
    var height = ((chartHeight + 10) * visible_layers);
    curr_height = height;
    var contextHeight = 50;
    var contextWidth = width;

    var svg = d3.select("#chart-container").append("svg")
            .attr("width", width + margin.left + margin.right)
            .attr("height", (height + margin.top + margin.bottom));

    $('#chart-container').css("width", (width + 150) + "px");
    $('#chart-container').css("height", height + "px");

    $('#SNV-boxes').empty();
    var snvBoxesSvg = d3.select("#SNV-boxes").append("svg")
                            .attr("id", "SNVBoxesSvg")
                            .attr("width", width + margin.left + margin.right)
                            .attr("height", height + margin.top);
    $('#SNV-boxes').css("width", (width + 150) + "px");
    $('#SNV-boxes').css("height", height + "px");
    $('#SNV-boxes').css("top", (margin.top - 20) + "px");


    charts = [];

    var layersCount = layers.length;

    coverage.forEach(function(d) {
        for (var prop in d) {
            if (d.hasOwnProperty(prop)) {
                d[prop] = parseFloat(d[prop]);
            }
        }
    });

    var max_coverage_values;
    var has_max_coverage_values = (typeof sessionStorage.max_coverage !== 'undefined');
    if (has_max_coverage_values) {
        max_coverage_values = JSON.parse(sessionStorage.max_coverage);
    }
    console.log("Plotting coverage (done: " + Math.round(Date.now()/1000) + ")");

    let gc_content_array = [];
    let gc_overlay_color = '#00FF00';
    let gc_content_window_size = 100;
    let gc_content_step_size = 10;

    if (typeof sessionStorage.gc_overlay_settings !== 'undefined') {
        let gc_overlay_settings = JSON.parse(sessionStorage.gc_overlay_settings);
        gc_content_window_size = parseInt(gc_overlay_settings['gc_window_size']);
        gc_content_step_size = parseInt(gc_overlay_settings['gc_step_size']);
        gc_content_array = computeGCContent(gc_content_window_size, gc_content_step_size);
        gc_overlay_color = gc_overlay_settings['gc_overlay_color'];
    }

    var j=0;
    for(var i = 0; i < layersCount; i++){
        var layer_index = layers.indexOf(layers_ordered[i]);

        if (parseFloat(state['layers'][layers_ordered[i]]['height']) == 0) {
            continue;
        }

        charts.push(new Chart({
                        name: layers[layer_index],
                        coverage: coverage[layer_index],
                        max_coverage: (has_max_coverage_values) ? max_coverage_values[j] : 0,
                        variability_a: variability[layer_index][0],
                        variability_b: variability[layer_index][1],
                        variability_c: variability[layer_index][2],
                        variability_d: variability[layer_index][3],
                        competing_nucleotides: competing_nucleotides[layer_index],
                        indels: indels[layer_index],
                        gc_content: gc_content_array,
                        'gc_content_window_size': gc_content_window_size,
                        'gc_content_step_size': gc_content_step_size,
                        'gc_overlay_color': gc_overlay_color,
                        id: j++,
                        width: width,
                        height: chartHeight,
                        maxVariability: maxVariability,
                        maxCountOverCoverage: maxCountOverCoverage,
                        svg: svg,
                        snv_svg: snvBoxesSvg,
                        margin: margin,
                        showBottomAxis: (j == visible_layers - 1),
                        color: state['layers'][layers[layer_index]]['color']
                }));

    }

    contextSvg = d3.select("#context-container").append("svg")
            .attr("width", width + margin.left + margin.right)
            .attr("height", 150);

    highlightBoxes = d3.select("#highlight-boxes").append("svg")
                                                  .attr("id", "highlightBoxesSvg")
                                                  .attr("width", width + margin.left + margin.right)
                                                  .attr("height", height);
    $('#highlight-boxes').css("width", (width + 150) + "px");

    var defs = contextSvg.append('svg:defs')
                         .attr('id', 'contextSvgDefs');

    contextSvg.append("rect")
       .attr("id", "gene-chart")
       .attr("width", width)
       .attr("height", "60px")
       .attr("fill", "black")
       .attr("fill-opacity", "0.2")
       .attr('transform', 'translate(50, 10)');

    // Define arrow markers
    if(cog_annotated) {
      ["none"].concat(Object.keys(COG_categories)).forEach(function(category){
        defs.append('svg:marker')
            .attr('id', 'arrow_' + category )
            .attr('markerHeight', 2)
            .attr('markerWidth', 2)
            .attr('orient', 'auto')
            .attr('refX', 0)
            .attr('refY', 0)
            .attr('viewBox', '-5 -5 10 10')
            .append('svg:path')
              .attr('d', 'M 0,0 m -5,-5 L 5,0 L -5,5 Z')
              .attr('fill', category != "none" ? $('#picker_' + category).attr('background-color') : "gray");
      });
    }
    if(kegg_annotated) {
      ["none"].concat(Object.keys(KEGG_categories)).forEach(function(category){
        defs.append('svg:marker')
            .attr('id', 'arrow_' + category )
            .attr('markerHeight', 2)
            .attr('markerWidth', 2)
            .attr('orient', 'auto')
            .attr('refX', 0)
            .attr('refY', 0)
            .attr('viewBox', '-5 -5 10 10')
            .append('svg:path')
              .attr('d', 'M 0,0 m -5,-5 L 5,0 L -5,5 Z')
              .attr('fill', category != "none" ? $('#picker_' + category).attr('background-color') : "gray");
      });
    }
    ["none"].concat(Object.keys(default_source_colors)).concat(Object.keys(state['highlight-genes'])).forEach(function(category){
      defs.append('svg:marker')
          .attr('id', 'arrow_' + category )
          .attr('markerHeight', 2)
          .attr('markerWidth', 2)
          .attr('orient', 'auto')
          .attr('refX', 0)
          .attr('refY', 0)
          .attr('viewBox', '-5 -5 10 10')
          .append('svg:path')
            .attr('d', 'M 0,0 m -5,-5 L 5,0 L -5,5 Z')
            .attr('fill', category != "none" ? $('#picker_' + category).attr('color') : "gray");
    });

    $('#context-container').css("width", (width + 150) + "px");

    /* Context down below */
    var contextXScale = d3.scale.linear().range([0, contextWidth]).domain(charts[0].xScale.domain());

    var contextAxis = d3.svg.axis()
                .scale(contextXScale)
                .tickSize(contextHeight);

    var contextArea = d3.svg.area()
                .interpolate("monotone")
                .x(function(d) { return contextXScale(d); })
                .y0(contextHeight)
                .y1(0);

    brush = d3.svg.brush()
                .x(contextXScale)
                .on("brushend", onBrush);

    var context = contextSvg.append("g")
                .attr("id", "context-chart")
                .attr("class","context")
                .attr("transform", "translate(" + (margin.left) + ", 80)");

    context.append("g")
                .attr("class", "x axis top noselect")
                .attr("transform", "translate(0,0)")
                .call(contextAxis)

    context.append("g")
                .attr("class", "x brush")
                .call(brush)
                .selectAll("rect")
                .attr("y", 0)
                .attr("height", contextHeight);

    $('#brush_start').val(contextXScale.domain()[0]);
    $('#brush_end').val(contextXScale.domain()[1]);

    if(show_nucleotides) display_nucleotides();
    if(state['show_highlights']) drawHighlightBoxes();

    function onBrush(){
        /* this will return a date range to pass into the chart object */
        var b = brush.empty() ? contextXScale.domain() : brush.extent();

        if (brush.empty()) {
            $('.btn-selection-sequence').addClass('disabled').prop('disabled', true);
        } else {
            $('.btn-selection-sequence').removeClass('disabled').prop('disabled', false);
        }

        b = [Math.floor(b[0]), Math.floor(b[1])];

        $('#brush_start').val(b[0]);
        $('#brush_end').val(b[1]);

        // rescale nucleotide display
        if(show_nucleotides) display_nucleotides();

        // highlight bars as % of page width without nucleotide display
        if($("#DNA_sequence").length == 0) drawHighlightBoxes();

        for(var i = 0; i < layersCount; i++){
            charts[i].showOnly(b);
        }
        resetArrowMarkers();
        drawArrows(b[0], b[1], $('#gene_color_order').val(), gene_offset_y, Object.keys(state['highlight-genes']));
    }

    drawArrows(0, charts[0].xScale.domain()[1], $('#gene_color_order').val(), gene_offset_y, Object.keys(state['highlight-genes']));
}


function Chart(options){
    this.coverage = options.coverage;
    this.max_coverage = options.max_coverage;
    this.variability_a = options.variability_a;
    this.variability_b = options.variability_b;
    this.variability_c = options.variability_c;
    this.variability_d = options.variability_d;
    this.competing_nucleotides = options.competing_nucleotides;
    this.indels = options.indels;
    this.gc_content = options.gc_content;
    this.gc_content_window_size = options.gc_content_window_size;
    this.gc_content_step_size = options.gc_content_step_size;
    this.gc_overlay_color = options.gc_overlay_color;
    this.width = options.width;
    this.height = options.height;
    this.maxVariability = options.maxVariability;
    this.maxCountOverCoverage = options.maxCountOverCoverage;
    this.svg = options.svg;
    this.snv_svg = options.snv_svg;
    this.id = options.id;
    this.name = options.name;
    this.margin = options.margin;
    this.showBottomAxis = options.showBottomAxis;
    this.color = options.color;

    var localName = this.name;
    var num_data_points = this.variability_a.length;

    this.xScale = d3.scale.linear()
                            .range([0, this.width])
                            .domain([0, this.coverage.length]);


    let cov_min_max = GetMaxMin(this.coverage);
    this.minCoverage = cov_min_max['Min'];
    this.maxCoverage = cov_min_max['Max'];

    // this.max_coverage comes from options, 0 means not available
    if (this.max_coverage == 0) {
        this.maxCoverageForyScale = Math.max(20, this.maxCoverage);
    } else {
        this.maxCoverageForyScale = this.max_coverage;
    }

    let gc_min_max = GetMaxMin(this.gc_content);
    this.maxGCContent = gc_min_max['Min'];
    this.minGCContent = gc_min_max['Max'];

    let yScaleMax = Math.max(this.maxVariability, this.maxCountOverCoverage);

    this.yScale = d3.scale.linear()
                            .range([this.height,0])
                            .domain([0,this.maxCoverageForyScale]);

    this.yScaleLine = d3.scale.linear()
                            .range([this.height, 0])
                            .domain([0, yScaleMax]);

    yScaleLineReverse = d3.scale.linear()
                            .range([0, this.height])
                            .domain([0, yScaleMax]);

    this.yScaleGC = d3.scale.linear()
                            .range([this.yScale(this.minCoverage),
                                   (this.maxCoverage < this.maxCoverageForyScale) ? this.yScale(this.maxCoverage) : 0])
                            .domain([this.minGCContent, this.maxGCContent]);

    var xS = this.xScale;
    var yS = this.yScale;
    var ySL = this.yScaleLine;
    var ySLR = yScaleLineReverse;
    var yGC = this.yScaleGC;

    this.area = d3.svg.area()
                            .x(function(d, i) { return xS(1+i); })
                            .y0(this.height)
                            .y1(function(d) { return (yS(d) < 0) ? 0 : yS(d); });
    if(indels_enabled) {
      this.line = d3.svg.line()
                              .x(function(d, i) { return xS(1+i)+4; })
                              .y(function(d, i) { if(i == 0) return ySL(0); if(i == num_data_points - 1) return ySL(0); return ySL(d); })
                              .interpolate('step-before');
    }

    this.reverseLine = d3.svg.line()
                            .x(function(d, i) { return xS(1+i); })
                            .y(function(d, i) { if(i == 0) return ySLR(0); if(i == num_data_points - 1) return ySLR(0); return ySLR(d); })
                            .interpolate('step-before');


    // .x() needs to stay as a arrow function, it has reference to scope.
    this.gc_line = d3.svg.line()
                            .x((d, i) => { return xS((this.gc_content_window_size / 2) + ((1+i) * this.gc_content_step_size)); })
                            .y(function(d) { return (yGC(d) < 0) ? 0 : yGC(d); });

    /*
        Assign it a class so we can assign a fill color
        And position it on the page
    */
    this.chartContainer = this.svg.append("g")
                        .attr('class',this.name.toLowerCase())
                        .attr("transform", "translate(" + this.margin.left + "," + (this.margin.top + (this.height * this.id) + (10 * this.id)) + ")");

    this.lineContainer = this.svg.append("g")
                        .attr('class',this.name.toLowerCase())
                        .attr("transform", "translate(" + this.margin.left + "," + (this.margin.top + (this.height * this.id) + (10 * this.id)) + ")");

    this.textContainer = this.snv_svg.append("g")
                        .attr('class',this.name.toLowerCase())
                        .attr("transform", "translate(" + this.margin.left + "," + (this.margin.top + (this.height * this.id) + (10 * this.id)) + ")");

    this.textContainerIndels = this.snv_svg.append("g")
                              .attr('class',this.name.toLowerCase())
                              .attr("transform", "translate(" + this.margin.left + "," + (this.margin.top + (this.height * this.id) + (10 * this.id)) + ")");

    this.gcContainer   = this.svg.append("g")
                        .attr('class',this.name.toLowerCase())
                        .attr("transform", "translate(" + this.margin.left + "," + (this.margin.top + (this.height * this.id) + (10 * this.id)) + ")");

    /* Add both into the page */
    this.chartContainer.append("path")
                              .data([this.coverage])
                              .attr("class", "chart")
                              .style("fill", this.color)
                              .style("fill-opacity", "0.5")
                              .attr("d", this.area);

    this.gcContainer.append("path")
                      .data([this.gc_content])
                      .attr("class", "line")
                      .style("stroke", this.gc_overlay_color)
                      .style("stroke-width", "1")
                      .style("fill", "none")
                      .attr("d", this.gc_line);

    if (state['show_snvs']) {
        this.lineContainer.append("path")
            .data([this.variability_b])
            .attr("class", "line")
            .attr("name", "first_pos")
            .style("fill", '#990000')
            .attr("d", this.reverseLine);

        this.lineContainer.append("path")
            .data([this.variability_c])
            .attr("class", "line")
            .attr("name", "second_pos")
            .style("fill", '#990000')
            .attr("d", this.reverseLine);

        this.lineContainer.append("path")
            .data([this.variability_d])
            .attr("class", "line")
            .attr("name", "third_pos")
            .style("fill", '#004400')
            .attr("d", this.reverseLine);

        this.lineContainer.append("path")
            .data([this.variability_a])
            .attr("class", "line")
            .attr("name", "outside_gene")
            .style("stroke", '#666666')
            .style("stroke-width", "0.2")
            .attr("d", this.reverseLine);

        this.textContainer.selectAll("text")
                                .data(d3.entries(this.competing_nucleotides))
                                .enter()
                                .append("text")
                                .attr("class", "SNV_text")
                                .attr("x", function (d) { return xS(0.5+parseInt(d.key)); })
                                .attr("y", function (d) { return 0; })
                                .attr("writing-mode", "tb")
                                .attr("font-size", "5px")
                                .attr("glyph-orientation-vertical", "0")
                                .attr("style", "cursor:pointer;")
                                .attr("paint-order", "stroke")
                                .attr("stroke", "#FFFFFF")
                                .attr("stroke-width", "1px")
                                .attr("fill", function (d){ return get_comp_nt_color(d.value['competing_nts']); })
                                .attr('data-content', function(d) {
                                    return '<span class="popover-close-button" onclick="$(this).closest(\'.popover\').popover(\'hide\');"></span> \
                                            <h3>Content</h3> \
                                            <table class="table table-striped" style="width: 100%; text-align: center; font-size: 12px;"> \
                                                <tr><td>Position in split</td><td>' + ((gene_mode) ? d.value['pos_in_split'] : d.value['pos']) +'</td></tr> \
                                                <tr><td>Position in contig</td><td>' + d.value['pos_in_contig'] +'</td></tr> \
                                                <tr><td>Reference</td><td>' + d.value['reference'] +'</td></tr> \
                                                <tr><td>Consensus</td><td>' + d.value['consensus'] +'</td></tr> \
                                                <tr><td>Departure from reference</td><td>' + d.value['departure_from_reference'].toFixed(4) +'</td></tr> \
                                                <tr><td>Departure from consensus</td><td>' + d.value['departure_from_consensus'].toFixed(4) +'</td></tr> \
                                                <tr><td>Competing nucleotides</td><td>' + d.value['competing_nts'] +'</td></tr> \
                                                <tr><td>Corresponding gene call</td><td>' + ((d.value['corresponding_gene_call'] == -1) ? 'No gene or in partial gene': d.value['corresponding_gene_call']) +'</td></tr> \
                                                <tr><td>Codon order in gene</td><td>' + ((d.value['codon_order_in_gene'] == -1) ? 'No gene or in noncoding gene': d.value['codon_order_in_gene']) +'</td></tr> \
                                                <tr><td>Base position in codon</td><td>' + ((d.value['base_pos_in_codon'] == 0) ? 'No gene or in noncoding gene': d.value['base_pos_in_codon']) +'</td></tr> \
                                                <tr><td>Coverage</td><td>' + d.value['coverage'] +'</td></tr> \
                                            </table> \
                                            <h3>Counts</h3> \
                                            <table class="table table-striped" style="width: 100%; text-align: center; font-size: 12px;"> \
                                                <tr><td>A</td><td>' + d.value['A'] +'</td></tr> \
                                                <tr><td>T</td><td>' + d.value['T'] +'</td></tr> \
                                                <tr><td>G</td><td>' + d.value['G'] +'</td></tr> \
                                                <tr><td>C</td><td>' + d.value['C'] +'</td></tr> \
                                                <tr><td>N</td><td>' + d.value['N'] +'</td></tr> \
                                            </table>';
                                })
                                .attr('data-toggle', 'popover')
                                .text(function (d) {
                                    return d.value['competing_nts'];
                                });
                                console.log("Drawing SNV markers (done: " + Math.round(Date.now()/1000) + ")");
    }

    if(state['show_indels']) {
      // save new dict where coverage vals < <VAL> are removed
      /*var f_indels = this.indels;
      f_indels.forEach((item, i) => { //foreachdoesnt work for dicts...
        Object.keys(item).forEach((j) => {
          if(item[j].coverage < 0) f_indels[i][j].coverage = 0;
        });

        //f_indels = Object.keys(item).reduce(function(f_indels,key){ if(item[key].coverage < 0) filtered[key] = 0; return f_indels; })
      });*/

      this.indel_coverage = [];
      this.indel_coverage.length = this.coverage.length;
      this.indel_coverage.fill(0);
      /* dictionary to hold data for positions with multiple indels */
      var mult_indels = {};
      d3.entries(this.indels).forEach(obj => {
        if(this.indel_coverage[obj.value['pos']] != 0) {
          this.indel_coverage[obj.value['pos']] += (obj.value['count']/obj.value['coverage']);

          var new_table = '<span class="popover-close-button" onclick="$(this).closest(\'.popover\').popover(\'hide\');"></span> \
                  <h3>' + ((obj.value['type'] == 'INS') ? 'Insertion' : 'Deletion') + '</h3> \
                  <table class="table table-striped" style="width: 100%; text-align: center; font-size: 12px;"> \
                      <tr><td>Position in split</td><td>' + obj.value['pos'] +'</td></tr> \
                      <tr><td>Position in contig</td><td>' + obj.value['pos_in_contig'] +'</td></tr> \
                      <tr><td>Reference</td><td>' + obj.value['reference'] +'</td></tr> \
                      <tr><td>Sequence</td><td>' + obj.value['sequence'] +'</td></tr> \
                      <tr><td>Type</td><td>' + ((obj.value['type'] == 'INS') ? 'Insertion' : 'Deletion') +'</td></tr> \
                      <tr><td>Length</td><td>' + obj.value['length'] +'</td></tr> \
                      <tr><td>Count</td><td>' + obj.value['count'] +'</td></tr> \
                      <tr><td>Corresponding gene call</td><td>' + ((obj.value['corresponding_gene_call'] == -1) ? 'No gene or in partial gene': obj.value['corresponding_gene_call']) +'</td></tr> \
                      <tr><td>Codon order in gene</td><td>' + ((obj.value['codon_order_in_gene'] == -1) ? 'No gene or in noncoding gene': obj.value['codon_order_in_gene']) +'</td></tr> \
                      <tr><td>Base position in codon</td><td>' + ((obj.value['base_pos_in_codon'] == 0) ? 'No gene or in noncoding gene': obj.value['base_pos_in_codon']) +'</td></tr> \
                      <tr><td>Coverage</td><td>' + obj.value['coverage'] +'</td></tr> \
                  </table>';

          if(mult_indels[obj.value['pos']]) {
            mult_indels[obj.value['pos']][0]++;
            mult_indels[obj.value['pos']][1] += new_table;
            mult_indels[obj.value['pos']][2] = mult_indels[obj.value['pos']][2] > obj.value['length'] ? mult_indels[obj.value['pos']][2] : obj.value['length'];
          } else {
            var initialKey = Object.keys(this.indels).find(key => this.indels[key]['pos'] === obj.value['pos']);
            var init_indel = this.indels[initialKey];
            init_table = '<span class="popover-close-button" onclick="$(this).closest(\'.popover\').popover(\'hide\');"></span> \
                    <h3>' + ((init_indel['type'] == 'INS') ? 'Insertion' : 'Deletion') + '</h3> \
                    <table class="table table-striped" style="width: 100%; text-align: center; font-size: 12px;"> \
                        <tr><td>Position in split</td><td>' + init_indel['pos'] +'</td></tr> \
                        <tr><td>Position in contig</td><td>' + init_indel['pos_in_contig'] +'</td></tr> \
                        <tr><td>Reference</td><td>' + init_indel['reference'] +'</td></tr> \
                        <tr><td>Sequence</td><td>' + init_indel['sequence'] +'</td></tr> \
                        <tr><td>Type</td><td>' + ((init_indel['type'] == 'INS') ? 'Insertion' : 'Deletion') +'</td></tr> \
                        <tr><td>Length</td><td>' + init_indel['length'] +'</td></tr> \
                        <tr><td>Count</td><td>' + init_indel['count'] +'</td></tr> \
                        <tr><td>Corresponding gene call</td><td>' + ((init_indel['corresponding_gene_call'] == -1) ? 'No gene or in partial gene': init_indel['corresponding_gene_call']) +'</td></tr> \
                        <tr><td>Codon order in gene</td><td>' + ((init_indel['codon_order_in_gene'] == -1) ? 'No gene or in noncoding gene': init_indel['codon_order_in_gene']) +'</td></tr> \
                        <tr><td>Base position in codon</td><td>' + ((init_indel['base_pos_in_codon'] == 0) ? 'No gene or in noncoding gene': init_indel['base_pos_in_codon']) +'</td></tr> \
                        <tr><td>Coverage</td><td>' + init_indel['coverage'] +'</td></tr> \
                    </table>';

            mult_indels[obj.value['pos']] = [];
            mult_indels[obj.value['pos']][0] = 2;
            mult_indels[obj.value['pos']][1] = init_table + new_table;
            mult_indels[obj.value['pos']][2] = init_indel['length'] > obj.value['length'] ? init_indel['length'] : obj.value['length'];
          }

        } else {
          this.indel_coverage[obj.value['pos']] = (obj.value['count']/obj.value['coverage'] < 1) ? (obj.value['count']/obj.value['coverage']) : 1;
        }
      });

      for (pos in mult_indels) {
        this.indel_coverage[pos] = this.indel_coverage[pos] / mult_indels[pos][0];
      }

      this.lineContainer.append("path")
          .data([this.indel_coverage])
          .attr("class", "reverseLine")
          .attr("name", "indel_1")
          .style("fill", '#800080')
          .attr("d", this.line);

      // add text to text container based on type, and data-content based on other variables
      this.textContainerIndels.selectAll("text")
                              .data(d3.entries(this.indels))
                              .enter()
                              .append("text")
                              //.filter(function(d){ return d.value['coverage'] >= state['min-indel-coverage']})
                              .attr("class", "indels_text")
                              .attr("x", function (d) { return xS(0.5+d.value['pos']); })
                              .attr("y", function (d) { return ySL(0); })
                              .attr("font-size", "14px")
                              .attr("style", "cursor:pointer;")
                              .attr("fill", function(d) { return (((d.value['pos'] in mult_indels ? mult_indels[d.value['pos']][2] : d.value['length']) > state['large-indel']) ? 'red' : '#CCCC00'); })
                              .attr('data-content', function(d) {
                                  if(d.value['pos'] in mult_indels) {
                                    return mult_indels[d.value['pos']][1];
                                  }

                                  return '<span class="popover-close-button" onclick="$(this).closest(\'.popover\').popover(\'hide\');"></span> \
                                          <h3>' + ((d.value['type'] == 'INS') ? 'Insertion' : 'Deletion') + '</h3> \
                                          <table class="table table-striped" style="width: 100%; text-align: center; font-size: 12px;"> \
                                              <tr><td>Position in split</td><td>' + d.value['pos'] +'</td></tr> \
                                              <tr><td>Position in contig</td><td>' + d.value['pos_in_contig'] +'</td></tr> \
                                              <tr><td>Reference</td><td>' + d.value['reference'] +'</td></tr> \
                                              <tr><td>Sequence</td><td>' + d.value['sequence'] +'</td></tr> \
                                              <tr><td>Type</td><td>' + ((d.value['type'] == 'INS') ? 'Insertion' : 'Deletion') +'</td></tr> \
                                              <tr><td>Length</td><td>' + d.value['length'] +'</td></tr> \
                                              <tr><td>Count</td><td>' + d.value['count'] +'</td></tr> \
                                              <tr><td>Corresponding gene call</td><td>' + ((d.value['corresponding_gene_call'] == -1) ? 'No gene or in partial gene': d.value['corresponding_gene_call']) +'</td></tr> \
                                              <tr><td>Codon order in gene</td><td>' + ((d.value['codon_order_in_gene'] == -1) ? 'No gene or in noncoding gene': d.value['codon_order_in_gene']) +'</td></tr> \
                                              <tr><td>Base position in codon</td><td>' + ((d.value['base_pos_in_codon'] == 0) ? 'No gene or in noncoding gene': d.value['base_pos_in_codon']) +'</td></tr> \
                                              <tr><td>Coverage</td><td>' + d.value['coverage'] +'</td></tr> \
                                          </table>';
                              })
                              .attr('data-toggle', 'popover')
                              .text(function (d) {
                                  if(d.value['pos'] in mult_indels) return 'x';
                                  return d.value['type'] == 'INS' ? '+' : '-';
                              });
                              console.log("Drawing indel markers (done: " + Math.round(Date.now()/1000) + ")");
    }



    this.xAxisTop = d3.svg.axis().scale(this.xScale).orient("top");

    if(this.id == 0){
        this.chartContainer.append("g")
                    .attr("class", "x axis top noselect")
                    .attr("transform", "translate(0,0)")
                    .call(this.xAxisTop);
    }


    this.yAxis = d3.svg.axis().scale(this.yScale).orient("left").ticks(5);
    this.yAxisLine = d3.svg.axis().scale(this.yScaleLine).orient("right").ticks(5);

    this.chartContainer.append("g")
                   .attr("class", "y axis noselect")
                   .attr("transform", "translate(-10,0)")
                   .call(this.yAxis);

    if(state['show_snvs'] || state['show_indels']) {
      this.lineContainer.append("g")
                     .attr("class", "y axis noselect")
                     .attr("transform", "translate(" + (this.width + 15) + ",0)")
                     .call(this.yAxisLine);
    }

    this.chartContainer.append("text")
                   .attr("class","sample-title")
                   .attr("transform", "translate(0,20)")
                   .text(this.name);

}

Chart.prototype.showOnly = function(b){
    this.xScale.domain(b); var xS = this.xScale;
    this.chartContainer.selectAll("path").data([this.coverage]).attr("d", this.area);
    this.gcContainer.selectAll("path").data([this.gc_content]).attr("d", this.gc_line);
    this.lineContainer.select("[name=outside_gene]").data([this.variability_a]).attr("d", this.reverseLine);
    this.lineContainer.select("[name=first_pos]").data([this.variability_b]).attr("d", this.reverseLine);
    this.lineContainer.select("[name=second_pos]").data([this.variability_c]).attr("d", this.reverseLine);
    this.lineContainer.select("[name=third_pos]").data([this.variability_d]).attr("d", this.reverseLine);
    if(indels_enabled) this.lineContainer.select("[name=indel_1]").data([this.indel_coverage]).attr("d", this.line);
    this.textContainer.selectAll(".SNV_text").data(d3.entries(this.competing_nucleotides)).attr("x", function (d) { return xS(0.5+parseInt(d.key)); });
    this.textContainerIndels.selectAll(".indels_text").data(d3.entries(this.indels)).attr("x", function (d) { return xS(0.5+d.value['pos']); });

    let numNucl = $('#brush_end').val()-$('#brush_start').val();
    let mk_font_size = 2000/numNucl;
    if(mk_font_size < 5) mk_font_size = 5;
    if(mk_font_size > 10) mk_font_size = 10;
    this.textContainer.selectAll(".SNV_text").data(d3.entries(this.competing_nucleotides)).attr("font-size", mk_font_size+"px");
    this.textContainerIndels.selectAll(".indels_text").data(d3.entries(this.indels)).attr("font-size", 2*mk_font_size+"px");

    this.chartContainer.select(".x.axis.top").call(this.xAxisTop);
}
