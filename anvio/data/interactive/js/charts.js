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

var layers;
var coverage;
var variability;
var maxVariability = 0;
var geneParser;
var contextSvg;
var state;
var layers_ordered;
var visible_layers;
var contig_id;
var page_header;
var highlight_gene;
var gene_mode;
var show_snvs;
var sequence;
var charts;
var brush;
var inspect_mode;
var show_nucleotides = true;


function loadAll() {
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
    show_snvs = getParameterByName('show_snvs') == 'true';

    if (typeof localStorage.state === 'undefined')
    {
        state = {}
    }
    else
    {
        state = JSON.parse(localStorage.state);
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

                for (var i=0; i<coverage.length; i++) {
                    variability[i] = [];
                    for (var l=0; l<4; l++) {
                        variability[i][l] = [];
                        for (var h=0; h<coverage[i].length; h++) {
                            if (contig_data.variability[i][l].hasOwnProperty(h)) {
                                variability[i][l].push(contig_data.variability[i][l][h]);
                                if (contig_data.variability[i][l][h] > maxVariability) {
                                    maxVariability = contig_data.variability[i][l][h];
                                }
                            } else {
                                variability[i][l].push(0);
                            }
                        }
                    }
                }

                competing_nucleotides = contig_data.competing_nucleotides;
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
                    next_str = '<a onclick="localStorage.state = JSON.stringify(state);" href="' + generate_inspect_link({'type': inspect_mode, 'item_name': next_contig_name, 'show_snvs': show_snvs}) +'" '+target_str+'> | next &gt;&gt;&gt;</a>';

                if(previous_contig_name)
                    prev_str = '<a onclick="localStorage.state = JSON.stringify(state);" href="' + generate_inspect_link({'type': inspect_mode, 'item_name': previous_contig_name, 'show_snvs': show_snvs}) + '" '+target_str+'>&lt;&lt;&lt; prev | </a>';

                $('#header').append("<strong>" + page_header + "</strong> detailed");
                $('#header').append("<p><small><small>" + prev_str + position + next_str + "</small></small></p>");
                $('#header').append("<p style='margin-top: -30px; margin-bottom: 15px;'><small><small><a href='#' onclick='showSearchItemsDialog();'>Select or Search Item</a></small></small></p>");

                $('.main').prepend(`<div style="float: right; text-align: right; padding-right: 60px; padding-bottom: 20px; display: inline-block;" class="form-inline"> \
                                        <b>Range:</b>
                                            <input class="form-control input-sm" id="brush_start" type="text" value="0" size="5">
                                        <b>:</b>
                                            <input class="form-control input-sm" id="brush_end" type="text" value="${sequence.length}" size="5">\
                                    </div>`);

                $('.main').prepend(`<div style="text-align: right; padding-left: 40px; padding-bottom: 20px; display: inline-block;"> \
                                        <button type="button" class="btn btn-primary btn-xs" onclick="show_sequence_modal('Sequence', page_header + '\\n' + sequence);">Get sequence</button> \
                                        <button type="button" class="btn btn-primary btn-xs disabled btn-selection-sequence"  onclick="show_selected_sequence();" disabled>Get sequence of selected area</button> \
                                    </div>`);

                $('.main').prepend(`<div style="text-align: right; padding-left: 40px; padding-bottom: 20px; display: inline-block;"> \
                                        <button type="button" class="btn btn-primary btn-xs" onclick="showOverlayGCContentDialog();" class="btn btn-outline-primary">Overlay GC Content</button> \
                                        <button type="button" class="btn btn-primary btn-xs" onclick="resetOverlayGCContent();" class="btn btn-outline-primary">Reset overlay</button> \
                                    </div>`);

                $('.main').prepend('<div style="text-align: left; padding-left: 40px; padding-bottom: 20px; display: inline-block;"> \
                                        <button type="button" class="btn btn-primary btn-xs" onclick="showSetMaxValuesDialog()" class="btn btn-outline-primary">Set maximum values</button> \
                                        <button type="button" class="btn btn-primary btn-xs" onclick="resetMaxValues()" class="btn btn-outline-primary">Reset maximum values</button> \
                                    </div>');

                createCharts(state);
                $('.loading-screen').hide();

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
            }
        });

}

function toggle_nucleotide_display() {
  show_nucleotides = !show_nucleotides;
  if(show_nucleotides) {
    display_nucleotides();
    $("div.nucl-activated").fadeIn(300).delay(1500).fadeOut(400);
  } else {
    contextSvg.select("#DNA_sequence").remove();
    contextSvg.select("#AA_sequence").remove();
    contextSvg.select("#solids").remove();
    contextSvg.attr("height", 150);
    $("div.nucl-deactivated").fadeIn(300).delay(1500).fadeOut(400);
  }
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

  if(end - start > 300 || end - start < 30) {
    contextSvg.attr("height", 150);
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

  var nucl_sequence = contextSvg.append("text")
                                .text(sequence.substring(start, end))
                                .attr("id", "DNA_sequence")
                                .attr("fill", "url(#solids)");
  /* width of monospaced character per font size */
  var nucl_text_font = width/((end-start)*.6002738402061856);
                   nucl_sequence.attr("font-size", nucl_text_font);
  var dna_seq_height = contextSvg.select("#DNA_sequence")[0][0].getBBox().height;
  var nucl_text_y = 140 + .75*dna_seq_height;
                   nucl_sequence.attr("y", nucl_text_y)
                                .attr("font-family", "monospace")
                                .attr("transform", "translate(" + (margin.left) + ", 0)");

  var show_AAs = false;
  geneParser["data"].forEach(function(gene){
    if(gene.start_in_split < end-2 && gene.stop_in_split > start+2) show_AAs = true;
  });

  if(show_AAs) {
    var aa_sequence = contextSvg.append("g")
                                .attr("id", "AA_sequence")
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
                     .attr("height", dna_seq_height)
                     .attr("width", 3*textWidth)
                     .attr("x", rect_x)
                     .attr("y", nucl_text_y + offset_y)
                     .attr("fill", i % 2 == 0 ? "rgb(144,137,250)" : "rgb(81,68,211)");
          rect_x += (3*textWidth);
        }

        aa_sequence.append("text")
                  .text(aa_string)
                  .attr('font-size', nucl_text_font)
                  .attr("font-family", "monospace")
                  .attr("fill", "white")
                  .attr("x", buffer*textWidth)
                  .attr("y", nucl_text_y + .67*dna_seq_height + offset_y);
      }
      prev_gene_stop = gene.stop_in_split;
    }
  }

  contextSvg.attr("height", 150 + dna_seq_height +
    (show_AAs? dna_seq_height : 0) +
    (overlapping_genes? 5+dna_seq_height : 0));
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
}


function resetOverlayGCContent() {
    delete sessionStorage.gc_overlay_settings;
    createCharts(state);
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
}


function resetMaxValues() {
    delete sessionStorage.max_coverage;
    createCharts(state);
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

                    let link = '<a onclick="localStorage.state = JSON.stringify(state);" href="' + generate_inspect_link({'type': inspect_mode, 'item_name': item_name, 'show_snvs': show_snvs}) +'" '+target_str+'>' + item_name_pretty + '</a>';
                    results_html += link + '<br />';

                }

                results_html += '<br /><br /><center>';

                if (results.length > 0) {
                    if (page + 1 < total_page) {
                        results_html += `<a href="#" onclick="search_items('${search_query}', ${page+1});">&lt;&lt;&lt; prev</a> | `;
                    }

                    results_html += " page " + (page + 1) + " of " + total_page;

                    if (page > 0) {
                        results_html += ` | <a href="#" onclick="search_items('${search_query}', ${page-1});"> next &gt;&gt;&gt;</a>`;
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

function createCharts(state){
    /* Adapted from Tyler Craft's Multiple area charts with D3.js article:
    http://tympanus.net/codrops/2012/08/29/multiple-area-charts-with-d3-js/  */
    $('#chart-container, #context-container').empty();

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

    geneParser = new GeneParser(genes);

    var margin = {top: 20, right: 50, bottom: 150, left: 50};
    var width = VIEWER_WIDTH * .80;
    var chartHeight = 200;
    var height = ((chartHeight + 10) * visible_layers);
    var contextHeight = 50;
    var contextWidth = width;

    var svg = d3.select("#chart-container").append("svg")
            .attr("width", width + margin.left + margin.right)
            .attr("height", (height + margin.top + margin.bottom));

    $('#chart-container').css("width", (width + 150) + "px");
    $('#chart-container').css("height", height + "px");


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
                        gc_content: gc_content_array,
                        'gc_content_window_size': gc_content_window_size,
                        'gc_content_step_size': gc_content_step_size,
                        'gc_overlay_color': gc_overlay_color,
                        id: j++,
                        width: width,
                        height: chartHeight,
                        maxVariability: maxVariability,
                        svg: svg,
                        margin: margin,
                        showBottomAxis: (j == visible_layers - 1),
                        color: state['layers'][layers[layer_index]]['color']
                }));

    }


    contextSvg = d3.select("#context-container").append("svg")
            .attr("width", width + margin.left + margin.right)
            .attr("height", 150);

    var defs = contextSvg.append('svg:defs');

    contextSvg.append("rect")
       .attr("width", width)
       .attr("height", "60px")
       .attr("fill", "black")
       .attr("fill-opacity", "0.2")
       .attr('transform', 'translate(50, 10)');

    // Define arrow markers
    ['green', 'gray', 'firebrick', '#226ab2'].forEach(function(color){
      defs.append('svg:marker')
          .attr('id', 'arrow_' + color )
          .attr('markerHeight', 2)
          .attr('markerWidth', 2)
          .attr('orient', 'auto')
          .attr('refX', 0)
          .attr('refY', 0)
          .attr('viewBox', '-5 -5 10 10')
          .append('svg:path')
            .attr('d', 'M 0,0 m -5,-5 L 5,0 L -5,5 Z')
            .attr('fill', color);
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
                .attr("class","context")
                .attr("transform", "translate(" + (margin.left) + ", 80)");

    context.append("g")
                .attr("class", "x axis top")
                .attr("transform", "translate(0,0)")
                .call(contextAxis)

    context.append("g")
                .attr("class", "x brush")
                .call(brush)
                .selectAll("rect")
                .attr("y", 0)
                .attr("height", contextHeight);

    display_nucleotides();

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

        for(var i = 0; i < layersCount; i++){
            charts[i].showOnly(b);
        }
        drawArrows(b[0], b[1]);
    }

    drawArrows(0, charts[0].xScale.domain()[1]);
}


function Chart(options){
    this.coverage = options.coverage;
    this.max_coverage = options.max_coverage;
    this.variability_a = options.variability_a;
    this.variability_b = options.variability_b;
    this.variability_c = options.variability_c;
    this.variability_d = options.variability_d;
    this.competing_nucleotides = options.competing_nucleotides;
    this.gc_content = options.gc_content;
    this.gc_content_window_size = options.gc_content_window_size;
    this.gc_content_step_size = options.gc_content_step_size;
    this.gc_overlay_color = options.gc_overlay_color;
    this.width = options.width;
    this.height = options.height;
    this.maxVariability = options.maxVariability;
    this.svg = options.svg;
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

    this.yScale = d3.scale.linear()
                            .range([this.height,0])
                            .domain([0,this.maxCoverageForyScale]);

    this.yScaleLine = d3.scale.linear()
                            .range([this.height, 0])
                            .domain([0, this.maxVariability]);

    this.yScaleGC = d3.scale.linear()
                            .range([this.yScale(this.minCoverage),
                                   (this.maxCoverage < this.maxCoverageForyScale) ? this.yScale(this.maxCoverage) : 0])
                            .domain([this.minGCContent, this.maxGCContent]);

    var xS = this.xScale;
    var yS = this.yScale;
    var ySL = this.yScaleLine;
    var yGC = this.yScaleGC;

    this.area = d3.svg.area()
                            .x(function(d, i) { return xS(i); })
                            .y0(this.height)
                            .y1(function(d) { return (yS(d) < 0) ? 0 : yS(d); });

    this.line = d3.svg.line()
                            .x(function(d, i) { return xS(i); })
                            .y(function(d, i) { if(i == 0) return ySL(0); if(i == num_data_points - 1) return ySL(0); return ySL(d); })
                            .interpolate('step-before');


    // .x() needs to stay as a arrow function, it has reference to scope.
    this.gc_line = d3.svg.line()
                            .x((d, i) => { return xS((this.gc_content_window_size / 2) + (i * this.gc_content_step_size)); })
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

    this.textContainer = this.svg.append("g")
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

    if (show_snvs) {
        this.lineContainer.append("path")
            .data([this.variability_b])
            .attr("class", "line")
            .attr("name", "first_pos")
            .style("fill", '#990000')
            .attr("d", this.line);

        this.lineContainer.append("path")
            .data([this.variability_c])
            .attr("class", "line")
            .attr("name", "second_pos")
            .style("fill", '#990000')
            .attr("d", this.line);

        this.lineContainer.append("path")
            .data([this.variability_d])
            .attr("class", "line")
            .attr("name", "third_pos")
            .style("fill", '#004400')
            .attr("d", this.line);

        this.lineContainer.append("path")
            .data([this.variability_a])
            .attr("class", "line")
            .attr("name", "outside_gene")
            .style("stroke", '#666666')
            .style("stroke-width", "0.2")
            .attr("d", this.line);

        this.textContainer.selectAll("text")
                                .data(d3.entries(this.competing_nucleotides))
                                .enter()
                                .append("text")
                                .attr("x", function (d) { return xS(d.key); })
                                .attr("y", function (d) { return 0; })
                                .attr("writing-mode", "tb")
                                .attr("font-size", "7px")
                                .attr("glyph-orientation-vertical", "0")
                                .attr("style", "cursor:pointer;")
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
    }



    this.xAxisTop = d3.svg.axis().scale(this.xScale).orient("top");

    if(this.id == 0){
        this.chartContainer.append("g")
                    .attr("class", "x axis top")
                    .attr("transform", "translate(0,0)")
                    .call(this.xAxisTop);
    }


    this.yAxis = d3.svg.axis().scale(this.yScale).orient("left").ticks(5);
    this.yAxisLine = d3.svg.axis().scale(this.yScaleLine).orient("right").ticks(5);

    this.chartContainer.append("g")
                   .attr("class", "y axis")
                   .attr("transform", "translate(-10,0)")
                   .call(this.yAxis);

    this.lineContainer.append("g")
                   .attr("class", "y axis")
                   .attr("transform", "translate(" + (this.width + 15) + ",0)")
                   .call(this.yAxisLine);

    this.chartContainer.append("text")
                   .attr("class","country-title")
                   .attr("transform", "translate(0,20)")
                   .text(this.name);

}

Chart.prototype.showOnly = function(b){
    this.xScale.domain(b); var xS = this.xScale;
    this.chartContainer.selectAll("path").data([this.coverage]).attr("d", this.area);
    this.gcContainer.selectAll("path").data([this.gc_content]).attr("d", this.gc_line);
    this.lineContainer.select("[name=outside_gene]").data([this.variability_a]).attr("d", this.line);
    this.lineContainer.select("[name=first_pos]").data([this.variability_b]).attr("d", this.line);
    this.lineContainer.select("[name=second_pos]").data([this.variability_c]).attr("d", this.line);
    this.lineContainer.select("[name=third_pos]").data([this.variability_d]).attr("d", this.line);
    this.textContainer.selectAll("text").data(d3.entries(this.competing_nucleotides)).attr("x", function (d) { return xS(d.key); });
    this.chartContainer.select(".x.axis.top").call(this.xAxisTop);
}
