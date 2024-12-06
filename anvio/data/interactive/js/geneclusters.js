/**
 * Javascript library to visualize anvi'o gene clusterss
 *
 *  Authors: A. Murat Eren <a.murat.eren@gmail.com>
 *           Özcan Esen
 *           Mahmoud Yousef
 *
 * Copyright 2016-2021, The anvi'o project (http://anvio.org)
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

var genomes;
var gene_caller_ids;
var gene_caller_ids_in_genomes;
var aa_sequences_in_gene_cluster;
var previous_gene_cluster_name;
var next_gene_cluster_name;
var index;
var total;

var state;
var gene_cluster_data;
var psgc_data = null;
var mode = null;

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

    gene_cluster_name = getUrlVars()["id"];
    document.title = gene_cluster_name + " detailed";
    state = JSON.parse(localStorage.state);

    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/geneclusters/' + state['order-by'] + '/' + gene_cluster_name,
        success: function(_gene_cluster_data) {
            gene_cluster_data = _gene_cluster_data;
            genomes = gene_cluster_data.genomes;
            gene_caller_ids = gene_cluster_data.gene_caller_ids;
            gene_caller_ids_in_genomes = gene_cluster_data.gene_caller_ids_in_genomes;
            aa_sequences_in_gene_cluster = gene_cluster_data.aa_sequences_in_gene_cluster;
            previous_gene_cluster_name = gene_cluster_data.previous_gene_cluster_name;
            next_gene_cluster_name = gene_cluster_data.next_gene_cluster_name;
            index = gene_cluster_data.index;
            total = gene_cluster_data.total;

            if(genomes.length == 0){
                console.log('Warning: no genomes returned')
            }

            next_str = " | next &gt;&gt;&gt;";
            prev_str = "&lt;&lt;&lt; prev | ";
            position = index + " of " + total;

            // anvi-server uses iframes for prettier urls, links need to be open _top
            var target_str = '';

            if (self != top) {
                target_str = 'target="_top"';
            }

            if(next_gene_cluster_name)
                next_str = '<a onclick="localStorage.state = JSON.stringify(state);" href="' + generate_inspect_link({'type': 'geneclusters', 'item_name': next_gene_cluster_name, 'show_snvs': false}) +'" '+target_str+'> | next &gt;&gt;&gt;</a>';

            if(previous_gene_cluster_name)
                prev_str = '<a onclick="localStorage.state = JSON.stringify(state);" href="' + generate_inspect_link({'type': 'geneclusters', 'item_name': previous_gene_cluster_name, 'show_snvs': false}) +'" '+target_str+'>&lt;&lt;&lt; prev | </a>';

            const headerContent = `<strong>${gene_cluster_name}</strong> with ${gene_caller_ids.length} genes detailed <br /><small><small>${prev_str} ${position} ${next_str}</small></small>`;
            document.getElementById("header").innerHTML = headerContent;

            if (typeof localStorage.state === 'undefined')
            {
                alert("Something went wrong, couldn't access to localStorage");
            }
            else
            {
                initializeCheckBoxes();
                createDisplay(true);
                $('.loading-screen').hide();
            }
        }
    });
}
async function loadGCAdditionalData(gc_id, gc_key, gc_key_short) {
    try {
        const response = await $.ajax({
            type: 'POST',
            cache: false,
            url: '/data/get_additional_gc_data/' + gc_id + '/' + gc_key
        });
        if (response['status'] === 0) {
            var newThHeader = $('<th>').text(gc_key_short);
            var newThData = $('<th>').text(gc_key === 'num_genes_in_gene_cluster' || gc_key === 'num_genomes_gene_cluster_has_hits' ? Math.round(response.gene_cluster_data): (response.gene_cluster_data).toFixed(2));

            var gc_title_list = {
                combined_homogeneity_index: "Combined Homogeneity Index",
                functional_homogeneity_index : "Functional Homogeneity Index",
                geometric_homogeneity_index : "Geometric Homogeneity Index",
                num_genes_in_gene_cluster : "Number of [g]enes in Gene Cluster",
                num_genomes_gene_cluster_has_hits : "Number of [G]enomes in Gene Cluster",
                max_num_paralogs : "Maximum Number of Paralogs",
                AAI_avg : "Amino Acid Identity Average",
                AAI_max : "Amino Acid Identity Maximum",
                AAI_min : "Amino Acid Identity Minimum",
                SCG : "Single-copy Core Gene"
            };

            $('#gc-acc-table-header').parent().append(newThHeader);
            $('#gc-acc-table-data').parent().append(newThData);

            newThHeader.prop('title', gc_title_list[gc_key]);
            $('#gc-acc-table').show();
        } else {
            console.log('Error:', response.message);
        }
    } catch (error) {
        console.error('AJAX Error:', error);
    }
}

async function createDisplay(display_table){
    var sequence_wrap_val = parseInt($('#wrap_length').val());
    var sequence_font_size_val = parseInt($('#font_size').val());

    var sequence_wrap = (isNumber(sequence_wrap_val) && sequence_wrap_val > 0) ? sequence_wrap_val : 1042;
    var sequence_font_size = (isNumber(sequence_font_size_val) && sequence_font_size_val > 0) ? sequence_font_size_val : 12;

    var svg = document.getElementById('svg');
    var table = document.getElementById('gc-acc-main');

    try {
        psgc_data =await loadPSGCData(gene_cluster_data.gene_cluster_name);
        
        if (psgc_data) {
            mode = 'structure';
            gc_type_in_psgc = await loadGCTypeInPSGCData(gene_cluster_data.gene_cluster_name);
        }
    } catch (error) {
        console.error('Error loading PSGC data:', error);
    }

    // clear content
    while (svg.firstChild) {
        svg.removeChild(svg.firstChild);
    }

    var y_cord = 0;
    var offset = 0;

    var acid_sequences = [];
    var order = {};
    var count = 0;
    var layer_id_list = {
        combined_homogeneity_index: "CHI",
        functional_homogeneity_index : "FHI",
        geometric_homogeneity_index : "GHI",
        num_genes_in_gene_cluster : "NgGC",
        num_genomes_gene_cluster_has_hits : "NGGC",
        max_num_paralogs : "MNP",
        AAI_avg : "AAI_avg",
        AAI_max : "AAI_max",
        AAI_min : "AAI_min",
        SCG : "SCG"
    };

    for (var layer_id = 0; layer_id < state['layer-order'].length; layer_id++)
    {
        var layer = state['layer-order'][layer_id];

        if (layer_id_list.hasOwnProperty(layer) && display_table) {
            await loadGCAdditionalData(gene_cluster_data.gene_cluster_name, layer,  layer_id_list[layer]);
        }

        if (gene_cluster_data.genomes.indexOf(layer) === -1)
            continue;

        gene_cluster_data.gene_caller_ids_in_genomes[layer].forEach(function(caller_id) {
            acid_sequences.push(gene_cluster_data.aa_sequences_in_gene_cluster[layer][caller_id]);
            order[caller_id] = count;
            count = count + 1;
        });
    }

    var max_length = maxLength(acid_sequences);
    var all_positions = [];

    for (var i=0; i < max_length; i++) {
        var new_item = [];
        for (var j=0; j < acid_sequences.length; j++) {
            new_item.push(acid_sequences[j][i]);
        }
        all_positions.push(new_item);
    }
    var coded_positions = determineColor(all_positions);

    // Ensure this loop is only called once
    if (!window.hasRunLayerLoop) {
        var fragment = document.createDocumentFragment();
        for (var layer_id = 0; layer_id < state['layer-order'].length; layer_id++) {
            var layer = state['layer-order'][layer_id];

            if (state['layers'][layer]['height'] == 0)
                continue;

            if (gene_cluster_data.genomes.indexOf(layer) === -1)
                continue;

            var rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
            rect.setAttribute('x', 0);
            rect.setAttribute('y', y_cord);
            rect.setAttribute('fill', state['layers'][layer]['color']);
            rect.setAttribute('opacity', 0.2);
            rect.setAttribute('height', (Math.max(gene_cluster_data.gene_caller_ids_in_genomes[layer].length,1) * sequence_font_size * 1.5) + 10);
            rect.setAttribute('width', 400);
            rect.setAttribute('class', 'sequenceBackground');
            rect.setAttribute('rx', 10);
            fragment.appendChild(rect);

            var text = document.createElementNS('http://www.w3.org/2000/svg', 'text');
            text.setAttribute('x', 0);
            text.setAttribute('y', parseFloat(rect.getAttribute('y')) + parseFloat(rect.getAttribute('height')) / 2);
            text.setAttribute('style', 'alignment-baseline:central');
            text.setAttribute('class', 'genomeTitle');
            text.appendChild(document.createTextNode(layer));
            fragment.appendChild(text);

            sub_y_cord = y_cord + 5;
            gene_cluster_data.gene_caller_ids_in_genomes[layer].forEach(function (caller_id) {
                sequence = gene_cluster_data.aa_sequences_in_gene_cluster[layer][caller_id];
                var text = document.createElementNS('http://www.w3.org/2000/svg', 'text');
                text.setAttribute('x', 0);
                text.setAttribute('y', sub_y_cord);
                text.setAttribute('font-size', sequence_font_size);
                text.setAttribute('style', 'alignment-baseline:text-before-edge; cursor: pointer;');
                text.setAttribute('class', 'callerTitle');
                text.setAttribute('gene-callers-id', caller_id);
                text.setAttribute('genome-name', layer);
                text.setAttribute('data-toggle', 'popover');
                text.onclick = function(event) {
                    var obj = event.target;
                    var caller_id = obj.getAttribute('gene-callers-id');
                    var layer = obj.getAttribute('genome-name');
                    if (!obj.getAttribute('data-content')) {
                        obj.setAttribute('data-content', get_gene_functions_table_html_for_pan(caller_id, layer) + '');
                    }
                };

                text.appendChild(document.createTextNode(caller_id));
                fragment.appendChild(text);

                if (mode === 'structure') {
                    let gc_id = '';

                    if (psgc_data) {
                        for (var psgc_id in psgc_data) {
                            var matchingGene = psgc_data[psgc_id].find(gene => 
                                gene.gene_callers_id === caller_id
                            );

                            gc_id = matchingGene.gene_cluster_id;
                            break;
                        }
                    }

                    // Get the type indicator [A] || [C] || [S] if type data is available
                    let type_indicator = '';
                    let type_class = '';
                    if (gc_type_in_psgc) {
                        for (let psgc_id in gc_type_in_psgc) {
                            const type = gc_type_in_psgc[psgc_id][gc_id];

                            if (type === 'core') {
                                type_indicator = 'C';
                                type_class = 'type-core';
                            }
                            else if (type === 'accessory') {
                                type_indicator = 'A';
                                type_class = 'type-accessory';
                            }
                            else if (type === 'singleton') {
                                type_indicator = 'S';
                                type_class = 'type-singleton';
                            }

                            break;
                        }
                    }

                    var text = document.createElementNS('http://www.w3.org/2000/svg', 'text');
                    text.setAttribute('x', 0);
                    text.setAttribute('y', sub_y_cord);
                    text.setAttribute('font-size', sequence_font_size);
                    text.setAttribute('style', 'alignment-baseline:text-before-edge; cursor: pointer;');
                    text.setAttribute('class', 'callerTitle gcName');
                    text.setAttribute('gc-caller-id', gc_id);
                    text.setAttribute('genome-name', layer);
                    text.setAttribute('data-toggle', 'popover');
                    text.onclick = function(event) {
                        var obj = event.target;
                        if (!obj.getAttribute('data-content')) {
                            obj.setAttribute('data-content', get_gene_functions_table_html_for_structure(psgc_data, gc_id) + '');
                        }
                    };

                    text.appendChild(document.createTextNode(gc_id));

                    if (type_indicator) {
                        var typeSpan = document.createElementNS('http://www.w3.org/2000/svg', 'tspan');
                        typeSpan.setAttribute('class', type_class);
                        typeSpan.setAttribute('style', 'alignment-baseline:text-before-edge; cursor: text;');
                        typeSpan.setAttribute('dy', '0');
                        typeSpan.appendChild(document.createTextNode(' [' + type_indicator + ']'));
                        text.appendChild(typeSpan);
                    }
                    fragment.appendChild(text);
                }

                var text = document.createElementNS('http://www.w3.org/2000/svg', 'text');
                text.setAttribute('x', 0);
                text.setAttribute('y', sub_y_cord);
                text.setAttribute('font-size', sequence_font_size);
                text.setAttribute('font-family',"monospace");
                text.setAttribute('font-weight', '100');
                text.setAttribute('style', 'alignment-baseline:text-before-edge');
                text.setAttribute('class', 'sequence');

                _sequence = sequence.substr(offset, sequence_wrap);
                var fragmentForTspans = document.createDocumentFragment();
                for (var _letter_index=0; _letter_index < _sequence.length; _letter_index++) {
                    var index = _letter_index + offset;
                    var num = order[caller_id];
                    var acid = _sequence[_letter_index];
                    var dict = coded_positions[index][num];

                    var tspan = document.createElementNS('http://www.w3.org/2000/svg', 'tspan');
                    tspan.setAttribute('fill', dict[acid]);
                    tspan.style.fontWeight = 'bold';
                    tspan.setAttribute('style', 'alignment-baseline:text-before-edge');
                    tspan.appendChild(document.createTextNode(acid));
                    fragmentForTspans.appendChild(tspan);
                }
                text.appendChild(fragmentForTspans);
                fragment.appendChild(text);

                sub_y_cord = sub_y_cord + sequence_font_size * 1.5;

            });

            y_cord = y_cord + parseFloat(rect.getAttribute('height')) + 5;

        }
        svg.appendChild(fragment);
        offset += sequence_wrap;
        y_cord = y_cord + 15;
    }

    calculateLayout();

    $('[data-toggle="popover"]').popover({"html": true, sanitize: false,"trigger": "click", "container": "body", "viewport": "body", "placement": "top"});

    $('[data-toggle="popover"]').on('shown.bs.popover', function (e) {
      var popover = $(e.target).data("bs.popover").tip;

      if ($(popover).css('top').charAt(0) === '-') {
        $(popover).css('top', '0px');
      }

      if ($(popover).css('left').charAt(0) === '-') {
        $(popover).css('left', '0px');
      }
    });
}

function calculateLayout() {
    var HORIZONTAL_PADDING = 10;
    var svg = document.getElementById('svg');

    var max_genome_title_width = 0;
    $('.genomeTitle').each(function(i, d) {
        dwidth = d.getBBox().width;
        if (dwidth > max_genome_title_width) {
            max_genome_title_width = dwidth;
        };
    });

    $('.sequenceBackground').attr('x', max_genome_title_width + HORIZONTAL_PADDING);
    $('.callerTitle').attr('x', max_genome_title_width + (HORIZONTAL_PADDING * 2));

    var max_caller_width = 0;
    $('.callerTitle:not(.gcName)').each(function(i, d) {
        dwidth = d.getBBox().width;
        if (dwidth > max_caller_width) {
            max_caller_width = dwidth;
        };
    });

    var max_gc_name_width = 0;
    var sequence_start_x;

    if (mode === 'structure') {
        $('.gcName').attr('x', max_genome_title_width + max_caller_width + (HORIZONTAL_PADDING * 3));
        $('.gcName').each(function(i, d) {
            dwidth = d.getBBox().width;
            if (dwidth > max_gc_name_width) {
                max_gc_name_width = dwidth;
            };
        });
        sequence_start_x = max_genome_title_width + max_caller_width + max_gc_name_width + (HORIZONTAL_PADDING * 4);
        $('.sequence').attr('x', sequence_start_x);
    } else {
        sequence_start_x = max_genome_title_width + max_caller_width + (HORIZONTAL_PADDING * 3);
        $('.sequence').attr('x', sequence_start_x);
    }

    var wrap_length = parseInt($('#wrap_length').val());
    var font_size = parseInt($('#font_size').val());

    var max_sequence_width = wrap_length * (font_size * 0.6);

    var total_width = max_caller_width + max_sequence_width;
    if (mode === 'structure') {
        total_width += max_gc_name_width + HORIZONTAL_PADDING;
    }

    $('.sequenceBackground').attr('width', total_width + (HORIZONTAL_PADDING * 3));

    bbox = svg.getBBox();
    svg.setAttribute('width', bbox.width);
    svg.setAttribute('height', bbox.height);
}


function removeGeneChart() {
  var node = document.getElementById("gene-arrow-chart");
  if (node && node.parentNode) {
    node.parentNode.removeChild(node);
  }
}

function scrollTableLeft() {
    document.querySelector('.gc-acc-main').scrollBy({
        left: -100,
        behavior: 'smooth'
    });
}

function scrollTableRight() {
    document.querySelector('.gc-acc-main').scrollBy({
        left: 100,
        behavior: 'smooth'
    });
}

async function loadPSGCData(psgc_name) {
    try {        
        var response = await $.ajax({
            type: 'GET',
            cache: false,
            url: '/data/get_psgc_data/' + psgc_name
        });

        if (typeof response === 'string') {
            response = JSON.parse(response);
        }

        if (response.status === 0 && response.data) {
            return response.data;
        } else {
            console.error('Error response:', response.message);
            return null;
        }
    } catch (error) {
        console.error('Error in loadPSGCData:', error);
    }
}

async function loadGCTypeInPSGCData(psgc_name) {
    try {
        const response = await $.ajax({
            type: 'GET',
            cache: false,
            url: '/data/get_psgc_type_data/' + psgc_name
        });

        if (response.status === 0 && response.data) {
            return response.data;
        } else {
            console.error('Error response from type data:', response.message);
            return null;
        }
    } catch (error) {
        console.error('Error in loadPSGCTypeData:', error);
        return null;
    }
}

function get_gene_functions_table_html_for_structure(psgc_data, selected_gc_id) {
    for (const psgc_id in psgc_data) {
        var genes = psgc_data[psgc_id].filter(gene => 
            gene.gene_cluster_id === selected_gc_id
        );

        if (genes && genes.length > 0) {
            let functions_table_html = '<span class="popover-close-button" onclick="$(this).closest(\'.popover\').popover(\'hide\');"></span>';

            let gc_type = 'Unknown';
            if (gc_type_in_psgc) {
                for (let psgc_id in gc_type_in_psgc) {
                    if (gc_type_in_psgc[psgc_id][selected_gc_id]) {
                        gc_type = gc_type_in_psgc[psgc_id][selected_gc_id].charAt(0).toUpperCase() + gc_type_in_psgc[psgc_id][selected_gc_id].slice(1);
                        break;
                    }
                }
            }

            functions_table_html += '<h2>Gene Cluster Information</h2>';
            functions_table_html += '<table class="table table-striped" style="width: 100%; text-align: center;">';
            functions_table_html += '<tbody>';
            functions_table_html += '<tr><th style="text-align: center;">Gene Cluster ID</th><td>' + selected_gc_id + '</td></tr>';
            functions_table_html += '<tr><th style="text-align: center;">Gene Cluster Type</th><td>' + gc_type + '</td></tr>';
            functions_table_html += '</tbody></table>';

            functions_table_html += '<h3>Genes in this cluster</h3>';
            functions_table_html += '<div style="max-height: 400px; overflow-y: auto;">';
            functions_table_html += '<table class="table table-striped" style="width: 100%;">';
            functions_table_html += '<thead><tr>' + '<th>Gene Caller ID</th>' + '<th>Genome Name</th>' + '<th>Alignment Summary</th>' + '</tr></thead>';
            functions_table_html += '<tbody>';

            genes.forEach(gene => {
                functions_table_html += '<tr>' +
                    '<td>' + gene.gene_callers_id + '</td>' +
                    '<td>' + gene.genome_name + '</td>' +
                    '<td><div style="max-width: 200px; overflow-x: auto;"><code>' + 
                        (gene.alignment_summary || 'No alignment data') + '</code></div></td>' +
                    '</tr>';
            });

            functions_table_html += '</tbody></table>';
            functions_table_html += '</div>';

            return functions_table_html;
        }
    }

    return '<p>No gene information available</p>';
}
