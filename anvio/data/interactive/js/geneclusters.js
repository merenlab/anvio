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

            document.getElementById("header").innerHTML = "<strong>" + gene_cluster_name + "</strong> with " + gene_caller_ids.length + " genes detailed <br /><small><small>" + prev_str + position + next_str + "</small></small>";

            if (typeof localStorage.state === 'undefined')
            {
                alert("Something went wrong, couldn't access to localStorage");
            }
            else
            {
                initializeCheckBoxes();
                createDisplay();
                $('.loading-screen').hide();
            }
        }
    });

    // gc_key = "AAI_avg"
    // gc_id = "Cluster_00000927"

    // loadGCAdditionalData(gc_id, gc_key);

}
function loadGCAdditionalData(gc_id, gc_key){

    $.ajax({
        type: 'POST',
        cache: false,
        url: '/data/get_additional_gc_data/' + gc_id + '/' + gc_key,
        success: function(data){
            if (data['status'] == 1){
                // Put in table element
                console.log(data)
            }
        }
    })
}

function createDisplay(){
    var sequence_wrap_val = parseInt($('#wrap_length').val());
    var sequence_font_size_val = parseInt($('#font_size').val());

    var sequence_wrap = (isNumber(sequence_wrap_val) && sequence_wrap_val > 0) ? sequence_wrap_val : 1042;
    var sequence_font_size = (isNumber(sequence_font_size_val) && sequence_font_size_val > 0) ? sequence_font_size_val : 12;

    var svg = document.getElementById('svg');

    // clear content
    while (svg.firstChild) {
        svg.removeChild(svg.firstChild);
    }

    var y_cord = 0;
    var offset = 0;

    var acid_sequences = [];
    var order = {};
    var count = 0;
    for (var layer_id = 0; layer_id < state['layer-order'].length; layer_id++)
    {
        var layer = state['layer-order'][layer_id];

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

    while (true)
    {
        ignoreTable = true;
        var fragment = document.createDocumentFragment();
        for (var layer_id = 0; layer_id < state['layer-order'].length; layer_id++)
        {
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

                var text = document.createElementNS('http://www.w3.org/2000/svg', 'text');
                text.setAttribute('x', 0);
                text.setAttribute('y', sub_y_cord);
                text.setAttribute('font-size', sequence_font_size);
                text.setAttribute('font-family',"monospace");
                text.setAttribute('font-weight', '100');
                text.setAttribute('style', 'alignment-baseline:text-before-edge');
                text.setAttribute('class', 'sequence');

                _sequence = sequence.substr(offset, sequence_wrap);
                for (var _letter_index=0; _letter_index < _sequence.length; _letter_index++) {
                    var tspan = document.createElementNS('http://www.w3.org/2000/svg', 'tspan');
                    var index = _letter_index+offset;
                    var num = order[caller_id];
                    var acid = _sequence[_letter_index];
                    var dict = coded_positions[index][num];
                    tspan.setAttribute('fill', dict[acid]);
                    tspan.style.fontWeight = 'bold';
                    tspan.appendChild(document.createTextNode(acid));
		            tspan.setAttribute('style', 'alignment-baseline:text-before-edge');
                    text.appendChild(tspan);
                }

                fragment.appendChild(text);

                sub_y_cord = sub_y_cord + sequence_font_size * 1.5;

                if (offset < sequence.length) {
                    ignoreTable = false;
                }
            });

            y_cord = y_cord + parseFloat(rect.getAttribute('height')) + 5;

        }
        if (ignoreTable) {
            break;
        } else {
            svg.appendChild(fragment);
        }
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
    $('.callerTitle').each(function(i, d) {
        dwidth = d.getBBox().width;
        if (dwidth > max_caller_width) {
            max_caller_width = dwidth;
        };
    });

    $('.sequence').attr('x', max_genome_title_width + max_caller_width + (HORIZONTAL_PADDING * 3));

    var max_sequence_width = 0;
    $('.sequence').each(function(i, d) {
        dwidth = d.getBBox().width;
        if (dwidth > max_sequence_width) {
            max_sequence_width = dwidth;
        };
    });

    $('.sequenceBackground').attr('width', max_caller_width + max_sequence_width + (HORIZONTAL_PADDING * 3));

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
