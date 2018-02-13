/**
 * Javascript library to visualize anvi'o gene clusterss
 *
 *  Author: A. Murat Eren <a.murat.eren@gmail.com>
 *  Credits: Özcan Esen
 *  Copyright 2016, The anvio Project
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


/****************************************************************************
 * Algorithm for determining color coding
 * Coded by: Mahmoud Yousef
 ****************************************************************************/

//Returns the length of the largest object in an array
function maxLength(arr){
    var length = 0;
    for (var i = 0; i < arr.length; i++)
    {
        if (arr[i].length > length){
            length = arr[i].length;
        }
    }
    return length;
}

//Use an array of dictionaries. The key is the amino acid, and the data is the color, or "null" if none
//Input: a 2D array of the amino acids
//Returns: a 2D array of dictionaries
function determineColor(sequences_array){
    var cols = sequences_array.length;
    var results = [];
    for (var n = 0; n < cols; n++){
        results.push([]);
    }

    var column = [];
    for(var c = 0; c < cols; c++){
        column = sequences_array[c];
        _positions = colorAlgorithm(column)
        results[c] = _positions
    }
    return results;
}


//compares all of the aa's and assigns colors or "null" to the respective amino acid
//input: the "column" array from determineColor
//output: the array with all of the color assignments\
function colorAlgorithm(positions){
    var _positions = []
    for (aa in positions){
        if (checked(positions[aa]) && aboveThreshold(positions, positions[aa])) {
            _positions[aa] = color(positions, positions[aa]); 	
        } else{
            var dict = {}
            dict[positions[aa]] = "BLACK";
            _positions[aa] = dict;
        }
    }
    return _positions;
}

function checked(letter){
    if (letter == undefined || letter == null){
        return false;
    }
    box =  document.getElementById(letter)
    if (box == null){
        return false;
    }
    return box.checked;
}


//does the actual comparisons
//This checks for amino acid conservation by common characteristics
function aboveThreshold(positions, aa) { 
    var number = 0;
    for (acid in positions) {
        if (acid === '') {
            continue;
        }
        number++;
    }
    var count = 0.0;
    var count2 = 0.0;
    var count3 = 0.0;
    
    var letter = aa
    switch (letter) {
        case "A":
        case "I":
        case "L":
        case "M":
        case "F":
        case "W":
        case "V":
            for (amino in positions){
                var acid = positions[amino]
                if (acid == "W" || acid == "L" || acid == "V" || acid == "I" || acid == "M" || acid == "A" || acid == "F" || acid == "C" || acid == "H" || acid == "P") {
                    count = count + 1;
                }
                }
              if ( (100 * count) / number >= 60) {
                return true;
              }
              break;
        case "P":
        case "G":
              return true;
              break; //P and G have a 0% threshold
        case "R":
        case "K":
          for (amino in positions){
              var acid = positions[amino]
                if (acid == "K" || acid == "R") {
                    count = count + 1;
                    count2 = count2 + 1;
                }
                else if (acid == "Q") {
                    count2 = count2 + 1;
                }
              }
              if ((100 * count) / number >= 60 || (100 * count2) / number >= 80) {
                return true;
              }
              break;
        case "N":
          for (amino in positions){
              var acid = positions[amino]
                if (acid == "N"){
                    count = count + 1;
                    count2 = count2 + 1;
                }
                if (acid == "Y"){
                    count2 = count2 + 1;
                }
              }
              if ((count * 100) / number >= 50 || (count2 * 100) / number >= 85) {
                    return true;
              }
              break;
        case "C":
          for (amino in positions){
              var acid = positions[amino]
                if (acid == "W" || acid == "L" || acid == "V" || acid == "I" || acid == "M" || acid == "A" || acid == "F" || acid == "C" || acid == "H" || acid == "P") {
                    count = count + 1;
                }
                else if (acid == "C"){
                    count2 = count2 + 1;
                }
              }
              if ((count * 100) / number >= 60 || count2 == number) {
                return true;
              }
              break;
        case "Q":
          for (amino in positions){
              var acid = positions[amino]
                if (acid == "K" || acid == "R"){
                    count++;
                    count3++;
                }else if (acid == "Q" || acid == "E"){
                    count2++;
                    count3++;
                }
              }
              if ((count * 100) / number >= 60 || (count2 * 100) / number >= 50 || (count3 * 100) / number >= 85){
                return true;
              }
              break;
        case "E":
          for (amino in positions){
              var acid = positions[amino]
                if (acid == "K" || acid == "R"){
                    count++;
                }else if (acid == "Q" || acid == "E"){
                    count2++;
                    count3++;
                }else if (acid == "D"){
                    count3++;
                }
              }
              if ((count * 100) / number >= 60 || (count2 * 100) / number >= 50 || (count3 * 100) / number >= 85){
                return true;
              }
              break;
        case "D":
          for (amino in positions){
              var acid = positions[amino]
                if (acid == "K" || acid == "R"){
                    count++;
                    count2++;
                }else if (acid == "Q"){
                    count2++;
                }else if (acid == "E" || acid == "D"){
                    count3++;
                }
              }
              if ((count * 100) / number >= 60 || (count2 * 100) / number >= 85 || (count3 * 100) / number >= 50){
                return true;
              }
              break;
        case "H":
        case "Y":
          for (amino in positions){
              var acid = positions[amino]
                if (acid == "W" || acid == "L" || acid == "V" || acid == "I" || acid == "M" || acid == "A" || acid == "F" || acid == "C" || acid == "H" || acid == "P") {
                    count++;
                    count2++;
                }else if (acid == "Y" || acid == "Q"){
                    count2++;
                }
              }
              if ((count * 100) / number >= 60 || (count2 * 100) / number >= 85) {
                return true;
              }
              break;
        case "S":
        case "T":
          for (amino in positions){
              var acid = positions[amino]
                if (acid == "W" || acid == "L" || acid == "V" || acid == "I" || acid == "M" || acid == "A" || acid == "F" || acid == "C" || acid == "H" || acid == "P") {
                    count++;
                }else if (acid == "S" || acid == "T"){
                    count2++;
                }
              }
              if ((count * 100) / number >= 60 || (count2 * 100) / number >= 50) {
                return true;
              }
              break;
        default: break;
    }
    return false;	

}

function color(positions, aa){
    var x = '';
    switch(aa){
        case "A":
        case "I":
        case "L":
        case "M":
        case "F":
        case "W":
        case "V":
              x = "BLUE";
              break;
        case "R":
        case "K":
              x = "RED";
              break;
        case "N":
        case "Q":
        case "S":
        case "T":
              x = "GREEN";
              break;
        case "E":
        case "D":
              x = "MAGENTA";
              break;
        case "G":
            x = "ORANGE";
              break;
        case "H":
        case "Y":
              x = "DARKTURQUOISE";
              break;
        case "P": 
            x = "YELLOW";
              break;
        case "C":
                check: {
                for (acid in positions){
                    if (positions[acid] != "C"){
                        x = "BLUE";
                        break check;
                      }
                  }
                x = "HOTPINK";
              }
              break;
        default: x = null;
    }
    var dict = {}
    dict[aa] = x;
    return dict;
}

/********************************************************************
 * End color coding algorithm
 ********************************************************************/


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

    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/geneclusters/' + gene_cluster_name,
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
                next_str = '<a onclick="sessionStorage.state = JSON.stringify(state, null, 4);" href="' + generate_inspect_link('geneclusters', next_gene_cluster_name) +'" '+target_str+'> | next &gt;&gt;&gt;</a>';

            if(previous_gene_cluster_name)
                prev_str = '<a onclick="sessionStorage.state = JSON.stringify(state, null, 4);" href="' + generate_inspect_link('geneclusters', previous_gene_cluster_name) +'" '+target_str+'>&lt;&lt;&lt; prev | </a>';

            document.getElementById("header").innerHTML = "<strong>" + gene_cluster_name + "</strong> with " + gene_caller_ids.length + " genes detailed <br /><small><small>" + prev_str + position + next_str + "</small></small>";

            if (typeof sessionStorage.state === 'undefined')
            {
                alert("Something went wrong, couldn't access to sessionStorage");
            }
            else
            {
                // backup the state, if user changes the page (prev, next) we are going to overwrite it.
                state = JSON.parse(sessionStorage.state);
                initializeCheckBoxes();
                createDisplay();
                $('.loading-screen').hide();
            }
        }
    });

}

function initializeCheckBoxes(){
    var svg = document.getElementById('svg');
    var container = document.getElementById('display-conservation-controls');

    var labels = [A = "A", C = "C", D = "D", E = "E", F = "F", G = "G", H = "H",
         I = "I", K = "K", L = "L", M = "M", N = "N", P = "P", Q = "Q", R = "R",
          S = "S", T = "T", V = "V", W = "W", Y = "Y"]

    for (i = 0; i < labels.length; i++){
        var word = String(labels[i])
        var box = document.createElement('input');
        box.type = "checkbox";
        box.name = word;
        box.value = word;
        box.id = word;
        if (i == 0){
            box.style = "margin-left:70px;"
        } else {
            box.style = "margin-left:5px;"
        }
        box.onclick = ( function() {
                        return createDisplay();
                         } )
        
        var label = document.createElement('label')
        label.htmlFor = word;
        label.appendChild(document.createTextNode(word));
        
        container.appendChild(box);
        container.appendChild(label);   
        box.checked = true 
    }

    var all = document.createElement("button")
    all.innerHTML = "check all"
    all.style = "margin-right:20px;"
    all.onclick  = (function() {
        for (i = 0; i < labels.length; i++){
            var letter = String(labels[i])
            document.getElementById(letter).checked = true
        }
        createDisplay();
    } )
    container.appendChild(all)

    var none = document.createElement("button")
    none.innerHTML = "uncheck all"
    none.onclick  = (function() {
        for (i = 0; i < labels.length; i++){
            var letter = String(labels[i])
            document.getElementById(letter).checked = false
        }
        createDisplay();
    } )
    container.appendChild(none)

}

function createDisplay(){
    var sequence_wrap_val = parseInt($('#wrap_length').val());
    var sequence_font_size_val = parseInt($('#font_size').val());

    var sequence_wrap = (isNumber(sequence_wrap_val) && sequence_wrap_val > 0) ? sequence_wrap_val : 140;
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
            order[layer] = count;
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
            text.setAttribute('font-size', "24px");
            text.setAttribute('font-family', "Lato, Arial");
            text.setAttribute('font-weight', '300');
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
                text.setAttribute('font-family', "Lato, Arial");
                text.setAttribute('font-weight', '300');
                text.setAttribute('style', 'alignment-baseline:text-before-edge; cursor: pointer;');
                text.setAttribute('class', 'callerTitle');
                text.setAttribute('gene-callers-id', caller_id);
                text.setAttribute('genome-name', layer);
                text.setAttribute('data-content', get_gene_functions_table_html_for_pan(caller_id, layer) + '')
                text.setAttribute('data-toggle', 'popover');

                text.appendChild(document.createTextNode(caller_id));
                fragment.appendChild(text);

                var text = document.createElementNS('http://www.w3.org/2000/svg', 'text');
                text.setAttribute('x', 0);
                text.setAttribute('y', sub_y_cord);
                text.setAttribute('font-size', sequence_font_size);
                text.setAttribute('font-family', "monospace");
                text.setAttribute('font-weight', '100');
                text.setAttribute('style', 'alignment-baseline:text-before-edge');
                text.setAttribute('class', 'sequence');

                _sequence = sequence.substr(offset, sequence_wrap);
                for (var _letter_index=0; _letter_index < _sequence.length; _letter_index++) {
                    var tspan = document.createElementNS('http://www.w3.org/2000/svg', 'tspan');
                    var index = _letter_index+offset;
                    var num = order[layer];
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

    $('[data-toggle="popover"]').popover({"html": true, "trigger": "click", "container": "body", "viewport": "body", "placement": "top"});

    // workaround for known popover bug
    // source: https://stackoverflow.com/questions/32581987/need-click-twice-after-hide-a-shown-bootstrap-popover
    $('body').on('hidden.bs.popover', function (e) {
      $(e.target).data("bs.popover").inState.click = false;
    });

    $('[data-toggle="popover"]').on('shown.bs.popover', function (e) {
      var popover = $(e.target).data("bs.popover").$tip;
      
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
