//ANCHOR - NODE FUNCTIONS
async function nodeinfo(e, data, group_dict, mapAS, gene_cluster_id='') {
  // console.log(gene_cluster_id)
  var id = e.id;
  var element = document.getElementById(id);

  if (element.getAttribute('class') == 'group') {
    gene_cluster_context = id;

    if (!gene_cluster_id) {
      gene_cluster_id = group_dict[id][0]
    }
  } else {
    gene_cluster_id = id;
    gene_cluster_context = null;
  }

  $('#InfoModalBody').empty()
  var bodyinfo = $('<div class="card-body overflow-scroll"></div>')
  // var closeBtn = $('<div class="row justify-content-end"><i class="btn d-flex col-1 justify-content-end text-danger fs-2x bi bi-x-circle-fill" data-dismiss="modal"></i></div>')
  // $('#InfoModalBody').append(closeBtn)
  $('#InfoModalBody').append(bodyinfo)

  var all_info = await get_gene_cluster_display_tables(gene_cluster_id, gene_cluster_context, data, group_dict, mapAS, add_align=1)

  bodyinfo.append(all_info)

  if (element.getAttribute('class') == 'group') {

    var container = document.querySelector("#InfoModalBody");
    var elem = container.querySelectorAll('.group_choice');
    
    elem.forEach(el => {
    
      el.onclick = function(){
        nodeinfo(e, data, group_dict, mapAS, gene_cluster_id=el.getAttribute("name_id"))
      }
    })
  }

  $('#InfoModal').modal('show');
}

async function get_gene_cluster_display_tables(gene_cluster_id, gene_cluster_context, data, group_dict, mapAS, add_align) {
    // The purpose of this function is to build HTML formatted tables to give access to
    // the details of a gene cluster. The parameters here are,
    //
    //   - `gene_cluster_id`: A singular gene cluster id to be detailed,
    //   - `gene_cluster_context`: A list of gene cluster ids that occur in the same context
    //      with `gene_cluster_id` (either they were binned together, or they were in the same
    //      group of gene clusters).
    //   - `data`: the primary data object from the JSON input
    //
    // If everything goes alright, this function will return a table that can be displayed in
    // any modal window.
    ///////////////////////////////////////////////////////////////////////////////////////////
    // BUILD CONTEXT
    // if this is a gene cluster that is a part of a context, then all others are going to be
    // shown here
    ///////////////////////////////////////////////////////////////////////////////////////////

    var gene_cluster_context_table = get_gene_cluster_context_table(gene_cluster_id, gene_cluster_context, data, group_dict, add_align);

    ///////////////////////////////////////////////////////////////////////////////////////////
    // BUILD BASIC INFORMATION TABLE
    ///////////////////////////////////////////////////////////////////////////////////////////

    var basic_info_table = get_gene_cluster_basics_table(gene_cluster_id, gene_cluster_context, data, add_align);

    ///////////////////////////////////////////////////////////////////////////////////////////
    // BUILD FUNCTIONS TABLE
    ///////////////////////////////////////////////////////////////////////////////////////////

    var basic_layer_table = get_layer_data(gene_cluster_id, data, add_align);

    var functions_table = await get_gene_cluter_functions_table(gene_cluster_id, data, add_align);

    ///////////////////////////////////////////////////////////////////////////////////////////
    // RETRIEVE AND BUILD SEQUENCE ALIGNMENTS
    ///////////////////////////////////////////////////////////////////////////////////////////
 
    if (add_align == 1) {

      var alignment = {}

      // if (gene_cluster_id != 'start' && gene_cluster_id != 'stop') {
      for (var genome of Object.keys(data['nodes'][gene_cluster_id]['gene_calls'])) {
        
        // if ($('#flex' + genome).prop('checked') == true){
        
        var genecall = data['nodes'][gene_cluster_id]['gene_calls'][genome]
        var name = data['nodes'][gene_cluster_id]['gene_cluster']
        alignment[genome] = [genecall, name]

        // }
      }
      // }

      gene_cluster_sequence_alignments_table = await appendalignment(gene_cluster_id, alignment, mapAS)

      ///////////////////////////////////////////////////////////////////////////////////////////
      // MERGE ALL AND RETURN
      ///////////////////////////////////////////////////////////////////////////////////////////

      return gene_cluster_context_table + basic_info_table + basic_layer_table + functions_table + gene_cluster_sequence_alignments_table

    } else {

      return gene_cluster_context_table + basic_info_table + basic_layer_table + functions_table

    }
}

function get_gene_cluster_consensus_functions(gene_cluster_name) {

  var func = new Object();
  func['genecluster'] = gene_cluster_name

  var d = $.ajax({
    url: "/pangraph/function",
    type: "POST",
    data: JSON.stringify(func),
    contentType: "application/json",
    dataType: "json",
    error: function(){
      console.log('Error while attempting to fetch function.')
    },
    success: function(){
      console.log('Successfully fetched function.')
    }
  })

  return d
}

function get_layer_data(gene_cluster_id, data, add_align) {

  var layers = Object.keys(data['meta']['layers'])
  var basic_info = {}

  for (var layer_name of layers) {
    if ($('#flex' + layer_name).prop('checked') == true){

      basic_info[layer_name] = data['nodes'][gene_cluster_id]['layer'][layer_name]

    }
  }

  if (add_align == 1) {
    basic_layer_table = `<p class="modal_header mt-0">Layers</p>`;
    basic_layer_table += `<table class="table table-striped table-bordered sortable" id="node_layers_table">`;
  } else {
    basic_layer_table = ''
    basic_layer_table += `<table class="table table-striped table-bordered sortable">`;
  }

  basic_layer_table += `<tbody>`;
  basic_layer_table += `<thead class="thead-light"><tr>`;
  for (const [key, value] of Object.entries(basic_info)) {
    basic_layer_table += `<th scope="row">` + key + `</th>`;
  }
  basic_layer_table += `</tr></thead><tbody>`;

  basic_layer_table += `<tbody><tr>`;
  for (const [key, value] of Object.entries(basic_info)) {
    basic_layer_table += `<td>` + value + `</td>`;
  }
  basic_layer_table += `</tbody></tr></table>`;

  return basic_layer_table;

}

function get_gene_cluster_basics_table(gene_cluster_id, gene_cluster_context, data, add_align) {
    // first, learn a few basics about the gene cluster to be displayed
    var x_pos = data['nodes'][gene_cluster_id]['position'][0]
    // var position_in_graph = x_pos + " / " + (data["infos"]["meta"]["global_x_offset"] - 1);
    var position_in_graph = x_pos;
    // var num_contributing_genomes = Object.keys(data['elements']['nodes'][gene_cluster_id]['genome']).length + " / " + (data['infos']['num_genomes']);
    var num_contributing_genomes = Object.keys(data['nodes'][gene_cluster_id]['gene_calls']).length;
    var gene_cluster_name = data['nodes'][gene_cluster_id]['gene_cluster']

    if (gene_cluster_context == null){
      var context = gene_cluster_id
    } else {
      var context = gene_cluster_context
    }

    basic_info = {'ID': gene_cluster_id, 'Gene Cluster': gene_cluster_name, 'Contributing Genomes': num_contributing_genomes, 'Position in Graph': position_in_graph}
    // build the basic information table
    if (add_align == 1) {
      basic_info_table = `<p class="modal_header mt-0">Basics</p>`;
      basic_info_table += `<table class="table table-striped table-bordered sortable" gc_context="` + context + `" gc_pos="` + x_pos + `" gc_id="` + gene_cluster_id + `" id="node_basics_table">`;
    } else {
      basic_info_table = ''
      basic_info_table += `<table class="table table-striped table-bordered sortable">`;
    }

    basic_info_table += `<tbody>`;
    basic_info_table += `<thead class="thead-light"><tr>`;
    for (const [key, value] of Object.entries(basic_info)) {
        basic_info_table += `<th scope="row">` + key + `</th>`;
    }
    basic_info_table += `</tr></thead><tbody>`;

    basic_info_table += `<tbody><tr>`;
    for (const [key, value] of Object.entries(basic_info)) {
        basic_info_table += `<td>` + value + `</td>`;
    }
    basic_info_table += `</tbody></tr></table>`;

    return basic_info_table;
}

async function get_gene_cluter_functions_table(gene_cluster_id, data, add_alig) {
    
    if (add_align == 1) {
      functions_table = `<p class="modal_header">Consensus functional annotations</p>`;
      functions_table += `<table class="table table-striped sortable" gc_id="` + gene_cluster_id + `" id="node_functions_table">`;
    } else {
      functions_table = ''
      functions_table += `<table class="table table-striped sortable">`;
    }
    functions_table += `<thead><tr>`;
    functions_table += `<th scope="col">Source</th>`;
    functions_table += `<th scope="col">Accession</th>`;
    functions_table += `<th scope="col">Function</th>`;
    functions_table += `</tr></thead><tbody>\n\n`;

    var gene_cluster_name = data['nodes'][gene_cluster_id]['gene_cluster']
    var d = await get_gene_cluster_consensus_functions(gene_cluster_name);

    // console.log(d)

    function_sources = Object.keys(d).sort();
    // console.log(function_sources);
    for (source of function_sources) {
        value = d[source];
        // console.log(source)
        var accession = value['accession']
        accession === undefined | accession === null ? accession = '' : accession = accession
        var func = value['function']
        func === undefined | func === null ? func = '' : func = func

        functions_table += `<tr>`;
        functions_table += `<td>` + source + `</td>`;
        functions_table += `<td>` + accession + `</td>`;
        functions_table += `<td>` + func + `</td>`;
        functions_table += `</tr>`;
    }
    functions_table += `</tbody></tr></table>\n\n`;

    return functions_table
}


function get_gene_cluster_context_table(gene_cluster_id_current, gene_cluster_context, data, group_dict, add_align) {
    
    if (gene_cluster_context == null){
        return '';
    } else {
      var group_context = []
      for (item in group_dict[gene_cluster_context]){
        group_context.push(group_dict[gene_cluster_context][item])
      }
    }

    // console.log(gene_cluster_id_current, gene_cluster_context, group_context)
    if (add_align == 1) {
      gene_cluster_context_table = `<p class="modal_header">Gene cluster context</p>`;
    }else {
      gene_cluster_context_table = ''
    }
    gene_cluster_context_table += `<div class="gene_cluster_context_items">`;
    for(index in group_context) {
        gene_cluster_id = group_context[index];
        gene_cluster_name = data['nodes'][gene_cluster_id]['gene_cluster']

        if (gene_cluster_id == gene_cluster_id_current){
            gene_cluster_context_table += `<span class="gene_cluster_id gene_cluster_id_current">` + gene_cluster_name + `</span>`;
        } else {
            gene_cluster_context_table += `<span class="gene_cluster_id"><a class="btn border-0 m-0 p-0 align-baseline group_choice" context="` + add_align + `" group="` + gene_cluster_context + `" name_id="` + gene_cluster_id + `">` + gene_cluster_name + `</a></span>`;
        }
    }
    gene_cluster_context_table += `</div>`;

    return gene_cluster_context_table;

}

function fetchalignment(alignment) {

  var d = $.ajax({
    url: "/pangraph/alignment",
    type: "POST",
    data: JSON.stringify(alignment),
    contentType: "application/json",
    dataType: "json",
    error: function(){
      console.log('Error while attempting to fetch alignment.')
    },
    success: function(){
      console.log('Successfully fetched alignment.')
    }
  })

  return d

};

async function appendalignment(gene_cluster_id, alignment, mapAS) {
    if (Object.keys(alignment).length == 0) {
        return ''
    }

    alignments_table = `<p class="modal_header">Sequence alignments</p>`;
    alignments_table += `<div class="scroll-wrapper"><table class="table sortable" gc_id="` + gene_cluster_id + `" id="node_sequence_alignments_table">`;
    alignments_table += `<thead class="thead-dark gc-table-header"><tr>`;
    alignments_table += `<th class="position-sticky" style="left:0px; z-index:2;" scope="col">Genome</th>`;
    alignments_table += `<th scope="col">Gene Call</th>`;
    alignments_table += `<th scope="col"><span id="th-sequence">Sequence</span></th>`;
    alignments_table += `</tr></thead><tbody>\n\n`;

    // console.log(alignment)

    var d = await fetchalignment(alignment)

    // console.log(d)

    for (var [genome, value] of Object.entries(d)) {
      var colored = value[1].replace(/A|R|N|D|C|Q|E|G|H|I|L|K|M|F|P|S|T|W|Y|V|-/gi, function(matched){return mapAS[matched];});

      alignments_table += `<tr>`
      alignments_table += `<td id="td-genome-cell">` + genome + `</td>`
      alignments_table += `<td id="td-value-cell">` + value[0] + `</a></td>`
      alignments_table += `<td id="gc-alignment-font"><div class="scrollable-content">` + colored + `</div></td>`
      alignments_table += `</tr>`
    }

    alignments_table += `</tbody></table></div>\n\n`

    // I guess this shouldn't return anything but insert the results into a context?
    // lol sorry if I screwed these thigns up badly :) IDK what I'm doing :p
    return alignments_table;
}






















//ANCHOR - BIN FUNCTIONS
function marknode(e, data, binid, bins, genome_size, group_dict){

  var bincolor = document.getElementById(binid + 'color').value
  var id = e.id;
  var current = ''

  var binkeys = Object.keys(bins)
    for (var key of binkeys) {
      if (bins[key].includes(id)) {
        current = key
        break;
    }
  }

  var core_color = $('#core_color')[0].value;
  var paralog_color = $('#paralog_color')[0].value;
  var singleton_color = $('#singleton_color')[0].value;
  var accessory_color = $('#accessory_color')[0].value;
  var rearranged_color = $('#rearranged_color')[0].value;
  var trna_color = $('#trna_color')[0].value;

  if ($('#flexsaturation').prop('checked') == true){
    var saturation = 1
  } else {
    var saturation = 0
  }

  if (e.getAttribute('class') == 'group') {

    var group = group_dict[id]
    var node_name = group[0]
    var node = data['nodes'][node_name];
    var node_type = node['type']
    var genome = Object.keys(node['gene_calls']).length;

  } else if (e.getAttribute('class') == 'node') {

    var node = data['nodes'][id];
    var node_type = node['type']
    var genome = Object.keys(node['gene_calls']).length;

  }

  if (node_type == 'core'){
    var node_color = core_color
  } else if (node_type == 'rearrangement') {
    var node_color = rearranged_color
  } else if (node_type == 'accessory') {
    var node_color = accessory_color
  } else if (node_type == 'duplication') {
    var node_color = paralog_color
  } else if (node_type == 'singleton') {
    var node_color = singleton_color
  } else if (node_type == 'rna') {
    var node_color = trna_color
  }

  if (current === binid) {

    if (saturation == 1){
      e.setAttribute("fill", lighter_color('#ffffff', node_color, genome / genome_size))
    } else {
      e.setAttribute("fill", node_color)
    }
    
    bins[binid] = bins[binid].filter(item => item !== id)
    $('#' + binid + 'value')[0].value = bins[binid].length

  } else if (current === '') {

    if (saturation == 1){
      e.setAttribute("fill", lighter_color('#ffffff', bincolor, genome / genome_size))
    } else {
      e.setAttribute("fill", bincolor)
    }

    bins[binid].push(id)
    $('#' + binid + 'value')[0].value = bins[binid].length

  } else {

    if (saturation == 1){
      e.setAttribute("fill", lighter_color('#ffffff', bincolor, genome / genome_size))
    } else {
      e.setAttribute("fill", bincolor)
    }

    bins[current] = bins[current].filter(item => item !== id)
    bins[binid].push(id)
    $('#' + binid + 'value')[0].value = bins[binid].length
    $('#' + current + 'value')[0].value = bins[current].length

  }

  return bins
}


//ANCHOR - SVG FUNCTIONS
function deg2rad(degrees)
{
  return degrees * Math.PI/180;
}

function downloadBlob(blob, name) {

  var blobUrl = URL.createObjectURL(blob);
  var link = document.createElement("a");

  link.href = blobUrl;
  link.download = name;

  document.body.appendChild(link);

  link.dispatchEvent(
    new MouseEvent('click', {
      bubbles: true,
      cancelable: true,
      view: window
    })
  );

  document.body.removeChild(link);
}

// NOTE - From https://coderwall.com/p/z8uxzw/javascript-color-blender
function int_to_hex(num) {
    var hex = Math.round(num).toString(16);
    if (hex.length == 1)
        hex = '0' + hex;
    return hex;
}

function lighter_color(color1, color2, percentage, threshold=0.25) {

  // console.log(color1, color2, percentage)

  percentage = threshold + (1 - threshold) * percentage

  // var color3 = $.xcolor.gradientlevel(color1, color2, percentage, 1);

  color1 = color1.substring(1);
  color2 = color2.substring(1);

  color1 = [parseInt(color1[0] + color1[1], 16), parseInt(color1[2] + color1[3], 16), parseInt(color1[4] + color1[5], 16)];
  color2 = [parseInt(color2[0] + color2[1], 16), parseInt(color2[2] + color2[3], 16), parseInt(color2[4] + color2[5], 16)];

  var color3 = [
      (1 - percentage) * color1[0] + percentage * color2[0],
      (1 - percentage) * color1[1] + percentage * color2[1],
      (1 - percentage) * color1[2] + percentage * color2[2]
  ];

  color3 = '#' + int_to_hex(color3[0]) + int_to_hex(color3[1]) + int_to_hex(color3[2]);

  // console.log(color3)

  return color3
}

function create_rectangle(i_x, i_y, j_x, j_y, theta, node_distance_x, linear, color, id='') {

  if (id != '') {
    var extra = '" id="' + id
  } else {
    var extra = ''
  }

  if (linear == 0) {
    var [a_x, a_y] = transform(i_x, i_y, theta)
    var [b_x, b_y] = transform(j_x, i_y, theta)
    var [c_x, c_y] = transform(i_x, j_y, theta)
    var [d_x, d_y] = transform(j_x, j_y, theta)

    var path = $('<path d="' +
      'M ' + a_x + ' ' + a_y + ' ' +
      'A ' + (i_y) + ' ' + (i_y) + ' 0 0 0 ' + b_x + ' ' + b_y + ' ' +
      'L' + d_x + ' ' + d_y +  ' ' +
      'A ' + (j_y) + ' ' + (j_y) + ' 0 0 1 ' + c_x + ' ' + c_y + ' ' +
      'L' + a_x + ' ' + a_y +
      '" fill="' + color + extra + '" stroke="" stroke-width="0"/>')
  } else {
    var [a_x, a_y] = [i_x * node_distance_x, -i_y]
    var [b_x, b_y] = [j_x * node_distance_x, -i_y]
    var [c_x, c_y] = [i_x * node_distance_x, -j_y]
    var [d_x, d_y] = [j_x * node_distance_x, -j_y]
  
    var path = $('<path d="' +
      'M ' + a_x + ' ' + a_y + ' ' +
      'L ' + b_x + ' ' + b_y + ' ' +
      'L ' + d_x + ' ' + d_y +  ' ' +
      'L ' + c_x + ' ' + c_y + ' ' +
      'L' + a_x + ' ' + a_y +
      '" fill="' + color + extra + '" stroke="" stroke-width="0"/>')
  }

  return(path)
}

function transform(x, y, theta) {
  var circle_x = y * Math.sin(deg2rad(theta * (x+0.5)))
  var circle_y = y * Math.cos(deg2rad(theta * (x+0.5)))
  return [circle_x, circle_y]
}

function pickcolor (edgecoloring, genomes) {

  var array = []

  for (var name of genomes) {
    array.push(edgecoloring[name])
    var sortedArray = array.sort(function(a, b) {
      return a[0] - b[0];
    });
  }

  return sortedArray[0][1]
}

function draw_newick(order, item_dist, max_dist, offset, max_size, line_thickness) {

  var output = ''
  var saving_positions = {}
  var saving_ids = {}
  
  var i = 0
  
  for (var [item, start, end] of order){
  
    var start_fraction = (max_size * (start / max_dist) - max_size) - offset
    var end_fraction = (max_size * (end / max_dist) - max_size) - offset
    
    if (item != 'branching') {
      var y_value = item_dist[item]
      output += '<path d="M ' + start_fraction + ' ' + y_value + ' L ' + end_fraction + ' ' + y_value + '" stroke-width="' + line_thickness + '" stroke="black"></path>'

      if (Object.values(saving_ids).includes(start_fraction)){
        var saving_id = Object.keys(saving_ids)[Object.values(saving_ids).indexOf(start_fraction)]
        // console.log(saving_id)
        saving_positions[saving_id].push(y_value)
      } else {
        saving_ids[i] = start_fraction
        saving_positions[i] = [y_value]
        i = i + 1
      }
      if (end_fraction != max_size){
        output += '<path d="M ' + end_fraction + ' ' + y_value + ' L ' + (0 - offset) + ' ' + y_value + '" stroke-dasharray="' + line_thickness * 5 + ',' + line_thickness * 5 + '" stroke-width="' + line_thickness + '" stroke="lightgray"></path>'
      }
    } else {

      if (Object.values(saving_ids).includes(end_fraction)){
        var saving_id_main = Object.keys(saving_ids)[Object.values(saving_ids).indexOf(end_fraction)]
        sorted_positions = saving_positions[saving_id_main].sort()
        // console.log(sorted_positions)
        
        // for (var j in sorted_positions) {
        for (var j = 0; j < sorted_positions.length -1; j++) {
          
          var y_value_i = sorted_positions[j]
          var y_value_j = sorted_positions[j+1]
        
          output += '<path d="M ' + end_fraction + ' ' + y_value_i + ' L ' + end_fraction + ' ' + y_value_j + '" stroke-width="' + line_thickness + '" stroke="black"></path>'

        }
      } else {
        var saving_id_main = -1
      }

      y_value = Math.min(...sorted_positions) + (Math.max(...sorted_positions) - Math.min(...sorted_positions)) / 2
      output += '<path d="M ' + start_fraction + ' ' + y_value + ' L ' + end_fraction + ' ' + y_value + '" stroke-width="' + line_thickness + '" stroke="black"></path>'

      if (Object.values(saving_ids).includes(start_fraction)) {
        var saving_id = Object.keys(saving_ids)[Object.values(saving_ids).indexOf(start_fraction)]
        saving_positions[saving_id].push(y_value)
      } else {
        saving_ids[i] = start_fraction
        saving_positions[i] = [y_value]
        i = i + 1
      }

      if (saving_id_main != -1) {
        delete saving_ids[saving_id_main]
        delete saving_positions[saving_id_main]
      }
    }
  }
  return(output)

}

function newick_to_order(string, prior = 0) {
  var result = []
  var newick = string.replace(' ', '')

  if (newick[0] == '('){
    var bracket_open = 0
    var bracket_closed = 0

    // for (var i in newick) {
    for (var i = 0; i < newick.length; i++) {

      var sub_letter = newick[i]
      if (sub_letter == '(') {
        bracket_open += 1 
      } else if (sub_letter == ')') {
        bracket_closed += 1
      }
          
      if (bracket_open == bracket_closed) {
        
        var sub_newick = newick.slice(1, i)
        var rest = newick.slice(i)
        
        if (rest.includes(',')) {
            var parts = rest.split(',')
            if (parts[0].includes(':')) {
                var value = parts[0].split(':')[1]
            } else {
                var value = 0 
            }

            result.push(...[['branching', prior, prior + parseFloat(value)]])

            // console.log(sub_newick)
            result.push(...newick_to_order(sub_newick, prior + parseFloat(value)))
            
            var next_iter = parts.slice(1).join(',')
            // console.log(next_iter)
            result.push(...newick_to_order(next_iter, prior))
        } else {
            if (rest.includes(':')) {
                var value = rest.split(':').slice(1)
            } else {
                var value = 0 
            }
            result.push(...[['branching', prior, prior + parseFloat(value)]])
            // console.log(sub_newick)
            result.push(...newick_to_order(sub_newick, prior + parseFloat(value)))
        }
        break;

      }

    }

  } else {

    if (newick.includes(',')) {
      var parts = newick.split(',')
      if (parts[0].includes(':')){
        var value = parts[0].split(':')[1]
        var branch = parts[0].split(':')[0]
      } else {
        var value = 0 
        var branch = parts[0]
      }

      var next_iter = parts.slice(1).join(',')
      result.push(...[[branch, prior, prior + parseFloat(value)]])
      // console.log(next_iter)
      result.push(...newick_to_order(next_iter, prior))

    } else{ 
      if (newick.includes(':')) {
        var value = newick.split(':')[1]
        var branch = newick.split(':')[0]
      } else {
        var value = 0 
        var branch =  newick
      }
      result.push(...[[branch, prior, prior + parseFloat(value)]])
    }
  }

  return(result)
}

function generate_svg(data, nodes, genomes, global_x, global_y, edges, layers, layers_min, layers_max, group_dict, scale=1) {

  // access all the relevant variables from UI
  var start = new Date().getTime();
  var svg_search = [];
  var svg_backbone = [];
  var svg_text = [];
  var svg_heatmaps = [];
  var svg_edges = [];
  var svg_nodes = [];
  var svg_groups = [];
  var svg_genome_tracks = {}
  for (var genome of genomes){
    svg_genome_tracks[genome] = []
  }

  var edgecoloring = {}
  $("#genomecolors :input[type='color']").each((index, element) => {
    edgecoloring[element.id] = [index, element.value]
  })

  if ($('#flexlinear').prop('checked') == true){
    var linear = 1
  } else {
    var linear = 0
  }

  if ($('#flexsaturation').prop('checked') == true){
    var saturation = 1
  } else {
    var saturation = 0
  }

  var angle = parseFloat($('#angle')[0].value);
  var outer_margin = parseFloat($('#outer_margin')[0].value);
  var inner_margin = parseFloat($('#inner_margin')[0].value);
  var node_size = parseFloat($('#size')[0].value);
  var node_thickness = parseFloat($('#circ')[0].value);
  var edge_thickness = parseFloat($('#edge')[0].value);
  var line_thickness = parseFloat($('#line')[0].value);
  var node_distance_x = parseFloat($('#distx')[0].value);
  var node_distance_y = parseFloat($('#disty')[0].value);
  var tree_length = parseFloat($('#tree_length')[0].value);
  var offset = parseFloat($('#tree_offset')[0].value);
  var tree_thickness = parseFloat($('#tree_thickness')[0].value);

  var core_color = $('#core_color')[0].value;
  var paralog_color = $('#paralog_color')[0].value;
  var singleton_color = $('#singleton_color')[0].value;
  var accessory_color = $('#accessory_color')[0].value;
  var rearranged_color = $('#rearranged_color')[0].value;
  var trna_color = $('#trna_color')[0].value;
  var layer_color = $('#layer_color')[0].value;
  var back_color = $('#back_color')[0].value;
  var non_back_color = $('#non_back_color')[0].value;

  var theta = angle / (global_x+1)

  if (linear == 0){
    var start_offset = parseFloat($('#inner')[0].value);
  } else {
    var start_offset = 0
  }
  
  var middle_layers = new Object();
  var outer_layers = new Object();

  var search_size = parseFloat($('#search_hit')[0].value);
  middle_layers['search'] = [search_size, start_offset, search_size + start_offset]

  if ($('#flexarrow').prop('checked') == true){
    var arrow_size = parseFloat($('#arrow')[0].value)
    middle_layers['arrow'] = [arrow_size, search_size + start_offset, arrow_size + search_size + start_offset]
  } else {
    var arrow_size = 0
  }

  var graph_size = node_size * 2 + node_thickness
  outer_layers['graph'] = [graph_size, inner_margin, inner_margin + graph_size]

  var sum_middle_layer = start_offset + search_size + arrow_size
  var sum_outer_layer = graph_size + inner_margin

  var current_middle_stop = sum_middle_layer
  var current_outer_stop = sum_outer_layer

  var enabled = []
  var genome_size = genomes.length;

  var order = newick_to_order(data['meta']['newick']).reverse()

  max_dist = 0
  item_order = []
  for (var item of order) {
    var [name, item_start, item_end] = item
    if (name != 'branching') {
      item_order.push(name)
    }
    if (item_end > max_dist) {
      max_dist = item_end
    }
  }

// layer size calculations

  for (var genome of item_order) {

    if ($('#flex' + genome).prop('checked') == true){
      enabled.push(genome)
    }

    var layer_name = genome + 'layer'
    if ($('#flex' + layer_name).prop('checked') == true){
      var layer_width = parseFloat($('#' + layer_name)[0].value)

      var layer_middle_start = current_middle_stop + inner_margin
      var layer_middle_stop = layer_middle_start + layer_width

      current_middle_stop = layer_middle_stop
      sum_middle_layer += layer_width + inner_margin
      
      middle_layers[layer_name] = [layer_width, layer_middle_start, layer_middle_stop]
    }
  }

  if ($('#flexbackbone').prop('checked') == true){
    //TEST BACKBONE NON BACKBONE LAYER
    var back_width = parseFloat($('#backbone')[0].value);

    var layer_width = back_width
    var layer_middle_start = current_middle_stop + inner_margin
    var layer_middle_stop = layer_middle_start + layer_width

    current_middle_stop = layer_middle_stop
    sum_middle_layer += layer_width + inner_margin
    
    middle_layers['back_vs_non_back'] = [layer_width, layer_middle_start, layer_middle_stop]
  }

  for (var layer_name of layers) {
    if ($('#flex' + layer_name).prop('checked') == true){

      var layer_width = parseFloat($('#' + layer_name)[0].value)
      var layer_outer_start = current_outer_stop + outer_margin
      var layer_outer_stop = layer_outer_start + layer_width

      current_outer_stop = layer_outer_stop 
      sum_outer_layer += layer_width + outer_margin
      
      outer_layers[layer_name] = [layer_width, layer_outer_start, layer_outer_stop]
    }
  }

  if (linear == 0){
    var radius = 0.5 * (node_distance_x / Math.sin(deg2rad(theta * (1/2))))
    var circle_dist = sum_middle_layer + graph_size * 0.5
    var extra_offset = 0

    // if (circle_dist < radius) {
    //   extra_offset = radius - circle_dist
    // }

    sum_middle_layer += extra_offset
    for (var layer in middle_layers) {
      var [layer_width, layer_start, layer_stop] = middle_layers[layer]
      middle_layers[layer] = [layer_width, layer_start + extra_offset, layer_stop + extra_offset]
    }

    var y_size = (sum_middle_layer + (global_y * node_distance_y) + sum_outer_layer);
    var x_size = (sum_middle_layer + (global_y * node_distance_y) + sum_outer_layer);
    if (scale == 0){
      var svg_core = $('<svg id="result" width="' + x_size*2 + 'px" height="' + y_size*2 + 'px" version="1.1" viewBox="-' + x_size + ' -' + y_size + ' ' + x_size*2 + ' ' + y_size*2 + '" position="absolute" xmlns="http://www.w3.org/2000/svg"></svg>')
    } else {
      var svg_core = $('<svg id="result" width="100%" height="100%" version="1.1" viewBox="-' + x_size + ' -' + y_size + ' ' + x_size*2 + ' ' + y_size*2 + '" position="absolute" xmlns="http://www.w3.org/2000/svg"></svg>')
    }
  } else {
    var x_size = (global_x + 1) * node_distance_x * 0.5;
    var y_size = (sum_middle_layer + (global_y * node_distance_y) + sum_outer_layer) * 0.5;
    if (scale == 0){
      var svg_core = $('<svg id="result" width="' + x_size*2 + 'px" height="' + y_size*2 + 'px" version="1.1" viewBox="-' + 0.5 * node_distance_y + ' -' + y_size*2 + ' ' + x_size*2 + ' ' + y_size*2 + '" position="absolute" xmlns="http://www.w3.org/2000/svg"></svg>')
    } else {
      var svg_core = $('<svg id="result" width="100%" height="100%" version="1.1" viewBox="-' + 0.5 * node_distance_y + ' -' + y_size*2 + ' ' + x_size*2 + ' ' + y_size*2 + '" position="absolute" xmlns="http://www.w3.org/2000/svg"></svg>')
    }
  }

  for (var genome of genomes) {
    var layer_name = genome + 'layer'
    if (Object.keys(middle_layers).includes(layer_name)){

      var [layer_width, layer_start, layer_stop] = middle_layers[layer_name]
      // for (var x = 1; x < global_x+1; x++) {

      //   var add_start = 1
      //   var add_stop = 0

      //   var i_x = x-add_start-0.5
      //   var i_y = layer_start
      //   var j_x = x+add_stop+0.5
      //   var j_y = layer_stop

      //   // console.log(i_x, i_y, j_x, j_y)

      //   svg_genome_tracks[genome].push(
      //     create_rectangle(i_x, i_y, j_x, j_y, theta, node_distance_x, linear, layer_color)
      //   )
      // }

      if (linear == 0){
        var [circle_a_x, circle_a_y] = transform(0-0.5, layer_start, theta)
        var [circle_b_x, circle_b_y] = transform(0-0.5, layer_stop, theta)
        var [circle_c_x, circle_c_y] = transform(global_x + 0.5, layer_start, theta)
        var [circle_d_x, circle_d_y] = transform(global_x + 0.5, layer_stop, theta)

        if ((global_x) * theta > 180) {
          var arc_flag = 1
        } else {
          var arc_flag = 0
        }
        
        svg_genome_tracks[genome].push(
          $('<path d="M ' + circle_c_x + ' ' + circle_c_y +
          ' A ' + layer_start + ' ' + layer_start + ' 0 ' + arc_flag + ' 1 ' + circle_a_x + ' ' + circle_a_y +
          ' L ' + circle_b_x + ' ' + circle_b_y +
          ' A ' + layer_stop + ' ' + layer_stop + ' 0 ' + arc_flag + ' 0 ' + circle_d_x + ' ' + circle_d_y +
          ' Z" stroke-width="0" fill="' + layer_color + '"></path>')
        )
      } else {
        var [circle_a_x, circle_a_y] = [(0-0.5) * node_distance_x, -layer_start]
        var [circle_b_x, circle_b_y] = [(0-0.5) * node_distance_x, -layer_stop]
        var [circle_c_x, circle_c_y] = [(global_x + 0.5) * node_distance_x, -layer_start]
        var [circle_d_x, circle_d_y] = [(global_x + 0.5) * node_distance_x, -layer_stop]

        svg_genome_tracks[genome].push(
          $('<path d="M ' + circle_c_x + ' ' + circle_c_y +
          ' L ' + circle_a_x + ' ' + circle_a_y +
          ' L ' + circle_b_x + ' ' + circle_b_y +
          ' L ' + circle_d_x + ' ' + circle_d_y +
          ' Z" stroke-width="0" fill="' + layer_color + '"></path>')
        )
      }

    }
  }

  if ($('#flexarrow').prop('checked') == true){
    
    var [arrow_size, arrow_start, arrow_stop] = middle_layers['arrow']

    var pointer_length = 10/ theta
    var pointer_height = arrow_stop - arrow_start
    var arrow_thickness = pointer_height / 4
    var steps = Math.round((30 / theta))

    // console.log(middle_layers)

    if (steps < 1) {
      steps = 1
    }

    if (linear == 0){
      var [circle_c_x, circle_c_y] = transform(global_x + 0.5 - pointer_length, arrow_start + arrow_thickness, theta)
      var [circle_a_x, circle_a_y] = transform(0-0.5, arrow_start + arrow_thickness, theta)
      var [circle_b_x, circle_b_y] = transform(0-0.5, arrow_stop - arrow_thickness, theta)
      var [circle_d_x, circle_d_y] = transform(global_x + 0.5 - pointer_length, arrow_stop - arrow_thickness, theta)
      var [circle_f_x, circle_f_y] = transform(global_x + 0.5 - pointer_length, arrow_stop, theta)
      var [circle_g_x, circle_g_y] = transform(global_x + 0.5, arrow_start + arrow_thickness * 2, theta)
      var [circle_e_x, circle_e_y] = transform(global_x + 0.5 - pointer_length, arrow_start, theta)

      if ((global_x) * theta > 180) {
        var arc_flag = 1
      } else {
        var arc_flag = 0
      }

      svg_core.append(
        $('<path d="M ' + circle_c_x + ' ' + circle_c_y +
        ' A ' + (arrow_start + arrow_thickness) + ' ' + (arrow_start + arrow_thickness) + ' 0 ' + arc_flag + ' 1 ' + circle_a_x + ' ' + circle_a_y +
        ' L ' + circle_b_x + ' ' + circle_b_y +
        ' A ' + (arrow_stop - arrow_thickness) + ' ' + (arrow_stop - arrow_thickness) + ' 0 ' + arc_flag + ' 0 ' + circle_d_x + ' ' + circle_d_y +
        ' L ' + circle_f_x + ' ' + circle_f_y +
        ' L ' + circle_g_x + ' ' + circle_g_y +
        ' L ' + circle_e_x + ' ' + circle_e_y + 
        ' Z" stroke-width="0" fill="slateGrey"></path>')
      )

      var [circle_h_x, circle_h_y] = transform(0-0.5, arrow_start + arrow_thickness * 2, theta)

    } else {
      var [circle_c_x, circle_c_y] = [(global_x + 0.5 - pointer_length) * node_distance_x, -(arrow_start + arrow_thickness)]
      var [circle_a_x, circle_a_y] = [(0-0.5) * node_distance_x, -(arrow_start + arrow_thickness)]
      var [circle_b_x, circle_b_y] = [(0-0.5) * node_distance_x, -(arrow_stop - arrow_thickness)]
      var [circle_d_x, circle_d_y] = [(global_x + 0.5 - pointer_length) * node_distance_x , -(arrow_stop - arrow_thickness)]
      var [circle_f_x, circle_f_y] = [(global_x + 0.5 - pointer_length) * node_distance_x, -arrow_stop]
      var [circle_g_x, circle_g_y] = [(global_x + 0.5) * node_distance_x, -(arrow_start + arrow_thickness * 2)]
      var [circle_e_x, circle_e_y] = [(global_x + 0.5 - pointer_length) * node_distance_x, -arrow_start]

      svg_core.append(
        $('<path d="M ' + circle_c_x + ' ' + circle_c_y +
        ' L ' + circle_a_x + ' ' + circle_a_y +
        ' L ' + circle_b_x + ' ' + circle_b_y +
        ' L ' + circle_d_x + ' ' + circle_d_y +
        ' L ' + circle_f_x + ' ' + circle_f_y +
        ' L ' + circle_g_x + ' ' + circle_g_y +
        ' L ' + circle_e_x + ' ' + circle_e_y + 
        ' Z" stroke-width="0" fill="slateGrey"></path>')
      )

      var [circle_h_x, circle_h_y] = [(0-0.5) * node_distance_x, -(arrow_start + arrow_thickness * 2)]
    }
      
    svg_text.push(
      $('<text text-anchor="end" transform="translate (-10)" dominant-baseline="middle" x="' + circle_h_x + '" y="' + circle_h_y + '" dy="0" font-size="' + $('#label')[0].value + '" font-family="sans-serif" fill="black">Orientation</text>')
    )

    var l = steps
    while (l < global_x) {

      if (l+steps <= global_x){
        var k = steps;
      } else {
        var k = steps - ((l+steps) - global_x);
      }

      if (linear == 0){

        var [circle_l_x, circle_l_y] = transform(l, arrow_start + arrow_thickness * 2, theta)
        var rotate = theta * (l+0.5)
        if (rotate >= 90 && rotate <= 180) {
          rotate += 180;
        } else if (rotate >= 180 && rotate <= 270) {
          rotate -= 180;
        }
        svg_core.append(
          $('<text text-anchor="middle" dominant-baseline="middle" transform="rotate(-' + rotate + ' ' + circle_l_x + ' ' + circle_l_y +')" x="' + circle_l_x + '" y="' + circle_l_y + '" dy="0" font-size="' + $('#label')[0].value + '" font-family="sans-serif" fill="white">' + l + '</text>')
        )

      } else {
        var [circle_l_x, circle_l_y] = [(l) * node_distance_x, -(arrow_start + arrow_thickness * 2)]
        svg_core.append(
          $('<text text-anchor="middle" dominant-baseline="middle" x="' + circle_l_x + '" y="' + circle_l_y + '" dy="0" font-size="' + $('#label')[0].value + '" font-family="sans-serif" fill="white">' + l + '</text>')
        )
        
      }

      l += k
    };
  }

  var end = new Date().getTime();
  var time = end - start;
  console.log('SVG drawing core', time, 'ms.')

  var start = new Date().getTime();
  for(var i in edges) {

    var edge = data['edges'][i];
    var edge_genomes = Object.keys(edge['directions'])

    var intersection = edge_genomes.filter(x => enabled.includes(x));
    if (intersection.length > 0) {
      var edge_genomes_length = edge_genomes.length;

      var color = pickcolor (edgecoloring, Object.keys(edge['directions']))

      if (saturation == 1){
        var pick = lighter_color('#ffffff', color, edge_genomes_length / genome_size);
      } else {
        var pick = color;
      }

      var source = edge['source']
      var target = edge['target']

      if (source != 'start' && target != 'stop' && edge['active'] == true){

        var i_x = nodes[source]['position'][0]
        var i_y = nodes[source]['position'][1]
        var j_x = nodes[target]['position'][0]
        var j_y = nodes[target]['position'][1]

        for (let e = 0; e <= genomes.length; e++) {

          if (e == genomes.length || (genomes[e] + 'layer' in middle_layers && edge_genomes.includes(genomes[e])) ) {

            if (e == genomes.length) {

              var dir_set = Object.values(edge['directions'])

              if (dir_set.includes('L') && dir_set.includes('R')) {
                var stroke = ' stroke-dasharray="' + line_thickness * 4 + ' ' + line_thickness + '" '
              } else if (dir_set.includes('L')) {
                var stroke = ' stroke-dasharray="' + line_thickness + '" '
              } else {
                var stroke = ''
              }
              
              var [graph_size, graph_start, graph_stop] = outer_layers['graph']
              var i_y_size = sum_middle_layer + graph_start + graph_size * 0.5 + i_y * node_distance_y
              var j_y_size = sum_middle_layer + graph_start + graph_size * 0.5 + j_y * node_distance_y
              var draw = pick
              var thickness = edge_thickness
            } else {
              var [layer_width, layer_start, layer_stop] = middle_layers[genomes[e] + 'layer']

              if (layer_width < line_thickness) {
                var draw = ''
              } else {
                layer_width -= line_thickness
                layer_start += line_thickness * 0.5
                layer_stop -= line_thickness * 0.5

                var i_y_size = layer_start + i_y * (layer_width / global_y)
                var j_y_size = layer_start + j_y * (layer_width / global_y)
                var draw = edgecoloring[genomes[e]][1]
                var thickness = line_thickness
                var stroke = ''
              }
            }

            if (linear == 0){
              var [circle_i_x, circle_i_y] = transform(i_x, i_y_size, theta);
              var [circle_j_x, circle_j_y] = transform(j_x, j_y_size, theta);
            } else {
              var [circle_i_x, circle_i_y] = [(i_x) * node_distance_x, -i_y_size];
              var [circle_j_x, circle_j_y] = [(j_x) * node_distance_x, -j_y_size];
            }

            if (draw !== "") {

              var route_edge = '<path class="path" d="M ' + circle_i_x + ' ' + circle_i_y

              if (edge['route'].length == 0){
                if (linear == 0){
                  if (i_y == j_y) {
                    route_edge += ' A ' + i_y_size  + ' ' + j_y_size + ' 0 0 0 ' + circle_j_x + ' ' + circle_j_y + '"' + stroke + ' stroke="' + draw + '" stroke-width="' + thickness + '" fill="none"/>'
                  } else {
                    route_edge += ' L ' + circle_j_x + ' ' + circle_j_y + '"' + stroke + ' stroke="' + draw + '" stroke-width="' + thickness + '" fill="none"/>'
                  }
                } else {
                  route_edge += ' L ' + circle_j_x + ' ' + circle_j_y + '"' + stroke + ' stroke="' + draw + '" stroke-width="' + thickness + '" fill="none"/>'
                }

              } else {

                var o_y = i_y

                for(var n in edge['route']) {

                  var n_x = edge['route'][n][0]
                  var n_y = edge['route'][n][1]

                  if (e == genomes.length) {
                    var o_y_size = sum_middle_layer + graph_start + graph_size * 0.5 + o_y * node_distance_y
                    var n_y_size = sum_middle_layer + graph_start + graph_size * 0.5 + n_y * node_distance_y
                  } else {
                    var o_y_size = layer_start + o_y * (layer_width / global_y)
                    var n_y_size = layer_start + n_y * (layer_width / global_y)
                  }

                  if (linear == 0){
                    var [circle_n_x, circle_n_y] = transform(n_x, n_y_size, theta);
                  } else {
                    var [circle_n_x, circle_n_y] = [(n_x) * node_distance_x, -n_y_size];
                  }

                  if (o_y == n_y) {
                    if (linear == 0){
                      route_edge += 'A ' + o_y_size  + ' ' + n_y_size + ' 0 0 0 ' + circle_n_x + ' ' + circle_n_y
                    } else {
                      route_edge += 'L ' + circle_n_x + ' ' + circle_n_y
                    }
                  } else {
                    route_edge += 'L ' + circle_n_x + ' ' + circle_n_y
                  }
          
                  var o_y = n_y
                }

                if (o_y == j_y) {
                  if (linear == 0){
                    route_edge += 'A ' + o_y_size  + ' ' + j_y_size + ' 0 0 0 ' + circle_j_x + ' ' + circle_j_y + '"' + stroke + ' stroke="' + draw + '" stroke-width="' + thickness + '" fill="none"/>'
                  } else {
                    route_edge += 'L ' + circle_j_x + ' ' + circle_j_y + '"' + stroke + ' stroke="' + draw + '" stroke-width="' + thickness + '" fill="none"/>'
                  }
                } else {
                  route_edge += 'L ' + circle_j_x + ' ' + circle_j_y + '"' + stroke + ' stroke="' + draw + '" stroke-width="' + thickness + '" fill="none"/>'
                }
              }

              if (e == genomes.length) {
                svg_edges.push(
                  $(route_edge)
                )
              } else {
                svg_genome_tracks[genomes[e]].push(
                  $(route_edge)
                )
              }
            }
          }
        }
      }
    }
  };

  var end = new Date().getTime();
  var time = end - start;
  console.log('SVG drawing lines', time, 'ms.')

  var backbone_pos = [];

  var start = new Date().getTime();
  var global_values = []
  for(var k in nodes) {

    var node = data['nodes'][k];
    var node_genomes = Object.keys(node['gene_calls']);
    var intersection = node_genomes.filter(x => enabled.includes(x));
    if (intersection.length > 0) {

      var node_genomes_length = node_genomes.length
        
      var k_x = node['position'][0]
      var k_y = node['position'][1]
      var node_group = node['group']
      var node_type = node['type']

      if (node['layer']['backbone'] == true)
        backbone_pos.push(k_x)

      if (!node_group) {
        var node_class = 'class="node'
        var x_value_start = k_x - 0.5
        var x_value_stop = k_x + 0.5

      } else {
        var node_class = 'stroke-opacity="0" fill-opacity="0" class="pseudo'
        var group = group_dict[node_group]
        var group_size = group.length
        var group_compress = $('#groupcompress')[0].value
        var group_size_compressed = Math.round(group_size * group_compress)

        if (group_size_compressed == 0) {
          group_size_compressed = 1
        }
        
        var z_x = data['nodes'][group[0]]['position'][0]

        var fraction = group_size_compressed / (group_size)
        var group_id = group.findIndex(x => x === k)

        var x_value_start = z_x - (1 - group_id * fraction) + 0.5
        var x_value_stop = z_x - (1 - (group_id + 1) * fraction) + 0.5

      }

      var color = pickcolor (edgecoloring, Object.keys(node['gene_calls']))

      if (node_type == 'core'){
        var node_color = core_color
      } else if (node_type == 'rearrangement') {
        var node_color = rearranged_color
      } else if (node_type == 'accessory') {
        var node_color = accessory_color
      } else if (node_type == 'duplication') {
        var node_color = paralog_color
      } else if (node_type == 'singleton') {
        var node_color = singleton_color
      } else if (node_type == 'rna') {
        var node_color = trna_color
      } else {
        console.log(node_type)
      }

      if (saturation == 1) {
        var draw = lighter_color('#ffffff', color, node_genomes_length / genome_size);
        var draw2 = lighter_color('#ffffff', node_color, node_genomes_length / genome_size)
      } else {
        var draw = color;
        var draw2 = node_color
      }
      
      var [graph_size, graph_start, graph_stop] = outer_layers['graph']
      var k_y_size = sum_middle_layer + graph_start + graph_size * 0.5 + k_y * node_distance_y
      
      if (linear == 0) {
        var [circle_k_x, circle_k_y] = transform(k_x, k_y_size, theta);
      } else {
        var [circle_k_x, circle_k_y] = [(k_x) * node_distance_x, -k_y_size];
      }

      svg_nodes.push(
        $('<circle ' + node_class + '" id="' + k + '" cx="' + circle_k_x + '" cy="' + circle_k_y + '" r="' + node_size + '" fill="' + draw2 + '" stroke="' + draw + '" stroke-width="' + node_thickness + '"/>')
      )

      // var add_start = 1
      // var add_stop = 0

      var [search_size, search_start, search_stop] = middle_layers['search']
      
      var i_x = k_x - 0.5
      var i_y = search_start
      var j_x = k_x + 0.5
      var j_y = search_stop

      if (!global_values.includes(k_x)) {
        svg_search.push(create_rectangle(i_x, i_y, j_x, j_y, theta, node_distance_x, linear, 'white', k_x))
      }

      for (var layer_name of layers) {

        if ($('#flex' + layer_name).prop('checked') == true){

          var value = node['layer'][layer_name]
          var max = layers_max[layer_name]
          var min = layers_min[layer_name]

          var [layer_width, layer_start, layer_stop] = outer_layers[layer_name]
          var k_y_size = sum_middle_layer + k_y * node_distance_y

          var i_x = x_value_start
          var i_y = layer_start + k_y_size
          var j_x = x_value_stop
          var j_y = layer_stop + k_y_size
          var color = lighter_color('#00ff00', '#ff0000', (value-min) / (max-min))

          svg_heatmaps.push(create_rectangle(i_x, i_y, j_x, j_y, theta, node_distance_x, linear, color))
        }
      }
    }
  }

  global_values.push(k_x)

  if ($('#flexbackbone').prop('checked') == true){
    var k_x = 0
    while (k_x <= global_x) {
    // for(var k_x of backbone_pos) {

      var [backbone_size, backbone_start, backbone_stop] = middle_layers['back_vs_non_back']
          
      var i_x = k_x - 0.5
      var i_y = backbone_start
      var j_x = k_x + 0.5
      var j_y = backbone_stop

      if (backbone_pos.includes(k_x)) {
        var color = back_color
      } else {
        var color = non_back_color
      }

      svg_backbone.push(create_rectangle(i_x, i_y, j_x, j_y, theta, node_distance_x, linear, color, k_x))
      k_x = k_x + 1
    }
  }
  // };

  for(var [l, group] of Object.entries(group_dict)) {

    var group_length = group.length
    var group_x = group.map((a) => (data['nodes'][a]['position'][0]));
    
    ind_max = group_x.indexOf(Math.max.apply(Math, group_x))
    ind_min = group_x.indexOf(Math.min.apply(Math, group_x))
    
    var left_node_name = group[ind_min]
    var right_node_name = group[ind_max]

    var left_node = data['nodes'][left_node_name];
    var right_node = data['nodes'][right_node_name];

    var group_genomes = Object.keys(left_node['gene_calls']);
    var group_type = left_node['type']

    if (group_type == 'core'){
      var group_color = core_color
    } else if (group_type == 'rearrangement') {
      var group_color = rearranged_color
    } else if (group_type == 'accessory') {
      var group_color = accessory_color
    } else if (group_type == 'duplication') {
      var group_color = paralog_color
    } else if (group_type == 'singleton') {
      var group_color = singleton_color
    } else if (group_type == 'rna') {
      var group_color = trna_color
    }

    var intersection = group_genomes.filter(x => enabled.includes(x));
    if (intersection.length > 0) {

      var group_genomes_length = group_genomes.length;
      var color = pickcolor (edgecoloring, Object.keys(left_node['gene_calls']))
      
      if (saturation == 1) {
        var draw = lighter_color('#ffffff', color, group_genomes_length / genome_size);
      } else {
        var draw = color;
      }
      
      var l_x = left_node['position'][0]
      var l_y = left_node['position'][1]

      var m_x = right_node['position'][0]
      var m_y = right_node['position'][1]

      var [graph_size, graph_start, graph_stop] = outer_layers['graph']
      var l_y_size = sum_middle_layer + graph_start + graph_size * 0.5 + l_y * node_distance_y
      var m_y_size = sum_middle_layer + graph_start + graph_size * 0.5 + m_y * node_distance_y

      var i_x = l_x
      var i_y = l_y_size + node_size
      var j_x = m_x
      var j_y = m_y_size - node_size

      if (saturation == 1) {
        var color = lighter_color('#ffffff', group_color, group_genomes_length / genome_size)
      } else {
        var color = group_color
      }

      if ((l_x - m_x) * theta >= 180) {
        var arc_flag = 1
      } else {
        var arc_flag = 0
      }
    
      if (linear == 0) {
        var [a_x, a_y] = transform(i_x, i_y, theta)
        var [b_x, b_y] = transform(i_x, j_y, theta)
        var [c_x, c_y] = transform(j_x, i_y, theta)
        var [d_x, d_y] = transform(j_x, j_y, theta)
    
        var path = $('<path class="group" id="' + l + '" d="' +
        'M ' + a_x + ' ' + a_y + ' ' +
        'A ' + (l_y_size + node_size) + ' ' + (l_y_size + node_size) + ' 0 ' + arc_flag + ' 0 ' + c_x + ' ' + c_y + ' ' +
        'A ' + node_size + ' ' + node_size + ' 0 0 0 ' + d_x + ' ' + d_y + ' ' +
        'A ' + (m_y_size - node_size) + ' ' + (m_y_size - node_size) + ' 0 ' + arc_flag + ' 1 ' + b_x + ' ' + b_y + ' ' +
        'A ' + node_size + ' ' + node_size + ' 0 0 0 ' + a_x + ' ' + a_y +
        '" fill="' + color + '" stroke="' + draw + '" stroke-width="' + node_thickness + '"/>')
    
      } else {
        var [a_x, a_y] = [i_x * node_distance_x, -i_y]
        var [b_x, b_y] = [i_x * node_distance_x, -j_y]
        var [c_x, c_y] = [j_x * node_distance_x, -i_y]
        var [d_x, d_y] = [j_x * node_distance_x, -j_y]
    
        var path = $('<path class="group" id="' + l + '" d="' +
        'M ' + a_x + ' ' + a_y + ' ' +
        'L ' + c_x + ' ' + c_y + ' ' +
        'A ' + node_size + ' ' + node_size + ' 0 0 1 ' + d_x + ' ' + d_y + ' ' +
        'L ' + b_x + ' ' + b_y + ' ' +
        'A ' + node_size + ' ' + node_size + ' 0 0 1 ' + a_x + ' ' + a_y +
        '" fill="' + color + '" stroke="' + draw + '" stroke-width="' + node_thickness + '"/>')
      }

      svg_groups.push(path)
    }
  };

  for (var layer_name of layers) {

    if ($('#flex' + layer_name).prop('checked') == true){
      var [layer_width, layer_start, layer_stop] = outer_layers[layer_name]
      var y_size = sum_middle_layer + layer_width * 0.5
      
      if (linear == 0){
        var [circle_x, circle_y] = transform(0-0.5, (layer_start + y_size), theta)
      } else {
        var [circle_x, circle_y] = [(0-0.5) * node_distance_x, -(layer_start + y_size)]
      }

      svg_text.push(
        $('<text text-anchor="end" transform="translate (-10)" dominant-baseline="middle" x="' + circle_x + '" y="' + circle_y + '" dy="0" font-size="' + $('#label')[0].value + '" font-family="sans-serif" fill="black">' + layer_name + '</text>')
      )
    }
  }

  item_dist = {}
  for (var genome_name of item_order) {
    if ($('#flex' + genome_name + 'layer').prop('checked') == true){
      var [layer_width, layer_start, layer_stop] = middle_layers[genome_name + 'layer']
      
      if (layer_width >= edge_thickness) {
      
        var y_size = layer_start + layer_width * 0.5

        if (linear == 0){
          var [circle_x, circle_y] = transform(0-0.5, y_size, theta)
        } else {
          var [circle_x, circle_y] = [(0-0.5) * node_distance_x, -y_size]
        }

        item_dist[genome_name] = circle_y

        svg_text.push(
          $('<text text-anchor="end" transform="translate (-10)" dominant-baseline="middle" x="' + circle_x + '" y="' + circle_y + '" dy="0" font-size="' + $('#label')[0].value + '" font-family="sans-serif" fill="black">' + genome_name + '</text>')
        )
      }
    }
  }

  // tree_length
  // console.log(Object.keys(item_dist).length, item_order.length)

  if ($('#flextree').prop('checked') == true){
    svg_core.append(draw_newick(order, item_dist, max_dist, offset, tree_length, tree_thickness))
  }

  var end = new Date().getTime();
  var time = end - start;
  console.log('SVG drawing remaining elements', time, 'ms.')

  var svg_group = $('<g></g>')
  for (var item of svg_search) svg_group.append(item);
  svg_core.append(svg_group);

  var svg_group = $('<g></g>')
  for (var item of svg_heatmaps) svg_group.append(item);
  svg_core.append(svg_group);

  for (var [genome, svg_genome_track] of Object.entries(svg_genome_tracks)) {
    var svg_group = $('<g></g>')
    for (var item of svg_genome_track) svg_group.append(item);
    svg_core.append(svg_group);
  }

  var svg_group = $('<g></g>')
  for (var item of svg_backbone) svg_group.append(item);
  svg_core.append(svg_group);

  var svg_group = $('<g></g>')
  for (var item of svg_edges) svg_group.append(item);
  svg_core.append(svg_group);

  var svg_group = $('<g></g>')
  for (var item of svg_nodes) svg_group.append(item);
  svg_core.append(svg_group);

  var svg_group = $('<g></g>')
  for (var item of svg_groups) svg_group.append(item);
  svg_core.append(svg_group);

  var svg_group = $('<g></g>')
  for (var item of svg_text) svg_group.append(item);
  svg_core.append(svg_group);

  return svg_core
}

function defineVariables(data) {
  var all_nodes = data['nodes'];
  var all_edges = data['edges'];

  // var layers = new Set();
  var layers = data['meta']['layers']
  layers = layers.filter(item => item !== 'backbone')

  var layers_min = new Object();
  var layers_max = new Object();
  // var groups = new Set();
  var global_x = 0;
  var global_y = 0;
  var genomes = data['meta']['genome_names']
  var group_dict = {}

  for(var n in all_nodes) {
    var node = all_nodes[n];
    for(var [layer, value] of Object.entries(node["layer"])) {
      if (layer in layers_min) {
        layers_min[layer] = value < layers_min[layer] ? value : layers_min[layer]
      } else {
        layers_min[layer] = value
      }

      if (layer in layers_max) {
        layers_max[layer] = value > layers_max[layer] ? value : layers_max[layer]
      } else {
        layers_max[layer] = value
      }

      // layers.add(layer)
    }

    var group = node["group"]
    var x = node["position"][0]
    var y = node["position"][1]
    
    if (group) {
      if (group in group_dict) {
        group_dict[group].push(n)
      } else {
        group_dict[group] = [n]
      } 
    }

    // groups.add(group)
    global_x = x < global_x ? global_x : x
    global_y = y < global_y ? global_y : y

  }

  for(var e in all_edges) {
    var edge = all_edges[e];
    var route = edge['route']
    if (route.length > 0 && edge['active'] == true) {
      for (var b in route) {
        var x = route[b][0]
        var y = route[b][1]
        global_y = y < global_y ? global_y : y
      }
    }
  }

  return [all_nodes, all_edges, layers, layers_min, layers_max, global_x, global_y, genomes, group_dict]
}

//ANCHOR - MAIN
$(document).ready(function() {

  var mapAS = {
    'A': '<span style="color: #000000;">A</span>', 'R': '<span style="color: #ff0000;">R</span>',
    'N': '<span style="color: #000000;">N</span>', 'D': '<span style="color: #000000;">D</span>',
    'C': '<span style="color: #000000;">C</span>', 'Q': '<span style="color: #000000;">Q</span>',
    'E': '<span style="color: #000000;">E</span>', 'G': '<span style="color: #ffa500;">G</span>',
    'H': '<span style="color: #ff0000;">H</span>', 'I': '<span style="color: #00ff00;">I</span>',
    'L': '<span style="color: #00ff00;">L</span>', 'K': '<span style="color: #ff0000;">K</span>',
    'M': '<span style="color: #00ff00;">M</span>', 'F': '<span style="color: #00ff00;">F</span>',
    'P': '<span style="color: #ffa500;">P</span>', 'S': '<span style="color: #ffa500;">S</span>',
    'T': '<span style="color: #ffa500;">T</span>', 'W': '<span style="color: #00ffff;">W</span>',
    'Y': '<span style="color: #00ffff;">Y</span>', 'V': '<span style="color: #00ff00;">V</span>',
    '-': '<span style="color: #000000;">-</span>'
  };
  var state = 'default'
  var functional_annotation_sources_available = [];
  var bins = {"bin1": []};
  var binnum = 1;
  var old_data ={}
  var new_data = {}
  var current_groups = {}

  $.ajax({
    url: "/pangraph/get_json",
    type: "POST",
    cache: false,
    contentType: "application/json",
    dataType: "json",
    error: function(){
      console.log('Error while attempting to load JSON data.')
    },
    success: function(data){
      console.log('Successfully load JSON data.')

      console.log(data)

      //ANCHOR - UI FUNCTIONS

      var start = new Date().getTime();
      var [all_nodes, all_edges, layers, layers_min, layers_max, global_x, global_y, genomes, group_dict] = defineVariables(data)

      for (var layer of layers) {
    
        if ($('#flex' + layer + '').length == 0) {
          var element = $('<div class="col-12 d-flex mb-1"></div>').append(
            $('<div class="col-2 d-flex align-items-center"></div>').append(
              $('<div class="form-switch d-flex"></div>').append(
                $('<input class="" type="checkbox" id="flex' + layer + '" name="flex' + layer + '" aria-label="..." data-toggle="tooltip" data-placement="top" title="Tooltip">')
              )
            )
          ).append(
              $('<div class="col-8 d-flex align-items-center"></div>').append(
                layer
              )
          ).append(
            $('<div class="d-flex col-2"></div>').append(
              $('<input type="text" class="form-control float-end text-end flex-fill p-0 border-0" style= "background-color: #e9ecef;" id="' + layer + '" name="' + layer + '" value=0 aria-label="..." data-toggle="tooltip" data-placement="top" title="Choose your color">')
            )
          )

          $('#layers').append(element)  
        }
      }

      $('#title-panel-first-line').text(data['meta']['project_name']);
      $('#title-panel-second-line').text('Pangraph Detail');

      if (!$('#genomecolors').children().length) {

        for (var genome of genomes) {

          $('#genomecolors').append(
            $('<div class="col-12 d-flex mb-1">').append(
              $('<div class="col-2 d-flex align-items-center">').append(
                $('<div class="form-switch d-flex">').append(
                  $('<input class="" type="checkbox" id="flex' + genome + '" name="' + genome + '" aria-label="..." data-bs-toggle="tooltip" data-bs-placement="top" title="Tooltip on top">')
                )
              )
            ).append(
              $('<div class="col-7 d-flex align-items-center">').append(
                genome
              )
            ).append(
              $('<div class="col-1 d-flex align-items-center">').append(
                $('<i class="user-handle bi bi-arrows-expand"></i>')
              )
            ).append(
              $('<div class="d-flex col-2">').append(
                $('<input type="color" class="form-control form-control-color flex-fill p-0 border-0" id="' + genome + '" name="' + genome + '" value="#000000" aria-label="..." data-bs-toggle="tooltip" data-bs-placement="top" title="Choose your color">')
              )
            )
          )

          $('#RightOffcanvasBodyTop').append(
            $('<tr>').append(
              $('<td class="col-8">').append(
                genome
              )
            ).append(
              $('<td class="col-4 text-end" id="number_' + genome + '">').append(
                0
              )
            )
          )

          if ($('#flex' + genome + 'layer').length == 0) {
            var element = $('<div class="col-12 d-flex mb-1"></div>').append(
              $('<div class="col-2 d-flex align-items-center"></div>').append(
                $('<div class="form-switch d-flex"></div>').append(
                  $('<input class="" type="checkbox" id="flex' + genome + 'layer" name="flex' + genome + 'layer" aria-label="..." data-toggle="tooltip" data-placement="top" title="Tooltip on top">')
                )
              )
            ).append(
              $('<div class="col-8 d-flex align-items-center"></div>').append(
                genome
              )
            ).append(
              $('<div class="d-flex col-2"></div>').append(
                $('<input type="text" class="form-control float-end text-end flex-fill p-0 border-0" style= "background-color: #e9ecef;" id="' + genome + 'layer" name="' + genome + 'layer" value=0 aria-label="..." data-toggle="tooltip" data-placement="top" title="Choose your color">')
              )
            )

            $('#layers').append(element)
          }
        }
      }

      for (var [setting, value] of Object.entries(data['states'])) {
        if (typeof value === 'number') {
          $('#' + setting)[0].value = value 
        } else if (value == true || value == false) {
          $('#' + setting).prop('checked', value);
        } else {
          $('#' + setting)[0].value = value 
        }
      }

      functional_annotation_sources_available = data['meta']['gene_function_sources'];
      $('#searchSources').empty()
      for (var annotation_source of functional_annotation_sources_available){
        $('#searchSources').append(
          $('<div class="col-12"></div>').append(
            $('<div class="row align-items-center"></div>').append(
              $('<div class="col-2 mb-1"></div>').append(
                $('<input class="" type="checkbox" id="flex' + annotation_source + '" value="" aria-label="..." data-toggle="tooltip" data-placement="top" title="Tooltip on top">')
              )
            ).append(  
              $('<div class="col-8 mb-1"></div>').append(
                annotation_source
              )
            ).append(
              $('<div class="col-2 mb-1"></div>')
            )
          )
        )
      }

      $('#expressiondrop').empty()
      $('#expressiondrop').append($('<option value="Choose item">Choose item</option>'))
      $('#expressiondrop').append($('<option value="Name">Name</option>'))
      $('#expressiondrop').append($('<option value="Position">Position</option>'))

      $('#filter').empty()
      $('#filter').append(
        $('<div class="col-12"></div>').append(
          $('<div class="row align-items-center"></div>').append(
            $('<div class="col-2 mb-1"></div>').append(
              $('<input class="" type="checkbox" id="minposition" value="" aria-label="..." data-toggle="tooltip" data-placement="top" title="Tooltip on top">')
            )
          ).append(  
            $('<div class="col-8 mb-1"></div>').append(
              'Min graph position'
            )
          ).append(
            $('<div class="col-2 mb-1"></div>').append(
              $('<input type="text" class="form-control flex-fill p-0 border-0" style= "background-color: #e9ecef;" value="" id="minpositiontext" aria-describedby="">')
            )  
          ).append(
            $('<div class="col-2 mb-1"></div>').append(
              $('<input class="" type="checkbox" id="maxposition" value="" aria-label="..." data-toggle="tooltip" data-placement="top" title="Tooltip on top">')
            )
          ).append( 
            $('<div class="col-8 mb-1"></div>').append(
              'Max graph position' 
            )
          ).append(
            $('<div class="col-2 mb-1"></div>').append(
              $('<input type="text" class="form-control flex-fill p-0 border-0" style= "background-color: #e9ecef;" value="" id="maxpositiontext" aria-describedby="">')
            )
          )
        )
      )

      for (var layer of layers) {
        if ($('#flex' + layer).prop('checked') == true){

          $('#expressiondrop').append($('<option value="' + layer + '">' + layer + '</option>'))
          $('#filter').append(
            $('<div class="col-12"></div>').append(
              $('<div class="row align-items-center"></div>').append(
                $('<div class="col-2 mb-1"></div>').append(
                  $('<input class="" type="checkbox" id="min' + layer + '" value="" aria-label="..." data-toggle="tooltip" data-placement="top" title="Tooltip on top">')
                )
              ).append(  
                $('<div class="col-8 mb-1"></div>').append(
                  'Min ' + layer
                )
              ).append(
                $('<div class="col-2 mb-1"></div>').append(
                  $('<input type="text" class="form-control flex-fill p-0 border-0" style= "background-color: #e9ecef;" value="" id="min' + layer + 'text" aria-describedby="">')
                )  
              ).append(
                $('<div class="col-2 mb-1"></div>').append(
                  $('<input class="" type="checkbox" id="max' + layer + '" value="" aria-label="..." data-toggle="tooltip" data-placement="top" title="Tooltip on top">')
                )
              ).append( 
                $('<div class="col-8 mb-1"></div>').append(
                  'Max ' + layer
                )
              ).append(
                $('<div class="col-2 mb-1"></div>').append(
                  $('<input type="text" class="form-control flex-fill p-0 border-0" style= "background-color: #e9ecef;" value="" id="max' + layer + 'text" aria-describedby="">')
                )
              )
            )
          )
        }
      }

      if ($('#flexlinear').prop('checked') == false) {
        $('#distx').prop('disabled', true)
        $('#inner').prop('disabled', false)
      } else {
        $('#distx').prop('disabled', false)
        $('#inner').prop('disabled', true)
      }

      if ($('#flexarrow').prop('checked') == false) {
        $('#arrow').prop('disabled', true)
      } else {
        $('#arrow').prop('disabled', false)
      }

      if ($('#flexbackbone').prop('checked') == false) {
        $('#backbone').prop('disabled', true)
      } else {
        $('#backbone').prop('disabled', false)
      }

      if ($('#flexcondtr').prop('checked') == false) {
        $('#condtr').prop('disabled', true)
      } else {
        $('#condtr').prop('disabled', false)
      }

      if ($('#flexmaxlength').prop('checked') == false) {
        $('#maxlength').prop('disabled', true)
      } else {
        $('#maxlength').prop('disabled', false)
      }

      if ($('#flexgroupcompress').prop('checked') == false) {
        $('#groupcompress').prop('disabled', true)
      } else {
        $('#groupcompress').prop('disabled', false)
      }

      $('#flexlinear').change(function() {
        if ($(this).prop('checked') == true){
          $('#distx').prop('disabled', false)
          $('#inner').prop('disabled', true)
        } else {
          $('#distx').prop('disabled', true)
          $('#inner').prop('disabled', false)
        }
      })

      $('#flexbackbone').change(function() {
        if ($(this).prop('checked') == true){
          $('#backbone')[0].value = 100;
          $('#backbone').prop('disabled', false)
        } else {
          $('#backbone')[0].value = 0;
          $('#backbone').prop('disabled', true)
        }
      })

      $('#flexarrow').change(function() {
        if ($(this).prop('checked') == true){
          $('#arrow')[0].value = 100;
          $('#arrow').prop('disabled', false)
        } else {
          $('#arrow')[0].value = 0;
          $('#arrow').prop('disabled', true)
        }
      })
  
      $('#flexcondtr').change(function() {
        if ($(this).prop('checked') == true){
          $('#condtr')[0].value = 2;
          $('#condtr').prop('disabled', false)
        } else {
          $('#condtr')[0].value = -1;
          $('#condtr').prop('disabled', true)
        }
      })
  
      $('#flexmaxlength').change(function() {
        if ($(this).prop('checked') == true){
          $('#maxlength')[0].value = 1;
          $('#maxlength').prop('disabled', false)
        } else {
          $('#maxlength')[0].value = -1;
          $('#maxlength').prop('disabled', true)
        }
      })
  
      $('#flexgroupcompress').change(function() {
        if ($(this).prop('checked') == true){
          $('#groupcompress')[0].value = 0.0;
          $('#groupcompress').prop('disabled', false)
        } else {
          $('#groupcompress')[0].value = 1.0;
          $('#groupcompress').prop('disabled', true)
        }
      })
  
      $('#flextree').change(function() {
        if ($(this).prop('checked') == true){
          for (var genome of genomes) {
            if ($('#flex' + genome + 'layer').prop('checked') == false){
              $('#' + genome + 'layer')[0].value = 50
              $('#flex' + genome + 'layer').prop('checked', true);
            }
            $('#flex' + genome + 'layer').prop('disabled', true);
          }
        } else {
          for (var genome of genomes) {
            $('#flex' + genome + 'layer').prop('disabled', false);
          }
        }
      })

      // ARE YOU SERIOUS JAVASCRIPT?
      // We have to get create a variable based on the functions
      // entities ID called "entry" because all other variables
      // would change dynamically forcing always a change in 
      // the function pointing to the last iteration.
      for (var layer_name of layers) {
        $('#flex' + layer_name).change(function() {
          var entry = $(this)[0].id.replace("flex", "")
          if ($(this).prop('checked') == true){
            $('#' + entry)[0].value = 25;
            $('#' + entry).prop('disabled', false)
          } else {
            $('#' + entry)[0].value = 0;
            $('#' + entry).prop('disabled', true)
          }
        })
      }

      for (var genome of genomes) {
        $('#flex' + genome + 'layer').change(function() {
          var entry = $(this)[0].id.replace("flex", "")
          if ($(this).prop('checked') == true){
            $('#' + entry)[0].value = 50;
            $('#' + entry).prop('disabled', false)
          } else {
            $('#' + entry)[0].value = 0;
            $('#' + entry).prop('disabled', true)
          }
        })
      }

      sortable('#genomecolors', {
        forcePlaceholderSize: true,
        handle: '.user-handle',
        items: 'div'
      });

      $('#InfoDownload').on('click', function() {

        // Variable to store the final csv data
        var csv_data = [];
        var basics = $('#node_basics_table')
        var title = basics[0].getAttribute("gc_id")
        var layers = $('#node_layers_table')
        var functions = $('#node_functions_table')

        var basics_rows = basics[0].getElementsByTagName('tr');
        var layers_rows = layers[0].getElementsByTagName('tr');
        
        var function_rows = functions[0].getElementsByTagName('tr');

        for (let i = 0; i < function_rows.length; i++) {
      
            if (i >= basics_rows.length) {
              var basics_cols = []
              basics_cols = Array.prototype.concat.apply(basics_cols, basics_rows[1].querySelectorAll('td,th'));
              basics_cols = Array.prototype.concat.apply(basics_cols, layers_rows[1].querySelectorAll('td,th'));
              
              var function_cols = function_rows[i].querySelectorAll('td,th');
            } else { 
              var basics_cols = []
              basics_cols = Array.prototype.concat.apply(basics_cols, basics_rows[i].querySelectorAll('td,th'));
              basics_cols = Array.prototype.concat.apply(basics_cols, layers_rows[i].querySelectorAll('td,th'));
              
              var function_cols = function_rows[i].querySelectorAll('td,th');
            }
              
            let csvrow = [];
            for (let j = 0; j < basics_cols.length; j++) {
      
              var info = basics_cols[j].innerHTML
              csvrow.push(info);
            }
            
            for (let k = 0; k < function_cols.length; k++) {
      
              var info = function_cols[k].innerHTML
              csvrow.push(info);
            }
      
            // console.log(csvrow)

            // Combine each column value with comma
            csv_data.push(csvrow.join(","));
        }
        // Combine each row data with new line character
        csv_data = csv_data.join('\n');
      
        var blob = new Blob([csv_data]);
        // var title = data['infos']['meta']['title']
        downloadBlob(blob, title + ".csv");
      });

      $('#AlignmentDownload').on('click', function() {

        // Variable to store the final csv data
        var csv_data = '';
        var alignment = $('#node_sequence_alignments_table')
        var basics = $('#node_basics_table')
        var title = alignment[0].getAttribute("gc_id")
        var xpos = basics[0].getAttribute("gc_pos")
        
        var alignment_rows = alignment[0].getElementsByTagName('tr');

        for (let i = 1; i < alignment_rows.length; i++) {
      
            var alignment_cols = alignment_rows[i].querySelectorAll('td,th');
      
            csv_data += ">" + title + "|Genome:" + alignment_cols[0].innerHTML +"|Genecall:" + alignment_cols[1].innerHTML + "|Position:" + xpos + '\n';
            var genome = ''
            
            var alignment_nucs = alignment_cols[2].getElementsByTagName('span');
            // console.log(alignment_nucs)

            for (let k = 0; k < alignment_nucs.length; k++) {

              genome += alignment_nucs[k].innerHTML
            }

            csv_data += genome.match(/.{1,60}/g).join("\r\n") + "\n"
        }
        // Combine each row data with new line character
        
        var blob = new Blob([csv_data]);
        // var title = data['infos']['meta']['title']
        downloadBlob(blob, title + ".fa");
      });

      old_data['condtr'] = data['states']['condtr']
      old_data['maxlength'] = data['states']['maxlength']
      old_data['groupcompress'] = data['states']['groupcompress']
      // old_data['ungroupfrom'] = data['states']['ungroupfrom']
      // old_data['ungroupto'] = data['states']['ungroupto']
      old_data['state'] = data['meta']['state']

      var end = new Date().getTime();
      var time = end - start;
      console.log('Document creation', time, 'ms.')

    }
  })  

  // ANCHOR - DRAW/REDRAW
  $('#redraw').on('click', function() {

    var state = 'default'
    new_data['condtr'] = parseInt($('#condtr')[0].value)
    new_data['maxlength'] = parseInt($('#maxlength')[0].value)
    new_data['groupcompress'] = parseFloat($('#groupcompress')[0].value)
    // new_data['ungroupfrom'] = $('#ungroupfrom')[0].value
    // new_data['ungroupto'] = $('#ungroupto')[0].value
    new_data['state'] = state

    var reiterate = false

    if (new_data['condtr'] != old_data['condtr']) {
      old_data['condtr'] = parseInt($('#condtr')[0].value)
      reiterate = true
    }

    if (new_data['maxlength'] != old_data['maxlength']) {
      old_data['maxlength'] = parseInt($('#maxlength')[0].value)
      reiterate = true
    }

    if (new_data['groupcompress'] != old_data['groupcompress']) {
      old_data['groupcompress'] = parseFloat($('#groupcompress')[0].value)
      reiterate = true
    }

    // if (new_data['ungroupfrom'] != old_data['ungroupfrom']) {
    //   old_data['ungroupfrom'] = parseInt($('#ungroupfrom')[0].value)
    //   reiterate = true
    // }

    // if (new_data['ungroupto'] != old_data['ungroupto']) {
    //   old_data['ungroupto'] = parseInt($('#ungroupto')[0].value)
    //   reiterate = true
    // }

    if (new_data['state'] != old_data['state']) {
      old_data['state'] = state
      reiterate = true
    }

    if (reiterate == true) {

      $.ajax({
        url: "/pangraph/settings",
        type: "POST",
        async: false,
        data: JSON.stringify(new_data),
        contentType: "application/json",
        dataType: "json",
        error: function(){
          console.log('Error while attempting to update JSON data.')
        },
        success: function(){
          console.log('Successfully updated JSON data.')
        }
      });
    }

    $.ajax({
      url: "/pangraph/get_json",
      type: "POST",
      cache: false,
      contentType: "application/json",
      dataType: "json",
      error: function(){
        console.log('Error while attempting to load JSON data.')
      },
      success: function(data){
        console.log('Successfully load JSON data.')

        var body = $('#svgbox')
        body.empty()

        var [all_nodes, all_edges, layers, layers_min, layers_max, global_x, global_y, genomes, group_dict] = defineVariables(data)
        var groups = Object.keys(group_dict)
        var svg_core = generate_svg(data, all_nodes, genomes, global_x, global_y, all_edges, layers, layers_min, layers_max, group_dict);
        var genome_size = genomes.length

        body.append(svg_core)
        body.html(body.html());

        if (typeof window.zoomSVG !== 'undefined') {
          window.zoomSVG.destroy()
        }
        window.zoomSVG = svgPanZoom('#result', {
          zoomEnabled: true,
          panEnabled: false,
          controlIconsEnabled: false,
          minZoom: 0.1,
          maxZoom: 100
        });

        $('#fit').off('click')
        $('#fit').on('click', function() {
          window.zoomSVG.resize();
          window.zoomSVG.fit();
          window.zoomSVG.center();
        })

        $('#svgDownload').off('click')
        $('#svgDownload').on('click', function() {
          var svg_download = generate_svg(data, all_nodes, genomes, global_x, global_y, all_edges, layers, layers_min, layers_max, group_dict, scale=0)
          var blob = new Blob([svg_download[0].outerHTML]);
          var title = data['meta']['project_name']
          downloadBlob(blob, title + ".svg");
        });

        var gc_nodes = document.querySelectorAll(".node")
        var divs = document.querySelectorAll(".node, .group");
        for (var el of divs) {

          if (el.getAttribute("class") == 'group'){
            var id = group_dict[el.getAttribute("id")][0]
            var name = el.getAttribute("id")
          } else {
            var id = el.getAttribute("id")
            var name = data['nodes'][el.getAttribute("id")]["gene_cluster"]
          }

          tippy(el, {
            content: '<strong>' + name + '</strong>' + '<br />',
            allowHTML: true,
            onHide(instance) {
              if (instance.reference.id.startsWith('GCG_')){
                var id = group_dict[instance.reference.id][0]
              } else {
                var id = instance.reference.id
              }
              var elements = Object.keys(data['nodes'][id]['gene_calls'])

              for (var element of elements) {
                $('#number_' + element)[0].innerText = '0';
              }
            },
            onShow(instance) {
              if (instance.reference.id.startsWith('GCG_')){
                var id = group_dict[instance.reference.id][0]
              } else {
                var id = instance.reference.id
              }
              var elements = Object.keys(data['nodes'][id]['gene_calls'])

              for (var element of elements) {
                $('#number_' + element)[0].innerText = '1';
              }
            },
            arrow: false,
            duration: 0,
            followCursor: true,
            theme: "light",
          });
        };

        var isDown = false
        var diff = 0

        var old_xpos = 0
        var old_ypos = 0

        var xpos = 0
        var ypos = 0

        var new_xpos = 0
        var new_ypos = 0

        body.off('mousedown')
        body.on('mousedown', function(e) {
          old_xpos = e.offsetX
          old_ypos = e.offsetY

          xpos = old_xpos
          ypos = old_ypos

          isDown = true
          diff = 0
        })

        body.off('mousemove')
        body.on('mousemove', function(e) {
          if (isDown === true) {
            new_xpos = e.offsetX
            new_ypos = e.offsetY

            diff += Math.sqrt((new_xpos-xpos)^2+(new_ypos-ypos)^2)

            if (!e.shiftKey) {
              window.zoomSVG.panBy({x: new_xpos-xpos, y: new_ypos-ypos})
            }

            xpos = new_xpos
            ypos = new_ypos
          }
        })

        body.off('mouseup')
        body.on('mouseup', function(e) {
          if (isDown === true) {

            var selection = document.querySelector('input[name="binradio"]:checked')

            isDown = false

            if (diff < 10) {
              if (e.target.getAttribute('class') === 'group' || e.target.getAttribute('class') === 'node') {

                if (e.shiftKey && selection !== null) {

                  var binid = selection.value
                  bins = marknode(e.target, data, binid, bins, genome_size, group_dict);

                } else {
                  //show_node_details_modal(e.target, data);
                  // console.log(e.target)
                  nodeinfo(e.target, data, group_dict, mapAS);
                }

              } else {
              }

            } else {
              if (e.shiftKey && selection !== null) {

                var binid = selection.value

                var max_xpos = Math.max(old_xpos, xpos)
                var min_xpos = Math.min(old_xpos, xpos)

                var max_ypos = Math.max(old_ypos, ypos)
                var min_ypos = Math.min(old_ypos, ypos)

                for (var n of gc_nodes) {

                  var bounding = n.getBoundingClientRect();
                  var left = bounding.left
                  var right = bounding.right
                  var bottom = bounding.bottom
                  var top = bounding.top

                  if (
                    min_xpos < left &&
                    max_xpos > right &&
                    min_ypos < bottom &&
                    max_ypos > top
                    ) {
                    bins = marknode(n, data, binid, bins, genome_size, group_dict);

                  }
                }

                for(var g of groups) {
                  // var inside = true;
                  var group = group_dict[g]
                  for (var k of group) {
                    var node = document.getElementById(k);

                    var bounding = node.getBoundingClientRect();
                    var left = bounding.left
                    var right = bounding.right
                    var bottom = bounding.bottom
                    var top = bounding.top

                    if (
                      min_xpos < left &&
                      max_xpos > right &&
                      min_ypos < bottom &&
                      max_ypos > top
                      ) {
                      var name = document.getElementById(g);
                      bins = marknode(name, data, binid, bins, genome_size, group_dict);
                      break
                    }
                    // if (
                    //   min_xpos < left &&
                    //   max_xpos > right &&
                    //   min_ypos < bottom &&
                    //   max_ypos > top
                    //   ) {
                    // } else {
                    //   inside = false
                    // }
                  }

                  // if (inside === true){
                  //   var name = document.getElementById(g);
                  //   bins = marknode(name, data, binid, bins, genome_size, group_dict);
                  // }
                }
              }
            }
          }
        })


        for (var binid of Object.keys(bins)) {
          var nodes = bins[binid];
          var updated_nodes = []
        
          for (var node of nodes) {
        
            var name = document.getElementById(node);
        
            if(name) {
              if (name.getAttribute('class') == 'pseudo') {
                for(var g in groups) {
                  var group = group_dict[g]
                  if (group.includes(node)) {
                    if (!updated_nodes.includes(g)){
                      updated_nodes.push(g)
                    }
                  }
                }
              } else {
                updated_nodes.push(node)
              }
            } else {
              updated_nodes.push(...group_dict[node])
            };
        
          }
        
          bins[binid] = updated_nodes
          // console.log(nodes)
          // console.log(updated_nodes)
          for (var node of bins[binid]) {
        
            bins[binid] = bins[binid].filter(item => item !== node);
            var name = document.getElementById(node);
            bins = marknode(name, data, binid, bins, genome_size, group_dict);
        
          }
        }


        for (var binid of Object.keys(bins)) {
      
          $('#' + binid + 'color').off("change")
          $('#' + binid + 'color').on("change", function() {

            var binid = this.name;
            var nodes = bins[binid];
      
            for (var node of nodes) {
      
              bins[binid] = bins[binid].filter(item => item !== node);
              var name = document.getElementById(node);
              bins = marknode(name, data, binid, bins, genome_size, group_dict);
      
            }
          })
        }

        $('#AddBin').on('click')
        $('#AddBin').on('click', function() {
        
          var selection = document.querySelector('input[name="binradio"]:checked');
          var binid = selection.value;
        
          var basics = $('#node_basics_table')
          var node = basics[0].getAttribute("gc_context")
          var name = document.getElementById(node);
        
          bins = marknode(name, data, binid, bins, genome_size, group_dict);
        
        });

        $('#binadd').off('click')
        $('#binadd').on('click', function() {
          binnum += 1;
          var binid = "bin" + binnum
        
          $('#bingrid').append(
              $('<div class="col-12" id="' + binid + '"></div>').append(
                  $('<div class="row gy-1 align-items-center" id="row' + binnum + '"></div>').append(
                      $('<div class="col-2"></div>').append(
                          $('<input type="radio" name="binradio" id="bin' + binnum + 'radio" value="bin' + binnum + '" checked></input>')
                      )
                  ).append(
                      $('<div class="col-6"></div>').append(
                          $('<input class="form-control flex-fill p-0 border-0" style="background-color: #e9ecef;" type="text" id="bin' + binnum + 'text" value="Bin_' + binnum + '"></input>')
                      )
                  ).append(
                      $('<div class="col-2"></div>').append(
                          $('<input type="text" class="form-control float-end text-end flex-fill p-0 border-0" id="bin' + binnum + 'value" name="condtr" value="0" readonly>')
                      )
                  ).append(
                      $('<div class="d-flex col-2"></div>').append(
                          $('<input class="form-control form-control-color flex-fill p-0 border-0 colorchange" type="color" name="bin' + binnum + '" id="bin' + binnum + 'color" value="#000000"></input>')
                      )
                  )
              )
          );
        
          $('#' + binid + 'color').on("change", function() {

            var binid = this.name;
            var nodes = bins[binid];
      
            for (var node of nodes) {
      
              bins[binid] = bins[binid].filter(item => item !== node);
              var name = document.getElementById(node);
              bins = marknode(name, data, binid, bins, genome_size, group_dict);
      
            }
          })

          bins['bin' + binnum] = [];
        });
        
        $('#binremove').off('click')
        $('#binremove').on('click', function() {
        
          var selection = document.querySelector('input[name="binradio"]:checked');
        
          if (selection !== null) {
            var binid = selection.value;
        
            for (var node of bins[binid]) {
              var name = document.getElementById(node);
              bins = marknode(name, data, binid, bins, genome_size, group_dict);
            }
        
            $("#" + binid).remove();
            delete bins[binid]
            var nextbin = document.querySelector('input[name="binradio"]')
        
            if (nextbin) {
              nextbin.checked = true;
            } else {
              $('#binadd').click();
            }
          }
        })

        $('#binanalyse').off('click')
        $('#binanalyse').on('click', async function() {
          var set = {}
          var bin_keys = Object.keys(bins)

          var selection = 'COG24_FUNCTION'
          set['selection'] = selection

          for (var binid of bin_keys) {
            set[binid] = []
            for (var id of bins[binid]) {
              if (Object.keys(group_dict).includes(id)) {
                for (var item of group_dict[id]) {
                  set[binid].push(data['nodes'][item]['gene_cluster'])  
                }
              } else {
                set[binid].push(data['nodes'][id]['gene_cluster'])
              }
            }
          }

          $.ajax({
            url: "/pangraph/analyse",
            type: "POST",
            data: JSON.stringify(set),
            contentType: "application/json",
            dataType: "json",
            success: function(result){
              
              $('#BinAnalyseBody').empty();
              var bin_analyses_table = ''

              bin_analyses_table += '<table class="table'
              bin_analyses_table += '<tbody><thead class="thead-light"><tr>'
              
              bin_analyses_table += '<th scope="row">' + 'Function' + '</th>'
              for (bin_key of bin_keys) {
                bin_analyses_table += '<th scope="row">' + bin_key + '</th>'
              }
              
              bin_analyses_table += '</tr></thead><tbody>'
              bin_analyses_table += '<tbody>'

              for (func of Object.keys(result)) {
                bin_analyses_table += '<tr><td>' + func + '</td>'
                for (bin_key of Object.keys(result[func])) {
                  var value = result[func][bin_key]
                  if (value == 1) {
                    bin_analyses_table += '<td class=td-analysis-cell-1></td>'
                  } else {
                    bin_analyses_table += '<td class=td-analysis-cell-0></td>'
                  }
                }
                bin_analyses_table += '</tr>'
              }

              bin_analyses_table += '</tbody></table>'

              $('#BinAnalyseBody')[0].innerHTML = bin_analyses_table

              $('#BinAnalyse').modal('show');

            }
          })
        })

        $('#bininfo').off('click')
        $('#bininfo').on('click', async function() {

          var selection = document.querySelector('input[name="binradio"]:checked');

          // if (selection !== null) {
          var binid = selection.value;
          var appendlist = [];

          $('#BinModalBody').empty();
          for (var id of bins[binid]) {
            var element = document.getElementById(id);

            if (element.getAttribute('class') == 'group') {
              appendlist.push(...group_dict[id])
            } else {
              appendlist.push(id);
            }
          }

          for (var gene_cluster_id of appendlist) {

            gene_cluster_context = null

            var body = $('<div></div>').append(
              $('<div class="card-body overflow-scroll" id="' + gene_cluster_id + 'div"></div>').append(
                await get_gene_cluster_display_tables(gene_cluster_id, gene_cluster_context, data, group_dict, mapAS, add_align=0)
              )
            )
            $('#BinModalBody').append(
              $('<div class="card border mb-3 w-100"></div>').append(
                body
              )
            )
          }
          $('#BinModal').modal('show');
        })


        // function get_passed_gene_clusters(searchfunction) {
  
        //   var d = $.ajax({
        //     url: "/pangraph/filter",
        //     type: "POST",
        //     data: JSON.stringify(searchfunction),
        //     contentType: "application/json",
        //     dataType: "json"
        //   })
        
        //   return d
        // }
        
        // async function checknode(searchfunction) {
        
        //     var keys = Object.keys(searchfunction)  
        //     if (keys.length > 0) {
        //         var d = await get_passed_gene_clusters(searchfunction)
        //         var result = d['gene_clusters']
        //     } else {
        //         var result = []
        //     }
        
        //     return result
        // }
        
        // $('#searchadd').off('click')
        // $('#searchadd').on('click', function() {
        
        //     var selection = document.querySelector('input[name="binradio"]:checked')
        //     var binid = selection.value
        
        //     for (var [id, members] of Object.entries(searched)) {
        //     if (!(bins[binid].includes(id))) {
        
        //         var e = document.getElementById(id);
        //         bins = marknode(e, data, binid, bins);
        
        //     }
        //     }
        // })
        
        // $('#searchremove').off('click')
        // $('#searchremove').on('click', function() {
        
        //     var selection = document.querySelector('input[name="binradio"]:checked')
        //     var binid = selection.value
        
        //     for (var [id, members] of Object.entries(searched)) {
        //     if ((bins[binid].includes(id))) {
        
        //         var e = document.getElementById(id);
        //         bins = marknode(e, data, binid, bins);
        
        //     }
        //     }
        // })
        
        // $('#searcherase').off('click')
        // $('#searcherase').on('click', function() {
        
        //     for (var [id, members] of Object.entries(searched)) {
        //     for (var mem of members) {
        //         var xpos = data['elements']['nodes'][mem]['position']['x_offset']
        //         var e = document.getElementById(xpos);
        //         e.setAttribute("fill", "white")
        //     }
        //     }
        
        //     searched = {};
        // })
        
        // $('#searchcolor').off('click')
        // $('#searchcolor').on('click', function() {
        
        //     for (var [id, members] of Object.entries(searched)) {
        //     for (var mem of members) {
        //         var xpos = data['elements']['nodes'][mem]['position']['x_offset']
        //         var e = document.getElementById(xpos);
        //         e.setAttribute("fill", "#ff0000")
        //     }
        //     }
        // })
        
        // const searchEraseButton = document.getElementById('searcherase');
        // const searchColorButton = document.getElementById('searchcolor');
        // const searchAppendBin = document.getElementById('searchadd');
        // const searchRemoveBin = document.getElementById('searchremove');
        
        // var searched = {}
        // $('#search').off('click')
        // $('#search').on('click', async function() {
        
        //     var searchpos = false
        
        //     //Filter Search
        //     var layers_filter = new Object
        //     var mingenomes = ($("#mingenomes").prop('checked') == true && !isNaN(parseFloat($("#mingenomestext")[0].value))) ? parseFloat($("#mingenomestext")[0].value) : -1
        //     var maxgenomes = ($("#maxgenomes").prop('checked') == true && !isNaN(parseFloat($("#maxgenomestext")[0].value))) ? parseFloat($("#maxgenomestext")[0].value) : -1
        //     var minposition = ($("#minposition").prop('checked') == true && !isNaN(parseFloat($("#minpositiontext")[0].value))) ? parseFloat($("#minpositiontext")[0].value) : -1
        //     var maxposition = ($("#maxposition").prop('checked') == true && !isNaN(parseFloat($("#maxpositiontext")[0].value))) ? parseFloat($("#maxpositiontext")[0].value) : -1
            
        //     layers_filter['genomes'] = {'min': mingenomes, 'max': maxgenomes}
        //     layers_filter['position'] = {'min': minposition, 'max': maxposition}
        
        //     var layers = Object.keys(data['infos']['layers_data'])
        //     for (var layer_name of layers) {
        
        //     // console.log("#min" + layer_name + "text")
        //     var minlayer = ($("#min" + layer_name).prop('checked') == true && !isNaN(parseFloat($("#min" + layer_name + "text")[0].value))) ? parseFloat($("#min" + layer_name + "text")[0].value) : -1
        //     var maxlayer = ($("#max" + layer_name).prop('checked') == true && !isNaN(parseFloat($("#max" + layer_name + "text")[0].value))) ? parseFloat($("#max" + layer_name + "text")[0].value) : -1
            
        //     layers_filter[layer_name] = {'min': minlayer, 'max': maxlayer}
        //     }
        
        //     //Expression Search
        //     var expressioncomparison = ''
        //     var expressiondrop = $('#expressiondrop')[0].value
        //     var expressionrel = $('#expressionrel')[0].value
        //     var expressiontext = $('#expressiontext')[0].value
        
        //     // console.log(expressiondrop, expressionrel, expressiontext)
        
        //     if (expressiondrop != "Choose item" && expressionrel != "Choose operator" && expressiontext != '') {
        //     if (expressionrel == '=') {
        //         if (!isNaN(expressiontext)) {
        //         expressioncomparison = '== ' + expressiontext
        //         } else {
        //         expressioncomparison = '== "' + expressiontext + '"'
        //         }
        //     } else if (expressionrel == '\u{2260}') {
        //         if (!isNaN(expressiontext)) {
        //         expressioncomparison = '!= ' + expressiontext
        //         } else {
        //         expressioncomparison = '!= "' + expressiontext + '"'
        //         }
        //     } else if (expressionrel == '\u{2264}' && !isNaN(expressiontext)) {
        //         expressioncomparison = '<= ' + expressiontext
        //     } else if (expressionrel == '\u{2265}' && !isNaN(expressiontext)) {
        //         expressioncomparison = '>= ' + expressiontext
        //     } else if (expressionrel == '\u{003C}' && !isNaN(expressiontext)) {
        //         expressioncomparison = '< ' + expressiontext
        //     } else if (expressionrel == '\u{003E}' && !isNaN(expressiontext)) {
        //         expressioncomparison = '> ' + expressiontext
        //     } else if (expressionrel == '\u{25C2}\u{25AA}\u{25B8}') {
        //         expressioncomparison = '.includes("' + expressiontext + '")'
        //     } else if (expressionrel == '\u{25C2}\u{25AA}') {
        //         expressioncomparison = '.endsWith("' + expressiontext + '")'
        //     } else if (expressionrel == '\u{25AA}\u{25B8}') {
        //         expressioncomparison = '.startsWith("' + expressiontext + '")'
        //     }
        //     } 
        
        //     // console.log(expressioncomparison)
        
        //     //Function Search
        //     var searchfunction = {}
        //     var searchterms = $('#searchFunctionsValue')[0].value.split(",");
        //     for (var source of functional_annotation_sources_available) {
        //     if ($("#flex" + source).prop('checked') == true) {
        //         searchfunction[source] = searchterms
        //     }
        //     }
        
        //     console.log(searchfunction)
        
        //     var layers_positions = {}
        //     for (var layer_name of Object.keys(layers_filter)) {
            
        //     var minlayer = layers_filter[layer_name]['min']
        //     var maxlayer = layers_filter[layer_name]['max']
        
        //     if (minlayer != -1 || maxlayer != -1 || (minlayer != -1 && maxlayer != -1)) {
        
        //         layers_positions[layer_name] = []
        //         searchpos = true
        //         if (layer_name == 'position') {
        //         var global_x = data["infos"]["meta"]["global_x_offset"] -1;
        //         var layerobjects = new Array(global_x - 1).fill().map((d, i) => i + 1);
        //         // console.log(layerobjects)
        //         } else {
        //         var layerobjects = document.querySelectorAll("." + layer_name)
        //         }
        
        //         for (var o of layerobjects) {
        
        //         var append = true
                
        //         if (layer_name == 'position') {
        //             var value = o
        //             var result = o
        //         } else {
        //             var value = o.getAttribute("name")
        //             var result = o.getAttribute("xpos")
        //         }
        
        //         if (minlayer != '-1'){
        //             if (eval(value + '<' + minlayer)){
        //             append = false
        //             }
        //         }
        
        //         if (maxlayer != '-1'){
        //             if (eval(value + '>' + maxlayer)){
        //             append = false
        //             }
        //         }
        
        //         if (append == true) {
        //             layers_positions[layer_name].push(parseInt(result))
        //         }
        
        //         }
        //     }
        //     }
        
        //     // console.log(layers_positions)
        
        //     var positions = []
        //     var keys = Object.keys(layers_positions)
        
        //     if (keys.length > 0) {
        //     for (var pos of layers_positions[keys[0]]) {
        
        //         if (keys.length > 0) {
        //         var add = true
        //         for (let k = 1; k < keys.length; k++) {
        //             if (!layers_positions[keys[k]].includes(pos)) {
        //             add = false
        //             }
        //         }
        
        //         if (add == true) {
        //             positions.push(pos)
        //         }
        //         } else {
        //         positions.push(pos)
        //         }
        //     }
        //     }
            
        //     // console.log(positions)
        
        //     if (expressioncomparison == '' && Object.keys(searchfunction).length == 0 && searchpos == false) {
        //     var message = "Please enter valid search parameters."
        //     } else {
        //     // console.log(expressioncomparison, searchfunction, searchpos)
        
        //     var passed_gcs = await checknode(searchfunction)
            
        //     var nodes = document.querySelectorAll(".node")
        //     for (var node of nodes) {
        
        //         var id = node.getAttribute("id")
        //         var node = data['elements']['nodes'][id]
                
        //         if (passed_gcs.includes(node['name'])) {
        //         var append = true
        //         } else {
        //         var append = false
        //         }
        
        //         if (searchpos == true) {
        //         if (!positions.includes(parseInt(node['position']['x_offset']))) {
        //             append = false
        //         }
        //         }
        
        //         if (append == true) {
        //         if (expressiondrop == "Name") {
        //             if (eval('"' + node["name"] + '"' + expressioncomparison)) {
        //             if (!(id in searched)) {
        //                 searched[id] = [id]
        //             }
        //             }
        //         } else {
        //             if (!(id in searched)) {
        //             searched[id] = [id]
        //             }
        //         }
        //         }
        //     }
        
        //     var groups = document.querySelectorAll(".group")
        //     for (var group of groups) {
        //         var groupid = group.getAttribute("id")
        //         var members = data["infos"]["groups"][groupid]
        //         for (var id of members) {
        
        //         node = data['elements']['nodes'][id]
        //         var append = true
        
        //         if (passed_gcs.includes(node['name'])) {
        //             var append = true
        //         } else {
        //             var append = false
        //         }              
        
        //         if (searchpos == true) {
        //             if (!positions.includes(parseInt(node['position']['x_offset']))) {
        //             append = false
        //             }
        //         }
        
        //         if (append == true) {
        //             if (expressiondrop == "Name") {
        //             if (eval('"' + node["name"] + '"' + expressioncomparison) || eval('"' + group + '"' + expressioncomparison)) {
        
        //                 if (!(groupid in searched)) {
        //                 searched[groupid] = [id]
        //                 } else {
        //                 searched[groupid].push(id)
        //                 }
        //             }
        //             } else {
        //             if (!(groupid in searched)) {
        //                 searched[groupid] = [id]
        //             } else {
        //                 searched[groupid].push(id)
        //             }
        //             }
        //         }
        //         }
        //     }
            
        //     var table = $('<div class="form-horizontal"></div>')
        
        //     for (var key of Object.keys(searched)) {
                
        //         table.append(
        //         $('<div class="col-12"></div>').append(
        //             $('<div class="row align-items-center"></div>').append(
        //             $('<div class="col-12 mb-1">' + searched[key] + '</td>')
        //             )
        //         )
        //         )
        //     }
        
        //     $('#searchtable').empty()
        //     $('#searchtable').append(
        //         $('<div class="pt-3"></div>').append(
        //         $('<div class="shadow-box pb-3 pr-3 pl-3 mb-3 rounded search-box-filter"></div>').append(
        //             table
        //         )
        //         )
        //     )
        
        //     var message = 'You have ' + Object.keys(searched).length + ' item(s) in your queue.'
        
        //     if (Object.keys(searched).length != 0) {
        //         searchEraseButton.disabled = false;
        //         searchColorButton.disabled = false;
        //         searchAppendBin.disabled = false;
        //         searchRemoveBin.disabled = false;
        //     }
        //     }
        
        //     var toastbody = $('#searchtoastbody')
        //     toastbody.empty()
        //     toastbody.append(
        //     message
        //     )
        //     // var searchtoast = bootstrap.Toast.getOrCreateInstance($('#searchtoast'))
        //     $('#searchtoast').toast('show')
        // })
        
        

      }
    })
  })

});