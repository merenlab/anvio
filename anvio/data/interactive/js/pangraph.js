//ANCHOR - Constants
const mapAS = {
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

// we will fill this later from the incoming data
var functional_annotation_sources_available = [];
//ANCHOR - Bin creation
var bins = {"bin1": []};
var binnum = 1;

var current_groups = {}

// NOTE - From https://stackoverflow.com/questions/1053843/get-the-element-with-the-highest-occurrence-in-an-array
function modeString(array) {
  if (array.length == 0) return null;

  var modeMap = {},
    maxEl = array[0],
    maxCount = 1;

  for (var i = 0; i < array.length; i++) {
    var el = array[i];

    if (modeMap[el] == null) modeMap[el] = 1;
    else modeMap[el]++;

    if (modeMap[el] > maxCount) {
      maxEl = el;
      maxCount = modeMap[el];
    } else if (modeMap[el] == maxCount && maxCount == 1) {
      maxEl = el;
      maxCount = modeMap[el];
    } else if (modeMap[el] == maxCount) {
      maxEl += ' & ' + el;
      maxCount = modeMap[el];
    }
  }
  return [maxEl, Math.round(maxCount/i * 100)];
}

//ANCHOR - Fetch GC consensus functions
function get_gene_cluster_consensus_functions(gene_cluster_data) {
  var d = new Object();
  if (functional_annotation_sources_available.length > 0) {
    for (var source of functional_annotation_sources_available) {

      var id = []
      var func = []

      if (gene_cluster_data != '') {
        for (var element of Object.keys(gene_cluster_data)) {
          var entry = gene_cluster_data[element][source]

          if (entry === 'None' || entry === undefined) {
            entry = ['-', '-', '-']
          }

          id.push(entry[0].split('!!!')[0]);
          func.push(entry[1].split('!!!')[0]);
        }

        var [id_maxEl, id_maxCount] = modeString(id);
        var [func_maxEl, func_maxCount] = modeString(func);
      } else {

        var id_maxEl = ''
        var id_maxCount = ''
        var func_maxEl = ''
        var func_maxCount = ''
      }
      d[source] = [id_maxEl, id_maxCount, func_maxEl, func_maxCount]
    }
  }

  return(d);
}


function get_gene_cluster_basics_table(gene_cluster_id, data, add_align) {
    // first, learn a few basics about the gene cluster to be displayed
    var position_in_graph = data['elements']['nodes'][gene_cluster_id]['position']['x'] + " / " + (data["infos"]["meta"]["global_x"] - 1);
    var num_contributing_genomes = Object.keys(data['elements']['nodes'][gene_cluster_id]['genome']).length + " / " + (data['infos']['num_genomes']);
    var gene_cluster_name = data['elements']['nodes'][gene_cluster_id]['name']

    basic_info = {'Gene Cluster': gene_cluster_name, 'Contributing Genomes': num_contributing_genomes, 'Position in Graph': position_in_graph}
    // build the basic information table
    if (add_align == 1) {
      basic_info_table = `<p class="modal_header mt-0">Basics</p>`;
      basic_info_table += `<table class="table table-striped table-bordered sortable" gc_id="` + gene_cluster_id + `" id="node_basics_table">`;
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


function get_gene_cluter_functions_table(gene_cluster_id, data, add_alig) {
    
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
    functions_table += `<th scope="col">Agreement</th>`;
    functions_table += `</tr></thead><tbody>\n\n`;

    var gene_cluster_data = data['elements']['nodes'][gene_cluster_id]['genome']
    var d = get_gene_cluster_consensus_functions(gene_cluster_data);

    function_sources = Object.keys(d).sort();
    // console.log(function_sources);
    for (index in function_sources) {
        source = function_sources[index];
        value = d[source];

        var [id_maxEl, id_maxCount, func_maxEl, func_maxCount] = value

        functions_table += `<tr>`;
        functions_table += `<td>` + source + `</td>`;
        functions_table += `<td>` + id_maxEl + `</td>`;
        functions_table += `<td>` + func_maxEl + `</td>`;
        functions_table += `<td>` + func_maxCount + `%</td>`;
        functions_table += `</tr>`;
    }
    functions_table += `</tbody></tr></table>\n\n`;

    return functions_table
}


function get_gene_cluster_context_table(gene_cluster_id_current, gene_cluster_context, data, add_align) {
    
    if (gene_cluster_context == null){
        return '';
    } else {
      var group_context = []
      for (item in data['infos']['groups'][gene_cluster_context]){
        group_context.push(data['infos']['groups'][gene_cluster_context][item])
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
        gene_cluster_name = data['elements']['nodes'][gene_cluster_id]['name']

        if (gene_cluster_id == gene_cluster_id_current){
            gene_cluster_context_table += `<span class="gene_cluster_id gene_cluster_id_current">` + gene_cluster_name + `</span>`;
        } else {
            gene_cluster_context_table += `<span class="gene_cluster_id"><a class="btn border-0 m-0 p-0 align-baseline group_choice" context="` + add_align + `" group="` + gene_cluster_context + `" name_id="` + gene_cluster_id + `">` + gene_cluster_name + `</a></span>`;
        }
    }
    gene_cluster_context_table += `</div>`;

    return gene_cluster_context_table;

}


//ANCHOR - Get tables for GC basic info and functions
async function get_gene_cluster_display_tables(gene_cluster_id, gene_cluster_context, data, add_align) {
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

    gene_cluster_context_table = get_gene_cluster_context_table(gene_cluster_id, gene_cluster_context, data, add_align);

    ///////////////////////////////////////////////////////////////////////////////////////////
    // BUILD BASIC INFORMATION TABLE
    ///////////////////////////////////////////////////////////////////////////////////////////

    basic_info_table = get_gene_cluster_basics_table(gene_cluster_id, data, add_align);

    ///////////////////////////////////////////////////////////////////////////////////////////
    // BUILD FUNCTIONS TABLE
    ///////////////////////////////////////////////////////////////////////////////////////////

    functions_table = get_gene_cluter_functions_table(gene_cluster_id, data, add_align);

    ///////////////////////////////////////////////////////////////////////////////////////////
    // RETRIEVE AND BUILD SEQUENCE ALIGNMENTS
    ///////////////////////////////////////////////////////////////////////////////////////////
 
    if (add_align == 1) {

      var alignment = {}

      // if (gene_cluster_id != 'start' && gene_cluster_id != 'stop') {
      for (var genome of Object.keys(data['elements']['nodes'][gene_cluster_id]['genome'])) {
        var genecall = data['elements']['nodes'][gene_cluster_id]['genome'][genome]['gene_call']
        var contig = data['elements']['nodes'][gene_cluster_id]['genome'][genome]['contig']
        var direction = data['elements']['nodes'][gene_cluster_id]['genome'][genome]['direction']
        var name = data['elements']['nodes'][gene_cluster_id]['name']
        alignment[genome] = [genecall, contig, direction, name]
      }
      // }

      gene_cluster_sequence_alignments_table = await appendalignment(gene_cluster_id, alignment)

      ///////////////////////////////////////////////////////////////////////////////////////////
      // MERGE ALL AND RETURN
      ///////////////////////////////////////////////////////////////////////////////////////////

      return gene_cluster_context_table + basic_info_table + functions_table + gene_cluster_sequence_alignments_table

    } else {

      return gene_cluster_context_table + basic_info_table + functions_table

    }
}

//ANCHOR - Fetch GC alignment
function fetchalignment(alignment) {

  var d = $.ajax({
    url: "/pangraph/alignment",
    type: "POST",
    data: JSON.stringify(alignment),
    contentType: "application/json",
    dataType: "json"
  })

  return d

};

//ANCHOR - Append GC alignment
async function appendalignment(gene_cluster_id, alignment) {
    // FIXME: We need to come up with more accurate and explicit function names.
    if (Object.keys(alignment).length == 0) {
        return ''
    }

    alignments_table = `<p class="modal_header">Sequence alignments</p>`;
    alignments_table += `<div class="scroll-wrapper"><table class="table sortable" gc_id="` + gene_cluster_id + `" id="node_sequence_alignments_table">`;
    alignments_table += `<thead class="thead-dark"><tr>`;
    alignments_table += `<th class="position-sticky" style="left:0px; z-index:2;" scope="col">Genome</th>`;
    alignments_table += `<th scope="col">Gene Call</th>`;
    alignments_table += `<th scope="col">Contig</th>`;
    alignments_table += `<th scope="col">Direction</th>`;
    alignments_table += `<th scope="col"><span id="th-sequence">Sequence</span></th>`;
    alignments_table += `</tr></thead><tbody>\n\n`;


    var d = await fetchalignment(alignment)

    // console.log(d)

    for (var [genome, value] of Object.entries(d)) {
      var colored = value[3].replace(/A|R|N|D|C|Q|E|G|H|I|L|K|M|F|P|S|T|W|Y|V|-/gi, function(matched){return mapAS[matched];});

      alignments_table += `<tr>`
      alignments_table += `<td id="td-genome-cell">` + genome + `</td>`
      alignments_table += `<td id="td-value-cell">` + value[0] + `</a></td>`
      alignments_table += `<td id="td-contig-cell">` + value[1] + `</td>`
      alignments_table += `<td id="td-direction-cell">` + value[2] + `</a></td>`
      alignments_table += `<td id="gc-alignment-font"><div class="scrollable-content">` + colored + `</div></td>`
      alignments_table += `</tr>`
    }

    alignments_table += `</tbody></table></div>\n\n`

    // I guess this shouldn't return anything but insert the results into a context?
    // lol sorry if I screwed these thigns up badly :) IDK what I'm doing :p
    return alignments_table;
}

//ANCHOR - Color node and add/remove from bin
function marknode(e, data, binid, bins){

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

  if (e.getAttribute('class') == 'group') {

    var node_color = $('#groups')[0].value;
    var genome_size = Object.keys(data["infos"]["genomes"]).length;

    var group = data['infos']['groups'][id]
    var node_name = group[0]

    var node = data['elements']['nodes'][node_name];
    var genome = Object.keys(node['genome']).length;


  } else if (e.getAttribute('class') == 'node') {

    var node_color = $('#nodes')[0].value;
    var genome_size = Object.keys(data["infos"]["genomes"]).length;

    var node = data['elements']['nodes'][id];
    var genome = Object.keys(node['genome']).length;

  }

  if (current === binid) {

    e.setAttribute("fill", lighter_color('#ffffff', node_color, genome / genome_size))
    bins[binid] = bins[binid].filter(item => item !== id)
    $('#' + binid + 'value')[0].value = bins[binid].length

  } else if (current === '') {

    e.setAttribute("fill", lighter_color('#ffffff', bincolor, genome / genome_size))
    bins[binid].push(id)
    $('#' + binid + 'value')[0].value = bins[binid].length

  } else {

    e.setAttribute("fill", lighter_color('#ffffff', bincolor, genome / genome_size))
    bins[current] = bins[current].filter(item => item !== id)
    bins[binid].push(id)
    $('#' + binid + 'value')[0].value = bins[binid].length
    $('#' + current + 'value')[0].value = bins[current].length

  }

  return bins
}

// //ANCHOR - Information for the GC
async function nodeinfo(e, data) {
  var id = e.id;
  var element = document.getElementById(id);

  if (element.getAttribute('class') == 'group') {
    gene_cluster_context = id;
    gene_cluster_id = data['infos']['groups'][id][0]
  } else {
    gene_cluster_id = id;
    gene_cluster_context = null;
  }

  $('#InfoModalBody').empty()
  var bodyinfo = $('<div class="card-body overflow-scroll"></div>')
  // var closeBtn = $('<div class="row justify-content-end"><i class="btn d-flex col-1 justify-content-end text-danger fs-2x bi bi-x-circle-fill" data-dismiss="modal"></i></div>')
  // $('#InfoModalBody').append(closeBtn)
  $('#InfoModalBody').append(bodyinfo)

  // console.log(gene_cluster_id, gene_cluster_context)

  var all_info = await get_gene_cluster_display_tables(gene_cluster_id, gene_cluster_context, data, add_align=1)

  // console.log(all_info)

  bodyinfo.append(all_info)

  $('#InfoModal').modal('show');
}

//ANCHOR - Degree to rad calculation
function deg2rad(degrees)
{
  return degrees * Math.PI/180;
}

//ANCHOR - General download function
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

//ANCHOR - Function to mix two colors
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

//ANCHOR - Transformation function
function transform(x, y, theta) {
  var circle_x = y * Math.sin(deg2rad(theta * x))
  var circle_y = y * Math.cos(deg2rad(theta * x))
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

//ANCHOR - All in one SVG creation function
//ANCHOR - All in one SVG creation function
async function generate_svg(body, data) {

  var svg_search = [];
  var svg_heatmaps = [];
  var svg_edges = [];
  var svg_nodes = [];
  var svg_groups = [];

  var edgecoloring = {}
  $("#genomecolors :input[type='color']").each((index, element) => {
    edgecoloring[element.id] = [index, element.value]
  })

  var outer_margin = parseInt($('#outer_margin')[0].value);
  var inner_margin = parseInt($('#inner_margin')[0].value);
  var node_size = parseInt($('#size')[0].value);
  var node_thickness = parseInt($('#circ')[0].value);
  var edge_thickness = parseInt($('#edge')[0].value);
  var line_thickness = parseInt($('#line')[0].value);
  // var node_distance_x = parseInt($('#distx')[0].value);
  // var node_distance_y = parseInt($('#disty')[0].value);

  var groups = data['infos']['groups']
  var edges = data['elements']['edges']
  var nodes = data['elements']['nodes']

  var global_x = data["infos"]["meta"]["global_x"] -1;
  var global_y = data["infos"]["meta"]["global_y"] +1;

  var group_color = $('#groups')[0].value;
  var node_color = $('#nodes')[0].value;

  var theta = 270 / (global_x)
  // var radius = 0.5 * (node_distance_x / Math.sin(deg2rad(theta * (1/2))))  
  var start_offset = parseInt($('#inner')[0].value)
  
  var middle_layers = new Object();
  var outer_layers = new Object();

  var search_size = 200
  middle_layers['search'] = [search_size, start_offset, search_size + start_offset]

  var arrow_size = parseInt($('#arrow')[0].value)
  middle_layers['arrow'] = [arrow_size, search_size + start_offset, arrow_size + search_size + start_offset]

  var graph_size = node_size * 2 + node_thickness
  outer_layers['graph'] = [graph_size, inner_margin, inner_margin + graph_size]

  var sum_middle_layer = start_offset + search_size + arrow_size
  var sum_outer_layer = graph_size + inner_margin

  var current_middle_stop = sum_middle_layer
  var current_outer_stop = sum_outer_layer

  var genomes = Object.keys(data["infos"]["genomes"])
  var genome_size = genomes.length;

  for (var genome of genomes) {

    var layer_name = genome + 'layer'
    if ($('#flex' + layer_name).prop('checked') == true){
      var layer_width = parseInt($('#' + layer_name)[0].value)

      var layer_middle_start = current_middle_stop + inner_margin
      var layer_middle_stop = layer_middle_start + layer_width

      current_middle_stop = layer_middle_stop
      sum_middle_layer += layer_width + inner_margin
      
      middle_layers[layer_name] = [layer_width, layer_middle_start, layer_middle_stop]
    }
  }

  var layers = data['infos']['layers_names']
  for (var layer_name of layers) {
    if ($('#flex' + layer_name).prop('checked') == true){

      var layer_width = parseInt($('#' + layer_name)[0].value)
      var layer_scale = data['infos']['layers_data'][layer_name]['scale']

      if (layer_scale == 'global') {
        var layer_middle_start = current_middle_stop + inner_margin
        var layer_middle_stop = layer_middle_start + layer_width

        current_middle_stop = layer_middle_stop
        sum_middle_layer += layer_width + inner_margin
        
        middle_layers[layer_name] = [layer_width, layer_middle_start, layer_middle_stop]
      } else {
        var layer_outer_start = current_outer_stop + outer_margin
        var layer_outer_stop = layer_outer_start + layer_width

        current_outer_stop = layer_outer_stop 
        sum_outer_layer += layer_width + outer_margin
        
        outer_layers[layer_name] = [layer_width, layer_outer_start, layer_outer_stop]
      }
    }
  }

  var size = (sum_middle_layer + global_y * sum_outer_layer);
  var svg_core = $('<svg id="result" width="100%" height="100%" version="1.1" viewBox="-' + size + ' -' + size + ' ' + size*2 + ' ' + size*2 + '" position="absolute" xmlns="http://www.w3.org/2000/svg">')

  if ($('#flexarrow').prop('checked') == true){
    
    var [arrow_size, arrow_start, arrow_stop] = middle_layers['arrow']

    var pointer_length = 5/ theta
    var pointer_height = arrow_stop - arrow_start
    var arrow_thickness = pointer_height / 4
    var steps = Math.round((30 / theta))

    if (steps < 1) {
      steps = 1
    }

    var [circle_c_x, circle_c_y] = transform(global_x - pointer_length, arrow_start + arrow_thickness, theta)
    var [circle_a_x, circle_a_y] = transform(0, arrow_start + arrow_thickness, theta)
    var [circle_b_x, circle_b_y] = transform(0, arrow_stop - arrow_thickness, theta)
    var [circle_d_x, circle_d_y] = transform(global_x - pointer_length, arrow_stop - arrow_thickness, theta)
    var [circle_f_x, circle_f_y] = transform(global_x - pointer_length, arrow_stop, theta)
    var [circle_g_x, circle_g_y] = transform(global_x, arrow_start + arrow_thickness * 2, theta)
    var [circle_e_x, circle_e_y] = transform(global_x - pointer_length, arrow_start, theta)

    if ((global_x) * theta >= 180) {
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

    var [circle_h_x, circle_h_y] = transform(0, arrow_start + arrow_thickness * 2, theta)
    svg_core.append(
      $('<text text-anchor="end" transform="translate (-10)" dominant-baseline="middle" x="' + circle_h_x + '" y="' + circle_h_y + '" dy="0" font-size="' + $('#label')[0].value + '" font-family="sans-serif" fill="black">Orientation</text>')
    )

    var l = 1
    while (l < global_x) {

      if (l+steps <= global_x){
        var k = steps;
      } else {
        var k = steps - ((l+steps) - global_x);
      }
      var [circle_l_x, circle_l_y] = transform(l-0.5, arrow_start + arrow_thickness * 2, theta)
      var rotate = theta * (l-0.5)
      if (rotate >= 90 && rotate <= 180) {
        rotate += 180;
      } else if (rotate >= 180 && rotate <= 270) {
        rotate -= 180;
      }
      svg_core.append(
        $('<text text-anchor="middle" dominant-baseline="middle" transform="rotate(-' + rotate + ' ' + circle_l_x + ' ' + circle_l_y +')" x="' + circle_l_x + '" y="' + circle_l_y + '" dy="0" font-size="' + $('#label')[0].value + '" font-family="sans-serif" fill="white">' + l + '</text>')
      )
      l += k
    };
  }

  for(var i in edges) {

    var edge = data['elements']['edges'][i];
    var edge_genomes = Object.keys(edge['genome'])
    var edge_genomes_length = edge_genomes.length;

    var color = pickcolor (edgecoloring, Object.keys(edge['genome']))
    var pick = lighter_color('#ffffff', color, edge_genomes_length / genome_size);

    var source = edge['source']
    var target = edge['target']

    if (source != 'start' && target != 'stop' && edge['shown'] == 1){

      var i_x = nodes[source]['position']['x']
      var i_y = nodes[source]['position']['y']
      var j_x = nodes[target]['position']['x']
      var j_y = nodes[target]['position']['y']

      for (let e = 0; e <= genomes.length; e++) {

        if (e == genomes.length || (genomes[e] + 'layer' in middle_layers && edge_genomes.includes(genomes[e])) ) {

          if (e == genomes.length) {

            if (edge['direction'] == 'L') {
              var stroke = ' stroke-dasharray="5,5" '
            } else if (edge['direction'] == 'B') {
              var stroke = ' stroke-dasharray="15,5" '
            } else {
              var stroke = ''
            }

            var [graph_size, graph_start, graph_stop] = outer_layers['graph']
            var i_y_size = sum_middle_layer + graph_start + graph_size * 0.5 + i_y * sum_outer_layer
            var j_y_size = sum_middle_layer + graph_start + graph_size * 0.5 + j_y * sum_outer_layer
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

              var i_y_size = layer_start + (i_y + 0.5) * (layer_width / global_y)
              var j_y_size = layer_start + (j_y + 0.5) * (layer_width / global_y)
              var draw = edgecoloring[genomes[e]][1]
              var thickness = line_thickness
              var stroke = ''
            }
          }

          var [circle_i_x, circle_i_y] = transform(i_x-0.5, i_y_size, theta);
          var [circle_j_x, circle_j_y] = transform(j_x-0.5, j_y_size, theta);

          if (draw !== "") {

            if (edge['bended'] == ""){

              if (i_y == j_y) {
                svg_edges.push(
                  $('<path class="path" d="M ' + circle_i_x + ' ' + circle_i_y + ' A ' + i_y_size  + ' ' + j_y_size + ' 0 0 0 ' + circle_j_x + ' ' + circle_j_y + '"' + stroke + ' stroke="' + draw + '" stroke-width="' + thickness + '" fill="none"/>')
                )  
              } else {
                svg_edges.push(
                  $('<path class="path" d="M ' + circle_i_x + ' ' + circle_i_y + ' L ' + circle_j_x + ' ' + circle_j_y + '"' + stroke + ' stroke="' + draw + '" stroke-width="' + thickness + '" fill="none"/>')
                )
              }

            } else {

              var bended_edge = '<path class="path" d="M ' + circle_i_x + ' ' + circle_i_y
              var o_y = i_y

              for(var n in edge['bended']) {

                var n_x = edge['bended'][n]['x']
                var n_y = edge['bended'][n]['y']

                if (e == genomes.length) {
                  var o_y_size = sum_middle_layer + graph_start + graph_size * 0.5 + o_y * sum_outer_layer
                  var n_y_size = sum_middle_layer + graph_start + graph_size * 0.5 + n_y * sum_outer_layer
                } else {
                  var o_y_size = layer_start + o_y * (layer_width / global_y)
                  var n_y_size = layer_start + n_y * (layer_width / global_y)
                }

                var [circle_n_x, circle_n_y] = transform(n_x-0.5, n_y_size, theta);

                if (o_y == n_y) {
                  bended_edge += 'A ' + o_y_size  + ' ' + n_y_size + ' 0 0 0 ' + circle_n_x + ' ' + circle_n_y
                } else {
                  bended_edge += 'L ' + circle_n_x + ' ' + circle_n_y
                }
        
                var o_y = n_y
              }

              if (o_y == j_y) {
                bended_edge += 'A ' + o_y_size  + ' ' + j_y_size + ' 0 0 0 ' + circle_j_x + ' ' + circle_j_y + '"' + stroke + ' stroke="' + draw + '" stroke-width="' + thickness + '" fill="none"/>'
              } else {
                bended_edge += 'L ' + circle_j_x + ' ' + circle_j_y + '"' + stroke + ' stroke="' + draw + '" stroke-width="' + thickness + '" fill="none"/>'
              }

              svg_edges.push(
                $(bended_edge)
              )
            }
          }
        }
      }
    }
  };

  var group_nodes = [];
  for(var g in groups) {
    var group = data['infos']['groups'][g]
    group_nodes = group_nodes.concat(group)
  };

  var global_values = []
  for(var k in nodes) {

    if (k !=  'start' && k != 'stop') {

      var node = data['elements']['nodes'][k];
      var node_genomes = Object.keys(node['genome']);
      var node_genomes_length = node_genomes.length
      
      var k_x = parseInt(node['position']['x'])
      var k_y = parseInt(node['position']['y'])

      if (!group_nodes.includes(k)) {
        var node_class = 'class="node'
      } else {
        var node_class = 'stroke-opacity="0" fill-opacity="0" class="pseudo'
      }

      var color = pickcolor (edgecoloring, Object.keys(node['genome']))
      var draw = lighter_color('#ffffff', color, node_genomes_length / genome_size);
      
      var [graph_size, graph_start, graph_stop] = outer_layers['graph']
      var k_y_size = sum_middle_layer + graph_start + graph_size * 0.5 + k_y * sum_outer_layer
      
      var [circle_k_x, circle_k_y] = transform(k_x-0.5, k_y_size, theta);

      svg_nodes.push(
        $('<circle ' + node_class + '" id="' + k + '" cx="' + circle_k_x + '" cy="' + circle_k_y + '" r="' + node_size + '" fill="' + lighter_color('#ffffff', node_color, node_genomes_length / genome_size) + '" stroke="' + draw + '" stroke-width="' + node_thickness + '"/>')
      )

      var add_start = 1
      var add_stop = 0

      var [search_size, search_start, search_stop] = middle_layers['search']
      
      var [e_x, e_y] = transform(parseInt(k_x)-add_start, search_start, theta)
      var [f_x, f_y] = transform(parseInt(k_x)+add_stop, search_start, theta)
      var [g_x, g_y] = transform(parseInt(k_x)-add_start, search_stop, theta)
      var [h_x, h_y] = transform(parseInt(k_x)+add_stop, search_stop, theta)

      if (!global_values.includes(k_x)) {
        svg_search.push(
          $('<path class="marker" id="' + k_x + '" d="' +
          'M ' + e_x + ' ' + e_y + ' ' +
          'A ' + search_start  + ' ' + search_start + ' 0 0 0 ' + f_x + ' ' + f_y + ' ' +
          'L' + h_x + ' ' + h_y +  ' ' +
          'A ' + search_stop  + ' ' + search_stop + ' 0 0 1 ' + g_x + ' ' + g_y + ' ' +
          'L' + e_x + ' ' + e_y +
          '" fill="white" stroke="" stroke-width="0"/>')
        )
      }

      for (var genome of genomes) {

        var layer_name = genome + 'layer'
        var [layer_width, layer_start, layer_stop] = middle_layers[layer_name]

        var [a_x, a_y] = transform(parseInt(k_x)-add_start, layer_start, theta)
        var [b_x, b_y] = transform(parseInt(k_x)+add_stop, layer_start, theta)
        var [c_x, c_y] = transform(parseInt(k_x)-add_start, layer_stop, theta)
        var [d_x, d_y] = transform(parseInt(k_x)+add_stop, layer_stop, theta)

        if (!global_values.includes(k_x)) {
          svg_heatmaps.push(
            $('<path d="' +
            'M ' + a_x + ' ' + a_y + ' ' +
            'A ' + (layer_start) + ' ' + (layer_start) + ' 0 0 0 ' + b_x + ' ' + b_y + ' ' +
            'L' + d_x + ' ' + d_y +  ' ' +
            'A ' + (layer_stop) + ' ' + (layer_stop) + ' 0 0 1 ' + c_x + ' ' + c_y + ' ' +
            'L' + a_x + ' ' + a_y +
            // '" fill="' + lighter_color('#00ff00', '#ff0000', mean_entropy[key] / max) + '" stroke="" stroke-width="2"/>')
            '" fill="lightgray" stroke="" stroke-width="0"/>')
          )
        }
      }

      for (var layer_name of layers) {

        if ($('#flex' + layer_name).prop('checked') == true){

          var layer_scale = data['infos']['layers_data'][layer_name]['scale']
          var value = data['infos']['layers_data'][layer_name][k]
          var max = data['infos']['layers_data'][layer_name]['max']

          if (layer_scale == 'local'){
            var [layer_width, layer_start, layer_stop] = outer_layers[layer_name]
            var k_y_size = sum_middle_layer + k_y * sum_outer_layer
          } else {
            var [layer_width, layer_start, layer_stop] = middle_layers[layer_name]
            var k_y_size = 0
          }

          var [a_x, a_y] = transform(parseInt(k_x)-add_start, layer_start + k_y_size, theta)
          var [b_x, b_y] = transform(parseInt(k_x)+add_stop, layer_start + k_y_size, theta)
          var [c_x, c_y] = transform(parseInt(k_x)-add_start, layer_stop + k_y_size, theta)
          var [d_x, d_y] = transform(parseInt(k_x)+add_stop, layer_stop + k_y_size, theta)

          if (layer_scale == 'local' || !global_values.includes(k_x) || (layer_scale == 'local' && !global_values.includes(k_x))) {
            svg_heatmaps.push(
              $('<path class="' + layer_name + '" xpos="' + k_x + '" name="' + (value / max).toFixed(3) + '" d="' +
              'M ' + a_x + ' ' + a_y + ' ' +
              'A ' + (layer_start + k_y_size) + ' ' + (layer_start + k_y_size) + ' 0 0 0 ' + b_x + ' ' + b_y + ' ' +
              'L' + d_x + ' ' + d_y +  ' ' +
              'A ' + (layer_stop + k_y_size) + ' ' + (layer_stop + k_y_size) + ' 0 0 1 ' + c_x + ' ' + c_y + ' ' +
              'L' + a_x + ' ' + a_y +
              // '" fill="' + lighter_color('#00ff00', '#ff0000', mean_entropy[key] / max) + '" stroke="" stroke-width="2"/>')
              '" fill="' + lighter_color('#00ff00', '#ff0000', value / max) + '" stroke="" stroke-width="0"/>')
            )
          }
        }
      }

    }

    global_values.push(k_x)

  };

  for(var l in groups) {

    var group = data['infos']['groups'][l]
    // console.log(group)

    var group_length = group.length
    var left_node_name = group[0]
    var right_node_name = group[group_length-1]

    var left_node = data['elements']['nodes'][left_node_name];
    var right_node = data['elements']['nodes'][right_node_name];

    var group_genomes = Object.keys(left_node['genome']);
    var group_genomes_length = group_genomes.length;

    var color = pickcolor (edgecoloring, Object.keys(left_node['genome']))
    var draw = lighter_color('#ffffff', color, group_genomes_length / genome_size);
    
    var l_x = parseInt(left_node['position']['x'])
    var l_y = parseInt(left_node['position']['y'])

    var m_x = parseInt(right_node['position']['x'])
    var m_y = parseInt(right_node['position']['y'])

    var [graph_size, graph_start, graph_stop] = outer_layers['graph']
    var l_y_size = sum_middle_layer + graph_start + graph_size * 0.5 + l_y * sum_outer_layer
    var m_y_size = sum_middle_layer + graph_start + graph_size * 0.5 + m_y * sum_outer_layer

    var [circle_t_x, circle_t_y] = transform(l_x-0.5, l_y_size + node_size, theta);
    var [circle_u_x, circle_u_y] = transform(l_x-0.5, l_y_size - node_size, theta);
    var [circle_v_x, circle_v_y] = transform(m_x-0.5, m_y_size + node_size, theta);
    var [circle_w_x, circle_w_y] = transform(m_x-0.5, m_y_size - node_size, theta);

    if ((m_x - l_x) * theta >= 180) {
      var arc_flag = 1
    } else {
      var arc_flag = 0
    }

    svg_groups.push(
      $('<path class="group" id="' + l + '" d="' +
      'M ' + circle_t_x + ' ' + circle_t_y + ' ' +
      'A ' + (l_y_size + node_size) + ' ' + (l_y_size + node_size) + ' 0 ' + arc_flag + ' 0 ' + circle_v_x + ' ' + circle_v_y + ' ' +
      'A ' + node_size + ' ' + node_size + ' 0 0 0 ' + circle_w_x + ' ' + circle_w_y + ' ' +
      'A ' + (m_y_size - node_size) + ' ' + (m_y_size - node_size) + ' 0 ' + arc_flag + ' 1 ' + circle_u_x + ' ' + circle_u_y + ' ' +
      'A ' + node_size + ' ' + node_size + ' 0 0 0 ' + circle_t_x + ' ' + circle_t_y +
      '" fill="' + lighter_color('#ffffff', group_color, group_genomes_length / genome_size) + '" stroke="' + draw + '" stroke-width="' + node_thickness + '"/>')
    )
  };

  for (var layer_name of layers) {

    if ($('#flex' + layer_name).prop('checked') == true){
      var layer_scale = data['infos']['layers_data'][layer_name]['scale']
      
      if (layer_scale == 'local'){
        var [layer_width, layer_start, layer_stop] = outer_layers[layer_name]
        var y_size = sum_middle_layer + layer_width * 0.5
      } else {
        var [layer_width, layer_start, layer_stop] = middle_layers[layer_name]
        var y_size = layer_width * 0.5
      }

      var [circle_x, circle_y] = transform(0, (layer_start + y_size), theta)
      svg_heatmaps.push(
        $('<text text-anchor="end" transform="translate (-10)" dominant-baseline="middle" x="' + circle_x + '" y="' + circle_y + '" dy="0" font-size="' + $('#label')[0].value + '" font-family="sans-serif" fill="black">' + layer_name + '</text>')
      )
    }
  }

  for (var genome_name of genomes) {
    if ($('#flex' + genome_name + 'layer').prop('checked') == true){
      var [layer_width, layer_start, layer_stop] = middle_layers[genome_name + 'layer']
      
      if (layer_width >= edge_thickness) {
      
        var y_size = layer_start + layer_width * 0.5

        var [circle_x, circle_y] = transform(0, y_size, theta)
        svg_edges.push(
          $('<text text-anchor="end" transform="translate (-10)" dominant-baseline="middle" x="' + circle_x + '" y="' + circle_y + '" dy="0" font-size="' + $('#label')[0].value + '" font-family="sans-serif" fill="black">' + genome_name + '</text>')
        )
      }
    }
  }

  for (var item of svg_search) svg_core.append(item);
  for (var item of svg_heatmaps) svg_core.append(item);
  for (var item of svg_edges) svg_core.append(item);
  for (var item of svg_nodes) svg_core.append(item);
  for (var item of svg_groups) svg_core.append(item);

  body.append(svg_core)
  body.html(body.html());
}

//ANCHOR - Check node
function checknode(searchpos, positions, node, searchfunction, expressiondrop, expressioncomparison) {

  var append = true

  if (searchpos == true) {

    if (!positions.includes(parseInt(node['position']['x']))) {
      append = false
    }
  }

  var keys = Object.keys(searchfunction)
  if (keys.length > 0) {
    var t = false
    var d = get_gene_cluster_consensus_functions(node['genome'])
    for (var source of keys) {
      for (var s of searchfunction[source]){

        if ((d[source][2].toLowerCase().includes(s.trim().toLowerCase()))) {
          t = true
        }
      }
    }
    if (t == false){
      append = false
    }
  }
    
  // if (append == true){
  //   console.log(node)
  // }

  return append
}

function main () {

  $.ajax({
    url: "/pangraph/get_json",
    type: "POST",
    cache: false,
    //data: JSON.stringify(data),
    contentType: "application/json",
    dataType: "json",
    success: function(data){

    var layers = data['infos']['layers_names']
  
    for (var layer_name of layers) {
    
      // var layer_name = data['infos']['layers_names'][layer]
      var layer_scale = data['infos']['layers_data'][layer_name]['scale']

      if (layer_scale == 'global') {
        var value = '50'
      } else {
        var value = '10'
      }
    
      if ($('#flex' + layer_name + '').length > 0) {
        // console.log(layer_name)
      } else {
        var element = $('<div class="col-12 d-flex mb-1"></div>').append(
          $('<div class="col-2 d-flex align-items-center"></div>').append(
            $('<div class="form-switch d-flex"></div>').append(
              $('<input class="" type="checkbox" id="flex' + layer_name + '" name="flex' + layer_name + '" aria-label="..." data-toggle="tooltip" data-placement="top" title="Tooltip on top" checked>')
            )
          )
        ).append(
            $('<div class="col-8 d-flex align-items-center"></div>').append(
              layer_name
            )
        ).append(
          $('<div class="d-flex col-2"></div>').append(
            $('<input type="text" class="form-control float-end text-end flex-fill p-0 border-0" style= "background-color: #e9ecef;" id="' + layer_name + '" name="' + layer_name + '" value=' + value + ' aria-label="..." data-toggle="tooltip" data-placement="top" title="Choose your color">')
          )
        )

        $('#layers').append(element)
      }

    }

    // This title must update with Generic data from JSON
    $('#title-panel-first-line').text(data['infos']['meta']['title']);
    $('#title-panel-second-line').text('Pangraph Detail');

    // It seems unused function after UI changes
    $('#redraw').on('click', function() {

      var data = new Object;

      data['condtr'] = $('#condtr')[0].value
      data['maxlength'] = $('#maxlength')[0].value

      $("#genomecolors :input[type='checkbox']").each((index, element) => {

        var genome = $(element).attr('name')
        if ($(element).prop('checked') == true){
          data[genome] = 'on'
        } else {
          data[genome] = 'off'
        }
      })

      $.ajax({
          url: "/pangraph/settings",
          type: "POST",
          async: false,
          data: JSON.stringify(data),
          contentType: "application/json",
          dataType: "json",
          success: function(){
            $(document).off("click", ".group_choice");
            $(document).off("click", ".binremove");
            $(document).off("click", ".binchoice li a");
            $(document).off("change", ".colorchange");
            $(document).find("*").off();
            main()
          }
      });
    })

    // console.log(data)

    if (!$('#genomecolors').children().length) {

      $("#condtr")[0].value = data['infos']['gene_cluster_grouping_threshold']
      $("#maxlength")[0].value = data['infos']['max_edge_length_filter']

      for (var [genome, value] of Object.entries(data['infos']['genomes'])) {

        var state = ''
        if (value == 'on') {
          state = ' checked'
        }

        $('#genomecolors').append(
          $('<div class="col-12 d-flex mb-1">').append(
              $('<div class="col-2 d-flex align-items-center">').append(
                $('<div class="form-switch d-flex">').append(
                  $('<input class="" type="checkbox" id="flex' + genome + '" name="' + genome + '" aria-label="..." data-bs-toggle="tooltip" data-bs-placement="top" title="Tooltip on top"' + state + ' disabled>')
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
                $('<input type="color" class="form-control form-control-color flex-fill p-0 border-0" id="' + genome + '" name="' + genome + '" aria-label="..." data-bs-toggle="tooltip" data-bs-placement="top" title="Choose your color">')
              )
            )
          )
        // )

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

        if ($('#flex' + genome + 'layer').length > 0) {
          // console.log(layer_name)
        } else {
          var element = $('<div class="col-12 d-flex mb-1"></div>').append(
            $('<div class="col-2 d-flex align-items-center"></div>').append(
              $('<div class="form-switch d-flex"></div>').append(
                $('<input class="" type="checkbox" id="flex' + genome + 'layer" name="flex' + genome + 'layer" aria-label="..." data-toggle="tooltip" data-placement="top" title="Tooltip on top" checked>')
              )
            )
          ).append(
              $('<div class="col-8 d-flex align-items-center"></div>').append(
                genome
              )
          ).append(
            $('<div class="d-flex col-2"></div>').append(
              $('<input type="text" class="form-control float-end text-end flex-fill p-0 border-0" style= "background-color: #e9ecef;" id="' + genome + 'layer" name="' + genome + 'layer" value=50 aria-label="..." data-toggle="tooltip" data-placement="top" title="Choose your color">')
            )
          )

          $('#layers').append(element)
        }
      }
    }

    var body = $('#svgbox')
    body.empty()

    // learn about the functional annotation sources
    functional_annotation_sources_available = data['infos']['functional_annotation_sources_available'];
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

    var layers = data['infos']['layers_names']

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

    for (var layer_name of layers) {
      if ($('#flex' + layer_name).prop('checked') == true){

        $('#expressiondrop').append($('<option value="' + layer_name + '">' + layer_name + '</option>'))
        $('#filter').append(
          $('<div class="col-12"></div>').append(
            $('<div class="row align-items-center"></div>').append(
              $('<div class="col-2 mb-1"></div>').append(
                $('<input class="" type="checkbox" id="min' + layer_name + '" value="" aria-label="..." data-toggle="tooltip" data-placement="top" title="Tooltip on top">')
              )
            ).append(  
              $('<div class="col-8 mb-1"></div>').append(
                'Min ' + layer_name
              )
            ).append(
              $('<div class="col-2 mb-1"></div>').append(
                $('<input type="text" class="form-control flex-fill p-0 border-0" style= "background-color: #e9ecef;" value="" id="min' + layer_name + 'text" aria-describedby="">')
              )  
            ).append(
              $('<div class="col-2 mb-1"></div>').append(
                $('<input class="" type="checkbox" id="max' + layer_name + '" value="" aria-label="..." data-toggle="tooltip" data-placement="top" title="Tooltip on top">')
              )
            ).append( 
              $('<div class="col-8 mb-1"></div>').append(
                'Max ' + layer_name
              )
            ).append(
              $('<div class="col-2 mb-1"></div>').append(
                $('<input type="text" class="form-control flex-fill p-0 border-0" style= "background-color: #e9ecef;" value="" id="max' + layer_name + 'text" aria-describedby="">')
              )
            )
          )
        )
      }
    }

    generate_svg(body, data);

    window.zoomSVG = svgPanZoom('#result', {
      zoomEnabled: true,
      panEnabled: false,
      controlIconsEnabled: false,
      minZoom: 0.1,
      maxZoom: 100
    });

    //ANCHOR - Zoom pan functions
    // $('#plus').on('click', function() {
    //   window.zoomSVG.zoomIn();
    // })

    // $('#minus').on('click', function() {
    //   window.zoomSVG.zoomOut();
    // })

    $('#fit').on('click', function() {
      window.zoomSVG.resize();
      window.zoomSVG.fit();
      window.zoomSVG.center();
    })

    //ANCHOR - Main panel response functions
    
    if ($('#condtr')[0].value == -1){
      // $('#customRange2').prop('disabled', true);
      $('#flexcondtr').prop('checked', false);
    }

    if ($('#maxlength')[0].value == -1){
      // $('#customRange3').prop('disabled', true);
      $('#flexmaxlength').prop('checked', false);
    }

    if ($('#groupcompress')[0].value == -1){
      $('#flexgroupcompress').prop('checked', false);
    }

    $('#flexcondtr').change(function() {
      if ($(this).prop('checked') == true){
        $('#condtr')[0].value = 2;
        // $('#customRange2')[0].value = 2;
        // $('#customRange2').prop('disabled', false);
      } else {
        $('#condtr')[0].value = -1;
        // $('#customRange2')[0].value = 2;
        // $('#customRange2').prop('disabled', true);
      }
    })

    $('#flexmaxlength').change(function() {
      if ($(this).prop('checked') == true){
        $('#maxlength')[0].value = 1;
        // $('#customRange3')[0].value = 1;
        // $('#customRange3').prop('disabled', false);
      } else {
        $('#maxlength')[0].value = -1;
        // $('#customRange3')[0].value = 1;
        // $('#customRange3').prop('disabled', true);
      }
    })

    $('#flexgroupcompress').change(function() {
      if ($(this).prop('checked') == true){
        $('#groupcompress')[0].value = 0;
        // $('#customRange4')[0].value = 0;
        // $('#customRange4').prop('disabled', false);
      } else {
        $('#groupcompress')[0].value = -1;
        // $('#customRange4')[0].value = 0;
        // $('#customRange4').prop('disabled', true);
      }
    })

    sortable('#genomecolors', {
      forcePlaceholderSize: true,
      handle: '.user-handle',
      items: 'div'
    });

    $(document).on("click", ".group_choice", async function() {

      // console.log(this)

      var gene_cluster_id = this.getAttribute("name_id");
      var gene_cluster_context = this.getAttribute("group");
      var add_align = parseInt(this.getAttribute("context"))

      var body = $(this).parent().parent().parent().parent()

      body.empty()
      var bodyinfo = $('<div class="card-body overflow-scroll"></div>')
      body.append(bodyinfo)

      // console.log(gene_cluster_id, gene_cluster_context)

      var all_info = await get_gene_cluster_display_tables(gene_cluster_id, gene_cluster_context, data, add_align)

      // console.log(all_info)

      bodyinfo.append(all_info)
      // nodeinfo(e.target, data);
    })

    //ANCHOR - Recolor nodes after redraw
    var groups = data['infos']['groups']
    for (var binid of Object.keys(bins)) {
      var nodes = bins[binid];
      var updated_nodes = []

      for (var node of nodes) {

        var name = document.getElementById(node);

        if(name) {
          if (name.getAttribute('class') == 'pseudo') {
            for(var g in groups) {
              var group = data['infos']['groups'][g]
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
          updated_nodes.push(...current_groups[node])
        };

      }

      bins[binid] = updated_nodes
      // console.log(nodes)
      // console.log(updated_nodes)
      for (var node of bins[binid]) {

        bins[binid] = bins[binid].filter(item => item !== node);
        var name = document.getElementById(node);
        bins = marknode(name, data, binid, bins);

      }
    }

    // console.log(bins)

    current_groups = data['infos']['groups']

    // console.log(current_groups)

    //ANCHOR - Bin dropdown choice function
    $(document).on("click", ".binremove", function() {

      var id = $(this).attr('name_id');
      var binid = $(this).attr('bin');

      var name = document.getElementById(id);
      // console.log(name)
      bins = marknode(name, data, binid, bins);

      $(this).parent().parent().parent().remove();

      // console.log(bins)

    })
    
    // ANCHOR - Add bin
    $('#binadd').on('click', function() {
        binnum += 1;
    
        $('#bingrid').append(
            $('<div class="col-12" id="bin' + binnum + '"></div>').append(
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
    
        bins['bin' + binnum] = [];
    });
    

    //ANCHOR - Remove bin
    $('#binremove').on('click', function() {

      var selection = document.querySelector('input[name="binradio"]:checked');

      if (selection !== null) {
        var binid = selection.value;

        for (var node of bins[binid]) {
          var name = document.getElementById(node);
          bins = marknode(name, data, binid, bins);
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

    //ANCHOR - Change bin
    $(document).on("change", ".colorchange", function() {

      var binid = this.name;
      var nodes = bins[binid];

      for (var node of nodes) {

        bins[binid] = bins[binid].filter(item => item !== node);
        var name = document.getElementById(node);
        bins = marknode(name, data, binid, bins);

      }

    });

    //ANCHOR - Info bin
    $('#bininfo').on('click', async function() {

      var selection = document.querySelector('input[name="binradio"]:checked');

      // if (selection !== null) {
      var binid = selection.value;
      var appendlist = [];

      $('#BinModalBody').empty();
      for (var id of bins[binid]) {
        var element = document.getElementById(id);

        if (element.getAttribute('class') == 'group') {
          gene_cluster_context = id;
          gene_cluster_id = data['infos']['groups'][id][0]
        } else {
          gene_cluster_id = id;
          gene_cluster_context = null;
        }

        appendlist.push([gene_cluster_id, gene_cluster_context]);

      }

      for (var [gene_cluster_id, gene_cluster_context] of appendlist) {

        var body = $('<div></div>').append(
          $('<div class="card-body overflow-scroll" id="' + gene_cluster_id + 'div"></div>').append(
            await get_gene_cluster_display_tables(gene_cluster_id, gene_cluster_context, data, add_align=0)
          )
        )

        if (gene_cluster_context == null) {
          var element_id = gene_cluster_id
        } else {
          var element_id = gene_cluster_context
        }

        $('#BinModalBody').append(
          $('<div class="card border mb-3 w-100"></div>').append(
            $('<div class="card-header"></div>').append(
              $('<div class="d-flex justify-content-end align-items-center"></div>').append(
                $('<button class="btn-close binremove" name_id="' + element_id + '" bin="' + binid + '"></button>')
              )
            )
          ).append(
            body
          )
        )

        // basic_info = {'Gene Cluster': id, 'Genomes': genomes, 'Position': position}
        // body.append(get_gene_cluster_display_tables(id, basic_info, gene_cluster_data))
      }

      $('#BinModal').modal('show');
      // }
    })

    $('#InfoDownload').on('click', function() {

      // Variable to store the final csv data
      var csv_data = [];
      var basics = $('#node_basics_table')
      var title = basics[0].getAttribute("gc_id")

      var functions = $('#node_functions_table')

      var basics_rows = basics[0].getElementsByTagName('tr');
      var function_rows = functions[0].getElementsByTagName('tr');

      for (let i = 0; i < function_rows.length; i++) {
    
          if (i >= basics_rows.length) {
            var basics_cols = basics_rows[1].querySelectorAll('td,th');
            var function_cols = function_rows[i].querySelectorAll('td,th');
          } else { 
            var basics_cols = basics_rows[i].querySelectorAll('td,th');
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
      var title = alignment[0].getAttribute("gc_id")
      
      var alignment_rows = alignment[0].getElementsByTagName('tr');

      for (let i = 1; i < alignment_rows.length; i++) {
    
          var alignment_cols = alignment_rows[i].querySelectorAll('td,th');
    
          csv_data += ">Genome:" + alignment_cols[0].innerHTML +"|Genecall:" + alignment_cols[1].innerHTML + "|Contig:" + alignment_cols[2].innerHTML + "|Direction:" + alignment_cols[3].innerHTML + '\n';
          var genome = ''
          
          var alignment_nucs = alignment_cols[4].getElementsByTagName('span');
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

    $('#svgDownload').on('click', function() {
      var blob = new Blob([$('#svgbox')[0].innerHTML]);
      var title = data['infos']['meta']['title']
      downloadBlob(blob, title + ".svg");
    });

    $('#searchadd').on('click', function() {

      var selection = document.querySelector('input[name="binradio"]:checked')
      var binid = selection.value

      for (var [id, members] of Object.entries(searched)) {
        if (!(bins[binid].includes(id))) {

          var e = document.getElementById(id);
          bins = marknode(e, data, binid, bins);

        }
      }
    })

    $('#searchremove').on('click', function() {

      var selection = document.querySelector('input[name="binradio"]:checked')
      var binid = selection.value

      for (var [id, members] of Object.entries(searched)) {
        if ((bins[binid].includes(id))) {

          var e = document.getElementById(id);
          bins = marknode(e, data, binid, bins);

        }
      }
    })

    $('#searcherase').on('click', function() {

      for (var [id, members] of Object.entries(searched)) {
        for (var mem of members) {
          var xpos = data['elements']['nodes'][mem]['position']['x']
          var e = document.getElementById(xpos);
          e.setAttribute("fill", "white")
        }
      }

      searched = {};
    })

    $('#searchcolor').on('click', function() {

      for (var [id, members] of Object.entries(searched)) {
        for (var mem of members) {
          var xpos = data['elements']['nodes'][mem]['position']['x']
          var e = document.getElementById(xpos);
          e.setAttribute("fill", "#ff0000")
        }
      }
    })

    const searchEraseButton = document.getElementById('searcherase');
    const searchColorButton = document.getElementById('searchcolor');
    const searchAppendBin = document.getElementById('searchadd');
    const searchRemoveBin = document.getElementById('searchremove');
    
    var searched = {}
    $('#search').on('click', function() {

      var searchpos = false

      //Filter Search
      var layers_filter = new Object
      var mingenomes = ($("#mingenomes").prop('checked') == true && !isNaN(parseFloat($("#mingenomestext")[0].value))) ? parseFloat($("#mingenomestext")[0].value) : -1
      var maxgenomes = ($("#maxgenomes").prop('checked') == true && !isNaN(parseFloat($("#maxgenomestext")[0].value))) ? parseFloat($("#maxgenomestext")[0].value) : -1
      var minposition = ($("#minposition").prop('checked') == true && !isNaN(parseFloat($("#minpositiontext")[0].value))) ? parseFloat($("#minpositiontext")[0].value) : -1
      var maxposition = ($("#maxposition").prop('checked') == true && !isNaN(parseFloat($("#maxpositiontext")[0].value))) ? parseFloat($("#maxpositiontext")[0].value) : -1
      
      layers_filter['genomes'] = {'min': mingenomes, 'max': maxgenomes}
      layers_filter['position'] = {'min': minposition, 'max': maxposition}

      var layers = data['infos']['layers_names']
      for (var layer_name of layers) {

        // console.log("#min" + layer_name + "text")
        var minlayer = ($("#min" + layer_name).prop('checked') == true && !isNaN(parseFloat($("#min" + layer_name + "text")[0].value))) ? parseFloat($("#min" + layer_name + "text")[0].value) : -1
        var maxlayer = ($("#max" + layer_name).prop('checked') == true && !isNaN(parseFloat($("#max" + layer_name + "text")[0].value))) ? parseFloat($("#max" + layer_name + "text")[0].value) : -1
        
        layers_filter[layer_name] = {'min': minlayer, 'max': maxlayer}
      }

      //Expression Search
      var expressioncomparison = ''
      var expressiondrop = $('#expressiondrop')[0].value
      var expressionrel = $('#expressionrel')[0].value
      var expressiontext = $('#expressiontext')[0].value

      console.log(expressiondrop, expressionrel, expressiontext)

      if (expressiondrop != "Choose item" && expressionrel != "Choose operator" && expressiontext != '') {
        if (expressionrel == '=') {
          if (!isNaN(expressiontext)) {
            expressioncomparison = '== ' + expressiontext
          } else {
            expressioncomparison = '== "' + expressiontext + '"'
          }
        } else if (expressionrel == '\u{2260}') {
          if (!isNaN(expressiontext)) {
            expressioncomparison = '!= ' + expressiontext
          } else {
            expressioncomparison = '!= "' + expressiontext + '"'
          }
        } else if (expressionrel == '\u{2264}' && !isNaN(expressiontext)) {
          expressioncomparison = '<= ' + expressiontext
        } else if (expressionrel == '\u{2265}' && !isNaN(expressiontext)) {
          expressioncomparison = '>= ' + expressiontext
        } else if (expressionrel == '\u{003C}' && !isNaN(expressiontext)) {
          expressioncomparison = '< ' + expressiontext
        } else if (expressionrel == '\u{003E}' && !isNaN(expressiontext)) {
          expressioncomparison = '> ' + expressiontext
        } else if (expressionrel == '\u{25C2}\u{25AA}\u{25B8}') {
          expressioncomparison = '.includes("' + expressiontext + '")'
        } else if (expressionrel == '\u{25C2}\u{25AA}') {
          expressioncomparison = '.endsWith("' + expressiontext + '")'
        } else if (expressionrel == '\u{25AA}\u{25B8}') {
          expressioncomparison = '.startsWith("' + expressiontext + '")'
        }
      } 

      console.log(expressioncomparison)

      //Function Search
      var searchfunction = {}
      var searchterms = $('#searchFunctionsValue')[0].value.split(",");
      for (var source of functional_annotation_sources_available) {
        if ($("#flex" + source).prop('checked') == true) {
          searchfunction[source] = searchterms
        }
      }

      var layers_positions = {}
      for (var layer_name of Object.keys(layers_filter)) {
        
        var minlayer = layers_filter[layer_name]['min']
        var maxlayer = layers_filter[layer_name]['max']

        if (minlayer != -1 || maxlayer != -1 || (minlayer != -1 && maxlayer != -1)) {

          layers_positions[layer_name] = []
          searchpos = true
          if (layer_name == 'position') {
            var global_x = data["infos"]["meta"]["global_x"] -1;
            var layerobjects = new Array(global_x - 1).fill().map((d, i) => i + 1);
            // console.log(layerobjects)
          } else {
            var layerobjects = document.querySelectorAll("." + layer_name)
          }

          for (var o of layerobjects) {

            var append = true
            
            if (layer_name == 'position') {
              var value = o
              var result = o
            } else {
              var value = o.getAttribute("name")
              var result = o.getAttribute("xpos")
            }

            if (minlayer != '-1'){
              if (eval(value + '<' + minlayer)){
                append = false
              }
            }

            if (maxlayer != '-1'){
              if (eval(value + '>' + maxlayer)){
                append = false
              }
            }

            if (append == true) {
              layers_positions[layer_name].push(parseInt(result))
            }

          }
        }
      }

      // console.log(layers_positions)

      var positions = []
      var keys = Object.keys(layers_positions)

      if (keys.length > 0) {
        for (var pos of layers_positions[keys[0]]) {

          if (keys.length > 0) {
            var add = true
            for (let k = 1; k < keys.length; k++) {
              if (!layers_positions[keys[k]].includes(pos)) {
                add = false
              }
            }

            if (add == true) {
              positions.push(pos)
            }
          } else {
            positions.push(pos)
          }
        }
      }
        
      // console.log(positions)

      if (expressioncomparison == '' && Object.keys(searchfunction).length == 0 && searchpos == false) {
        var message = "Please enter valid search parameters."
      } else {
        // console.log(expressioncomparison, searchfunction, searchpos)

        var nodes = document.querySelectorAll(".node")
        for (var node of nodes) {

          var id = node.getAttribute("id")
          var node = data['elements']['nodes'][id]

          if (checknode(searchpos, positions, node, searchfunction, expressiondrop, expressioncomparison) == true) {
            if (expressiondrop == "Name") {
              if (eval('"' + node["name"] + '"' + expressioncomparison)) {
                if (!(id in searched)) {
                  searched[id] = [id]
                }
              }
            } else {
              if (!(id in searched)) {
                searched[id] = [id]
              }
            }
          }
        }

        var groups = document.querySelectorAll(".group")
        for (var group of groups) {
          var groupid = group.getAttribute("id")
          var members = data["infos"]["groups"][groupid]
          for (var id of members) {

            node = data['elements']['nodes'][id]

            if (checknode(searchpos, positions, node, searchfunction, expressiondrop, expressioncomparison) == true) {
              
              if (expressiondrop == "Name") {
                if (eval('"' + node["name"] + '"' + expressioncomparison) || eval('"' + group + '"' + expressioncomparison)) {

                  if (!(groupid in searched)) {
                    searched[groupid] = [id]
                  } else {
                    searched[groupid].push(id)
                  }
                }
              } else {
                if (!(groupid in searched)) {
                  searched[groupid] = [id]
                } else {
                  searched[groupid].push(id)
                }
              }
            }
          }
        }
        
        var table = $('<div class="form-horizontal"></div>')

        for (var key of Object.keys(searched)) {
          
          table.append(
            $('<div class="col-12"></div>').append(
              $('<div class="row align-items-center"></div>').append(
                $('<div class="col-12 mb-1">' + searched[key] + '</td>')
              )
            )
          )
        }

        $('#searchtable').empty()
        $('#searchtable').append(
          $('<div class="pt-3"></div>').append(
            $('<div class="shadow-box pb-3 pr-3 pl-3 mb-3 rounded search-box-filter"></div>').append(
              table
            )
          )
        )

        var message = 'You have ' + Object.keys(searched).length + ' item(s) in your queue.'

        if (Object.keys(searched).length != 0) {
          searchEraseButton.disabled = false;
          searchColorButton.disabled = false;
          searchAppendBin.disabled = false;
          searchRemoveBin.disabled = false;
        }
      }

      var toastbody = $('#searchtoastbody')
      toastbody.empty()
      toastbody.append(
        message
      )
      // var searchtoast = bootstrap.Toast.getOrCreateInstance($('#searchtoast'))
      $('#searchtoast').toast('show')
    })

    var nodes = document.querySelectorAll(".node")
    var divs = document.querySelectorAll(".node, .group");
    for (var el of divs) {

      if (el.getAttribute("class") == 'group'){
        var id = data["infos"]["groups"][el.getAttribute("id")][0]
        var name = el.getAttribute("id")
      } else {
        var id = el.getAttribute("id")
        var name = data['elements']['nodes'][el.getAttribute("id")]['name']
      }

      tippy(el, {
        content: '<strong>' + name + '</strong>' + '<br />',
        allowHTML: true,
        onHide(instance) {
          if (instance.reference.id.startsWith('GCG_')){
            var id = data["infos"]["groups"][instance.reference.id][0]
          } else {
            var id = instance.reference.id
          }
          var elements = Object.keys(data['elements']['nodes'][id]['genome'])

          for (var element of elements) {
            $('#number_' + element)[0].innerText = '0';
          }
        },
        onShow(instance) {
          if (instance.reference.id.startsWith('GCG_')){
            var id = data["infos"]["groups"][instance.reference.id][0]
          } else {
            var id = instance.reference.id
          }
          var elements = Object.keys(data['elements']['nodes'][id]['genome'])

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

    $("#svgbox").on('mousedown', function(e) {
      old_xpos = e.offsetX
      old_ypos = e.offsetY

      xpos = old_xpos
      ypos = old_ypos

      isDown = true
      diff = 0
    })

    $("#svgbox").on('mousemove', function(e) {
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

    $("#svgbox").on('mouseup', function(e) {
      if (isDown === true) {

        var selection = document.querySelector('input[name="binradio"]:checked')

        isDown = false

        if (diff < 10) {
          if (e.target.getAttribute('class') === 'group' || e.target.getAttribute('class') === 'node') {

            if (e.shiftKey && selection !== null) {

              var binid = selection.value
              bins = marknode(e.target, data, binid, bins);

            } else {
              //show_node_details_modal(e.target, data);
              // console.log(e.target)
              nodeinfo(e.target, data);
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

            for (var n of nodes) {

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
                bins = marknode(n, data, binid, bins);

              }
            }

            var groups = data['infos']['groups']
            for(var g in groups) {
              // var inside = true;
              var group = data['infos']['groups'][g]
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
                  bins = marknode(name, data, binid, bins);
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
              //   bins = marknode(name, data, binid, bins);
              // }
            }
          }
        }
      }
    })
    document.body.addEventListener('keydown', function(ev) {
      if ((/^(?:input|select|textarea|button)$/i).test(ev.target.nodeName))
      {
          // shortcuts should not work if user is entering text to input.
          return false;
      }
      if (ev.keyCode === 83) { // S = 83
          $('#toggle-panel-left').trigger('click');
      }
      if (ev.keyCode === 84) { // T = 84
          $('#title-panel').toggle();
          $('#toggle-panel-top').toggleClass('invisible visible');
      }
    });
  }})

}

//ANCHOR - Main function after loading DOM
$(document).ready(function() {

  // console.log('start pangrah layout creation!')
  main()

});
