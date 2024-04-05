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

  return(d);
}


function get_gene_cluster_basics_table(gene_cluster_id, data) {
    // first, learn a few basics about the gene cluster to be displayed
    var position_in_graph = data['elements']['nodes'][gene_cluster_id]['position']['x'] + " / " + (data["infos"]["meta"]["global_x"] - 1);
    var num_contributing_genomes = Object.keys(data['elements']['nodes'][gene_cluster_id]['genome']).length + " / " + (data['infos']['num_genomes']);

    basic_info = {'Gene Cluster': gene_cluster_id, 'Contributing Genomes': num_contributing_genomes, 'Position in Graph': position_in_graph}
    // build the basic information table
    basic_info_table = `<p class="modal_header">Basics</p>`;
    basic_info_table += `<table class="table table-striped table-bordered sortable" id="node_basics_table">`;
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


function get_gene_cluter_functions_table(gene_cluster_id, data) {
    functions_table = `<p class="modal_header">Consensus functional annotations</p>`;
    functions_table += `<table class="table table-striped sortable" id="node_functions_table">`;
    functions_table += `<thead class="thead-dark"><tr>`;
    functions_table += `<th scope="col">Source</th>`;
    functions_table += `<th scope="col">Accession</th>`;
    functions_table += `<th scope="col">Function</th>`;
    functions_table += `<th scope="col">Agreement</th>`;
    functions_table += `</tr></thead><tbody>\n\n`;

    var gene_cluster_data = data['elements']['nodes'][gene_cluster_id]['genome']
    var d = get_gene_cluster_consensus_functions(gene_cluster_data);

    function_sources = Object.keys(d).sort();
    console.log(function_sources);
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


function get_gene_cluster_context_table(gene_cluster_id_current, gene_cluster_context) {
    if (gene_cluster_context == null)
        return '';

    gene_cluster_context_table = `<p class="modal_header">Gene cluster context</p>`;
    gene_cluster_context_table += `<div class="gene_cluster_context_items">`;
    for(index in gene_cluster_context) {
        gene_cluster_id = gene_cluster_context[index];
        if (gene_cluster_id == gene_cluster_id_current){
            gene_cluster_context_table += `<span class="gene_cluster_id gene_cluster_id_current">` + gene_cluster_id + `</span>`;
        } else {
            // FIXME: we will need to find a way to fix this <a> tag below
            gene_cluster_context_table += `<span class="gene_cluster_id"><a href="#">` + gene_cluster_id + `</a></span>`;
        }
    }
    gene_cluster_context_table += `</div>`;

    return gene_cluster_context_table;

}


//ANCHOR - Get tables for GC basic info and functions
function get_gene_cluster_display_tables(gene_cluster_id, gene_cluster_context, data) {
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

    gene_cluster_context_table = get_gene_cluster_context_table(gene_cluster_id, gene_cluster_context);

    ///////////////////////////////////////////////////////////////////////////////////////////
    // BUILD BASIC INFORMATION TABLE
    ///////////////////////////////////////////////////////////////////////////////////////////

    basic_info_table = get_gene_cluster_basics_table(gene_cluster_id, data);

    ///////////////////////////////////////////////////////////////////////////////////////////
    // BUILD FUNCTIONS TABLE
    ///////////////////////////////////////////////////////////////////////////////////////////

    functions_table = get_gene_cluter_functions_table(gene_cluster_id, data);

    ///////////////////////////////////////////////////////////////////////////////////////////
    // RETRIEVE AND BUILD SEQUENCE ALIGNMENTS
    ///////////////////////////////////////////////////////////////////////////////////////////
//
    var alignment = {}

    if (gene_cluster_id != 'start' && gene_cluster_id != 'stop') {
      for (var genome of Object.keys(data['elements']['nodes'][gene_cluster_id]['genome'])) {
        alignment[genome] = [data['elements']['nodes'][gene_cluster_id]['genome'][genome]['gene_call'], data['elements']['nodes'][gene_cluster_id]['name']]
      }
    }

    gene_cluster_sequence_alignments_table = appendalignment(alignment)

    ///////////////////////////////////////////////////////////////////////////////////////////
    // MERGE ALL AND RETURN
    ///////////////////////////////////////////////////////////////////////////////////////////

    return gene_cluster_context_table + basic_info_table + functions_table + gene_cluster_sequence_alignments_table
}

//ANCHOR - Fetch genecall functions
function fetchgenome(info) {

  var d = new Object();
  for (var source of functional_annotation_sources_available) {

    if (info[source] == 'None') {
      var func_id = 'None'
      var func = 'None'
      var evalue = 0
    } else {

      var func_id = info[source][0]
      var func = info[source][1]
      var evalue = info[source][2]
    }

    d[source] = [func_id, func, evalue]
  }
  return(d);
}

//ANCHOR - Append genecall functions
async function appendgenome(bodygenome, name, call, genome, length, direction, paralog, partial, info) {

  bodygenome.append(
    $('<div class="row gx-2"></div>').append(
      $('<div class="col-2"></div>').append(
        $('<b>Name</b>')
      )
    ).append(
      $('<div class="col-10"></div>').append(
        name
      )
    ).append(
      $('<div class="col-2"></div>').append(
        $('<b>Genecall</b>')
      )
    ).append(
      $('<div class="col-10" id="genecall"></div>').append(
        call
      )
    ).append(
      $('<div class="col-2"></div>').append(
        $('<b>Genome</b>')
      )
    ).append(
      $('<div class="col-10" id="genome"></div>').append(
        genome
      )
    ).append(
      $('<div class="col-2"></div>').append(
        $('<b>Length</b>')
      )
    ).append(
      $('<div class="col-10" id="length"></div>').append(
        length
      )
    ).append(
      $('<div class="col-2"></div>').append(
        $('<b>Direction</b>')
      )
    ).append(
      $('<div class="col-10" id="direction"></div>').append(
        direction
      )
    ).append(
      $('<div class="col-2"></div>').append(
        $('<b>Paralog</b>')
      )
    ).append(
      $('<div class="col-10" id="paralog"></div>').append(
        paralog
      )
    ).append(
      $('<div class="col-2"></div>').append(
        $('<b>Partial</b>')
      )
    ).append(
      $('<div class="col-10" id="partial"></div>').append(
        partial
      )
    )
  )

  var node = $('<div class="row gx-2"></div>').append(
    $('<div class="col-2"></div>').append(
      $('<b>Source</b>')
    )
  ).append(
    $('<div class="col-9"></div>').append(
      $('<div class="row g-0"></div>').append(
        $('<div class="col-2"></div>').append(
          $('<b>Accession</b>')
        )
      ).append(
        $('<div class="col-10"></div>').append(
          $('<b>Function</b>')
        )
      )
    )
  ).append(
    $('<div class="col-1"></div>').append(
      $('<b>e-value</b>')
    )
  )

  var d = fetchgenome(info)
  for (var [source, value]  of Object.entries(d)) {

    var [func_id, func, evalue] = value

    node.append(
      $('<div class="col-2"></div>').append(
        source
      )
    )

    var acc = $('<div class="row g-0"></div>')
    var list_func_id = func_id.split('!!!')
    var list_func = func.split('!!!')

    for (var i = 0; i < list_func_id.length; i++) {
      acc.append(
        $('<div class="col-2"></div>').append(
          list_func_id[i]
        )
      )
      acc.append(
        $('<div class="col-10"></div>').append(
          list_func[i]
        )
      )
    }

    node.append(
      $('<div/>', {
        class: 'col-9'
      }).append(acc)
    )

    node.append(
      $('<div/>', {
        class: 'col-1 text-end'
      }).append(evalue === '' ? '' : evalue.toFixed(3))
    )

  }

  bodygenome.append(
    $('<hr>')
  ).append(
    node
  )

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
async function appendalignment(alignment) {
    // FIXME: We need to come up with more accurate and explicit function names.
    if (Object.keys(alignment).length == 0) {
        return ''
    }

    alignments_table = `<p class="modal_header">Sequence alignments</p>`;
    alignments_table += `<table class="table table-striped sortable" id="node_sequence_alignments_table">`;
    alignments_table += `<thead><tr>`;
    alignments_table += `<th scope="col">Genome</th>`;
    alignments_table += `<th scope="col">Gene Call</th>`;
    alignments_table += `<th scope="col">Sequence</th>`;
    alignments_table += `</tr></thead><tbody>\n\n`;


    var d = await fetchalignment(alignment)

    for (var [genome, value] of Object.entries(d)) {
      var colored = value[1].replace(/A|R|N|D|C|Q|E|G|H|I|L|K|M|F|P|S|T|W|Y|V|-/gi, function(matched){return mapAS[matched];});

      alignments_table += `<tr>`
      alignments_table += `<td>` + genome + `</td>`
      alignments_table += `<td><a class="btn border-0 m-0 p-0 align-baseline genome">` + value[0] + `</a></td>`
      alignments_table += `<td>` + colored + `</td>`
      alignments_table += `</tr>`
    }

    alignments_table += `</tbody></table>\n\n`;

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

function get_gene_cluster_functions_table(){

}

function show_node_details_modal(e, data){
  var id = e.id;
  var element = document.getElementById(id);

  if (element.getAttribute('class') == 'group') {
      // this is a group of gene clusters
  } else {
      // this is a single gene cluster
    $('#InfoModal').modal('show');
  }
}

// //ANCHOR - Information for the GC
function nodeinfo(e, data) {
  var id = e.id;
  var element = document.getElementById(id);

  if (element.getAttribute('class') == 'group') {
    gene_cluster_context = [];
    for (item in data['infos']['groups'][id]){
      gene_cluster_context.push(data['infos']['groups'][id][item])
    }
    gene_cluster_id = gene_cluster_context[0]
  } else {
    gene_cluster_id = id;
    gene_cluster_context = null;
  }

  $('#InfoModalBody').empty()
  var bodyinfo = $('<div class="card-body overflow-scroll"></div>')
  $('#InfoModalBody').append(bodyinfo)

  bodyinfo.append(get_gene_cluster_display_tables(gene_cluster_id, gene_cluster_context, data))

  $('#AlignmentModalBody').empty()
  var bodyalign = $('<div class="card-body overflow-scroll"></div>')
  $('#AlignmentModalBody').append(
    bodyalign
  )

  $('#GenomeModalBody').empty()

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

  var edgecoloring = {}

  $("#genomecolors :input[type='color']").each((index, element) => {
    edgecoloring[element.id] = [index, element.value]
  })

  var node_size = parseInt($('#size')[0].value);
  var node_thickness = parseInt($('#circ')[0].value);
  var node_distance_x = parseInt($('#distx')[0].value) + node_size + node_thickness;
  var node_distance_y = parseInt($('#disty')[0].value) + node_size;
  var global_x = data["infos"]["meta"]["global_x"] -1;
  var genome_size = Object.keys(data["infos"]["genomes"]).length;
  var global_y = data["infos"]["meta"]["global_y"];
  var group_color = $('#groups')[0].value;
  var node_color = $('#nodes')[0].value;

  var startstop_color = $('#startstop')[0].value;
  var theta = 270 / (global_x)
  
  var radius = 0.5 * (node_distance_x / Math.sin(deg2rad(theta * (1/2))))  
  var start_offset = 0
  
  // Link to Frontend
  var marker_width = 200
  var marker_start = start_offset / node_distance_y
  var marker_stop = marker_start + marker_width / node_distance_y

  // Link to Frontend
  var arrow_width = 50
  var arrow_start = marker_stop
  var arrow_stop = arrow_start + arrow_width / node_distance_y

  var sum_layer = arrow_width + marker_width
  var current_stop = arrow_stop

  var layers = data['infos']['layers_names']
  var layers_sizes = new Object();

  for (var layer in layers) {
    
    // Link to Frontend
    var layer_width = 50
    var layer_start = current_stop
    var layer_stop = layer_start + layer_width / node_distance_y

    current_stop = layer_stop
    sum_layer += layer_width

    var layer_name = data['infos']['layers_names'][layer]
    layers_sizes[layer_name] = [layer_width, layer_start, layer_stop]

  }

  // console.log(arrow_start, arrow_stop, marker_start, marker_stop, entropie_start, entropie_stop)
  console.log(layers_sizes)

  if (radius < sum_layer + start_offset + node_distance_y) {
    radius = sum_layer + start_offset + node_distance_y
  }

  graph_start = radius / node_distance_y

  var size = (global_y * node_distance_y + radius);
  var svg = $('<svg id="result" width="100%" height="100%" version="1.1" viewBox="-' + size + ' -' + size + ' ' + size*2 + ' ' + size*2 + '" position="absolute" xmlns="http://www.w3.org/2000/svg">');

  // var arrow_start = 0
  // var arrow_stop = global_x
  var pointer_length = global_x / 50
  var pointer_height = arrow_stop - arrow_start
  var arrow_thickness = pointer_height / 4
  var steps = 10

  var [circle_c_x, circle_c_y] = transform(global_x - pointer_length, (arrow_start + arrow_thickness) * node_distance_y, theta)
  var [circle_a_x, circle_a_y] = transform(0, (arrow_start + arrow_thickness) * node_distance_y, theta)
  var [circle_b_x, circle_b_y] = transform(0, (arrow_stop - arrow_thickness) * node_distance_y, theta)
  var [circle_d_x, circle_d_y] = transform(global_x - pointer_length, (arrow_stop - arrow_thickness) * node_distance_y, theta)
  var [circle_f_x, circle_f_y] = transform(global_x - pointer_length, (arrow_stop) * node_distance_y, theta)
  var [circle_g_x, circle_g_y] = transform(global_x, (arrow_start + arrow_thickness * 2) * node_distance_y, theta)
  var [circle_e_x, circle_e_y] = transform(global_x - pointer_length, (arrow_start) * node_distance_y, theta)

  if ((global_x) * theta >= 180) {
    var arc_flag = 1
  } else {
    var arc_flag = 0
  }

  svg.append(
    $('<path d="M ' + circle_c_x + ' ' + circle_c_y +
    ' A ' + ((arrow_start + arrow_thickness) * node_distance_y) + ' ' + ((arrow_start + arrow_thickness) * node_distance_y) + ' 0 ' + arc_flag + ' 1 ' + circle_a_x + ' ' + circle_a_y +
    ' L ' + circle_b_x + ' ' + circle_b_y +
    ' A ' + ((arrow_stop - arrow_thickness) * node_distance_y) + ' ' + ((arrow_stop - arrow_thickness) * node_distance_y) + ' 0 ' + arc_flag + ' 0 ' + circle_d_x + ' ' + circle_d_y +
    ' L ' + circle_f_x + ' ' + circle_f_y +
    ' L ' + circle_g_x + ' ' + circle_g_y +
    ' L ' + circle_e_x + ' ' + circle_e_y + 
    ' Z" stroke-width="0" fill="slateGrey"></path>')
  )

  var [circle_h_x, circle_h_y] = transform(0, (arrow_start + arrow_thickness * 2) * node_distance_y, theta)
  svg.append(
    $('<text text-anchor="end" transform="translate (-10)" dominant-baseline="middle" x="' + circle_h_x + '" y="' + circle_h_y + '" dy="0" font-family="sans-serif" fill="black">orientation</text>')
  )

  var l = 1
  while (l < global_x) {

    if (l+steps <= global_x){
      var k = steps;
    } else {
      var k = steps - ((l+steps) - global_x);
    }
    var [circle_l_x, circle_l_y] = transform(l-0.5, (arrow_start + arrow_thickness * 2) * node_distance_y, theta)
    var rotate = theta * (l-0.5)
    if (rotate >= 90 && rotate <= 180) {
      rotate += 180;
    } else if (rotate >= 180 && rotate <= 270) {
      rotate -= 180;
    }
    svg.append(
      $('<text text-anchor="middle" dominant-baseline="middle" transform="rotate(-' + rotate + ' ' + circle_l_x + ' ' + circle_l_y +')" x="' + circle_l_x + '" y="' + circle_l_y + '" dy="0" font-family="sans-serif" fill="white">' + l + '</text>')
    )
    l += k
  };

  var heatmap_max = 1
  for (var key = 0; key < global_x; key++ ) {

    // if (key == 1) {
    //   var add_start = 0
    //   var add_stop = 0.5
    // } else if (key == global_x-1) {
    //   var add_start = 0.5
    //   var add_stop = 0
    // } else {
    //   var add_start = 0.5
    //   var add_stop = 0.5
    // }

    var add_start = 0
    var add_stop = 1

    var [e_x, e_y] = transform(parseInt(key)-add_start, (marker_start * node_distance_y), theta)
    var [f_x, f_y] = transform(parseInt(key)+add_stop, (marker_start * node_distance_y), theta)
    var [g_x, g_y] = transform(parseInt(key)-add_start, (marker_stop * node_distance_y), theta)
    var [h_x, h_y] = transform(parseInt(key)+add_stop, (marker_stop * node_distance_y), theta)

    svg.append(
      $('<path class="marker" id="' + key + '" d="' +
      'M ' + e_x + ' ' + e_y + ' ' +
      'A ' + (marker_start * node_distance_y)  + ' ' + (marker_start * node_distance_y) + ' 0 0 0 ' + f_x + ' ' + f_y + ' ' +
      'L' + h_x + ' ' + h_y +  ' ' +
      'A ' + (marker_stop * node_distance_y)  + ' ' + (marker_stop * node_distance_y) + ' 0 0 1 ' + g_x + ' ' + g_y + ' ' +
      'L' + e_x + ' ' + e_y +
      '" fill="white" stroke="" stroke-width="0"/>')
    )

    for (var layer in layers) {

      var layer_name = data['infos']['layers_names'][layer]
      var [layer_width, layer_start, layer_stop] = layers_sizes[layer_name]
      
      console.log(layer_name, layer_width, layer_start, layer_stop)
      // console.log(layer_name, ll)

      var [a_x, a_y] = transform(parseInt(key)-add_start, (layer_start * node_distance_y), theta)
      var [b_x, b_y] = transform(parseInt(key)+add_stop, (layer_start * node_distance_y), theta)
      var [c_x, c_y] = transform(parseInt(key)-add_start, (layer_stop * node_distance_y), theta)
      var [d_x, d_y] = transform(parseInt(key)+add_stop, (layer_stop * node_distance_y), theta)

      svg.append(
        // $('<path class="entropy" xpos="' + key + '" name="' + (mean_entropy[key] / max).toFixed(3) + '" d="' +
        $('<path class="' + layer_name + '" xpos="' + key + '" name="' + (0 / heatmap_max).toFixed(3) + '" d="' +
        'M ' + a_x + ' ' + a_y + ' ' +
        'A ' + (layer_start * node_distance_y)  + ' ' + (layer_start * node_distance_y) + ' 0 0 0 ' + b_x + ' ' + b_y + ' ' +
        'L' + d_x + ' ' + d_y +  ' ' +
        'A ' + (layer_stop * node_distance_y)  + ' ' + (layer_stop * node_distance_y) + ' 0 0 1 ' + c_x + ' ' + c_y + ' ' +
        'L' + a_x + ' ' + a_y +
        // '" fill="' + lighter_color('#00ff00', '#ff0000', mean_entropy[key] / max) + '" stroke="" stroke-width="2"/>')
        '" fill="' + lighter_color('#00ff00', '#ff0000', 0 / heatmap_max) + '" stroke="" stroke-width="0"/>')
      )

      var [circle_h_x, circle_h_y] = transform(0, (layer_start + (layer_stop-layer_start)/2) * node_distance_y, theta)
      svg.append(
        $('<text text-anchor="end" transform="translate (-10)" dominant-baseline="middle" x="' + circle_h_x + '" y="' + circle_h_y + '" dy="0" font-family="sans-serif" fill="black">' + layer_name + '</text>')
      )
    }

  };
  
  var groups = data['infos']['groups']
  var edges = data['elements']['edges']
  var nodes = data['elements']['nodes']

  for(var i in edges) {

    var edge = data['elements']['edges'][i];
    var genome = Object.keys(edge['genome']).length;

    if ($('#coloredge').prop('checked') == true) {
      var color = pickcolor (edgecoloring, Object.keys(edge['genome']))
      var draw = lighter_color('#ffffff', color, genome / genome_size);
    } else {
      var draw = lighter_color('#ffffff', '#000000', genome / genome_size);
    }

    if (edge['direction'] == 'L') {
      var stroke = ' stroke-dasharray="5,5" '
    } else if (edge['direction'] == 'B') {
      var stroke = ' stroke-dasharray="15,5" '
    } else {
      var stroke = ''
    }

    var source = edge['source']
    var target = edge['target']

    if (source !=  'start' && target != 'stop' && edge['shown'] == 1){

      var i_x = nodes[source]['position']['x']
      var i_y = nodes[source]['position']['y']
      var j_x = nodes[target]['position']['x']
      var j_y = nodes[target]['position']['y']

      var [circle_i_x, circle_i_y] = transform(i_x-0.5, (i_y + graph_start) * node_distance_y, theta);
      var [circle_j_x, circle_j_y] = transform(j_x-0.5, (j_y + graph_start) * node_distance_y, theta);

      if (edge['bended'] == ""){

        svg.append(
          $('<path class="path" d="M ' + circle_i_x + ' ' + circle_i_y + ' A ' + ((i_y + graph_start) * node_distance_y)  + ' ' + ((j_y + graph_start) * node_distance_y) + ' 0 0 0 ' + circle_j_x + ' ' + circle_j_y + '"' + stroke + ' stroke="' + draw + '" stroke-width="2" fill="none"/>')
        )

      } else {

        var bended_edge = '<path class="path" d="M ' + circle_i_x + ' ' + circle_i_y
        var o_y = i_y

        for(var n in edge['bended']) {

          var n_x = edge['bended'][n]['x']
          var n_y = edge['bended'][n]['y']

          var [circle_n_x, circle_n_y] = transform(n_x-0.5, (n_y + graph_start) * node_distance_y, theta);

          bended_edge += 'A ' + ((o_y + graph_start) * node_distance_y)  + ' ' + ((n_y + graph_start) * node_distance_y) + ' 0 0 0 ' + circle_n_x + ' ' + circle_n_y
          var o_y = n_y
        }

        bended_edge += 'A ' + ((o_y + graph_start) * node_distance_y)  + ' ' + ((j_y + graph_start) * node_distance_y) + ' 0 0 0 ' + circle_j_x + ' ' + circle_j_y + '"' + stroke + ' stroke="' + draw + '" stroke-width="2" fill="none"/>'

        svg.append(
          $(bended_edge)
        )
      }
    }
  };

  var group_nodes = [];
  for(var g in groups) {
    var group = data['infos']['groups'][g]
    group_nodes = group_nodes.concat(group)
  };

  // var shannon = {}
  for(var k in nodes) {

    if (k !=  'start' && k != 'stop') {

      var node = data['elements']['nodes'][k];
      var k_x = node['position']['x']
      var k_y = node['position']['y']
      var genome = Object.keys(node['genome']).length;
      var name = node['name']

      // if (pos_x in shannon){
      //   if (name in shannon[pos_x]) {
      //     shannon[pos_x][name] += genome
      //   } else {
      //     shannon[pos_x][name] = genome
      //   }
      // } else {
      //   shannon[pos_x] = {}
      //   shannon[pos_x][name] = genome
      // }

      if (!group_nodes.includes(k)) {
        var node_class = 'class="node'
      } else {
        var node_class = 'stroke-opacity="0" fill-opacity="0" class="pseudo'
      }

      if ($('#colorgenecluster').prop('checked') == true) {
        var color = pickcolor (edgecoloring, Object.keys(node['genome']))
        var draw = lighter_color('#ffffff', color, genome / genome_size);
      } else {
        var draw = lighter_color('#ffffff', '#000000', genome / genome_size);
      }

      var [circle_k_x, circle_k_y] = transform(k_x-0.5, (k_y + graph_start) * node_distance_y, theta);

      if (k ==  'stop' || k == 'start'){
        svg.append(
          $('<circle ' + node_class + '" id="' + k + '" cx="' + circle_k_x + '" cy="' + circle_k_y + '" r="' + node_size + '" fill="' + lighter_color('#ffffff', startstop_color, genome / genome_size) + '" stroke="' + draw + '" stroke-width="' + node_thickness + '"/>')
        )
      } else {
        svg.append(
          $('<circle ' + node_class + '" id="' + k + '" cx="' + circle_k_x + '" cy="' + circle_k_y + '" r="' + node_size + '" fill="' + lighter_color('#ffffff', node_color, genome / genome_size) + '" stroke="' + draw + '" stroke-width="' + node_thickness + '"/>')
        )
      }
    }
  };

  for(var l in groups) {
    var group = data['infos']['groups'][l]

    var group_length = group.length
    var left_node_name = group[0]
    var right_node_name = group[group_length-1]

    var left_node = data['elements']['nodes'][left_node_name];
    var right_node = data['elements']['nodes'][right_node_name];

    var genome = Object.keys(left_node['genome']).length;

    if ($('#colorgenecluster').prop('checked') == true) {
      var color = pickcolor (edgecoloring, Object.keys(left_node['genome']))
      var draw = lighter_color('#ffffff', color, genome / genome_size);
    } else {
      var draw = lighter_color('#ffffff', '#000000', genome / genome_size);
    }

    var l_x = left_node['position']['x']
    var l_y = left_node['position']['y']

    var m_x = right_node['position']['x']
    var m_y = right_node['position']['y']

    var group_offset = node_size / node_distance_y

    var [circle_t_x, circle_t_y] = transform(l_x-0.5, (l_y + graph_start + group_offset) * node_distance_y, theta);
    var [circle_u_x, circle_u_y] = transform(l_x-0.5, (l_y + graph_start - group_offset) * node_distance_y, theta);
    var [circle_v_x, circle_v_y] = transform(m_x-0.5, (m_y + graph_start + group_offset) * node_distance_y, theta);
    var [circle_w_x, circle_w_y] = transform(m_x-0.5, (m_y + graph_start - group_offset) * node_distance_y, theta);

    if ((m_x - l_x) * theta >= 180) {
      var arc_flag = 1
    } else {
      var arc_flag = 0
    }

    svg.append(
      $('<path class="group" id="' + l + '" d="' +
      'M ' + circle_t_x + ' ' + circle_t_y + ' ' +
      'A ' + ((l_y + graph_start + group_offset) * node_distance_y) + ' ' + ((m_y + graph_start + group_offset) * node_distance_y) + ' 0 ' + arc_flag + ' 0 ' + circle_v_x + ' ' + circle_v_y + ' ' +
      'A ' + node_size + ' ' + node_size + ' 0 0 0 ' + circle_w_x + ' ' + circle_w_y + ' ' +
      'A ' + ((l_y + graph_start - group_offset) * node_distance_y) + ' ' + ((m_y + graph_start - group_offset) * node_distance_y) + ' 0 ' + arc_flag + ' 1 ' + circle_u_x + ' ' + circle_u_y + ' ' +
      'A ' + node_size + ' ' + node_size + ' 0 0 0 ' + circle_t_x + ' ' + circle_t_y +
      '" fill="' + lighter_color('#ffffff', group_color, genome / genome_size) + '" stroke="' + draw + '" stroke-width="' + node_thickness + '"/>')
    )
  };

  body.append(svg)

  body.html(body.html());
}

//ANCHOR - Check node
function checknode(searchpos, positions, node, mingenomes, maxgenomes, minposition, maxposition, searchfunction, expressiondrop, expressioncomparison) {

  var append = true

  if (searchpos == true) {

    if (!positions.includes(node['position']['x'].toString())) {
      append = false
    }
  }

  if (mingenomes != '-1'){
    if (eval(Object.keys(node['genome']).length + '<=' + mingenomes)){
      append = false
    }
  }

  if (maxgenomes != '-1'){
    if (eval(Object.keys(node['genome']).length + '>=' + maxgenomes)){
      append = false
    }
  }

  if (minposition != '-1'){
    if (eval(node['position']['x'] + '<=' + minposition)){
      append = false
    }
  }

  if (maxposition != '-1'){
    if (eval(node['position']['x'] + '>=' + maxposition)){
      append = false
    }
  }

  var d = get_gene_cluster_consensus_functions(node['genome'])
  for (var source of Object.keys(searchfunction)) {
    if (!(d[source][2].includes(searchfunction[source]))) {
      append = false
    }
  }

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

  // $.getJSON('static/json/result.json?' + new Date().getTime(), function(data) {

    // $('#redraw').on('click', function() {

    //   // [...document.querySelectorAll('*')].forEach(node => {
    //   //   if (node._tippy) {
    //   //     node._tippy.destroy();
    //   //   }
    //   // });
    //   $(document).off().find("*").off();
    //   main()

    // })


    // This title must update with Generic data from JSON
    $('#title-panel-first-line').text('TITLE');
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
          success: function(data){
            console.log(data)
          }
      });

      // [...document.querySelectorAll('*')].forEach(node => {
      //   if (node._tippy) {
      //     node._tippy.destroy();
      //   }
      // });

      $(document).off().find("*").off();
      main()
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
            // $('<div class="row gy-0 align-items-center">').append(
              $('<div class="col-2 d-flex align-items-center">').append(
                $('<div class="form-switch">').append(
                  $('<input class="form-check-input" type="checkbox" id="flex' + genome + '" name="' + genome + '" aria-label="..." data-bs-toggle="tooltip" data-bs-placement="top" title="Tooltip on top"' + state + '>')
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
          $('<div class="col-8">').append(
            genome
          )
        ).append(
          $('<div class="col-4 text-end" id="number_' + genome + '">').append(
            0
          )
        )

      }
    }

    var body = $('#svgbox')
    body.empty()

    // learn about the functional annotation sources
    functional_annotation_sources_available = data['infos']['functional_annotation_sources_available'];

    generate_svg(body, data);

    window.zoomSVG = svgPanZoom('#result', {
      zoomEnabled: true,
      panEnabled: false,
      controlIconsEnabled: false,
      minZoom: 0.1,
      maxZoom: 100
    });

    //ANCHOR - Zoom pan functions
    $('#plus').on('click', function() {
      window.zoomSVG.zoomIn();
    })

    $('#minus').on('click', function() {
      window.zoomSVG.zoomOut();
    })

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

    //ANCHOR - Choose genecall
    $(document).on("click", ".genome", function() {

      var id = $('#name').attr('name')
      var call = this.innerText;
      var genome = $(this).parent().parent().children(":first").text();

      var name = data['elements']['nodes'][id]['name'];
      var length = data['elements']['nodes'][id]['genome'][genome]['length'];
      var direction = data['elements']['nodes'][id]['genome'][genome]['direction'];
      var paralog = data['elements']['nodes'][id]['genome'][genome]['max_paralog'];
      var partial = data['elements']['nodes'][id]['genome'][genome]['partial'];

      var info = data['elements']['nodes'][id]['genome'][genome];

      $('#GenomeModalBody').empty()
      var bodygenome = $('<div class="card-body overflow-scroll"></div>')
      $('#GenomeModalBody').append(
        bodygenome
      )
      appendgenome(bodygenome, name, call, genome, length, direction, paralog, partial, info)

      $('#AlignmentModal').modal('hide');
      $('#GenomeModal').modal('show');
    });

    $(document).on("click", ".choice li a", function() {
      var name = $(this).attr('name');
      var dropitem = $(this).parent().parent().parent().children(":first");

      dropitem[0].name = name
      dropitem.empty()
      dropitem.append(
        $('<span class="caret"></span>').append(
          name
        )
      )
    })

    //ANCHOR - Bin dropdown choice function
    $(document).on("click", ".binremove", function() {

      var id = $(this).attr('name');
      var binid = $(this).attr('bin');

      var name = document.getElementById(id);
      bins = marknode(name, data, binid, bins);

      $(this).parent().parent().parent().remove();

    })

    //ANCHOR - Change GC in group window
    $(document).on("click", ".gcchoice li a", function() {

      var id = $(this).attr('name');
      var drop = $('#drop');
      var group = drop.attr('name');
      var name = data['elements']['nodes'][id]['name']

      var dropitem = $('#name');
      dropitem[0].name = id;
      dropitem.empty();
      dropitem.append(
        $('<span class="caret"></span>').append(name)
      );

      var position = data['elements']['nodes'][id]['position']['x'] + " / " + (data["infos"]["meta"]["global_x"] - 1);
      var genomes = Object.keys(data['elements']['nodes'][id]['genome']).length + " / " + (data['infos']['num_genomes']);
      var gene_cluster_data = data['elements']['nodes'][id]['genome'];

      $('#InfoModalBody').empty();
      var bodyinfo = $('<div class="card-body overflow-scroll"></div>');
      $('#InfoModalBody').append(bodyinfo);

      basic_info = {'Name': id, 'Genomes': genomes, 'Position': position};
      bodyinfo.append(get_gene_cluster_display_tables('', basic_info, gene_cluster_data));

      var alignment = {}

      if (id != 'start' && id != 'stop') {
        for (var genome of Object.keys(data['elements']['nodes'][id]['genome'])) {
          alignment[genome] = [data['elements']['nodes'][id]['genome'][genome]['gene_call'], data['elements']['nodes'][id]['name']];
        }
      }

      $('#AlignmentModalBody').empty()
      var bodyalign = $('<div class="card-body overflow-scroll"></div>');
      $('#AlignmentModalBody').append(bodyalign);
      appendalignment(bodyalign, alignment)
    });

    //ANCHOR - Bin dropdown choice function
    $(document).on("click", ".binchoice li a", function() {

      var id = $(this).attr('name');
      var drop = $(this).parent().parent().parent();
      var group = drop.attr('name');
      var name = data['elements']['nodes'][id]['name'];
  
      var dropitem = $('#' + group + 'name');
      dropitem[0].name = id;
      dropitem.empty();
      dropitem.append(
          $('<span class="caret"></span>').append(name)
      );
  
      var position = data['elements']['nodes'][id]['position']['x'] + "/" + (data["infos"]["meta"]["global_x"] - 1);
      var genomes = Object.keys(data['elements']['nodes'][id]['genome']).length + "/" + (data['infos']['genomes'].length);
      var gene_cluster_data = data['elements']['nodes'][id]['genome'];
  
      var body = $('#' + group + 'div');
      body.empty();
  
      var basic_info = {'Name': id, 'Genomes': genomes, 'Position': position};
      body.append(get_gene_cluster_display_tables(group, basic_info, gene_cluster_data));
  });
  

    //ANCHOR - Bin creation
    var bins = {"bin1": []};
    var binnum = 1;
    
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
    $('#bininfo').on('click', function() {

      var selection = document.querySelector('input[name="binradio"]:checked');

      // if (selection !== null) {
      var binid = selection.value;
      var appendlist = [];

      $('#BinModalBody').empty();
      for (var id of bins[binid]) {
        var element = document.getElementById(id);
        if (element.getAttribute('class') == 'group') {

          var drop = $('<div class="dropdown" id="' + id + 'drop" name="' + id + '"></div>');
          drop.append(
            $('<a class= "btn border-1 m-0 p-0 dropdown-toggle" id="' + id + 'name" name="Choose GC" data-toggle="dropdown" href="#"></a>').append(
              $('<span class="caret"></span>').append('Choose GC')
            )
          );

          var dropitem = $('<ul class="dropdown-menu pre-scrollable binchoice"></ul>')
          var grouplist = data['infos']['groups'][id]

          for (var listitem of grouplist) {
            dropitem.append(
              $('<li class="caret"></li>').append(
                $('<a class="dropdown-item" name="' + listitem + '" href="#">' + data['elements']['nodes'][listitem]['name'] + '</a>')
              )
            );
          }
          drop.append(dropitem)

          var position = '';
          var genomes = '';
          var group = id;
          var gene_cluster_data = '';

          appendlist.push([id, drop, position, genomes, group, gene_cluster_data]);
        } else {

          var drop = $('<div id="' + id + 'name" name="' + id + '"></div>').append(
            data['elements']['nodes'][id]['name']
          )
          var position = data['elements']['nodes'][id]['position']['x'] + "/" + (data["infos"]["meta"]["global_x"] - 1);
          var genomes = Object.keys(data['elements']['nodes'][id]['genome']).length + "/" + (data['infos']['genomes'].length);
          var group = 'None'
          var gene_cluster_data = data['elements']['nodes'][id]['genome']
          appendlist.push([id, drop, position, genomes, group, gene_cluster_data])
        }
      }

      for (var [id, drop, position, genomes, group, gene_cluster_data] of appendlist) {

        var body = $('<div class="card-body overflow-scroll" id="' + id + 'div"></div>')

        $('#BinModalBody').append(
          $('<div class="card border mb-3 w-100" style="height:400px;"></div>').append(
            $('<div class="card-header"></div>').append(
              $('<div class="d-flex justify-content-end align-items-center"></div>').append(
                $('<button class="btn-close binremove" name="' + id + '" bin="' + binid + '"></button>')
              )
            )
          ).append(
            body
          )
        )

        basic_info = {'Gene Cluster': id, 'Genomes': genomes, 'Position': position}
        body.append(get_gene_cluster_display_tables(id, basic_info, gene_cluster_data))
      }

      $('#BinModal').modal('show');
      // }
    })

    $('#InfoDownload').on('click', async function() {

      var id = $('#name').attr('name');
      var name = $('#name').text();

      if (id !=  'Choose GC' && id != 'start' && id != 'stop'){

        var group = $('#group').text();
        var genomes = $('#genomes').text();
        var position = $('#position').text();

        var csv = "Name\t" + name + "\nGroup\t" + group + "\nGenomes\t" + genomes + "\nPosition\t" + position + "\nSource\tAccession\tFunction\tConfidence";

        var func = get_gene_cluster_consensus_functions(data['elements']['nodes'][id]['genome']);

        for (var [key, value] of Object.entries(func)) {
          csv += "\n" + key + "\t" + value[0] + "\t" + value[1] + "\t" + value[2];
        };

        var blob = new Blob([csv]);
        downloadBlob(blob, name + ".csv");
      }
    });

    $('#GenomeDownload').on('click', async function() {


      var id = $('#name').attr('name')
      var name = $('#name').text();

      if (id !=  'Choose GC' && id != 'start' && id != 'stop'){

        var genome = $('#genome').text();
        var genecall = $('#genecall').text();
        var length = $('#length').text();
        var partial = $('#partial').text();
        var paralog = $('#paralog').text();
        var direction = $('#direction').text();

        var csv = "Genome\t" + genome + "\nGenecall\t" + genecall + "\nLength\t" + length + "\Partial\t" + partial + "\ParalogS\t" + paralog + "\Direction\t" + direction + "\nSource\tAccession\tFunction\te-Value";

        var func = fetchgenome(data['elements']['nodes'][id]['genome'][genome]);

        for (var [key, value] of Object.entries(func)) {
          csv += "\n" + key + "\t" + value[0] + "\t" + value[1] + "\t" + value[2];
        };

        var blob = new Blob([csv]);
        downloadBlob(blob, name + ".csv");
      }
    });

    $('#AlignmentDownload').on('click', async function() {

      var id = $('#name').attr('name')
      var name = $('#name').text();

      if (id !=  'Choose GC' && id != 'start' && id != 'stop'){

        var al = {}
        for (var genome of Object.keys(data['elements']['nodes'][id]['genome'])) {
          al[genome] = [data['elements']['nodes'][id]['genome'][genome]['gene_call'], data['elements']['nodes'][id]['name']]
        }

        var align = await fetchalignment(al);
        var csv = "";

        for (var [genome, value] of Object.entries(align)) {

          csv += ">" + name + "|Genome:" + genome +"|Genecall:" + value[0] + "\n";
          csv += value[1].match(/.{1,60}/g).join("\r\n") + "\n";

        }

        var blob = new Blob([csv]);
        downloadBlob(blob, name + ".fa");
      }
    });

    $('#jsonDownload').on('click', async function() {
      var csv = JSON.stringify(data);
      var blob = new Blob([csv]);
      downloadBlob(blob,  pass_project_name + ".json");
    });

    $('#svgDownload').on('click', async function() {
      var blob = new Blob([$('#svgbox')[0].innerHTML]);
      downloadBlob(blob, pass_project_name + ".svg");
    });

    $('#expressiontext, #searchFunctionsValue').on('input', function(){
      // Check the length of #expressiontext and #searchFunctionsValue values
      if ($('#expressiontext').val().trim() !== '' || $('#searchFunctionsValue').val().trim() !== '') {
          // Enable the search button if at least one of them is not empty
          $('#search').prop('disabled', false);
      } else {
          // Disable the search button if both of them are empty
          $('#search').prop('disabled', true);
      }
    });


    $('#searchadd').on('click', function() {

      var selection = document.querySelector('input[name="binradio"]:checked')
      var binid = selection.value

      for (var [id, members] of Object.entries(searched)) {
        if (!(id in bins[binid])) {

          var e = document.getElementById(id);
          bins = marknode(e, data, binid, bins);

        }
      }

    })

    $('#searchremove').on('click', function() {

      var selection = document.querySelector('input[name="binradio"]:checked')
      var binid = selection.value

      for (var [id, members] of Object.entries(searched)) {
        if (id in bins[binid]) {

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

      searchEraseButton.disabled = false;
      searchColorButton.disabled = false;
      searchAppendBin.disabled = false;
      searchRemoveBin.disabled = false;

      var mingenomes = ($("#mingenomes").prop('checked') == true && !isNaN($("#mingenomestext")[0].value)) ? $("#mingenomestext")[0].value : '-1'
      var maxgenomes = ($("#maxgenomes").prop('checked') == true && !isNaN($("#maxgenomestext")[0].value)) ? $("#maxgenomestext")[0].value : '-1'
      var minentropy = ($("#minentropy").prop('checked') == true && !isNaN($("#minentropytext")[0].value)) ? $("#minentropytext")[0].value : '-1'
      var maxentropy = ($("#maxentropy").prop('checked') == true && !isNaN($("#maxentropytext")[0].value)) ? $("#maxentropytext")[0].value : '-1'
      var minposition = ($("#minposition").prop('checked') == true && !isNaN($("#minpositiontext")[0].value)) ? $("#minpositiontext")[0].value : '-1'
      var maxposition = ($("#maxposition").prop('checked') == true && !isNaN($("#maxpositiontext")[0].value)) ? $("#maxpositiontext")[0].value : '-1'
      var searchfunction = {}
      var expressioncomparison = ''

      var expressiondrop = $('#expressiondrop').attr('value')
      var expressionrel = $('#expressionrel').attr('value')
      var expressiontext = $('#expressiontext')[0].value

      if (expressionrel != "Choose operator" && expressiontext != '') {
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

      for (var source of functional_annotation_sources_available) {
        if ($("#" + source).prop('checked') == true) {
          searchfunction[source] = $('#functiontext')[0].value
        }
      }

      var positions = []
      var searchpos = false

      if ((expressiondrop == "Entropy" && expressionrel != "Choose operator") || minentropy != '-1' || maxentropy != '-1') {
        var searchpos = true
        var entropy = document.querySelectorAll(".entropy")
        for (var en of entropy) {

          var append = true
          var value = en.getAttribute("name")

          if (minentropy != '-1'){
            if (eval(value + '<' + minentropy)){
              append = false
            }
          }

          if (maxentropy != '-1'){
            if (eval(value + '>' + maxentropy)){
              append = false
            }
          }

          if ((expressiondrop == "Entropy" && expressionrel != "Choose operator") && expressionrel != '\u{25C2}\u{25AA}\u{25B8}' && expressionrel != '\u{25AA}\u{25B8}' && expressionrel != '\u{25C2}\u{25AA}') {
            if (!eval(value + expressioncomparison)) {
              append = false
            }
          }

          if (append == true) {
            positions.push(en.getAttribute("xpos"))
          }
        }
      }

      // console.log(positions)

      if ((expressioncomparison != "Choose operator") || Object.keys(searchfunction).length != 0 || mingenomes != '-1' || maxgenomes != '-1' || minposition != '-1' || maxposition != '-1' || searchpos == true) {

        var nodes = document.querySelectorAll(".node")
        for (var node of nodes) {

          var id = node.getAttribute("id")
          var node = data['elements']['nodes'][id]

          if (checknode(searchpos, positions, node, mingenomes, maxgenomes, minposition, maxposition, searchfunction, expressiondrop, expressioncomparison) == true) {
            if (expressiondrop == "Name" && expressionrel != '\u{2264}' && expressionrel != '\u{2265}' && expressionrel != '\u{003C}' && expressionrel != '\u{003E}') {
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

            if (checknode(searchpos, positions, node, mingenomes, maxgenomes, minposition, maxposition, searchfunction, expressiondrop, expressioncomparison) == true) {
              if (expressiondrop == "Name" && expressionrel != '\u{2264}' && expressionrel != '\u{2265}' && expressionrel != '\u{003C}' && expressionrel != '\u{003E}') {
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
      }

      // console.log(searched)

      var toastbody = $('#toastbody')
      toastbody.empty()
      toastbody.append(
        'You have ' + Object.keys(searched).length + ' item(s) in your queue.'
      )
      var searchtoast = bootstrap.Toast.getOrCreateInstance($('#searchtoast'))
      searchtoast.show()
    })

    var entropy = document.querySelectorAll(".entropy")
    // for (var en of entropy) {
    //   tippy(en, {
    //     content: '<strong>' + en.getAttribute("name") + '</strong>' + '<br />',
    //     allowHTML: true,
    //     arrow: true,
    //     duration: 0,
    //     followCursor: true,
    //     theme: "light",
    //   });
    // }
    // }

    var nodes = document.querySelectorAll(".node")
    var divs = document.querySelectorAll(".node, .group");
    for (var el of divs) {

      if (el.getAttribute("id").startsWith('GCG_')){
        var id = data["infos"]["groups"][el.getAttribute("id")][0]
        var name = el.getAttribute("id")
      } else {
        var id = el.getAttribute("id")
        var name = data['elements']['nodes'][el.getAttribute("id")]['name']
      }

      // tippy(el, {
      //   content: '<strong>' + name + '</strong>' + '<br />',
      //   allowHTML: true,
      //   onHide(instance) {
      //     if (instance.reference.id.startsWith('GCG_')){
      //       var id = data["infos"]["groups"][instance.reference.id][0]
      //     } else {
      //       var id = instance.reference.id
      //     }
      //     var elements = Object.keys(data['elements']['nodes'][id]['genome'])

      //     for (var element of elements) {
      //       $('#number_' + element)[0].innerText = '0';
      //     }
      //   },
      //   onShow(instance) {
      //     if (instance.reference.id.startsWith('GCG_')){
      //       var id = data["infos"]["groups"][instance.reference.id][0]
      //     } else {
      //       var id = instance.reference.id
      //     }
      //     var elements = Object.keys(data['elements']['nodes'][id]['genome'])

      //     for (var element of elements) {
      //       $('#number_' + element)[0].innerText = '1';
      //     }
      //   },
      //   arrow: true,
      //   duration: 0,
      //   followCursor: true,
      //   theme: "light",
      // });
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
  }})

}

//ANCHOR - Main function after loading DOM
$(document).ready(function() {

  // console.log('start pangrah layout creation!')
  main()

});
