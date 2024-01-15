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
const notation = ["COG20_PATHWAY", "KEGG_Class", "Transfer_RNAs", "KOfam", "KEGG_Module", "COG20_CATEGORY", "COG20_FUNCTION"];

for (var source of notation) {

  $('#functiondiv').append(
    $('<div class="col-2"></div').append(
      $('<input class="form-check-input" type="checkbox" id="' + source + '" value="" data-bs-toggle="tooltip" data-bs-placement="top"></input>')
    )
  ).append(
    $('<div class="col-10"></div').append(
      source
    )
  )
}

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
  for (var source of notation) {

    var id = []
    var func = []

    if (gene_cluster_data != '') {
      for (var element of Object.keys(gene_cluster_data)) {
        var entry = gene_cluster_data[element][source]

        if (entry === 'None' || entry === undefined) {
          entry = ['None', 'None', 'None']
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
    basic_info_table += `<table class="table table-striped table-bordered" id="node_basics_table">`;
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
    functions_table += `<table class="table table-striped" id="node_functions_table">`;
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
  for (var source of notation) {

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
    alignments_table += `<table class="table table-striped" id="node_sequence_alignments_table">`;
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
function transform(x, y, node_distance_y, radius, theta) {
  var circle_x = ((y * node_distance_y) + radius) * Math.sin(deg2rad(theta * x))
  var circle_y = ((y * node_distance_y) + radius) * Math.cos(deg2rad(theta * x))
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
async function generate_svg(body, data) {

  var edgecoloring = {}

  $("#genomecolors :input[type='color']").each((index, element) => {
    edgecoloring[element.id] = [index, element.value]
  })

  var node_size = parseInt($('#size')[0].value);
  var node_thickness = parseInt($('#circ')[0].value);
  var node_distance_x = parseInt($('#distx')[0].value) + node_size + node_thickness;
  var node_distance_y = parseInt($('#disty')[0].value) + node_size;
  var global_x = data["infos"]["meta"]["global_x"];
  var genome_size = Object.keys(data["infos"]["genomes"]).length;
  var global_y = data["infos"]["meta"]["global_y"];
  var group_color = $('#groups')[0].value;
  var node_color = $('#nodes')[0].value;

  var startstop_color = $('#startstop')[0].value;
  var theta = 360 / global_x
  var radius = 0.5 * (node_distance_x / Math.sin(deg2rad(theta * (1/2))))
  if (radius < 1000) {
    radius = 1000
  }
  var size = (global_y * node_distance_y + radius);
  var svg = $('<svg id="result" width="100%" height="100%" version="1.1" viewBox="-' + size + ' -' + size + ' ' + size*2 + ' ' + size*2 + '" position="absolute" xmlns="http://www.w3.org/2000/svg">');

  // svg.append(
  //   $('<text text-anchor="middle" alignment-baseline="central" id="info" x="0" y="0" dy="0" font-family="sans-serif" fill="slateGrey"></text>').append(
  //     $('<tspan x="0" dy="-400" font-size="400">' + pass_project_name + ' PanGraph</tspan>')
  //   ).append(
  //     $('<tspan x="0" dy="800" font-size="250">Preserved nodes: ' + (data['infos']['visualization']['instances'] / data['infos']['original']['instances'] *100).toFixed(3) + '%</tspan>')
  //   )
  // )

  var arrow_start = 0
  var arrow_stop = global_x
  var pointer_length = 2
  var pointer_height = 4
  var arow_thickness = 2
  var arrow_y_offset = 21
  var steps = 50

  var [a_x, a_y] = transform(arrow_start, - arrow_y_offset + arow_thickness / 2, node_distance_y, radius, theta)
  var [b_x, b_y] = transform(arrow_start, - arrow_y_offset - arow_thickness / 2, node_distance_y, radius, theta)
  var [c_x, c_y] = transform(arrow_stop - pointer_length, - arrow_y_offset + arow_thickness / 2, node_distance_y, radius, theta)
  var [d_x, d_y] = transform(arrow_stop - pointer_length, - arrow_y_offset - arow_thickness / 2, node_distance_y, radius, theta)
  var [e_x, e_y] = transform(arrow_stop - pointer_length, - arrow_y_offset + pointer_height / 2, node_distance_y, radius, theta)
  var [f_x, f_y] = transform(arrow_stop - pointer_length, - arrow_y_offset - pointer_height / 2, node_distance_y, radius, theta)
  var [g_x, g_y] = transform(arrow_stop, - arrow_y_offset, node_distance_y, radius, theta)

  svg.append(
    $('<path d="M ' + c_x + ' ' + c_y +
    ' A ' + (((- arrow_y_offset + (arow_thickness / 2)) * node_distance_y) + radius) + ' ' + (((- arrow_y_offset + (arow_thickness / 2)) * node_distance_y) + radius) + ' 0 1 1 ' + a_x + ' ' + a_y +
    ' L ' + b_x + ' ' + b_y +
    ' A ' + (((- arrow_y_offset - (arow_thickness / 2)) * node_distance_y) + radius) + ' ' + (((- arrow_y_offset - (arow_thickness / 2)) * node_distance_y) + radius) + ' 0 1 0 ' + d_x + ' ' + d_y +
    ' L ' + f_x + ' ' + f_y +
    ' L ' + g_x + ' ' + g_y +
    ' L ' + e_x + ' ' + e_y + ' Z" stroke="slateGrey" stroke-width="2" fill="slateGrey"></path>')
  )

  var l = 1
  while (l < global_x) {

    if (l+steps <= global_x){
      var k = steps;
    } else {
      var k = steps - ((l+steps) - global_x);
    }
    var [l_x, l_y] = transform(l, - arrow_y_offset, node_distance_y, radius, theta)
    var rotate = theta * l
    if (rotate >= 90 && rotate <= 180) {
      rotate += 180;
    } else if (rotate >= 180 && rotate <= 270) {
      rotate -= 180;
    }
    svg.append(
      $('<text text-anchor="middle" dominant-baseline="middle" transform="rotate(-' + rotate + ' ' + l_x + ' ' + l_y +')" x="' + l_x + '" y="' + l_y + '" dy="0" font-family="sans-serif" fill="white">' + l + '</text>')
    )
    l += k
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
    } else {
      var stroke = ''
    }

    var source = edge['source']
    var target = edge['target']

    if (edge['bended'] == ""){

      var [circle_i_x, circle_i_y] = transform(nodes[source]['position']['x'], nodes[source]['position']['y'], node_distance_y, radius, theta);
      var [circle_j_x, circle_j_y] = transform(nodes[target]['position']['x'], nodes[target]['position']['y'], node_distance_y, radius, theta);

      svg.append(
        $('<path class="path" d="M ' + circle_i_x + ' ' + circle_i_y + ' L ' + circle_j_x + ' ' + circle_j_y + '"' + stroke + ' stroke="' + draw + '" stroke-width="2" fill="none"/>')
      )

    } else {

      var [circle_i_x, circle_i_y] = transform(nodes[source]['position']['x'], nodes[source]['position']['y'], node_distance_y, radius, theta);
      var [circle_j_x, circle_j_y] = transform(nodes[target]['position']['x'], nodes[target]['position']['y'], node_distance_y, radius, theta);
      var circle_svg = '<path class="path" d="M ' + circle_i_x + ' ' + circle_i_y

      for(var j in edge['bended']) {

        var bend_n_x = edge['bended'][j]['x']
        var bend_n_y = edge['bended'][j]['y']
        var [circle_n_x, circle_n_y] = transform(bend_n_x, bend_n_y, node_distance_y, radius, theta);

        circle_svg += 'L ' + circle_n_x + ' ' + circle_n_y

      }

      circle_svg += 'L ' + circle_j_x + ' ' + circle_j_y + '"' + stroke + ' stroke="' + draw + '" stroke-width="2" fill="none"/>'

      svg.append(
        $(circle_svg)
      )
    }
  };

  var group_nodes = [];
  for(var g in groups) {
    var group = data['infos']['groups'][g]
    group_nodes = group_nodes.concat(group)
  };

  var shannon = {}
  for(var k in nodes) {

    var node = data['elements']['nodes'][k];
    var pos_x = node['position']['x']
    var pos_y = node['position']['y']
    var genome = Object.keys(node['genome']).length;
    var name = node['name']

    if (pos_x in shannon){
      if (name in shannon[pos_x]) {
        shannon[pos_x][name] += genome
      } else {
        shannon[pos_x][name] = genome
      }
    } else {
      shannon[pos_x] = {}
      shannon[pos_x][name] = genome
    }

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

    var [circle_x, circle_y] = transform(pos_x, pos_y, node_distance_y, radius, theta);

    if (k ==  'stop' || k == 'start'){
      svg.append(
        $('<circle ' + node_class + '" id="' + k + '" cx="' + circle_x + '" cy="' + circle_y + '" r="' + node_size + '" fill="' + lighter_color('#ffffff', startstop_color, genome / genome_size) + '" stroke="' + draw + '" stroke-width="2"/>')
      )
    } else {
      svg.append(
        $('<circle ' + node_class + '" id="' + k + '" cx="' + circle_x + '" cy="' + circle_y + '" r="' + node_size + '" fill="' + lighter_color('#ffffff', node_color, genome / genome_size) + '" stroke="' + draw + '" stroke-width="2"/>')
      )
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

    var left_pos_x = left_node['position']['x']
    var left_pos_y = left_node['position']['y']

    var right_pos_x = right_node['position']['x']
    var right_pos_y = right_node['position']['y']

    var [circle_t_x, circle_t_y] = transform(left_pos_x, left_pos_y, node_distance_y, radius + 10, theta);
    var [circle_u_x, circle_u_y] = transform(left_pos_x, left_pos_y, node_distance_y, radius - 10, theta);
    var [circle_v_x, circle_v_y] = transform(right_pos_x, right_pos_y, node_distance_y, radius + 10, theta);
    var [circle_w_x, circle_w_y] = transform(right_pos_x, right_pos_y, node_distance_y, radius - 10, theta);

    svg.append(
      $('<path class="group" id="' + l + '" d="' +
      'M ' + circle_t_x + ' ' + circle_t_y + ' ' +
      'A ' + (left_pos_y * node_distance_y + radius + 10)  + ' ' + (right_pos_y * node_distance_y + radius + 10) + ' 0 0 0 ' + circle_v_x + ' ' + circle_v_y + ' ' +
      'A 10 10 0 0 0 ' + circle_w_x + ' ' + circle_w_y + ' ' +
      'A ' + (right_pos_y * node_distance_y + radius - 10) + ' ' + (left_pos_y * node_distance_y + radius - 10) + ' 0 0 1 ' + circle_u_x + ' ' + circle_u_y + ' ' +
      'A 10 10 0 0 0 ' + circle_t_x + ' ' + circle_t_y +
      '" fill="' + lighter_color('#ffffff', group_color, genome / genome_size) + '" stroke="' + draw + '" stroke-width="2"/>')
    )
  };

  var entropy = {}
  for (var position of Object.keys(shannon)) {

    var res = 0
    var sum = Object.values(shannon[position]).reduce((partialSum, a) => partialSum + a, 0);
    for (var genome of Object.keys(shannon[position])) {
        var ref_freq = shannon[position][genome] / sum

        res += -(ref_freq * Math.log2(ref_freq))
    }

    entropy[position] = res
  }

  var window = 2
  var mean_entropy = {}
  for (let i = 0; i < Object.keys(entropy).length; i += 1) {
    if (i in entropy){
      var values = []
      values.push(entropy[i])

      for (let j = 1; j <= window; j += 1) {
        if (i-j in entropy) {values.push(entropy[i-j])}
        if (i+j in entropy) {values.push(entropy[i+j])}
      }

      mean_entropy[i] = values.reduce((partialSum, a) => partialSum + a, 0) / values.length
    }
  }

  var max = Math.max(...Object.values(mean_entropy))

  for (var key of Object.keys(mean_entropy)) {

    var [a_x, a_y] = transform(parseInt(key)-0.5, -2, node_distance_y, radius, theta)
    var [b_x, b_y] = transform(parseInt(key)+0.5, -2, node_distance_y, radius, theta)
    var [c_x, c_y] = transform(parseInt(key)-0.5, -14, node_distance_y, radius, theta)
    var [d_x, d_y] = transform(parseInt(key)+0.5, -14, node_distance_y, radius, theta)

    svg.append(
      $('<path class="entropy" xpos="' + key + '" name="' + (mean_entropy[key] / max).toFixed(3) + '" d="' +
      'M ' + a_x + ' ' + a_y + ' ' +
      'A ' + (radius - 2)  + ' ' + (radius - 2) + ' 0 0 0 ' + b_x + ' ' + b_y + ' ' +
      'L' + d_x + ' ' + d_y +  ' ' +
      'A ' + (radius - 14)  + ' ' + (radius - 14) + ' 0 0 1 ' + c_x + ' ' + c_y + ' ' +
      'L' + a_x + ' ' + a_y +
      '" fill="' + lighter_color('#00ff00', '#ff0000', mean_entropy[key] / max) + '" stroke="" stroke-width="2"/>')
    )

    var [e_x, e_y] = transform(parseInt(key)-0.5, -25, node_distance_y, radius, theta)
    var [f_x, f_y] = transform(parseInt(key)+0.5, -25, node_distance_y, radius, theta)
    var [g_x, g_y] = transform(parseInt(key)-0.5, -40, node_distance_y, radius, theta)
    var [h_x, h_y] = transform(parseInt(key)+0.5, -40, node_distance_y, radius, theta)

    svg.append(
      $('<path class="marker" id="' + key + '" d="' +
      'M ' + e_x + ' ' + e_y + ' ' +
      'A ' + (radius - 25)  + ' ' + (radius - 25) + ' 0 0 0 ' + f_x + ' ' + f_y + ' ' +
      'L' + h_x + ' ' + h_y +  ' ' +
      'A ' + (radius - 40)  + ' ' + (radius - 40) + ' 0 0 1 ' + g_x + ' ' + g_y + ' ' +
      'L' + e_x + ' ' + e_y +
      '" fill="white" stroke="" stroke-width="2"/>')
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

    $('#redraw').on('click', function() {

      // [...document.querySelectorAll('*')].forEach(node => {
      //   if (node._tippy) {
      //     node._tippy.destroy();
      //   }
      // });

      $(document).off().find("*").off();
      main()

    })

    $('#settings').on('click', function() {

      var data = new Object;

      data['conntr'] = $('#conntr')[0].value
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
          dataType: "json"
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

      $("#conntr")[0].value = -1;
      $("#condtr")[0].value = data['infos']['max_edge_length_filter']
      $("#maxlength")[0].value = data['infos']['gene_cluster_grouping_threshold']

      for (var [genome, value] of Object.entries(data['infos']['genomes'])) {

        var state = ''
        if (value == 'on') {
          state = ' checked'
        }

        $('#genomecolors').append(
          $('<div class="col-12">').append(
            $('<div class="row gy-0 align-items-center">').append(
              $('<div class="col-2">').append(
                $('<div class="form-switch">').append(
                  $('<input class="form-check-input" type="checkbox" id="flex' + genome + '" name="' + genome + '" aria-label="..." data-bs-toggle="tooltip" data-bs-placement="top" title="Tooltip on top"' + state + '>')
                )
              )
            ).append(
              $('<div class="col-7">').append(
                genome
              )
            ).append(
              $('<div class="col-1">').append(
                $('<i class="user-handle bi bi-arrows-expand"></i>')
              )
            ).append(
              $('<div class="d-flex col-2">').append(
                $('<input type="color" class="form-control form-control-color flex-fill p-0 border-0" id="' + genome + '" name="' + genome + '" aria-label="..." data-bs-toggle="tooltip" data-bs-placement="top" title="Choose your color">')
              )
            )
          )
        )

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
    if ($('#conntr')[0].value == -1){
      // $('#customRange1').prop('disabled', true);
      $('#flexconntr').prop('checked', false);
    }

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

    $('#flexconntr').change(function() {
      if ($(this).prop('checked') == true){
        $('#conntr')[0].value = 0;
        // $('#customRange1')[0].value = 0;
        // $('#customRange1').prop('disabled', false);
      } else {
        $('#conntr')[0].value = -1;
        // $('#customRange1')[0].value = 0;
        // $('#customRange1').prop('disabled', true);
      }
    })

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
        $('#customRange4')[0].value = 0;
        $('#customRange4').prop('disabled', false);
      } else {
        $('#groupcompress')[0].value = -1;
        $('#customRange4')[0].value = 0;
        $('#customRange4').prop('disabled', true);
      }
    })

    $(function () {
      $(".grid").sortable({
        tolerance: 'pointer',
        revert: 'invalid',
        handle: ".user-handle",
        forceHelperSize: true
      });
    });

    //ANCHOR - Choose genecall
    $(document).on("click", ".genome", function() {

      var id = $('#name').attr('name')
      var call = this.innerText;
      var genome = this.parentNode.parentNode.firstChild.innerText;

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
      var name = $(this).attr('name')
      var dropitem = $(this).parent().parent().parent().children(":first")

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

      var id = $(this).attr('name')
      var binid = $(this).attr('bin')

      var name = document.getElementById(id);
      bins = marknode(name, data, binid, bins);

      $(this).parent().parent().parent().remove()

    })

    //ANCHOR - Change GC in group window
    $(document).on("click", ".gcchoice li a", function() {

      var id = $(this).attr('name')
      var drop = $('#drop')
      var group = drop.attr('name')
      var name = data['elements']['nodes'][id]['name']

      var dropitem = $('#name')
      dropitem[0].name = id
      dropitem.empty()
      dropitem.append(
        $('<span class="caret"></span>').append(
          name
        )
      )

      var position = data['elements']['nodes'][id]['position']['x'] + " / " + (data["infos"]["meta"]["global_x"] - 1);
      var genomes = Object.keys(data['elements']['nodes'][id]['genome']).length + " / " + (data['infos']['num_genomes']);
      var gene_cluster_data = data['elements']['nodes'][id]['genome'];

      $('#InfoModalBody').empty()
      var bodyinfo = $('<div class="card-body overflow-scroll"></div>')
      $('#InfoModalBody').append(
        bodyinfo
      )

      basic_info = {'Name': id, 'Genomes': genomes, 'Position': position}
      bodyinfo.append(get_gene_cluster_display_tables('', basic_info, gene_cluster_data))

      var alignment = {}

      if (id != 'start' && id != 'stop') {
        for (var genome of Object.keys(data['elements']['nodes'][id]['genome'])) {
          alignment[genome] = [data['elements']['nodes'][id]['genome'][genome]['gene_call'], data['elements']['nodes'][id]['name']]
        }
      }

      $('#AlignmentModalBody').empty()
      var bodyalign = $('<div class="card-body overflow-scroll"></div>')
      $('#AlignmentModalBody').append(
        bodyalign
      )
      appendalignment(bodyalign, alignment)
    });

    //ANCHOR - Bin dropdown choice function
    $(document).on("click", ".binchoice li a", function() {

      var id = $(this).attr('name')
      var drop = $(this).parent().parent().parent()
      var group = drop.attr('name')
      var name = data['elements']['nodes'][id]['name']

      var dropitem = $('#' + group + 'name')
      dropitem[0].name = id
      dropitem.empty()
      dropitem.append(
        $('<span class="caret"></span>').append(
          name
        )
      )

      var position = data['elements']['nodes'][id]['position']['x'] + "/" + (data["infos"]["meta"]["global_x"] - 1);
      var genomes = Object.keys(data['elements']['nodes'][id]['genome']).length + "/" + (data['infos']['genomes'].length);
      var gene_cluster_data = data['elements']['nodes'][id]['genome'];

      var body = $('#' + group + 'div')
      body.empty()

      basic_info = {'Name': id, 'Genomes': genomes, 'Position': position}
      body.append(get_gene_cluster_display_tables(group, basic_info, gene_cluster_data))
    })

    //ANCHOR - Bin creation
    var bins = {"bin1": []}
    var binnum = 1

    //ANCHOR - Add bin
    $('#binadd').on('click', function() {
      binnum += 1

      $('#bingrid').append(
        $('<div class="col-12" id="bin' + binnum + '"></div>').append(
          $('<div class="row gy-1 align-items-center" id="row' + binnum + '"></div>').append(
            $('<div class="col-2"></div>').append(
              $('<input class="form-check-input" type="radio" name="binradio" id="bin' + binnum + 'radio" value="bin' + binnum + '" checked></input>')
            )
          ).append(
            $('<div class="col-6"></div>').append(
              $('<input class="form-control flex-fill p-0 border-0" style= "background-color: #e9ecef;" type="text" id="bin' + binnum + 'text" value="Bin_' + binnum + '"></input>')
            )
          ).append(
            $('<div class="col-2"></div>').append(
              $('<input type="text" class="form-control float-end text-end flex-fill p-0 border-0" id="bin' + binnum + 'value" name="condtr" value=0 readonly>')
            )
          ).append(
            $('<div class="d-flex col-2"></div>').append(
              $('<input class="form-control form-control-color flex-fill p-0 border-0 colorchange" type="color" name="bin' + binnum + '" id="bin' + binnum + 'color" value="#000000"></input>')
            )
          )
        )
      )

      bins['bin' + binnum] = []
    })

    //ANCHOR - Remove bin
    $('#binremove').on('click', function() {

      var selection = document.querySelector('input[name="binradio"]:checked')

      if (selection !== null) {
        var binid = selection.value

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

      var binid = this.name
      var nodes = bins[binid]

      for (var node of nodes) {

        bins[binid] = bins[binid].filter(item => item !== node)
        var name = document.getElementById(node);
        bins = marknode(name, data, binid, bins);

      }

    });

    //ANCHOR - Info bin
    $('#bininfo').on('click', function() {

      var selection = document.querySelector('input[name="binradio"]:checked')

      // if (selection !== null) {
      var binid = selection.value
      var appendlist = []

      $('#BinModalBody').empty()
      for (var id of bins[binid]) {
        var element = document.getElementById(id);
        if (element.getAttribute('class') == 'group') {

          var drop = $('<div class="dropdown" id="' + id + 'drop" name="' + id + '"></div>')

          drop.append(
            $('<a class= "btn border-1 m-0 p-0 dropdown-toggle" id="' + id + 'name" name="Choose GC" data-bs-toggle="dropdown" href="#"></a>').append(
              $('<span class="caret"></span>').append('Choose GC')
            )
          )

          var dropitem = $('<ul class="dropdown-menu pre-scrollable binchoice"></ul>')
          var grouplist = data['infos']['groups'][id]

          for (var listitem of grouplist) {
            dropitem.append(
              $('<li class="caret"></li>').append(
                $('<a class="dropdown-item" name="' + listitem + '" href="#">' + data['elements']['nodes'][listitem]['name'] + '</a>')
              )
            )
          }
          drop.append(dropitem)

          var position = ''
          var genomes = ''
          var group = id
          var gene_cluster_data = ''

          appendlist.push([id, drop, position, genomes, group, gene_cluster_data])
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

      var id = $('#name').attr('name')
      var name = $('#name')[0].innerText;

      if (id !=  'Choose GC' && id != 'start' && id != 'stop'){

        var group = $('#group')[0].innerText;
        var genomes = $('#genomes')[0].innerText;
        var position = $('#position')[0].innerText;

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
      var name = $('#name')[0].innerText;

      if (id !=  'Choose GC' && id != 'start' && id != 'stop'){

        var genome = $('#genome')[0].innerText;
        var genecall = $('#genecall')[0].innerText;
        var length = $('#length')[0].innerText;
        var partial = $('#partial')[0].innerText;
        var paralog = $('#paralog')[0].innerText;
        var direction = $('#direction')[0].innerText;

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
      var name = $('#name')[0].innerText;

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

    var searched = {}
    $('#search').on('click', function() {

      var mingenomes = ($("#mingenomes").prop('checked') == true && !isNaN($("#mingenomestext")[0].value)) ? $("#mingenomestext")[0].value : '-1'
      var maxgenomes = ($("#maxgenomes").prop('checked') == true && !isNaN($("#maxgenomestext")[0].value)) ? $("#maxgenomestext")[0].value : '-1'
      var minentropy = ($("#minentropy").prop('checked') == true && !isNaN($("#minentropytext")[0].value)) ? $("#minentropytext")[0].value : '-1'
      var maxentropy = ($("#maxentropy").prop('checked') == true && !isNaN($("#maxentropytext")[0].value)) ? $("#maxentropytext")[0].value : '-1'
      var minposition = ($("#minposition").prop('checked') == true && !isNaN($("#minpositiontext")[0].value)) ? $("#minpositiontext")[0].value : '-1'
      var maxposition = ($("#maxposition").prop('checked') == true && !isNaN($("#maxpositiontext")[0].value)) ? $("#maxpositiontext")[0].value : '-1'
      var searchfunction = {}
      var expressioncomparison = ''

      var expressiondrop = $('#expressiondrop').attr('name')
      var expressionrel = $('#expressionrel').attr('name')
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

      for (var source of notation) {
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

  console.log('start pangrah layout creation!')
  main()

});
