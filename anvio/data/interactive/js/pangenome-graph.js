class PangenomeGraphUserInterface {
    constructor() {
        this.bin_dict = {'bin_1': []};
        this.bin_group_dict = {};
        this.current_bin_id = 'bin_1';
        this.current_bin_number = 1;
        this.data = null;
        this.settings_dict = {};
        this.state = 'default';
        this.server_offline = false;

        this.selection = null;
        this.panZoomInstance = null;
        this._region_label_size_initialized = false;

        this.isDown = false
        this.diff = 0
        this.cur_xpos = 0
        this.cur_ypos = 0
        this.new_xpos = 0
        this.new_ypos = 0
        this.drawing_zoom = false
        this.zoomBox = {}

        this.svg_cur_xpos = 0
        this.svg_cur_ypos = 0
        this.svg_new_xpos = 0
        this.svg_new_ypos = 0

        this.synteny = {}
        this.layers_min = {};
        this.layers_max = {};
        this.global_x = 0;
        this.global_y = 0;
        this.nodes = [];
        this.edges = [];
        this.layers = {};
        this.genomes = [];
        this.group_dict = {};
        this.functional_annotation_sources_available = [];

        this.add_bin = this.add_bin.bind(this);
        this.delete_bin = this.delete_bin.bind(this);
        this.delete_all_bins = this.delete_all_bins.bind(this);
        this.update_bin = this.update_bin.bind(this);
        this.switch_bin = this.switch_bin.bind(this);

        this.marknode = this.marknode.bind(this);

        this.start_draw = this.start_draw.bind(this);
        this.hide_tippy = this.hide_tippy.bind(this);
        this.show_tippy = this.show_tippy.bind(this);
        this.press_down = this.press_down.bind(this);
        this.press_move = this.press_move.bind(this);
        this.press_up = this.press_up.bind(this);
        this.svg_download = this.svg_download.bind(this);
        this.fit_aspect = this.fit_aspect.bind(this);
        this.nodeinfo = this.nodeinfo.bind(this);
        this.get_gene_cluster_display_tables = this.get_gene_cluster_display_tables.bind(this);
        this.get_gene_cluster_context_table = this.get_gene_cluster_context_table.bind(this);
        this.get_gene_cluster_basics_table = this.get_gene_cluster_basics_table.bind(this);
        this.get_layer_data = this.get_layer_data.bind(this);
        this.get_gene_cluster_functions_table = this.get_gene_cluster_functions_table.bind(this);
        this.get_region_data = this.get_region_data.bind(this)
        this.appendalignment = this.appendalignment.bind(this);
        this.recolor_alignment = this.recolor_alignment.bind(this);
        this.alignment_download = this.alignment_download.bind(this);
        this.info_download = this.info_download.bind(this);
        this.add_info_to_bin = this.add_info_to_bin.bind(this);
        this.flextree_change = this.flextree_change.bind(this);
        this.save_bin = this.save_bin.bind(this);
        this.load_bin = this.load_bin.bind(this);
        this.save_state = this.save_state.bind(this);
        this.load_state = this.load_state.bind(this);
        this.set_UI_settings = this.set_UI_settings.bind(this);
        this.draw_bin_rings = this.draw_bin_rings.bind(this);
        this._render_bin_visuals = this._render_bin_visuals.bind(this);

        this.draw_newick = this.draw_newick.bind(this);
        this.generate_svg = this.generate_svg.bind(this);

        this.initialize_JSON();
    }

    generate_svg(scale=1) {

        var start = new Date().getTime();
        var svg_arrow = [];
        var svg_search = [];
        var svg_backbone = [];
        var svg_text = [];
        var svg_heatmaps = [];
        var svg_edges = [];
        var svg_nodes = [];
        var svg_tree = [];
        var svg_groups = [];
        var svg_region_labels = [];
        var svg_genome_tracks = {}
        for (var genome of this.genomes){
            svg_genome_tracks[genome] = []
        }
        
        var edgecoloring = {}
        $("#genomecolors .pangraph-colorpicker").each((index, element) => {
            edgecoloring[element.id] = [index, $(element).attr('color')]
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

        var start_angle = parseFloat($('#start_angle')[0].value);
        var end_angle = parseFloat($('#end_angle')[0].value);
        var outer_margin = parseFloat($('#outer_margin')[0].value);
        var inner_margin = parseFloat($('#inner_margin')[0].value);
        var node_size = parseFloat($('#size')[0].value);
        var node_thickness = parseFloat($('#circ')[0].value);
        var edge_thickness = parseFloat($('#edge')[0].value);
        var line_thickness = parseFloat($('#line')[0].value);
        var track_line_width = parseFloat($('#track_line_width')[0].value);
        var node_distance_x = parseFloat($('#distx')[0].value);
        var node_distance_y = parseFloat($('#disty')[0].value);
        var num_position = parseFloat($('#num_position')[0].value);
        var tree_length = parseFloat($('#tree_length')[0].value);
        var offset = parseFloat($('#tree_offset')[0].value);
        var tree_thickness = parseFloat($('#tree_thickness')[0].value);
        var text_offset = parseFloat($('#label_offset')[0].value);
        
        var core_color = $('#core_color').attr('color');
        var paralog_color = $('#paralog_color').attr('color');
        var singleton_color = $('#singleton_color').attr('color');
        var accessory_color = $('#accessory_color').attr('color');
        var rearranged_color = $('#rearranged_color').attr('color');
        var trna_color = $('#trna_color').attr('color');
        var layer_color = $('#layer_color').attr('color');
        var back_color = $('#back_color').attr('color');
        var non_back_color = $('#non_back_color').attr('color');
        var genome_size = this.genomes.length
        
        var theta = (end_angle - start_angle) / (this.global_x+1)
        
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
        var order = []
        var max_dist = 0
        var item_order = []

        // Only parse newick tree if it exists (may be empty for identical genomes)
        if (this.data['meta']['newick'] && this.data['meta']['newick'].length > 0 && $('#flextree').prop('checked') == true) {
            order = this.newick_to_order(this.data['meta']['newick']).reverse()
            for (var item of order) {
                var [name, item_start, item_end] = item
                if (name != 'branching') {
                    item_order.push(name)
                }
                if (item_end > max_dist) {
                    max_dist = item_end
                }
            }
        } else {
            // No newick tree - just use genome names in their original order
            // item_order = this.genomes.slice()

            // Comment from Alex: maybe let's use the current order of layers as
            // per color order instead, this opens more options for users if 
            // there is no current tree loaded or drawn.

            var array_order = []
            for (var [key, value] of Object.entries(edgecoloring)) {
                 array_order.push([value[0], key])
            }
            
            var sorted_array_order = array_order.sort(function(a, b) {
              return b[0] - a[0];
            });
            
            item_order = sorted_array_order.map(arr => arr[1]);
        }
        
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
        
        if ($('#flexglobalbackbone').prop('checked') == true){
        //TEST BACKBONE NON BACKBONE LAYER
            var back_width = parseFloat($('#globalbackbone')[0].value);
            
            var layer_width = back_width
            var layer_middle_start = current_middle_stop + inner_margin
            var layer_middle_stop = layer_middle_start + layer_width
            
            current_middle_stop = layer_middle_stop
            sum_middle_layer += layer_width + inner_margin
            
            middle_layers['back_vs_non_back'] = [layer_width, layer_middle_start, layer_middle_stop]
        }
        
        for (var layer_name of this.layers) {
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
            var radius = 0.5 * (node_distance_x / Math.sin(this.deg2rad(theta * (1/2))))
            var circle_dist = sum_middle_layer + graph_size * 0.5
            var extra_offset = 0
            
            sum_middle_layer += extra_offset
            for (var layer in middle_layers) {
                var [layer_width, layer_start, layer_stop] = middle_layers[layer]
                middle_layers[layer] = [layer_width, layer_start + extra_offset, layer_stop + extra_offset]
            }
        
            var y_size = (sum_middle_layer + (this.global_y * node_distance_y) + sum_outer_layer);
            var x_size = (sum_middle_layer + (this.global_y * node_distance_y) + sum_outer_layer);
            if (scale == 0){
                var svg_core = $('<svg id="result" width="' + x_size*2 + 'px" height="' + y_size*2 + 'px" version="1.1" viewBox="-' + x_size + ' -' + y_size + ' ' + x_size*2 + ' ' + y_size*2 + '" xmlns="http://www.w3.org/2000/svg"></svg>')
            } else {
                var svg_core = $('<svg id="result" width="100%" height="100%" version="1.1" viewBox="-' + x_size + ' -' + y_size + ' ' + x_size*2 + ' ' + y_size*2 + '" xmlns="http://www.w3.org/2000/svg"></svg>')
            }
        } else {
            var x_size = (this.global_x + 1) * node_distance_x * 0.5;
            var y_size = (sum_middle_layer + (this.global_y * node_distance_y) + sum_outer_layer) * 0.5;
            if (scale == 0){
                var svg_core = $('<svg id="result" width="' + x_size*2 + 'px" height="' + y_size*2 + 'px" version="1.1" viewBox="-' + 0.5 * node_distance_y + ' -' + y_size*2 + ' ' + x_size*2 + ' ' + y_size*2 + '" xmlns="http://www.w3.org/2000/svg"></svg>')
            } else {
                var svg_core = $('<svg id="result" width="100%" height="100%" version="1.1" viewBox="-' + 0.5 * node_distance_y + ' -' + y_size*2 + ' ' + x_size*2 + ' ' + y_size*2 + '" xmlns="http://www.w3.org/2000/svg"></svg>')
            }
        }
        
        if ($('#flexarrow').prop('checked') == true){
        
            var [arrow_size, arrow_start, arrow_stop] = middle_layers['arrow']
            var pointer_height = arrow_stop - arrow_start
            var pointer_length = pointer_height / 20
            var arrow_thickness = pointer_height / 4
            var steps = Math.round(this.global_x / (num_position + 1))
            
            if (steps < 1) {
                steps = 1
            }
            
            if (linear == 0){
                var [circle_c_x, circle_c_y] = this.circle_transform(this.global_x + 0.5 - pointer_length, arrow_start + arrow_thickness, theta, start_angle)
                var [circle_a_x, circle_a_y] = this.circle_transform(0-0.5, arrow_start + arrow_thickness, theta, start_angle)
                var [circle_b_x, circle_b_y] = this.circle_transform(0-0.5, arrow_stop - arrow_thickness, theta, start_angle)
                var [circle_d_x, circle_d_y] = this.circle_transform(this.global_x + 0.5 - pointer_length, arrow_stop - arrow_thickness, theta, start_angle)
                var [circle_f_x, circle_f_y] = this.circle_transform(this.global_x + 0.5 - pointer_length, arrow_stop, theta, start_angle)
                var [circle_g_x, circle_g_y] = this.circle_transform(this.global_x + 0.5, arrow_start + arrow_thickness * 2, theta, start_angle)
                var [circle_e_x, circle_e_y] = this.circle_transform(this.global_x + 0.5 - pointer_length, arrow_start, theta, start_angle)
                    
                if ((this.global_x) * theta > 180) {
                    var arc_flag = 1
                } else {
                    var arc_flag = 0
                }

                svg_arrow.push(
                    $('<path d="M ' + circle_c_x + ' ' + circle_c_y +
                    ' A ' + (arrow_start + arrow_thickness) + ' ' + (arrow_start + arrow_thickness) + ' 0 ' + arc_flag + ' 1 ' + circle_a_x + ' ' + circle_a_y +
                    ' L ' + circle_b_x + ' ' + circle_b_y +
                    ' A ' + (arrow_stop - arrow_thickness) + ' ' + (arrow_stop - arrow_thickness) + ' 0 ' + arc_flag + ' 0 ' + circle_d_x + ' ' + circle_d_y +
                    ' L ' + circle_f_x + ' ' + circle_f_y +
                    ' L ' + circle_g_x + ' ' + circle_g_y +
                    ' L ' + circle_e_x + ' ' + circle_e_y + 
                    ' Z" stroke-width="0" fill="slateGrey"></path>')
                )
                
                var [circle_h_x, circle_h_y] = this.circle_transform(0-0.5, arrow_start + arrow_thickness * 2, theta, start_angle)

                var rotate = start_angle
                if (rotate >= 90 && rotate <= 180) {
                    rotate += 180;
                    var anchor = "start"
                } else if (rotate >= 180 && rotate <= 270) {
                    rotate -= 180;
                    var anchor = "start"
                } else {
                    var anchor = "end"
                }

                var [angle_h_x, angle_h_y] = this.angle_transform(circle_h_x, circle_h_y, 180 - start_angle, text_offset)
                
                svg_text.push(
                    $('<text text-anchor="' + anchor + '" transform="rotate(-' + rotate + ' ' + angle_h_x + ' ' + angle_h_y +')" dominant-baseline="middle" x="' + angle_h_x + '" y="' + angle_h_y + '" dy="0" font-size="' + $('#label')[0].value + '" font-family="sans-serif" fill="black">Orientation</text>')
                )
                
            } else {
                var [circle_c_x, circle_c_y] = [(this.global_x + 0.5 - pointer_length) * node_distance_x, -(arrow_start + arrow_thickness)]
                var [circle_a_x, circle_a_y] = [(0-0.5) * node_distance_x, -(arrow_start + arrow_thickness)]
                var [circle_b_x, circle_b_y] = [(0-0.5) * node_distance_x, -(arrow_stop - arrow_thickness)]
                var [circle_d_x, circle_d_y] = [(this.global_x + 0.5 - pointer_length) * node_distance_x , -(arrow_stop - arrow_thickness)]
                var [circle_f_x, circle_f_y] = [(this.global_x + 0.5 - pointer_length) * node_distance_x, -arrow_stop]
                var [circle_g_x, circle_g_y] = [(this.global_x + 0.5) * node_distance_x, -(arrow_start + arrow_thickness * 2)]
                var [circle_e_x, circle_e_y] = [(this.global_x + 0.5 - pointer_length) * node_distance_x, -arrow_start]
                
                svg_arrow.push(
                    $('<path d="M ' + circle_c_x + ' ' + circle_c_y +
                    ' L ' + circle_a_x + ' ' + circle_a_y +
                    ' L ' + circle_b_x + ' ' + circle_b_y +
                    ' L ' + circle_d_x + ' ' + circle_d_y +
                    ' L ' + circle_f_x + ' ' + circle_f_y +
                    ' L ' + circle_g_x + ' ' + circle_g_y +
                    ' L ' + circle_e_x + ' ' + circle_e_y + 
                    ' Z" stroke-width="0" fill="slateGrey"></path>')
                )
                
                var [circle_h_x, circle_h_y] = [(0-0.5) * node_distance_x - text_offset, -(arrow_start + arrow_thickness * 2)]

                svg_text.push(
                    $('<text text-anchor="end" dominant-baseline="middle" x="' + circle_h_x + '" y="' + circle_h_y + '" dy="0" font-size="' + $('#label')[0].value + '" font-family="sans-serif" fill="black">Orientation</text>')
                )
            }
                
            var l = steps
            var k = 1
            while (k <= num_position) {
                
                if (linear == 0){
                    var [circle_l_x, circle_l_y] = this.circle_transform(l, arrow_start + arrow_thickness * 2, theta, start_angle)
                    var rotate = theta * (l+0.5) + start_angle
                    if (rotate >= 90 && rotate <= 180) {
                        rotate += 180;
                    } else if (rotate >= 180 && rotate <= 270) {
                        rotate -= 180;
                    }

                    svg_arrow.push(
                        $('<text text-anchor="middle" dominant-baseline="middle" transform="rotate(-' + rotate + ' ' + circle_l_x + ' ' + circle_l_y +')" x="' + circle_l_x + '" y="' + circle_l_y + '" dy="0" font-size="' + $('#label')[0].value + '" font-family="sans-serif" fill="white">' + l + '</text>')
                    )
                } else {
                    var [circle_l_x, circle_l_y] = [(l) * node_distance_x, -(arrow_start + arrow_thickness * 2)]
                    svg_arrow.push(
                        $('<text text-anchor="middle" dominant-baseline="middle" x="' + circle_l_x + '" y="' + circle_l_y + '" dy="0" font-size="' + $('#label')[0].value + '" font-family="sans-serif" fill="white">' + l + '</text>')
                    )
                }
                l += steps
                k += 1
            };
        }

        var edge_synteny = {}
        
        var end = new Date().getTime();
        var time = end - start;
        console.log('SVG drawing core', time, 'ms.')
        
        var start = new Date().getTime();
        for(var i in this.edges) {

            var edge = this.data['edges'][i];
            // console.log(edge)
            var edge_genomes = Object.keys(edge['directions'])
            
            var intersection = edge_genomes.filter(x => enabled.includes(x));
            if (intersection.length > 0) {
                var edge_genomes_length = edge_genomes.length;
                var color = this.pickcolor(edgecoloring, Object.keys(edge['directions']))
        
                if (saturation == 1){
                    var pick = this.lighter_color('#ffffff', color, edge_genomes_length / genome_size);
                } else {
                    var pick = color;
                }
            
                var source = edge['source']
                var target = edge['target']
                
                if (source != 'start' && target != 'stop' && edge['active'] == true){
            
                    var i_x = this.nodes[source]['position'][0]
                    var i_y = this.nodes[source]['position'][1]
                    var j_x = this.nodes[target]['position'][0]
                    var j_y = this.nodes[target]['position'][1]
                    
                    var dir_set = Object.values(edge['directions'])
    
                    if (dir_set.includes('L') && dir_set.includes('R')) {
                        var stroke = ' stroke-dasharray="' + line_thickness * 4 + ' ' + line_thickness + '" '
                    } else if (dir_set.includes('L')) {
                        var stroke = ' stroke-dasharray="' + line_thickness + '" '
                    } else {
                        var stroke = ''
                    }

                    if (dir_set.includes('R')) {
                        if (source in edge_synteny) {
                        } else {
                            edge_synteny[source] = {}
                        }
    
                        edge_synteny[source][target] = edge['route']
                    } else {
                        if (target in edge_synteny) {
                        } else {
                            edge_synteny[target] = {}
                        }

                        edge_synteny[target][source] = [...edge['route']].reverse()
                    }
    
                    var [graph_size, graph_start, graph_stop] = outer_layers['graph']
                    var i_y_size = sum_middle_layer + graph_start + graph_size * 0.5 + i_y * node_distance_y
                    var j_y_size = sum_middle_layer + graph_start + graph_size * 0.5 + j_y * node_distance_y
                    var draw = pick
                    var thickness = edge_thickness
                
                    if (linear == 0){
                        var [circle_i_x, circle_i_y] = this.circle_transform(i_x, i_y_size, theta, start_angle);
                        var [circle_j_x, circle_j_y] = this.circle_transform(j_x, j_y_size, theta, start_angle);
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
                                var o_y_size = sum_middle_layer + graph_start + graph_size * 0.5 + o_y * node_distance_y
                                var n_y_size = sum_middle_layer + graph_start + graph_size * 0.5 + n_y * node_distance_y
        
                                if (linear == 0){
                                    var [circle_n_x, circle_n_y] = this.circle_transform(n_x, n_y_size, theta, start_angle);
                                } else {
                                    var [circle_n_x, circle_n_y] = [(n_x) * node_distance_x, -n_y_size];
                                }
        
                                if (o_y == n_y) {
                                    if (linear == 0){
                                        route_edge += ' A ' + o_y_size  + ' ' + n_y_size + ' 0 0 0 ' + circle_n_x + ' ' + circle_n_y
                                    } else {
                                        route_edge += ' L ' + circle_n_x + ' ' + circle_n_y
                                    }
                                } else {
                                    route_edge += ' L ' + circle_n_x + ' ' + circle_n_y
                                }
        
                                var o_y = n_y
                            }
        
                            if (o_y == j_y) {
                                if (linear == 0){
                                    route_edge += ' A ' + o_y_size  + ' ' + j_y_size + ' 0 0 0 ' + circle_j_x + ' ' + circle_j_y + '"' + stroke + ' stroke="' + draw + '" stroke-width="' + thickness + '" fill="none"/>'
                                } else {
                                    route_edge += ' L ' + circle_j_x + ' ' + circle_j_y + '"' + stroke + ' stroke="' + draw + '" stroke-width="' + thickness + '" fill="none"/>'
                                }
                            } else {
                                route_edge += ' L ' + circle_j_x + ' ' + circle_j_y + '"' + stroke + ' stroke="' + draw + '" stroke-width="' + thickness + '" fill="none"/>'
                            }
                        }

                        svg_edges.push($(route_edge))
                    }
                }
            }
        };

        for (var genome of this.genomes) {
            var layer_name = genome + 'layer'
            if (Object.keys(middle_layers).includes(layer_name)){
        
                var [layer_width, layer_start, layer_stop] = middle_layers[layer_name]
            
                if (linear == 0){
                    var [circle_a_x, circle_a_y] = this.circle_transform(0-0.5, layer_start, theta, start_angle)
                    var [circle_b_x, circle_b_y] = this.circle_transform(0-0.5, layer_stop, theta, start_angle)
                    var [circle_c_x, circle_c_y] = this.circle_transform(this.global_x + 0.5, layer_start, theta, start_angle)
                    var [circle_d_x, circle_d_y] = this.circle_transform(this.global_x + 0.5, layer_stop, theta, start_angle)
                    
                    if ((this.global_x) * theta > 180) {
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
                    var [circle_c_x, circle_c_y] = [(this.global_x + 0.5) * node_distance_x, -layer_start]
                    var [circle_d_x, circle_d_y] = [(this.global_x + 0.5) * node_distance_x, -layer_stop]
                    
                    svg_genome_tracks[genome].push(
                        $('<path d="M ' + circle_c_x + ' ' + circle_c_y +
                        ' L ' + circle_a_x + ' ' + circle_a_y +
                        ' L ' + circle_b_x + ' ' + circle_b_y +
                        ' L ' + circle_d_x + ' ' + circle_d_y +
                        ' Z" stroke-width="0" fill="' + layer_color + '"></path>')
                    )
                }

                var sorted_keys = Object.keys(this.synteny[genome]).sort(function (a, b) {return parseInt(a) - parseInt(b);});
                if (layer_width >= track_line_width) {

                    layer_width -= track_line_width
                    layer_start += track_line_width * 0.5
                    layer_stop -= track_line_width * 0.5

                    var draw = edgecoloring[genome][1]
                    var thickness = track_line_width
                    var stroke = ''

                    var edge_chain = []
                    for(var i of sorted_keys){
    
                        var j = (parseInt(i) + 1).toString()
                        var node_i = this.synteny[genome][i]
                        var node_j = this.synteny[genome][j]

                        var pass_on = 1
                        
                        if (typeof(node_i) != 'undefined' && typeof(node_j) != 'undefined'){

                            var i_x = this.nodes[node_i]['position'][0]
                            var i_y = this.nodes[node_i]['position'][1]
                            var j_x = this.nodes[node_j]['position'][0]
                            var j_y = this.nodes[node_j]['position'][1]

                            if (node_i in edge_synteny) {
                                if (node_j in edge_synteny[node_i]) {
                                    if (edge_chain.length == 0) {
                                        edge_chain.push([i_x, i_y])
                                    }

                                    edge_chain.push(...edge_synteny[node_i][node_j])
                                    edge_chain.push([j_x, j_y])
                                } else {pass_on = 0}
                            } else {pass_on = 0}
                        } else {pass_on = 0}

                        if (pass_on == 0) {
                            var route_edge = ''

                            for (var p=0; p < edge_chain.length - 1; p++){
                                var p_x = edge_chain[p][0]
                                var p_y = edge_chain[p][1]
                                var p_y_size = layer_start + p_y * (layer_width / this.global_y)

                                if (linear == 0){
                                    var [circle_p_x, circle_p_y] = this.circle_transform(p_x, p_y_size, theta, start_angle);
                                } else {
                                    var [circle_p_x, circle_p_y] = [(p_x) * node_distance_x, -p_y_size];
                                }

                                if (route_edge == "") {
                                    route_edge = '<path class="path" d="M ' + circle_p_x + ' ' + circle_p_y
                                }
                                
                                var q_x = edge_chain[p+1][0]
                                var q_y = edge_chain[p+1][1]
                                var q_y_size = layer_start + q_y * (layer_width / this.global_y)

                                if (linear == 0){
                                    var [circle_q_x, circle_q_y] = this.circle_transform(q_x, q_y_size, theta, start_angle);
                                } else {
                                    var [circle_q_x, circle_q_y] = [(q_x) * node_distance_x, -q_y_size];
                                }

                                if (p_y == q_y) {
                                    if (linear == 0){
                                        route_edge += ' A ' + p_y_size  + ' ' + q_y_size + ' 0 0 0 ' + circle_q_x + ' ' + circle_q_y
                                    } else {
                                        route_edge += ' L ' + circle_q_x + ' ' + circle_q_y
                                    }
                                } else {
                                    route_edge += ' L ' + circle_q_x + ' ' + circle_q_y
                                }
                            }

                            // console.log(route_edge)
                            edge_chain = []

                            if (route_edge != "") {
                                svg_genome_tracks[genome].push($(route_edge + '" stroke="' + draw + '" stroke-width="' + thickness + '" fill="none"/>'))
                            }
                        }                            
                    }
                }
            }
        }
        
        var end = new Date().getTime();
        var time = end - start;
        console.log('SVG drawing lines', time, 'ms.')
        
        var globalbackbone_pos = [];
        
        var start = new Date().getTime();
        var global_values = []
        for(var k in this.nodes) {
        
            var node = this.data['nodes'][k];
            var node_genomes = Object.keys(node['gene_calls']);
            var intersection = node_genomes.filter(x => enabled.includes(x));
            if (intersection.length > 0) {
            
                var node_genomes_length = node_genomes.length
            
                var k_x = node['position'][0]
                var k_y = node['position'][1]
                var node_group = node['group']
                var node_type = node['type']
                
                if (node['layer']['backbone'] == 1) {
                    globalbackbone_pos.push(k_x)
                }
            
                if (!node_group) {
                    var node_class = 'class="node'
                    var x_value_start = k_x - 0.5
                    var x_value_stop = k_x + 0.5
        
                } else {
                    var node_class = 'stroke-opacity="0" fill-opacity="0" class="pseudo'
                    var group = this.group_dict[node_group]
                    var group_size = group.length
                    var group_compress = $('#groupcompress')[0].value
                    var group_size_compressed = Math.round(group_size * group_compress)
        
                    if (group_size_compressed == 0) {
                        group_size_compressed = 1
                    }
    
                    var z_x = this.data['nodes'][group[0]]['position'][0]
    
                    var fraction = group_size_compressed / (group_size)
                    var group_id = group.findIndex(x => x === k)
                    
                    var x_value_start = z_x - (1 - group_id * fraction) + 0.5
                    var x_value_stop = z_x - (1 - (group_id + 1) * fraction) + 0.5
    
                }
    
                var color = this.pickcolor(edgecoloring, Object.keys(node['gene_calls']))
    
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
                    console.warn('Unexpected node type:', node_type)
                }
    
                if (saturation == 1) {
                    var draw = this.lighter_color('#ffffff', color, node_genomes_length / genome_size);
                    var draw2 = this.lighter_color('#ffffff', node_color, node_genomes_length / genome_size)
                } else {
                    var draw = color;
                    var draw2 = node_color
                }
    
                var [graph_size, graph_start, graph_stop] = outer_layers['graph']
                var k_y_size = sum_middle_layer + graph_start + graph_size * 0.5 + k_y * node_distance_y
    
                if (linear == 0) {
                    var [circle_k_x, circle_k_y] = this.circle_transform(k_x, k_y_size, theta, start_angle);
                } else {
                    var [circle_k_x, circle_k_y] = [(k_x) * node_distance_x, -k_y_size];
                }
    
                svg_nodes.push(
                    $('<circle ' + node_class + '" id="' + k + '" cx="' + circle_k_x + '" cy="' + circle_k_y + '" r="' + node_size + '" fill="' + draw2 + '" stroke="' + draw + '" stroke-width="' + node_thickness + '"/>')
                )
    
                var [search_size, search_start, search_stop] = middle_layers['search']
                
                var i_x = k_x - 0.5
                var i_y = search_start
                var j_x = k_x + 0.5
                var j_y = search_stop
                
                if (!global_values.includes(k_x)) {
                    svg_search.push(this.create_rectangle(i_x, i_y, j_x, j_y, theta, start_angle, node_distance_x, linear, 'white', k_x))
                }
    
                for (var layer_name of this.layers) {
                
                    if ($('#flex' + layer_name).prop('checked') == true){
                
                        var value = node['layer'][layer_name] 
                        var max = parseFloat($('#' + layer_name + '_max')[0].value);
                        var min = parseFloat($('#' + layer_name + '_min')[0].value);

                        if (value < min) {
                            value = min
                        }
                        
                        if (value > max) {
                            value = max
                        }
                        
                        var [layer_width, layer_start, layer_stop] = outer_layers[layer_name]
                        var k_y_size = sum_middle_layer + k_y * node_distance_y
                        
                        var i_x = x_value_start
                        var i_y = layer_start + k_y_size
                        var j_x = x_value_stop
                        var j_y = layer_stop + k_y_size
                        var color = this.lighter_color('#00ff00', '#ff0000', (value-min) / (max-min))
                        
                        svg_heatmaps.push(this.create_rectangle(i_x, i_y, j_x, j_y, theta, start_angle, node_distance_x, linear, color))
                    }
                }
            }
        }
        
        global_values.push(k_x)
        
        if ($('#flexglobalbackbone').prop('checked') == true){
            var k_x = 0
            while (k_x <= this.global_x) {
        
                var [globalbackbone_size, globalbackbone_start, globalbackbone_stop] = middle_layers['back_vs_non_back']
        
                var i_x = k_x - 0.5
                var i_y = globalbackbone_start
                var j_x = k_x + 0.5
                var j_y = globalbackbone_stop
                
                if (globalbackbone_pos.includes(k_x)) {
                    var color = back_color
                } else {
                    var color = non_back_color
                }
        
                svg_backbone.push(this.create_rectangle(i_x, i_y, j_x, j_y, theta, start_angle, node_distance_x, linear, color, k_x))
                k_x = k_x + 1
            }
        }
        
        for(var [l, group] of Object.entries(this.group_dict)) {
        
            var group_length = group.length
            var group_x = group.map((a) => (this.data['nodes'][a]['position'][0]));
            
            var ind_max = group_x.indexOf(Math.max.apply(Math, group_x))
            var ind_min = group_x.indexOf(Math.min.apply(Math, group_x))
            
            var left_node_name = group[ind_min]
            var right_node_name = group[ind_max]
        
            var left_node = this.data['nodes'][left_node_name];
            var right_node = this.data['nodes'][right_node_name];
        
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
                var color = this.pickcolor(edgecoloring, Object.keys(left_node['gene_calls']))
        
                if (saturation == 1) {
                    var draw = this.lighter_color('#ffffff', color, group_genomes_length / genome_size);
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
                    var color = this.lighter_color('#ffffff', group_color, group_genomes_length / genome_size)
                } else {
                    var color = group_color
                }
        
                if ((l_x - m_x) * theta >= 180) {
                    var arc_flag = 1
                } else {
                    var arc_flag = 0
                }
        
                if (linear == 0) {
                    var [a_x, a_y] = this.circle_transform(i_x, i_y, theta, start_angle)
                    var [b_x, b_y] = this.circle_transform(i_x, j_y, theta, start_angle)
                    var [c_x, c_y] = this.circle_transform(j_x, i_y, theta, start_angle)
                    var [d_x, d_y] = this.circle_transform(j_x, j_y, theta, start_angle)
            
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
        
        for (var layer_name of this.layers) {
        
            if ($('#flex' + layer_name).prop('checked') == true){
                var [layer_width, layer_start, layer_stop] = outer_layers[layer_name]
                var y_size = sum_middle_layer + layer_width * 0.5
            
                if (linear == 0){
                    var [circle_x, circle_y] = this.circle_transform(0-0.5, (layer_start + y_size), theta, start_angle)

                    var rotate = start_angle
                    if (rotate >= 90 && rotate <= 180) {
                        rotate += 180;
                        var anchor = "start"
                    } else if (rotate >= 180 && rotate <= 270) {
                        rotate -= 180;
                        var anchor = "start"
                    } else {
                        var anchor = "end"
                    }

                    svg_text.push(
                        $('<text text-anchor="' + anchor + '" transform="rotate(-' + rotate + ' ' + circle_x + ' ' + circle_y +')" dominant-baseline="middle" x="' + circle_x + '" y="' + circle_y + '" dy="0" font-size="' + $('#label')[0].value + '" font-family="sans-serif" fill="black">' + layer_name + '</text>')
                    )
    
                } else {
                    var [circle_x, circle_y] = [(0-0.5) * node_distance_x, -(layer_start + y_size)]

                    svg_text.push(
                        $('<text text-anchor="end" dominant-baseline="middle" x="' + circle_x + '" y="' + circle_y + '" dy="0" font-size="' + $('#label')[0].value + '" font-family="sans-serif" fill="black">' + layer_name + '</text>')
                    )
                }
            }
        }
        
        var item_dist = {}
        for (var genome_name of item_order) {
            if ($('#flex' + genome_name + 'layer').prop('checked') == true){
                var [layer_width, layer_start, layer_stop] = middle_layers[genome_name + 'layer']
        
                if (layer_width >= edge_thickness) {
                    var y_size = layer_start + layer_width * 0.5
                    if (linear == 0){
                        var [circle_x, circle_y] = this.circle_transform(0-0.5, y_size, theta, start_angle)

                        var rotate = start_angle
                        if (rotate >= 90 && rotate <= 180) {
                            rotate += 180;
                            var anchor = "start"
                        } else if (rotate >= 180 && rotate <= 270) {
                            rotate -= 180;
                            var anchor = "start"
                        } else {
                            var anchor = "end"
                        }

                        var [angle_x, angle_y] = this.angle_transform(circle_x, circle_y, 180 - start_angle, text_offset)
                        
                        svg_text.push(
                            $('<text text-anchor="' + anchor + '" transform="rotate(-' + rotate + ' ' + angle_x + ' ' + angle_y +')" dominant-baseline="middle" x="' + angle_x + '" y="' + angle_y + '" dy="0" font-size="' + $('#label')[0].value + '" font-family="sans-serif" fill="black">' + genome_name + '</text>')
                        )
                    } else {
                        var [circle_x, circle_y] = [(0-0.5) * node_distance_x - text_offset, -y_size]

                        svg_text.push(
                            $('<text text-anchor="end" dominant-baseline="middle" x="' + circle_x + '" y="' + circle_y + '" dy="0" font-size="' + $('#label')[0].value + '" font-family="sans-serif" fill="black">' + genome_name + '</text>')
                        )
                    }
        
                    item_dist[genome_name] = y_size
                }
            }
        }

        // Only draw newick tree if we have one (order will be empty if no newick data)
        if ($('#flextree').prop('checked') == true && order.length > 0){
            svg_tree = this.draw_newick(order, item_dist, max_dist, offset, tree_length, tree_thickness, theta, start_angle, end_angle, linear, node_distance_x)
        }
        
        if ($('#flexregionlabels').prop('checked') && this.data['regions']) {
            // Base radial distance without any node layers.
            var base_outer_r = current_outer_stop + sum_middle_layer;

            for (var [rid, rinfo] of Object.entries(this.data['regions'])) {
                var x_min = rinfo['x_min'];
                var x_max = rinfo['x_max'];
                var x_span = x_max - x_min;

                // Skip single-position (backbone-only) regions
                if (x_span === 0) continue;

                // Find the tallest node in this region so the label clears it.
                var region_max_y = 0;
                for (var [nid, ndata] of Object.entries(this.nodes)) {
                    var nx = ndata['position'][0];
                    if (nx >= x_min && nx <= x_max) {
                        var ny = ndata['position'][1];
                        if (ny > region_max_y) region_max_y = ny;
                    }
                }
                var outer_content_r = base_outer_r + region_max_y * node_distance_y;

                var x_mid = (x_min + x_max) / 2;
                var svg_region_width = x_span * node_distance_x;

                // Store geometry as data attributes; position and font-size are set
                // dynamically by the zoom handler so the label always clears the graph.
                if (linear == 0) {
                    var label_angle = theta * (x_mid + 0.5) + start_angle;
                    var angle_start = theta * x_min + start_angle;
                    var angle_end   = theta * (x_max + 1) + start_angle;
                    svg_region_labels.push($(
                        `<text class="region-label"
                            data-region-id="${rid}"
                            data-svg-width="${svg_region_width}"
                            data-layout="circular"
                            data-outer-r="${outer_content_r}"
                            data-inner-r="${start_offset}"
                            data-angle="${label_angle}"
                            data-angle-start="${angle_start}"
                            data-angle-end="${angle_end}"
                            text-anchor="middle" dominant-baseline="middle"
                            x="0" y="0" font-size="1" font-family="sans-serif"
                            fill="#555555" opacity="0.85"
                            style="display:none">#${rid}</text>`
                    ));
                } else {
                    var label_x   = x_mid * node_distance_x;
                    var x_start_svg = x_min * node_distance_x;
                    var x_end_svg   = (x_max + 1) * node_distance_x;
                    svg_region_labels.push($(
                        `<text class="region-label"
                            data-region-id="${rid}"
                            data-svg-width="${svg_region_width}"
                            data-layout="linear"
                            data-outer-r="${outer_content_r}"
                            data-label-x="${label_x}"
                            data-x-start="${x_start_svg}"
                            data-x-end="${x_end_svg}"
                            text-anchor="middle" dominant-baseline="middle"
                            x="${label_x}" y="0" font-size="1" font-family="sans-serif"
                            fill="#555555" opacity="0.85"
                            style="display:none">#${rid}</text>`
                    ));
                }
            }
        }

        var end = new Date().getTime();
        var time = end - start;
        console.log('SVG drawing remaining elements', time, 'ms.')

        var svg_group = $('<g></g>')
        for (var item of svg_arrow) svg_group.append(item);
        svg_core.append(svg_group);
        
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
        for (var item of svg_tree) svg_group.append(item);
        svg_core.append(svg_group);

        svg_core.append($('<g id="bin-rings-group"></g>'));

        var svg_group = $('<g id="region-labels-group"></g>')
        for (var item of svg_region_labels) svg_group.append(item);
        svg_core.append(svg_group);

        var svg_group = $('<g></g>')
        for (var item of svg_text) svg_group.append(item);
        svg_core.append(svg_group);

        // Store layout parameters for dynamic bin ring drawing
        this._layout = {
            theta: theta,
            start_angle: start_angle,
            end_angle: end_angle,
            node_distance_x: node_distance_x,
            node_distance_y: node_distance_y,
            sum_middle_layer: sum_middle_layer,
            outer_layers: outer_layers,
            linear: linear,
            current_outer_stop: current_outer_stop,
            start_offset: start_offset,
        };

        return svg_core
    }

    deg2rad(degrees) {
        return degrees * Math.PI/180;
    }

    download_blob(blob, name) {

        var blob_url = URL.createObjectURL(blob);
        var link = document.createElement("a");
        
        link.href = blob_url;
        link.download = name;
        
        document.body.appendChild(link);
        
        link.dispatchEvent(new MouseEvent('click', {
            bubbles: true,
            cancelable: true,
            view: window
        }));
        
        document.body.removeChild(link);
    }

    // NOTE - From https://stackoverflow.com/questions/19799777/how-to-add-transparency-information-to-a-hex-color-code
    addAlpha(color, opacity) {
        var _opacity = Math.round(Math.min(Math.max(opacity ?? 1, 0), 1) * 255);
        return color + _opacity.toString(16).toUpperCase();
    }

    // NOTE - From https://coderwall.com/p/z8uxzw/javascript-color-blender
    int_to_hex(num) {
        var hex = Math.round(num).toString(16);
        if (hex.length == 1)
            hex = '0' + hex;
        return hex;
    }

    lighter_color(color1, color2, percentage, threshold=0.25) {
        percentage = threshold + (1 - threshold) * percentage
        
        color1 = color1.substring(1);
        color2 = color2.substring(1);
        
        color1 = [parseInt(color1[0] + color1[1], 16), parseInt(color1[2] + color1[3], 16), parseInt(color1[4] + color1[5], 16)];
        color2 = [parseInt(color2[0] + color2[1], 16), parseInt(color2[2] + color2[3], 16), parseInt(color2[4] + color2[5], 16)];
        
        var color3 = [
        (1 - percentage) * color1[0] + percentage * color2[0],
        (1 - percentage) * color1[1] + percentage * color2[1],
        (1 - percentage) * color1[2] + percentage * color2[2]
        ];
        
        var color4 = '#' + this.int_to_hex(color3[0]) + this.int_to_hex(color3[1]) + this.int_to_hex(color3[2]);

        return color4
    }

    create_rectangle(i_x, i_y, j_x, j_y, theta, start_angle, node_distance_x, linear, color, id='') {
    
        if (id != '') {
            var extra = '" id="' + id
        } else {
            var extra = ''
        }
        
        if (linear == 0) {
            var [a_x, a_y] = this.circle_transform(i_x, i_y, theta, start_angle)
            var [b_x, b_y] = this.circle_transform(j_x, i_y, theta, start_angle)
            var [c_x, c_y] = this.circle_transform(i_x, j_y, theta, start_angle)
            var [d_x, d_y] = this.circle_transform(j_x, j_y, theta, start_angle)
            
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

    circle_transform(x, y, theta, start_angle) {
        var circle_x = y * Math.sin(this.deg2rad(theta * (x+0.5) + start_angle))
        var circle_y = y * Math.cos(this.deg2rad(theta * (x+0.5) + start_angle))
        return [circle_x, circle_y]
    }

    pickcolor(edgecoloring, genomes) {
    
        var array = []
        for (var name of genomes) {
            array.push(edgecoloring[name])
        }

        var sortedArray = array.sort(function(a, b) {
            return a[0] - b[0];
        });
        
        return sortedArray[0][1]
    }

    angle_transform(x, y, angle, dist) {
        var angle_x = x + dist * Math.cos(this.deg2rad(angle))
        var angle_y = y + dist * Math.sin(this.deg2rad(angle))
        return [angle_x, angle_y]
    }

    draw_newick(order, item_dist, max_dist, offset, max_size, line_thickness, theta, start_angle, end_angle, linear, node_distance_x) {

        var output = []
        var saving_positions = {}
        var saving_ids = {}
        
        var i = 0
        
        for (var [item, start, end] of order){
        
            var start_fraction = ((max_size * (start / max_dist) - max_size) - offset) * (-1)
            var end_fraction = ((max_size * (end / max_dist) - max_size) - offset) * (-1)
            
            if (item != 'branching') {
                var y_value = item_dist[item]

                if (linear == 0){
                    var [x, y] = this.circle_transform(this.global_x + 0.5, y_value, theta, start_angle)
                    var [a_x, a_y] = this.angle_transform(x, y, 360 - end_angle, end_fraction)
                    var [b_x, b_y] = this.angle_transform(x, y, 360 - end_angle, start_fraction)
                } else {
                    var [a_x, a_y] = [(this.global_x + 0.5) * node_distance_x + end_fraction, -y_value]
                    var [b_x, b_y] = [(this.global_x + 0.5) * node_distance_x + start_fraction, -y_value]
                }

                output.push('<path d="M ' + a_x + ' ' + a_y + ' L ' + b_x + ' ' + b_y + '" stroke-width="' + line_thickness + '" stroke="black"></path>')

                if (Object.values(saving_ids).includes(start_fraction)){
                    var saving_id = Object.keys(saving_ids)[Object.values(saving_ids).indexOf(start_fraction)]
                    saving_positions[saving_id].push(y_value)
                } else {
                    saving_ids[i] = start_fraction
                    saving_positions[i] = [y_value]
                    i = i + 1
                }
                if (end_fraction != max_size){
                    if (linear == 0){
                        var [c_x, c_y] = this.angle_transform(x, y, 360 - end_angle, offset)
                    } else {
                        var [c_x, c_y] = [(this.global_x + 0.5) * node_distance_x + offset, -y_value]
                    }
                        
                    output.push('<path d="M ' + a_x + ' ' + a_y + ' L ' + c_x + ' ' + c_y + '" stroke-dasharray="' + line_thickness + '" stroke-width="' + line_thickness + '" stroke="lightgray"></path>')
                }
            } else {
        
                if (Object.values(saving_ids).includes(end_fraction)){
                    var saving_id_main = Object.keys(saving_ids)[Object.values(saving_ids).indexOf(end_fraction)]
                    var sorted_positions = saving_positions[saving_id_main].sort()
            
                    for (var j = 0; j < sorted_positions.length -1; j++) {
            
                        var y_value_i = sorted_positions[j]
                        var y_value_j = sorted_positions[j+1]

                        if (linear == 0){
                            var [i_x, i_y] = this.circle_transform(this.global_x + 0.5, y_value_i, theta, start_angle)
                            var [m_x, m_y] = this.angle_transform(i_x, i_y, 360 - end_angle, end_fraction)
    
                            var [j_x, j_y] = this.circle_transform(this.global_x + 0.5, y_value_j, theta, start_angle)
                            var [n_x, n_y] = this.angle_transform(j_x, j_y, 360 - end_angle, end_fraction)
                        } else {
                            var [m_x, m_y] = [(this.global_x + 0.5) * node_distance_x + end_fraction, -y_value_i]
                            var [n_x, n_y] = [(this.global_x + 0.5) * node_distance_x + end_fraction, -y_value_j]
                        }
                            
                        output.push('<path d="M ' + m_x + ' ' + m_y + ' L ' + n_x + ' ' + n_y + '" stroke-width="' + line_thickness + '" stroke="black"></path>')
                    }
                } else {
                    var saving_id_main = -1
                }
                
                y_value = Math.min(...sorted_positions) + (Math.max(...sorted_positions) - Math.min(...sorted_positions)) / 2

                if (linear == 0){
                    var [x, y] = this.circle_transform(this.global_x + 0.5, y_value, theta, start_angle)
                    var [a_x, a_y] = this.angle_transform(x, y, 360 - end_angle, end_fraction)
                    var [b_x, b_y] = this.angle_transform(x, y, 360 - end_angle, start_fraction)
                } else {
                    var [a_x, a_y] = [(this.global_x + 0.5) * node_distance_x + end_fraction, -y_value]
                    var [b_x, b_y] = [(this.global_x + 0.5) * node_distance_x + start_fraction, -y_value]
                }
                
                output.push('<path d="M ' + a_x + ' ' + a_y + ' L ' + b_x + ' ' + b_y + '" stroke-width="' + line_thickness + '" stroke="black"></path>')

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

    newick_to_order(string, prior = 0) {
        var result = []
        var newick = string.replace(' ', '')
        
        if (newick[0] == '('){
            var bracket_open = 0
            var bracket_closed = 0
            
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
                        result.push(...this.newick_to_order(sub_newick, prior + parseFloat(value)))
        
                        var next_iter = parts.slice(1).join(',')
                        result.push(...this.newick_to_order(next_iter, prior))
                    } else {
                        if (rest.includes(':')) {
                            var value = rest.split(':').slice(1)
                        } else {
                            var value = 0 
                        }
                        result.push(...[['branching', prior, prior + parseFloat(value)]])
                        result.push(...this.newick_to_order(sub_newick, prior + parseFloat(value)))
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
                result.push(...this.newick_to_order(next_iter, prior))    
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
    
    start_draw() {
        var new_settings_dict = {};
        
        new_settings_dict['condtr'] = parseInt($('#condtr')[0].value);
        new_settings_dict['maxlength'] = parseInt($('#maxlength')[0].value);
        new_settings_dict['groupcompress'] = parseFloat($('#groupcompress')[0].value);
        
        if (JSON.stringify(this.settings_dict) !== JSON.stringify(new_settings_dict)) {
            this.rerun_JSON(new_settings_dict);
        }
        
        $('#svgbox').css('opacity', '0.2');

        $.ajax({
            url: "/pangraph/get_pangraph_json_data",
            type: "GET",
            cache: false,
            contentType: "application/json",
            dataType: "json",
            success: (data) => {
                this.data = data['data'];
                this.initialize_variables();
                this.settings_dict = JSON.parse(JSON.stringify(new_settings_dict));
                this.main_draw();
                $('#svgbox').css('opacity', '');
            },
            error: (err) => {
                $('#svgbox').css('opacity', '');
                toastr.error('Failed to reload pan-graph data.', 'Error');
            }
        })
    }

    hide_tippy(instance) {
        
        if (instance.reference.id.startsWith('GCG_')){
            var element_id = this.group_dict[instance.reference.id][0]
        } else {
            var element_id = instance.reference.id
        }

        $('#number_sgc')[0].innerText = '';
        $('#number_gc')[0].innerText = '';
        $('#number_type')[0].innerText = '';
        $('#number_position')[0].innerText = '0';
        
        for (var source of this.functional_annotation_sources_available) {
            $('#number_' + source)[0].innerText = '';
        }
        
        for (var genome of this.genomes) {
            $('#number_' + genome)[0].innerText = '0';
        }

        for (var layer of this.layers) {
            $('#number_' + layer)[0].innerText = '0';
        }
    }
    
    async show_tippy(instance) {

        if (instance.reference.id.startsWith('GCG_')){
            var element_id = this.group_dict[instance.reference.id][0]
        } else {
            var element_id = instance.reference.id
        }
        
        var d = await this.get_gene_cluster_consensus_functions([element_id]);
        for (var source of this.functional_annotation_sources_available) {

            var gene_cluster_functions = d['data'] && d['data'][element_id];
            if (gene_cluster_functions && source in gene_cluster_functions) {
                var func = gene_cluster_functions[source]['function']
                func === undefined | func === null ? func = '' : func = func
                $('#number_' + source)[0].innerText = func;
            } else {
                $('#number_' + source)[0].innerText = '';
            }
        }
        
        for (var genome of this.genomes) {
            if (genome in this.data['nodes'][element_id]['gene_calls']) {
                var gene_call = this.data['nodes'][element_id]['gene_calls'][genome]
                $('#number_' + genome)[0].innerText = '1';
            } else {
                $('#number_' + genome)[0].innerText = '0';
            }
        }

        for (var layer of this.layers) {
            if (layer in this.data['nodes'][element_id]['layer']) {
                var value = this.data['nodes'][element_id]['layer'][layer]
                $('#number_' + layer)[0].innerText = value;
            } else {
                $('#number_' + layer)[0].innerText = '';
            }
        }

        $('#number_sgc')[0].innerText = element_id;
        $('#number_gc')[0].innerText = this.data['nodes'][element_id]['gene_cluster'];
        $('#number_type')[0].innerText = this.data['nodes'][element_id]['type'];
        $('#number_position')[0].innerText = parseInt(this.data['nodes'][element_id]['position']);

        var raw_type = this.data['nodes'][element_id]['type'];
        var gene_cluster = this.data['nodes'][element_id]['gene_cluster'];
        var num_genomes = Object.keys(this.data['nodes'][element_id]['gene_calls']).length;
        var total_genomes = this.genomes.length;

        const type_labels = {
            'core':          'Core',
            'accessory':     'Accessory',
            'singleton':     'Singleton',
            'duplication':   'Multi-Copy',
            'rearrangement': 'Rearrangement',
            'rna':           'tRNA',
        };
        var type_label = type_labels[raw_type] || raw_type;
        var type_display = (raw_type === 'core')
            ? type_label
            : type_label + ' (' + num_genomes + ' of ' + total_genomes + ' genomes)';

        if (raw_type === 'duplication') {
            var copies_in_graph = Object.values(this.data['nodes']).filter(n => n['gene_cluster'] === gene_cluster).length;
            type_display += ', ' + copies_in_graph + ' cop' + (copies_in_graph === 1 ? 'y' : 'ies') + ' in graph';
        }

        instance.setContent('<strong>' + gene_cluster + '</strong><br />' + type_display);
        $('#number_type')[0].innerText = type_display;
    }

    press_down(instance) {
        if (instance.button !== 0) return;

        this.cur_xpos = instance.offsetX;
        this.cur_ypos = instance.offsetY;

        this.isDown = true;
        this.diff = 0;

        var svgEl = document.getElementById('result');
        if (instance.altKey && typeof(svgEl) != 'undefined' && svgEl != null) {
            var svgEl = document.getElementById('result');
            var viewport = svgEl.querySelector('.svg-pan-zoom_viewport');
            var pt = svgEl.createSVGPoint();
            pt.x = instance.clientX;
            pt.y = instance.clientY;

            var svgP = pt.matrixTransform(viewport.getScreenCTM().inverse());
            this.svg_cur_xpos = svgP.x
            this.svg_cur_ypos = svgP.y

            this.selection = document.createElementNS("http://www.w3.org/2000/svg", "rect");
            this.selection.setAttribute("x", 0);
            this.selection.setAttribute("y", 0);
            this.selection.setAttribute("width", 0);
            this.selection.setAttribute("height", 0);

            var bin_color = $('#' + this.current_bin_id + '_color').attr('color')
            var fill_color = this.addAlpha(bin_color, 0.1);

            this.selection.setAttribute("fill", fill_color);
            viewport.appendChild(this.selection);
        } else if (instance.shiftKey) {
            this.drawing_zoom = true;
            this.zoomBox = {start_x: instance.clientX, start_y: instance.clientY};
            $('#divzoom').css({"top": 0, "left": 0, "width": 0, "height": 0});
            $('#divzoom').show();
        }
    }

    press_move(instance) {

        if (this.isDown === true) {
            this.new_xpos = instance.offsetX
            this.new_ypos = instance.offsetY
            this.diff += Math.sqrt((this.new_xpos-this.cur_xpos)^2+(this.new_ypos-this.cur_ypos)^2)

            var svgEl = document.getElementById('result');
            if (instance.shiftKey && this.drawing_zoom) {
                var _top = this.zoomBox.start_y > instance.clientY ? instance.clientY : this.zoomBox.start_y;
                var _left = this.zoomBox.start_x > instance.clientX ? instance.clientX : this.zoomBox.start_x;
                var _height = Math.abs(this.zoomBox.start_y - instance.clientY);
                var _width = Math.abs(this.zoomBox.start_x - instance.clientX);

                var divzoom = document.getElementById('divzoom');
                divzoom.style.top = _top + "px";
                divzoom.style.left = _left + "px";
                divzoom.style.width = _width + "px";
                divzoom.style.height = _height + "px";

                clearTextSelection();
            } else if (instance.altKey && typeof(svgEl) != 'undefined' && svgEl != null) {
                var viewport = svgEl.querySelector('.svg-pan-zoom_viewport');
                var pt = svgEl.createSVGPoint();
                pt.x = instance.clientX;
                pt.y = instance.clientY;

                var svgP = pt.matrixTransform(viewport.getScreenCTM().inverse());
                this.svg_new_xpos = svgP.x
                this.svg_new_ypos = svgP.y

                var min_svg_xpos = Math.min(this.svg_cur_xpos, this.svg_new_xpos);
                var min_svg_ypos = Math.min(this.svg_cur_ypos, this.svg_new_ypos);
                var width = Math.abs(this.svg_cur_xpos - this.svg_new_xpos);
                var height = Math.abs(this.svg_cur_ypos - this.svg_new_ypos);

                this.selection.setAttribute("x", min_svg_xpos);
                this.selection.setAttribute("y", min_svg_ypos);
                this.selection.setAttribute("width", width);
                this.selection.setAttribute("height", height);
            } else if (!instance.altKey && !instance.shiftKey && typeof(this.panZoomInstance) != 'undefined' && this.panZoomInstance != null) {
                this.panZoomInstance.panBy({x: this.new_xpos-this.cur_xpos, y: this.new_ypos-this.cur_ypos})
            }

            this.cur_xpos = this.new_xpos;
            this.cur_ypos = this.new_ypos;
        }
    }

    press_up(instance) {
        if (instance.button !== 0) return;
        if (this.isDown === true) {
            this.isDown = false

            if (this.drawing_zoom) {
                var zoom_rect = document.getElementById('divzoom').getBoundingClientRect();
                var svgbox = document.getElementById('svgbox');
                var svgbox_rect = svgbox.getBoundingClientRect();

                if (zoom_rect.width > 2 && zoom_rect.height > 2 && this.panZoomInstance) {
                    var sizes = this.panZoomInstance.getSizes();
                    var currentZoom = this.panZoomInstance.getZoom();
                    var currentPan = this.panZoomInstance.getPan();

                    var zoomCenterX = zoom_rect.left + zoom_rect.width / 2 - svgbox_rect.left;
                    var zoomCenterY = zoom_rect.top + zoom_rect.height / 2 - svgbox_rect.top;

                    var svgCenterX = (zoomCenterX - currentPan.x) / currentZoom;
                    var svgCenterY = (zoomCenterY - currentPan.y) / currentZoom;

                    var zoomFactor = Math.min(svgbox_rect.width / zoom_rect.width, svgbox_rect.height / zoom_rect.height);
                    var newZoom = currentZoom * zoomFactor;

                    this.panZoomInstance.zoom(newZoom);

                    var newPanX = svgbox_rect.width / 2 - svgCenterX * newZoom;
                    var newPanY = svgbox_rect.height / 2 - svgCenterY * newZoom;
                    this.panZoomInstance.pan({x: newPanX, y: newPanY});
                }

                this.drawing_zoom = false;
                this.zoomBox = {};
                $('#divzoom').hide();
                return;
            }

            if (this.diff < 10) {
                if (instance.target.getAttribute('class') === 'group' || instance.target.getAttribute('class') === 'node') {
                    if (instance.altKey) {
                        var bin_id = this.current_bin_id;
                        this.marknode(instance.target, bin_id);
                    } else {
                        this.nodeinfo(instance.target);
                    }
                }
            } else {
                if (instance.altKey) {
                    var bin_id = this.current_bin_id;

                    var min_svg_xpos = Math.min(this.svg_cur_xpos, this.svg_new_xpos);
                    var min_svg_ypos = Math.min(this.svg_cur_ypos, this.svg_new_ypos);
                    var width = Math.abs(this.svg_cur_xpos - this.svg_new_xpos);
                    var height = Math.abs(this.svg_cur_ypos - this.svg_new_ypos);

                    var nodes = document.querySelectorAll(".node")
                    for (var n of nodes) {
                        var bbox = n.getBBox();

                        if (
                            bbox.x + bbox.width >= min_svg_xpos &&
                            bbox.x <= min_svg_xpos + width &&
                            bbox.y + bbox.height >= min_svg_ypos &&
                            bbox.y <= min_svg_ypos + height
                        ) {
                            this.marknode(n, bin_id);
                        }
                    }

                    var groups = Object.keys(this.group_dict);
                    for(var g of groups) {
                        var group = this.group_dict[g]
                        for (var k of group) {
                            var node = document.getElementById(k);
                            var bbox = node.getBBox();
                            if (
                                bbox.x + bbox.width >= min_svg_xpos &&
                                bbox.x <= min_svg_xpos + width &&
                                bbox.y + bbox.height >= min_svg_ypos &&
                                bbox.y <= min_svg_ypos + height
                            ) {
                                var name = document.getElementById(g);
                                this.marknode(name, bin_id);
                                break
                            }
                        }
                    }
                }
            }

            if (this.selection) {
                this.selection.remove();
                this.selection = null;
            }
        }
    }

    update_bin() {
        this._suppress_bin_ring_draw = true;
        for (var bin_id of Object.keys(this.bin_dict)) {
            var nodes = this.bin_dict[bin_id];
            var updated_nodes = []
            for (var node of nodes) {
                var name = document.getElementById(node);

                if(name) {
                    if (name.getAttribute('class') == 'pseudo') {
                        var groups = Object.keys(this.group_dict);
                        for(var g of groups) {
                            var group = this.group_dict[g]
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
                    updated_nodes.push(...this.bin_group_dict[node])
                    delete this.bin_group_dict[node];
                };
            }

            this.bin_dict[bin_id] = updated_nodes
            for (var node of this.bin_dict[bin_id]) {

                this.bin_dict[bin_id] = this.bin_dict[bin_id].filter(item => item !== node);
                var name = document.getElementById(node);
                this.marknode(name, bin_id);

            }
        }
        this._suppress_bin_ring_draw = false;
    }

    draw_bin_rings() {
        if (!this._layout) return;
        const L = this._layout;

        // Collect per-bin data: runs, color, name, label midpoint
        this._bin_ring_data = [];

        for (const [bin_id, node_ids] of Object.entries(this.bin_dict)) {
            if (!node_ids || node_ids.length === 0) continue;

            const bin_color = $('#' + bin_id + '_color').attr('color') || '#000000';
            const raw_name = $('#' + bin_id + '_text').val() || bin_id;
            const bin_name = raw_name.replace(/_/g, ' ').replace(/^./, c => c.toUpperCase());

            // Gather all x-positions for nodes in this bin
            const x_positions = new Set();
            for (const nid of node_ids) {
                if (this.group_dict[nid]) {
                    for (const member of this.group_dict[nid]) {
                        if (this.data['nodes'][member]) {
                            x_positions.add(this.data['nodes'][member]['position'][0]);
                        }
                    }
                } else if (this.data['nodes'][nid]) {
                    x_positions.add(this.data['nodes'][nid]['position'][0]);
                }
            }
            if (x_positions.size === 0) continue;

            // Sort and group into contiguous runs
            const sorted_x = [...x_positions].sort((a, b) => a - b);
            const runs = [];
            let run_start = sorted_x[0], run_end = sorted_x[0];
            for (let i = 1; i < sorted_x.length; i++) {
                if (sorted_x[i] <= run_end + 1) {
                    run_end = sorted_x[i];
                } else {
                    runs.push([run_start, run_end]);
                    run_start = sorted_x[i];
                    run_end = sorted_x[i];
                }
            }
            runs.push([run_start, run_end]);

            this._bin_ring_data.push({ bin_color, bin_name, runs });
        }

        // Now render everything with zoom-aware positioning
        this._render_bin_visuals();
    }

    _render_bin_visuals() {
        // Locate (or create) the SVG group for bin visuals
        const svg_el = document.getElementById('result');
        if (!svg_el) return;
        const vp = svg_el.querySelector('.svg-pan-zoom_viewport') || svg_el;
        let grp = vp.querySelector('#bin-rings-group');
        if (!grp) {
            grp = document.createElementNS('http://www.w3.org/2000/svg', 'g');
            grp.setAttribute('id', 'bin-rings-group');
            const rlg = vp.querySelector('#region-labels-group');
            if (rlg) vp.insertBefore(grp, rlg);
            else vp.appendChild(grp);
        }
        while (grp.firstChild) grp.removeChild(grp.firstChild);

        if (!this._layout) return;
        const L = this._layout;

        // --- Compute the max region-label radius at the current zoom ---
        const realZoom = this.panZoomInstance ? this.panZoomInstance.getSizes().realZoom : 1;
        const DISTANCE = parseFloat($('#region_label_distance')[0].value) || 2;
        const region_font_svg = (parseFloat($('#region_label_size')[0].value) || 13) / Math.pow(realZoom, 0.6);

        // Compute the true radial extent of each visible region label
        // analytically.  Labels are horizontal (non-rotated) text centred
        // at radius r_center.  For such text at angle θ the outward radial
        // projection of the bounding-box corners is:
        //   r_center + (text_width/2)*|sin θ| + (font_size/2)*|cos θ|
        // This avoids getBBox() whose axis-aligned boxes overestimate the
        // radial extent for text placed at angles around the circle.
        let max_region_label_edge = 0;
        document.querySelectorAll('.region-label').forEach(el => {
            if (el.style.display === 'none') return;
            const outer_r = parseFloat(el.dataset.outerR);
            const r_center = outer_r + Math.max(region_font_svg * DISTANCE, outer_r * 0.005);
            if (el.dataset.layout === 'circular') {
                const angle_rad = parseFloat(el.dataset.angle) * Math.PI / 180;
                let text_w;
                try { text_w = el.getComputedTextLength(); } catch (_) { text_w = region_font_svg * 3; }
                const edge = r_center
                    + (text_w / 2) * Math.abs(Math.sin(angle_rad))
                    + (region_font_svg / 2) * Math.abs(Math.cos(angle_rad));
                if (edge > max_region_label_edge) max_region_label_edge = edge;
            } else {
                // Linear layout: radial extent is just the centre plus half the font height
                const edge = r_center + region_font_svg * 0.5;
                if (edge > max_region_label_edge) max_region_label_edge = edge;
            }
        });

        // If no region labels are visible, fall back to the content edge
        if (max_region_label_edge === 0) {
            let max_y = 0;
            for (const [nid, ndata] of Object.entries(this.nodes)) {
                if (ndata['position'][1] > max_y) max_y = ndata['position'][1];
            }
            max_region_label_edge = L.current_outer_stop + L.sum_middle_layer + max_y * L.node_distance_y;
        }

        // Ring sits beyond the outermost region label edge with a clear gap
        const ring_gap = region_font_svg * 0.8;
        const ring_thickness = L.node_distance_y * (parseFloat($('#bin_ring_height')[0].value) || 4);
        const ring_inner = max_region_label_edge + ring_gap;
        const ring_outer = ring_inner + ring_thickness;

        // Background fill starts from the innermost graph content
        const bg_inner_r = L.start_offset;

        // Cache radii so the region-label hover highlight can match bin pie size
        this._hover_bg_inner_r = bg_inner_r;
        this._hover_ring_outer = ring_outer;

        if (!this._bin_ring_data || this._bin_ring_data.length === 0) return;

        // Bin label font
        const bin_font_svg = (parseFloat($('#bin_label_size')[0].value) || 19.5) / Math.pow(realZoom, 0.6);
        // Label sits just beyond the ring
        const label_r = ring_outer + bin_font_svg * 0.4;

        // Ring and edge appearance
        const ring_opacity = parseFloat($('#bin_ring_opacity')[0].value) || 0.8;
        const show_edges = $('#flexbinedges').prop('checked');
        const edge_thickness = parseFloat($('#bin_edge_thickness')[0].value) || 4;
        const edge_color = $('#bin_edge_color').attr('color') || '#FFFFFF';
        const edge_opacity = parseFloat($('#bin_edge_opacity')[0].value) || 1;

        let clip_counter = 0;
        const deferred_ring_els = [];
        for (const { bin_color, bin_name, runs } of this._bin_ring_data) {
            for (const [rx_min, rx_max] of runs) {
                const run_x_mid = (rx_min + rx_max + 1) / 2;

                if (L.linear == 0) {
                    const a1_deg = L.theta * rx_min + L.start_angle;
                    const a2_deg = L.theta * (rx_max + 1) + L.start_angle;
                    const a1 = this.deg2rad(a1_deg);
                    const a2 = this.deg2rad(a2_deg);
                    const span_deg = a2_deg - a1_deg;
                    const large_arc = span_deg > 180 ? 1 : 0;

                    // Background fill wedge
                    const bg_d = this._arc_path(bg_inner_r, ring_outer, a1, a2, large_arc);
                    const bg_el = document.createElementNS('http://www.w3.org/2000/svg', 'path');
                    bg_el.setAttribute('d', bg_d);
                    bg_el.setAttribute('fill', bin_color);
                    bg_el.setAttribute('fill-opacity', '0.08');
                    if (show_edges) {
                        // Double stroke-width + clip to shape = inside-only stroke
                        const clip_id = 'bin-clip-' + (clip_counter++);
                        const clipPath = document.createElementNS('http://www.w3.org/2000/svg', 'clipPath');
                        clipPath.setAttribute('id', clip_id);
                        const clipUse = document.createElementNS('http://www.w3.org/2000/svg', 'path');
                        clipUse.setAttribute('d', bg_d);
                        clipPath.appendChild(clipUse);
                        grp.appendChild(clipPath);
                        bg_el.setAttribute('clip-path', 'url(#' + clip_id + ')');
                        bg_el.setAttribute('stroke', edge_color);
                        bg_el.setAttribute('stroke-width', edge_thickness * 2);
                        bg_el.setAttribute('stroke-opacity', edge_opacity);
                    } else {
                        bg_el.setAttribute('stroke-width', '0');
                    }
                    bg_el.setAttribute('pointer-events', 'none');
                    grp.appendChild(bg_el);

                    // Solid outer ring arc (deferred so it paints on top of all edges)
                    const ring_d = this._arc_path(ring_inner, ring_outer, a1, a2, large_arc);
                    const ring_el = document.createElementNS('http://www.w3.org/2000/svg', 'path');
                    ring_el.setAttribute('d', ring_d);
                    ring_el.setAttribute('fill', bin_color);
                    ring_el.setAttribute('fill-opacity', ring_opacity);
                    ring_el.setAttribute('stroke-width', '0');
                    ring_el.setAttribute('pointer-events', 'none');
                    deferred_ring_els.push(ring_el);
                } else {
                    const x1 = rx_min * L.node_distance_x;
                    const x2 = (rx_max + 1) * L.node_distance_x;

                    // Background fill
                    const bg_rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
                    bg_rect.setAttribute('x', Math.min(x1, x2));
                    bg_rect.setAttribute('y', -(ring_outer));
                    bg_rect.setAttribute('width', Math.abs(x2 - x1));
                    bg_rect.setAttribute('height', ring_outer - bg_inner_r);
                    bg_rect.setAttribute('fill', bin_color);
                    bg_rect.setAttribute('fill-opacity', '0.08');
                    if (show_edges) {
                        const clip_id = 'bin-clip-' + (clip_counter++);
                        const clipPath = document.createElementNS('http://www.w3.org/2000/svg', 'clipPath');
                        clipPath.setAttribute('id', clip_id);
                        const clipRect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
                        clipRect.setAttribute('x', Math.min(x1, x2));
                        clipRect.setAttribute('y', -(ring_outer));
                        clipRect.setAttribute('width', Math.abs(x2 - x1));
                        clipRect.setAttribute('height', ring_outer - bg_inner_r);
                        clipPath.appendChild(clipRect);
                        grp.appendChild(clipPath);
                        bg_rect.setAttribute('clip-path', 'url(#' + clip_id + ')');
                        bg_rect.setAttribute('stroke', edge_color);
                        bg_rect.setAttribute('stroke-width', edge_thickness * 2);
                        bg_rect.setAttribute('stroke-opacity', edge_opacity);
                    } else {
                        bg_rect.setAttribute('stroke-width', '0');
                    }
                    bg_rect.setAttribute('pointer-events', 'none');
                    grp.appendChild(bg_rect);

                    // Solid outer ring strip (deferred so it paints on top of all edges)
                    const rect_el = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
                    rect_el.setAttribute('x', Math.min(x1, x2));
                    rect_el.setAttribute('y', -(ring_outer));
                    rect_el.setAttribute('width', Math.abs(x2 - x1));
                    rect_el.setAttribute('height', ring_thickness);
                    rect_el.setAttribute('fill', bin_color);
                    rect_el.setAttribute('fill-opacity', ring_opacity);
                    rect_el.setAttribute('stroke-width', '0');
                    rect_el.setAttribute('pointer-events', 'none');
                    deferred_ring_els.push(rect_el);
                }

                // Bin label — one per run so split bins each get their own label
                if ($('#flexbinlabels').prop('checked')) {
                    const label_el = document.createElementNS('http://www.w3.org/2000/svg', 'text');
                    label_el.textContent = bin_name;
                    label_el.setAttribute('text-anchor', 'middle');
                    label_el.setAttribute('dominant-baseline', 'middle');
                    label_el.setAttribute('fill', bin_color);
                    label_el.setAttribute('opacity', '0.9');
                    label_el.setAttribute('pointer-events', 'none');
                    label_el.setAttribute('font-family', 'HelveticaNeue-CondensedBold, Helvetica Neue, Helvetica, sans-serif');
                    label_el.setAttribute('font-weight', 'bold');
                    label_el.setAttribute('font-size', bin_font_svg);

                    const orientation = $('#bin_label_orientation')[0].value || 'natural';

                    if (L.linear == 0) {
                        const angle_rad = this.deg2rad(L.theta * (run_x_mid + 0.5) + L.start_angle);
                        const sin_a = Math.sin(angle_rad);
                        const cos_a = Math.cos(angle_rad);
                        const lx = label_r * sin_a;
                        const ly = label_r * cos_a;

                        if (orientation === 'radial') {
                            // Rotate text along the radial direction; flip on the
                            // left side so it never reads upside-down.
                            let rot = Math.atan2(cos_a, sin_a) * 180 / Math.PI;
                            if (sin_a >= 0) {
                                label_el.setAttribute('text-anchor', 'start');
                            } else {
                                label_el.setAttribute('text-anchor', 'end');
                                rot += 180;
                            }
                            label_el.setAttribute('x', lx);
                            label_el.setAttribute('y', ly);
                            label_el.setAttribute('transform', `rotate(${rot}, ${lx}, ${ly})`);
                        } else {
                            // Natural: smart text-anchor based on position
                            if (sin_a > 0.25) label_el.setAttribute('text-anchor', 'start');
                            else if (sin_a < -0.25) label_el.setAttribute('text-anchor', 'end');
                            label_el.setAttribute('x', lx);
                            label_el.setAttribute('y', ly);
                        }
                    } else {
                        const lx = run_x_mid * L.node_distance_x;
                        const ly = -label_r;

                        if (orientation === 'radial') {
                            label_el.setAttribute('text-anchor', 'start');
                            label_el.setAttribute('x', lx);
                            label_el.setAttribute('y', ly);
                            label_el.setAttribute('transform', `rotate(-90, ${lx}, ${ly})`);
                        } else {
                            label_el.setAttribute('x', lx);
                            label_el.setAttribute('y', ly);
                        }
                    }
                    grp.appendChild(label_el);
                }
            }
        }

        // Append solid ring elements on top of all backgrounds and edges
        for (const el of deferred_ring_els) grp.appendChild(el);
    }

    // Helper: build an SVG arc-ring path between two radii at two angles
    _arc_path(r_inner, r_outer, a1, a2, large_arc) {
        const ox1 = r_outer * Math.sin(a1), oy1 = r_outer * Math.cos(a1);
        const ox2 = r_outer * Math.sin(a2), oy2 = r_outer * Math.cos(a2);
        const ix1 = r_inner * Math.sin(a1), iy1 = r_inner * Math.cos(a1);
        const ix2 = r_inner * Math.sin(a2), iy2 = r_inner * Math.cos(a2);
        return `M ${ix1} ${iy1} L ${ox1} ${oy1}
                A ${r_outer} ${r_outer} 0 ${large_arc} 0 ${ox2} ${oy2}
                L ${ix2} ${iy2}
                A ${r_inner} ${r_inner} 0 ${large_arc} 1 ${ix1} ${iy1} Z`;
    }

    fit_aspect() {
        this.panZoomInstance.resize();
        this.panZoomInstance.fit();
        this.panZoomInstance.center();
    }

    svg_download() {
        var svg_download = this.generate_svg(0);
        var blob = new Blob([svg_download[0].outerHTML]);
        var title = this.data['meta']['project_name']
        this.download_blob(blob, title + ".svg");
    }
    
    main_draw() {
        var svg_core = this.generate_svg();
        var genome_size = this.genomes.length;

        // TODO this is just a temporary fix the whole generate_svg() has to be rewritten in either
        // createElementNS or raw html string fashion. FML.
        // Save current pan/zoom so a redraw (e.g. color change) doesn't snap back.
        // Only restore if the drawing type hasn't changed — circular↔linear have incompatible coordinate systems.
        const currentLinear = $('#flexlinear').prop('checked');
        var savedPan = null, savedZoom = null;
        if (this.panZoomInstance !== null) {
            if (currentLinear === this._lastLinear) {
                savedPan  = this.panZoomInstance.getPan();
                savedZoom = this.panZoomInstance.getZoom();
            }
            this.panZoomInstance.destroy();
        }
        this._lastLinear = currentLinear;

        $('#svgbox').empty().html(svg_core[0].outerHTML);

        this.refresh_region_label_visibility = () => {
            if (!$('#flexregionlabels').prop('checked')) return;
            const realZoom = this.panZoomInstance.getSizes().realZoom;
            const TARGET_PX = parseFloat($('#region_label_size')[0].value) || 13;
            const MIN_WIDTH_PX = parseFloat($('#region_label_min_width')[0].value) || 80;
            // Exponent < 1 makes screen size grow with zoom (exponent 1 = constant screen
            // size; exponent 0 = scales with graph).  0.6 gives moderate growth: at 10x
            // zoom the label is ~2.5x its original screen size, feeling part of the graph.
            const font_size_svg = TARGET_PX / Math.pow(realZoom, 0.6);
            // r = content_edge + font_size * distance, where distance is user-controlled.
            // distance=1 means the label baseline sits flush; higher values push it further out.
            const DISTANCE = parseFloat($('#region_label_distance')[0].value) || 2;
            document.querySelectorAll('.region-label').forEach(el => {
                const screen_width = parseFloat(el.dataset.svgWidth) * realZoom;
                if (screen_width >= MIN_WIDTH_PX) {
                    const outer_r = parseFloat(el.dataset.outerR);
                    // The font-based term shrinks as realZoom grows, so at extreme zoom it
                    // can become smaller than the arc/stroke overshoot beyond outer_r, making
                    // the label appear inside the graph.  A 0.5% floor ensures the label
                    // clears the graph at high zoom without pushing it far off-screen.
                    const r = outer_r + Math.max(font_size_svg * DISTANCE, outer_r * 0.005);
                    if (el.dataset.layout === 'circular') {
                        const angle_rad = parseFloat(el.dataset.angle) * Math.PI / 180;
                        el.setAttribute('x', r * Math.sin(angle_rad));
                        el.setAttribute('y', r * Math.cos(angle_rad));
                    } else {
                        el.setAttribute('x', el.dataset.labelX);
                        el.setAttribute('y', -r);
                    }
                    el.setAttribute('font-size', font_size_svg);
                    el.style.display = '';
                } else {
                    el.style.display = 'none';
                }
            });
        };

        this.panZoomInstance = svgPanZoom('#result', {
            zoomEnabled: true,
            panEnabled: false,
            controlIconsEnabled: false,
            minZoom: 0.1,
            maxZoom: 100,
            onZoom: () => { this.refresh_region_label_visibility(); this._render_bin_visuals(); }
        });

        // Restore pan/zoom from before the redraw (e.g. after a color change).
        if (savedZoom !== null) {
            this.panZoomInstance.zoom(savedZoom);
            this.panZoomInstance.pan(savedPan);
        }

        // Auto-size region labels once on first load.
        //
        // Goal: labels should appear at ~40 screen-pixels when the whole graph fills the
        // viewport (the typical initial view after panning/zooming to fit).
        //
        // With font_size_svg = TARGET_PX / realZoom^0.6, screen size = TARGET_PX * realZoom^0.4.
        // At "graph fills viewport" zoom for a circular layout:
        //   realZoom_fill ≈ container_width / (2 * max_outer_r)
        //   screen_px = TARGET_PX * (container_width / (2*max_outer_r))^0.4 = 16
        //   TARGET_PX = 16 * (2*max_outer_r / container_width)^0.4
        //
        // Using max_outer_r (the graph's total SVG extent) rather than per-region widths
        // ensures a sensible default even for graphs with many narrow regions (e.g. SAR11).
        if (!this._region_label_size_initialized) {
            this._region_label_size_initialized = true;
            const labels = document.querySelectorAll('.region-label');
            if (labels.length > 0) {
                const max_outer_r = Math.max(...Array.from(labels)
                    .map(el => parseFloat(el.dataset.outerR))
                    .filter(r => r > 0));
                const container_width = document.getElementById('result').clientWidth || window.innerWidth;
                const auto_px = Math.round(16 * Math.pow(2 * max_outer_r / container_width, 0.4));
                $('#region_label_size')[0].value = Math.max(8, Math.min(200, auto_px));
                $('#bin_label_size')[0].value = Math.max(12, Math.min(300, Math.round(auto_px * 1.5)));
            }
        }

        this.refresh_region_label_visibility();

        // Hover highlight: draw a filled wedge/rectangle showing the region extent
        const svg_el = document.getElementById('result');
        const get_viewport = () => svg_el.querySelector('.svg-pan-zoom_viewport') || svg_el;

        const remove_highlight = () => {
            const el = document.getElementById('region-hover-highlight');
            if (el) el.remove();
        };

        document.querySelectorAll('.region-label').forEach(el => {
            el.style.cursor = 'pointer';

            el.addEventListener('click', (event) => {
                event.stopPropagation();
                const rid = el.dataset.regionId;
                this.show_region_context_menu(event, rid);
            });

            el.addEventListener('mouseover', () => {
                remove_highlight();

                const layout = el.dataset.layout;
                // Use the same radii as the bin pie so hover preview matches bin visuals
                const outer_r = this._hover_ring_outer || parseFloat(el.dataset.outerR) || 0;
                const inner_r = this._hover_bg_inner_r || parseFloat(el.dataset.innerR) || 0;
                // Active bin color
                const bin_color = $('#' + this.current_bin_id + '_color').attr('color') || '#000000';

                let path_d;
                if (layout === 'circular') {
                    const a1 = parseFloat(el.dataset.angleStart) * Math.PI / 180;
                    const a2 = parseFloat(el.dataset.angleEnd)   * Math.PI / 180;
                    const span = parseFloat(el.dataset.angleEnd) - parseFloat(el.dataset.angleStart);
                    const large_arc = span > 180 ? 1 : 0;
                    if (inner_r <= 0) {
                        const ox1 = outer_r * Math.sin(a1), oy1 = outer_r * Math.cos(a1);
                        const ox2 = outer_r * Math.sin(a2), oy2 = outer_r * Math.cos(a2);
                        path_d = `M 0 0 L ${ox1} ${oy1} A ${outer_r} ${outer_r} 0 ${large_arc} 0 ${ox2} ${oy2} Z`;
                    } else {
                        path_d = this._arc_path(inner_r, outer_r, a1, a2, large_arc);
                    }
                } else {
                    const xs = parseFloat(el.dataset.xStart);
                    const xe = parseFloat(el.dataset.xEnd);
                    path_d = `M ${xs} ${-inner_r} L ${xe} ${-inner_r} L ${xe} ${-outer_r} L ${xs} ${-outer_r} Z`;
                }

                const highlight = document.createElementNS('http://www.w3.org/2000/svg', 'path');
                highlight.setAttribute('id', 'region-hover-highlight');
                highlight.setAttribute('d', path_d);
                highlight.setAttribute('fill', bin_color);
                highlight.setAttribute('fill-opacity', '0.08');
                highlight.setAttribute('stroke-width', '0');
                highlight.setAttribute('pointer-events', 'none');
                // Insert at the bottom of the viewport so it sits behind all graph elements
                const vp = get_viewport();
                vp.insertBefore(highlight, vp.firstChild);
            });

            el.addEventListener('mouseout', remove_highlight);
        });

        var elements = document.querySelectorAll(".node, .group");
        for (var element of elements) {

            if (element.getAttribute("class") == 'group'){
                var element_id = this.group_dict[element.getAttribute("id")][0]
                var node_id = element.getAttribute("id")
            } else {
                var element_id = element.getAttribute("id")
                var node_id = this.data['nodes'][element.getAttribute("id")]["gene_cluster"]
            }

            element.addEventListener('contextmenu', (ev) => {
                this.show_node_context_menu(ev, ev.currentTarget);
            });

            tippy(element, {
                content: '<strong>' + node_id + '</strong>' + '<br />',
                allowHTML: true,
                onHide: this.hide_tippy,
                onShow: this.show_tippy,
                arrow: false,
                duration: 0,
                followCursor: true,
                theme: "light",
            });
        };

        this.update_bin();
        this.draw_bin_rings();
    }

    rerun_JSON(new_data) {
        $.ajax({
            url: "/pangraph/rerun_pangraph_json_data",
            type: "POST",
            async: false,
            data: JSON.stringify(new_data),
            contentType: "application/json",
            dataType: "json",
            error: function(){
                toastr.error('Failed to update pan-graph data.', 'Error');
            },
            success: function(json){}
        });
    }

    initialize_JSON() {
        $.ajax({
            url: "/pangraph/initial_pangraph_json_data",
            type: "GET",
            cache: false,
            contentType: "application/json",
            dataType: "json",
            success: (data) => {
                if (data['status'] != 0) {
                    toastr.error(data['message'], 'Initialization error', { 'timeOut': '0', 'extendedTimeOut': '0' });
                    return;
                }
                this.data = data['data'];
                this.session_id = data['session_id'];
                setTimeout(() => this.checkBackgroundProcess(), 5000);
                this.initialize_variables();
                this.initialize_user_interface();
                this.set_UI_settings();
            },
            error: (err) => {
                toastr.error('Could not reach the server during initialization.', 'Initialization error', { 'timeOut': '0', 'extendedTimeOut': '0' });
            }
        });
    }

    checkBackgroundProcess() {
        const errorMessage = "It seems the server that's been serving this page is no longer accessible. " +
                             "You may lose your unsaved changes in this window.";
        $.ajax({
            type: 'GET',
            cache: false,
            url: '/pangraph/session_id',
            success: (data) => {
                if (data != this.session_id) {
                    this.server_offline = true;
                    toastr.error(errorMessage, "", { 'timeOut': '0', 'extendedTimeOut': '0' });
                } else {
                    setTimeout(() => this.checkBackgroundProcess(), 5000);
                }
            },
            error: () => {
                this.server_offline = true;
                toastr.error(errorMessage, "", { 'timeOut': '0', 'extendedTimeOut': '0' });
            }
        });
    }

    marknode(element, bin_id) {
        
        var bin_color = $('#' + bin_id + '_color').attr('color')
        var id = element.id;
        var current = ''
        
        var bin_keys = Object.keys(this.bin_dict)
        for (var key of bin_keys) {
            if (this.bin_dict[key].includes(id)) {
                current = key;
                break;
            }
        }
        
        var core_color = $('#core_color').attr('color');
        var paralog_color = $('#paralog_color').attr('color');
        var singleton_color = $('#singleton_color').attr('color');
        var accessory_color = $('#accessory_color').attr('color');
        var rearranged_color = $('#rearranged_color').attr('color');
        var trna_color = $('#trna_color').attr('color');

        var genome_size = this.genomes.length
        
        if ($('#flexsaturation').prop('checked') == true){
            var saturation = 1
        } else {
            var saturation = 0
        }
        
        if (element.getAttribute('class') == 'group') {
            var group = this.group_dict[id]
            var node_name = group[0]
            var node = this.data['nodes'][node_name];
            var node_type = node['type']
            var genome = Object.keys(node['gene_calls']).length;
        } else if (element.getAttribute('class') == 'node') {
            var node = this.data['nodes'][id];
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
        
        if (current === bin_id) {
        
            if (saturation == 1){
                element.setAttribute("fill", this.lighter_color('#ffffff', node_color, genome / genome_size))
            } else {
                element.setAttribute("fill", node_color)
            }
            this.bin_dict[bin_id] = this.bin_dict[bin_id].filter(item => item !== id)
            $('#' + bin_id + '_value')[0].value = this.bin_dict[bin_id].length

            if (element.getAttribute('class') == 'group') {
                delete this.bin_group_dict[id];
            }
            
        } else if (current === '') {
        
            if (saturation == 1){
                element.setAttribute("fill", this.lighter_color('#ffffff', bin_color, genome / genome_size))
            } else {
                element.setAttribute("fill", bin_color)
            }
            this.bin_dict[bin_id].push(id)
            $('#' + bin_id + '_value')[0].value = this.bin_dict[bin_id].length

            if (element.getAttribute('class') == 'group') {
                this.bin_group_dict[id] = group
            }
        
        } else {
        
            if (saturation == 1){
                element.setAttribute("fill", this.lighter_color('#ffffff', bin_color, genome / genome_size))
            } else {
                element.setAttribute("fill", bin_color)
            }
            
            this.bin_dict[current] = this.bin_dict[current].filter(item => item !== id)
            this.bin_dict[bin_id].push(id)
            
            $('#' + bin_id + '_value')[0].value = this.bin_dict[bin_id].length
            $('#' + current + '_value')[0].value = this.bin_dict[current].length
        }

        if (!this._suppress_bin_ring_draw) this.draw_bin_rings();
    }

    
    switch_bin(instance) {
        var bin_id = instance.currentTarget.id.replace("_radio", "");
        this.current_bin_id = bin_id
    }

    add_info_to_bin() {
    
        var bin_id = this.current_bin_id;
        var target = document.getElementById('node_basics_table');
        var gc_context = target.getAttribute('gc_context');
        var element = document.getElementById(gc_context)
        
        this.marknode(element, bin_id);
    }
    
    delete_bin(bin_id) {
        this._suppress_bin_ring_draw = true;
        for (var node of this.bin_dict[bin_id]) {
            var element = document.getElementById(node);
            this.marknode(element, bin_id);
        }
        this._suppress_bin_ring_draw = false;

        $("#" + bin_id + "_grid").remove();
        delete this.bin_dict[bin_id];
        this.draw_bin_rings();

        if (Object.keys(this.bin_dict).length !== 0) {
            var next_bin_id = Object.keys(this.bin_dict)[0];
            $('#' + next_bin_id + '_radio').click();
        } else {
            this.add_bin();
        }
    }

    delete_all_bins() {
        this._suppress_bin_ring_draw = true;
        for (var bin_id of Object.keys(this.bin_dict)) {
            for (var node of this.bin_dict[bin_id]) {
                var element = document.getElementById(node);
                this.marknode(element, bin_id);
            }
        }
        this._suppress_bin_ring_draw = false;

        $('#bingrid').empty();
        this.bin_dict = {};
        this.bin_group_dict = {};
        this.current_bin_number = 0;
        this.add_bin();
        this.draw_bin_rings();
    }

    add_bin() {
        this.current_bin_number += 1;
        this.current_bin_id = "bin_" + this.current_bin_number;
        const n = this.current_bin_number;
        const new_color = randomColor({luminosity: 'dark'});

        const $row = $(`<tr class="bin-row" id="bin_${n}_grid">
            <td><input type="radio" name="binradio" id="bin_${n}_radio" bin_id="bin_${n}" checked></td>
            <td><div class="pangraph-colorpicker" id="bin_${n}_color" color="${new_color}" style="background-color: ${new_color}; width: 30px; height: 22px; cursor: pointer; border: 1px solid #ccc;"></div></td>
            <td><input type="text" class="form-control form-control-sm p-0 border-0" style="background-color: #e9ecef;" value="Bin_${n}" id="bin_${n}_text"></td>
            <td><input type="button" class="form-control form-control-sm p-0 border-0 bin-count-btn" id="bin_${n}_value" value=0 title="Click for functions summary"></td>
            <td><center><span class="default-bin-icon bi bi-trash-fill fa-lg" aria-hidden="true" title="Delete this bin" onclick="pgui.delete_bin('bin_${n}');"></span></center></td>
        </tr>`);

        $row.find(`#bin_${n}_value`).on('click', ((bin_id) => () => { this.show_bin_functions(bin_id); })(`bin_${n}`));
        $('#bingrid').append($row);

        $(`#bin_${n}_radio`).on("click", this.switch_bin);
        this._init_bin_colorpicker(`bin_${n}`);
        this.bin_dict[`bin_${n}`] = [];
    }
    
    initialize_variables() {
        this.nodes = this.data['nodes'];
        this.edges = this.data['edges'];
        // this.layers = this.data['meta']['layers'].filter(item => item !== 'backbone');
        this.layers = this.data['meta']['layers'];
        this.genomes = this.data['meta']['genome_names'];
        this.functional_annotation_sources_available = this.data['meta']['gene_function_sources'];

        this.group_dict = {}
        this.synteny = {}

        this.global_x = 0
        this.global_y = 0
        
        for(var g in this.genomes) {
            this.synteny[this.genomes[g]] = {}
        }
        
        for(var n in this.nodes) {
            var node = this.nodes[n];
            for(var [layer, value] of Object.entries(node["layer"])) {
                if (layer in this.layers_min) {
                    this.layers_min[layer] = value < this.layers_min[layer] ? value : this.layers_min[layer];
                } else {
                    this.layers_min[layer] = value;
                }

                if (layer in this.layers_max) {
                    this.layers_max[layer] = value > this.layers_max[layer] ? value : this.layers_max[layer];
                } else {
                    this.layers_max[layer] = value;
                }        
            }

            var group = node["group"];
            var x = node["position"][0];
            var y = node["position"][1];

            if (group) {
                if (group in this.group_dict) {
                    this.group_dict[group].push(n);
                } else {
                    this.group_dict[group] = [n];
                } 
            }
            this.global_x = x < this.global_x ? this.global_x : x;
            this.global_y = y < this.global_y ? this.global_y : y;

            for(var [genome, synteny_position] of Object.entries(node["synteny"])) {
                this.synteny[genome][synteny_position] = n
            }
        }
        
        for(var e in this.edges) {
            var edge = this.edges[e];
            var route = edge['route'];
            if (route.length > 0 && edge['active'] == true) {
                for (var b in route) {
                    var x = route[b][0];
                    var y = route[b][1];
                    this.global_y = y < this.global_y ? this.global_y : y;
                }
            }
        }
    }

    // Initialize a colpick color picker on the given selector (e.g. '#core_color').
    // Updates the swatch color live while the picker is open, then redraws when closed.
    _init_colorpicker(selector) {
        $(selector).colpick({
            layout: 'hex',
            submit: 0,
            colorScheme: 'light',
            onChange: (hsb, hex, rgb, el, bySetColor) => {
                $(el).css('background-color', '#' + hex);
                $(el).attr('color', '#' + hex);
            },
            onHide: () => {
                this.main_draw();
            }
        });
    }

    // Initialize a colpick color picker for a bin swatch.
    // Updates the swatch live, then re-marks all nodes in the bin when closed.
    _init_bin_colorpicker(bin_id) {
        const selector = '#' + bin_id + '_color';
        $(selector).colpick({
            layout: 'hex',
            submit: 0,
            colorScheme: 'light',
            onChange: (hsb, hex, rgb, el, bySetColor) => {
                $(el).css('background-color', '#' + hex);
                $(el).attr('color', '#' + hex);
            },
            onHide: () => {
                this._suppress_bin_ring_draw = true;
                const nodes = [...(this.bin_dict[bin_id] || [])];
                for (var node of nodes) {
                    this.bin_dict[bin_id] = this.bin_dict[bin_id].filter(item => item !== node);
                    var element = document.getElementById(node);
                    this.marknode(element, bin_id);
                }
                this._suppress_bin_ring_draw = false;
                this.draw_bin_rings();
            }
        });
    }

    // Initialize colpick on static (HTML-declared) color pickers.
    // Per-genome pickers are initialized individually in the genome loop.
    initialize_colorpickers() {
        const genomeColorIds = new Set(this.genomes);
        document.querySelectorAll('.pangraph-colorpicker').forEach(el => {
            if (!genomeColorIds.has(el.id)) {
                this._init_colorpicker('#' + el.id);
            }
        });
    }

    set_UI_settings() {

        var genome_order = []

        for (var [setting, value] of Object.entries(this.data['states'])) {
            const el = $('#' + setting);
            if (el.hasClass('pangraph-colorpicker')) {
                // colpick element: set background and color attribute
                const hex = value.replace('#', '');
                el.css('background-color', value).attr('color', value);
                el.colpickSetColor(hex);
            } else if (typeof value === 'number') {
                el[0].value = value;
            } else if (value == true || value == false) {
                el.prop('checked', value);
            } else {
                el[0].value = value;
            }

            if (this.genomes.includes(setting)) {
                genome_order.push(setting)
            }
        }

        var container = document.getElementById('genomecolors');

        genome_order.forEach(id => {
            var element = document.getElementById(id + '_row');
            if (element) {
                container.appendChild(element);
            }
        });

        // If max edge length filter was never set (legacy -1 default), apply the new default.
        if (parseInt($('#maxlength')[0].value) === -1) {
            $('#maxlength')[0].value = 1000;
            $('#flexmaxlength').prop('checked', true);
        }

        for(var [layer, max_value] of Object.entries(this.layers_max)) {
            $('#' + layer + '_max')[0].value = max_value;
        }

        for(var [layer, min_value] of Object.entries(this.layers_min)) {
            $('#' + layer + '_min')[0].value = min_value;
        }

        const isLinear = $('#flexlinear').prop('checked');
        $('#drawing_type_select').val(isLinear ? 'linear' : 'circular');
        $('#radius_row').toggle(!isLinear);
    }

    initialize_user_interface() {
        
        $('#RightOffcanvasBodyTop').append(
            $('<tr>').append(
                $('<td class="col-4">').append(
                    'Synteny gene cluster'
                )
            ).append(
                $('<td class="col-8 text-end" id="number_sgc">').append(
                    ''
                )
            )
        );

        $('#RightOffcanvasBodyTop').append(
            $('<tr>').append(
                $('<td class="col-4">').append(
                    'Gene cluster'
                )
            ).append(
                $('<td class="col-8 text-end" id="number_gc">').append(
                    ''
                )
            )
        );

        $('#RightOffcanvasBodyTop').append(
            $('<tr>').append(
                $('<td class="col-4">').append(
                    'SynGC type'
                )
            ).append(
                $('<td class="col-8 text-end" id="number_type">').append(
                    ''
                )
            )
        );

        $('#RightOffcanvasBodyTop').append(
            $('<tr>').append(
                $('<td class="col-4">').append(
                    'Position'
                )
            ).append(
                $('<td class="col-8 text-end" id="number_position">').append(
                    ''
                )
            )
        );
        
        for (var layer of this.layers) {
            // if ($('#flex' + layer + '').length == 0) {
                var element = $('<div class="col-12 d-flex mb-1"></div>').append(
                    $('<div class="col-1 d-flex align-items-center"></div>').append(
                        $('<div class="form-switch d-flex"></div>').append(
                            $('<input class="" type="checkbox" id="flex' + layer + '" name="flex' + layer + '" aria-label="..." data-toggle="tooltip" data-placement="top" title="Tooltip">')
                        )
                    )
                ).append(
                    $('<div class="col-5 d-flex align-items-center"></div>').append(
                        layer
                    )
                ).append(
                    $('<div class="d-flex col-2"></div>').append(
                        $('<input type="text" class="form-control float-end text-end flex-fill p-0 border-0" style= "background-color: #e9ecef;" id="' + layer + '_min" name="' + layer + '_min" value=0 aria-label="..." data-toggle="tooltip" data-placement="top" title="Choose your color">')
                    )
                ).append(
                    $('<div class="d-flex col-2"></div>').append(
                        $('<input type="text" class="form-control float-end text-end flex-fill p-0 border-0" style= "background-color: #e9ecef;" id="' + layer + '_max" name="' + layer + '_max" value=0 aria-label="..." data-toggle="tooltip" data-placement="top" title="Choose your color">')
                    )
                ).append(
                    $('<div class="d-flex col-2"></div>').append(
                        $('<input type="text" class="form-control float-end text-end flex-fill p-0 border-0" style= "background-color: #e9ecef;" id="' + layer + '" name="' + layer + '" value=0 aria-label="..." data-toggle="tooltip" data-placement="top" title="Choose your color">')
                    )
                );

                $('#RightOffcanvasBodyTop').append(
                    $('<tr>').append(
                        $('<td class="col-4">').append(
                            layer
                        )
                    ).append(
                        $('<td class="col-8 text-end" id="number_' + layer + '">').append(
                            0
                        )
                    )
                );
                
                $('#local_layers').append(element);
            // }
        }

        $('#title-panel-first-line').text(this.data['meta']['project_name']);
        $('#title-panel-second-line').text('Pangraph Detail');
        
        // if (!$('#genomecolors').children().length) {
        for (var genome of this.genomes) {  
            $('#genomecolors').append(
                $('<div class="col-12 d-flex mb-1" id="' + genome + '_row">').append(
                    $('<div class="col-1 d-flex align-items-center">').append(
                        $('<div class="form-switch d-flex">').append(
                            $('<input class="" type="checkbox" id="flex' + genome + '" name="' + genome + '" aria-label="..." data-bs-toggle="tooltip" data-bs-placement="top" title="Tooltip on top">')
                        )
                    )
                ).append(
                    $('<div class="col-8 d-flex align-items-center">').append(
                        genome
                    )
                ).append(
                    $('<div class="col-1 d-flex align-items-center">').append(
                        $('<i class="user-handle bi bi-arrows-expand"></i>')
                    )
                ).append(
                    $('<div class="d-flex col-2 align-items-center">').append(
                        $('<div class="pangraph-colorpicker" id="' + genome + '" color="#000000" style="background-color: #000000; width: 100%; height: 22px; cursor: pointer; border: 1px solid #ccc;"></div>')
                    )
                )
            );
            this._init_colorpicker('#' + genome);
    
            $('#RightOffcanvasBodyTop').append(
                $('<tr>').append(
                    $('<td class="col-4">').append(
                        genome
                    )
                ).append(
                    $('<td class="col-8 text-end" id="number_' + genome + '">').append(
                        0
                    )
                )
            );
    
            // if ($('#flex' + genome + 'layer').length == 0) {
            var element = $('<div class="col-12 d-flex mb-1"></div>').append(
                $('<div class="col-1 d-flex align-items-center"></div>').append(
                    $('<div class="form-switch d-flex"></div>').append(
                        $('<input class="" type="checkbox" id="flex' + genome + 'layer" name="flex' + genome + 'layer" aria-label="..." data-toggle="tooltip" data-placement="top" title="Tooltip on top">')
                    )
                )
            ).append(
                $('<div class="col-9 d-flex align-items-center"></div>').append(
                    genome
                )
            ).append(
                $('<div class="d-flex col-2"></div>').append(
                    $('<input type="text" class="form-control float-end text-end flex-fill p-0 border-0" style= "background-color: #e9ecef;" id="' + genome + 'layer" name="' + genome + 'layer" value=0 aria-label="..." data-toggle="tooltip" data-placement="top" title="Choose your color">')
                )
            );

            $('#genome_tracks').append(element);
        }
        // }
        // }
        
        // $('#searchSources').empty();
        for (var annotation_source of this.functional_annotation_sources_available){
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
            );

            $('#RightOffcanvasBodyTop').append(
                $('<tr>').append(
                    $('<td class="col-4">').append(
                        annotation_source
                    )
                ).append(
                    $('<td class="col-8 text-end" id="number_' + annotation_source + '">').append(
                        ''
                    )
                )
            );
        }
        
        // $('#expressiondrop').empty();
        $('#expressiondrop').append($('<option value="Choose item">Choose item</option>'));
        $('#expressiondrop').append($('<option value="Name">Name</option>'));
        $('#expressiondrop').append($('<option value="Position">Position</option>'));
        
        // $('#filter').empty();
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
        );
        
        for (var layer of this.layers) {
            $('#expressiondrop').append($('<option value="' + layer + '">' + layer + '</option>'));
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
                );
        }

        // Initialize main values based on standard settings
        
        // if ($('#flexlinear').prop('checked') == false) {
        //     $('#distx').prop('disabled', true);
        //     $('#inner').prop('disabled', false);
        // } else {
        //     $('#distx').prop('disabled', false);
        //     $('#inner').prop('disabled', true);
        // }
        
        // if ($('#flexarrow').prop('checked') == false) {
        //     $('#arrow').prop('disabled', true);
        // } else {
        //     $('#arrow').prop('disabled', false);
        // }
        
        // if ($('#flexbackbone').prop('checked') == false) {
        //     $('#backbone').prop('disabled', true);
        // } else {
        //     $('#backbone').prop('disabled', false);
        // }
        
        // if ($('#flexcondtr').prop('checked') == false) {
        //     $('#condtr').prop('disabled', true);
        // } else {
        //     $('#condtr').prop('disabled', false);
        // }
        
        // if ($('#flexmaxlength').prop('checked') == false) {
        //     $('#maxlength').prop('disabled', true);
        // } else {
        //     $('#maxlength').prop('disabled', false);
        // }
        
        // if ($('#flexgroupcompress').prop('checked') == false) {
        //     $('#groupcompress').prop('disabled', true);
        // } else {
        //     $('#groupcompress').prop('disabled', false);
        // }

        // Initialize main buttons with anonymous functions
        
        $('#flexlinear').change(function() {
            if ($(this).prop('checked') == true){
                $('#distx').prop('disabled', false);
                $('#inner').prop('disabled', true);
            } else {
                $('#distx').prop('disabled', true);
                $('#inner').prop('disabled', false);
            }
        })

        $('#drawing_type_select').on('change', () => {
            const isLinear = $('#drawing_type_select').val() === 'linear';
            $('#flexlinear').prop('checked', isLinear).trigger('change');
            $('#radius_row').toggle(!isLinear);
        });

        $('#flexglobalbackbone').change(function() {
            if ($(this).prop('checked') == true){
                $('#globalbackbone')[0].value = 100;
                $('#globalbackbone').prop('disabled', false);
            } else {
                $('#globalbackbone')[0].value = 0;
                $('#globalbackbone').prop('disabled', true);
            }
        })
        
        $('#flexarrow').change(function() {
            if ($(this).prop('checked') == true){
                $('#arrow')[0].value = 100;
                $('#arrow').prop('disabled', false);
            } else {
                $('#arrow')[0].value = 0;
                $('#arrow').prop('disabled', true);
            }
        })
        
        $('#flexcondtr').change(function() {
            if ($(this).prop('checked') == true){
                $('#condtr')[0].value = 2;
                $('#condtr').prop('disabled', false);
            } else {
                $('#condtr')[0].value = -1;
                $('#condtr').prop('disabled', true);
            }
        })
        
        $('#flexmaxlength').change(function() {
            if ($(this).prop('checked') == true){
                $('#maxlength')[0].value = 1000;
                $('#maxlength').prop('disabled', false);
            } else {
                $('#maxlength')[0].value = -1;
                $('#maxlength').prop('disabled', true);
            }
        })
        
        $('#flexgroupcompress').change(function() {
            if ($(this).prop('checked') == true){
                $('#groupcompress')[0].value = 0.0;
                $('#groupcompress').prop('disabled', false);
            } else {
                $('#groupcompress')[0].value = 1.0;
                $('#groupcompress').prop('disabled', true);
            }
        })
        
        for (var layer_name of this.layers) {
            $('#flex' + layer_name).change(function() {
                var entry = $(this)[0].id.replace("flex", "");
                if ($(this).prop('checked') == true){
                    $('#' + entry)[0].value = 25;
                    $('#' + entry).prop('disabled', false);
                } else {
                    $('#' + entry)[0].value = 0;
                    $('#' + entry).prop('disabled', true);
                }
            })
        }
        
        for (var genome of this.genomes) {
            $('#flex' + genome + 'layer').change(function() {
                var entry = $(this)[0].id.replace("flex", "");
                if ($(this).prop('checked') == true){
                    $('#' + entry)[0].value = 50;
                    $('#' + entry).prop('disabled', false);
                } else {
                    $('#' + entry)[0].value = 0;
                    $('#' + entry).prop('disabled', true);
                }
            })
        }

        // Initialize main buttons with "this" bound functions
        
        // Assign a random color to bin_1 and initialize its colorpicker
        const bin1_color = randomColor({luminosity: 'dark'});
        $('#bin_1_color').css('background-color', bin1_color).attr('color', bin1_color);
        this._init_bin_colorpicker('bin_1');
        $('#bin_1_radio').on("click", this.switch_bin)
        $('#bin_1_value').on("click", () => this.show_bin_functions('bin_1'))

        // Replace spaces with underscores in bin name inputs (delegated for dynamic bins too)
        $('#bingrid').on('input', 'input[id$="_text"]', function () {
            const el = this;
            const pos = el.selectionStart;
            el.value = el.value.replace(/ /g, '_');
            el.setSelectionRange(pos, pos);
        });
        // Update bin labels live when a bin name is changed
        $('#bingrid').on('change', 'input[id$="_text"]', () => this.draw_bin_rings());

        $('#flextree').on("change", this.flextree_change)
        $('#flexsaturation').on("change", () => this.main_draw())
        $('#flexregionlabels').on("change", () => this.main_draw())
        $('#region_label_size, #region_label_min_width, #region_label_distance').on("change", () => {
            if (this.panZoomInstance && $('#flexregionlabels').prop('checked')) {
                this.refresh_region_label_visibility();
            }
        })
        $('#flexbinlabels, #bin_label_orientation, #bin_label_size, #bin_ring_height, #bin_ring_opacity, #flexbinedges, #bin_edge_thickness, #bin_edge_opacity').on("change", () => this._render_bin_visuals())

        $('#binadd').on("click", this.add_bin);
        $('#binremove').on("click", this.delete_all_bins);
        $('#redraw').on("click", this.start_draw);
        $('#fit').on('click', this.fit_aspect);
        $('#svgDownload').on('click', this.svg_download);
        $('#genome_tracks_select_all').on('click', () => {
            $('#genome_tracks input[type="checkbox"]').prop('checked', true).trigger('change');
        })
        $('#genome_tracks_unselect_all').on('click', () => {
            $('#genome_tracks input[type="checkbox"]').prop('checked', false).trigger('change');
        })
        $('#genomes_select_all').on('click', () => {
            $('#genomecolors input[type="checkbox"]').prop('checked', true).trigger('change');
        })
        $('#genomes_unselect_all').on('click', () => {
            $('#genomecolors input[type="checkbox"]').prop('checked', false).trigger('change');
        })
        $('#svgbox').on('mousedown', this.press_down)
        $('#svgbox').on('mousemove', this.press_move)
        $('#svgbox').on('mouseup', this.press_up)
        $('#svgbox').on('mouseleave', this.press_up)
        $('#AddBin').on("click", this.add_info_to_bin);
        $('#AlignmentDownload').on("click", this.alignment_download);
        $('#InfoDownload').on("click", this.info_download);

        $('#stateload').on("click", () => this.state_modal('load'));
        $('#statesave').on("click", () => this.state_modal('save'));
        $('#binload').on("click", this.load_bin);
        $('#binsave').on("click", this.save_bin);

        $('#savestatebutton').on("click", this.save_state);
        $('#loadstatebutton').on("click", this.load_state);
        
        sortable('#genomecolors', {
            forcePlaceholderSize: true,
            handle: '.user-handle',
            items: 'div'
        });
        this.settings_dict['condtr'] = JSON.parse(JSON.stringify(this.data['states']['condtr']))
        this.settings_dict['maxlength'] = JSON.parse(JSON.stringify(this.data['states']['maxlength']))
        this.settings_dict['groupcompress'] = JSON.parse(JSON.stringify(this.data['states']['groupcompress']))
        this.settings_dict['state'] = JSON.parse(JSON.stringify(this.data['meta']['state']))

        // Delegated handlers for amino acid conservation checkboxes in the alignment modal
        $(document).on('change', '.pangraph-aa-checkbox', () => this.recolor_alignment());
        $(document).on('click', '#pangraph-aa-check-all', () => {
            document.querySelectorAll('.pangraph-aa-checkbox').forEach(cb => { cb.checked = true; });
            this.recolor_alignment();
        });
        $(document).on('click', '#pangraph-aa-uncheck-all', () => {
            document.querySelectorAll('.pangraph-aa-checkbox').forEach(cb => { cb.checked = false; });
            this.recolor_alignment();
        });

        // Keyboard shortcuts: D = Draw, S = toggle settings panel
        document.body.addEventListener('keydown', (ev) => {
            if ((/^(?:input|select|textarea|button)$/i).test(ev.target.nodeName)) return;
            if (ev.keyCode === 68) this.start_draw();          // D
            if (ev.keyCode === 83) toggleLeftPanel();           // S
        });

        // Initialize colpick on all static color pickers (per-genome pickers are
        // initialized individually as they are created in the genome loop above).
        this.initialize_colorpickers();

        $("#redraw").removeClass("disabled");
    }

    flextree_change(instance) {
        if ($(instance.currentTarget).prop('checked') == true){
            for (var genome of this.genomes) {
                if ($('#flex' + genome + 'layer').prop('checked') == false){
                    $('#' + genome + 'layer')[0].value = 50;
                    $('#flex' + genome + 'layer').prop('checked', true);
                }
                $('#flex' + genome + 'layer').prop('disabled', true);
            }
        } else {
            for (var genome of this.genomes) {
                $('#flex' + genome + 'layer').prop('disabled', false);
            }
        }
    }
    
    _resolve_node_ids(e, gene_cluster_id='') {
        const id = e.id;
        const element = document.getElementById(id);
        let gene_cluster_context;
        if (element.getAttribute('class') == 'group') {
            gene_cluster_context = id;
            if (!gene_cluster_id) gene_cluster_id = this.group_dict[id][0];
        } else {
            gene_cluster_id = id;
            gene_cluster_context = null;
        }
        return { gene_cluster_id, gene_cluster_context };
    }

    async nodeinfo(e, gene_cluster_id='') {
        const { gene_cluster_id: gcid, gene_cluster_context } = this._resolve_node_ids(e, gene_cluster_id);

        const all_info = await this.get_gene_cluster_display_tables(gcid, gene_cluster_context, 1, false);
        const title = `Synteny gene cluster: ${gcid}`;
        showPangraphFunctionsSummaryTableDialog(title, all_info);

        setTimeout(() => {
            document.querySelectorAll('.group_choice').forEach(el => {
                el.addEventListener('click', () => { this.nodeinfo(e, el.getAttribute("name_id")); });
            });
            setupItemTableFiltering(this._last_gene_clusters);
        }, 100);
    }

    async nodeinfo_with_functions(e, gene_cluster_id='') {
        const { gene_cluster_id: gcid, gene_cluster_context } = this._resolve_node_ids(e, gene_cluster_id);

        if (this.server_offline) {
            toastr.error('The server is no longer accessible.', 'Request failed');
            return;
        }
        waitingDialog.show('Fetching functions and metabolism data...', { dialogSize: 'sm' });
        let all_info;
        try {
            all_info = await this.get_gene_cluster_display_tables(gcid, gene_cluster_context, 1, true);
        } catch(err) {
            waitingDialog.hide();
            toastr.error('Could not load data.', "Request failed");
            return;
        }
        waitingDialog.hide();

        const title = `Synteny gene cluster: ${gcid}`;
        showPangraphFunctionsSummaryTableDialog(title, all_info);

        setTimeout(() => {
            document.querySelectorAll('.group_choice').forEach(el => {
                el.addEventListener('click', () => { this.nodeinfo_with_functions(e, el.getAttribute("name_id")); });
            });
            setupItemTableFiltering(this._last_gene_clusters);
        }, 100);
    }

    show_node_context_menu(event, node_el) {
        event.preventDefault();
        document.querySelectorAll('.context-menu').forEach(m => m.remove());

        const menu = document.createElement('ul');
        menu.setAttribute('class', 'dropdown-menu context-menu');
        menu.setAttribute('role', 'menu');
        menu.style.display = 'block';
        menu.style.visibility = 'hidden';

        const items = [
            { title: 'Show functions for this SynGC', action: () => this.nodeinfo_with_functions(node_el) },
            { title: 'Add SynGC as a new bin', action: () => {
                this.add_bin();
                this.marknode(node_el, this.current_bin_id);
            }},
            { title: 'Append SynGC into the active bin', action: () => {
                this.marknode(node_el, this.current_bin_id);
            }},
        ];

        for (const item of items) {
            const li = document.createElement('li');
            const a = document.createElement('a');
            a.setAttribute('class', 'dropdown-item');
            a.setAttribute('href', '#');
            a.textContent = item.title;
            a.addEventListener('click', (e) => { e.preventDefault(); item.action(); menu.remove(); });
            li.appendChild(a);
            menu.appendChild(li);
        }

        const hr = document.createElement('li');
        hr.innerHTML = '<hr class="dropdown-divider">';
        menu.appendChild(hr);

        const copyLi = document.createElement('li');
        const copyA = document.createElement('a');
        copyA.setAttribute('class', 'dropdown-item');
        copyA.setAttribute('href', '#');
        copyA.textContent = 'Copy SynGC name into clipboard';
        copyA.addEventListener('click', (e) => {
            e.preventDefault();
            const { gene_cluster_id: gcid } = this._resolve_node_ids(node_el);
            navigator.clipboard.writeText(gcid).then(() => {
                toastr.info(`'${gcid}' copied to clipboard.`);
            });
            menu.remove();
        });
        copyLi.appendChild(copyA);
        menu.appendChild(copyLi);

        document.body.appendChild(menu);

        const maxLeft = window.innerWidth - menu.clientWidth;
        const maxTop = window.innerHeight - menu.clientHeight;
        menu.style.left = Math.min(maxLeft, event.clientX) + 'px';
        menu.style.top = Math.min(maxTop, event.clientY) + 'px';
        menu.style.visibility = '';

        document.addEventListener('click', () => menu.remove(), { once: true });
    }

    async get_gene_cluster_display_tables(gene_cluster_id, gene_cluster_context, add_align, show_functions=false) {
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
        
        var gene_cluster_context_table = this.get_gene_cluster_context_table(gene_cluster_id, gene_cluster_context, add_align);
        
        ///////////////////////////////////////////////////////////////////////////////////////////
        // BUILD BASIC INFORMATION TABLE
        ///////////////////////////////////////////////////////////////////////////////////////////
        
        var basic_info_table = this.get_gene_cluster_basics_table(gene_cluster_id, gene_cluster_context, add_align);
        
        ///////////////////////////////////////////////////////////////////////////////////////////
        // BUILD FUNCTIONS TABLE
        ///////////////////////////////////////////////////////////////////////////////////////////
        
        var basic_layer_table = this.get_layer_data(gene_cluster_id, add_align);

        var functions_table = show_functions ? await this.get_gene_cluster_functions_table(gene_cluster_id, add_align) : '';

        var regions_table = await this.get_region_data(gene_cluster_id, add_align);
        
        ///////////////////////////////////////////////////////////////////////////////////////////
        // RETRIEVE AND BUILD SEQUENCE ALIGNMENTS
        ///////////////////////////////////////////////////////////////////////////////////////////
        
        ///////////////////////////////////////////////////////////////////////////////////////////
        // MERGE ALL AND RETURN
        // Order depends on whether this is a left-click (no functions) or right-click (with
        // functions). Right-click order: Basics → Metabolism+Functions → Alignments → the rest.
        ///////////////////////////////////////////////////////////////////////////////////////////

        if (add_align == 1) {
            var gene_cluster_sequence_alignments_table = await this.appendalignment(gene_cluster_id)

            if (show_functions) {
                return gene_cluster_context_table + basic_info_table + functions_table + gene_cluster_sequence_alignments_table + basic_layer_table + regions_table;
            } else {
                return gene_cluster_context_table + basic_info_table + gene_cluster_sequence_alignments_table + basic_layer_table + regions_table;
            }
        } else {
            if (show_functions) {
                return gene_cluster_context_table + basic_info_table + functions_table + basic_layer_table + regions_table;
            } else {
                return gene_cluster_context_table + basic_info_table + basic_layer_table + regions_table;
            }
        }
    }
    
    get_layer_data(gene_cluster_id, add_align) {
    
        var layers = Object.keys(this.data['meta']['layers'])
        var basic_info = {}
        
        for (var layer_name of this.layers) {
            // if ($('#flex' + layer_name).prop('checked') == true){
            basic_info[layer_name] = this.data['nodes'][gene_cluster_id]['layer'][layer_name]
            // }
        }
        
        if (add_align == 1) {
            var basic_layer_table = `<p class="bin-modal-header" style="background: #c8e6c978;">Layers</p>`;
            basic_layer_table += `<table class="table table-striped table-bordered sortable" id="node_layers_table">`;
        } else {
            var basic_layer_table = ''
            basic_layer_table += `<table class="table table-striped table-bordered sortable">`;
        }

        basic_layer_table += `<thead><tr>`;
        basic_layer_table += `<th scope="col" style="width: 25%;">Item</th>`;
        basic_layer_table += `<th scope="col">Value</th>`;
        
        basic_layer_table += `</tr></thead>`;
        basic_layer_table += `<tbody>`;
        for (const [key, value] of Object.entries(basic_info)) {
            basic_layer_table += `<tr>`
            basic_layer_table += `<td>` + key + `</td>`;
            basic_layer_table += `<td>` + value + `</td>`;
            basic_layer_table += `</tr>`;
        }
        basic_layer_table += `</tbody></table>`;
        
        return basic_layer_table;

    }
    
    get_gene_cluster_basics_table(gene_cluster_id, gene_cluster_context, add_align) {
        // first, learn a few basics about the gene cluster to be displayed
        var x_pos = this.data['nodes'][gene_cluster_id]['position'][0]
        // var position_in_graph = x_pos + " / " + (data["infos"]["meta"]["global_x_offset"] - 1);
        var position_in_graph = x_pos;
        // var num_contributing_genomes = Object.keys(data['elements']['nodes'][gene_cluster_id]['genome']).length + " / " + (data['infos']['num_genomes']);
        var num_contributing_genomes = Object.keys(this.data['nodes'][gene_cluster_id]['gene_calls']).length;
        var gene_cluster_name = this.data['nodes'][gene_cluster_id]['gene_cluster']
        
        if (gene_cluster_context == null){
            var context = gene_cluster_id
        } else {
            var context = gene_cluster_context
        }
        
        var basic_info = {'Synteny Gene Cluster': gene_cluster_id, 'Gene Cluster': gene_cluster_name, 'Contributing Genomes': num_contributing_genomes, 'Position in Graph': position_in_graph}
        // build the basic information table
        if (add_align == 1) {
            var basic_info_table = `<p class="bin-modal-header" style="background: #d0e4f578;">Basics</p>`;
            basic_info_table += `<table class="table table-striped table-bordered sortable" gc_context="` + context + `" gc_pos="` + x_pos + `" gc_id="` + gene_cluster_id + `" id="node_basics_table">`;
        } else {
            var basic_info_table = ''
            basic_info_table += `<table class="table table-striped table-bordered sortable">`;
        }
        
        // basic_info_table += `<thead class="thead-light"><tr>`;
        basic_info_table += `<thead><tr>`;
        for (var [key, value] of Object.entries(basic_info)) {
            basic_info_table += `<th scope="row">` + key + `</th>`;
        }
        basic_info_table += `</tr></thead>`;
        basic_info_table += `<tbody>`;
        basic_info_table += `<tr>`;
        for (var [key, value] of Object.entries(basic_info)) {
            basic_info_table += `<td>` + value + `</td>`;
        }
        basic_info_table += `</tr></tbody></table>`;
        
        return basic_info_table;
    }

    async get_region_data(gene_cluster_id, add_align) {

      var d = await this.get_gene_cluster_region_data([gene_cluster_id]);
      var region_info = d['data'][gene_cluster_id]
        
      if (add_align == 1) {
            var region_table = `<p class="bin-modal-header" style="background: #b2dfdb78;">Details from graph</p>`;
            region_table += `<table class="table table-striped table-bordered sortable" id="node_regions_table">`;
        } else {
            var region_table = ''
            region_table += `<table class="table table-striped table-bordered sortable">`;
        }

        region_table += `<thead><tr>`;
        region_table += `<th scope="col" style="width: 25%;">Metric</th>`;
        region_table += `<th scope="col">Value</th>`;
        
        region_table += `</tr></thead>`;
        region_table += `<tbody>`;
        for (const [key, value] of Object.entries(region_info)) {
            region_table += `<tr>`
            region_table += `<td>` + key + `</td>`;
            region_table += `<td>` + value + `</td>`;
            region_table += `</tr>`;
        }
        region_table += `</tbody></table>`;
        
        return region_table;

    }
    
    get_pangraph_gc_config() {
        return {
            itemLabel: 'synteny gene clusters',
            itemIdLabel: 'SynGC',
            metabolismDescription: 'Metabolic modules these synteny gene clusters are involved in. Completeness scores show what fraction of each module is represented by the gene clusters in view.',
            functionsDescription: 'Functional annotations for each synteny gene cluster by each annotation source available in the database.',
            dialogFunction: 'showPangraphFunctionsSummaryTableDialog',
            getAccessionString: (d, source) => {
                const result = d[source]?.accession;
                return (!result || result === '-') ? 'N/A' : result;
            },
            getFunctionString: (d, source) => {
                const result = d[source]?.function;
                return (!result || result === '-') ? 'N/A' : result;
            }
        };
    }

    async fetch_functions_and_metabolism(sgc_ids) {
        const response = await $.ajax({
            url: '/pangraph/get_pangraph_synteny_gc_functions_and_metabolism',
            type: 'POST',
            data: JSON.stringify({ synteny_gene_clusters: sgc_ids }),
            contentType: 'application/json',
            dataType: 'json',
            timeout: 10000,
        });
        return response;
    }

    async get_gene_cluster_functions_table(gene_cluster_id, add_align) {
        let response;
        try {
            response = await this.fetch_functions_and_metabolism([gene_cluster_id]);
        } catch(err) {
            return `<p class="text-muted" style="margin:10px 0">Could not reach the functions endpoint.</p>`;
        }

        if (!response || response.status !== 0) {
            return `<p class="text-muted" style="margin:10px 0">${(response && response.message) || 'Could not load functional annotations.'}</p>`;
        }

        this._last_gene_clusters = response['gene_clusters'] || {};
        return buildFunctionsContent(response, this.get_pangraph_gc_config());
    }

    async show_bin_functions(bin_id) {
        const items = this.bin_dict[bin_id] || [];
        if (!items.length) {
            toastr.warning('There are no synteny gene clusters in this bin yet.', "Nothing to show");
            return;
        }

        // Expand group nodes to their constituent synteny gene clusters
        const sgc_ids = [];
        for (const id of items) {
            const members = id.startsWith('GCG_') ? (this.group_dict[id] || []) : [id];
            for (const gc of members) {
                if (!sgc_ids.includes(gc)) sgc_ids.push(gc);
            }
        }

        const raw_name = document.getElementById(bin_id + '_text')?.value || bin_id;
        const bin_name = raw_name.replace(/_/g, ' ').replace(/^./, c => c.toUpperCase());

        if (this.server_offline) {
            toastr.error('The server is no longer accessible.', 'Request failed');
            return;
        }
        waitingDialog.show('Fetching functions and metabolism data...', { dialogSize: 'sm' });

        let response;
        try {
            response = await this.fetch_functions_and_metabolism(sgc_ids);
        } catch(err) {
            waitingDialog.hide();
            toastr.error('Could not reach the functions endpoint.', "Request failed");
            return;
        }

        waitingDialog.hide();

        if (!response || response.status !== 0) {
            toastr.error((response && response.message) || 'Could not load functional annotations.', "Server error");
            return;
        }

        this._last_gene_clusters = response['gene_clusters'] || {};
        const title = `A summary of functions for ${sgc_ids.length} synteny gene clusters in "${bin_name}"`;
        showPangraphFunctionsSummaryTableDialog(title, buildFunctionsContent(response, this.get_pangraph_gc_config()));
        setTimeout(() => setupItemTableFiltering(this._last_gene_clusters), 100);
    }

    get_region_svg_node_ids(rid) {
        const rinfo = this.data['regions'][rid];
        if (!rinfo) return [];
        const x_min = rinfo['x_min'];
        const x_max = rinfo['x_max'];

        // Build reverse mapping: node_id -> group_id
        const node_to_group = {};
        for (const [group_id, node_ids] of Object.entries(this.group_dict)) {
            for (const node_id of node_ids) {
                node_to_group[node_id] = group_id;
            }
        }

        const result = new Set();
        for (const [node_id, node_data] of Object.entries(this.data['nodes'])) {
            const x = node_data['position'][0];
            if (x >= x_min && x <= x_max) {
                result.add(node_id in node_to_group ? node_to_group[node_id] : node_id);
            }
        }
        return Array.from(result);
    }

    show_region_context_menu(event, rid) {
        document.querySelectorAll('.context-menu').forEach(m => m.remove());

        const menu = document.createElement('ul');
        menu.setAttribute('class', 'dropdown-menu context-menu');
        menu.setAttribute('role', 'menu');
        menu.style.display = 'block';
        menu.style.visibility = 'hidden';

        const items = [
            { title: 'Show summary of functions in region', action: () => this.show_region_functions(rid) },
            { title: 'Add SynGCs in region as a new bin',   action: () => this.add_region_as_new_bin(rid) },
            { title: 'Append SynGCs in region into the active bin', action: () => this.append_region_to_active_bin(rid) },
        ];

        for (const item of items) {
            const li = document.createElement('li');
            const a = document.createElement('a');
            a.setAttribute('class', 'dropdown-item');
            a.setAttribute('href', '#');
            a.textContent = item.title;
            a.addEventListener('click', (e) => { e.preventDefault(); item.action(); menu.remove(); });
            li.appendChild(a);
            menu.appendChild(li);
        }

        document.body.appendChild(menu);

        const maxLeft = window.innerWidth - menu.clientWidth;
        const maxTop = window.innerHeight - menu.clientHeight;
        menu.style.left = Math.min(maxLeft, event.clientX) + 'px';
        menu.style.top = Math.min(maxTop, event.clientY) + 'px';
        menu.style.visibility = '';

        document.addEventListener('click', () => menu.remove(), { once: true });
    }

    async show_region_functions(rid) {
        const svg_ids = this.get_region_svg_node_ids(rid);
        if (!svg_ids.length) {
            toastr.warning('There are no synteny gene clusters in this region.', "Nothing to show");
            return;
        }

        const sgc_ids = [];
        for (const id of svg_ids) {
            const members = id.startsWith('GCG_') ? (this.group_dict[id] || []) : [id];
            for (const gc of members) {
                if (!sgc_ids.includes(gc)) sgc_ids.push(gc);
            }
        }

        if (this.server_offline) {
            toastr.error('The server is no longer accessible.', 'Request failed');
            return;
        }
        waitingDialog.show('Fetching functions and metabolism data...', { dialogSize: 'sm' });

        let response;
        try {
            response = await this.fetch_functions_and_metabolism(sgc_ids);
        } catch(err) {
            waitingDialog.hide();
            toastr.error('Could not reach the functions endpoint.', "Request failed");
            return;
        }

        waitingDialog.hide();

        if (!response || response.status !== 0) {
            toastr.error((response && response.message) || 'Could not load functional annotations.', "Server error");
            return;
        }

        this._last_gene_clusters = response['gene_clusters'] || {};
        const title = `A summary of functions for ${sgc_ids.length} synteny gene clusters in region #${rid}`;
        showPangraphFunctionsSummaryTableDialog(title, buildFunctionsContent(response, this.get_pangraph_gc_config()));
        setTimeout(() => setupItemTableFiltering(this._last_gene_clusters), 100);
    }

    add_region_as_new_bin(rid) {
        const svg_ids = this.get_region_svg_node_ids(rid);
        if (!svg_ids.length) {
            toastr.warning('There are no synteny gene clusters in this region.', "Nothing to show");
            return;
        }
        this.add_bin();
        const bin_id = this.current_bin_id;
        for (const svg_id of svg_ids) {
            const el = document.getElementById(svg_id);
            if (el) this.marknode(el, bin_id);
        }
        toastr.success(`Added ${svg_ids.length} items to new bin "${bin_id}".`, "Bin created");
    }

    append_region_to_active_bin(rid) {
        const svg_ids = this.get_region_svg_node_ids(rid);
        if (!svg_ids.length) {
            toastr.warning('There are no synteny gene clusters in this region.', "Nothing to show");
            return;
        }
        const bin_id = this.current_bin_id;
        let added = 0;
        for (const svg_id of svg_ids) {
            if (!this.bin_dict[bin_id].includes(svg_id)) {
                const el = document.getElementById(svg_id);
                if (el) { this.marknode(el, bin_id); added++; }
            }
        }
        toastr.success(`Appended ${added} new item(s) to "${bin_id}".`, "Bin updated");
    }


    get_gene_cluster_context_table(gene_cluster_id_current, gene_cluster_context, add_align) {
        
        if (gene_cluster_context == null){
            return '';
        } else {
          var group_context = []
          for (var item in this.group_dict[gene_cluster_context]){
            group_context.push(this.group_dict[gene_cluster_context][item])
          }
        }
    
        // console.log(gene_cluster_id_current, gene_cluster_context, group_context)
        if (add_align == 1) {
          var gene_cluster_context_table = `<p class="bin-modal-header" style="background: #ddd8f578;">Gene cluster context</p>`;
        }else {
          var gene_cluster_context_table = ''
        }
        gene_cluster_context_table += `<div class="gene_cluster_context_items">`;
        for(var gene_cluster_id of group_context) {
            var gene_cluster_name = this.data['nodes'][gene_cluster_id]['gene_cluster']
    
            if (gene_cluster_id == gene_cluster_id_current){
                gene_cluster_context_table += `<span class="gene_cluster_id gene_cluster_id_current">` + gene_cluster_name + `</span>`;
            } else {
                gene_cluster_context_table += `<span class="gene_cluster_id"><a class="btn border-0 m-0 p-0 align-baseline group_choice" context="` + add_align + `" group="` + gene_cluster_context + `" name_id="` + gene_cluster_id + `">` + gene_cluster_name + `</a></span>`;
            }
        }
        gene_cluster_context_table += `</div>`;
    
        return gene_cluster_context_table;
    
    }

    async appendalignment(gene_cluster_id) {
        var d = await this.fetchalignment([gene_cluster_id]);
        var alignment = d['data'][gene_cluster_id];

        // Collect sequences in display order to build conservation-based per-column colors
        var sequence_entries = [];
        for (var genome of this.genomes) {
            if (genome in alignment) {
                for (var [gene_call, sequence] of Object.entries(alignment[genome])) {
                    sequence_entries.push({genome, gene_call, sequence});
                }
            }
        }

        // Build sequences_array[col][seq_idx] as required by determineColor; cache for recolor_alignment()
        var max_len = sequence_entries.reduce((m, e) => Math.max(m, e.sequence.length), 0);
        var sequences_array = [];
        for (var col = 0; col < max_len; col++) {
            sequences_array.push(sequence_entries.map(e => col < e.sequence.length ? e.sequence[col] : '-'));
        }
        this._alignment_sequence_entries = sequence_entries;
        this._alignment_sequences_array = sequences_array;

        var coded_positions = window.determineColor(sequences_array);

        // Build checkbox controls (IDs match AA letters so checked() in color-coding.js finds them)
        const aa_letters = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'];
        var controls_html = `<div style="margin: 6px 0 10px 0; font-size: 12px;">`;
        controls_html += `<span style="margin-right: 6px; font-weight: bold;">Color by conservation:</span>`;
        for (var aa of aa_letters) {
            controls_html += `<label style="margin-right: 5px; cursor: pointer;"><input type="checkbox" id="${aa}" class="pangraph-aa-checkbox" checked style="margin-right: 2px;">${aa}</label>`;
        }
        controls_html += ` <button type="button" id="pangraph-aa-check-all" class="btn btn-sm btn-outline-secondary" style="margin-left: 8px; font-size: 11px; padding: 0 6px;">All</button>`;
        controls_html += ` <button type="button" id="pangraph-aa-uncheck-all" class="btn btn-sm btn-outline-secondary" style="font-size: 11px; padding: 0 6px;">None</button>`;
        controls_html += `</div>`;

        var alignments_table = `<p class="bin-modal-header" style="background: #f8bbd078;">Sequence alignments</p>`;
        alignments_table += controls_html;
        alignments_table += `<div class="scroll-wrapper"><table class="table sortable" gc_id="${gene_cluster_id}" id="node_sequence_alignments_table">`;
        alignments_table += `<thead class="gc-table-header"><tr>`;
        alignments_table += `<th class="position-sticky" style="left:0px; z-index:2;" scope="col">Genome</th>`;
        alignments_table += `<th scope="col">Gene Call</th>`;
        alignments_table += `<th scope="col"><span id="th-sequence">Sequence</span></th>`;
        alignments_table += `</tr></thead><tbody>`;

        for (var seq_idx = 0; seq_idx < sequence_entries.length; seq_idx++) {
            var entry = sequence_entries[seq_idx];
            var colored_sequence = '';
            for (var col = 0; col < entry.sequence.length; col++) {
                var aa = entry.sequence[col];
                var color_name = 'black';
                if (coded_positions[col] && coded_positions[col][seq_idx]) {
                    color_name = coded_positions[col][seq_idx][aa] || 'black';
                }
                colored_sequence += `<span style="color: ${color_name};">${aa}</span>`;
            }
            alignments_table += `<tr>`;
            alignments_table += `<td id="td-genome-cell">${entry.genome}</td>`;
            alignments_table += `<td id="td-value-cell">${entry.gene_call}</td>`;
            alignments_table += `<td id="gc-alignment-font"><div class="scrollable-content">${colored_sequence}</div></td>`;
            alignments_table += `</tr>`;
        }

        alignments_table += `</tbody></table></div>`;
        return alignments_table;
    }

    recolor_alignment() {
        if (!this._alignment_sequence_entries || !this._alignment_sequences_array) return;
        var coded_positions = window.determineColor(this._alignment_sequences_array);
        var rows = document.querySelectorAll('#node_sequence_alignments_table tbody tr');
        for (var seq_idx = 0; seq_idx < rows.length; seq_idx++) {
            var cell = rows[seq_idx].querySelector('#gc-alignment-font .scrollable-content');
            if (!cell) continue;
            var entry = this._alignment_sequence_entries[seq_idx];
            if (!entry) continue;
            var colored_sequence = '';
            for (var col = 0; col < entry.sequence.length; col++) {
                var aa = entry.sequence[col];
                var color_name = 'black';
                if (coded_positions[col] && coded_positions[col][seq_idx]) {
                    color_name = coded_positions[col][seq_idx][aa] || 'black';
                }
                colored_sequence += `<span style="color: ${color_name};">${aa}</span>`;
            }
            cell.innerHTML = colored_sequence;
        }
    }

    get_gene_cluster_region_data(gene_cluster_names) {

        var func = {};
        func['synteny_gene_clusters'] = gene_cluster_names

        var d = $.ajax({
            url: "/pangraph/get_pangraph_synteny_gene_cluster_region",
            type: "POST",
            data: JSON.stringify(func),
            contentType: "application/json",
            dataType: "json",
            timeout: 10000,
            error: function(){
                toastr.error('Failed to fetch region data.', 'Error');
            },
            success: function(){}
        })

        return d
    }

    get_gene_cluster_consensus_functions(gene_cluster_names) {

        var func = {};
        func['synteny_gene_clusters'] = gene_cluster_names

        var d = $.ajax({
            url: "/pangraph/get_pangraph_synteny_gene_cluster_function",
            type: "POST",
            data: JSON.stringify(func),
            contentType: "application/json",
            dataType: "json",
            timeout: 10000,
            error: function(){
                toastr.error('Failed to fetch function data.', 'Error');
            },
            success: function(){}
        })

        return d
        }

    fetchalignment(gene_cluster_names) {

        var func = {};
        func['synteny_gene_clusters'] = gene_cluster_names

        var d = $.ajax({
            url: "/pangraph/get_pangraph_synteny_gene_cluster_alignment",
            type: "POST",
            data: JSON.stringify(func),
            contentType: "application/json",
            dataType: "json",
            timeout: 10000,
            error: function(){
                toastr.error('Failed to fetch alignment data.', 'Error');
            },
            success: function(){}
        })
        
        return d
    }

    info_download () {
    
        var csv_data = [];
        var basics = $('#node_basics_table')
        var title = basics[0].getAttribute("gc_id")
        var layers = $('#node_layers_table')
        var regions = $('#node_regions_table')
        var functions = $('#node_functions_table')
        
        var basics_rows = basics[0].getElementsByTagName('tr');
        var layers_rows = layers[0].getElementsByTagName('tr');
        var regions_rows = regions[0].getElementsByTagName('tr');
        
        var function_rows = functions[0].getElementsByTagName('tr');
        
        for (var i = 0; i < function_rows.length; i++) {
            
            if (i >= basics_rows.length) {
                var basics_cols = []
                basics_cols = Array.prototype.concat.apply(basics_cols, basics_rows[1].querySelectorAll('td,th'));
                basics_cols = Array.prototype.concat.apply(basics_cols, layers_rows[1].querySelectorAll('td,th'));
                basics_cols = Array.prototype.concat.apply(basics_cols, regions_rows[1].querySelectorAll('td,th'));
                var function_cols = function_rows[i].querySelectorAll('td,th');
            } else { 
                var basics_cols = []
                basics_cols = Array.prototype.concat.apply(basics_cols, basics_rows[i].querySelectorAll('td,th'));
                basics_cols = Array.prototype.concat.apply(basics_cols, layers_rows[i].querySelectorAll('td,th'));
                basics_cols = Array.prototype.concat.apply(basics_cols, regions_rows[i].querySelectorAll('td,th'));
                var function_cols = function_rows[i].querySelectorAll('td,th');
            }
            
            let csvrow = [];
            for (var j = 0; j < basics_cols.length; j++) {
                var info = basics_cols[j].innerHTML
                csvrow.push(info);
            }
            
            for (var k = 0; k < function_cols.length; k++) {
                var info = function_cols[k].innerHTML
                csvrow.push(info);
            }
            csv_data.push(csvrow.join(","));
        }
        csv_data = csv_data.join('\n');
        
        var blob = new Blob([csv_data]);
        var title = this.data['meta']['project_name']
        this.download_blob(blob, title + ".csv");
    }

    alignment_download () {

        var csv_data = '';
        var alignment = $('#node_sequence_alignments_table')
        var basics = $('#node_basics_table')
        var title = alignment[0].getAttribute("gc_id")
        var xpos = basics[0].getAttribute("gc_pos")
        
        var alignment_rows = alignment[0].getElementsByTagName('tr');
        
        for (var i = 1; i < alignment_rows.length; i++) {
        
            var alignment_cols = alignment_rows[i].querySelectorAll('td,th');
            
            csv_data += ">" + title + "|Genome:" + alignment_cols[0].innerHTML +"|Genecall:" + alignment_cols[1].innerHTML + "|Position:" + xpos + '\n';
            var genome = ''
            
            var alignment_nucs = alignment_cols[2].getElementsByTagName('span');
            for (let k = 0; k < alignment_nucs.length; k++) {
                genome += alignment_nucs[k].innerHTML
            }
            
            csv_data += genome.match(/.{1,60}/g).join("\r\n") + "\n"
        }
        var blob = new Blob([csv_data]);
        this.download_blob(blob, title + ".fa");
    }

    save_bin() {
        const exportFn = () => {
            const data = {};
            const colors = {};
            for (const [bin_id, node_ids] of Object.entries(this.bin_dict)) {
                if (!node_ids.length) continue;
                const bin_name = $('#' + bin_id + '_text').val() || bin_id;
                data[bin_name] = node_ids;
                colors[bin_name] = $('#' + bin_id + '_color').attr('color') || '#000000';
            }
            return { data, colors };
        };

        new StoreCollectionDialog(exportFn).Show();
    }

    load_bin() {
        const importFn = (collection_data) => {
            const bin_data = collection_data['data'];
            const colors   = collection_data['colors'];

            this._suppress_bin_ring_draw = true;

            // Uncolor and clear all existing bins
            for (const [bin_id, node_ids] of Object.entries(this.bin_dict)) {
                for (const node of [...node_ids]) {
                    const el = document.getElementById(node);
                    if (el) this.marknode(el, bin_id);
                }
                $('#' + bin_id + '_grid').remove();
            }

            this.bin_dict = {};
            this.bin_group_dict = {};
            this.current_bin_number = 0;

            // Recreate bins from the loaded collection
            for (const [bin_name, node_ids] of Object.entries(bin_data)) {
                this.add_bin();
                const bin_id = this.current_bin_id;
                const color = colors[bin_name] || '#000000';

                $('#' + bin_id + '_text').val(bin_name);
                $('#' + bin_id + '_color').css('background-color', color).attr('color', color);

                for (const node_id of node_ids) {
                    const el = document.getElementById(node_id);
                    if (el) this.marknode(el, bin_id);
                }
            }

            // If nothing was loaded, add an empty default bin
            if (!Object.keys(this.bin_dict).length) {
                this.add_bin();
            }

            this._suppress_bin_ring_draw = false;
            this.draw_bin_rings();
        };

        new LoadCollectionDialog(importFn).Show();
    }

    load_state () {
        if ($('#loadstatename')[0].value != "") {
            this.state = $('#loadstatename')[0].value

            var state = {}
            state['state_name'] = this.state
    
            $.ajax({
                url: "/pangraph/load_pangraph_state",
                type: "POST",
                cache: false,
                contentType: "application/json",
                dataType: "json",
                data: JSON.stringify(state),
                success: (data) => {
                    this.data = data['data'];
                    this.initialize_variables();
                    this.set_UI_settings();
                    this.main_draw();
                    toastr.success(`State "${this.state}" has been loaded.`, 'State loaded');
                },
                error: (err) => {
                    toastr.error('Failed to load state.', 'Error');
                }
            })
            
            $('#loadstatemodal').modal('hide');
        }
    }
    
    save_state () {
        var new_state = {}
        
        for (var [setting, value] of Object.entries(this.data['states'])) {
            const el = $('#' + setting);
            if (el.hasClass('pangraph-colorpicker')) {
                new_state[setting] = el.attr('color');
            } else if (typeof value === 'number') {
                new_state[setting] = Number(el[0].value);
            } else if (value == true || value == false) {
                new_state[setting] = el.prop('checked');
            } else {
                new_state[setting] = el[0].value;
            }
        }

        var genome_order = [...document.getElementById("genomecolors").children]
            .filter(element => element.classList.contains("col-12"))
            .map(element => element.id.replace('_row', ''));

        genome_order.forEach(key => {
            if (key in new_state) {
              const value = new_state[key];
              delete new_state[key];
              new_state[key] = value;
            }
        });
         
        var result = {}
        result['state_name'] = $('#savestatename')[0].value
        result['state_values'] = new_state

        $.ajax({
            url: "/pangraph/save_pangraph_state",
            type: "POST",
            data: JSON.stringify(result),
            contentType: "application/json",
            dataType: "json",
            error: function(){
                toastr.error('Failed to save state.', 'Error');
            },
            success: function(){
                toastr.success(`State "${result['state_name']}" has been saved.`, 'State saved');
            }
        })
        
        $('#savestatemodal').modal('hide');
    }
    
    state_modal (save_load) {

        $.ajax({
            type: 'GET',
            cache: false,
            url: '/pangraph/get_pangraph_states',
            success: function(data) {

                $('#' + save_load + 'statelisttab').empty()
                $('#' + save_load + 'statename')[0].value = ""
                
                for (var state_name of data['data']) {
                    $('#' + save_load + 'statelisttab').append(
                        $('<a class="list-group-item list-group-item-action" id="list-profile-list" data-toggle="list" href="#list-profile" role="tab" aria-controls="profile">').append(
                            state_name
                        )
                    )
                }
    
                $('#' + save_load + 'statelisttab a').on('click', function (e) {
                    e.preventDefault()
                    $(this).tab('show')
                    $('#' + save_load + 'statename')[0].value = $(this)[0].innerText
                })
    
                $('#' + save_load + 'statemodal').modal('show');
            }
        });
    }
}

$(document).ready(function () {
    const pgui = new PangenomeGraphUserInterface();
    window.pgui = pgui;  // Make accessible to Selenium for SVG export
});
