/**
 * Javascript library to display phylogenetic trees and more
 *
 *  Authors: Özcan Esen <ozcanesen@gmail.com>
 *           A. Murat Eren
 *
 * Copyright 2015-2021, The anvi'o project (http://anvio.org)
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

var Drawer = function(settings) {
    this.timer = new BasicTimer('tree_draw');
    this.settings = settings;

    this.tree_svg_id = 'tree';
    this.fontHeight = 10;
    this.root_length = 0.1;

    this.has_tree = (clusteringData.constructor !== Array);

    this.layerdata_dict = new Array();
    this.layer_fonts = new Array();
};

Drawer.prototype.iterate_layers = function(callback) {
    for (var i = 0; i < this.settings['layer-order'].length; i++) {

        // order 0 belongs to tree itself, so layer indexes start from 1
        var order = i + 1;
        var index = parseInt(this.settings['layer-order'][i]);

        var ret = callback.call(this, {
            order: order,
            index: index,
            is_parent:      (layer_types[index] == 0),
            is_stackbar:    (layer_types[index] == 1),
            is_categorical: (layer_types[index] == 2),
            is_numerical:   (layer_types[index] == 3),

            get_view_attribute: (function(key) {
                return this.settings['views'][current_view][index][key];
            }).bind(this),

            get_visual_attribute: (function(key) {
                return this.settings['layers'][index][key];
            }).bind(this),
        });

        if (ret == -1)
            break;
    }
};

Drawer.prototype.draw = function() {
    this.initialize_screen();

    if (parseFloat(this.settings['angle-min']) <= 0) {
        this.settings['angle-min'] = 0;
    }

    if (parseFloat(this.settings['angle-max']) >= 360) {
        this.settings['angle-max'] = 359.9999999;
    }

    this.initialize_tree();

    if (this.has_tree) {
        this.collapse_nodes();
        this.generate_mock_data_for_collapsed_nodes();
    }

    this.assign_leaf_order();
    this.generate_tooltips();
    this.normalize_values();
    this.calculate_bar_sizes();
    this.calculate_leaf_sizes();

    createBin('svg', 'viewport');
    createBin('viewport', 'tree_bin');
    createBin('tree_bin', 'tree');
    drawLine('tree', {'id': '_origin'}, {'x': 0, 'y': 0}, {'x': 0, 'y': 0}, false);
    createBin('tree_bin', 'guide_lines');

    this.calculate_tree_coordinates();
    this.draw_tree();
    this.calculate_layer_boundaries();

    total_radius = this.layer_boundaries[this.layer_boundaries.length - 1][1];
    beginning_of_layers = this.layer_boundaries[1][0];
    layer_boundaries = this.layer_boundaries;

    this.draw_layer_backgrounds();
    this.draw_guide_lines();
    this.draw_categorical_layers();
    this.draw_numerical_layers();
    this.draw_stack_bar_layers();
    this.draw_layer_names();
    this.draw_samples();
    this.overlay_collapsed_node_layers();
    this.draw_collapsed_nodes();

    bins.MigrateCollection();
    bins.RebuildIntersections();
    bins.RedrawLineColors();

    createBin('tree_bin', 'bin');

    this.initialize_tree_observer();

    // Scale to fit window
    bbox = svg.getBBox();
    zoom_reset();
    this.show_drawing_statistics();

    // pan and mouse zoom
    $('svg').svgPan(
        {
            'viewportId': 'viewport',
            'onlyPanOnMouseUp': (this.total_object_count > 10000 || this.tree.leaves.length > 5000),
        });

    this.bind_tree_events();
    initialize_area_zoom(); // area-zoom.js
    this.update_title_panel();

    ANIMATIONS_ENABLED = false;
};


Drawer.prototype.draw_collapsed_nodes = function() {
    for (var i=0; i < collapsedNodes.length; i++) {
        let collapse_attributes = collapsedNodes[i];
        let left_most = this.tree.label_to_leaves[collapse_attributes['left_most']];
        let right_most = this.tree.label_to_leaves[collapse_attributes['right_most']];
        var p = this.tree.FindLowestCommonAncestor(left_most, right_most);

        this.draw_collapsed_node(p, collapse_attributes);
    }
};


Drawer.prototype.assign_leaf_order = function() {
    if (this.has_tree) {
        var n = new NodeIterator(this.tree.root);
        var q = n.Begin();

        var order_counter = 0;
        while (q != null) {
            if (q.IsLeaf()) {
                q.order = order_counter++;
                this.tree.leaves[q.order] = q;
            }
            q=n.Next();
        }

        this.tree.num_leaves = order_counter;
        leaf_count = order_counter;
    }
    else
    {
        for (var i = 0; i < clusteringData.length; i++) {
            let q = this.tree.NewNode();
            q.order = this.tree.num_leaves++;
            q.label = clusteringData[i];
            this.tree.label_to_leaves[q.label] = q;
            this.tree.leaves[q.order] = q;
        }

        leaf_count = clusteringData.length;
    }
};


Drawer.prototype.generate_mock_data_for_collapsed_nodes = function(node_list) {
    if (!this.has_tree)
        return;

    for (var i=0; i < collapsedNodes.length; i++) {
        let collapse_attributes = collapsedNodes[i];
        let left_most = this.tree.label_to_leaves[collapse_attributes['left_most']];
        let right_most = this.tree.label_to_leaves[collapse_attributes['right_most']];
        var q = this.tree.FindLowestCommonAncestor(left_most, right_most);

        var mock_data = [q.label];
        for (var j = 1; j < parameter_count; j++) {
            if (layerdata[0][j].indexOf(';') > -1) {
                var pcount = layerdata[0][j].split(';').length;
                mock_data.push(Array(pcount).fill('1').join(';'));
            }
            else if (isNumber(layerdata[1][j])) {
                mock_data.push(0);
            }
            else {
                mock_data.push("None");
            }
        }

        this.layerdata_dict[q.label] = mock_data;
    }
};

Drawer.prototype.generate_tooltips = function() {
    empty_tooltip = '<tr><td>Item Name</td><td>n/a</td></tr>';
    empty_tooltip += '<tr><td>parent</td><td>n/a</td></tr>';

    for (var i = 1; i < this.settings['layer-order'].length; i++)
    {
        var pindex = this.settings['layer-order'][i];
        var layer_title = layerdata[0][pindex];
        if (layer_title.indexOf('!') > -1)
            layer_title = layer_title.split('!')[0];
        empty_tooltip += '<tr><td>' + layer_title + '</td><td>n/a</td></tr>';
    }

    $('#tooltip_content').html(empty_tooltip);
    $('#mouse_hover_scroll').css('top', '0');

    for (var index = 1; index < layerdata.length; index++)
    {
        var params = layerdata[index];
        this.layerdata_dict[params[0]] = params.slice(0);

        var title = [];
        for (var i = 0; i < this.settings['layer-order'].length; i++)
        {
            var pindex = this.settings['layer-order'][i];

            if (layer_types[pindex] == 0) // check if parent
            {
                if (layerdata[index][pindex] == '')
                {
                    title.push('<td>parent</td><td>n/a</td>');
                }
                else
                {
                    title.push('<td>parent</td><td>' + layerdata[index][pindex] + '</td>');
                }

            }
            else if (layer_types[pindex] == 1) {
                let stack_names = layerdata[0][pindex].split('!')[1].split(';');
                let stack_items = layerdata[index][pindex].split(';');
                let layer_title = layerdata[0][pindex];
                message = '<td>' + layer_title.split('!')[0] + '</td><td><table>';
                for (let j = stack_names.length - 1; j >= 0; j--) {
                    let bar_name = stack_names[j];
                    message += `<tr><td><div class="colorpicker" style="background-color: ${stack_bar_colors[pindex][bar_name]}"></div>${bar_name}</td><td>${stack_items[j]}</td></tr>`;
                }
                message += '</table></td>';

                title.push(message);
            }
            else
            {
                var layer_title = layerdata[0][pindex];
                title.push('<td>' + layer_title + '</td><td>' + layerdata[index][pindex] + '</td>');
            }
        }

        layerdata_title[params[0]] = title;
    }
};

Drawer.prototype.normalize_values = function() {
    // this object maps layer.index to max value in the layer
    // we will need this when calculating bar sizes
    this.param_max = {};

    for (var id in this.layerdata_dict)
    {
        this.iterate_layers(function(layer) {
            if (layer.is_categorical || layer.is_parent) {
                return;
            }
            else if (layer.is_stackbar) {
                // convert ";" string to array after normalization
                var stack_bar_items = this.layerdata_dict[id][layer.index].split(";");
                var normalization = layer.get_view_attribute('normalization');

                for (var j=0; j < stack_bar_items.length; j++) {
                    if (normalization == 'sqrt') {
                        stack_bar_items[j] = Math.sqrt(parseFloat(stack_bar_items[j]));
                    }
                    else if (normalization == 'log') {
                        stack_bar_items[j] = log10(parseFloat(stack_bar_items[j]) + 1);
                    }
                }

                this.layerdata_dict[id][layer.index] = stack_bar_items.slice(0);
            }
            else { // numerical
                if (isNumber(this.layerdata_dict[id][layer.index])) {
                    var normalization = layer.get_view_attribute('normalization');
                    if (normalization == 'sqrt')
                    {
                        this.layerdata_dict[id][layer.index] = Math.sqrt(parseFloat(this.layerdata_dict[id][layer.index]));
                    }
                    else if (normalization == 'log')
                    {
                        this.layerdata_dict[id][layer.index] = log10(parseFloat(this.layerdata_dict[id][layer.index]) + 1);
                    }

                    if (!this.param_max.hasOwnProperty(layer.index)|| parseFloat(this.layerdata_dict[id][layer.index]) > parseFloat(this.param_max[layer.index]))
                    {
                        this.param_max[layer.index] = parseFloat(this.layerdata_dict[id][layer.index]);
                    }
                }
            }
        });
    }
};

Drawer.prototype.calculate_bar_sizes = function() {
    this.iterate_layers(function(layer) {
        if (layer.is_parent || layer.is_categorical) {
            return;
        }
        else
        {
            var min_max_disabled = layer.get_view_attribute('min')['disabled'];
            var min = parseFloat(layer.get_view_attribute('min')['value']);
            var max = parseFloat(layer.get_view_attribute('max')['value']);
            var min_new = null;
            var max_new = null;

            for (var id in this.layerdata_dict)
            {
                if (layer.is_stackbar)
                {
                    var total = 0;

                    for (var j=0; j < this.layerdata_dict[id][layer.index].length; j++)
                    {
                        if (isNumber(this.layerdata_dict[id][layer.index][j])) {
                            total = total + parseFloat(this.layerdata_dict[id][layer.index][j]);
                        }
                    }

                    var multiplier = (total == 0) ? 0 : (parseFloat(layer.get_visual_attribute('height')) / total);

                    for (var j=0; j < this.layerdata_dict[id][layer.index].length; j++)
                    {
                        if (isNumber(this.layerdata_dict[id][layer.index][j])) {
                            this.layerdata_dict[id][layer.index][j] = this.layerdata_dict[id][layer.index][j] * multiplier;
                        }
                    }
                }
                else // numerical data
                {
                    var bar_size = isNumber(this.layerdata_dict[id][layer.index]) ? parseFloat(this.layerdata_dict[id][layer.index]) : 0;

                    if (!min_max_disabled)
                    {
                        if (bar_size > max) {
                            bar_size = max - min;
                        }
                        else if (bar_size < min) {
                            bar_size = 0;
                        }
                        else {
                            bar_size = bar_size - min;
                        }

                        if (bar_size == 0) {
                            this.layerdata_dict[id][layer.index] = 0;
                        } else {
                            this.layerdata_dict[id][layer.index] = bar_size * parseFloat(layer.get_visual_attribute('height')) / (max - min);
                        }
                    }
                    else
                    {
                        if ((min_new == null) || bar_size < min_new) {
                            min_new = bar_size;
                        }

                        if ((max_new == null) || bar_size > max_new) {
                            max_new = bar_size;
                        }

                        if (bar_size == 0) {
                            this.layerdata_dict[id][layer.index] = 0;
                        } else {
                            this.layerdata_dict[id][layer.index] = bar_size * parseFloat(layer.get_visual_attribute('height')) / this.param_max[layer.index];
                        }

                        var min_max_str = "Min: " + min_new + " - Max: " + max_new;
                        $('#min' + layer.index).attr('title', min_max_str);
                        $('#max' + layer.index).attr('title', min_max_str);
                    }
                }
            }

            if (min_max_disabled)
            {
                $('#min' + layer.index).prop('disabled', false);
                $('#max' + layer.index).val(max_new).prop('disabled', false);
            }
        }
    });
};

Drawer.prototype.initialize_tree = function() {
    this.tree = new Tree();
    if (this.has_tree)
    {
        this.tree.Parse(clusteringData.trim(), this.settings['edge-normalization']);

        if (this.tree.error != 0) {
            toastr.error('Error while parsing tree data.');
            return;
        }
        this.tree.ComputeDepths();

        var n = new NodeIterator(this.tree.root);
        var q = n.Begin();
        while (q != null)
        {
            if (!this.tree.has_edge_lengths) {
                q.edge_length = 1;
            }
            q=n.Next();
        }
    }
};

Drawer.prototype.collapse_nodes = function(node_list) {
    for (var i=0; i < collapsedNodes.length; i++) {
        let collapse_attributes = collapsedNodes[i];
        let left_most = this.tree.label_to_leaves[collapse_attributes['left_most']];
        let right_most = this.tree.label_to_leaves[collapse_attributes['right_most']];
        var cnode = this.tree.FindLowestCommonAncestor(left_most, right_most);

        var max_edge = 0;
        var sum_size = 0;
        var n = new PreorderIterator(cnode);
        var q = n.Begin();
        while (q != null) {
            if (q == cnode) {
                 q = n.Next();
                continue;
            }

            var d = q.edge_length;

            if (d < 0.00001) {
                d = 0.0;
            }

            if (q.ancestor == cnode) {
                q.edge_length = d;
            } else {
                q.edge_length = q.ancestor.edge_length + d;
            }

            max_edge = Math.max(max_edge, (q.collapsed) ? q.edge_length + q.max_child_path : q.edge_length);

            if (q.IsLeaf()) {
                sum_size += q.size;
            }
            q = n.Next();
        }

        cnode.label = collapse_attributes['label'];

        cnode.max_child_path = max_edge;
        cnode.size = Math.max(1, sum_size * collapse_attributes['size']);
        cnode.collapse_order = i;
        cnode.child = null;
        cnode.collapsed = true;
    }
};


Drawer.prototype.overlay_collapsed_node_layers = function() {
    for (var i=0; i < collapsedNodes.length; i++) {
        let collapse_attributes = collapsedNodes[i];
        let left_most = this.tree.label_to_leaves[collapse_attributes['left_most']];
        let right_most = this.tree.label_to_leaves[collapse_attributes['right_most']];
        var p = this.tree.FindLowestCommonAncestor(left_most, right_most);

        if (this.settings['tree-type'] == 'circlephylogram') {
            drawPie('tree_bin',
                'overlay_collapsed_' + p.id,
                p.angle - p.size / 2,
                p.angle + p.size / 2,
                beginning_of_layers - 10,
                total_radius + 10,
                (p.size > Math.PI) ? 1 : 0,
                '#F0F0F0',
                1,
                false);
        } else {
            drawPhylogramRectangle('tree_bin',
                'overlay_collapsed_' + p.id,
                beginning_of_layers - 10,
                p.xy.y,
                p.size,
                total_radius - beginning_of_layers + 20,
                '#F0F0F0',
                1,
                false);
        }
    }
}

Drawer.prototype.bind_tree_events = function() {
    var tree_bin = document.getElementById('tree_bin');
    tree_bin.addEventListener('click', lineClickHandler, false);
    tree_bin.addEventListener('contextmenu', lineContextMenuHandler, false);
    tree_bin.addEventListener('mouseover',lineMouseEnterHandler, false);
    tree_bin.addEventListener('mouseout', lineMouseLeaveHandler, false);
};

Drawer.prototype.calculate_leaf_sizes = function() {
    var total_size = 0;

    for (var i=0; i < this.tree.leaves.length; i++) {
        total_size += this.tree.leaves[i].size;
    }

    var total_width;

    if (this.settings['tree-type'] == 'circlephylogram')
    {
        var total_width = Math.toRadians(parseFloat(this.settings['angle-max']) - parseFloat(this.settings['angle-min']));
    }
    else
    {
        var total_width = this.width;
    }

    for (var i=0; i < this.tree.leaves.length; i++) {
        q = this.tree.leaves[i];
        q.size = q.size * (total_width / total_size);

        if (typeof this.smallest_leaf_size === 'undefined' || q.size < this.smallest_leaf_size) {
            this.smallest_leaf_size = q.size;
        }
    }
};

Drawer.prototype.calculate_tree_coordinates = function() {
    if (this.has_tree) {
        this.max_path_length = 0;
        this.tree.root.path_length = this.tree.root.edge_length;

        if (this.settings['tree-type'] == 'circlephylogram') {
            this.root_length = (this.radius / 2) * this.root_length;
        } else {
            this.root_length = (this.height) * this.root_length;
        }

        var n = new PreorderIterator(this.tree.root);
        var q = n.Begin();
        while (q != null) {
            var d = q.edge_length;

            if (d < 0.00001) {
                d = 0.0;
            }

            if (q != this.tree.root) {
                q.path_length = q.ancestor.path_length + d;
            }

            this.max_path_length = Math.max(this.max_path_length, (q.collapsed) ? q.path_length + q.max_child_path : q.path_length);
            q = n.Next();
        }

        var n = new NodeIterator(this.tree.root);
        var p = n.Begin();
        while (p != null) {
            if (p.IsLeaf())
            {
                this.calculate_leaf_coordinate(p);
            }
            else
            {
                this.calculate_internal_node_coordinate(p);
            }
            p = n.Next();
        }
    }
    else
    {
        for (var i=0; i < this.tree.leaves.length; i++) {
            this.calculate_leaf_coordinate(this.tree.leaves[i]);
        }
    }
};

Drawer.prototype.calculate_leaf_coordinate = function(p) {
    if (this.settings['tree-type'] == 'circlephylogram')
    {
        if (p.order == 0) {
            p.angle = Math.min(Math.toRadians(this.settings['angle-min']), Math.toRadians(this.settings['angle-max'])) + p.size / 2;
        } else {
            var prev_leaf = this.tree.leaves[p.order - 1];
            p.angle = prev_leaf.angle + prev_leaf.size / 2 + p.size / 2;
        }

        if (this.has_tree) {
            if (this.tree.has_edge_lengths) {
                p.radius = this.root_length + (p.path_length / this.max_path_length) * ((this.radius / 2) - this.root_length);
            }
            else
            {
                p.radius = this.root_length + ((this.tree.root.depth - p.depth) / this.tree.root.depth) * ((this.radius / 2) - this.root_length);
            }
        } else {
            p.radius = this.radius / 2;
        }

        var pt = [];
        pt['x'] = p.radius * Math.cos(p.angle);
        pt['y'] = p.radius * Math.sin(p.angle);

        p.xy['x'] = pt['x'];
        p.xy['y'] = pt['y'];
    }
    else // phylogram
    {
        var pt = [];

        if (p.order == 0) {
            pt['y'] = p.size / 2;
        } else {
            var prev_leaf = this.tree.leaves[p.order - 1];
            pt['y'] = prev_leaf.xy['y'] + prev_leaf.size / 2 + p.size / 2;
        }

        if (this.has_tree) {
            if (this.tree.has_edge_lengths) {
                pt['x'] = this.root_length + (p.path_length / this.max_path_length) * (this.height - this.root_length);
            } else {
                pt['x'] = this.root_length + ((this.tree.root.depth - p.depth) / this.tree.root.depth) * (this.height - this.root_length);
            }
        } else {
            pt['x'] = 0;
        }

        p.xy['x'] = pt['x'];
        p.xy['y'] = pt['y'];
    }
};

Drawer.prototype.calculate_internal_node_coordinate = function(p) {
    if (this.settings['tree-type'] == 'circlephylogram')
    {
        var left_angle = p.child.angle;
        var right_angle = p.child.GetRightMostSibling().angle;

        p.angle = left_angle + (right_angle - left_angle) / 2;

        if (this.tree.has_edge_lengths) {
            p.radius = this.root_length + (p.path_length / this.max_path_length) * ((this.radius / 2) - this.root_length);
        }
        else
        {
            p.radius = this.root_length + ((this.tree.root.depth - p.depth) / this.tree.root.depth) * ((this.radius / 2) - this.root_length);
        }

        var pt = [];
        pt['x'] = p.radius * Math.cos(p.angle);
        pt['y'] = p.radius * Math.sin(p.angle);

        p.xy['x'] = pt['x'];
        p.xy['y'] = pt['y'];

        var q = p.child;
        while (q) {
            pt = [];

            pt['x'] = p.radius * Math.cos(q.angle);
            pt['y'] = p.radius * Math.sin(q.angle);

            q.backarc = [];
            q.backarc['x'] = pt['x'];
            q.backarc['y'] = pt['y'];

            q = q.sibling;
        }
    }
    else
    {
        var pt = [];

        if (this.tree.has_edge_lengths) {
            pt['x'] = this.root_length + (p.path_length / this.max_path_length) * (this.height - this.root_length);
        } else {
            pt['x'] = this.root_length + ((this.tree.root.depth - p.depth) / this.tree.root.depth) * (this.height - this.root_length);
        }

        var pl = p.child.xy;
        var pr = p.child.GetRightMostSibling().xy;

        pt['y'] = pl['y'] + (pr['y'] - pl['y']) / 2;
        p.xy['x'] = pt['x'];
        p.xy['y'] = pt['y'];
    }
};

Drawer.prototype.draw_tree = function() {
    if (this.settings['tree-type'] == 'phylogram')
        $('#tree_bin').attr('transform', 'rotate(90)');

    if (this.has_tree) {
        this.draw_root();

        var n = new NodeIterator(this.tree.root);
        var p = n.Begin();
        while (p != null) {
            if (p.IsLeaf())
            {
                this.draw_leaf(p);
            }
            else
            {
                this.draw_internal_node(p);
            }

            p = n.Next();
        }
    } else {
        for (var i=0; i < this.tree.leaves.length; i++) {
            var p = this.tree.leaves[i];
            p.backarc = p.xy;
            this.draw_leaf(p);
        }
    }

    // a check is triggered if the user clicks `support_value_checkbox` AFTER
    // drawing the tree. This one is here for those who first click, and then
    // draw so the warning message is added retrospectively:

    if (this.settings['show-support-values']) {
        if (max_branch_support_value_seen == 0 || max_branch_support_value_seen == null){
            $('#max_branch_support_value_seen_is_zero_warning').show();
            $('#support_value_checkbox').prop("checked", false);
        }
    }
};

Drawer.prototype.draw_leaf = function(p) {
    if (this.settings['tree-type'] == 'circlephylogram') {
            var p0 = p.xy
            var p1 = p.backarc;
            drawLine(this.tree_svg_id, p, p0, p1);
    }
    else
    {
        var p0 = p.xy
        var p1 = [];
        var anc = p.ancestor;
        if (anc) {
            p1['x'] = anc.xy['x'];
            p1['y'] = p0['y'];

            drawLine(this.tree_svg_id, p, p0, p1);
        }
    }
};

Drawer.prototype.draw_internal_node = function(p) {
    let PADDING_STYLE = 'stroke:rgba(0,0,0,0);stroke-width:16;';
    const branch_support_values = [];

    let supportValueData = {
        numberRange : [ this.settings['support-range-low'], this.settings['support-range-high'] ],
        colorRange : [ this.settings['support-color-low'], this.settings['support-color-high']],
        showSymbol : this.settings['support-display-symbol'],
        showNumber : this.settings['support-display-number'],
        invertSymbol : this.settings['support-symbol-invert'],
        minRadius: this.settings['support-min-symbol-size'],
        maxRadius : this.settings['support-symbol-size'],
        symbolColor : this.settings['support-symbol-color'],
        secondSymbolColor : this.settings['second-support-symbol-color'],
        fontColor : this.settings['support-font-color'],
        secondFontColor : this.settings['second-support-font-color'],
        fontSize : this.settings['support-font-size'],
        textRotation : this.settings['support-text-rotation'],
        thresholdValue: this.settings['support-threshold'],
        thresholdOperator: this.settings['support-show-operator'],
        thresholdRange0 : [this.settings['support-bootstrap0-range-low'], this.settings['support-bootstrap0-range-high']],
        thresholdRange1: [this.settings['support-bootstrap1-range-low'], this.settings['support-bootstrap1-range-high']],
        symbolDataSource: this.settings['support-symbol-data']
    }

    if (this.settings['tree-type'] == 'circlephylogram')
    {
        var p0 = [];
        var p1 = [];

        p0['x'] = p.xy['x'];
        p0['y'] = p.xy['y'];

        var anc = p.ancestor;
        // if has ancestor we are going to draw vertical line
        if (anc) {
            p0 = p.xy;
            p1 = p.backarc;

            // we draw two lines because one is transparent and makes clicking easier.
            // you can see same thing happening multiple times in this function.
            // ---
            // also actual line drawing should be placed before padding line drawing funtion.
            // in an ideal world they should not share same element ID but they do to trigger
            // same mouse events. So when bin changes color first object with that ID changes color.
            // padding line should placed second to avoid that.

            drawLine(this.tree_svg_id, p, p0, p1);

            this.settings['show-support-values'] ? drawSupportValue(this.tree_svg_id, p, p0, p1, supportValueData) : null;

            let line = drawLine(this.tree_svg_id, p, p0, p1);
            line.setAttribute('style', 'stroke:rgba(0,0,0,0);stroke-width:16;');
            line.classList.add('clone');

        }

        p0 = p.child.backarc;
        p1 = p.child.GetRightMostSibling().backarc;

        var large_arc_flag = (Math.abs(p.child.GetRightMostSibling().angle - p.child.angle) > Math.PI) ? true : false;

        if (typeof p.branch_support === 'string' && p.branch_support.includes('/')) {
            multiple_support_value_seen = true;
            const [branch_support_value0, branch_support_value1] = p.branch_support.split('/').map(parseFloat);
            
            branch_support_values.push(branch_support_value0, branch_support_value1);

            if (branch_support_values.length > 0) {
                const min_support = Math.min(...branch_support_values);
                const max_support = Math.max(...branch_support_values);
            
                if (min_branch_support_value_seen === null || min_support < min_branch_support_value_seen) {
                    min_branch_support_value_seen = min_support;
                } else {
                    this.min_support = null;
                }
            
                if (max_branch_support_value_seen === null || max_support > max_branch_support_value_seen) {
                    max_branch_support_value_seen = max_support;
                } else {
                    this.max_support = null;
                }
            } 

        } else {
                // If there are no string branch support value like '100/100':
                min_branch_support_value_seen == null ? min_branch_support_value_seen = p.branch_support : null;
                max_branch_support_value_seen == null ? max_branch_support_value_seen = p.branch_support : null;
                p.branch_support > max_branch_support_value_seen ? max_branch_support_value_seen = p.branch_support : null;
                p.branch_support < min_branch_support_value_seen ? min_branch_support_value_seen = p.branch_support : null;
        }

        this.settings['show-support-values'] ? drawSupportValue(this.tree_svg_id, p, p0, p1, supportValueData) : null;
        drawCircleArc(this.tree_svg_id, p, p0, p1, p.radius, large_arc_flag);

        let arc = drawCircleArc(this.tree_svg_id, p, p0, p1, p.radius, large_arc_flag);
        arc.setAttribute('style', PADDING_STYLE);
        arc.classList.add('clone');

    }
    else
    {
        var p0 = [];
        var p1 = [];

        p0['x'] = p.xy['x'];
        p0['y'] = p.xy['y'];

        var anc = p.ancestor;
        if (anc) {
            p1['x'] = anc.xy['x'];
            p1['y'] = p0['y'];

            drawLine(this.tree_svg_id, p, p0, p1);

            let line = drawLine(this.tree_svg_id, p, p0, p1);
            line.setAttribute('style', PADDING_STYLE);
            line.classList.add('clone');
        }

        // vertical line
        var pl = p.child.xy;
        var pr = p.child.GetRightMostSibling().xy;

        p0['x'] = p0['x'];
        p0['y'] = pl['y'];
        p1['x'] = p0['x'];
        p1['y'] = pr['y'];

        drawLine(this.tree_svg_id, p, p0, p1, true);

        // support value business happens here:
        if (typeof p.branch_support === 'string' && p.branch_support.includes('/')) {
            multiple_support_value_seen = true;
            const [branch_support_value0, branch_support_value1] = p.branch_support.split('/').map(parseFloat);
            
            branch_support_values.push(branch_support_value0, branch_support_value1);

            if (branch_support_values.length > 0) {
                const min_support = Math.min(...branch_support_values);
                const max_support = Math.max(...branch_support_values);
            
                if (min_branch_support_value_seen === null || min_support < min_branch_support_value_seen) {
                    min_branch_support_value_seen = min_support;
                } else {
                    this.min_support = null;
                }
            
                if (max_branch_support_value_seen === null || max_support > max_branch_support_value_seen) {
                    max_branch_support_value_seen = max_support;
                } else {
                    this.max_support = null;
                }
            } 

        } else {
                // If there are no string branch support value like '100/100':
                min_branch_support_value_seen == null ? min_branch_support_value_seen = p.branch_support : null;
                max_branch_support_value_seen == null ? max_branch_support_value_seen = p.branch_support : null;
                p.branch_support > max_branch_support_value_seen ? max_branch_support_value_seen = p.branch_support : null;
                p.branch_support < min_branch_support_value_seen ? min_branch_support_value_seen = p.branch_support : null;
        }
        this.settings['show-support-values'] ? drawSupportValue(this.tree_svg_id, p, p0, p1, supportValueData) : null;

        let line = drawLine(this.tree_svg_id, p, p0, p1, true);
        line.setAttribute('style', PADDING_STYLE);
        line.classList.add('clone');
    }
};

Drawer.prototype.draw_collapsed_node = function(p, attributes) {
    var p0 = p.xy

    var triangle = document.createElementNS('http://www.w3.org/2000/svg', 'polygon');
    triangle.setAttribute('id', 'line' + p.id);
    triangle.setAttribute('vector-effect', 'non-scaling-stroke');
    triangle.setAttribute('fill', attributes['color']);
    triangle.setAttribute('style', 'stroke:' + LINE_COLOR + ';stroke-width:1;');

    var tp1_x, tp1_y, tp2_x, tp2_y;
    if (this.settings['tree-type'] == 'phylogram') {
        var tp1_x = this.root_length + (((p.path_length + p.max_child_path) / this.max_path_length) * ((this.height) - this.root_length));
        var tp2_x = tp1_x;

        var tp1_y = p0['y'] - p.size / 3;
        var tp2_y = p0['y'] + p.size / 3;

        drawRotatedText('tree_bin',
            {
                'x': tp1_x + 20,
                'y': (tp1_y + tp2_y) / 2
            },
            p.label,
            0,
            'left',
            ((attributes['font_size'] !== "0") ? parseFloat(attributes['font_size']) : (p.size / 2)) + 'px',
            'sans-serif',
            attributes['color'],
            0,
            'central');

    } else {
        var _radius = this.root_length + (((p.path_length + p.max_child_path) / this.max_path_length) * ((this.radius / 2) - this.root_length));
        let size = Math.min(p.size * 2/3, Math.PI / 3);
        var tp1_x = _radius * Math.cos(p.angle + size / 3);
        var tp1_y = _radius * Math.sin(p.angle + size / 3);

        var tp2_x = _radius * Math.cos(p.angle - size / 3);
        var tp2_y = _radius * Math.sin(p.angle - size / 3);

        let bottom_size = Math.hypot(tp1_x - tp2_x, tp1_y - tp2_y);

        let align = 'left';
        let angle = Math.toDegrees(p.angle);

        if ((p.angle > Math.PI / 2.0) && (p.angle < 1.5 * Math.PI)) {
            align = 'right';
            angle += 180.0;
        }

        drawRotatedText('tree_bin',
            {
                'x': (_radius + 20) * Math.cos(p.angle),
                'y': (_radius + 20) * Math.sin(p.angle)
            },
            p.label,
            angle,
            align,
            ((attributes['font_size'] !== "0") ? parseFloat(attributes['font_size']) : (bottom_size / 2)) + 'px',
            'sans-serif',
            attributes['color'],
            0,
            'central');
    }
    triangle.setAttribute('points', p0['x'] + ',' + p0['y'] + ' ' + tp1_x + ',' + tp1_y + ' ' + tp2_x + ',' + tp2_y);
    document.getElementById('tree').appendChild(triangle);
};

Drawer.prototype.draw_root = function() {
    var p0 = this.tree.root.xy
    var p1 = [];

    if (this.settings['tree-type'] == 'circlephylogram') {
        p1['x'] = 0;
        p1['y'] = 0;
    } else {
        p1['x'] = 0;
        p1['y'] = p0['y'];
    }

    drawLine(this.tree_svg_id, this.tree.root, p0, p1);
};

Drawer.prototype.initialize_screen = function() {
    this.width = parseFloat(this.settings['tree-width']);
    this.height = parseFloat(this.settings['tree-height']);
    this.radius = parseFloat(this.settings['tree-radius']);

    if (this.width == 0) {
        this.width = VIEWER_WIDTH;
    }

    if (this.height == 0) {
        this.height = VIEWER_HEIGHT;
    }

    if (this.radius == 0) {
        this.radius = (this.height > this.width) ? this.height : this.width;
    }

    // this global variable required for calculating the event position
    // look at getNodeFromEvents function in mouse-events.js
    original_width = this.width;
};

Drawer.prototype.draw_categorical_layers = function() {
    this.iterate_layers(function(layer) {
        if (layer.get_visual_attribute('height') == 0)
            return;

        if (!(layer.is_categorical || layer.is_parent))
            return;

        var layer_items = [];

        for (var node_i=0; node_i < this.tree.leaves.length; node_i++) {
            q = this.tree.leaves[node_i];

            if ((layer.is_categorical && layer.get_visual_attribute('type') == 'color') || layer.is_parent)
            {
                layer_items.push(this.layerdata_dict[q.label][layer.index]);
            }
            else
            {
                var _label = (this.layerdata_dict[q.label][layer.index] == null) ? '' : this.layerdata_dict[q.label][layer.index];

                if (this.settings['tree-type'] == 'circlephylogram')
                {
                    var align = 'left';
                    var new_angle = q.angle * 180.0 / Math.PI;
                    var offset_xy = [];
                    var _radius = this.layer_boundaries[layer.order][0] + this.layer_fonts[layer.order] * MONOSPACE_FONT_ASPECT_RATIO;
                    var font_gap = Math.atan(this.layer_fonts[layer.order] / _radius) / 3;

                    if ((q.angle > Math.PI / 2.0) && (q.angle < 1.5 * Math.PI)) {
                        align = 'right';
                        new_angle += 180.0;
                        offset_xy['x'] = Math.cos(q.angle - font_gap) * _radius;
                        offset_xy['y'] = Math.sin(q.angle - font_gap) * _radius;
                    } else {
                        offset_xy['x'] = Math.cos(q.angle + font_gap) * _radius;
                        offset_xy['y'] = Math.sin(q.angle + font_gap) * _radius;
                    }

                    drawRotatedText('layer_' + layer.order,
                                    offset_xy,
                                    _label,
                                    new_angle,
                                    align,
                                    this.layer_fonts[layer.order],
                                    "monospace",
                                    layer.get_visual_attribute('color'),
                                    layer.get_visual_attribute('height'),
                                    'center');
                }
                else
                {

                    var _offsetx = this.layer_boundaries[layer.order][0] + this.layer_fonts[layer.order] * MONOSPACE_FONT_ASPECT_RATIO;
                    var offset_xy = [];
                    offset_xy['x'] = _offsetx;
                    offset_xy['y'] = q.xy['y'] + this.layer_fonts[layer.order] / 4;

                    drawRotatedText('layer_' + layer.order,
                                    offset_xy,
                                    _label,
                                    0,
                                    'left',
                                    this.layer_fonts[layer.order],
                                    "monospace",
                                    layer.get_visual_attribute('color'),
                                    layer.get_visual_attribute('height'),
                                    'center');
                }
            }
        }

        if ((layer.is_categorical && layer.get_visual_attribute('type') == 'color') || layer.is_parent) {
            // This algorithm finds changing points in the categorical layer
            // to draw sequential items with single object, so we have to add
            // a dummy item to the end in order to detect last group.
            // I have used 'null', but if the last group is already 'null', algorithm skips
            // that group since there is no change, so we have to check the last item
            // and if it is 'null' we have to add something different.

            if (layer_items[layer_items.length-1] == null) {
                layer_items.push(-1);
            } else {
                layer_items.push(null)
            }

            var prev_value = layer_items[0];
            var prev_start = 0;

            var items_to_draw = new Array();

            for (var j=1; j < layer_items.length; j++)
            {
                if (prev_value != layer_items[j])
                {
                    if (layer.is_categorical || prev_value != '')
                        items_to_draw.push(new Array(prev_start, j - 1, prev_value)); // start, end, item;

                    prev_start = j;
                }
                prev_value = layer_items[j];
            }

            for (var j=0; j < items_to_draw.length; j++)
            {
                var categorical_item = items_to_draw[j];

                var color;

                if (layer.is_categorical)
                {
                    var _category_name = categorical_item[2];
                    if (_category_name == null || _category_name == '' || _category_name == 'null')
                        _category_name = 'None';
                    color = categorical_data_colors[layer.index][_category_name];
                }
                else // parent
                {
                    if (j % 2 == 1 && j == items_to_draw.length - 1)
                        color = '#AAAAAA';
                    else if (j % 2 == 1)
                        color = '#888888';
                    else
                        color = '#666666';
                }

                var start = this.tree.leaves[categorical_item[0]];
                var end = this.tree.leaves[categorical_item[1]];

                if (tree_type == 'circlephylogram')
                {
                    drawPie('layer_' + layer.order,
                        'categorical_' + layer.order + '_' + j, // path_<layer>_<id>
                        start.angle - start.size / 2,
                        end.angle + end.size / 2,
                        this.layer_boundaries[layer.order][0],
                        this.layer_boundaries[layer.order][1],
                        (end.angle - start.angle + (end.size / 2) + (start.size / 2) > Math.PI) ? 1 : 0,
                        color,
                        1,
                        false);
                }
                else
                {
                    var height = end.xy['y'] + end.size / 2 - start.xy['y'] + start.size / 2;

                    drawPhylogramRectangle('layer_' + layer.order,
                        'categorical_' + layer.order + '_' + j, // path_<layer>_<id>
                        this.layer_boundaries[layer.order][0],
                        start.xy['y'] - start.size / 2 + height / 2,
                        height,
                        this.layer_boundaries[layer.order][1] - this.layer_boundaries[layer.order][0],
                        color,
                        1,
                        false);
                }
            }
        }
    }.bind(this));
};

Drawer.prototype.calculate_layer_boundaries = function() {
    this.layer_boundaries = new Array();

    if (this.settings['tree-type'] == 'phylogram') {
        if (this.has_tree)
        {
            this.layer_boundaries.push([0, this.height]);
        }
        else
        {
            this.layer_boundaries.push([0, 0]);
        }
    }
    else
    {
        this.layer_boundaries.push([0, parseFloat(this.radius) / 2]);
    }

    this.iterate_layers(function(layer) {
        if (layer.is_categorical && layer.get_visual_attribute('type') == 'text' && layer.get_visual_attribute('height') == 0) {
            this.calculate_font_size_for_text_layer(layer);
        }

        var margin = parseFloat(this.settings['layer-margin']) ? parseFloat(layer.get_visual_attribute('margin')) : parseFloat(this.settings['layer-margin']);
        var height = parseFloat(layer.get_visual_attribute('height'));

        var ending_of_previous_layer = this.layer_boundaries[layer.order - 1][1];

        var layer_start = ending_of_previous_layer + ((height > 0) ? margin : 0);
        var layer_end   = layer_start + height;

        this.layer_boundaries.push([layer_start, layer_end]);
    }.bind(this));
};

Drawer.prototype.calculate_font_size_for_text_layer = function(layer) {
    var leaf_perimeter;
    if (this.settings['tree-type'] == 'circlephylogram') {
        var margin = (this.settings['custom-layer-margin']) ? parseFloat(layer.get_visual_attribute('margin')) : parseFloat(this.settings['layer-margin']);
        var leaf_perimeter = this.smallest_leaf_size * (this.layer_boundaries[layer.order - 1][1] + margin);
    } else {
        var leaf_perimeter = this.smallest_leaf_size;
    }
    var layer_font = Math.min(leaf_perimeter, parseFloat(this.settings['max-font-size']));

    this.layer_fonts[layer.order] = layer_font;

    if (layer.is_categorical && layer.get_visual_attribute('type') == 'text' && layer.get_visual_attribute('height') == 0)
    {
        // find longest item.
        var longest_text_len = 0;
        for (var _pos = 1; _pos < layerdata.length; _pos++)
        {
            if (layerdata[_pos][layer.index] != null && layerdata[_pos][layer.index].length > longest_text_len)
            {
                longest_text_len = layerdata[_pos][layer.index].length;
            }
        }
        // make layer height bit longer than text
        longest_text_len += 2;

        layers[layer.index]['height'] = Math.ceil(longest_text_len * MONOSPACE_FONT_ASPECT_RATIO * layer_font) + 1;
        $('#height' + layer.index).val(layers[layer.index]['height']);
    }
};

Drawer.prototype.draw_layer_backgrounds = function() {
    this.iterate_layers(function(layer) {
        createBin('tree_bin', 'layer_background_' + layer.order);
        createBin('tree_bin', 'layer_' + layer.order);
        createBin('tree_bin', 'event_catcher_' + layer.order);

        // draw event catcher of the layer
        if (this.settings['tree-type']=='phylogram')
        {
            drawPhylogramRectangle('event_catcher_' + layer.order,
                'event',
                this.layer_boundaries[layer.order][0],
                this.width / 2,
                this.width,
                this.layer_boundaries[layer.order][1] - this.layer_boundaries[layer.order][0],
                '#ffffff',
                0,
                true);
        }
        else
        {
            var _first = this.tree.leaves[0];
            var _last = this.tree.leaves[this.tree.leaves.length - 1];
            drawPie('event_catcher_' + layer.order,
                'event',
                _first.angle - _first.size / 2,
                _last.angle + _last.size / 2,
                this.layer_boundaries[layer.order][0],
                this.layer_boundaries[layer.order][1],
                (_last.angle - _first.angle + _first.size / 2 + _last.size / 2 > Math.PI) ? 1:0, // large arc flag
                '#ffffff',
                0,
                true);
        }

        var _bgcolor;
        var _opacity;

        if (layer.is_categorical && layer.get_visual_attribute('type') == 'text') {
            _bgcolor = layer.get_visual_attribute('color-start');
            _opacity = 1;
        } else {
            _bgcolor = layer.get_visual_attribute('color');
            _opacity = parseFloat(this.settings['background-opacity']);
        }

        // draw backgrounds
        if (this.settings['tree-type']=='phylogram' && ((layer.is_numerical && layer.get_visual_attribute('type') == 'bar') || (layer.is_categorical && layer.get_visual_attribute('type') == 'text')))
        {
            drawPhylogramRectangle('layer_background_' + layer.order,
                'all',
                this.layer_boundaries[layer.order][0],
                this.width / 2,
                this.width,
                this.layer_boundaries[layer.order][1] - this.layer_boundaries[layer.order][0],
                _bgcolor,
                _opacity,
                false);
        }


        if (this.settings['tree-type'] == 'circlephylogram' && ((layer.is_numerical && layer.get_visual_attribute('type') == 'bar') || (layer.is_categorical && layer.get_visual_attribute('type') == 'text')))
        {
            var _first = this.tree.leaves[0];
            var _last = this.tree.leaves[this.tree.leaves.length - 1];

            drawPie('layer_background_' + layer.order,
                'all',
                _first.angle - _first.size / 2,
                _last.angle + _last.size / 2,
                this.layer_boundaries[layer.order][0],
                this.layer_boundaries[layer.order][1],
                (Math.abs(_last.angle - _first.angle) + _first.size / 2 + _last.size / 2 > Math.PI) ? 1:0, // large arc flag
                _bgcolor,
                _opacity,
                false);
        }
    }.bind(this));
};

Drawer.prototype.draw_guide_lines = function() {
    if (!this.has_tree)
        return;

    if (this.settings['draw-guide-lines'] == 'no')
        return;

    var beginning_of_layers = this.layer_boundaries[0][1];
    var odd_even_flag = 1;

    var n = new NodeIterator(this.tree.root);
    var q = n.Begin();

    while (q != null) {
        if (q.IsLeaf() && !q.collapsed) {
            if (this.settings['draw-guide-lines'] == 'odd_even') {
                odd_even_flag = odd_even_flag * -1;
            }
            if (odd_even_flag > 0) {
                if (this.settings['tree-type'] == 'phylogram') {
                    drawStraightGuideLine('guide_lines', q.id, q.xy, beginning_of_layers);
                } else {
                    drawGuideLine('guide_lines', q.id, q.angle, q.radius, beginning_of_layers);
                }

            }
        }
        q = n.Next();
    }
};

Drawer.prototype.draw_numerical_layers = function() {
    this.iterate_layers(function(layer) {
        if (layer.get_visual_attribute('height') == 0)
            return;

        if (!layer.is_numerical)
            return;

        var numeric_cache = [];

        var previous_non_zero_order = 0;

        for (var i=0; i < this.tree.leaves.length; i++) {
            q = this.tree.leaves[i];

            if (this.settings['tree-type'] == 'phylogram') {
                if (layer.get_visual_attribute('type') == 'intensity') {
                     drawPhylogramRectangle('layer_' + layer.order,
                        q.id,
                        this.layer_boundaries[layer.order][0] ,
                        q.xy['y'],
                        q.size,
                        this.layer_boundaries[layer.order][1] - this.layer_boundaries[layer.order][0],
                        getGradientColor(
                            layer.get_visual_attribute('color-start'),
                            layer.get_visual_attribute('color'),
                            this.layerdata_dict[q.label][layer.index] / layer.get_visual_attribute('height')
                        ),
                        1,
                        false);
                }
                else
                {
                    if (this.settings['optimize-speed'] || layer.get_visual_attribute('type') == 'line') {
                        if (this.layerdata_dict[q.label][layer.index] > 0)
                        {
                            if (q.order == 0 || (q.order > 0 && this.layerdata_dict[this.tree.leaves[i-1].label][layer.index] == 0)) {
                                numeric_cache.push(
                                    "M",
                                    this.layer_boundaries[layer.order][1],
                                    q.xy['y'] - q.size / 2
                                    );
                            }

                            if (layer.get_visual_attribute('type') == 'line')
                            {
                                numeric_cache.push(
                                    "L",
                                    this.layer_boundaries[layer.order][1] - this.layerdata_dict[q.label][layer.index],
                                    q.xy['y']
                                    );
                            }
                            else
                            {
                                numeric_cache.push(
                                    "L",
                                    this.layer_boundaries[layer.order][1] - this.layerdata_dict[q.label][layer.index],
                                    q.xy['y'] - q.size / 2,
                                    "L",
                                    this.layer_boundaries[layer.order][1] - this.layerdata_dict[q.label][layer.index],
                                    q.xy['y'] + q.size / 2
                                    );
                            }

                            if ((q.order == this.tree.leaves.length - 1) || (q.order < (this.tree.leaves.length - 1) && this.layerdata_dict[this.tree.leaves[i+1].label][layer.index] == 0)) {
                                numeric_cache.push(
                                    "L",
                                    this.layer_boundaries[layer.order][1],
                                    q.xy['y'] + q.size / 2,
                                    "Z"
                                    );
                            }
                        }
                    }
                    else
                    {
                        if (this.layerdata_dict[q.label][layer.index] > 0) {
                             drawPhylogramRectangle('layer_' + layer.order,
                                q.id,
                                this.layer_boundaries[layer.order][1] - this.layerdata_dict[q.label][layer.index],
                                q.xy['y'],
                                q.size,
                                this.layerdata_dict[q.label][layer.index],
                                layer.get_visual_attribute('color'),
                                1,
                                false);
                        }
                    }
                }
            }
            else
            {
                if (layer.get_visual_attribute('type') == 'intensity') {
                    drawPie('layer_' + layer.order,
                        q.id,
                        q.angle - q.size / 2,
                        q.angle + q.size / 2,
                        this.layer_boundaries[layer.order][0],
                        this.layer_boundaries[layer.order][1],
                        (Math.abs(q.size) > Math.PI) ? 1 : 0,
                        getGradientColor(
                            layer.get_visual_attribute('color-start'),
                            layer.get_visual_attribute('color'),
                            this.layerdata_dict[q.label][layer.index] / layer.get_visual_attribute('height')
                        ),
                        1,
                        false);
                }
                else
                {
                    if (this.settings['optimize-speed'] || layer.get_visual_attribute('type') == 'line') {
                        var start_angle  = q.angle - q.size / 2;
                        var end_angle    = q.angle + q.size / 2;
                        var inner_radius = this.layer_boundaries[layer.order][0];
                        var outer_radius = this.layer_boundaries[layer.order][0] + this.layerdata_dict[q.label][layer.index];

                        if (this.layerdata_dict[q.label][layer.index] > 0) {
                            if (q.order == 0 || this.layerdata_dict[this.tree.leaves[i-1].label][layer.index] == 0)
                            {
                                var ax = Math.cos(start_angle) * inner_radius;
                                var ay = Math.sin(start_angle) * inner_radius;

                                numeric_cache.push("M", ax, ay);
                                previous_non_zero_order = q.order;
                            }

                            if (layer.get_visual_attribute('type') == 'line') {
                                var bx = Math.cos(q.angle) * outer_radius;
                                var by = Math.sin(q.angle) * outer_radius;

                                numeric_cache.push("L", bx, by);
                            }
                            else
                            {
                                var bx = Math.cos(start_angle) * outer_radius;
                                var by = Math.sin(start_angle) * outer_radius;

                                var cx = Math.cos(end_angle) * outer_radius;
                                var cy = Math.sin(end_angle) * outer_radius;

                                numeric_cache.push("L", bx, by, "A", outer_radius, outer_radius, 0, is_large_angle(start_angle, end_angle), 1, cx, cy);
                            }

                            if ((q.order == this.tree.leaves.length - 1) || this.layerdata_dict[this.tree.leaves[i+1].label][layer.index] == 0) {
                                var dx = Math.cos(end_angle) * inner_radius;
                                var dy = Math.sin(end_angle) * inner_radius;

                                numeric_cache.push("L", dx, dy);
                                var first_node = this.tree.leaves[previous_non_zero_order];
                                var first_node_start_angle = first_node.angle - first_node.size / 2;

                                var ex = Math.cos(first_node_start_angle) * inner_radius;
                                var ey = Math.sin(first_node_start_angle) * inner_radius;

                                numeric_cache.push("A", inner_radius, inner_radius, 1, is_large_angle(end_angle, first_node_start_angle), 0, ex, ey, "Z");
                            }
                        }
                    }
                    else
                    {
                        if (layer.get_visual_attribute('height') > 0) {
                            drawPie('layer_' + layer.order,
                                q.id,
                                q.angle - q.size / 2,
                                q.angle + q.size / 2,
                                this.layer_boundaries[layer.order][0],
                                this.layer_boundaries[layer.order][0] + this.layerdata_dict[q.label][layer.index],
                                0,
                                layer.get_visual_attribute('color'),
                                1,
                                false);
                        }
                    }
                }
            }
        }

        if (numeric_cache.length > 0) {
            var path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
            path.setAttribute('shape-rendering', 'auto');
            path.setAttribute('pointer-events', 'none');

            if (layer.get_visual_attribute('type') == 'line') {
                path.setAttribute('stroke', layer.get_visual_attribute('color-start'));
                path.setAttribute('stroke-width', '1px');
                path.setAttribute('vector-effect', 'non-scaling-stroke');
                path.setAttribute('fill', layer.get_visual_attribute('color'));
            }
            else {
                path.setAttribute('fill', layer.get_visual_attribute('color'));
                path.setAttribute('stroke-width', '0');
            }

            path.setAttribute('d', numeric_cache.join(' '));
            var layer_group = document.getElementById('layer_' + layer.order);
            layer_group.appendChild(path);
        }

    }.bind(this));
};

Drawer.prototype.draw_stack_bar_layers = function() {
    this.iterate_layers(function(layer) {
        if (layer.get_visual_attribute('height') == 0)
            return;

        if (!layer.is_stackbar)
            return;

        var path_cache = [];

        for (var i=0; i < this.tree.leaves.length; i++) {
            q = this.tree.leaves[i];
            var offset = 0;
            for (var j=0; j < this.layerdata_dict[q.label][layer.index].length; j++)
            {
                let stack_item_value = parseFloat(this.layerdata_dict[q.label][layer.index][j]);
                if (isNaN(stack_item_value)) {
                    stack_item_value = Number.EPSILON; // smallest positive number
                }

                if (this.settings['tree-type'] == 'phylogram') {
                    if (q.order == 0) {
                        path_cache[j] = [];
                        path_cache[j].push(
                            "M",
                            this.layer_boundaries[layer.order][1],
                            q.xy['y'] - q.size / 2
                            );
                    }

                    path_cache[j].push(
                        "L",
                        this.layer_boundaries[layer.order][1] - stack_item_value - offset,
                        q.xy['y'] - q.size / 2,
                        "L",
                        this.layer_boundaries[layer.order][1] - stack_item_value - offset,
                        q.xy['y'] + q.size / 2
                        );

                    if (q.order == (this.tree.leaves.length - 1)) {
                        path_cache[j].push(
                            "L",
                            this.layer_boundaries[layer.order][1],
                            q.xy['y'] + q.size / 2,
                            "Z"
                            );
                    }
                }
                else
                {
                    var start_angle  = q.angle - q.size / 2;
                    var end_angle    = q.angle + q.size / 2;
                    var inner_radius = this.layer_boundaries[layer.order][0];
                    var outer_radius = this.layer_boundaries[layer.order][0] + stack_item_value + offset;

                    if (q.order == 0)
                    {
                        path_cache[j] = [];

                        var ax = Math.cos(start_angle) * inner_radius;
                        var ay = Math.sin(start_angle) * inner_radius;

                        path_cache[j].push("M", ax, ay);
                    }

                    var bx = Math.cos(start_angle) * outer_radius;
                    var by = Math.sin(start_angle) * outer_radius;

                    var cx = Math.cos(end_angle) * outer_radius;
                    var cy = Math.sin(end_angle) * outer_radius;

                    path_cache[j].push("L", bx, by, "A", outer_radius, outer_radius, 0, is_large_angle(start_angle, end_angle), 1, cx, cy);

                    if (q.order == this.tree.leaves.length - 1) {
                        var bx = Math.cos(end_angle) * inner_radius;
                        var by = Math.sin(end_angle) * inner_radius;

                        var _min = Math.toRadians(this.settings['angle-min']);
                        var _max = Math.toRadians(this.settings['angle-max']);

                        path_cache[j].push("L", bx, by,
                            "A", inner_radius, inner_radius, 0, is_large_angle(_min, _max), 0, path_cache[j][1], path_cache[j][2],
                            "Z");
                    }
                }

                offset += stack_item_value;
            }
        }

        let layer_name = getLayerName(layer.index);
        let bars = (layer_name.indexOf('!') > -1) ? layer_name.split('!')[1].split(';') : layer_name.split(';');

        for (var j = path_cache.length - 1; j >= 0; j--) {
            var path = document.createElementNS('http://www.w3.org/2000/svg', 'path');

            path.setAttribute('shape-rendering', 'auto');
            path.setAttribute('pointer-events', 'none');
            path.setAttribute('fill', stack_bar_colors[layer.index][bars[j]]);
            path.setAttribute('stroke-width', '0');
            path.setAttribute('d', path_cache[j].join(' '));

            var layer_group = document.getElementById('layer_' + layer.order);
            layer_group.appendChild(path);
        }
    }.bind(this));
};

Drawer.prototype.initialize_tree_observer = function() {
    var observer = new MutationObserver(updateSingleBackgroundGlobals); // in mouse-events.js

    observer.observe(document.getElementById('viewport'), {
        attributes:    true,
        attributeFilter: ["transform"]
    });
};

Drawer.prototype.draw_layer_names = function() {

    createBin('tree_bin', 'layer_labels');

    this.iterate_layers(function(layer) {
        var height = parseFloat(layer.get_visual_attribute('height'));

        if (height == 0)
            return;

        var layer_title = " " + getPrettyLayerTitle(layerdata[0][layer.index]);
        var font_size = Math.min(height, parseFloat(this.settings['max-font-size-label']));

        if (this.settings['tree-type'] == 'circlephylogram')
        {
            if (parseFloat(this.settings['angle-max']) > 315)
                return;

            var angle_max = Math.toRadians(parseFloat(this.settings['angle-max']));
            var distance = this.layer_boundaries[layer.order][0] + (height / 2) - (font_size * 1/3);

            var cx = Math.cos(angle_max) * distance;
            var cy = Math.sin(angle_max) * distance;

            var textobj = drawFixedWidthText('layer_labels', {
                    'x': cx,
                    'y': cy
                },
                layer_title,
                font_size + 'px',
                layer.get_visual_attribute('color'),
                total_radius,
                height);

            textobj.setAttribute('transform', 'rotate(' + Math.toDegrees(angle_max + Math.PI / 2) + ' ' + cx + ' ' + cy + ')');
        }
        else
        {
            drawRotatedText('layer_labels', {
                'x': this.layer_boundaries[layer.order][1] - (height / 2) + (font_size * 1/3),
                'y': this.width + 40, 
                },
                layer_title,
                -90,
                'right',
                font_size + 'px',
                'sans-serif',
                layer.get_visual_attribute('color'),
                0,
                'baseline');
        }
    }.bind(this));
};

Drawer.prototype.draw_samples = function() {
    if (this.settings['tree-type'] == 'circlephylogram' && this.settings['angle-min'] == 0 && this.settings['angle-max'] <= 270)
    {
        createBin('viewport', 'samples');
        drawSamples(); // in sample.js
    }
    else if (this.settings['tree-type'] == 'phylogram')
    {
        createBin('viewport', 'samples');
        $('#samples').attr('transform', 'rotate(90)');
        drawSamples(); // in sample.js
    }
};

Drawer.prototype.update_title_panel = function() {
    var _sub_title = "Items order: <b>" + getClusteringPrettyName(this.settings['order-by']) + "</b> | ";
        _sub_title += "Current view: <b>" + this.settings['current-view'] + "</b> | ";
        _sub_title += "Sample order: <b>" + this.settings['samples-order'] + "</b>";
    $('#title-panel-second-line').html(_sub_title);
};

Drawer.prototype.show_drawing_statistics = function() {
    var tree_object_count = document.getElementById('tree').getElementsByTagName('*').length + document.getElementById('guide_lines').getElementsByTagName('*').length;
    var total_object_count = document.getElementById('svg').getElementsByTagName('*').length;

    $('#draw_delta_time').html(this.tree.leaves.length + ' splits and ' + total_object_count +' objects drawn in ' + this.timer.getDeltaSeconds('done')['deltaSecondsStart'] + ' seconds.');

    console.log('[info] Leaf count: ' + this.tree.leaves.length);
    console.log('[info] Object count in tree (with guide lines): ' + tree_object_count);
    console.log('[info] Total objects in SVG: ' + total_object_count);
};
