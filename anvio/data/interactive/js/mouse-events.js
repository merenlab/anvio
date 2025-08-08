/**
 * Javascript user input handler functions for anvi'o interactive interface
 *
 *  Authors: Ozcan Esen
 *           Matthew Klein <mtt.l.kln@gmail.com>
 *           A. Murat Eren <a.murat.eren@gmail.com>
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

$(document).ready(function() {
    document.body.addEventListener('mousemove', mouseMoveHandler, false); // for tooltip
});

function lineClickHandler(event) {
    if (dragging || drawing_zoom)
        return;

    if (event.target.parentNode && event.target.parentNode.id == 'samples_tree')
    {
        var id = event.target.id.match(/\d+/);
        var node = samples_id_to_node_map[id[0]];

        var _n = new NodeIterator(node);
        var _q = _n.Begin();

        $('#table_layers').find('.layer_selectors').prop('checked', false);
        while (_q != null)
        {
            if (_q.IsLeaf()) {

                if(_q.label){
                    $('#table_layers').find('.titles').each(
                        function(index, obj){
                            if (_q.label.toLowerCase() == obj.title.toLowerCase())
                            {
                                $(obj).parent().find('.layer_selectors').prop('checked','checked');
                            }
                        }
                    );
                }
            }
            _q = _n.Next();
        }
        return;
    }

    if (IsCtrlPressed(event)) {
        bins.NewBin();
    }

    var p = getNodeFromEvent(event);
    bins.AppendNode(p);
}


function lineContextMenuHandler(event) {

    if (event.preventDefault) {
        event.preventDefault();
    }
    if (event.target.parentNode.id == 'samples_tree') {
        var id = event.target.id.match(/\d+/);
        var p = samples_id_to_node_map[id[0]];
        let allLeaves = drawer.tree.leaves
        console.log(allLeaves)

        let menu = new ContextMenu({'container': document.body,
                                    'event': event,
                                    'node': p,
                                    'layer': null,
                                    'all' : allLeaves,
                                    'isSample': true});

        menu.Show();
    } else {
        let p = getNodeFromEvent(event);
        let layer_id = event.target.parentNode.id.match(/\d+/);
        let allLeaves = drawer.tree.leaves

        if (p.IsLeaf() || IsCtrlPressed(event)) {
            let menu = new ContextMenu({'container': document.body,
                                        'event': event,
                                        'node': p,
                                        'layer': layer_id,
                                        'all' : allLeaves,
                                        'isSample': false});

            menu.Show();
        }
        else
        {
            bins.RemoveNode(p);
        }
    }
}


function lineMouseEnterHandler(event) {
    if (drawing_zoom)
        return;

    var p = getNodeFromEvent(event);

    $('#path_hover').remove();

    if (typeof p === 'undefined' || p.id == 0)
        return; // skip root

    if (p.collapsed)
        return;

    var bin_color = bins.GetSelectedBinColor();
    var [p1, p2] = p.GetBorderNodes();

    if (tree_type == 'circlephylogram')
    {
        drawPie('tree_bin',
            'hover',
            p1.angle - p1.size / 2,
            p2.angle + p2.size / 2,
            distance(p.backarc, {
                'x': 0,
                'y': 0
            }),
            total_radius,
            (p2.angle - p1.angle + (p1.size / 2) + (p2.size / 2) > Math.PI) ? 1 : 0,
            bin_color,
            0.3,
            false);
    }
    else
    {
        var _x = (p.ancestor) ? p.ancestor.xy.x : p.xy.x;
        var begin = p1.xy.y - p1.size / 2;
        var end = p2.xy.y + p2.size / 2;

        drawPhylogramRectangle('tree_bin',
            'hover',
             _x,
            (begin + end) / 2,
            end - begin,
            Math.abs(total_radius - _x),
            bin_color,
            0.3,
            false);
    }
}

function lineMouseLeaveHandler(event) {
    if (drawing_zoom)
        return;

    $('#path_hover').remove();
}

function mouseMoveHandler(event) {

    if (drawing_zoom)
        return;

    if (event.target.parentNode.id == 'samples_tree' || samples_tree_hover)
    {
        var samples_tree = document.getElementById('samples_tree');
        for (var i=0; i < samples_tree.childNodes.length; i++)
        {
            var obj = samples_tree.childNodes[i];
            if (obj.className == 'clone') {
                continue;
            }
            obj.style['stroke'] = LINE_COLOR;
            obj.style['stroke-width'] = '1px';

        }
        samples_tree_hover = false;
    }

    if (event.target.parentNode.id == 'samples_tree')
    {

        var id = event.target.id.match(/\d+/);
        var node = samples_id_to_node_map[id[0]];
        var _n = new NodeIterator(node);
        var _q = _n.Begin();

        write_mouse_table(`<tr><td>Label</td><td>${node.label ? node.label : 'N/A'}</td></tr>
                           <tr><td>Support</td><td>${node.branch_support}</td></tr>
                           <tr><td>Edge length</td><td>${node.original_edge_length}</td></tr>`, 'Layers order branch', '', 0);

        while (_q != null)
        {
            var lineobj = document.getElementById('samples_line' + _q.id);
            if (lineobj)
            {
                lineobj.style['stroke-width'] = '3px';
                lineobj.style['stroke'] = '#FF0000';
            }
            var lineobj = document.getElementById('samples_arc' + _q.id);
            if (lineobj)
            {
                lineobj.style['stroke-width'] = '3px';
                lineobj.style['stroke'] = '#FF0000';
            }
            _q = _n.Next();
        }
        samples_tree_hover = true;
        return;
    }

    if (event.target.id == 'path_samples')
    {
        // samples tooltip
        var sample_name = event.target.getAttribute('sample-name');
        var sample_group = event.target.getAttribute('sample-group');
        var layer_name_hover = event.target.getAttribute('layer-name');
        var stack_bar_name_hover = event.target.getAttribute('bar-name');

        var message = "";
        var layer_pos = 0;
        var layer_counter=0;
        let layer_stackbar_pos = -1;

        for (var i=0; i < last_settings['samples-layer-order'].length; i++)
        {
            var layer_name = last_settings['samples-layer-order'][i]['layer_name'];
            var group      = last_settings['samples-layer-order'][i]['group'];

            if (!is_sample_group_visible(group)) {
                continue;
            }

            var pretty_name = (layer_name.indexOf('!') > -1) ? layer_name.split('!')[0] : layer_name;

            let highlight_row = false;
            if (layer_name == layer_name_hover && group == sample_group)
            {
                highlight_row = true;
                layer_pos = layer_counter;
            }

            if (layer_name.indexOf('!') > -1) { // stack bar
                message += `<tr style="${highlight_row ? 'background-color: #dddddd;' : ''}"><td>` + pretty_name + '</td><td>';

                let stack_names = layer_name.split('!')[1].split(';');
                let stack_items = samples_information_dict[group][sample_name][layer_name].split(';');
                message += '<table class="table table-striped">';
                for (let j = stack_names.length - 1; j >= 0; j--) {
                    let bar_name = stack_names[j];
                    let bar_pretty_name = bar_name.replace('Unknown_t_', '').replace('_', ' ');
                    message += `<tr>
                                    <td style="${(highlight_row && bar_name == stack_bar_name_hover) ? 'background-color: rgb(232, 202, 207);' : ''}">
                                        <div class="colorpicker" style="background-color: ${samples_stack_bar_colors[group][layer_name][bar_name]}"></div>&nbsp;${bar_pretty_name}
                                    </td>
                                    <td style="white-space: nowrap; ${(highlight_row && bar_name == stack_bar_name_hover) ? 'background-color: rgb(232, 202, 207);' : ''}">${stack_items[j]}</td>
                                </tr>`;

                    if (highlight_row && bar_name == stack_bar_name_hover) {
                        layer_stackbar_pos = stack_names.length - j;
                    }
                }
                message += '</table>';
            } else {
                message += `<tr style="${highlight_row ? 'background-color: rgb(232, 202, 207);' : ''}"><td>` + pretty_name + '</td><td>';
                message += samples_information_dict[group][sample_name][layer_name];
            }

            message += '</td></tr>';

            // since we skip hidden layer groups, we can not use 'i' to refer layer position.
            layer_counter++;
        }

        write_mouse_table(message,
            sample_name,
            (layer_name_hover.indexOf('!') > -1) ? layer_name_hover.split('!')[0] : layer_name_hover,
            layer_pos,
            layer_stackbar_pos);

        return;
    }

    var p = getNodeFromEvent(event);

    if (event.target.id && event.target.id == 'path_event')
        lineMouseEnterHandler(event);

    if (!p)
        return;

    if (!p.IsLeaf()) {
        write_mouse_table(`<tr><td>Label</td><td>${p.label ? p.label : 'N/A'}</td></tr>
                           <tr><td>Branch support</td><td>${p.branch_support}</td></tr>
                           <tr><td>Edge length</td><td>${p.original_edge_length}</td></tr>`, 'Branch', '', 0);
    }

    var layer_id_exp = event.target.parentNode.id.match(/\d+/);
    if (!layer_id_exp)
        return;
    var layer_id = layer_id_exp[0];
    var target_node = drawer.tree.nodes[p.id];

    if (target_node.collapsed) {
        write_mouse_table(target_node.label, 'CollapsedNode', '', 0);
        return;
    }

    var tooltip_arr = layerdata_title[target_node.label].slice(0);

    var message = "";
    for (var i=0; i < tooltip_arr.length; i++)
    {
        if (i == layer_id - 1)
        {
            message += '<tr style="background-color: rgb(232, 202, 207);">' + tooltip_arr[i] + '</tr>';
        }
        else
        {
            message += '<tr>' + tooltip_arr[i] + '</tr>';
        }
    }

    write_mouse_table(message,
        target_node.label,
        getPrettyLayerTitle(tooltip_arr[layer_id - 1].split('>')[1].split('<')[0]),
        layer_id);
}


function write_mouse_table(content, item_name, layer_name, layer_id, stackbar_layer_id) {
    $('#cell_item_name').html(item_name);

    if (layer_name && layer_name.length > 0) {
        $('#cell_layer_name').closest('tr').show();
        $('#cell_layer_name').html(layer_name);
    } else {
        $('#cell_layer_name').closest('tr').hide();
        $('#cell_layer_name').html('');
    }

    $('#tooltip_content').html(content);

    if ($('#tooltip_content').height() + 300 > $(window).height()) {
        let top = 0;
        if (stackbar_layer_id && stackbar_layer_id > -1) {
            top = $('#tooltip_content tr').eq(layer_id).find('tr').eq(stackbar_layer_id).position()['top'];
        } else {
            top = $('#tooltip_content tr').eq(layer_id).position()['top'];
        }

        $('#mouse_hover_scroll').css('top', Math.min(0, ($(window).height()-300) / 2 + -1 * top));
    } else {
        $('#mouse_hover_scroll').css('top', 0);
    }
}

// globals related single background
var rect_left;
var rect_width;
var original_width; // this will be set in Drawer.initialize_screen;

var origin_x;
var origin_y;

function updateSingleBackgroundGlobals()
{
    if (last_settings['tree-type'] == 'phylogram')
    {
        var path_event = document.getElementById('path_event');
        var rect = path_event.getBoundingClientRect();

        rect_left = rect.left;
        rect_width = rect.width;
    }
    else if (last_settings['tree-type'] == 'circlephylogram')
    {
        var root = document.getElementById('line_origin');
        var rect = root.getBoundingClientRect();

        origin_x = rect.left;
        origin_y = rect.top;
    }
}

function getNodeFromEvent(event)
{
    if (!drawer)
        return null;

    if (event.target.id == 'path_event' || event.target.parentNode.id == 'bin')
    {
        if (last_settings['tree-type'] == 'phylogram')
        {
            var _x = original_width - ((event.clientX - rect_left) * (window['original_width'] / rect_width));

            for (var i=0; i < drawer.tree.leaves.length; i++) {
                var node = drawer.tree.leaves[i];
                if ((_x > (node.xy['y'] - node.size / 2)) && (_x < (node.xy['y'] + node.size / 2))) {
                    return node;
                }
            }
        }
        else if (last_settings['tree-type'] == 'circlephylogram')
        {
            var _y = event.clientY - origin_y;
            var _x = event.clientX - origin_x;
            var angle = Math.atan2(_y, _x);// - angle_per_leaf / 2;
            if (angle < 0)
                angle = 2 * Math.PI + angle;

            for (var i=0; i < drawer.tree.leaves.length; i++) {
                var node = drawer.tree.leaves[i];
                if ((angle > (node.angle - node.size / 2)) && (angle < (node.angle + node.size / 2))) {
                    return node;
                }
            }
        }
    }
    else
    {
        var id = event.target.id.match(/\d+/);

        if (id)
            return drawer.tree.nodes[id[0]];
    }
}
