/**
 * Javascript user input handler functions for anvi'o interactive interface
 *
 *  Author: Özcan Esen <ozcanesen@gmail.com>
 *  Copyright 2015, The anvio Project
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


function getBinId() {
    var radios = document.getElementsByName('active_bin');
    for(var i=0; i < radios.length; i++)
    {
        if (radios[i].checked)
            return radios[i].value;
    }
}

function lineClickHandler(event) {
    if (dragging || drawing_zoom)
        return;
    
    var p = getNodeFromEvent(event);

    if (p.id == 0)
        return; // skip root

    if ((navigator.platform.toUpperCase().indexOf('MAC')>=0 && event.metaKey) || event.ctrlKey)
        newBin();

    var bin_id = getBinId();

    if (bin_id === 'undefined')
        return;

    var bin_color = document.getElementById('bin_color_' + bin_id).getAttribute('color');

    var bins_to_update = [];
    for (var i = 0; i < p.child_nodes.length; i++) {
        var pos = SELECTED[bin_id].indexOf(id_to_node_map[p.child_nodes[i]].label);
        if (pos == -1) {
            SELECTED[bin_id].push(id_to_node_map[p.child_nodes[i]].label);

            if (bins_to_update.indexOf(bin_id) == -1)
                bins_to_update.push(bin_id);
        }

        // remove nodes from other bins
        for (var bid = 1; bid <= bin_counter; bid++) {
            // don't remove nodes from current bin
            if (bid == bin_id)
                continue;

            var pos = SELECTED[bid].indexOf(id_to_node_map[p.child_nodes[i]].label);
            if (pos > -1) {
                SELECTED[bid].splice(pos, 1);

                if (bins_to_update.indexOf(bid) == -1)
                    bins_to_update.push(bid);
            }
        }
    }

    redrawBins();
    updateBinsWindow(bins_to_update);
}

function lineContextMenuHandler(event) {
    if (event.preventDefault) event.preventDefault();
    var bin_id = getBinId();

    if (event.target.id.indexOf('path_') > -1)
    {
        context_menu_target_id = getNodeFromEvent(event).id;

        $('#control_contextmenu').show();

        if (bin_id > 0)
        {
            var pos = SELECTED[bin_id].indexOf(id_to_node_map[parseInt(context_menu_target_id)].label);

            if (pos == -1) {
                $('#control_contextmenu #select').show();
                $('#control_contextmenu #remove').hide();
            }
            else
            {
                $('#control_contextmenu #select').hide();
                $('#control_contextmenu #remove').show();
            }
        }
        else
        {
            $('#control_contextmenu #select').hide();
            $('#control_contextmenu #remove').hide();
        }

        $('#control_contextmenu').offset({left:event.pageX-2,top:event.pageY-2});
        return false;
    }

    var p = getNodeFromEvent(event);

    if (p.id == 0)
        return; // skip root

    if (bin_id === 'undefined')
        return;

    var bins_to_update = [];
    for (var i = 0; i < p.child_nodes.length; i++) {
        // remove nodes from all bins
        for (var bin_id = 1; bin_id <= bin_counter; bin_id++) {
            var pos = SELECTED[bin_id].indexOf(id_to_node_map[p.child_nodes[i]].label);
            if (pos > -1) {
                SELECTED[bin_id].splice(pos, 1);

                if (bins_to_update.indexOf(bin_id) == -1)
                    bins_to_update.push(bin_id);
            }
        }
    }
    redrawBins();
    updateBinsWindow(bins_to_update);
    lineMouseLeaveHandler(event);
    return false;
}

function lineMouseEnterHandler(event) {
    if (drawing_zoom)
        return;

    var p = getNodeFromEvent(event);

    $('#path_hover').remove();

    if (typeof p === 'undefined' || p.id == 0)
        return; // skip root

    var bin_id = getBinId();

    if (bin_id === 'undefined')
        return;

    var bin_color = document.getElementById('bin_color_' + bin_id).getAttribute('color');

    var p1 = p;
    while (p1.child) {
        p1 = p1.child;
    }

    var p2 = p;

    while (p2.child) {
        p2 = p2.child.GetRightMostSibling();
    }

    if (tree_type == 'circlephylogram')
    {
        drawPie('tree_bin',
            'hover',
            p1.angle - angle_per_leaf / 2,
            p2.angle + angle_per_leaf / 2,
            distance(p.backarc, {
                'x': 0,
                'y': 0
            }),
            total_radius,
            (p2.angle - p1.angle + angle_per_leaf > Math.PI) ? 1 : 0,
            bin_color,
            0.3,
            false);
    }
    else
    {  
        drawPhylogramRectangle('tree_bin',
            'hover',
            p.ancestor.xy.x,
            (p1.xy.y + p2.xy.y) / 2,
            p2.xy.y - p1.xy.y + height_per_leaf,
            total_radius - p.ancestor.xy.x,
            bin_color,
            0.3,
            false);
   }

    for (var index = 0; index < p.child_nodes.length; index++) {
        var _line = document.getElementById('line' + p.child_nodes[index]);
        if (_line) {
            _line.style['stroke-width'] = '3';
            _line.style['stroke'] = bin_color;       
        }

        var _arc = document.getElementById('arc' + p.child_nodes[index]);
        if (_arc) {
            _arc.style['stroke-width'] = '3';
            _arc.style['stroke'] = bin_color;
        }
    }
}

function lineMouseLeaveHandler(event) {
    if (drawing_zoom)
        return;

    var p = getNodeFromEvent(event);

    $('#path_hover').remove();

    var bin_id = getBinId();

    if (bin_id === 'undefined') {
        document.focus();
        return;
    }

    if (!p)
        return;

    for (var index = 0; index < p.child_nodes.length; index++) {
        var _line = document.getElementById('line' + p.child_nodes[index]);
        if (_line) {
            _line.style['stroke-width'] = '1';       
        }

        var _arc = document.getElementById('arc' + p.child_nodes[index]);
        if (_arc) {
            _arc.style['stroke-width'] = '1';
        }
    }

    var node_stack = [];
    for (var bin_id = 1; bin_id <= bin_counter; bin_id++) {
        var color_picker = document.getElementById('bin_color_' + bin_id);

        if (!color_picker)
            continue;

        var bin_color = color_picker.getAttribute('color');

        for (var i = 0; i < SELECTED[bin_id].length; i++) {
            node_stack.push(label_to_node_map[SELECTED[bin_id][i]].id);

            var _line = document.getElementById('line' + label_to_node_map[SELECTED[bin_id][i]].id);
            if (_line) {
                _line.style['stroke-width'] = '2';
                _line.style['stroke'] = bin_color;       
            }

            var _arc = document.getElementById('arc' + label_to_node_map[SELECTED[bin_id][i]].id);
            if (_arc) {
                _arc.style['stroke-width'] = '2';
                _arc.style['stroke'] = bin_color;
            }
        }
    }

    for (var i = 0; i < p.child_nodes.length; i++) {
        if (node_stack.indexOf(p.child_nodes[i]) > -1)
            continue;

        var _line = document.getElementById('line' + p.child_nodes[i]);
        if (_line) {
            _line.style['stroke'] = LINE_COLOR;       
        }

        var _arc = document.getElementById('arc' + p.child_nodes[i]);
        if (_arc) {
            _arc.style['stroke'] = LINE_COLOR;
        }
    }
}

function mouseMoveHandler(event) {
    if (drawing_zoom)
        return;
    
    if (event.target.id == 'path_samples')
    {   
        // samples tooltip
        var sample_name = event.target.getAttribute('sample-name');
        var layer_name_hover = event.target.getAttribute('layer-name');

        var message = "";
        for (var i=0; i < last_settings['samples-layer-order'].length; i++)
        {
            var layer_name = last_settings['samples-layer-order'][i];

            if (layer_name == layer_name_hover)
            {
                message += '<tr style="background-color: rgb(232, 202, 207);"><td>' + layer_name + '</td><td>' + samples_information_dict[sample_name][layer_name] + '</td></tr>';
            }
            else
            {
                message += '<tr><td>' + layer_name + '</td><td>' + samples_information_dict[sample_name][layer_name] + '</td></tr>';
            }
        }

        $('#tooltip_content').html(message);
        return;
    }

    var p = getNodeFromEvent(event);

    if (event.target.id && event.target.id == 'path_event')
        lineMouseEnterHandler(event);

    if (!p)
        return;

    var layer_id_exp = event.target.parentNode.id.match(/\d+/);
    if (!layer_id_exp)
        return;
    var layer_id = layer_id_exp[0];

    var tooltip_arr = layerdata_title[id_to_node_map[p.id].label].slice(0);
    
    var message = "";
    for (var i=0; i < tooltip_arr.length; i++)
    {
        if (i == layer_id)
        {
            message += '<tr style="background-color: rgb(232, 202, 207);">' + tooltip_arr[i] + '</tr>';
        }
        else
        {
            message += '<tr>' + tooltip_arr[i] + '</tr>';
        }
    }

    var belongs = "n/a";
    var stop = false;
    var bin_color = '#FFFFFF';

    for (var bin_id = 1; !stop && bin_id <= bin_counter; bin_id++) 
    {
        for (var i = 0; !stop && i < SELECTED[bin_id].length; i++) {
            if (SELECTED[bin_id][i] == p.label) {
                belongs = $('#bin_name_' + bin_id).val();
                bin_color = $('#bin_color_'+ bin_id).attr('color');
                stop = true; // break nested loop
                break;
            }
        }
    }

    var tr_bin = '<tr><td class="tk">bin</td><td class="tv"><div class="colorpicker" style="margin-right: 5px; display: inline-block; background-color:' + bin_color + '"></div>' + belongs + '</td></tr>'

    $('#tooltip_content').html(message + tr_bin);
}


function menu_callback(action) {
    var contig_name = id_to_node_map[context_menu_target_id].label;

    switch (action) {

        case 'select':
            var fake_event = {'target': {'id': '#line' + context_menu_target_id}};
            lineClickHandler(fake_event);
            break;

        case 'remove':
            var fake_event = {'target': {'id': '#line' + context_menu_target_id}};
            lineContextMenuHandler(fake_event);
            break;

        case 'content':
            $.ajax({
                type: 'GET',
                cache: false,
                url: '/data/contig/' + contig_name + '?timestamp=' + new Date().getTime(),
                success: function(data) {
                    $('#splitSequence').val('>' + contig_name + '\n' + data);
                    $('#modSplitSequence').modal('show');
                }
            });
            break;

        case 'blastn_nr': fire_up_ncbi_blast(contig_name, 'blastn', 'nr'); break;
        case 'blastx_nr': fire_up_ncbi_blast(contig_name, 'blastx', 'nr'); break;
        case 'blastn_refseq_genomic': fire_up_ncbi_blast(contig_name, 'blastn', 'refseq_genomic'); break;
        case 'blastx_refseq_protein': fire_up_ncbi_blast(contig_name, 'blastx', 'refseq_genomic'); break;
        
        case 'inspect':
            $.ajax({
                type: 'POST',
                cache: false,
                async: false,
                url: "/data/charts/set_state?timestamp=" + new Date().getTime(), 
                data: {state: JSON.stringify(serializeSettings(true), null, 4)},
                success: function() {
                    window.open('charts.html?contig=' + contig_name, '_blank');
                }
            });
            break;
    }
}

// globals related single background
var rect_left;
var rect_width;

var origin_x;
var origin_y;

function updateSingleBackgroundGlobals()
{
    if (tree_type == 'phylogram')
    {
        var path_event = document.getElementById('path_event');
        var rect = path_event.getBoundingClientRect();

        rect_left = rect.left;
        rect_width = rect.width;
    }
    else // circlephylogram
    {
        var root = document.getElementById('line0');
        var rect = root.getBoundingClientRect();

        var angle = id_to_node_map[0].angle;

        var halfPI = Math.PI / 2;

        if (angle < halfPI)
        {
            origin_x = rect.left;
            origin_y = rect.top;
        }
        else if (angle < 2 * halfPI)
        {
            origin_x = rect.left + rect.width;
            origin_y = rect.top;
        }
        else if (angle < 3 * halfPI)
        {
            origin_x = rect.left + rect.width;
            origin_y = rect.top + rect.height;
        }
        else // 4 * halfPI
        {
            origin_x = rect.left;
            origin_y = rect.top + rect.height;
        }
    }
}

function getNodeFromEvent(event)
{
    if (event.target.id == 'path_event')
    {
        if (tree_type == 'phylogram')
        {
            return order_to_node_map[leaf_count - parseInt((event.clientX - rect_left) / (rect_width / leaf_count)) - 1];
        }
        else
        {
            var _y = event.clientY - origin_y;
            var _x = event.clientX - origin_x;

            var angle = Math.atan2(_y, _x) - angle_per_leaf / 2;
            if (angle < 0)
                angle = 2 * Math.PI + angle;

            var order = Math.ceil((angle - Math.toRadians(last_settings['angle-min']) - (angle_per_leaf / 2)) / angle_per_leaf);
            
            if (order < 1 || order > leaf_count)
                order = 0;

            return order_to_node_map[order]
        }
    }
    else
    {
        var id = event.target.id.match(/\d+/);

        if (id)
            return id_to_node_map[id[0]];
    }
}
