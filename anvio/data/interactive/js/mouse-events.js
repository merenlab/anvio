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

    if (event.target.parentNode && event.target.parentNode.id == 'samples_tree')
    {
        var id = event.target.id.match(/\d+/);
        var node = samples_id_to_node_map[id];

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

        if (mode == "collection") {
            $('#collection_mode_right_click_menu').show();
            $('#collection_mode_right_click_menu').offset({left:event.pageX-2,top:event.pageY-2});
        } else if (mode == "pan"){
            $('#pan_mode_right_click_menu').show();
            $('#pan_mode_right_click_menu').offset({left:event.pageX-2,top:event.pageY-2});
        } else {
            $('#default_right_click_menu').show();
            $('#default_right_click_menu').offset({left:event.pageX-2,top:event.pageY-2});
        }

        context_menu_layer_id = event.target.parentNode.id.match(/\d+/);
        var is_layer_selected = $('#tbody_layers tr:nth-child(' + context_menu_layer_id + ') input:checkbox').is(':checked')

        if (is_layer_selected) {
            $('#default_right_click_menu #unselect_layer').show();
            $('#default_right_click_menu #select_layer').hide();
            $('#collection_mode_right_click_menu #unselect_layer').show();
            $('#collection_mode_right_click_menu #select_layer').show();
            $('#pan_mode_right_click_menu #unselect_layer').show();
            $('#pan_mode_right_click_menu #select_layer').show();
        } else {
            $('#default_right_click_menu #select_layer').show();
            $('#default_right_click_menu #unselect_layer').hide();
            $('#collection_mode_right_click_menu #select_layer').show();
            $('#collection_mode_right_click_menu #unselect_layer').show();
            $('#pan_mode_right_click_menu #select_layer').show();
            $('#pan_mode_right_click_menu #unselect_layer').show();
        }

        if (bin_id > 0) {
            var pos = SELECTED[bin_id].indexOf(id_to_node_map[parseInt(context_menu_target_id)].label);

            if (pos == -1) {
                $('#default_right_click_menu #select').show();
                $('#default_right_click_menu #remove').hide();
                $('#collection_mode_right_click_menu #select').show();
                $('#collection_mode_right_click_menu #remove').hide();
                $('#pan_mode_right_click_menu #select').show();
                $('#pan_mode_right_click_menu #remove').hide();
            }
            else {
                $('#default_right_click_menu #select').hide();
                $('#default_right_click_menu #remove').show();
                $('#collection_mode_right_click_menu #select').hide();
                $('#collection_mode_right_click_menu #remove').show();
                $('#pan_mode_right_click_menu #select').hide();
                $('#pan_mode_right_click_menu #remove').show();
            }
        }
        else {
            $('#default_right_click_menu #select').hide();
            $('#default_right_click_menu #remove').hide();
            $('#collection_mode_right_click_menu #select').hide();
            $('#collection_mode_right_click_menu #remove').hide();
            $('#pan_mode_right_click_menu #select').hide();
            $('#pan_mode_right_click_menu #remove').hide();
        }

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
        var _x = (p.ancestor) ? p.ancestor.xy.x : p.xy.x;
        drawPhylogramRectangle('tree_bin',
            'hover',
             _x,
            (p1.xy.y + p2.xy.y) / 2,
            p2.xy.y - p1.xy.y + height_per_leaf,
            total_radius - _x,
            bin_color,
            0.3,
            false);
    }

    if (p.IsLeaf())
        return;

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
        var layer_name_hover = event.target.getAttribute('layer-name');

        var message = "";
        var layer_pos = 0;
        for (var i=0; i < last_settings['samples-layer-order'].length; i++)
        {
            var layer_name = last_settings['samples-layer-order'][i];
            var pretty_name = (layer_name.indexOf('!') > -1) ? layer_name.split('!')[0] : layer_name;

            if (layer_name == layer_name_hover)
            {
                message += '<tr style="background-color: rgb(232, 202, 207);"><td>' + pretty_name + '</td><td>' + samples_information_dict[sample_name][layer_name] + '</td></tr>';
                layer_pos = i;
            }
            else
            {
                message += '<tr><td>' + pretty_name + '</td><td>' + samples_information_dict[sample_name][layer_name] + '</td></tr>';
            }
        }

        $('#tooltip_content').html(message);
        if ($('#tooltip_content').height() + 300 > $(window).height()) {
            $('#mouse_hover_scroll').css('top', ($(window).height()-300) / 2 + -1 * $('#tooltip_content tr').eq(layer_pos).position()['top']);
        } else {
            $('#mouse_hover_scroll').css('top', 0);
        }
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
    if ($('#tooltip_content').height() + 300 > $(window).height()) {
        $('#mouse_hover_scroll').css('top', ($(window).height()-300) / 2 + -1 * $('#tooltip_content tr').eq(layer_id).position()['top']);
    } else {
        $('#mouse_hover_scroll').css('top', 0);
    }
}


function menu_callback(action, param) {
    var item_name = id_to_node_map[context_menu_target_id].label;

    switch (action) {

        case 'select':
            var fake_event = {'target': {'id': '#line' + context_menu_target_id}};
            lineClickHandler(fake_event);
            break;

        case 'remove':
            var fake_event = {'target': {'id': '#line' + context_menu_target_id}};
            lineContextMenuHandler(fake_event);
            break;

        case 'select_layer':
            $('#tbody_layers tr:nth-child(' + context_menu_layer_id + ') input:checkbox').prop('checked', true);
            break;
        case 'unselect_layer':
            $('#tbody_layers tr:nth-child(' + context_menu_layer_id + ') input:checkbox').prop('checked', false);
            break;

        case 'get_split_sequence':
            $.ajax({
                type: 'GET',
                cache: false,
                url: '/data/contig/' + item_name + '?timestamp=' + new Date().getTime(),
                success: function(data) {
                    $('#splitSequence').val('>' + data['header'] + '\n' + data['sequence']);
                    $('#modSplitSequence').modal('show');
                }
            });
            break;

        case 'blastn_nr': fire_up_ncbi_blast(item_name, 'blastn', 'nr', 'contig'); break;
        case 'blastx_nr': fire_up_ncbi_blast(item_name, 'blastx', 'nr', 'contig'); break;
        case 'blastn_refseq_genomic': fire_up_ncbi_blast(item_name, 'blastn', 'refseq_genomic', 'contig'); break;
        case 'blastx_refseq_protein': fire_up_ncbi_blast(item_name, 'blastx', 'refseq_genomic', 'contig'); break;

        // collection mode-specific:
        case 'refine_bin': toastr.error('Refine function from the interface is not currently implemented :/ ' +
                                        'Please use `anvi-refine` program for "' + item_name  +'"'); break;
        case 'get_hmm_sequence':
            $.ajax({
                type: 'GET',
                cache: false,
                url: '/data/hmm/' + item_name + '/' + param + '?timestamp=' + new Date().getTime(),
                success: function(data) {
                    if ('error' in data){
                        $('#modGenerateSummary').modal('hide');
                        waitingDialog.hide();
                        toastr.error(data['error'], "", { 'timeOut': '0', 'extendedTimeOut': '0' });
                    } else {
                        $('#splitSequence').val('>' + data['header'] + '\n' + data['sequence']);
                        $('#modSplitSequence').modal('show');
                    }
                }
            });
            break;

        case 'inspect_contig':
            $.ajax({
                type: 'POST',
                cache: false,
                async: false,
                url: "/data/charts/set_state?timestamp=" + new Date().getTime(), 
                data: {state: JSON.stringify(serializeSettings(true), null, 4)},
                success: function() {
                    window.open('charts.html?contig=' + item_name, '_blank');
                }
            });
            break;

        case 'inspect_protein_cluster':
            $.ajax({
                type: 'POST',
                cache: false,
                async: false,
                url: "/data/proteinclusters/set_state?timestamp=" + new Date().getTime(), 
                data: {state: JSON.stringify(serializeSettings(true), null, 4)},
                success: function() {
                    window.open('proteinclusters.html?id=' + item_name, '_blank');
                }
            });
            break;

        case 'get_AA_sequences_for_PC':
            $.ajax({
                type: 'GET',
                cache: false,
                url: '/data/get_AA_sequences_for_PC/' + item_name + '?timestamp=' + new Date().getTime(),
                success: function(data) {
                    data = JSON.parse(data);

                    var output = '';

                    for (var key in data)
                        output = output + ">" + key + "\n" + data[key] + "\n";

                    $('#splitSequence').val(output);
                    $('#modSplitSequence').modal('show');
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
        var root = document.getElementById('line_origin');
        var rect = root.getBoundingClientRect();

        origin_x = rect.left;
        origin_y = rect.top;
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
