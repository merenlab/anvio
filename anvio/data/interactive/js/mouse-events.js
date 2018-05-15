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

        let menu = new ContextMenu({'container': document.body,
                                    'event': event,
                                    'node': p,
                                    'layer': null,
                                    'isSample': true});

        menu.Show();
    } else {
        let p = getNodeFromEvent(event);
        let layer_id = event.target.parentNode.id.match(/\d+/);

        if (p.IsLeaf() || IsCtrlPressed(event)) {
            let menu = new ContextMenu({'container': document.body,
                                        'event': event,
                                        'node': p,
                                        'layer': layer_id,
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

        write_mouse_table(message, "Layers", layer_pos);
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
    var target_node = drawer.tree.nodes[p.id];

    if (target_node.collapsed) {
        $('#tooltip_content').html("Collapsed branch");
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

    write_mouse_table(message, target_node.label, layer_id);
}


function write_mouse_table(content, item_name, layer_id) {
    $('#cell_item_name').html(item_name);
    $('#tooltip_content').html(content);

    if ($('#tooltip_content').height() + 300 > $(window).height()) {
        $('#mouse_hover_scroll').css('top', Math.min(0, ($(window).height()-300) / 2 + -1 * $('#tooltip_content tr').eq(layer_id).position()['top']));
    } else {
        $('#mouse_hover_scroll').css('top', 0);
    } 
} 


function menu_callback(action, param) {
    var item_name = drawer.tree.nodes[context_menu_target_id].label;
    var target = (mode == 'gene') ? 'gene' : 'contig';
    var new_tree;

    switch (action) {
        case 'collapse':
            new_tree = new Tree();
            new_tree.Parse(clusteringData.trim(), false);
            new_tree.FindNode(item_name).collapsed = true;
            clusteringData = new_tree.Serialize();
            $('#tree_modified_warning').show();
            drawTree();
            break;

        case 'expand':
            new_tree = new Tree();
            new_tree.Parse(clusteringData.trim(), false);
            new_tree.FindNode(item_name).collapsed = false;
            clusteringData = new_tree.Serialize();
            $('#tree_modified_warning').show();
            drawTree();
            break;

        case 'rotate':
            new_tree = new Tree();
            new_tree.Parse(clusteringData.trim(), false);
            new_tree.FindNode(item_name).Rotate();
            clusteringData = new_tree.Serialize();
            $('#tree_modified_warning').show();
            drawTree();
            break;

        case 'reroot':
            $.ajax({
                type: 'POST',
                cache: false,
                url: '/data/reroot_tree',
                data: {
                    'newick': clusteringData,
                    'leaves': item_name  
                },
                success: function(data) {
                    clusteringData = data['newick'];
                    $('#tree_modified_warning').show();
                    drawTree();
                }
            });
            break;

        case 'select':
            var fake_event = {'target': {'id': '#line' + context_menu_target_id}};
            lineClickHandler(fake_event);
            break;

        case 'remove':
            bins.RemoveNode(drawer.tree.nodes[context_menu_target_id]);
            break;

        case 'select_layer':
            $('#tbody_layers tr:nth-child(' + context_menu_layer_id + ') input:checkbox').prop('checked', true);
            break;
        case 'unselect_layer':
            $('#tbody_layers tr:nth-child(' + context_menu_layer_id + ') input:checkbox').prop('checked', false);
            break;

        case 'get_gene_sequence':
            $.ajax({
                type: 'GET',
                cache: false,
                url: '/data/gene/' + item_name,
                success: function(data) {
                    $('#modSplitSequence .modal-title').html('Gene Sequence');
                    $('#splitSequence').val('>' + data['header'] + '\n' + data['sequence']);
                    $('#modSplitSequence').modal('show');
                }
            });
            break;

        case 'get_split_sequence':
            $.ajax({
                type: 'GET',
                cache: false,
                url: '/data/contig/' + item_name,
                success: function(data) {
                    $('#modSplitSequence .modal-title').html('Split Sequence');
                    $('#splitSequence').val('>' + data['header'] + '\n' + data['sequence']);
                    $('#modSplitSequence').modal('show');
                }
            });
            break;

        case 'blastn_nr': get_sequence_and_blast(item_name, 'blastn', 'nr', target); break;
        case 'blastx_nr': get_sequence_and_blast(item_name, 'blastx', 'nr', target); break;
        case 'blastn_refseq_genomic': get_sequence_and_blast(item_name, 'blastn', 'refseq_genomic', target); break;
        case 'blastx_refseq_protein': get_sequence_and_blast(item_name, 'blastx', 'refseq_genomic', target); break;

        // collection mode-specific:
        case 'refine_bin': toastr.error('Refine function from the interface is not currently implemented :/ ' +
                                        'Please use `anvi-refine` program for "' + item_name  +'"'); break;
        case 'get_hmm_sequence':
            $.ajax({
                type: 'GET',
                cache: false,
                url: '/data/hmm/' + item_name + '/' + param,
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
            localStorage.state = JSON.stringify(serializeSettings(true), null, 4);
            window.open(generate_inspect_link('inspect', item_name), '_blank');
            break;

        case 'inspect_gene':
            localStorage.state = JSON.stringify(serializeSettings(true), null, 4);
            window.open(generate_inspect_link('inspect_gene', item_name), '_blank');
            break;

        case 'inspect_context':
            localStorage.state = JSON.stringify(serializeSettings(true), null, 4);
            window.open(generate_inspect_link('inspect_context', item_name), '_blank');
            break;

        case 'inspect_gene_cluster':
            localStorage.state = JSON.stringify(serializeSettings(true), null, 4);
            window.open(generate_inspect_link('geneclusters', item_name), '_blank');
            break;

        case 'get_AA_sequences_for_gene_cluster':
            $.ajax({
                type: 'GET',
                cache: false,
                url: '/data/get_AA_sequences_for_gene_cluster/' + item_name,
                success: function(data) {
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
    if (event.target.id == 'path_event')
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
