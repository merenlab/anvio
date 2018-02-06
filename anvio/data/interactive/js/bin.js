/**
 * Draw bins, bin labels stuff.
 *
 *  Author: Özcan Esen <ozcanesen@gmail.com>
 *  Credits: A. Murat Eren
 *  Copyright 2017, The anvio Project
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


function redrawBins()
{
    // check if tree parsed, if not there is nothing to redraw.
    if ($.isEmptyObject(label_to_node_map)) 
        return;

    var leaf_list = Array.apply(null, new Array(leaf_count+1)).map(Number.prototype.valueOf,0);

    // put bin numbers of selected leaves to leaf list
    // maybe we should write directly into leaf_list in mouse events, instead of generate it everytime.
    for (var bin_id = 1; bin_id <= bin_counter; bin_id++) {
        for (var j = SELECTED[bin_id].length - 1; j >= 0; j--) {
            var node = label_to_node_map[SELECTED[bin_id][j]];
            if (typeof node === 'undefined')
            {
                SELECTED[bin_id].splice(j, 1);
                continue;
            }

            if (node.IsLeaf() && !node.collapsed) {
                leaf_list[node.order] = bin_id;
            }
        }
    }

    // cluster bins and put them into bins_to_draw array with (start, end, bin_id);
    var prev_value = leaf_list[0];
    var prev_start = 0;

    var bins_to_draw = new Array();

    for (var i=1; i < leaf_list.length; i++)
    {
        if (prev_value != leaf_list[i])
        {
            if (prev_value != 0) {
                bins_to_draw.push(new Array(prev_start, i - 1, prev_value)); // start, end, bin_id;
            }

            prev_start = i;
        }
        prev_value = leaf_list[i];
    }

    // remove exist bin drawings
    var bin = document.getElementById('bin');
    while (bin.hasChildNodes()) {
        bin.removeChild(bin.lastChild);
    }

    // draw new bins
    var show_grid = $('#show_grid_for_bins')[0].checked;
    var grid_color = document.getElementById('grid_color').getAttribute('color');
    var grid_width = $('#grid_width').val();
    var show_bin_labels = $('#show_bin_labels')[0].checked;
    var bin_labels_font_size = parseFloat($('#bin_labels_font_size').val());
    var autorotate_bin_labels = $('#autorotate_bin_labels')[0].checked;
    var bin_labels_angle = $('#bin_labels_angle').val();
    
    var outer_ring_size = parseFloat($('#outer-ring-height').val());
    var outer_ring_margin = parseFloat($('#outer-ring-margin').val());

    for (var i=0; i < bins_to_draw.length; i++) {
        var start = order_to_node_map[bins_to_draw[i][0]];
        var end = order_to_node_map[bins_to_draw[i][1]];

        var color = document.getElementById('bin_color_' + bins_to_draw[i][2]).getAttribute('color');

        if (tree_type == 'circlephylogram')
        {

            drawPie('bin',
                'bin_outer_' + i,
                start.angle - start.size / 2,
                end.angle + end.size / 2,
                total_radius + outer_ring_margin,
                total_radius + outer_ring_margin + outer_ring_size,
                (end.angle - start.angle + (start.size / 2) + (end.size / 2) > Math.PI) ? 1 : 0,
                color,
                1,
                false);

            var align = 'left';
            var angle = (end.angle + end.size / 2 + start.angle - start.size / 2) / 2;
            var new_angle = angle * 180.0 / Math.PI;
            if ((angle > Math.PI / 2.0) && (angle < 1.5 * Math.PI)) {
                align = 'right';
                new_angle += 180.0;
            }

            if (show_bin_labels)
            {
                // so much trigonometry, sorry :(
                var bin_label_radius = total_radius + outer_ring_margin * 1.5 + outer_ring_size * (highlighted_splits.length > 0 ? 2 : 1);
                var bin_label_angle = (end.angle + end.size / 2 + start.angle - start.size / 2) / 2;

                var bin_label_px = bin_label_radius * Math.cos(bin_label_angle);
                var bin_label_py = bin_label_radius * Math.sin(bin_label_angle);

                bin_label_px = bin_label_px - Math.cos(Math.PI / 2 + bin_label_angle) * (bin_labels_font_size / 3) * (align == 'right' ? 1 : -1);
                bin_label_py = bin_label_py - Math.sin(Math.PI / 2 + bin_label_angle) * (bin_labels_font_size / 3) * (align == 'right' ? 1 : -1);

                drawRotatedText(
                    'bin',
                    {
                        'x': bin_label_px, 
                        'y': bin_label_py, 
                    },
                    $('#bin_name_' + bins_to_draw[i][2]).val().replace("_", " "),
                    (autorotate_bin_labels) ? new_angle : bin_labels_angle,
                    align,
                    bin_labels_font_size + "px",
                    "HelveticaNeue-CondensedBold, Helvetica Neue, Helvetica, sans-serif",
                    color,
                    0,
                    'baseline'
                    );

            }

            var pie = drawPie('bin',
                'bin_background_' + i,
                start.angle - start.size / 2,
                end.angle + end.size / 2,
                beginning_of_layers,
                (show_grid) ? total_radius + outer_ring_margin + outer_ring_size : total_radius,
                (Math.abs(end.angle - start.angle) + start.size / 2 + end.size / 2 > Math.PI) ? 1 : 0,
                color,
                (show_grid) ? 0 : 0.1,
                false);

            if (show_grid) {
                pie.setAttribute('vector-effect', 'non-scaling-stroke');
                pie.setAttribute('stroke-opacity', '1');
                pie.setAttribute('stroke-width', grid_width);
                pie.setAttribute('stroke', grid_color);
            }


        }
        else
        {
            var height = end.xy['y'] + end.size / 2 - start.xy['y'] + start.size / 2;

            drawPhylogramRectangle('bin',
                'bin_outer_' + i,
                total_radius + outer_ring_margin,
                start.xy['y'] - start.size / 2 + height / 2,
                height,
                outer_ring_size,
                color,
                1,
                false);

            if (show_bin_labels)
            {
                drawRotatedText(
                    'bin',
                    {
                        'y':  (start.xy.y - start.size / 2 + end.xy.y + end.size / 2) / 2 + (bin_labels_font_size / 3), 
                        'x': (total_radius + outer_ring_margin * 1.5 + outer_ring_size * (highlighted_splits.length > 0 ? 2 : 1)), 
                    },
                    $('#bin_name_' + bins_to_draw[i][2]).val().replace("_", " "),
                    (autorotate_bin_labels) ? 0 : bin_labels_angle,
                    'left',
                    bin_labels_font_size + "px",
                    "HelveticaNeue-CondensedBold, Helvetica Neue, Helvetica, sans-serif",
                    color,
                    0,
                    'baseline'
                    );

            }

            var rect = drawPhylogramRectangle('bin',
                'bin_background_' + i,
                beginning_of_layers,
                start.xy['y'] - start.size / 2 + height / 2,
                height,
                (show_grid) ? total_radius + outer_ring_margin + outer_ring_size - beginning_of_layers : total_radius - beginning_of_layers,
                color,
                (show_grid) ? 0 : 0.1,
                false);

            if (show_grid) {
                rect.setAttribute('vector-effect', 'non-scaling-stroke');
                rect.setAttribute('stroke-opacity', '1');
                rect.setAttribute('stroke-width', grid_width);
                rect.setAttribute('stroke', grid_color);
            }
        }
    }


    // draw higlighted splits
    for (var i=0; i < highlighted_splits.length; i++) {
        // TO DO: more performance
        var start = label_to_node_map[highlighted_splits[i]];
        var end = start;

        var color = document.getElementById('picker_highlight').getAttribute('color');

        if (tree_type == 'circlephylogram')
        {
            drawPie('bin',
                'bin_outer_' + 1,
                start.angle - start.size / 2,
                end.angle + end.size / 2,
                total_radius + outer_ring_margin + outer_ring_size,
                total_radius + outer_ring_margin + outer_ring_size * 2,
                (end.angle - start.angle + start.size / 2 + end.size / 2 > Math.PI) ? 1 : 0,
                color,
                1,
                false);     
        }
        else
        {
            var height = end.xy['y'] + end.size / 2 - start.xy['y'] + start.size / 2;
            
            drawPhylogramRectangle('bin',
                'bin_outer_' + 1,
                total_radius + outer_ring_margin + outer_ring_size,
                start.xy['y'] - start.size / 2 + height / 2,
                height,
                outer_ring_size,
                color,
                1,
                false);
        }
    }

    try{
        var fake_event = {'target': {'id': '#line' + order_to_node_map[0].id}};
        lineMouseLeaveHandler(fake_event);
    }catch(err){
        console.log("Triggering mouseLeaveHandler failed.");
        console.log(err);
    }
}


function rebuildIntersections()
{
    for (var bin_id = 1; bin_id <= bin_counter; bin_id++) {
        // try to make new intersections

        var selected_set = new Set(SELECTED[bin_id]);

        var next_iteration = [].concat(SELECTED[bin_id]);
        var inserted;
        do {
            inserted = 0;
            var nodes = [].concat(next_iteration);
            var length = nodes.length;
            for (var cursor = 0; cursor < length; cursor++)
            {
                if (!label_to_node_map.hasOwnProperty(nodes[cursor])) {
                    continue;
                }

                var node = label_to_node_map[nodes[cursor]];
                var parent = node.ancestor;

                if (parent == null || parent.ancestor == null) 
                {
                    // skip root
                    continue;
                }

                if (selected_set.has(parent.label))
                {
                    // parent already in selected list
                    continue;
                }

                if (node.sibling != null && selected_set.has(node.sibling.label))
                {
                    selected_set.add(parent.label);
                    next_iteration.push(parent.label);
                    next_iteration.splice(next_iteration.indexOf(node.label), 1);
                    next_iteration.splice(next_iteration.indexOf(node.sibling.label), 1);
                    inserted++;
                }
            }

        } while (inserted > 0)

        SELECTED[bin_id] = Array.from(selected_set);
    }
}