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


function Bins(prefix, container) {
    this.selections = {}
    this.bin_counter = 0;
    this.prefix = prefix || "Bin_";
    this.higlighted_items = [];
    this.container = container || document.createElement("div");

    this.cache = {};
}


Bins.prototype.NewBin = function(id, binState) {
    if (typeof id === 'undefined')
    {
        var from_state = false;
        var id = this.bin_counter;
        var name = this.prefix + (id + 1);
        var color = randomColor({luminosity: 'dark'});
        var contig_count = 0;
        var contig_length = "N/A";
        var completeness = '---';
        var redundancy = '---';

        this.selections[id] = [];
        this.bin_counter++;
    }
    else
    {
        // we are adding bins from collection
        var from_state = true;
        var name = binState['name'];
        var color = binState['color'];
        var contig_count = 0;
        var contig_length = "N/A";
        var completeness = "---";
        var redundancy = "---";
    }

    var template = '<tr bin-id="{id}">' +
                   '    <td><input type="radio" name="active_bin" value="{id}" checked></td>' +
                   '    <td><div id="bin_color_{id}" class="colorpicker" color="{color}" style="background-color: {color}"></td>' +
                   '    <td data-value="{name}"><input type="text" onChange="redrawBins();" size="21" id="bin_name_{id}" value="{name}"></td>';

    if (mode != 'pan')
    {
        template +='    <td data-value="{count}" class="num-items"><input type="button" value="{count}" title="Click for contig names" onClick="showContigNames({id});"></td> ' +
                   '    <td data-value="{length}" class="length-sum"><span>{length}</span></td>';
    }

    template +=    '    <td data-value="{completeness}"><input id="completeness_{id}" type="button" value="{completeness}" title="Click for completeness table" onClick="showCompleteness({id});"></td> ' +
                   '    <td data-value="{redundancy}"><input id="redundancy_{id}" type="button" value="{redundancy}" title="Click for redundant hits" onClick="showRedundants({id});"></td> ' +
                   '    <td><center><span class="glyphicon glyphicon-trash" aria-hidden="true" alt="Delete this bin" title="Delete this bin" onClick="bins.DeleteBin({id});"></span></center></td>' +
                   '</tr>';

    template = template.replace(new RegExp('{id}', 'g'), id)
                       .replace(new RegExp('{name}', 'g'), name)
                       .replace(new RegExp('{color}', 'g'), color)
                       .replace(new RegExp('{count}', 'g'), contig_count)
                       .replace(new RegExp('{completeness}', 'g'), completeness)
                       .replace(new RegExp('{redundancy}', 'g'), redundancy)
                       .replace(new RegExp('{length}', 'g'), contig_length);

    this.container.insertAdjacentHTML('beforeend', template);

/*    if(!from_state){
        $('#completeness_' + id).attr("disabled", true);
        $('#redundancy_' + id).attr("disabled", true);
    }*/
/*
    $('#bin_color_' + id).colpick({
        layout: 'hex',
        submit: 0,
        colorScheme: 'light',
        onChange: function(hsb, hex, rgb, el, bySetColor) {
            $(el).css('background-color', '#' + hex);
            $(el).attr('color', '#' + hex);

            if (!bySetColor) $(el).val(hex);
        },
        onHide: function() {
            redrawBins();
        }
    }).keyup(function() {
        $(this).colpickSetColor(this.value);
    });*/
}


Bins.prototype.IsEmpty = function() {
    return !this.container.hasChildNodes();
};


Bins.prototype.GetSelectedBinId = function() {
    return this.container.querySelector('input[name=active_bin]:checked').value;
};


Bins.prototype.GetSelectedBinColor = function(id) {
    return this.container.querySelector('#bin_color_' + this.GetSelectedBinId()).getAttribute('color');
};


Bins.prototype.DeleteBin = function(id, show_confirm=true) {
    if (show_confirm && !confirm('Are you sure?')) {
        return;
    }

    this.container.querySelector('#bin_row_' + id).remove();
    this.container.querySelectorAll('input[type=radio]')[0].setAttribute('checked', true);

    for (var i = 0; i < SELECTED[id].length; i++) {
        var node = drawer.tree.nodes[SELECTED[id][i]];

        if (typeof node === 'undefined' || !node.hasOwnProperty('id')) {
            continue;
        }

        var node_id = node.id;

        let line = document.getElementById('line' + child.id);
        if (line) {
            line.style['stroke-width'] = '3';
            line.style['stroke'] = bin_color;       
        }

        let arc = document.getElementById('arc' + child.id);
        if (arc) {
            arc.style['stroke-width'] = '3';
            arc.style['stroke'] = bin_color;
        }
    }

    SELECTED[id] = [];

    if (this.IsEmpty())
    {
        this.NewBin();
    }

    this.RedrawBins();
};


Bins.prototype.DeleteAllBins = function() {
    if (!confirm('Are you sure you want to remove all bins?')) {
        return;
    }
    var bin_ids_to_delete = [];

    $('#tbody_bins tr').each(
        function(index, bin) {
            bin_ids_to_delete.push($(bin).attr('bin-id'));
        }
    );

    bin_ids_to_delete.map(function(bin_id) { 
        deleteBin(bin_id, false);
    });
};


Bins.prototype.AppendBranch = function(p) {
    var bin_id = this.GetSelectedBinId();
    var bin_color = this.GetSelectedBinColor();
    var bins_to_update = [];

    for (const child of p.IterateChildren()) {
        let pos = this.selections[bin_id].indexOf(child.id);
        if (pos == -1) {
            this.selections[bin_id].push(child.id);

            if (bins_to_update.indexOf(bin_id) == -1)
                bins_to_update.push(bin_id);
        }

        // remove nodes from other bins
        for (let bid in this.selections) {
            // don't remove nodes from current bin
            if (bid == bin_id)
                continue;

            let pos = this.selections[bin_id].indexOf(child.id)
            if (pos > -1) {
                this.selections[bin_id].splice(pos, 1);

                if (bins_to_update.indexOf(bid) == -1)
                    bins_to_update.push(bid);
            }
        }

        let line = document.getElementById('line' + child.id);
        if (line) {
            line.style['stroke-width'] = '3';
            line.style['stroke'] = bin_color;       
        }

        let arc = document.getElementById('arc' + child.id);
        if (arc) {
            arc.style['stroke-width'] = '3';
            arc.style['stroke'] = bin_color;
        }
    }

    this.RedrawBins();
    this.UpdateBinsWindow(bins_to_update);
};


Bins.prototype.RemoveBranch = function(p) {
    var bins_to_update = [];
    for (const child of p.IterateChildren()) {
        for (let bin_id in this.selections) {
            let pos = this.selections[bin_id].indexOf(child.id);
            if (pos > -1) {
                this.selections[bin_id].splice(pos, 1);

                if (bins_to_update.indexOf(bin_id) == -1)
                    bins_to_update.push(bin_id);
            }
        }

        let line = document.getElementById('line' + child.id);
        if (line) {
            line.style['stroke-width'] = '1';
            line.style['stroke'] = LINE_COLOR;       
        }

        let arc = document.getElementById('arc' + child.id);
        if (arc) {
            arc.style['stroke-width'] = '1';
            arc.style['stroke'] = LINE_COLOR;
        }
    }

    this.RedrawBins();
    this.UpdateBinsWindow(bins_to_update);
};


Bins.prototype.UpdateBinsWindow = function(bin_list) {
    if (typeof bin_list === 'undefined')
    {
        var bin_list = Object.keys(self.selections);
    }

    for (let i = 0; i < bin_list.length; i++) {
        let bin_id = bin_list[i];

        if (mode == 'pan') {

        } else {
            let num_items = 0;
            let length_sum = 0;

            for (let j = 0; j < this.selections[bin_id].length; j++) {
                let node_id = this.selections[bin_id][j];
                let node = drawer.tree.nodes[node_id];

                if (node.IsLeaf()) {
                    num_items++;
                    length_sum += parseInt(contig_lengths[node.label]);
                }

                let bin_row = this.container.querySelector(`tr[bin-id="${bin_id}"]`);

                bin_row.querySelector('td.num-items').setAttribute('data-value', num_items);
                bin_row.querySelector('td.num-items>input').value = num_items;

                if (isNaN(length_sum)) {
                    bin_row.querySelector('td.length-sum').setAttribute('data-value', 0);
                    bin_row.querySelector('td.length-sum>span').innerHTML = 'n/a';
                } else {
                    bin_row.querySelector('td.length-sum').setAttribute('data-value', length_sum);
                    bin_row.querySelector('td.length-sum>span').innerHTML = readableNumber(length_sum);
                }
            }        
        }
    }

    $('#bin_settings_tab:not(.active) a').css('color', "#ff0000");
};


Bins.prototype.RedrawBins = function() {
    if (!drawer)
        return;

    var leaf_list = [];
    for (var i=0; i < drawer.tree.leaves.length + 1; i++) {
        leaf_list.push(-1);
    }

    for (let bin_id in this.selections) {
        for (let j = this.selections[bin_id].length - 1; j >= 0; j--) {
            let node = drawer.tree.nodes[this.selections[bin_id][j]];

            if (typeof node === 'undefined')
            {
                this.selections[bin_id].splice(j, 1);
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
            if (prev_value != -1) {
                bins_to_draw.push(new Array(prev_start, i - 1, prev_value)); // start, end, bin_id;
            }

            prev_start = i;
        }
        prev_value = leaf_list[i];
    }

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
        var start = drawer.tree.leaves[bins_to_draw[i][0]];
        var end = drawer.tree.leaves[bins_to_draw[i][1]];

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
                var bin_label_radius = total_radius + outer_ring_margin * 1.5 + outer_ring_size * (this.higlighted_items.length > 0 ? 2 : 1);
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
                        'x': (total_radius + outer_ring_margin * 1.5 + outer_ring_size * (this.higlighted_items.length > 0 ? 2 : 1)), 
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


/*    // draw higlighted splits
    for (var i=0; i < highlighted_splits.length; i++) {
        // TO DO: more performance
        var start = drawer.tree.nodes[highlighted_splits[i]];
        var end = start;

        if (!start)
            continue;

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
    }*/
}


function rebuildIntersections()
{
    return;
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
                if (!drawer.tree.nodes.hasOwnProperty(nodes[cursor].id)) {
                    continue;
                }

                var node = drawer.tree.nodes[nodes[cursor].id];
                var parent = node.ancestor;

                if (parent == null || parent.ancestor == null) 
                {
                    // skip root
                    continue;
                }

                if (selected_set.has(parent.id))
                {
                    // parent already in selected list
                    continue;
                }

                if (node.sibling != null && selected_set.has(node.sibling.id))
                {
                    selected_set.add(parent.id);
                    next_iteration.push(parent.id);
                    next_iteration.splice(next_iteration.indexOf(node.id), 1);
                    next_iteration.splice(next_iteration.indexOf(node.sibling.id), 1);
                    inserted++;
                }
            }

        } while (inserted > 0)

        SELECTED[bin_id] = Array.from(selected_set);
    }
}
