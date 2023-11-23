/**
 * Javascript library to visualize additional layer info
 * (previously known as samples db)
 *
 *  Authors: Ozcan Esen
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

var samples_categorical_colors = {};
var samples_categorical_stats = {};
var samples_stack_bar_colors = {};
var samples_stack_bar_stats = {};
var samples_information_dict;
var samples_order_dict;

var samplesClusteringData = {'newick': '', 'basic': null};

function get_newick_leaf_order(newick)
{
    var _order_list = [];

    var t = new Tree();
    t.Parse(newick, false);

    var n = new NodeIterator(t.root);
    var q = n.Begin();

    order_counter = 0;
    while (q != null)
    {
        if (q.IsLeaf()) {
            _order_list.push(q.label);
        }
        q=n.Next();
    }

    return _order_list;
}

$(document).ready(function() {
    $('#samples_order').change(function() {
        $('#btn_redraw_samples').prop('disabled', true);
        $('#btn_redraw_samples_layer').prop('disabled', true);

        if (samples_order_dict.hasOwnProperty(this.value)) {
            samplesClusteringData = samples_order_dict[this.value];
        }

        var samples_new_order;

        // get new sample order
        if (samplesClusteringData['basic'] != null && samplesClusteringData['basic'] != "")
        {
            samples_new_order = samplesClusteringData['basic'].split(',');
        }
        else
        {
            samples_new_order = get_newick_leaf_order(samplesClusteringData['newick']);
        }
        $.map(samples_new_order, $.trim);

        var new_order = [];

        var i=0;
        $('#tbody_layers tr').each(
            function(index, layer) {
                var layer_id = $(layer).find('.input-height')[0].id.replace('height', '');
                var layer_name = getLayerName(layer_id);

                if (samples_new_order.indexOf(layer_name) > -1)
                {
                    new_order.push(samples_new_order[i++]);
                }
                else
                {
                    new_order.push(layer_name);
                }
            }
        );

        for(var i=0; i < new_order.length; i++)
        {
            // sort main layers with new order.
            var layer_id = getLayerId(new_order[i]);
            var detached_row = $('#height' + layer_id).closest('tr').detach();
            $('#tbody_layers').append(detached_row);
        }

        // sort sample layers with new order, but only if they are member of a group starting with "ANI_" or "SourMash_"
        for (let group in samples_information_dict) {
            if (group.startsWith('ANI_') || group.startsWith('SourMash_')) {
                for(var i=new_order.length - 1; i >= 0; i--) {
                    let detached_sample_row = $(`tr[samples-group-name='${group}'][samples-layer-name='${new_order[i]}']`);
                    $('#tbody_samples').append(detached_sample_row);
                }
            }
        }
    });


    $('#group_select_all').click(function() {
        let target_group = $('#group_list_for_select_all').val();

        $('#tbody_samples tr').each((index, tr) => {
            let group = $(tr).attr('samples-group-name');

            if (group == target_group) {
                $(tr).find('.layer_selectors').prop('checked', true);
            }
        });
    });

    $('#group_unselect_all').click(function() {
        let target_group = $('#group_list_for_select_all').val();

        $('#tbody_samples tr').each((index, tr) => {
            let group = $(tr).attr('samples-group-name');

            if (group == target_group) {
                $(tr).find('.layer_selectors').prop('checked', false);
            }
        });
    });
});


function samplesMergeStackbarLayers(data_dict, order_list) {
    let new_dict = {};
    let new_order_list = {};

    for (let group in data_dict) {
        new_dict[group] = {};

        for (let item in data_dict[group]) {

            new_dict[group][item] = {};
            for (let i=0; i < order_list[group].length; i++) {
                let data_key = order_list[group][i];
                if (data_key.indexOf('!') > -1) {
                    let stack_bar_name = data_key.split('!')[0];
                    let stack_bar_layer_name = data_key.split('!')[1];
                    let found_key_in_new_dict = "";

                    for (let data_key_in_new in new_dict[group][item]) {
                        if (data_key_in_new.startsWith(stack_bar_name + '!')) {
                            found_key_in_new_dict = data_key_in_new;
                            break;
                        }
                    }

                    if (found_key_in_new_dict === "") {
                        new_dict[group][item][data_key] = data_dict[group][item][data_key];
                    } else {
                        new_dict[group][item][found_key_in_new_dict + ';' + stack_bar_layer_name] = new_dict[group][item][found_key_in_new_dict] + ';' + data_dict[group][item][data_key];
                        delete new_dict[group][item][found_key_in_new_dict];
                    }
                } else {
                    new_dict[group][item][data_key] = data_dict[group][item][data_key];
                }
            }
        }

        let first_item = Object.keys(new_dict[group])[0];
        new_order_list[group] = Object.keys(new_dict[group][first_item]);
    }

    return {'dict': new_dict, 'default_order': new_order_list};
}

function convert_samples_order_to_array(input_dict) {
    let all_information_layers = [];

    for (let group in input_dict) {
        input_dict[group].forEach((layer_name) => {
            if (samples_information_dict.hasOwnProperty(group)) {
                let first_sample = Object.keys(samples_information_dict[group])[0];
                if (samples_information_dict[group][first_sample].hasOwnProperty(layer_name)) {
                    all_information_layers.push({
                        'group': group,
                        'layer_name': layer_name,
                    });
                }
            }
        });
    }

    return all_information_layers;
}

function is_sample_group_visible(group_name) {
    return $('input:checkbox#group_' + group_name).is(':checked');
}

function toggleSampleGroups() {
    let visible_groups = new Set([]);

    $('#tbody_samples tr').each((index, tr) => {
        let group = $(tr).attr('samples-group-name');

        if (is_sample_group_visible(group)) {
            $(tr).show();
            visible_groups.add(group);
        } else {
            $(tr).hide();
            $(tr).find('.layer_selectors').prop('checked', false);
        }
    });

    $('#group_list_for_select_all').empty();
    visible_groups.forEach((group) => {
        $('#group_list_for_select_all').append(`<option val="${group}">${group}</option>`);
    });
};


function update_samples_layer_min_max(select) {
    let tr = $(select).closest('tr');

    let group = $(tr).attr('samples-group-name');
    let sample_name = $(tr).attr('samples-layer-name');

    let norm = $(select).val();

    let max = null;

    for (let lname in samples_information_dict[group]) {
        if (norm == 'none') {
            if (max === null || parseFloat(samples_information_dict[group][lname][sample_name]) > max) {
                max = parseFloat(samples_information_dict[group][lname][sample_name]);
            }
        }
        else if (norm == 'sqrt') {
            if (max === null || Math.sqrt(parseFloat(samples_information_dict[group][lname][sample_name])) > max) {
                max = Math.sqrt(parseFloat(samples_information_dict[group][lname][sample_name]));
            }
        }
        else if (norm == 'log') {
            if (max === null || log10(parseFloat(samples_information_dict[group][lname][sample_name]) + 1) > max) {
                max = log10(parseFloat(samples_information_dict[group][lname][sample_name]) + 1);
            }
        }
    }

    $(tr).find('.input-min').val('0');
    $(tr).find('.input-max').val(max);
};


function buildSamplesTable(samples_layer_order, samples_layers) {
    var all_information_layers = samples_layer_order;

    for (let group in samples_information_dict) {
        let first_sample = Object.keys(samples_information_dict[group])[0];
        for (let layer_name in samples_information_dict[group][first_sample]) {
            let found = false;

            for (const [index, entry] of all_information_layers.entries()) {
                if (entry['group'] == group && entry['layer_name'] == layer_name) {
                    found = true;
                    break;
                }
            };

            if (!found) {
                all_information_layers.push({
                    'group': group,
                    'layer_name': layer_name,
                });
            }
        }
    }

    if (all_information_layers.length == 0) {
        return;
    }

    $('#tbody_samples').empty();

    for (var i=0; i < all_information_layers.length; i++)
    {
        var layer_name = all_information_layers[i]['layer_name'];
        var group = all_information_layers[i]['group'];

        if (!samples_information_dict.hasOwnProperty(group) || Object.keys(samples_information_dict[group]) < 1) {
            continue;
        }

        var first_sample = Object.keys(samples_information_dict[group])[0];

        if (!samples_information_dict[group][first_sample].hasOwnProperty(layer_name)) {
            continue;
        }

        var pretty_name = getNamedLayerDefaults(layer_name, 'pretty_name', layer_name);
        pretty_name = (pretty_name.indexOf('!') > -1) ? pretty_name.split('!')[0] : pretty_name;

        var short_name = (pretty_name.length > 10) ? pretty_name.slice(0,10) + "..." : pretty_name;

        var hasSettings = false;
        if (typeof(samples_layers) !== 'undefined' &&
            typeof(samples_layers[group]) !== 'undefined' &&
            typeof(samples_layers[group][layer_name]) !== 'undefined') {
                hasSettings = true;
                layer_settings = samples_layers[group][layer_name];
        }

        if (isNumber(samples_information_dict[group][first_sample][layer_name]))
        {
            var data_type = "numeric";

            if (hasSettings)
            {
                var norm         = layer_settings['normalization'];
                var min          = layer_settings['min']['value'];
                var max          = layer_settings['max']['value'];
                var height       = layer_settings['height'];
                var color        = layer_settings['color'];
                var margin       = layer_settings['margin'];
                var color_start  = layer_settings['color-start'];
                var type         = layer_settings['type'];
            }
            else
            {
                var norm         = getNamedLayerDefaults(layer_name, 'norm', 'none', group);
                var min          = getNamedLayerDefaults(layer_name, 'min', null, group);
                var max          = getNamedLayerDefaults(layer_name, 'max', null, group);
                var height       = getNamedLayerDefaults(layer_name, 'height', 500, group);
                var color        = getNamedLayerDefaults(layer_name, 'color', '#919191', group);
                var color_start  = getNamedLayerDefaults(layer_name, 'color-start', '#EFEFEF', group);
                var type         = getNamedLayerDefaults(layer_name, 'type', 'bar', group);
                var margin       = 15;
            }

            var template = '<tr samples-group-name="{group}" samples-layer-name="{name}" data-type="{data-type}">' +
                '<td><img class="drag-icon" src="images/drag.gif" /></td>' +
                '<td title="{name}" class="titles">{short-name}</td>' +
                '<td><div class="colorpicker picker_start" color="{color-start}" style="background-color: {color-start}; {color-start-hide}"></div><div class="colorpicker" color="{color}" style="background-color: {color}"></div></td>' +
                '<td style="width: 50px;">' +
                '    <select style="width: 50px;" class="type type_multiple form-control form-control-sm select-sm" onChange="togglePickerStart(this);">' +
                '        <option value="bar"{option-type-bar}>Bar</option>' +
                '        <option value="intensity"{option-type-intensity}>Intensity</option>' +
                '    </select>' +
                '</td>' +
                '<td>' +
                '    <select onChange="update_samples_layer_min_max(this);" class="normalization type type_multiple form-control form-control-sm col-12 select-sm">' +
                '        <option value="none"{option-none}>none</option>' +
                '        <option value="sqrt"{option-sqrt}>sqrt</option>' +
                '        <option value="log"{option-log}>log</option>' +
                '    </select>' +
                '</td>' +
                '<td><input class="input-height form-control form-control-sm" type="text" size="3" value="{height}"></input></td>' +
                '<td><input class="input-margin form-control form-control-sm" type="text" size="3" value="{margin}"></input></td>' +
                '<td><input class="input-min form-control form-control-sm" type="text" size="4" value="{min}"></input></td>' +
                '<td><input class="input-max form-control form-control-sm" type="text" size="4" value="{max}"></input></td>' +
                '<td><input type="checkbox" class="layer_selectors"></input></td>' +
                '</tr>';

            template = template.replace(new RegExp('{name}', 'g'), layer_name)
                               .replace(new RegExp('{group}', 'g'), group)
                               .replace(new RegExp('{data-type}', 'g'), data_type)
                               .replace(new RegExp('{short-name}', 'g'), short_name)
                               .replace(new RegExp('{option-' + norm + '}', 'g'), ' selected')
                               .replace(new RegExp('{option-([a-z]*)}', 'g'), '')
                               .replace(new RegExp('{option-type-' + type + '}', 'g'), ' selected')
                               .replace(new RegExp('{option-type-([a-z]*)}', 'g'), '')
                               .replace(new RegExp('{color}', 'g'), color)
                               .replace(new RegExp('{color-start}', 'g'), color_start)
                               .replace(new RegExp('{color-start-hide}', 'g'), (type!='intensity') ? '; visibility: hidden;' : '')
                               .replace(new RegExp('{height}', 'g'), height)
                               .replace(new RegExp('{min}', 'g'), min)
                               .replace(new RegExp('{max}', 'g'), max)
                               .replace(new RegExp('{margin}', 'g'), margin);
        }
        else if (layer_name.indexOf(';') > -1)
        {
            var data_type = "stack-bar";

            if (hasSettings)
            {
                var norm   = layer_settings['normalization'];
                var height = layer_settings['height'];
                var margin = layer_settings['margin'];
            }
            else
            {
                var norm   = getNamedLayerDefaults(layer_name, 'norm', 'none', group);
                var height = getNamedLayerDefaults(layer_name, 'height', 500, group);
                var margin = 15;
            }

            var template = '<tr samples-group-name="{group}" samples-layer-name="{name}" data-type="{data-type}">' +
                '<td><img class="drag-icon" src="images/drag.gif" /></td>' +
                '<td title="{pretty-name}" class="titles">{short-name}</td>' +
                '<td style="width: 50px;">n/a</td>' +
                '<td style="width: 50px;">n/a</td>' +
                '<td>' +
                '    <select class="normalization type type_multiple form-control form-control-sm col-12 select-sm">' +
                '        <option value="none"{option-none}>none</option>' +
                '        <option value="sqrt"{option-sqrt}>sqrt</option>' +
                '        <option value="log"{option-log}>log</option>' +
                '    </select>' +
                '</td>' +
                '<td><input class="input-height form-control form-control-sm" type="text" size="3" value="{height}"></input></td>' +
                '<td><input class="input-margin form-control form-control-sm" type="text" size="3" value="{margin}"></input></td>' +
                '<td style="width: 50px;">n/a</td>' +
                '<td style="width: 50px;">n/a</td>' +
                '<td><input type="checkbox" class="layer_selectors"></input></td>' +
                '</tr>';

            template = template.replace(new RegExp('{name}', 'g'), layer_name)
                               .replace(new RegExp('{group}', 'g'), group)
                               .replace(new RegExp('{data-type}', 'g'), data_type)
                               .replace(new RegExp('{short-name}', 'g'), short_name)
                               .replace(new RegExp('{pretty-name}', 'g'), pretty_name)
                               .replace(new RegExp('{option-' + norm + '}', 'g'), ' selected')
                               .replace(new RegExp('{option-([a-z]*)}', 'g'), '')
                               .replace(new RegExp('{height}', 'g'), height)
                               .replace(new RegExp('{margin}', 'g'), margin);
        }
        else
        {
            var data_type = "categorical";

            if (hasSettings)
            {
                var height = layer_settings['height'];
                var margin = layer_settings['margin'];
            }
            else
            {
                var height = getNamedLayerDefaults(layer_name, 'height', 80, group);
                var margin = 15;
            }

            var template = '<tr samples-group-name="{group}" samples-layer-name="{name}" data-type="{data-type}">' +
                '<td><img class="drag-icon" src="images/drag.gif" /></td>' +
                '<td title="{name}" class="titles">{short-name}</td>' +
                '<td style="width: 50px;">n/a</td>' +
                '<td style="width: 50px;">n/a</td>' +
                '<td style="width: 50px;">n/a</td>' +
                '<td><input class="input-height form-control form-control-sm" type="text" size="3" value="{height}"></input></td>' +
                '<td><input class="input-margin form-control form-control-sm" type="text" size="3" value="{margin}"></input></td>' +
                '<td style="width: 50px;">n/a</td>' +
                '<td style="width: 50px;">n/a</td>' +
                '<td><input type="checkbox" class="layer_selectors"></input></td>' +
                '</tr>';

            template = template.replace(new RegExp('{name}', 'g'), layer_name)
                               .replace(new RegExp('{group}', 'g'), group)
                               .replace(new RegExp('{data-type}', 'g'), data_type)
                               .replace(new RegExp('{short-name}', 'g'), short_name)
                               .replace(new RegExp('{height}', 'g'), height)
                               .replace(new RegExp('{margin}', 'g'), margin);
        }
        $('#tbody_samples').append(template);
    }

    $('#tbody_samples .normalization').each((index, select) => {
        if ($(select).closest('tr').find('.input-min').val() == 'null' && $(select).closest('tr').find('.input-max').val() == 'null') {
            $(select).trigger('change');
        }
    });

    $('.colorpicker').colpick({
        layout: 'hex',
        submit: 0,
        colorScheme: 'light',
        onChange: function(hsb, hex, rgb, el, bySetColor) {
            $(el).css('background-color', '#' + hex);
            $(el).attr('color', '#' + hex);

            if (!bySetColor) $(el).val(hex);
        }
    }).keyup(function() {
        $(this).colpickSetColor(this.value);
    });
}

function drawSamples() {
    $('#samples').empty();
    drawSamplesLayers(serializeSettings());
    var samples_tree = document.getElementById('samples_tree');
    if (samples_tree)
    {
        samples_tree.addEventListener('click', lineClickHandler, false);
        samples_tree.addEventListener('contextmenu', lineContextMenuHandler, false);
    }
}


function drawSamplesLayers(settings) {
    for (var i = settings['samples-layer-order'].length - 1; i >= 0; i--) {
        if (!is_sample_group_visible(settings['samples-layer-order'][i]['group'])) {
            settings['samples-layer-order'].splice(i, 1);
        }
    }

    var _samples_information_dict = jQuery.extend(true, {}, samples_information_dict); // keep original

    // calculate sample information layer boundaries
    var samples_layer_boundaries = [];

    for (var i=0; i < settings['samples-layer-order'].length; i++)
    {
        var samples_layer_name     = settings['samples-layer-order'][i]['layer_name'];
        var group                  = settings['samples-layer-order'][i]['group'];
        var samples_layer_settings = settings['samples-layers'][group][samples_layer_name];

        var start = (samples_layer_settings['height'] == 0) ? 0 : samples_layer_settings['margin'];
        var end   = start + samples_layer_settings['height'];

        if (i > 0)
        {
            start += samples_layer_boundaries[i-1][1];
            end   += samples_layer_boundaries[i-1][1];
        }

        samples_layer_boundaries.push([start,end]);
    }

    var gradient_done = false;

    if (settings['tree-type'] == 'phylogram') {
        gradient_done = true;
        // there is no need for bacnground gradient for phylogram tree
    }

    var sample_xy = {};

    var samples_start = -1;
    var samples_end = -1;

    for (var j = 0; j < settings['layer-order'].length; j++) {
        var layer_index = j+1;
        var pindex = settings['layer-order'][j];
        var sample_name = getLayerName(pindex);

        if (settings['layers'][pindex]['height'] == 0)
            continue;

        sample_xy[sample_name] = {
            'x': layer_boundaries[layer_index][0] + (layer_boundaries[layer_index][1] - layer_boundaries[layer_index][0]) / 2,
            'y': (samples_layer_boundaries.length > 0) ? 0 - samples_layer_boundaries[samples_layer_boundaries.length-1][1] : 0,
        }

        if(!gradient_done)
        {
            drawGradientBackground(layer_boundaries[layer_index][0]);
            gradient_done = true;
        }

        for (var i=0; i < settings['samples-layer-order'].length; i++)
        {
            var samples_layer_name     = settings['samples-layer-order'][i]['layer_name'];
            var group                  = settings['samples-layer-order'][i]['group'];
            var samples_layer_settings = settings['samples-layers'][group][samples_layer_name];
            var samples_pretty_name    = (samples_layer_name.indexOf('!') > -1) ? samples_layer_name.split('!')[0] : samples_layer_name;

            if (samples_layer_settings['height'] == 0) {
                continue;
            }

            if (!_samples_information_dict[group].hasOwnProperty(sample_name) || !_samples_information_dict[group][sample_name].hasOwnProperty(samples_layer_name)) {
                continue;
            }

            if (samples_start == -1) {
                samples_start = layer_index;
            } else {
                samples_start = Math.min(samples_start, layer_index);
            }
            samples_end = Math.max(layer_index);

            if (samples_layer_settings['data-type'] == 'numeric')
            {
                var value = _samples_information_dict[group][sample_name][samples_layer_name];
                var min = parseFloat(samples_layer_settings['min']['value']);
                var max = parseFloat(samples_layer_settings['max']['value']);
                var ratio;
                if (value > max) {
                    ratio = 1;
                }
                else if (value < min || (max - min) == 0)  {
                    ratio = 0;
                }
                else {
                    ratio = (value - min) / (max - min);
                }

                var size;
                var color;
                if (samples_layer_settings['type'] == 'intensity')
                {
                    var size = samples_layer_settings['height'];
                    var color = getGradientColor(samples_layer_settings['color-start'], samples_layer_settings['color'],  ratio);
                }
                else
                {
                    // bar
                    var size = ratio * samples_layer_settings['height'];
                    var color = samples_layer_settings['color'];
                }

                var start = layer_boundaries[layer_index][0];
                var width = layer_boundaries[layer_index][1] - layer_boundaries[layer_index][0];

                var rect = drawPhylogramRectangle('samples',
                    'samples',
                    start,
                    0 - samples_layer_boundaries[i][0] - (size / 2),
                    size,
                    width,
                    color,
                    1,
                    true);

                rect.setAttribute('sample-name', sample_name);
                rect.setAttribute('sample-group', group);
                rect.setAttribute('layer-name', samples_layer_name);
            }
            else if (samples_layer_settings['data-type'] == 'stack-bar')
            {
                var norm = samples_layer_settings['normalization'];
                var stack_bar_items = _samples_information_dict[group][sample_name][samples_layer_name].split(';');
                var _sum = 0;
                for (var _j=0; _j < stack_bar_items.length; _j++)
                {
                    if (norm == 'sqrt')
                    {
                        stack_bar_items[_j] = Math.sqrt(parseFloat(stack_bar_items[_j]));
                    }
                    else if (norm == 'log')
                    {
                        stack_bar_items[_j] = log10(parseFloat(stack_bar_items[_j]) + 1);
                    }
                    else
                    {
                        stack_bar_items[_j] = parseFloat(stack_bar_items[_j]);
                    }

                    _sum = _sum + stack_bar_items[_j];
                }

                if (_sum > 0) {
                    for (var _j=0; _j < stack_bar_items.length; _j++)
                    {
                        stack_bar_items[_j] = stack_bar_items[_j] / _sum;
                    }
                }

                var offset = 0;
                for (var _i=0; _i < stack_bar_items.length; _i++)
                {
                    if (!isNumber(stack_bar_items[_i]))
                        continue;

                    let bar_name = samples_layer_name.split('!')[1].split(';')[_i];
                    var color = samples_stack_bar_colors[group][samples_layer_name][bar_name];
                    var size  = (samples_layer_boundaries[i][1] - samples_layer_boundaries[i][0]) * stack_bar_items[_i];

                    var rect = drawPhylogramRectangle('samples',
                        'samples',
                        layer_boundaries[layer_index][0],
                        0 - samples_layer_boundaries[i][0] - offset - (size / 2),
                        size,
                        layer_boundaries[layer_index][1] - layer_boundaries[layer_index][0],
                        color,
                        1,
                        true);

                    rect.setAttribute('sample-name', sample_name);
                    rect.setAttribute('sample-group', group);
                    rect.setAttribute('layer-name', samples_layer_name);
                    rect.setAttribute('bar-name', bar_name);

                    offset = offset + size;
                }
            }
            else
            {
                // categorical
                var value = _samples_information_dict[group][sample_name][samples_layer_name];

                if (value == null || value == 'null' || value == '') {
                    value == 'None';
                }

                if (typeof samples_categorical_colors[group][samples_layer_name] === 'undefined') {
                    samples_categorical_colors[group][samples_layer_name] = {};
                }

                if (typeof samples_categorical_colors[group][samples_layer_name][value] === 'undefined')
                {
                    samples_categorical_colors[group][samples_layer_name][value] = randomColor({luminosity: 'dark'});
                }

                var color = samples_categorical_colors[group][samples_layer_name][value];
                var size  = samples_layer_boundaries[i][1] - samples_layer_boundaries[i][0];

                var rect = drawPhylogramRectangle('samples',
                    'samples',
                    layer_boundaries[layer_index][0],
                    0 - samples_layer_boundaries[i][0] - (size / 2),
                    size,
                    layer_boundaries[layer_index][1] - layer_boundaries[layer_index][0],
                    color,
                    1,
                    true);

                rect.setAttribute('sample-name', sample_name);
                rect.setAttribute('sample-group', group);
                rect.setAttribute('layer-name', samples_layer_name);
            }
        }
    }

    if (!(samples_start == -1 || samples_end == -1))
    {
        // draw sample backgrounds and titles.
        for (var i=0; i < settings['samples-layer-order'].length; i++)
        {
            var samples_layer_name     = settings['samples-layer-order'][i]['layer_name'];
            var group                  = settings['samples-layer-order'][i]['group'];
            var samples_layer_settings = settings['samples-layers'][group][samples_layer_name];
            var samples_pretty_name    = (samples_layer_name.indexOf('!') > -1) ? samples_layer_name.split('!')[0] : samples_layer_name;
            var min = parseFloat(samples_layer_settings['min']['value']);
            var max = parseFloat(samples_layer_settings['max']['value']);

            if (samples_layer_settings['height'] == 0) {
                continue;
            }

            if (samples_layer_settings['data-type'] == 'numeric')
            {
                if (samples_layer_settings['type'] != 'intensity')
                {
                    var start = samples_layer_boundaries[i][0];
                    var end   = samples_layer_boundaries[i][1];

                    drawPhylogramRectangle('samples',
                        'samples_background',
                        layer_boundaries[samples_start][0],
                        0 - end + (end - start) / 2,
                        end - start,
                        layer_boundaries[samples_end][1] - layer_boundaries[samples_start][0],
                        samples_layer_settings['color'],
                        0.2,
                        false);
                }

                let font_size = Math.min(samples_layer_settings['height'] / 3, parseFloat(settings['max-font-size-label']));
                drawText('samples', {
                        'x': layer_boundaries[samples_end][1] + 20,
                        'y': 0 - (samples_layer_boundaries[i][0] + samples_layer_boundaries[i][1]) / 2 + font_size / 6
                    },
                    getNamedLayerDefaults(samples_layer_name, 'pretty_name', samples_layer_name),
                    font_size + 'px',
                    'left',
                    samples_layer_settings['color'],
                    'baseline');

                drawText('samples', {
                    'x': layer_boundaries[samples_end][1] + 20,
                    'y': 0 - samples_layer_boundaries[i][1] + samples_layer_settings['height'] / 6
                }, max , samples_layer_settings['height'] / 6 + 'px', 'left', '#000000', 'baseline');

                drawText('samples', {
                    'x': layer_boundaries[samples_end][1] + 20,
                    'y': 0 - samples_layer_boundaries[i][0]
                }, min , samples_layer_settings['height'] / 6 + 'px', 'left', '#000000', 'baseline');

            }
            else if (samples_layer_settings['data-type'] == 'stack-bar')
            {
                let font_size = Math.min(samples_layer_settings['height'] / 3, parseFloat(settings['max-font-size-label']));
                drawText('samples', {
                        'x': layer_boundaries[samples_end][1] + 20,
                        'y': 0 - (samples_layer_boundaries[i][0] + samples_layer_boundaries[i][1]) / 2 + font_size / 6
                    },
                    getNamedLayerDefaults(samples_pretty_name, 'pretty_name', samples_pretty_name),
                    font_size + 'px',
                    'left',
                    '#919191',
                    'baseline');
            }
            else
            {
                let font_size = Math.min(samples_layer_settings['height'], parseFloat(settings['max-font-size-label']));
                drawText('samples', {
                        'x': layer_boundaries[samples_end][1] + 20,
                        'y': 0 - (samples_layer_boundaries[i][0] + samples_layer_boundaries[i][1]) / 2 + font_size / 2
                    },
                    getNamedLayerDefaults(samples_layer_name, 'pretty_name', samples_layer_name),
                    font_size + 'px',
                    'left',
                    samples_layer_settings['color'],
                    'baseline');
            }
        }
    }
    drawSamplesTree(settings, sample_xy);
}

function drawSamplesTree(settings, sample_xy)
{
    createBin('samples', 'samples_tree');
    var samples_order = settings['samples-order'];

    if (samplesClusteringData['newick'] == null || samplesClusteringData['newick'] == '')
        return;

    var newick = samplesClusteringData['newick'];
    var t = new Tree();
    t.Parse(newick, settings['samples-edge-length-normalization']);
    t.ComputeDepths();

    var use_edge_lengths = t.has_edge_lengths;
    if (settings['samples-ignore-branch-length'])
    {
        use_edge_lengths = false;
    }

    if (use_edge_lengths)
    {
        max_path_length = 0;
        var n = new PreorderIterator(t.root);
        var q = n.Begin();
        while (q != null) {
            var d = q.edge_length;
            if (d < 0.00001) {
                d = 0.0;
            }
            if (q != t.root) {
                q.path_length = q.ancestor.path_length + d;
            }

            max_path_length = Math.max(max_path_length, q.path_length);
            q = n.Next();
        }

        var n = new NodeIterator(t.root);
        var q = n.Begin();
    } else {
        var max_depth = t.root.depth;
        max_path_length = max_depth;

        for (var i=0; i < t.nodes.length; i++)
        {
            var n = t.nodes[i];
            if (n != t.root)
            {
                n.path_length = max_depth - n.depth;
            }
        }
    }

    var samples_top = -1;
    var samples_height = parseInt(settings['samples-tree-height']);

    samples_id_to_node_map = new Array();

    var n = new NodeIterator(t.root);
    var q = n.Begin();
    while (q != null)
    {
        samples_id_to_node_map[q.id] = q;

        if (q.IsLeaf())
        {
            if (!sample_xy.hasOwnProperty(q.label))
            {
                var _message = "Error: Layer additional order '" + samples_order + "' has leaf named '" + q.label + "' but it is not a valid layer name. Falling back to custom order. Maybe that layer is hidden??";
                toastr.warning(_message, "", { 'timeOut': '0', 'extendedTimeOut': '0' });
                samplesClusteringData = {'newick': '', 'basic': null};
                $('#samples_order').val('custom').trigger('change');
                return;
            }

            if (samples_top == -1)
            {
                samples_top = sample_xy[q.label]['y'] - samples_height;
            }
            q.xy['x'] = sample_xy[q.label]['x'];
            q.xy['y'] = samples_top + (samples_height * q.path_length/max_path_length);
        }
        else
        {
            var pl = q.child.xy;
            var pr = q.child.GetRightMostSibling().xy;
            q.xy['x'] = pl['x'] + (pr['x'] - pl['x']) / 2;
            q.xy['y'] = samples_top + (samples_height * q.path_length/max_path_length);
        }
        q=n.Next();
    }

    var n = new NodeIterator(t.root);
    var q = n.Begin();

    var _lines = [];
    while (q != null)
    {
        if (q.IsLeaf())
        {
            var p0 = q.xy
            var p1 = [];
            var anc = q.ancestor;
            if (anc) {
                p1['y'] = anc.xy['y'];
                p1['x'] = p0['x'];

                _lines.push(drawLine('samples_tree', q, p0, p1));
            }
        }
        else
        {
            var p0 = [];
            var p1 = [];

            p0['y'] = q.xy['y'];
            p0['x'] = q.xy['x'];

            var anc = q.ancestor;
            if (anc) {
                p1['y'] = anc.xy['y'];
                p1['x'] = p0['x'];
                _lines.push(drawLine('samples_tree', q, p0, p1));
            }

            // vertical line
            var pl = q.child.xy;
            var pr = q.child.GetRightMostSibling().xy;

            p0['y'] = p0['y'];
            p0['x'] = pl['x'];
            p1['y'] = p0['y'];
            p1['x'] = pr['x'];

            _lines.push(drawLine('samples_tree', q, p0, p1, true));
        }
        q=n.Next();
    }

    _lines.forEach(function(_line) {
        _line.setAttribute('id', 'samples_' + _line.getAttribute('id'));
        var new_line = drawLine('samples_tree', {'id': 0}, {'x': 0, 'y': 0}, {'x': 0, 'y': 0});
        new_line.setAttribute('id', _line.getAttribute('id') + '_clone');
        new_line.setAttribute('d', _line.getAttribute('d'));
        new_line.classList.add('clone');
        new_line.style['stroke-width'] = '20px';
        new_line.style['stroke-opacity'] = '0';
        new_line.setAttribute('pointer-events', 'all');
    });

    // draw guide lines for samples tree
    if (use_edge_lengths) {
        for (var i=0; i < t.nodes.length; i++) {
            q = t.nodes[i];
            if (q.IsLeaf())
            {
                var _line = drawLine('samples_tree', q, q.xy, {'x': q.xy['x'], 'y': samples_top + samples_height});
                _line.setAttribute('id', _line.getAttribute + '_guide');
                _line.style['stroke-opacity'] = '0.2';
                _line.setAttribute('stroke-dasharray', '1,1');
            }
        }
    }
}

function drawGradientBackground(start)
{
    // draw gradient over labels
    createGradient(document.getElementById('svg'),'gradient1',[
        {offset:'0%', style:'stop-color:rgb(255,255,255);stop-opacity:0'},
        {offset:'100%', style:'stop-color:rgb(255,255,255);stop-opacity:1'},
    ]);

    var grect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
    grect.setAttribute('id', 'label_gradient1');
    grect.setAttribute('fill', 'url(#gradient1)');
    grect.setAttribute('x', start - 400);
    grect.setAttribute('y', 0 - total_radius);
    grect.setAttribute('width', 410);
    grect.setAttribute('height', total_radius);
    grect.setAttribute('stroke-width', '0px');
    document.getElementById('samples').appendChild(grect);

    var rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
    rect.setAttribute('id', 'label_gradient2');
    rect.setAttribute('fill', '#ffffff');
    rect.setAttribute('x', start);
    rect.setAttribute('y', 0 - total_radius);
    rect.setAttribute('width', total_radius - start);
    rect.setAttribute('height', total_radius);
    rect.setAttribute('stroke-width', '0px');
    rect.setAttribute('pointer-events', 'none');

    document.getElementById('samples').appendChild(rect);
}
