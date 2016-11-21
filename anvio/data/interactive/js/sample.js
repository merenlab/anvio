var samples_categorical_colors = {};
var samples_stack_bar_colors = {};
var samples_information_dict;
var samples_order_dict;

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

        if (this.value == 'custom') 
            return;

        var organization = samples_order_dict[this.value];
        var samples_new_order;

        // get new sample order
        if (organization['basic'] != null && organization['basic'] != "")
        {
            samples_new_order = organization['basic'].split(',');
        }
        else
        {
            samples_new_order = get_newick_leaf_order(organization['newick']);
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
            var layer_id = getLayerId(new_order[i]);
            var detached_row = $('#height' + layer_id).closest('tr').detach();
            $('#tbody_layers').append(detached_row);
        }
    });
});

function buildSamplesTable(samples_layer_order, samples_layers) {
    var first_sample = Object.keys(samples_information_dict)[0];

    if (typeof first_sample === 'undefined')
    {
        return;
    }
    
    var order_from_state = true;
    var all_information_layers = Object.keys(samples_information_dict[first_sample]);

    if (typeof(samples_layer_order) === 'undefined') {
        samples_layer_order = all_information_layers;
        order_from_state = false;
    } else {
        // add missing layers to state order list
        all_information_layers.forEach(function(_lname) {
            if (samples_layer_order.indexOf(_lname) == -1) {
                samples_layer_order.push(_lname);
            }
        });
    }
    
    $('#tbody_samples').empty();

    for (var i=0; i < samples_layer_order.length; i++)
    {
        var layer_name  = samples_layer_order[i];
        var pretty_name = getNamedLayerDefaults(layer_name, 'pretty_name', layer_name);
        pretty_name = (pretty_name.indexOf('!') > -1) ? pretty_name.split('!')[0] : pretty_name;
        
        var short_name = (pretty_name.length > 10) ? pretty_name.slice(0,10) + "..." : pretty_name;

        var hasSettings = false;
        if (typeof(samples_layers) !== 'undefined' && typeof(samples_layers[layer_name]) !== 'undefined') {
            hasSettings = true;
            layer_settings = samples_layers[layer_name];
        }

        if (isNumber(samples_information_dict[first_sample][layer_name]))
        {
            var data_type = "numeric";
            
            if (hasSettings) 
            {
                var norm         = layer_settings['normalization'];
                var min          = layer_settings['min']['value'];
                var max          = layer_settings['max']['value'];
                var min_disabled = layer_settings['min']['disabled'];
                var max_disabled = layer_settings['max']['disabled'];
                var height       = layer_settings['height'];
                var color        = layer_settings['color'];
                var margin       = layer_settings['margin'];
                var color_start  = layer_settings['color-start'];
                var type         = layer_settings['type'];
            }
            else
            {
                var norm         = getNamedLayerDefaults(layer_name, 'norm', 'none');
                var min          = 0;
                var max          = 0;
                var min_disabled = true;
                var max_disabled = true;
                var height       = getNamedLayerDefaults(layer_name, 'height', 500);
                var color        = getNamedLayerDefaults(layer_name, 'color', '#919191');
                var margin       = 15;
                var color_start  = "#FFFFFF";
                var type         = "bar";
            }

            var template = '<tr samples-layer-name="{name}" data-type="{data-type}">' +
                '<td><img class="drag-icon" src="images/drag.gif" /></td>' +
                '<td title="{name}" class="titles">{short-name}</td>' +
                '<td><div class="colorpicker picker_start" color="{color-start}" style="background-color: {color-start}; {color-start-hide}"></div><div class="colorpicker" color="{color}" style="background-color: {color}"></div></td>' +
                '<td style="width: 50px;">' +
                '    <select style="width: 50px;" class="type" onChange="togglePickerStart(this);">' +
                '        <option value="bar"{option-type-bar}>Bar</option>' +
                '        <option value="intensity"{option-type-intensity}>Intensity</option>' +
                '    </select>' +
                '</td>' +
                '<td>' +
                '    <select onChange="clearMinMax(this);" class="normalization">' +
                '        <option value="none"{option-none}>none</option>' +
                '        <option value="sqrt"{option-sqrt}>sqrt</option>' +
                '        <option value="log"{option-log}>log</option>' +
                '    </select>' +
                '</td>' +
                '<td><input class="input-height" type="text" size="3" value="{height}"></input></td>' +
                '<td><input class="input-margin" type="text" size="3" value="{margin}"></input></td>' +
                '<td><input class="input-min" type="text" size="4" value="{min}"{min-disabled}></input></td>' +
                '<td><input class="input-max" type="text" size="4" value="{max}"{min-disabled}></input></td>' +
                '<td><input type="checkbox" class="layer_selectors"></input></td>' +
                '</tr>';

            template = template.replace(new RegExp('{name}', 'g'), layer_name)
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
                               .replace(new RegExp('{min-disabled}', 'g'), (min_disabled) ? ' disabled': '')
                               .replace(new RegExp('{max-disabled}', 'g'), (max_disabled) ? ' disabled': '')
                               .replace(new RegExp('{margin}', 'g'), margin);
            
            $('#tbody_samples').append(template); 
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
                var norm   = getNamedLayerDefaults(layer_name, 'norm', 'none');
                var height = getNamedLayerDefaults(layer_name, 'height', 500);
                var margin = 15;

                // pick random color for stack bar items
                if (!(layer_name in samples_stack_bar_colors))
                {
                    samples_stack_bar_colors[layer_name] = new Array();
                    for (var j=0; j < layer_name.split(";").length; j++)
                    {
                        samples_stack_bar_colors[layer_name].push(randomColor());
                    } 
                }  
            }

            var template = '<tr samples-layer-name="{name}" data-type="{data-type}">' +
                '<td><img class="drag-icon" src="images/drag.gif" /></td>' +
                '<td title="{pretty-name}" class="titles">{short-name}</td>' +
                '<td>n/a</td>' +
                '<td style="width: 50px;">n/a</td>' +
                '<td>' +
                '    <select onChange="clearMinMax(this);" class="normalization">' +
                '        <option value="none"{option-none}>none</option>' +
                '        <option value="sqrt"{option-sqrt}>sqrt</option>' +
                '        <option value="log"{option-log}>log</option>' +
                '    </select>' +
                '</td>' +
                '<td><input class="input-height" type="text" size="3" value="{height}"></input></td>' +
                '<td><input class="input-margin" type="text" size="3" value="{margin}"></input></td>' +
                '<td>n/a</td>' +
                '<td>n/a</input></td>' +
                '<td><input type="checkbox" class="layer_selectors"></input></td>' +
                '</tr>';

            template = template.replace(new RegExp('{name}', 'g'), layer_name)
                               .replace(new RegExp('{data-type}', 'g'), data_type)
                               .replace(new RegExp('{short-name}', 'g'), short_name)
                               .replace(new RegExp('{pretty-name}', 'g'), pretty_name)
                               .replace(new RegExp('{option-' + norm + '}', 'g'), ' selected')
                               .replace(new RegExp('{option-([a-z]*)}', 'g'), '')
                               .replace(new RegExp('{height}', 'g'), height)
                               .replace(new RegExp('{margin}', 'g'), margin);
        
            if (order_from_state)
            {
                $('#tbody_samples').append(template);
            }
            else
            {
                $('#tbody_samples').prepend(template);
            }
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
                var height = getNamedLayerDefaults(layer_name, 'height', 80);
                var margin = 15;
            }

            var template = '<tr samples-layer-name="{name}" data-type="{data-type}">' +
                '<td><img class="drag-icon" src="images/drag.gif" /></td>' +
                '<td title="{name}" class="titles">{short-name}</td>' +
                '<td>n/a</td>' +
                '<td style="width: 50px;">n/a</td>' +
                '<td>n/a</td>' +
                '<td><input class="input-height" type="text" size="3" value="{height}"></input></td>' +
                '<td><input class="input-margin" type="text" size="3" value="{margin}"></input></td>' +
                '<td>n/a</td>' +
                '<td>n/a</input></td>' +
                '<td><input type="checkbox" class="layer_selectors"></input></td>' +
                '</tr>';

            template = template.replace(new RegExp('{name}', 'g'), layer_name)
                               .replace(new RegExp('{data-type}', 'g'), data_type)
                               .replace(new RegExp('{short-name}', 'g'), short_name)
                               .replace(new RegExp('{height}', 'g'), height)
                               .replace(new RegExp('{margin}', 'g'), margin);
        
            if (order_from_state)
            {
                $('#tbody_samples').append(template);
            }
            else
            {
                $('#tbody_samples').prepend(template);
            }
        }  
    }

    $('.colorpicker').colpick({
        layout: 'hex',
        submit: 0,
        colorScheme: 'dark',
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
    }
}

function drawSamplesLayers(settings) {
    var samples_layer_max = {};
    var samples_layer_min = {};

    var _samples_information_dict = jQuery.extend(true, {}, samples_information_dict); // keep original

    for (sample in _samples_information_dict)
    {
        for (layer in _samples_information_dict[sample])
        {
            if (settings['samples-layers'][layer]['data-type'] == 'numeric') 
            {
                var norm = settings['samples-layers'][layer]['normalization'];

                if (norm == 'sqrt')
                {
                    _samples_information_dict[sample][layer] = Math.sqrt(parseFloat(_samples_information_dict[sample][layer]));
                }
                else if (norm == 'log')
                {
                    _samples_information_dict[sample][layer] = log10(parseFloat(_samples_information_dict[sample][layer]) + 1);
                }

                if (typeof samples_layer_max[layer] === 'undefined' || parseFloat(_samples_information_dict[sample][layer]) > samples_layer_max[layer])
                {
                    samples_layer_max[layer] = parseFloat(_samples_information_dict[sample][layer]);
                }
            }
            else if (settings['samples-layers'][layer]['data-type'] == 'stack-bar') 
            {
                var norm = settings['samples-layers'][layer]['normalization'];

                var stack_bar_items = _samples_information_dict[sample][layer].split(';');
                var _sum = 0;
                for (var j=0; j < stack_bar_items.length; j++)
                {
                    if (norm == 'sqrt')
                    {
                        stack_bar_items[j] = Math.sqrt(parseFloat(stack_bar_items[j]));
                    }
                    else if (norm == 'log')
                    {
                        stack_bar_items[j] = log10(parseFloat(stack_bar_items[j]) + 1);
                    }
                    else
                    {
                        stack_bar_items[j] = parseFloat(stack_bar_items[j]);
                    }

                    _sum = _sum + stack_bar_items[j];
                }

                for (var j=0; j < stack_bar_items.length; j++)
                {
                    stack_bar_items[j] = stack_bar_items[j] / _sum;
                }

                _samples_information_dict[sample][layer] = stack_bar_items;
            }
            else
            {
                //categorical
                if (typeof samples_categorical_colors[layer] === 'undefined')
                    samples_categorical_colors[layer] = {};
            }
        }
    }

    // calculate sample information layer boundaries
    var samples_layer_boundaries = [];

    for (var i=0; i < settings['samples-layer-order'].length; i++)
    {
        var samples_layer_name     = settings['samples-layer-order'][i];
        var samples_layer_settings = settings['samples-layers'][samples_layer_name];
        
        if (samples_layer_settings['min']['disabled'])
        {
            $('#tbody_samples [samples-layer-name=' + samples_layer_name + '] .input-min').prop('disabled', false);
            $('#tbody_samples [samples-layer-name=' + samples_layer_name + '] .input-max').prop('disabled', false).val(samples_layer_max[samples_layer_name])
            samples_layer_min[samples_layer_name] = 0;
        }
        else
        {
            samples_layer_max[samples_layer_name] = samples_layer_settings['max']['value'];
            samples_layer_min[samples_layer_name] = samples_layer_settings['min']['value'];
        }

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

    var sample_xy = {};

    var samples_start = -1;
    var samples_end = -1;

    for (var j = 0; j < settings['layer-order'].length; j++) {
        var layer_index = j+1;
        var pindex = settings['layer-order'][j];
        var sample_name = getLayerName(pindex);

        sample_xy[sample_name] = {
            'x': layer_boundaries[layer_index][0] + (layer_boundaries[layer_index][1] - layer_boundaries[layer_index][0]) / 2,
            'y': (samples_layer_boundaries.length > 0) ? 0 - samples_layer_boundaries[samples_layer_boundaries.length-1][1] : 0,
        }

        if (!(sample_name in samples_information_dict)) // skip if not sample
            continue;

        // update start once
        if (samples_start == -1)
            samples_start = layer_index;
        samples_end = layer_index;

        if(!gradient_done)
        {
            drawGradientBackground(layer_boundaries[layer_index][0]);
            gradient_done = true;
        }

        for (var i=0; i < settings['samples-layer-order'].length; i++)
        {
            var samples_layer_name     = settings['samples-layer-order'][i];
            var samples_layer_settings = settings['samples-layers'][samples_layer_name];
            var samples_pretty_name    = (samples_layer_name.indexOf('!') > -1) ? samples_layer_name.split('!')[0] : samples_layer_name;

            if (samples_layer_settings['height'] == 0) {
                continue;
            }

            if (samples_layer_settings['data-type'] == 'numeric') 
            {
                var value = _samples_information_dict[sample_name][samples_layer_name];
                var min = samples_layer_min[samples_layer_name];
                var max = samples_layer_max[samples_layer_name];
                
                var ratio;
                if (value > max) {
                    ratio = 1;
                }
                else if (value < min) {
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
                rect.setAttribute('layer-name', samples_layer_name);
            }
            else if (samples_layer_settings['data-type'] == 'stack-bar') 
            {
                var stack_bar_items = _samples_information_dict[sample_name][samples_layer_name];

                var offset = 0;
                for (var _i=0; _i < stack_bar_items.length; _i++)
                {
                    var color = samples_stack_bar_colors[samples_layer_name][_i];
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
                    rect.setAttribute('layer-name', samples_layer_name);

                    offset = offset + size;
                }
            }
            else
            {
                // categorical
                var value = _samples_information_dict[sample_name][samples_layer_name];

                if (typeof samples_categorical_colors[samples_layer_name][value] === 'undefined')
                {
                    samples_categorical_colors[samples_layer_name][value] = randomColor({luminosity: 'dark'});
                }

                var color = samples_categorical_colors[samples_layer_name][value];
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
                rect.setAttribute('layer-name', samples_layer_name);
            }
        }
    }

    // draw sample backgrounds and titles.
    for (var i=0; i < settings['samples-layer-order'].length; i++)
    {
        var samples_layer_name     = settings['samples-layer-order'][i];
        var samples_layer_settings = settings['samples-layers'][samples_layer_name];
        var samples_pretty_name    = (samples_layer_name.indexOf('!') > -1) ? samples_layer_name.split('!')[0] : samples_layer_name;
        var min = samples_layer_min[samples_layer_name];
        var max = samples_layer_max[samples_layer_name];

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

            drawText('samples', {
                'x': layer_boundaries[samples_end][1] + 20,
                'y': 0 - (samples_layer_boundaries[i][0] + samples_layer_boundaries[i][1]) / 2 + samples_layer_settings['height'] / 6
            }, getNamedLayerDefaults(samples_layer_name, 'pretty_name', samples_layer_name) , samples_layer_settings['height'] / 3 + 'px', 'left', samples_layer_settings['color'], 'baseline');

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
            drawText('samples', {
                'x': layer_boundaries[samples_end][1] + 20,
                'y': 0 - (samples_layer_boundaries[i][0] + samples_layer_boundaries[i][1]) / 2 + samples_layer_settings['height'] / 6
            }, getNamedLayerDefaults(samples_pretty_name, 'pretty_name', samples_pretty_name), samples_layer_settings['height'] / 3 + 'px', 'left', '#919191', 'baseline');
        }
        else
        {
            drawText('samples', {
                'x': layer_boundaries[samples_end][1] + 20,
                'y': 0 - (samples_layer_boundaries[i][0] + samples_layer_boundaries[i][1]) / 2 + samples_layer_settings['height'] / 2
            }, getNamedLayerDefaults(samples_layer_name, 'pretty_name', samples_layer_name), samples_layer_settings['height'] + 'px', 'left', samples_layer_settings['color'], 'baseline');   
        }
    }


    drawSamplesTree(settings, sample_xy);
}

function drawSamplesTree(settings, sample_xy)
{
    createBin('samples', 'samples_tree');
    var samples_order = settings['samples-order'];

    if (!samples_order_dict.hasOwnProperty(samples_order) || samples_order_dict[samples_order]['newick'] == null || samples_order_dict[samples_order]['newick'] == '')
        return;

    var newick = samples_order_dict[samples_order]['newick'];
    var t = new Tree();
    t.Parse(newick, settings['samples-edge-length-normalization']);
    t.ComputeDepths();
    t.ComputeWeights(t.root);

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
                var _message = "Error: Sample order '" + samples_order + "' has leaf named '" + q.label + "' but it is not a valid layer name.";
                console.log(_message);
                toastr.error(_message, "", { 'timeOut': '0', 'extendedTimeOut': '0' });
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
