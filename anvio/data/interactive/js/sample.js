var metadata_categorical_colors = {};
var metadata;
var organizations;

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
    $('#sample_organization').change(function() {
        if (this.value == 'custom') 
            return;

        $('#btn_redraw_metadata').prop('disabled', true);

        var organization = organizations[this.value];
        var sample_order;

        // get new sample order
        if (organization['order'] != "")
        {
            sample_order = organization['order'].split(',');
        }
        else
        {
            sample_order = get_newick_leaf_order(organization['newick']);
        }
        $.map(sample_order, $.trim);

        for(var i=0; i < sample_order.length; i++)
        {
            var layer_id = getLayerId(sample_order[i]);
            var detached_row = $('#height' + layer_id).closest('tr').detach();
            $('#tbody_layers').append(detached_row);
        }
    });
});

function buildMetadataTable(metadata_layer_order, metadata_layers) {
    var first_sample = Object.keys(metadata)[0];
    
    if (typeof(metadata_layer_order) === 'undefined') {
        metadata_layer_order = Object.keys(metadata[first_sample]); // get layer order from first sample's metadata
    }
    
    $('#tbody_metadata').empty();

    for (var i=0; i < metadata_layer_order.length; i++)
    {
        var layer_name = metadata_layer_order[i];
        var short_name = (layer_name.length > 10) ? layer_name.slice(0,10) + "..." : layer_name;

        var hasSettings = false;
        if (typeof(metadata_layers) !== 'undefined' && typeof(metadata_layers[layer_name]) !== 'undefined') {
            hasSettings = true;
            layer_settings = metadata_layers[layer_name];
        }

        if (isNumber(metadata[first_sample][layer_name]))
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
                var norm = "none";
                var min    = 0;
                var max    = 0;
                var min_disabled = true;
                var max_disabled = true;
                var height = 500;
                var color  = '#919191';
                var margin = 15;
                var color_start = "#FFFFFF";
                var type = "bar";
            }

            var template = '<tr metadata-layer-name="{name}" data-type="{data-type}">' +
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
                var height = 80;
                var margin = 15;
            }

            var template = '<tr metadata-layer-name="{name}" data-type="{data-type}">' +
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
        }

        $('#tbody_metadata').append(template);   
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

function drawMetadataLayers(settings) {
    var metadata_layer_max = {};
    var metadata_layer_min = {};

    var _metadata = jQuery.extend(true, {}, metadata); // keep original

    for (sample in _metadata)
    {
        for (layer in _metadata[sample])
        {
            if (settings['metadata-layers'][layer]['data-type'] == 'numeric') 
            {
                var norm = settings['metadata-layers'][layer]['normalization'];

                if (norm == 'sqrt')
                {
                    _metadata[sample][layer] = Math.sqrt(parseFloat(_metadata[sample][layer]));
                }
                else if (norm == 'log')
                {
                    _metadata[sample][layer] = log10(parseFloat(_metadata[sample][layer]) + 1);
                }

                if (typeof metadata_layer_max[layer] === 'undefined' || parseFloat(_metadata[sample][layer]) > metadata_layer_max[layer])
                {
                    metadata_layer_max[layer] = parseFloat(_metadata[sample][layer]);
                }
            }
            else
            {
                //categorical
                if (typeof metadata_categorical_colors[layer] === 'undefined')
                    metadata_categorical_colors[layer] = {};
            }
        }
    }

    // calculate metadata layer boundaries
    var metadata_layer_boundaries = [];

    for (var i=0; i < settings['metadata-layer-order'].length; i++)
    {
        var metadata_layer_name     = settings['metadata-layer-order'][i];
        var metadata_layer_settings = settings['metadata-layers'][metadata_layer_name];
        
        if (metadata_layer_settings['min']['disabled'])
        {
            $('#tbody_metadata [metadata-layer-name=' + metadata_layer_name + '] .input-min').prop('disabled', false);
            $('#tbody_metadata [metadata-layer-name=' + metadata_layer_name + '] .input-max').prop('disabled', false).val(metadata_layer_max[metadata_layer_name])
            metadata_layer_min[metadata_layer_name] = 0;
        }
        else
        {
            metadata_layer_max[metadata_layer_name] = metadata_layer_settings['max']['value'];
            metadata_layer_min[metadata_layer_name] = metadata_layer_settings['min']['value'];
        }

        var start = metadata_layer_settings['margin'];
        var end   = start + metadata_layer_settings['height'];
        
        if (i > 0)
        {
            start += metadata_layer_boundaries[i-1][1];
            end   += metadata_layer_boundaries[i-1][1];
        }

        metadata_layer_boundaries.push([start,end]);
    }

    var backgrounds_done = false;
    var gradient_done = false;

    var sample_xy = {};

    for (var j = 0; j < settings['layer-order'].length; j++) {
        var layer_index = j+1;
        var pindex = settings['layer-order'][j];
        var sample_name = getLayerName(pindex);

        if (!(sample_name in metadata)) // skip if not sample
            continue;

        if(!gradient_done)
        {
            drawGradientBackground(layer_boundaries[layer_index][0]);
            gradient_done = true;
        }

        sample_xy[sample_name] = {
            'x': layer_boundaries[layer_index][0] + (layer_boundaries[layer_index][1] - layer_boundaries[layer_index][0]) / 2,
            'y': 0 - metadata_layer_boundaries[metadata_layer_boundaries.length-1][1],
        }

        for (var i=0; i < settings['metadata-layer-order'].length; i++)
        {
            var metadata_layer_name     = settings['metadata-layer-order'][i];
            var metadata_layer_settings = settings['metadata-layers'][metadata_layer_name];

            if (metadata_layer_settings['data-type'] == 'numeric') 
            {
                var value = _metadata[sample_name][metadata_layer_name];
                var min = metadata_layer_min[metadata_layer_name];
                var max = metadata_layer_max[metadata_layer_name];
                
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
                if (metadata_layer_settings['type'] == 'intensity')
                {
                    var size = metadata_layer_settings['height'];
                    var color = getGradientColor(metadata_layer_settings['color-start'], metadata_layer_settings['color'],  ratio);
                }
                else
                {
                    // bar
                    var size = ratio * metadata_layer_settings['height'];
                    var color = metadata_layer_settings['color'];

                    if (!backgrounds_done)
                    {
                        var start = metadata_layer_boundaries[i][0];
                        var end   = metadata_layer_boundaries[i][1];

                        drawPhylogramRectangle('metadata',
                            'metadata_background',
                            layer_boundaries[layer_index][0],
                            0 - end + (end - start) / 2,
                            end - start,
                            total_radius - layer_boundaries[layer_index][0],
                            metadata_layer_settings['color'],
                            0.2,
                            false);
                    }
                }
                
                if (!backgrounds_done)
                {
                    drawText('metadata', {
                        'x': total_radius + 20,
                        'y': 0 - (metadata_layer_boundaries[i][0] + metadata_layer_boundaries[i][1]) / 2
                    }, metadata_layer_name , metadata_layer_settings['height'] / 3 + 'px', 'left', metadata_layer_settings['color']);
                    
                    drawText('metadata', {
                        'x': total_radius + 10,
                        'y': 0 - metadata_layer_boundaries[i][1]
                    }, max , metadata_layer_settings['height'] / 6 + 'px', 'left', '#000000', 'text-before-edge');

                    drawText('metadata', {
                        'x': total_radius + 10,
                        'y': 0 - metadata_layer_boundaries[i][0]
                    }, min , metadata_layer_settings['height'] / 6 + 'px', 'left', '#000000', 'text-after-edge');

                }

                var rect = drawPhylogramRectangle('metadata',
                    'metadata',
                    layer_boundaries[layer_index][0],
                    0 - metadata_layer_boundaries[i][0] - (size / 2),
                    size,
                    layer_boundaries[layer_index][1] - layer_boundaries[layer_index][0],
                    color,
                    1,
                    true);

                rect.setAttribute('sample-name', sample_name);
                rect.setAttribute('layer-name', metadata_layer_name);
            }
            else
            {

                // categorical
                var value = _metadata[sample_name][metadata_layer_name];

                if (typeof metadata_categorical_colors[metadata_layer_name][value] === 'undefined')
                {
                    metadata_categorical_colors[metadata_layer_name][value] = randomColor({luminosity: 'dark'});
                }

                var color = metadata_categorical_colors[metadata_layer_name][value];
                var size  = metadata_layer_boundaries[i][1] - metadata_layer_boundaries[i][0];

                var rect = drawPhylogramRectangle('metadata',
                    'metadata',
                    layer_boundaries[layer_index][0],
                    0 - metadata_layer_boundaries[i][0] - (size / 2),
                    size,
                    layer_boundaries[layer_index][1] - layer_boundaries[layer_index][0],
                    color,
                    1,
                    true);

                rect.setAttribute('sample-name', sample_name);
                rect.setAttribute('layer-name', metadata_layer_name);
                
                if (!backgrounds_done)
                {
                    drawText('metadata', {
                        'x': total_radius + 20,
                        'y': 0 - (metadata_layer_boundaries[i][0] + metadata_layer_boundaries[i][1]) / 2
                    }, metadata_layer_name , metadata_layer_settings['height'] + 'px', 'left', metadata_layer_settings['color']);
                }
            }
        }

        backgrounds_done = true;
    }

    drawMetadataTree(settings, sample_xy);
}

function drawMetadataTree(settings, sample_xy)
{
    createBin('metadata', 'metadata_tree');
    var organization_name = settings['organization-name'];

    if (!organizations.hasOwnProperty(organization_name) || organizations[organization_name]['newick'] == '')
        return;

    var newick = organizations[organization_name]['newick'];
    var t = new Tree();
    t.Parse(newick, false);
    t.ComputeDepths();
    t.ComputeWeights();

    var n = new NodeIterator(t.root);
    var q = n.Begin();

    var sample_y = -1;

    while (q != null)
    {
        if (q.IsLeaf())
        {
            q.xy = sample_xy[q.label];

            if (sample_y == -1)
            {
                sample_y = q.xy['y'];
            }
        }
        else
        {
            var pl = q.child.xy;
            var pr = q.child.GetRightMostSibling().xy;

            q.xy['x'] = pl['x'] + (pr['x'] - pl['x']) / 2;
            q.xy['y'] = sample_y - (q.depth * 100);
        }
        q=n.Next();
    }

    var n = new NodeIterator(t.root);
    var q = n.Begin();

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

                drawLine('metadata_tree', q, p0, p1);
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
                drawLine('metadata_tree', q, p0, p1);
            }

            // vertical line
            var pl = q.child.xy;
            var pr = q.child.GetRightMostSibling().xy;

            p0['y'] = p0['y'];
            p0['x'] = pl['x'];
            p1['y'] = p0['y'];
            p1['x'] = pr['x'];

            drawLine('metadata_tree', q, p0, p1);
        }
        q=n.Next();
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
    document.getElementById('metadata').appendChild(grect);

    var rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
    rect.setAttribute('id', 'label_gradient2');
    rect.setAttribute('fill', '#ffffff');
    rect.setAttribute('x', start);
    rect.setAttribute('y', 0 - total_radius);
    rect.setAttribute('width', total_radius - start);
    rect.setAttribute('height', total_radius);
    rect.setAttribute('stroke-width', '0px');
    document.getElementById('metadata').appendChild(rect);
}
