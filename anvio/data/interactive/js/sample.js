sampleOrganizationResponse = [JSON.parse('{"num_reads":{"order":"","newick":"((DAY_19:3.8927e-12,(DAY_15A:1.17876e-12,DAY_17A:1.17876e-12)Int16:3.8927e-12)Int18:1.58727e-11,((DAY_15B:1.14795e-12,(DAY_17B:4.17546e-14,DAY_23:4.17546e-14)Int12:1.14795e-12)Int15:5.55138e-12,(DAY_22B:1.20141e-12,((DAY_18:3.08426e-14,DAY_22A:3.08426e-14)Int11:7.76164e-13,(DAY_16:3.99589e-13,DAY_24:3.99589e-13)Int13:7.76164e-13)Int14:1.20141e-12)Int17:5.55138e-12)Int19:1.58727e-11);"},"even_odd":{"order":"","newick":"((DAY_22B:0.000332086,((DAY_18:8.54156e-06,DAY_22A:8.54156e-06)Int11:0.00021462,(DAY_16:0.000110463,DAY_24:0.000110463)Int13:0.00021462)Int14:0.000332086)Int15:0.0376259,((DAY_15B:0.000506205,(DAY_17B:1.8431e-05,DAY_23:1.8431e-05)Int12:0.000506205)Int16:0.00441781,(DAY_19:0.00170499,(DAY_15A:0.00051519,DAY_17A:0.00051519)Int17:0.00170499)Int18:0.00441781)Int19:0.0376259);"},"basic":{"order":"DAY_15A,DAY_15B,DAY_16,DAY_17A,DAY_17B,DAY_18,DAY_19,DAY_22A,DAY_22B,DAY_23,DAY_24","newick":""}, "mini_test": {"order":"s204_6M,s204_7M,s204_9M", "newick":""}}')];
sampleMetadataResponse = [JSON.parse('{"DAY_17A":{"days_after_birth":"17","num_reads":"2409083","percent_mapped_reads":"84.4","days":"even","sections":"A"},"DAY_17B":{"days_after_birth":"17","num_reads":"14467205","percent_mapped_reads":"92.96","days":"even","sections":"A"},"DAY_18":{"days_after_birth":"18","num_reads":"11157806","percent_mapped_reads":"96.65","days":"odd","sections":"A"},"DAY_19":{"days_after_birth":"19","num_reads":"12605189","percent_mapped_reads":"88.5","days":"even","sections":"A"},"DAY_15B":{"days_after_birth":"15","num_reads":"8924116","percent_mapped_reads":"91.76","days":"even","sections":"A"},"DAY_22A":{"days_after_birth":"22","num_reads":"12585100","percent_mapped_reads":"96.69","days":"odd","sections":"B"},"DAY_16":{"days_after_birth":"16","num_reads":"10505839","percent_mapped_reads":"95.71","days":"odd","sections":"A"},"DAY_15A":{"days_after_birth":"15","num_reads":"5639445","percent_mapped_reads":"85.58","days":"even","sections":"A"},"DAY_23":{"days_after_birth":"23","num_reads":"6318127","percent_mapped_reads":"92.91","days":"even","sections":"B"},"DAY_22B":{"days_after_birth":"22","num_reads":"7386372","percent_mapped_reads":"95.11","days":"odd","sections":"B"},"DAY_24":{"days_after_birth":"24","num_reads":"10756936","percent_mapped_reads":"96.22","days":"odd","sections":"B"}}')];
var metadata = sampleMetadataResponse[0];
var metadata_categorical_colors = {};

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
        if (this.value == 'none') 
            return;

        var organization = sampleOrganizationResponse[0][this.value];
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

function buildMetadataTable(state) {
    var first_sample = Object.keys(metadata)[0];
    var layers = metadata[first_sample]; // get layers from first sample's metadata
    $('#tbody_metadata').empty();

    for (layer in layers)
    {
        var layer_name = layer;
        var short_name = (layer_name.length > 10) ? layer_name.slice(0,10) + "..." : layer_name;

        if (isNumber(metadata[first_sample][layer]))
        {
            var data_type = "numeric";
            var norm = "none";
            var min    = 0;
            var max    = 0;
            var min_disabled = true;
            var max_disabled = true;
            var height = 500;
            var color  = '#000000';
            var margin = 15;
            var color_start = "#FFFFFF";
            var type = "bar";

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
            var height = 80;
            var margin = 15;

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

        $('#tbody_metadata').prepend(template);   
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

    var _metadata = metadata; // keep original

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

                if (typeof metadata_categorical_colors[layer][value] === 'undefined')
                {
                    metadata_categorical_colors[layer][value] = randomColor();
                }

                var color = metadata_categorical_colors[layer][value];
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

            }
        }

        backgrounds_done = true;
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
