sampleOrganizationResponse = [JSON.parse('{"num_reads":{"order":"","newick":"((DAY_19:3.8927e-12,(DAY_15A:1.17876e-12,DAY_17A:1.17876e-12)Int16:3.8927e-12)Int18:1.58727e-11,((DAY_15B:1.14795e-12,(DAY_17B:4.17546e-14,DAY_23:4.17546e-14)Int12:1.14795e-12)Int15:5.55138e-12,(DAY_22B:1.20141e-12,((DAY_18:3.08426e-14,DAY_22A:3.08426e-14)Int11:7.76164e-13,(DAY_16:3.99589e-13,DAY_24:3.99589e-13)Int13:7.76164e-13)Int14:1.20141e-12)Int17:5.55138e-12)Int19:1.58727e-11);"},"even_odd":{"order":"","newick":"((DAY_22B:0.000332086,((DAY_18:8.54156e-06,DAY_22A:8.54156e-06)Int11:0.00021462,(DAY_16:0.000110463,DAY_24:0.000110463)Int13:0.00021462)Int14:0.000332086)Int15:0.0376259,((DAY_15B:0.000506205,(DAY_17B:1.8431e-05,DAY_23:1.8431e-05)Int12:0.000506205)Int16:0.00441781,(DAY_19:0.00170499,(DAY_15A:0.00051519,DAY_17A:0.00051519)Int17:0.00170499)Int18:0.00441781)Int19:0.0376259);"},"basic":{"order":"DAY_15A,DAY_15B,DAY_16,DAY_17A,DAY_17B,DAY_18,DAY_19,DAY_22A,DAY_22B,DAY_23,DAY_24","newick":""}, "mini_test": {"order":"s204_6M,s204_7M,s204_9M", "newick":""}}')];
sampleMetadataResponse = [JSON.parse('{"DAY_17A":{"days_after_birth":"17","num_reads":"17","percent_mapped_reads":"17","days":"17","bar3":"17","sections":"17"},"DAY_17B":{"days_after_birth":"17","num_reads":"17","percent_mapped_reads":"17","days":"17","bar3":"17","sections":"17"},"DAY_18":{"days_after_birth":"17","num_reads":"17","percent_mapped_reads":"17","days":"17","bar3":"17","sections":"17"},"DAY_19":{"days_after_birth":"17","num_reads":"17","percent_mapped_reads":"17","days":"17","bar3":"17","sections":"17"},"DAY_15B":{"days_after_birth":"17","num_reads":"17","percent_mapped_reads":"17","days":"17","bar3":"17","sections":"17"},"DAY_22A":{"days_after_birth":"17","num_reads":"17","percent_mapped_reads":"17","days":"17","bar3":"17","sections":"17"},"DAY_16":{"days_after_birth":"17","num_reads":"17","percent_mapped_reads":"17","days":"17","bar3":"17","sections":"17"},"DAY_15A":{"days_after_birth":"17","num_reads":"17","percent_mapped_reads":"17","days":"17","bar3":"17","sections":"17"},"DAY_23":{"days_after_birth":"17","num_reads":"17","percent_mapped_reads":"17","days":"17","bar3":"17","sections":"17"},"DAY_22B":{"days_after_birth":"17","num_reads":"17","percent_mapped_reads":"17","days":"17","bar3":"17","sections":"17"},"DAY_24":{"days_after_birth":"17","num_reads":"17","percent_mapped_reads":"17","days":"17","bar3":"17","sections":"17"}}')];

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

function buildMetadataTable() {
    var layers = sampleMetadataResponse[0][Object.keys(sampleMetadataResponse[0])[0]];
    
    for (layer in layers)
    {
        var layer_name = layer;
        var short_name = (layer_name.length > 10) ? layer_name.slice(0,10) + "..." : layer_name;


        var norm = "none";
        var min    = 0;
        var max    = 0;
        var min_disabled = true;
        var max_disabled = true;
        var height = 180;
        var color  = '#000000';
        var margin = '15';
        var color_start = "#FFFFFF";
        var type = "bar";

        var template = '<tr metadata-layer-name="{name}">' +
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
            '        <option value="sqrt"{option-sqrt}>Square root</option>' +
            '        <option value="log"{option-log}>Logarithm</option>' +
            '    </select>' +
            '</td>' +
            '<td><input class="input-height" type="text" size="3" value="50"></input></td>' +
            '<td><input class="input-margin" type="text" size="3" value="15"></input></td>' +
            '<td><input class="input-min" type="text" size="4" value="{min}"{min-disabled}></input></td>' +
            '<td><input class="input-max" type="text" size="4" value="{max}"{min-disabled}></input></td>' +
            '<td><input type="checkbox" class="layer_selectors"></input></td>' +
            '</tr>';

        template = template.replace(new RegExp('{name}', 'g'), layer_name)
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

function drawMetadataLayers() {

}
