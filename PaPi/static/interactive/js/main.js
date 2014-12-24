//--------------------------------------------------------------------------------------------------
//  Globals
//--------------------------------------------------------------------------------------------------

var VIEWER_WIDTH = window.innerWidth || document.documentElement.clientWidth || document.getElementsByTagName('body')[0].clientWidth;
var VIEWER_HEIGHT = window.innerHeight || document.documentElement.clientHeight || document.getElementsByTagName('body')[0].clientHeight;

var LINE_COLOR='#888888';

var scale = 0;

var id_to_node_map = new Array();
var label_to_node_map = {};

var angle_per_leaf;
var height_per_leaf;
var tree_type;

var total_radius = 0;

var SELECTED = new Array();

var newick;

var metadata;
var contig_lengths;
var parameter_count;

var group_counter = 0; // for id
var group_count = 0;

var layer_types;

var categorical_data_colors = {};
var stack_bar_colors = {};

var context_menu_target_id = 0;

var metadata_title = {};
var metadata_dict;

//---------------------------------------------------------
//  Init
//---------------------------------------------------------

$(document).ready(function() {

    $('.dialogs').hide();

    var timestamp = new Date().getTime(); 

    $.when(    
        $.ajax({
            type: 'GET',
            cache: false,
            url: '/data/title?timestamp=' + timestamp,
        }),
        $.ajax({
            type: 'GET',
            cache: false,
            url: '/data/clusterings?timestamp=' + timestamp,
        }),
        $.ajax({
            type: 'GET',
            cache: false,
            url: '/data/contig_lengths?timestamp=' + timestamp,
        }),
        $.ajax({
            type: 'GET',
            cache: false,
            url: '/data/state?timestamp=' + timestamp,
        }),
        $.ajax({
            type: 'GET',
            cache: false,
            url: '/data/meta?timestamp=' + timestamp,
        }))
    .then(
        function (titleResponse, clusteringsResponse, contigLengthsResponse, stateResponse, metaResponse) 
        {
            var state = eval(stateResponse[0]);
            var hasState = !$.isEmptyObject(state);
            document.title = titleResponse[0];
            contig_lengths = eval(contigLengthsResponse[0]);

            /*
                Get metadata and create layers table
            */
            metadata = eval(metaResponse[0]);
            parameter_count = metadata[0].length;

            // since we are painting parent layers odd-even, 
            // we should remove single parents (single means no parent)
            removeSingleParents(); // in utils.js

            var layer_order;

            if (hasState) {
                layer_order = state['layer-order'];

                categorical_data_colors = state['categorical_data_colors'];
                stack_bar_colors = state['stack_bar_colors'];

                $('#tree_type').val(state['tree-type']);
                $('#angle-min').val(state['angle-min']);
                $('#angle-max').val(state['angle-max']);
                $('#layer-margin').val(state['layer-margin']);
                $('#')
            }
            else {
                // range(1, prameter_count), we skipped column 0 because its not a layer, its name column.
                layer_order = Array.apply(null, Array(parameter_count-1)).map(function (_, i) {return i+1;}); 
            }
            //
            //  add layers to table
            //
            layer_types = {};

            for (var i = 0; i < layer_order.length; i++) 
            {
                // common layer variables
                var layer_id = layer_order[i];
                var layer_name = metadata[0][layer_id];
                var short_name = (layer_name.length > 10) ? layer_name.slice(0,10) + "..." : layer_name;

                if (hasState)
                    var layer_settings = state['layers'][layer_id];

                //
                //  parent layer
                //
                if (layer_name == '__parent__')
                {
                   layer_types[layer_id] = 0;

                    if (hasState)
                        var height = layer_settings['height'];
                    else
                        var height = '50';

                    var template = '<tr>' +
                        '<td><img src="images/drag.gif" /></td>' +
                        '<td>Parent</td>' +
                        '<td>n/a</td>' +
                        '<td>n/a</td>' +
                        '<td><input class="input-height" type="text" size="3" id="height{id}" value="{height}"></input></td>' +
                        '<td>n/a</td>' +
                        '<td>n/a</td>' +
                        '</tr>';

                    template = template.replace(new RegExp('{id}', 'g'), layer_id)
                                       .replace(new RegExp('{height}', 'g'), height);

                    $('#tbody_layers').prepend(template);
                }
                //
                // stack bar layer
                //
                else if (layer_name.indexOf(';') > -1) 
                {
                    layer_types[layer_id] = 1;

                    if (hasState)
                    {
                        var norm   = layer_settings['normalization'];
                        var height = layer_settings['height'];
                    }
                    else
                    {
                        var norm   = 'log';
                        var height = '30';  

                        // pick random color for stack bar items
                        stack_bar_colors[layer_id] = new Array();
                        for (var j=0; j < layer_name.split(";").length; j++)
                        {
                            stack_bar_colors[layer_id].push(randomColor());
                        }              
                    }

                    var template = '<tr>' +
                        '<td><img class="drag-icon" src="images/drag.gif" /></td>' +
                        '<td title="{name}">{short-name}</td>' +
                        '<td>n/a</td>' +
                        '<td>' +
                        '    <select id="normalization{id}" onChange="clearMinMax(this);">' +
                        '        <option value="none"{option-none}>none</option>' +
                        '        <option value="sqrt"{option-sqrt}>Square root</option>' +
                        '        <option value="log"{option-log}>Logarithm</option>' +
                        '    </select>' +
                        '</td>' +
                        '<td><input class="input-height" type="text" size="3" id="height{id}" value="{height}"></input></td>' +
                        '<td>n/a</td>' +
                        '<td>n/a</td>' +
                        '</tr>';

                    template = template.replace(new RegExp('{id}', 'g'), layer_id)
                                       .replace(new RegExp('{name}', 'g'), layer_name)
                                       .replace(new RegExp('{short-name}', 'g'), short_name)
                                       .replace(new RegExp('{option-' + norm + '}', 'g'), ' selected')
                                       .replace(new RegExp('{option-([a-z]*)}', 'g'), '')
                                       .replace(new RegExp('{height}', 'g'), height);

                    $('#tbody_layers').append(template);
                }
                //
                // categorical layer
                //
                else if (metadata[1][layer_id] === '' || !isNumber(metadata[1][layer_id]))
                { 
                    layer_types[layer_id] = 2;

                    if (hasState)
                    {
                        var height = layer_settings['height'];
                    }
                    else
                    {
                        var height = 30;

                        categorical_data_colors[layer_id] = {};
                    }
                    
                    var template = '<tr>' +
                        '<td><img class="drag-icon" src="images/drag.gif" /></td>' +
                        '<td title="{name}">{short-name}</td>' +
                        '<td>n/a</td>' +
                        '<td>n/a</td>' +
                        '<td><input class="input-height" type="text" size="3" id="height{id}" value="{height}"></input></td>' +
                        '<td>n/a</td>' +
                        '<td>n/a</td>' +
                        '</tr>';

                    template = template.replace(new RegExp('{id}', 'g'), layer_id)
                                       .replace(new RegExp('{name}', 'g'), layer_name)
                                       .replace(new RegExp('{short-name}', 'g'), short_name)
                                       .replace(new RegExp('{height}', 'g'), height);

                    $('#tbody_layers').append(template);
                } 
                //
                // numerical layer
                //
                else
                {
                    layer_types[layer_id] = 3;

                    if (hasState)
                    {
                        var height = layer_settings['height'];
                        var norm   = layer_settings['normalization'];
                        var color  = layer_settings['color'];
                        var min    = layer_settings['min']['value'];
                        var max    = layer_settings['max']['value'];
                        var min_disabled = layer_settings['min']['disabled'];
                        var max_disabled = layer_settings['max']['disabled'];
                    }
                    else
                    {
                        var height = getNamedLayerDefaults(layer_name, 'height', '180');
                        var norm   = getNamedLayerDefaults(layer_name, 'norm', 'log');
                        var color  = getNamedLayerDefaults(layer_name, 'color', '#000000');
                        var min    = 0;
                        var max    = 0;
                        var min_disabled = true;
                        var max_disabled = true;
                    }

                    var template = '<tr>' +
                        '<td><img class="drag-icon" src="images/drag.gif" /></td>' +
                        '<td title="{name}">{short-name}</td>' +
                        '<td><div id="picker{id}" class="colorpicker" color="{color}" style="background-color: {color}"></td>' +
                        '<td>' +
                        '    <select id="normalization{id}" onChange="clearMinMax(this);">' +
                        '        <option value="none"{option-none}>none</option>' +
                        '        <option value="sqrt"{option-sqrt}>Square root</option>' +
                        '        <option value="log"{option-log}>Logarithm</option>' +
                        '    </select>' +
                        '</td>' +
                        '<td><input class="input-height" type="text" size="3" id="height{id}" value="150"></input></td>' +
                        '<td><input class="input-min" type="text" size="4" id="min{id}" value="{min}"{min-disabled}></input></td>' +
                        '<td><input class="input-max" type="text" size="4" id="max{id}" value="{max}"{min-disabled}></input></td>' +
                        '</tr>';

                    template = template.replace(new RegExp('{id}', 'g'), layer_id)
                                       .replace(new RegExp('{name}', 'g'), layer_name)
                                       .replace(new RegExp('{short-name}', 'g'), short_name)
                                       .replace(new RegExp('{option-' + norm + '}', 'g'), ' selected')
                                       .replace(new RegExp('{option-([a-z]*)}', 'g'), '')
                                       .replace(new RegExp('{color}', 'g'), color)
                                       .replace(new RegExp('{height}', 'g'), height)
                                       .replace(new RegExp('{min}', 'g'), min)
                                       .replace(new RegExp('{max}', 'g'), max)
                                       .replace(new RegExp('{min-disabled}', 'g'), (min_disabled) ? ' disabled': '')
                                       .replace(new RegExp('{max-disabled}', 'g'), (max_disabled) ? ' disabled': '');

                    $('#tbody_layers').append(template);
                }

                $('#picker'+ layer_id).colpick({
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
                
            } // layer loop

            // make layers table sortable
            $("#tbody_layers").sortable({helper: fixHelperModified, handle: '.drag-icon', items: "> tr:not(:first)"}).disableSelection(); 

            /* 
            //  Clusterings
            */
            var default_tree = (hasState) ? state['order-by'] : clusteringsResponse[0][0];
            var available_trees = clusteringsResponse[0][1];

            var available_trees_combo = '';
            var available_trees_combo_item = '<option value="{val}"{sel}>{text}</option>';
            
            $.each(available_trees, function(index, value) {
                if(index == default_tree)
                {
                    available_trees_combo += available_trees_combo_item
                                .replace('{val}', index)
                                .replace('{sel}', ' selected')
                                .replace('{text}', index);
                }
                else
                {
                    available_trees_combo += available_trees_combo_item
                                .replace('{val}', index)
                                .replace('{sel}', '')
                                .replace('{text}', index);
                }
            }); 
            
            $('#trees_container').append(available_trees_combo);

            $('#trees_container').change(function() {

                $('#trees_container').prop('disabled', true);
                $('#btn_draw_tree').prop('disabled', true);

                $.ajax({
                    type: 'GET',
                    cache: false,
                    url: '/tree/' + $('#trees_container').val() + '?timestamp=' + new Date().getTime(),
                    success: function(data) {
                        newick = data;
                        $('#trees_container').attr('disabled', false);
                        $('#btn_draw_tree').attr('disabled', false); 
                    }
                });
            });

            $('#trees_container').trigger('change'); // load default newick tree

            /*
            //  Add groups
            */
            if (hasState)
            {
                SELECTED = state['SELECTED'];
                group_counter = state['group-counter'];

                for (gid in state['groups'])
                {
                    newGroup(gid, state['groups'][gid]);
                }
            }
            else
            {
                newGroup();
            }

            initializeDialogs();
        } // response callback
    ); // promise

    document.body.addEventListener('click', function() {
        $('#control_contextmenu').hide();
    }, false);
}); // document ready


//---------------------------------------------------------
//  ui callbacks
//---------------------------------------------------------

function serializeSettings() {
    var state = {};

    state['group-counter'] = group_counter;
    state['tree-type'] = $('#tree_type').val();
    state['order-by'] = $('#trees_container').val();
    state['angle-min'] = $('#angle-min').val();
    state['angle-max'] = $('#angle-max').val();
    state['layer-margin'] = $('#layer-margin').val();
    state['edge-normalization'] = $('#edge_length_normalization').is(':checked');

    state['layers'] = {};
    state['layer-order'] = new Array();
    $('#tbody_layers tr').each(
        function(index, layer) {
            var layer_id = $(layer).find('.input-height')[0].id.replace('height', '');

            state['layer-order'].push(layer_id);

            state['layers'][layer_id] = {};
            state['layers'][layer_id]["normalization"] = $(layer).find('select').val();
            state['layers'][layer_id]["color"] = $(layer).find('.colorpicker').attr('color');
            state['layers'][layer_id]["height"] = $(layer).find('.input-height').val();
            state['layers'][layer_id]["min"] = {'value': $(layer).find('.input-min').val(), 'disabled': $(layer).find('.input-min').is(':disabled') }; 
            state['layers'][layer_id]["max"] = {'value': $(layer).find('.input-max').val(), 'disabled': $(layer).find('.input-max').is(':disabled') };
        }
    );

    state['SELECTED'] = SELECTED;

    state['groups'] = {};
    $('#tbody_groups tr').each(
        function(index, group) {
            var gid = $(group).attr('group-id');

            state['groups'][gid] = {};
            state['groups'][gid]['name'] = $('#group_name_' + gid).val();
            state['groups'][gid]['color'] = $('#group_color_' + gid).attr('color');
            state['groups'][gid]['contig-length'] = $('#contig_length_' + gid).text();
            state['groups'][gid]['contig-count'] = $('#contig_count_' + gid).val();

        }
    );

    state['categorical_data_colors'] = categorical_data_colors;
    state['stack_bar_colors'] = stack_bar_colors;

    return state;
}

function saveCurrentState() {
    $.post("/save_state", {
        state: JSON.stringify(serializeSettings(), null, 4),
    });
}

function drawTree() {

    var settings = serializeSettings();

    tree_type = settings['tree-type'];

    $('#img_loading').show();
    $('#btn_draw_tree').prop('disabled', true);

    setTimeout(function () 
        { 
            draw_tree(settings); // call treelib.js where the magic happens

            $('#img_loading').hide();
            $('#btn_draw_tree').prop('disabled', false);

        }, 1); 
}

function showContigNames(gid) {
    var names = new Array();

    for (var j = 0; j < SELECTED[gid].length; j++) {
        if (label_to_node_map[SELECTED[gid][j]].IsLeaf()) {
            names.push(SELECTED[gid][j]);
        }
    }

    if (names.length == 0)
        return;

    messagePopupShow('Contig Names', names.join('\n'));
    /*
    $("#contig_name_dialog").dialog("option", "title", "Contig Names");
    $('#contig_names').val(names.join("\n"));
    $('#contig_name_dialog').dialog('open');
    $('#contig_names').click(); // focus & select all*/
}

function newGroup(id, groupState) {

    group_count++;

    if (typeof id === 'undefined')
    {
        group_counter++;
        var id = group_counter;
        var name = "Group_" + id;
        var color = '#000000';
        var contig_count = 0;
        var contig_length = 0;

        SELECTED[group_counter] = [];
    }
    else
    {
        // we are adding groups from state file
        var name = groupState['name'];
        var color = groupState['color'];
        var contig_count = groupState['contig-count'];
        var contig_length = groupState['contig-length'];
    }

    var template = '<tr group-id="{id}" id="group_row_{id}">' +
                   '    <td><input type="radio" name="active_group" value="{id}" checked></td>' +
                   '    <td><div id="group_color_{id}" class="colorpicker" color="{color}" style="background-color: {color}"></td>' +
                   '    <td><input type="text" size="12" id="group_name_{id}" value="{name}"></td>' +
                   '    <td><input id="contig_count_{id}" type="button" value="{count}" title="Click for contig names" onClick="showContigNames({id});"></td> ' +
                   '    <td><span id="contig_length_{id}">{length}</span></td>' +
                   '    <td><span class="delete-button ui-icon ui-icon-trash" alt="Delete this group" title="Delete this group" onClick="deleteGroup({id});"></span></td>' +
                   '</tr>';

    template = template.replace(new RegExp('{id}', 'g'), id)
                       .replace(new RegExp('{name}', 'g'), name)
                       .replace(new RegExp('{color}', 'g'), color)
                       .replace(new RegExp('{count}', 'g'), contig_count)
                       .replace(new RegExp('{length}', 'g'), contig_length);

    $('#tbody_groups').append(template);

    $('#group_color_' + id).colpick({
        layout: 'hex',
        submit: 0,
        colorScheme: 'dark',
        onChange: function(hsb, hex, rgb, el, bySetColor) {
            $(el).css('background-color', '#' + hex);
            $(el).attr('color', '#' + hex);

            if (!bySetColor) $(el).val(hex);
        },
        onHide: function() {
            redrawGroupColors(id);
        }
    }).keyup(function() {
        $(this).colpickSetColor(this.value);
    });
}

function deleteGroup(id) {
    if (confirm('Are you sure?')) {

        $('#group_row_' + id).remove();
        $('#tbody_groups input[type=radio]').last().prop('checked', true);
        group_count--;

        for (var i = 0; i < SELECTED[id].length; i++) {
            var node_id = label_to_node_map[SELECTED[id][i]].id;
            $("#line" + node_id).css('stroke-width', '1');
            $("#arc" + node_id).css('stroke-width', '1');
            $("#line" + node_id).css('stroke', LINE_COLOR);
            $("#arc" + node_id).css('stroke', LINE_COLOR);

            if (label_to_node_map[SELECTED[id][i]].IsLeaf())
            {
                $('.path_' + node_id + "_background").css({'fill': '#FFFFFF', 'fill-opacity': '0.0'});
                $('.path_' + node_id + "_outer_ring").css('fill', '#FFFFFF');
            }
        }

        SELECTED[id] = [];

        if (group_count==0)
        {
            newGroup();
        }
    }
}

function submitGroups() {
    if ($.isEmptyObject(label_to_node_map)) {
        alert('You should draw tree before submit.');
        return;
    }

    var output = {};
    var msg_group_count = 0;
    var msg_contig_count = 0;

    for (var gid = 1; gid <= group_counter; gid++) {
        if (SELECTED[gid].length > 0) {
            msg_group_count++;
            var group_name = $('#group_name_' + gid).val();

            output[group_name] = new Array();
            for (var i = 0; i < SELECTED[gid].length; i++) {
                if (label_to_node_map[SELECTED[gid][i]].IsLeaf()) {
                    output[group_name].push(SELECTED[gid][i]);
                    msg_contig_count++;
                }
            }
        }
    }

    if (!confirm('You\'ve selected ' + msg_contig_count + ' contigs in ' + msg_group_count + ' group. You won\'t able to select more contigs after submit. Do you want to continue?')) {
        return;
    }

    $('#tbody_groups input[type=radio]').prop('checked', false).attr("disabled", true);
    $('#submit-groups').attr("disabled", true);
    $('#btn_new_group').attr("disabled", true);
    window.deleteGroup = function() {};

    $.post("/submit", {
        groups: JSON.stringify(output),
        svg: null
    });
}

function updateGroupWindow() { 
    // count contigs and update group labels
    for (var gid = 1; gid <= group_counter; gid++) {
        var contigs = 0;
        var length_sum = 0;

        for (var j = 0; j < SELECTED[gid].length; j++) {
            if (label_to_node_map[SELECTED[gid][j]].IsLeaf())
            {
                contigs++;
                length_sum += parseInt(contig_lengths[SELECTED[gid][j]]);
            }
        }

        $('#contig_count_' + gid).val(contigs);
        $('#contig_length_' + gid).html(readableNumber(length_sum));
    }
}

function exportSvg() {
    // check if tree parsed, which means there is a tree on the screen
    if ($.isEmptyObject(label_to_node_map)) 
        return;

    // draw group list to output svg
    drawGroupLegend();

    // move group highlights to new svg groups
    for (var gid = 1; gid <= group_counter; gid++) {

        createGroup('tree_group', 'selected_group_' + gid);

        for (var j = 0; j < SELECTED[gid].length; j++) {
            if (label_to_node_map[SELECTED[gid][j]].IsLeaf()) {
                $('.path_' + label_to_node_map[SELECTED[gid][j]].id + "_background").detach().appendTo('#selected_group_' + gid);
                $('.path_' + label_to_node_map[SELECTED[gid][j]].id + "_outer_ring").detach().appendTo('#selected_group_' + gid);
            }
        }
    }

    // remove ungrouped backgrounds.
    if (tree_type == 'circlephylogram')
    {
        var detached_paths = $('#tree_group > path').detach();        
    }
    else
    {
        var detached_paths = $('#tree_group > rect').detach();   
    }

    svgCrowbar();

    // add removed ungrouped backgrounds back
    $(detached_paths).appendTo('#tree_group');

    $('#group_legend').remove();
}
