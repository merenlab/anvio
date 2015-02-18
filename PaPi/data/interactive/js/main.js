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

var last_settings;

var search_column;
var search_results = [];
var highlight_backup = {};

var views = {};
var current_view = '';
var layer_order;

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
            url: '/data/views?timestamp=' + timestamp,
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
            url: '/data/default_view?timestamp=' + timestamp,
        }))
    .then(
        function (titleResponse, clusteringsResponse, viewsResponse, contigLengthsResponse, stateResponse, defaultViewResponse) 
        {
            var state = eval(stateResponse[0]);
            var hasState = !$.isEmptyObject(state);
            document.title = titleResponse[0];
            contig_lengths = eval(contigLengthsResponse[0]);

            /*
                Get metadata and create layers table
            */

            metadata = eval(defaultViewResponse[0]);
            parameter_count = metadata[0].length;

            // since we are painting parent layers odd-even, 
            // we should remove single parents (single means no parent)
            removeSingleParents(); // in utils.js

            if (hasState) {
                if (state.hasOwnProperty('layer-order')) {
                    layer_order = state['layer-order'];
                } else {
                    layer_order = Array.apply(null, Array(parameter_count-1)).map(function (_, i) {return i+1;}); // range(1, parameter_count)
                }

                if (state.hasOwnProperty('views'))
                    views = state['views'];

                if (state.hasOwnProperty('categorical_data_colors'))
                    categorical_data_colors = state['categorical_data_colors'];

                if (state.hasOwnProperty('stack_bar_colors'))
                    stack_bar_colors = state['stack_bar_colors'];

                if (state.hasOwnProperty('tree-type'))
                    $('#tree_type').val(state['tree-type']);
                if (state.hasOwnProperty('angle-min'))
                    $('#angle-min').val(state['angle-min']);
                if (state.hasOwnProperty('angle-max'))
                    $('#angle-max').val(state['angle-max']);
                if (state.hasOwnProperty('layer-margin'))
                    $('#layer-margin').val(state['layer-margin']);
            }
            else {
                layer_order = Array.apply(null, Array(parameter_count-1)).map(function (_, i) {return i+1;}); // range(1, parameter_count)
            }

            layer_types = {};

            //  Edit Attributes For Multiple Layers

            $('#select_layer').on('change', function() {
                var layer_name = $('#select_layer').val();

                // clean prior selections
                $('.layer_selectors').each(
                        function() {
                            this.checked = false;
                        }
                    );

                if(layer_name){ // if layer_name is empty, there is nothing to select, move on.
                    $('.titles').each(
                        function(){
                            if (this.title.indexOf(layer_name) > -1)
                            {
                                $('#' + 'select_this_' + getNumericPart(this.id)).attr('checked','checked');
                            }
                        }
                    );
                }
            });

            $('#picker_multiple').colpick({
                layout: 'hex',
                submit: 0,
                colorScheme: 'light',
                onChange: function(hsb, hex, rgb, el, bySetColor) {
                    $(el).css('background-color', '#' + hex);
                    $(el).attr('color', '#' + hex);
                    if (!bySetColor) $(el).val(hex);

                    $('.layer_selectors:checked').each(
                        function(){
                            $('#' + 'picker' + getNumericPart(this.id)).attr('color', '#' + hex);
                            $('#' + 'picker' + getNumericPart(this.id)).css('background-color', '#' + hex);
                        }
                    );
                }
            }).keyup(function() {
                    $(this).colpickSetColor(this.value);
            });

            $('#min_multiple').on('change', function(){
                var intend_value = $('#min_multiple').val();
                $('.layer_selectors:checked').each(
                    function(){
                        $('#' + 'min' + getNumericPart(this.id)).attr('value', intend_value);
                    }
                );
            });

            $('#max_multiple').on('change', function(){
                var intend_value = $('#max_multiple').val();
                $('.layer_selectors:checked').each(
                    function(){
                        $('#' + 'max' + getNumericPart(this.id)).attr('value', intend_value);
                    }
                );
            });

            $('#height_multiple').on('change', function(){
                var intend_value = $('#height_multiple').val();
                $('.layer_selectors:checked').each(
                    function(){
                        $('#' + 'height' + getNumericPart(this.id)).attr('value', intend_value);
                    }
                );
            });

            $('#normalization_multiple').on('change', function(){
                var intend_value = $('#normalization_multiple option:selected').val();
                $('.layer_selectors:checked').each(
                    function(){
                        var picker = $('#' + 'normalization' + getNumericPart(this.id));
                        $(picker).attr('value', intend_value);
                        clearMinMax(picker);
                    }
                );
            });

            $('#select_all').on('click', function() {
                if(this.checked) {
                    $('.layer_selectors').each(
                        function() {
                            this.checked = true;
                        }
                    );

                }else{
                    $('.layer_selectors').each(
                        function() {
                            this.checked = false;
                        }
                    );
                }
            });


            /* 
            //  Clusterings
            */
            var default_tree = (hasState && state.hasOwnProperty('order-by')) ? state['order-by'] : clusteringsResponse[0][0];
            var available_trees = clusteringsResponse[0][1];
            var available_trees_combo = getComboBoxContent(default_tree, available_trees);

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
            //  Views
            */
            var default_view = (hasState && state.hasOwnProperty('current-view')) ? state['current-view'] : viewsResponse[0][0];
            var available_views = viewsResponse[0][1];
            var available_views_combo = getComboBoxContent(default_view, available_views);

            $('#views_container').append(available_views_combo);

            $('#views_container').change(function() {

                $('#views_container').prop('disabled', false);
                $('#btn_draw_tree').prop('disabled', true);

                $.ajax({
                    type: 'GET',
                    cache: false,
                    url: '/data/view/' + $('#views_container').val() + '?timestamp=' + new Date().getTime(),
                    success: function(data) {
                        metadata = eval(data);
                        removeSingleParents(); // in utils.js
                        
                        $('#views_container').attr('disabled', false);
                        $('#btn_draw_tree').attr('disabled', false);

                        if (current_view != '') {
                            // backup current layer order and layers table to global views object
                            syncViews();
                        }
                        current_view = $('#views_container').val();

                        $("#tbody_layers").empty();

                        buildLayersTable(layer_order, views[current_view]);

                        // make layers table sortable
                        $("#tbody_layers").sortable({helper: fixHelperModified, handle: '.drag-icon', items: "> tr:not(:first)"}).disableSelection(); 

                        
                    }
                });
            });

            $('#views_container').trigger('change'); // load default view


            /*
            //  Add groups
            */
            if (hasState)
            {
                if (state.hasOwnProperty('SELECTED'))
                    SELECTED = state['SELECTED'];

                if (state.hasOwnProperty('group-counter'))
                    group_counter = state['group-counter'];

                if (state.hasOwnProperty('groups'))
                {
                    for (gid in state['groups'])
                    {
                        newGroup(gid, state['groups'][gid]);
                    }              
                }
            }
            else
            {
                newGroup();
            }

            initializeDialogs();

            // add metadata columns to search window
            for (var i=0; i < metadata[0].length; i++)
            {
                $('#searchLayerList').append(new Option(metadata[0][i],i));
            }
        } // response callback
    ); // promise

    document.body.addEventListener('click', function() {
        $('#control_contextmenu').hide();
    }, false);
}); // document ready

function syncViews() {
    views[current_view] = {};
    layer_order = new Array();
    $('#tbody_layers tr').each(
        function(index, layer) {
            var layer_id = $(layer).find('.input-height')[0].id.replace('height', '');

            layer_order.push(layer_id);

            views[current_view][layer_id] = {};
            views[current_view][layer_id]["normalization"] = $(layer).find('select').val();
            views[current_view][layer_id]["color"] = $(layer).find('.colorpicker').attr('color');
            views[current_view][layer_id]["height"] = $(layer).find('.input-height').val();
            views[current_view][layer_id]["min"] = {'value': $(layer).find('.input-min').val(), 'disabled': $(layer).find('.input-min').is(':disabled') }; 
            views[current_view][layer_id]["max"] = {'value': $(layer).find('.input-max').val(), 'disabled': $(layer).find('.input-max').is(':disabled') };
        }
    );    
}


function getComboBoxContent(default_item, available_items){
    var combo = '';
    var combo_item = '<option value="{val}"{sel}>{text}</option>';

    $.each(available_items, function(index, value) {
        if(index == default_item)
        {
            combo += combo_item
                        .replace('{val}', index)
                        .replace('{sel}', ' selected')
                        .replace('{text}', index);
        }
        else
        {
            combo += combo_item
                        .replace('{val}', index)
                        .replace('{sel}', '')
                        .replace('{text}', index);
        }
    });

    return combo;
}


// get numeric part from id
function getNumericPart(id){
    var $num = id.replace(/[^\d]+/, '');

    return $num;
}

function buildLayersTable(order, settings)
{
    for (var i = 0; i < order.length; i++) 
    {
        // common layer variables
        var layer_id = order[i];
        var layer_name = metadata[0][layer_id];
        var short_name = (layer_name.length > 10) ? layer_name.slice(0,10) + "..." : layer_name;

        var hasSettings = false;
        if (typeof settings !== 'undefined') {
            var layer_settings = settings[layer_id];
            var hasSettings = true;
        }

        //
        //  parent layer
        //
        if (layer_name == '__parent__')
        {
           layer_types[layer_id] = 0;

            if (hasSettings)
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
                '<td><input type="checkbox" id="select_this_{id}" class="layer_selectors"></input></td>' +
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

            if (hasSettings)
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
                '<td title="{name}" class="titles" id="title{id}">{short-name}</td>' +
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
                '<td><input type="checkbox" id="select_this_{id}" class="layer_selectors"></input></td>' +
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
        else if (metadata[1][layer_id] === null || !isNumber(metadata[1][layer_id]))
        { 
            layer_types[layer_id] = 2;

            if (hasSettings)
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
                '<td title="{name}" class="titles" id="title{id}">{short-name}</td>' +
                '<td>n/a</td>' +
                '<td>n/a</td>' +
                '<td><input class="input-height" type="text" size="3" id="height{id}" value="{height}"></input></td>' +
                '<td>n/a</td>' +
                '<td>n/a</td>' +
                '<td><input type="checkbox" id="select_this_{id}" class="layer_selectors"></input></td>' +
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

            if (hasSettings)
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
                '<td title="{name}" class="titles" id="title{id}">{short-name}</td>' +
                '<td><div id="picker{id}" class="colorpicker" color="{color}" style="background-color: {color}"></td>' +
                '<td>' +
                '    <select id="normalization{id}" onChange="clearMinMax(this);">' +
                '        <option value="none"{option-none}>none</option>' +
                '        <option value="sqrt"{option-sqrt}>Square root</option>' +
                '        <option value="log"{option-log}>Logarithm</option>' +
                '    </select>' +
                '</td>' +
                '<td><input class="input-height" type="text" size="3" id="height{id}" value="{height}"></input></td>' +
                '<td><input class="input-min" type="text" size="4" id="min{id}" value="{min}"{min-disabled}></input></td>' +
                '<td><input class="input-max" type="text" size="4" id="max{id}" value="{max}"{min-disabled}></input></td>' +
                '<td><input type="checkbox" id="select_this_{id}" class="layer_selectors"></input></td>' +
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
    }
}

//---------------------------------------------------------
//  ui callbacks
//---------------------------------------------------------

function serializeSettings() {
    var state = {};

    state['group-counter'] = group_counter;
    state['tree-type'] = $('#tree_type').val();
    state['order-by'] = $('#trees_container').val();
    state['current-view'] = $('#views_container').val();
    state['angle-min'] = $('#angle-min').val();
    state['angle-max'] = $('#angle-max').val();
    state['layer-margin'] = $('#layer-margin').val();
    state['edge-normalization'] = $('#edge_length_normalization').is(':checked');

    state['SELECTED'] = SELECTED;

    // sync views object and layers table
    syncViews();

    state['views'] = views;
    state['layer-order'] = layer_order;

    state['groups'] = {};
    $('#tbody_groups tr').each(
        function(index, group) {
            var gid = $(group).attr('group-id');

            state['groups'][gid] = {};
            state['groups'][gid]['name'] = $('#group_name_' + gid).val();
            state['groups'][gid]['color'] = $('#group_color_' + gid).attr('color');
            state['groups'][gid]['contig-length'] = $('#contig_length_' + gid).text();
            state['groups'][gid]['contig-count'] = $('#contig_count_' + gid).val();
            state['groups'][gid]['completeness'] = $('#completeness_' + gid).val();
            state['groups'][gid]['contamination'] = $('#contamination_' + gid).val();
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

            // last_settings used in export svg for layer information,
            // we didn't use "settings" sent to draw_tree because draw_tree updates layer's min&max
            // running serializeSettings() twice costs extra time but we can ignore it to keep code simple.
            last_settings = serializeSettings(); 

            $('#img_loading').hide();
            $('#btn_draw_tree').prop('disabled', false);

        }, 1); 
}


function getContigNames(gid) {
    var names = new Array();

    for (var j = 0; j < SELECTED[gid].length; j++) {
        if (label_to_node_map[SELECTED[gid][j]].IsLeaf()) {
            names.push(SELECTED[gid][j]);
        }
    }

    return names
}


function showContigNames(gid) {
    names = getContigNames(gid)
    messagePopupShow('Contig Names', names.join('\n'));
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
        var completeness = 0;
        var contamination = 0;

        SELECTED[group_counter] = [];
    }
    else
    {
        // we are adding groups from state file
        var name = groupState['name'];
        var color = groupState['color'];
        var contig_count = groupState['contig-count'];
        var contig_length = groupState['contig-length'];
        var completeness = groupState['completeness'];
        var contamination = groupState['contamination'];
    }

    var template = '<tr group-id="{id}" id="group_row_{id}">' +
                   '    <td><input type="radio" name="active_group" value="{id}" checked></td>' +
                   '    <td><div id="group_color_{id}" class="colorpicker" color="{color}" style="background-color: {color}"></td>' +
                   '    <td><input type="text" size="12" id="group_name_{id}" value="{name}"></td>' +
                   '    <td><input id="contig_count_{id}" type="button" value="{count}" title="Click for contig names" onClick="showContigNames({id});"></td> ' +
                   '    <td><span id="contig_length_{id}">{length}</span></td>' +
                   '    <td><span id="completeness_{id}">{completeness}%</span></td>' +
                   '    <td><input id="contamination_{id}" type="button" value="{contamination}%" title="Click for contaminants" onClick="showContaminants({id});"></td> ' +
                   '    <td><center><span class="delete-button ui-icon ui-icon-trash" alt="Delete this group" title="Delete this group" onClick="deleteGroup({id});"></span></center></td>' +
                   '</tr>';

    template = template.replace(new RegExp('{id}', 'g'), id)
                       .replace(new RegExp('{name}', 'g'), name)
                       .replace(new RegExp('{color}', 'g'), color)
                       .replace(new RegExp('{count}', 'g'), contig_count)
                       .replace(new RegExp('{completeness}', 'g'), completeness)
                       .replace(new RegExp('{contamination}', 'g'), contamination)
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
    // FIXME: this loop starts from 1 to group_counter, but things go bad when
    //        a group is removed :) so "group_counter" should not be an integer, but a
    //        list that contains all group ids, and ids should be removed from it when
    //        a group is removed. all thes for loops can then go through ids in that 
    //        list.
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

        // FIXME:
        // updateGroupWindow goes through every group (whether a group is recently updated or not). 
        // calling updateGroupWindow for every selection action is somewhat manageable
        // since what this function is doing is pretty starightforward at this point.
        // however, the following addition will add a lot to the computational budget
        // to it. I think this funciton should take a parameter (group_id) that needs
        // to be updated. when necessary, it would be called for all group_id's.
        split_names = getContigNames(gid);
        group_name = $('#group_name_' + gid).val();

        updateComplateness(gid);
    }
}

function updateComplateness(gid) {
    $.ajax({
        type: "POST",
        url: "/data/completeness",
        cache: false,
        data: {split_names: JSON.stringify(split_names), group_name: JSON.stringify(group_name)},
        success: function(data){
            completeness_info_dict = JSON.parse(data);
            console.log(completeness_info_dict);

            default_source = completeness_info_dict['default_source'];
            stats = completeness_info_dict['stats'];

            $('#completeness_' + gid).html(stats[default_source]['percent_complete']);
            $('#contamination_' + gid).val(stats[default_source]['percent_contamination']);
        },
    });
}

function exportSvg() {
    // check if tree parsed, which means there is a tree on the screen
    if ($.isEmptyObject(label_to_node_map)) 
        return;

    // draw group and layer legend to output svg
    var settings = serializeSettings();

    var groups_to_draw = new Array();
    for (var gid = 1; gid <= group_counter; gid++) {
        if(SELECTED[gid].length > 0) {
            groups_to_draw.push(settings['groups'][gid]);
        }
    }

    var left = 0 - total_radius - 400; // draw on the left top
    var top = 20 - total_radius;

    if (groups_to_draw.length > 0) {
        drawGroupLegend(groups_to_draw, top, left);
        top = top + 100 + (groups_to_draw.length + 2.5) * 20
    }

    // important,
    // we used current settings because we want current group information.
    // now we are going to use "last_settings" which updated by draw button.
    var settings = {};
    settings = last_settings; 

    drawLayerLegend(settings['views'][current_view], settings['layer-order'], top, left);

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
    $('#layer_legend').remove();
}

function searchContigs() 
{
    var svalue = $('#searchValue').val();

    if (svalue == "")
    {
        alert("Search value shouldn't be empty.");
        return;
    }
    var column = $('#searchLayerList').val();
    search_column = column;
    var operator = $('#searchOperator').val();
    
    if (operator < 6)
    {
        var operator_text = $('#searchOperator option:selected').text();

        // logical operator
        var _pre = "metadata[";
        var _post = "][" + column + "] " + operator_text + " \"" + svalue.trim() + "\"";

    }
    else if (operator == 6)
    {
        // contains
        var _pre = "metadata[";
        var _post = "][" + column + "].toString().indexOf(\"" + svalue + "\") != -1";
    }

    var _len = metadata.length;
    var _counter = 0;
    search_results = [];

    $('#search_result_message').html("Searching...");

    for (var row=1; row < _len; row++)
    {
        if (eval(_pre + row + _post)){
            search_results.push(row);
            _counter++;
        }
    }
    $('#search_result_message').html(_counter + " contigs found.");
}

function showSearchResult() {
    var msg = "Line\t\tContig Name\t\t" + metadata[0][search_column] + "\n";

    var _len = search_results.length;
    for (var i=0; i < _len; i++)
    {
        msg = msg + search_results[i] + "\t\t" + metadata[search_results[i]][0] + "\t\t" + metadata[search_results[i]][search_column] + "\n";
    }
    messagePopupShow('Search Results ('+_len+" items)", msg);
}

function highlightResult() {
    var HIGHLIGHT_COLOR= "#FFC000";

    // check if tree exists
    if ($.isEmptyObject(label_to_node_map)) {
        alert('Draw tree first.');
        return;
    }

    // clear previous highlight
    clearHighlight();

    var _len = search_results.length;
    for (var i=0; i < _len; i++) {
        var _contig_name = metadata[search_results[i]][0];
        var _id = label_to_node_map[_contig_name].id;

        var _path_background = document.getElementsByClassName('path_' + _id + '_background');
        for (var _i=0; _i < _path_background.length; _i++) {
            _path_background[_i].style['fill'] = HIGHLIGHT_COLOR;
            _path_background[_i].style['fill-opacity'] = '0.1';  
        }
        var _path_outer_ring = document.getElementsByClassName('path_' + _id + '_outer_ring');
        for (var _i=0; _i < _path_outer_ring.length; _i++) {
            if (_i==0) {
                highlight_backup[_contig_name] = _path_outer_ring[_i].style['fill'];
            }
            _path_outer_ring[_i].style['fill'] = HIGHLIGHT_COLOR;
        }
    }
}

function clearHighlight() {
    for (_contig_name in highlight_backup)
    {
        var _id = label_to_node_map[_contig_name].id;
        var _color = highlight_backup[_contig_name];
        var _opacity = (_color == "#FFFFFF" || _color == "") ? '0.0' : '0.1';

        var _path_background = document.getElementsByClassName('path_' + _id + '_background');
        for (var _i=0; _i < _path_background.length; _i++) {
            _path_background[_i].style['fill'] = _color;
            _path_background[_i].style['fill-opacity'] = _opacity;  
        }
        var _path_outer_ring = document.getElementsByClassName('path_' + _id + '_outer_ring');
        for (var _i=0; _i < _path_outer_ring.length; _i++) {
            _path_outer_ring[_i].style['fill'] = _color;
        }
    }

    highlight_backup = {};

}

function appendResult() {
    // check if tree exists
    if ($.isEmptyObject(label_to_node_map)) {
        alert('Draw tree first.');
        return;
    }

    var group_id = getGroupId();

    if (group_id === 'undefined')
        return;

    clearHighlight();

    var _len = search_results.length;
    for (var i=0; i < _len; i++) {
        _contig_name = metadata[search_results[i]][0];
        if (SELECTED[group_id].indexOf(_contig_name) == -1)
            SELECTED[group_id].push(_contig_name);

        for (var gid = 1; gid <= group_counter; gid++) {
            // don't remove nodes from current group
            if (gid == group_id)
                continue;

            var pos = SELECTED[gid].indexOf(_contig_name);
            if (pos > -1) {
                SELECTED[gid].splice(pos, 1);
            }
        }
    }

    updateGroupWindow();
    redrawGroupColors(group_id);
}

function removeResult() {
    // check if tree exists
    if ($.isEmptyObject(label_to_node_map)) {
        alert('Draw tree first.');
        return;
    }

    var group_id = getGroupId();

    if (group_id === 'undefined')
        return;

    var _len = search_results.length;
    for (var i=0; i < _len; i++) {
        _contig_name = metadata[search_results[i]][0];
        var _id = label_to_node_map[_contig_name].id;

        var pos = SELECTED[group_id].indexOf(_contig_name);
        if (pos > -1) {
            SELECTED[group_id].splice(pos, 1);
        }

        var _path_background = document.getElementsByClassName('path_' + _id + '_background');
        for (var _i=0; _i < _path_background.length; _i++) {
            _path_background[_i].style['fill'] = '#FFFFFF';
            _path_background[_i].style['fill-opacity'] = '0.0';  
        }
        var _path_outer_ring = document.getElementsByClassName('path_' + _id + '_outer_ring');
        for (var _i=0; _i < _path_outer_ring.length; _i++) {
            _path_outer_ring[_i].style['fill'] = '#FFFFFF';
        }
    }

    updateGroupWindow();
}
