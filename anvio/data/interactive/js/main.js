//--------------------------------------------------------------------------------------------------
//  Globals
//--------------------------------------------------------------------------------------------------

var VIEWER_WIDTH;
var VIEWER_HEIGHT;
var dragging = false;
var shiftPressed = false;

var LINE_COLOR='#888888';

var scale = 0;

var id_to_node_map = new Array();
var label_to_node_map = {};
var order_to_node_map = {};
var leaf_count;

var angle_per_leaf;
var height_per_leaf;
var tree_type;
var margin;
var order_counter;

var total_radius = 0;
var layer_boundaries;

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
var empty_tooltip = "";

var last_settings;

var search_column;
var search_results = [];
var highlight_backup = {};

var views = {};
var current_view = '';
var layer_order;

var completeness_dict = {};

//---------------------------------------------------------
//  Init
//---------------------------------------------------------

$(document).ready(function() {

    $('.dialogs').hide();
    $('.dialogs2').hide();

    $('#tree_type').change(function() {
        if ($('#tree_type').val()=='circlephylogram') 
        {
            $('.phylogram_settings').hide();
            $('.circlephylogram_settings').show();
        }
        else
        {
            $('.phylogram_settings').show();
            $('.circlephylogram_settings').hide();
        }
    });

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
                if (state.hasOwnProperty('tree-height'))
                    $('#tree_height').val(state['tree-height']);
                if (state.hasOwnProperty('tree-width'))
                    $('#tree_width').val(state['tree-width']);
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
            newGroup();

            initializeDialogs();

            // add metadata columns to search window
            for (var i=0; i < metadata[0].length; i++)
            {
                $('#searchLayerList').append(new Option(metadata[0][i],i));
            }
        } // response callback
    ); // promise

    // initialize colorpicker for search result highlight color.
    $('#picker_highlight').colpick({
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

    document.body.addEventListener('click', function() {
        $('#control_contextmenu').hide();
    }, false);

    $(window).keydown(function(evt) {
      if (evt.which == 16) { // shift
        shiftPressed = true;
      }
    }).keyup(function(evt) {
      if (evt.which == 16) { // shift
        shiftPressed = false;
      }
    });
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
                if (!(layer_id in stack_bar_colors))
                {
                    stack_bar_colors[layer_id] = new Array();
                    for (var j=0; j < layer_name.split(";").length; j++)
                    {
                        stack_bar_colors[layer_id].push(randomColor());
                    } 
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

                if (!(layer_id in categorical_data_colors))
                {
                    categorical_data_colors[layer_id] = {};
                }
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
    state['tree-height'] = $('#tree_height').val();
    state['tree-width'] = $('#tree_width').val();
    state['layer-margin'] = $('#layer-margin').val();
    state['edge-normalization'] = $('#edge_length_normalization').is(':checked');

    state['SELECTED'] = SELECTED;

    // sync views object and layers table
    syncViews();

    state['views'] = views;
    state['layer-order'] = layer_order;

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
    // get current client size
    VIEWER_WIDTH = window.innerWidth || document.documentElement.clientWidth || document.getElementsByTagName('body')[0].clientWidth;
    VIEWER_HEIGHT = window.innerHeight || document.documentElement.clientHeight || document.getElementsByTagName('body')[0].clientHeight;

    var settings = serializeSettings();
    tree_type = settings['tree-type'];

    $('#img_loading').show();
    $('#draw_delta_time').html('');
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
        var from_state = false;
        var id = group_counter;
        var name = "Group_" + id;
        var color = '#000000';
        var contig_count = 0;
        var contig_length = 0;
        var completeness = '---';
        var contamination = '---';

        SELECTED[group_counter] = [];
    }
    else
    {
        // we are adding groups from collection
        var from_state = true;
        var name = groupState['name'];
        var color = groupState['color'];
        var contig_count = 0;
        var contig_length = 0;
        var completeness = "---";
        var contamination = "---";
    }

    var template = '<tr group-id="{id}" id="group_row_{id}">' +
                   '    <td><input type="radio" name="active_group" value="{id}" checked></td>' +
                   '    <td><div id="group_color_{id}" class="colorpicker" color="{color}" style="background-color: {color}"></td>' +
                   '    <td><input type="text" size="12" id="group_name_{id}" value="{name}"></td>' +
                   '    <td><input id="contig_count_{id}" type="button" value="{count}" title="Click for contig names" onClick="showContigNames({id});"></td> ' +
                   '    <td><span id="contig_length_{id}">{length}</span></td>' +
                   '    <td><input id="completeness_{id}" type="button" value="{completeness}" title="Click for completeness table" onClick="showCompleteness({id});"></td> ' +
                   '    <td><input id="contamination_{id}" type="button" value="{contamination}" title="Click for contaminants" onClick="showContaminants({id});"></td> ' +
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

    if(!from_state){
        $('#completeness_' + id).attr("disabled", true);
        $('#contamination_' + id).attr("disabled", true);
    }

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
            redrawGroups();
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
        }

        SELECTED[id] = [];
        delete completeness_dict[id];

        if (group_count==0)
        {
            newGroup();
        }

        redrawGroups();
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

function updateGroupWindow(group_list) {
    if (typeof group_list === 'undefined')
    {
        var group_list = [];
        $('#tbody_groups tr').each(
        function(index, group) {
            group_list.push(parseInt($(group).attr('group-id')));
        });
    }

    for (var _i = 0; _i < group_list.length; _i++) {
        var gid = group_list[_i];
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

            stats = completeness_info_dict['stats'];
            refs = completeness_info_dict['refs'];

            completeness_dict[gid] = completeness_info_dict;

            showCompleteness(gid, true);
            showContaminants(gid, true);

            // it is sad that one liner with a list comprehension takes this much code in JS:
            sum_completeness = 0.0;
            sum_contamination = 0.0;
            for(source in stats){
                sum_completeness += stats[source]['percent_complete'];
                sum_contamination += stats[source]['percent_contamination'];
            }
            average_completeness = sum_completeness / Object.keys(stats).length;
            average_contamination = sum_contamination / Object.keys(stats).length;

            $('#completeness_' + gid).val(average_completeness.toFixed(1) + '%');
            $('#contamination_' + gid).val(average_contamination.toFixed(1) + '%');

            $('#completeness_' + gid).attr("disabled", false);
            $('#contamination_' + gid).attr("disabled", false);
        },
    });
}

function showCompleteness(gid, updateOnly) {
    if (!completeness_dict.hasOwnProperty(gid))
        return;

    if ($('#completeness_content_' + gid).length)
    {
        $('#completeness_content_' + gid).html(buildCompletenessTable(completeness_dict[gid]));
        return;
    }

    if (!updateOnly)
    {
        $('<div> \
           <div id="completeness_content_' + gid + '">' + buildCompletenessTable(completeness_dict[gid]) + '</div> \
           </div>').dialog({
                resizable: false,
                collapseEnabled: false,
                width: 'auto',
                title: 'Completeness Info',
                position: {
                    my: "center",
                    at: "center",
                    of: window
                },
                close: function(ev, ui) {
                    $(this).remove();
                }});
    }
}

function showContaminants(gid, updateOnly) {
    if (!completeness_dict.hasOwnProperty(gid))
        return;

    if ($('#contaminants_content_' + gid).length)
    {
        $('#contaminants_content_' + gid).html(buildContaminantsTable(completeness_dict[gid]));
        return;
    }

    if (!updateOnly)
    {
        $('<div> \
           <div id="contaminants_content_' + gid + '">' + buildContaminantsTable(completeness_dict[gid]) + '</div> \
           </div>').dialog({
                resizable: false,
                collapseEnabled: false,
                width: 'auto',
                title: 'Contaminants Info',
                position: {
                    my: "center",
                    at: "center",
                    of: window
                },
                close: function(ev, ui) {
                    $(this).remove();
                }});
    }
}

function buildContaminantsTable(info_dict) {
    var stats = info_dict['stats'];

    var output = '<div style="max-height: 400px; width: 480px; margin; 5px; overflow: auto;"><table style="width: 100%;">';
    var header = '<tr>';
    var body = '<tr>';

    for(var source in stats) {
        header += '<td><b>' + source + ' (' + Object.keys(stats[source]['contaminants']).length + ')</br></td>';

        var contaminants_html = '';

        for (var contaminant in stats[source]['contaminants']) {
            var title = '';
            var order_array = '[';
            for (var i = 0; i < stats[source]['contaminants'][contaminant].length; i++)
            {
                var contig = stats[source]['contaminants'][contaminant][i];

                for (var j = 0; j < contig.length; j++) {
                    // splits
                    title += contig[j] + '&#xA;'; // it's a hack for using newline in tooltips
                    order_array += label_to_node_map[contig[j]].order + ', ';
                }
            }
            order_array += ']';

            contaminants_html += '<span style="cursor:pointer;" \
                                    title="' + title + '" \
                                    onclick="redrawGroups(' + order_array + ');"> \
                                    ' + contaminant + ' (' + stats[source]['contaminants'][contaminant].length + ') \
                                  </span><br />';
        }

        body += '<td valign="top">' + contaminants_html + '</td>';
    }

    output += header + '</tr><tr><td colspan="' + Object.keys(stats).length + '"><hr /></td>' + body + '</tr></table></div>';
    return output;
}

function buildCompletenessTable(info_dict) {
    var stats = info_dict['stats'];
    var refs = info_dict['refs'];

    var output = '<table style="margin: 5px;">';

    header_template = '<td style="font-variant: small-caps; font-size: 125%; padding: 0px 10px;">{TXT}</td>';
    percentages_template = '<td style="text-align: center; font-weight: bold;">{VAL}%</td>';
    source_template = '<td style="text-align: left;"><a href="{REF}" target="_blank" style="text-decoration:none; outline: 0;">{VAL}</a></td>';

    headers = ['source', 'completeness', 'contamination'];
    entry_keys = ['percent_complete', 'percent_contamination'];

    // header line:
    output += '<tr>';
    headers.forEach(function(header) {
        output += header_template.replace('{TXT}', header);
    });
    output += '</tr><tr><td colspan="' + headers.length + '"><hr /></td>';

    for (var source in stats) {
        output += '<tr>';

        output += source_template.replace('{VAL}', source).replace('{REF}', refs[source]);

        entry_keys.forEach(function(key) {
            output += percentages_template.replace('{VAL}', stats[source][key].toFixed(1));
        })
        output += '</tr>';
    }

    output += '</tr><tr><td colspan="' + headers.length + '"><hr /></td></table>';
    return output;
}

function exportSvg() {
    // check if tree parsed, which means there is a tree on the screen
    if ($.isEmptyObject(label_to_node_map)) 
        return;

    // draw group and layer legend to output svg
    var settings = serializeSettings();

    var groups_to_draw = new Array();
    $('#tbody_groups tr').each(
        function(index, group) {
            var gid = $(group).attr('group-id');
            groups_to_draw.push({
                'name': $('#group_name_' + gid).val(),
                'color': $('#group_color_' + gid).attr('color'),
                'contig-length': $('#contig_length_' + gid).html(),
                'contig-count': $('#contig_count_' + gid).val(),
            });
        }
    );

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

    svgCrowbar();

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
        if (metadata[row][column]==null)
            continue;

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
    // check if tree exists
    if ($.isEmptyObject(label_to_node_map)) {
        alert('Draw tree first.');
        return;
    }

    var order_list = new Array();

    for (var i=0; i < search_results.length; i++) {
        var _contig_name = metadata[search_results[i]][0];
        var _order = label_to_node_map[_contig_name].order;

        order_list.push(_order);
    }

    redrawGroups(order_list); 
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

    var groups_to_update = [];
    var _len = search_results.length;
    for (var i=0; i < _len; i++) {
        _contig_name = metadata[search_results[i]][0];
        if (SELECTED[group_id].indexOf(_contig_name) == -1) {
            SELECTED[group_id].push(_contig_name);

            if (groups_to_update.indexOf(group_id) == -1)
                groups_to_update.push(group_id);
        }

        for (var gid = 1; gid <= group_counter; gid++) {
            // don't remove nodes from current group
            if (gid == group_id)
                continue;

            var pos = SELECTED[gid].indexOf(_contig_name);
            if (pos > -1) {
                SELECTED[gid].splice(pos, 1);

                if (groups_to_update.indexOf(gid) == -1)
                    groups_to_update.push(gid);
            }
        }
    }

    updateGroupWindow(groups_to_update);
    redrawGroups();
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

    var groups_to_update = [];
    var _len = search_results.length;
    for (var i=0; i < _len; i++) {
        _contig_name = metadata[search_results[i]][0];
        var _id = label_to_node_map[_contig_name].id;

        var pos = SELECTED[group_id].indexOf(_contig_name);
        if (pos > -1) {
            SELECTED[group_id].splice(pos, 1);
            
            if (groups_to_update.indexOf(group_id) == -1)
                groups_to_update.push(group_id);
        }
    }

    updateGroupWindow(groups_to_update);
    redrawGroups();
}

function showStoreCollectionWindow() {
    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/collections?timestamp=' + new Date().getTime(),
        success: function(data) {
            $('#storeCollection_list').empty();

            for (source in data) {
                var read_only = data[source]["read_only"];

                if (read_only) {
                    var _disabled = ' disabled="true"';
                    var _name = source + ' (read only)';
                }
                else
                {
                    var _disabled = '';
                    var _name = source;
                }

                $('#storeCollection_list').append('<option value="' + source + '"' + _disabled + '>' + _name + '</option>');
            }

            $('#storeCollectionWindow').dialog('open');
        }
    });
}

function storeCollection() {
    var collection_name = $('#storeCollection_name').val();

    data = {};
    colors = {};

    $('#tbody_groups tr').each(
        function(index, group) {
            var gid = $(group).attr('group-id');
            var gname = $('#group_name_' + gid).val();

            colors[gname] = $('#group_color_' + gid).attr('color');
            data[gname] = new Array();

            for (var i=0; i < SELECTED[gid].length; i++)
            {
                if (label_to_node_map[SELECTED[gid][i]].IsLeaf())
                {
                    data[gname].push(SELECTED[gid][i]);
                }
            }
        }
    );

    $.post("/store_collections", {
        source: collection_name,
        data: JSON.stringify(data, null, 4),
        colors: JSON.stringify(colors, null, 4),
    },
    function(server_response, status){
          alert("Server: " + server_response + "\n\n(status: " + status + ")");
    });

    $('#storeCollectionWindow').dialog('close');    
}

function showLoadCollectionWindow() {
    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/collections?timestamp=' + new Date().getTime(),
        success: function(data) {
            $('#loadCollection_list').empty();

            for (source in data) {
                var read_only = data[source]["read_only"];

                if (read_only) {
                    var _name = source + ' (read only)';
                }
                else
                {
                    var _name = source;
                }

                $('#loadCollection_list').append('<option value="' + source + '">' + _name + '</option>');
            }

            $('#loadCollectionWindow').dialog('open');
        }
    });
}

function showCollectionDetails() {
    var cname = $('#loadCollection_list>option:selected').val();
    
    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/collections?timestamp=' + new Date().getTime(),
        success: function(data) {
            var tbl = "<table>" +
                "<tr><td>Read Only: </td><td>" + data[cname]['read_only'] + "</td></tr>" + 
                "<tr><td>Source DB Path: </td><td>" + data[cname]['source_db_path'] + "</td></tr>" +
                "<tr><td>Source DB Version: </td><td>" + data[cname]['source_db_version'] + "</td></tr>" +
                "<tr><td colspan='2'><hr></td></tr>" +
                "<tr><td>Number of Splits: </td><td>" + data[cname]['num_splits'] + "</td></tr>" +
                "<tr><td>Number of Clusters: </td><td>" + data[cname]['num_clusters'] + "</td></tr>" +
                "</table>";

            $('#loadCollection_details').html(tbl);
        }
    });
}

function loadCollection() {
    if ($.isEmptyObject(label_to_node_map)) {
        alert('You should draw tree before load collection.');
        return;
    }

    var collection = $('#loadCollection_list').val();
    if (collection === null)
        return;

    if (!confirm("You will lost current groups, please be sure you stored current groups. Do you want to continue?"))
        return;

    var coldata 
    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/collection/' + collection + '?timestamp=' + new Date().getTime(),
        success: function(data) {
            // empty group window
            $('#tbody_groups').empty();
            SELECTED = new Array();
            group_count = 0;
            group_counter = 0;

            // count groups.
            var skip_empty = false;
            var num_groups = Object.keys(data['data']).length;

            if (num_groups > 16) {
                skip_empty = confirm("This collection has " + num_groups + " groups, Do you want to skip empty groups to improve user interface performance?");
            }

            // load new groups
            var gid=0;
            for (group in data['data'])
            {
                // collection may be contain unknown splits/contigs, we should clear them.
                var contigs = new Array();

                for (index in data['data'][group])
                {
                    if (typeof contig_lengths[data['data'][group][index]] !== 'undefined') {
                        contigs.push(data['data'][group][index]);
                    }
                    
                }

                if (!(skip_empty && contigs.length == 0))
                {
                    gid++;
                    group_counter++;
                    SELECTED[gid] = contigs;

                    var _color =  (data['colors'][group]) ? data['colors'][group] : randomColor();

                    newGroup(gid, {'name': group, 'color': _color});
                }
            }

            rebuildIntersections();
            updateGroupWindow();
            redrawGroups();
            $('#loadCollectionWindow').dialog('close');
        }
    });
}
