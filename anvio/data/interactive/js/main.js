/**
 * Javascript library for anvi'o interactive interface
 *
 *  Author: Özcan Esen <ozcanesen@gmail.com>
 *  Credits: A. Murat Eren, Doğan Can Kilment
 *  Copyright 2015, The anvio Project
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

//--------------------------------------------------------------------------------------------------
//  Globals
//--------------------------------------------------------------------------------------------------
var VERSION = '0.2.0';

var VIEWER_WIDTH;
var VIEWER_HEIGHT;
var dragging = false;
var windowsMinimized = false;

var zoomBox = {};
var drawing_zoom = false;

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

var bin_counter = 0; // for id
var bin_count = 0;

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
var layers = {};
var current_view = '';
var layer_order;

var completeness_dict = {};

var sort_column;
var sort_order;

var readOnly = false;
var refineMode = false;
var bin_prefix;

var current_state_name = "";
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
            url: '/data/default_view?timestamp=' + timestamp,
        }),
        $.ajax({
            type: 'GET',
            cache: false,
            url: '/data/mode?timestamp=' + timestamp,
        }),
        $.ajax({
            type: 'GET',
            cache: false,
            url: '/data/read_only?timestamp=' + timestamp,
        }),
        $.ajax({
            type: 'GET',
            cache: false,
            url: '/data/bin_prefix?timestamp=' + timestamp,
        }))
    .then(
        function (titleResponse, clusteringsResponse, viewsResponse, contigLengthsResponse, defaultViewResponse, modeResponse, readOnlyResponse, prefixResponse) 
        {
            if (modeResponse[0] == 'refine')
            {
                refineMode = true;
                $('.full-mode').hide();
                $('.refine-mode').show();
            }

            if (readOnlyResponse[0] == true)
            {
                readOnly = true;
                alert('It seems that this is a read-only instance, therefore the database-writing functions will be inaccessible.');

                $('#btn_save_current_state').prop('disabled', 'disabled');
                $('#btn_store_collection').prop('disabled', 'disabled');
                $('#generate-summary').prop('disabled', 'disabled');
                $('#btn_store_refined_bins').prop('disabled', 'disabled');
            }

            bin_prefix = prefixResponse[0];

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

            layer_order = Array.apply(null, Array(parameter_count-1)).map(function (_, i) {return i+1;}); // range(1, parameter_count)
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

            $('#picker_start_multiple').colpick({
                layout: 'hex',
                submit: 0,
                colorScheme: 'light',
                onChange: function(hsb, hex, rgb, el, bySetColor) {
                    $(el).css('background-color', '#' + hex);
                    $(el).attr('color', '#' + hex);
                    if (!bySetColor) $(el).val(hex);

                    $('.layer_selectors:checked').each(
                        function(){
                            $('#' + 'picker_start' + getNumericPart(this.id)).attr('color', '#' + hex);
                            $('#' + 'picker_start' + getNumericPart(this.id)).css('background-color', '#' + hex);
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

            $('#margin_multiple').on('change', function(){
                var intend_value = $('#margin_multiple').val();
                $('.layer_selectors:checked').each(
                    function(){
                        $('#' + 'margin' + getNumericPart(this.id)).attr('value', intend_value);
                    }
                );
            });

            $('#type_multiple').on('change', function(){
                var intend_value = $('#type_multiple').val();
                $('.layer_selectors:checked').each(
                    function(){
                        $('#' + 'type' + getNumericPart(this.id)).attr('value', intend_value).trigger('change');
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
            var default_tree = clusteringsResponse[0][0];
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
            var default_view = viewsResponse[0][0];
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
            //  Add bins
            */
            newBin();

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

    $(window).keyup(function(evt) {
        if (evt.altKey && (evt.which == 84)) // alt+t for toggle windows
        {
            if (windowsMinimized)
            {
                $("button[title='collapse']").trigger('click');
            }
            else
            {
                $("button[title='restore']").trigger('click');
            }
            windowsMinimized = !windowsMinimized;
        }
    });
}); // document ready

function syncViews() {
    views[current_view] = {};
    layer_order = new Array();
    $('#tbody_layers tr').each(
        function(index, layer) {
            var layer_id = $(layer).find('.input-height')[0].id.replace('height', '');
            layers[layer_id] = {};
            layer_order.push(layer_id);

            views[current_view][layer_id] = {};
            views[current_view][layer_id]["normalization"] = $(layer).find('.normalization').val();
            views[current_view][layer_id]["min"] = {'value': $(layer).find('.input-min').val(), 'disabled': $(layer).find('.input-min').is(':disabled') }; 
            views[current_view][layer_id]["max"] = {'value': $(layer).find('.input-max').val(), 'disabled': $(layer).find('.input-max').is(':disabled') };

            layers[layer_id]["color"] = $(layer).find('.colorpicker:last').attr('color');
            layers[layer_id]["height"] = $(layer).find('.input-height').val();
            layers[layer_id]["margin"] = $(layer).find('.input-margin').val();
            layers[layer_id]["type"] = $(layer).find('.type').val();
            layers[layer_id]["color-start"] = $(layer).find('.colorpicker:first').attr('color');


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

        var hasViewSettings = false;
        if (typeof settings !== 'undefined') {
            var view_settings = settings[layer_id];
            var hasViewSettings = true;
        }

        var hasLayerSettings = false;
        if (typeof layers[layer_id] !== 'undefined')
        {
            var layer_settings = layers[layer_id];
            hasLayerSettings = true;
        }

        //
        //  parent layer
        //
        if (layer_name == '__parent__')
        {
           layer_types[layer_id] = 0;

            if (hasLayerSettings) 
            {
                var height = layer_settings['height'];
                var margin = layer_settings['margin'];
            }
            else 
            {
                var height = '50';
                var margin = '15';
            }

            var template = '<tr>' +
                '<td><img src="images/drag.gif" /></td>' +
                '<td>Parent</td>' +
                '<td>n/a</td>' +
                '<td>n/a</td>' +
                '<td>n/a</td>' +
                '<td><input class="input-height" type="text" size="3" id="height{id}" value="{height}"></input></td>' +
                '<td class="column-margin"><input class="input-margin" type="text" size="3" id="margin{id}" value="{margin}"></input></td>' +
                '<td>n/a</td>' +
                '<td>n/a</td>' +
                '<td><input type="checkbox" id="select_this_{id}" class="layer_selectors"></input></td>' +
                '</tr>';

            template = template.replace(new RegExp('{id}', 'g'), layer_id)
                               .replace(new RegExp('{height}', 'g'), height)
                               .replace(new RegExp('{margin}', 'g'), margin);

            $('#tbody_layers').prepend(template);
        }
        //
        // stack bar layer
        //
        else if (layer_name.indexOf(';') > -1) 
        {
            layer_types[layer_id] = 1;

            if (hasLayerSettings)
            {
                var height = layer_settings['height'];
                var margin = layer_settings['margin'];
            }
            else
            {
                var height = '30';
                var margin = '15';

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

            if (hasViewSettings)
            {
                var norm = view_settings['normalization'];
            }
            else
            {
                var norm = 'log';
            }

            var template = '<tr>' +
                '<td><img class="drag-icon" src="images/drag.gif" /></td>' +
                '<td title="{name}" class="titles" id="title{id}">{short-name}</td>' +
                '<td>n/a</td>' +
                '<td>n/a</td>' +
                '<td>' +
                '    <select id="normalization{id}" onChange="clearMinMax(this);" class="normalization">' +
                '        <option value="none"{option-none}>none</option>' +
                '        <option value="sqrt"{option-sqrt}>Square root</option>' +
                '        <option value="log"{option-log}>Logarithm</option>' +
                '    </select>' +
                '</td>' +
                '<td><input class="input-height" type="text" size="3" id="height{id}" value="{height}"></input></td>' +
                '<td class="column-margin"><input class="input-margin" type="text" size="3" id="margin{id}" value="{margin}"></input></td>' +
                '<td>n/a</td>' +
                '<td>n/a</td>' +
                '<td><input type="checkbox" id="select_this_{id}" class="layer_selectors"></input></td>' +
                '</tr>';

            template = template.replace(new RegExp('{id}', 'g'), layer_id)
                               .replace(new RegExp('{name}', 'g'), layer_name)
                               .replace(new RegExp('{short-name}', 'g'), short_name)
                               .replace(new RegExp('{option-' + norm + '}', 'g'), ' selected')
                               .replace(new RegExp('{option-([a-z]*)}', 'g'), '')
                               .replace(new RegExp('{height}', 'g'), height)
                               .replace(new RegExp('{margin}', 'g'), margin);

            $('#tbody_layers').append(template);
        }
        //
        // categorical layer
        //
        else if (metadata[1][layer_id] === null || !isNumber(metadata[1][layer_id]))
        { 
            layer_types[layer_id] = 2;

            if (hasLayerSettings)
            {
                var height = layer_settings['height'];
                var margin = layer_settings['margin'];
            }
            else
            {
                var height = '30';
                var margin = '15';

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
                '<td>n/a</td>' +
                '<td><input class="input-height" type="text" size="3" id="height{id}" value="{height}"></input></td>' +
                '<td class="column-margin"><input class="input-margin" type="text" size="3" id="margin{id}" value="{margin}"></input></td>' +
                '<td>n/a</td>' +
                '<td>n/a</td>' +
                '<td><input type="checkbox" id="select_this_{id}" class="layer_selectors"></input></td>' +
                '</tr>';

            template = template.replace(new RegExp('{id}', 'g'), layer_id)
                               .replace(new RegExp('{name}', 'g'), layer_name)
                               .replace(new RegExp('{short-name}', 'g'), short_name)
                               .replace(new RegExp('{height}', 'g'), height)
                               .replace(new RegExp('{margin}', 'g'), margin);

            $('#tbody_layers').append(template);
        } 
        //
        // numerical layer
        //
        else
        {
            layer_types[layer_id] = 3;

            if (hasViewSettings)
            {
                var norm   = view_settings['normalization'];
                var min    = view_settings['min']['value'];
                var max    = view_settings['max']['value'];
                var min_disabled = view_settings['min']['disabled'];
                var max_disabled = view_settings['max']['disabled'];
            }
            else
            {
                var norm   = getNamedLayerDefaults(layer_name, 'norm', 'log');
                var min    = 0;
                var max    = 0;
                var min_disabled = true;
                var max_disabled = true;
            }

            if (hasLayerSettings)
            {
                var height = layer_settings['height'];
                var color  = layer_settings['color'];
                var margin = layer_settings['margin'];
                var color_start = layer_settings['color-start'];
                var type = layer_settings['type'];
            }
            else
            {
                var height = getNamedLayerDefaults(layer_name, 'height', '180');
                var color  = getNamedLayerDefaults(layer_name, 'color', '#000000');
                var margin = '15';
                var color_start = "#FFFFFF";
                var type = "bar";
            }

            /* Some ad-hoc manipulation of special hmmx_ split hmm layers */ 
            if (layer_name.substring(0, 5) == "hmmx_"){
                var height = '30';
                var norm   = 'none';
                var color  = '#882222';
            }

            /* Some ad-hoc manipulation of special hmms_ single hmm layers */ 
            if (layer_name.substring(0, 5) == "hmms_"){
                var height = '150';
                var norm   = 'sqrt';
                var color  = '#882222';
            }

            var template = '<tr>' +
                '<td><img class="drag-icon" src="images/drag.gif" /></td>' +
                '<td title="{name}" class="titles" id="title{id}">{short-name}</td>' +
                '<td><div id="picker_start{id}" class="colorpicker picker_start" color="{color-start}" style="background-color: {color-start}; {color-start-hide}"></div><div id="picker{id}" class="colorpicker" color="{color}" style="background-color: {color}"></div></td>' +
                '<td style="width: 50px;">' +
                '    <select id="type{id}" style="width: 50px;" class="type" onChange="if(this.value==\'intensity\') {  $(\'#picker_start\' + getNumericPart(this.id)).css(\'visibility\', \'visible\'); } else { $(\'#picker_start\' + getNumericPart(this.id)).css(\'visibility\', \'hidden\'); }  ">' +
                '        <option value="bar"{option-type-bar}>Bar</option>' +
                '        <option value="intensity"{option-type-intensity}>Intensity</option>' +
                '    </select>' +
                '</td>' +
                '<td>' +
                '    <select id="normalization{id}" onChange="clearMinMax(this);" class="normalization">' +
                '        <option value="none"{option-none}>none</option>' +
                '        <option value="sqrt"{option-sqrt}>Square root</option>' +
                '        <option value="log"{option-log}>Logarithm</option>' +
                '    </select>' +
                '</td>' +
                '<td><input class="input-height" type="text" size="3" id="height{id}" value="{height}"></input></td>' +
                '<td class="column-margin"><input class="input-margin" type="text" size="3" id="margin{id}" value="{margin}"></input></td>' +
                '<td><input class="input-min" type="text" size="4" id="min{id}" value="{min}"{min-disabled}></input></td>' +
                '<td><input class="input-max" type="text" size="4" id="max{id}" value="{max}"{min-disabled}></input></td>' +
                '<td><input type="checkbox" id="select_this_{id}" class="layer_selectors"></input></td>' +
                '</tr>';

            template = template.replace(new RegExp('{id}', 'g'), layer_id)
                               .replace(new RegExp('{name}', 'g'), layer_name)
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
                               .replace(new RegExp('{margin}', 'g'), margin)


            $('#tbody_layers').append(template);
        }

        if($('#custom_layer_margin').is(':checked'))
        {
            $('.column-margin').show();
        }
        else
        {
            $('.column-margin').hide();
        }

        $('#picker'+ layer_id + ', #picker_start' + layer_id).colpick({
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
function getLayerName(layer_id)
{
    return metadata[0][layer_id];
}

function getLayerId(layer_name) 
{
    for (var i=0; i < parameter_count; i++)
    {
        if (layer_name == metadata[0][i])
            return i;
    }
    return -1;
}

function serializeSettings(use_layer_names) {

    if (typeof use_layer_names === 'undefined')
        use_layer_names = false;

    var state = {};
    state['version'] = VERSION;
    state['tree-type'] = $('#tree_type').val();
    state['order-by'] = $('#trees_container').val();
    state['current-view'] = $('#views_container').val();
    state['angle-min'] = $('#angle-min').val();
    state['angle-max'] = $('#angle-max').val();
    state['tree-radius'] = $('#tree-radius').val();
    state['tree-height'] = $('#tree_height').val();
    state['tree-width'] = $('#tree_width').val();
    state['layer-margin'] = $('#layer-margin').val();
    state['outer-ring-height'] = $('#outer-ring-height').val();
    state['edge-normalization'] = $('#edge_length_normalization').is(':checked');
    state['custom-layer-margin'] = $('#custom_layer_margin').is(':checked');
    state['show-grid-for-bins'] = $('#show_grid_for_bins').is(':checked');

    // sync views object and layers table
    syncViews();

    if (use_layer_names)
    {
        // save state file with layer name instead of id.

        state['layer-order'] = [];
        for (var i=0; i < layer_order.length; i++)
        {
            state['layer-order'].push(getLayerName(layer_order[i]));
        }

        state['categorical_data_colors'] = {};
        for (var key in categorical_data_colors)
        {
            state['categorical_data_colors'][getLayerName(key)] = categorical_data_colors[key];
        }

        state['stack_bar_colors'] = {};
        for (var key in stack_bar_colors)
        {
            state['stack_bar_colors'][getLayerName(key)] = stack_bar_colors[key];
        }

        state['layers'] = {};
        for (var key in layers)
        {
            state['layers'][getLayerName(key)] = layers[key];
        }

        state['views'] = {};
        for (var view_key in views)
        {
            state['views'][view_key] = {};
            for (var key in views[view_key])
            {
                state['views'][view_key][getLayerName(key)] = views[view_key][key];
            }
        }


    }
    else
    {
        state['views'] = views;
        state['layer-order'] = layer_order;
        state['layers'] = layers;

        state['categorical_data_colors'] = categorical_data_colors;
        state['stack_bar_colors'] = stack_bar_colors;
    }

    return state;
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
    $('#bins').dialog('open');

    setTimeout(function () 
        { 
            draw_tree(settings); // call treelib.js where the magic happens

            // last_settings used in export svg for layer information,
            // we didn't use "settings" sent to draw_tree because draw_tree updates layer's min&max
            // running serializeSettings() twice costs extra time but we can ignore it to keep code simple.
            last_settings = serializeSettings();

            redrawBins();

            $('#img_loading').hide();
            $('#btn_draw_tree').prop('disabled', false);

            if (settings['tree-radius'] == 0)
            {
                $('#tree-radius-container').show();
                $('#tree-radius').val(Math.max(VIEWER_HEIGHT, VIEWER_WIDTH));
            }
        }, 1); 
}


function getContigNames(bin_id) {
    var names = new Array();

    for (var j = 0; j < SELECTED[bin_id].length; j++) {
        if (label_to_node_map[SELECTED[bin_id][j]].IsLeaf()) {
            names.push(SELECTED[bin_id][j]);
        }
    }

    return names
}


function showContigNames(bin_id) {
    names = getContigNames(bin_id)
    messagePopupShow('Contig Names', names.join('<br />'));
}

function newBin(id, binState) {

    bin_count++;

    if (typeof id === 'undefined')
    {
        bin_counter++;
        var from_state = false;
        var id = bin_counter;
        var name = bin_prefix + id;
        var color = '#000000';
        var contig_count = 0;
        var contig_length = 0;
        var completeness = '---';
        var contamination = '---';

        SELECTED[bin_counter] = [];
    }
    else
    {
        // we are adding bins from collection
        var from_state = true;
        var name = binState['name'];
        var color = binState['color'];
        var contig_count = 0;
        var contig_length = 0;
        var completeness = "---";
        var contamination = "---";
    }

    var template = '<tr bin-id="{id}" id="bin_row_{id}">' +
                   '    <td><input type="radio" name="active_bin" value="{id}" checked></td>' +
                   '    <td><div id="bin_color_{id}" class="colorpicker" color="{color}" style="background-color: {color}"></td>' +
                   '    <td><input type="text" size="12" id="bin_name_{id}" value="{name}"></td>' +
                   '    <td><input id="contig_count_{id}" type="button" value="{count}" title="Click for contig names" onClick="showContigNames({id});"></td> ' +
                   '    <td><span id="contig_length_{id}">{length}</span></td>' +
                   '    <td><input id="completeness_{id}" type="button" value="{completeness}" title="Click for completeness table" onClick="showCompleteness({id});"></td> ' +
                   '    <td><input id="contamination_{id}" type="button" value="{contamination}" title="Click for contaminants" onClick="showContaminants({id});"></td> ' +
                   '    <td><center><span class="delete-button ui-icon ui-icon-trash" alt="Delete this bin" title="Delete this bin" onClick="deleteBin({id});"></span></center></td>' +
                   '</tr>';

    template = template.replace(new RegExp('{id}', 'g'), id)
                       .replace(new RegExp('{name}', 'g'), name)
                       .replace(new RegExp('{color}', 'g'), color)
                       .replace(new RegExp('{count}', 'g'), contig_count)
                       .replace(new RegExp('{completeness}', 'g'), completeness)
                       .replace(new RegExp('{contamination}', 'g'), contamination)
                       .replace(new RegExp('{length}', 'g'), contig_length);


    $('#tbody_bins').append(template);

    if(!from_state){
        $('#completeness_' + id).attr("disabled", true);
        $('#contamination_' + id).attr("disabled", true);
    }

    $('#bin_color_' + id).colpick({
        layout: 'hex',
        submit: 0,
        colorScheme: 'dark',
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
    });
}

function deleteBin(id) {
    if (confirm('Are you sure?')) {

        $('#bin_row_' + id).remove();
        $('#tbody_bins input[type=radio]').last().prop('checked', true);
        bin_count--;

        for (var i = 0; i < SELECTED[id].length; i++) {
            var node_id = label_to_node_map[SELECTED[id][i]].id;
            $("#line" + node_id).css('stroke-width', '1');
            $("#arc" + node_id).css('stroke-width', '1');
            $("#line" + node_id).css('stroke', LINE_COLOR);
            $("#arc" + node_id).css('stroke', LINE_COLOR);
        }

        SELECTED[id] = [];
        delete completeness_dict[id];

        if (bin_count==0)
        {
            newBin();
        }

        redrawBins();
    }
}

function showGenSummaryWindow() {
    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/collections?timestamp=' + new Date().getTime(),
        success: function(data) {
            $('#summaryCollection_list').empty();

            for (source in data) {
                var read_only = data[source]["read_only"];

                if (read_only) {
                    var _name = source + ' (read only)';
                }
                else
                {
                    var _name = source;
                }

                $('#summaryCollection_list').append('<option value="' + source + '">' + _name + '</option>');
            }

            $('#summarizeCollectionWindow').dialog('open');
        }
    });
}


function updateBinsWindow(bin_list) {
    if (typeof bin_list === 'undefined')
    {
        var bin_list = [];
        $('#tbody_bins tr').each(
        function(index, bin) {
            bin_list.push(parseInt($(bin).attr('bin-id')));
        });
    }

    for (var _i = 0; _i < bin_list.length; _i++) {
        var bin_id = bin_list[_i];
        var contigs = 0;
        var length_sum = 0;

        for (var j = 0; j < SELECTED[bin_id].length; j++) {
            if (label_to_node_map[SELECTED[bin_id][j]].IsLeaf())
            {
                contigs++;
                length_sum += parseInt(contig_lengths[SELECTED[bin_id][j]]);
            }
        }

        $('#contig_count_' + bin_id).val(contigs);
        $('#contig_length_' + bin_id).html(readableNumber(length_sum));

        split_names = getContigNames(bin_id);
        bin_name = $('#bin_name_' + bin_id).val();

        updateComplateness(bin_id);
    }

    sortBins();
}

function updateComplateness(bin_id) {
    $.ajax({
        type: "POST",
        url: "/data/completeness",
        cache: false,
        data: {split_names: JSON.stringify(split_names), bin_name: JSON.stringify(bin_name)},
        success: function(data){
            completeness_info_dict = JSON.parse(data);

            stats = completeness_info_dict['stats'];
            refs = completeness_info_dict['refs'];

            completeness_dict[bin_id] = completeness_info_dict;

            showCompleteness(bin_id, true);
            showContaminants(bin_id, true);

            // it is sad that one liner with a list comprehension takes this much code in JS:
            sum_completeness = 0.0;
            sum_contamination = 0.0;
            for(source in stats){
                sum_completeness += stats[source]['percent_complete'];
                sum_contamination += stats[source]['percent_contamination'];
            }
            average_completeness = sum_completeness / Object.keys(stats).length;
            average_contamination = sum_contamination / Object.keys(stats).length;

            $('#completeness_' + bin_id).val(average_completeness.toFixed(1) + '%');
            $('#contamination_' + bin_id).val(average_contamination.toFixed(1) + '%');

            $('#completeness_' + bin_id).attr("disabled", false);
            $('#contamination_' + bin_id).attr("disabled", false);
        },
    });
}

function showCompleteness(bin_id, updateOnly) {
    if (!completeness_dict.hasOwnProperty(bin_id))
        return;

    if ($('#completeness_content_' + bin_id).length)
    {
        $('#completeness_content_' + bin_id).html(buildCompletenessTable(completeness_dict[bin_id]));
        return;
    }

    if (!updateOnly)
    {
        $('<div> \
           <div id="completeness_content_' + bin_id + '">' + buildCompletenessTable(completeness_dict[bin_id]) + '</div> \
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

function showContaminants(bin_id, updateOnly) {
    if (!completeness_dict.hasOwnProperty(bin_id))
        return;

    if ($('#contaminants_content_' + bin_id).length)
    {
        $('#contaminants_content_' + bin_id).html(buildContaminantsTable(completeness_dict[bin_id]));
        return;
    }

    if (!updateOnly)
    {
        $('<div> \
           <div id="contaminants_content_' + bin_id + '">' + buildContaminantsTable(completeness_dict[bin_id]) + '</div> \
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
                                    onclick="redrawBins(' + order_array + ');"> \
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

    // draw bin and layer legend to output svg
    var settings = serializeSettings();

    var bins_to_draw = new Array();
    $('#tbody_bins tr').each(
        function(index, bin) {
            var bin_id = $(bin).attr('bin-id');
            bins_to_draw.push({
                'name': $('#bin_name_' + bin_id).val(),
                'color': $('#bin_color_' + bin_id).attr('color'),
                'contig-length': $('#contig_length_' + bin_id).html(),
                'contig-count': $('#contig_count_' + bin_id).val(),
            });
        }
    );

    var left = 0 - total_radius - 400; // draw on the left top
    var top = 20 - total_radius;

    if (bins_to_draw.length > 0) {
        drawBinLegend(bins_to_draw, top, left);
        top = top + 100 + (bins_to_draw.length + 2.5) * 20
    }

    // important,
    // we used current settings because we want current bin information.
    // now we are going to use "last_settings" which updated by draw button.
    var settings = {};
    settings = last_settings; 

    drawLayerLegend(settings['views'][current_view], settings['layer-order'], top, left);

    svgCrowbar();

    $('#bin_legend').remove();
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

    redrawBins(order_list); 
}

function appendResult() {
    // check if tree exists
    if ($.isEmptyObject(label_to_node_map)) {
        alert('Draw tree first.');
        return;
    }

    var bin_id = getBinId();

    if (bin_id === 'undefined')
        return;

    var bins_to_update = [];
    var _len = search_results.length;
    for (var i=0; i < _len; i++) {
        _contig_name = metadata[search_results[i]][0];
        if (SELECTED[bin_id].indexOf(_contig_name) == -1) {
            SELECTED[bin_id].push(_contig_name);

            if (bins_to_update.indexOf(bin_id) == -1)
                bins_to_update.push(bin_id);
        }

        for (var bid = 1; bid <= bin_counter; bid++) {
            // don't remove nodes from current bin
            if (bid == bin_id)
                continue;

            var pos = SELECTED[bid].indexOf(_contig_name);
            if (pos > -1) {
                SELECTED[bid].splice(pos, 1);

                if (bins_to_update.indexOf(bid) == -1)
                    bins_to_update.push(bid);
            }
        }
    }

    updateBinsWindow(bins_to_update);
    redrawBins();
}

function removeResult() {
    // check if tree exists
    if ($.isEmptyObject(label_to_node_map)) {
        alert('Draw tree first.');
        return;
    }

    var bin_id = getBinId();

    if (bin_id === 'undefined')
        return;

    var bins_to_update = [];
    var _len = search_results.length;
    for (var i=0; i < _len; i++) {
        _contig_name = metadata[search_results[i]][0];
        var _id = label_to_node_map[_contig_name].id;

        var pos = SELECTED[bin_id].indexOf(_contig_name);
        if (pos > -1) {
            SELECTED[bin_id].splice(pos, 1);
            
            if (bins_to_update.indexOf(bin_id) == -1)
                bins_to_update.push(bin_id);
        }
    }

    updateBinsWindow(bins_to_update);
    redrawBins();
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


function storeRefinedBins() {
    data = {};
    colors = {};

    $('#tbody_bins tr').each(
        function(index, bin) {
            var bin_id = $(bin).attr('bin-id');
            var bin_name = $('#bin_name_' + bin_id).val();

            colors[bin_name] = $('#bin_color_' + bin_id).attr('color');
            data[bin_name] = new Array();

            for (var i=0; i < SELECTED[bin_id].length; i++)
            {
                if (label_to_node_map[SELECTED[bin_id][i]].IsLeaf())
                {
                    data[bin_name].push(SELECTED[bin_id][i]);
                }
            }
        }
    );

    $.post("/store_refined_bins", {
        data: JSON.stringify(data, null, 4),
        colors: JSON.stringify(colors, null, 4),
    },
    function(server_response, status){
          alert("Server: " + server_response + "\n\n(status: " + status + ")");
    });
}

function storeCollection() {
    var collection_name = $('#storeCollection_name').val();

    data = {};
    colors = {};

    $('#tbody_bins tr').each(
        function(index, bin) {
            var bin_id = $(bin).attr('bin-id');
            var bin_name = $('#bin_name_' + bin_id).val();

            colors[bin_name] = $('#bin_color_' + bin_id).attr('color');
            data[bin_name] = new Array();

            for (var i=0; i < SELECTED[bin_id].length; i++)
            {
                if (label_to_node_map[SELECTED[bin_id][i]].IsLeaf())
                {
                    data[bin_name].push(SELECTED[bin_id][i]);
                }
            }
        }
    );

    $.post("/store_collection", {
        source: collection_name,
        data: JSON.stringify(data, null, 4),
        colors: JSON.stringify(colors, null, 4),
    },
    function(server_response, status){
          alert("Server: " + server_response + "\n\n(status: " + status + ")");
    });

    $('#storeCollectionWindow').dialog('close');    
}

function generateSummary() {
    var collection = $('#summaryCollection_list').val();

    if (collection === null)
        return;

    alert('This may take a while. Please click OK and wait while anvio generates the summary page..');

    $.ajax({
        type: 'GET',
        cache: false,
        url: '/summarize/' + collection + '?timestamp=' + new Date().getTime(),
        success: function(data) {
            $('#summarizeCollectionWindow').dialog('close');
            if ('error' in data){
                alert(data['error']);
            } else {
                $('#summaryResultWindow').dialog('open');
                $('#summaryDestination').html('<a href="file://' + data['success'] + '">here</a>.');
            }
        }
    });
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

function showCollectionDetails(target_id) {
    var cname = $('#' + target_id + 'Collection_list>option:selected').val();
    
    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/collections?timestamp=' + new Date().getTime(),
        success: function(data) {
            var tbl = "<table>" +
                "<tr><td>Read Only: </td><td>" + data[cname]['read_only'] + "</td></tr>" + 
                "<tr><td colspan='2'><hr></td></tr>" +
                "<tr><td>Number of Splits: </td><td>" + data[cname]['num_splits'] + "</td></tr>" +
                "<tr><td>Number of Clusters: </td><td>" + data[cname]['num_clusters'] + "</td></tr>" +
                "</table>";

            $('#' + target_id + 'Collection_details').html(tbl);
        }
    });
}

function loadCollection() {
    if ($.isEmptyObject(label_to_node_map)) {
        alert('You should draw tree before load collection.');
        return;
    }

    var collection = $('#loadCollection_list').val();
    if (collection === null) {
        alert('Please select a collection.');
        return;
    }

    if (!confirm("You will lost current bins, please be sure you stored current bins. Do you want to continue?"))
        return;

    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/collection/' + collection + '?timestamp=' + new Date().getTime(),
        success: function(data) {
            // empty bin window
            var bins_cleared = false;
            SELECTED = new Array();
            bin_count = 0;
            bin_counter = 0;

            // calculate treshold.
            var threshold = parseFloat($('#loadCollection_threshold').val()) * $('#loadCollection_threshold_base').val();

            // load new bins
            var bin_id=0;
            for (bin in data['data'])
            {
                // collection may be contain unknown splits/contigs, we should clear them.
                var contigs = new Array();
                var sum_contig_length = 0;

                for (index in data['data'][bin])
                {
                    if (typeof contig_lengths[data['data'][bin][index]] !== 'undefined') {
                        contigs.push(data['data'][bin][index]);
                        sum_contig_length += contig_lengths[data['data'][bin][index]];
                    }
                    
                }

                if (sum_contig_length >= threshold)
                {
                    if (!bins_cleared)
                    {
                        $('#tbody_bins').empty();
                        bins_cleared = true;
                    }
                    bin_id++;
                    bin_counter++;
                    SELECTED[bin_id] = contigs;

                    var _color =  (data['colors'][bin]) ? data['colors'][bin] : randomColor();

                    newBin(bin_id, {'name': bin, 'color': _color});
                }
            }

            rebuildIntersections();
            updateBinsWindow();
            redrawBins();
            $('#loadCollectionWindow').dialog('close');
        }
    });
}

function showSaveStateWindow()
{
    $.ajax({
        type: 'GET',
        cache: false,
        url: '/state/all?timestamp=' + new Date().getTime(),
        success: function(state_list) {
            $('#saveState_list').empty();

            for (state_name in state_list) {
                var _select = "";
                if (state_name == current_state_name)
                {
                    _select = ' selected="selected"'; 
                }
                $('#saveState_list').append('<option ' + _select + '>' + state_name + '</option>');
            }

            $('#saveStateWindow').dialog('open');
            $('#saveState_list').trigger('change');
        }
    });
}

function saveState() 
{
    var name = $('#saveState_name').val();

    $.ajax({
        type: 'GET',
        cache: false,
        url: '/state/all?timestamp=' + new Date().getTime(),
        success: function(state_list) {
            for (state_name in state_list) {
                if (state_name == name)
                {
                    if (!confirm('"' + name + '" already exist, do you want to overwrite it?'))
                        return;
                }
            }

        }
    });

    $.ajax({
        type: 'POST',
        cache: false,
        url: '/state/save?timestamp=' + new Date().getTime(),
        data: {
            'name': name,
            'content': JSON.stringify(serializeSettings(true), null, 4)
        },
        success: function(response) {
            response = JSON.parse(response);

            if (response['status_code']==0)
            {
                alert("Failed, Interface running in read only mode.");
            }
            else if (response['status_code']==1)
            {
                // successfull
                $('#saveStateWindow').dialog('close');

                current_state_name = name;
                $('#working_on').html('Working on: ' + current_state_name);
            }
        }
    });
}


function showLoadStateWindow()
{
    $.ajax({
        type: 'GET',
        cache: false,
        url: '/state/all?timestamp=' + new Date().getTime(),
        success: function(state_list) {
            $('#loadState_list').empty();

            for (state_name in state_list) {
                $('#loadState_list').append('<option lastmodified="' + state_list[state_name]['last_modified'] + '">' + state_name + '</option>');
            }

            $('#loadStateWindow').dialog('open');
        }
    });
}

function loadState()
{
    $.ajax({
        type: 'POST',
        cache: false,
        url: '/state/get?timestamp=' + new Date().getTime(),
        data: {'name': $('#loadState_list').val() },
        success: function(response) {
            try{
                var state = JSON.parse(response);
            }catch(e){
                alert('Failed to parse state data, ' + e);
                return;
            }

            if ((state['version'] !== VERSION) || !state.hasOwnProperty('version'))
            {
                alert("Version of the given state file doesn't match with version of the interactive tree, ignoring state file.");
                return;
            }

            if (state.hasOwnProperty('layer-order')) {
                layer_order = [];
                for (var i = 0; i < state['layer-order'].length; i++)
                {
                    // remove non-exists layers.
                    var layer_id = getLayerId(state['layer-order'][i]);

                    if (layer_id != -1)
                    {
                        layer_order.push(layer_id);
                    }
                }

                // add layers that not exist in state and exist in metadata
                for (var i=1; i < parameter_count; i++)
                {
                    if ($.inArray(i, layer_order) === -1)
                    {
                        layer_order.push(i);
                    }
                }

            } else {
                layer_order = Array.apply(null, Array(parameter_count-1)).map(function (_, i) {return i+1;}); // range(1, parameter_count)
            }

            if (state.hasOwnProperty('views')) {
                views = {};
                for (var view_key in state['views'])
                {
                    views[view_key] = {};
                    for (var key in state['views'][view_key])
                    {
                        var layer_id = getLayerId(key);
                        if (layer_id != -1)
                        {
                            views[view_key][layer_id] = state['views'][view_key][key];
                        }
                    }
                }
            }

            if (state.hasOwnProperty('layers')) {
                layers = {};
                for (var key in state['layers'])
                {
                    var layer_id = getLayerId(key);
                    if (layer_id != -1)
                    {
                        layers[layer_id] = state['layers'][key];
                    }
                }
            }

            if (state.hasOwnProperty('categorical_data_colors')) {
                categorical_data_colors = {};
                for (var key in state['categorical_data_colors'])
                {
                    var layer_id = getLayerId(key);
                    if (layer_id != -1)
                    {
                        categorical_data_colors[layer_id] = state['categorical_data_colors'][key];
                    }
                }
            }

            if (state.hasOwnProperty('stack_bar_colors')) {
                stack_bar_colors = {};
                for (var key in state['stack_bar_colors'])
                {
                    var layer_id = getLayerId(key);
                    if (layer_id != -1)
                    {
                        stack_bar_colors[layer_id] = state['stack_bar_colors'][key];
                    }
                }
            }

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
            if (state.hasOwnProperty('tree-radius')) {
                $('#tree-radius-container').show();
                $('#tree-radius').val(state['tree-radius']);
            }
            if (state.hasOwnProperty('layer-margin'))
                $('#layer-margin').val(state['layer-margin']);
            if (state.hasOwnProperty('outer-ring-height'))
                $('#outer-ring-height').val(state['outer-ring-height']);
            if (state.hasOwnProperty('edge-normalization'))
                $('#edge_length_normalization').prop('checked', state['edge-normalization']);
            if (state.hasOwnProperty('custom-layer-margin')) {
                $('#custom_layer_margin').prop('checked', state['custom-layer-margin']).trigger('change');
            }
            if (state.hasOwnProperty('show-grid-for-bins')) {
                $('#show_grid_for_bins').prop('checked', state['show-grid-for-bins']).trigger('change');
            }

            // reload layers
            var current_view = $('#views_container').val();
            $("#tbody_layers").empty();
            buildLayersTable(layer_order, views[current_view]);

            current_state_name = $('#loadState_list').val();
            $('#working_on').html('Working on: ' + current_state_name);
            $('#loadStateWindow').dialog('close');
        }
    });
}


function sortBins(column)
{
    if (typeof column !== 'undefined')
    {
        // clear exist sort icons.
        $('.sort-icon').attr('class', 'sort-icon');

        if (sort_column != column)
        {
            sort_column = column;
            sort_order = 'asc';
        }
        else
        {
            if (sort_order != 'desc')
            {
                sort_order = 'desc';
            }
            else
            {
                sort_order = 'asc';
            }
        }

        var sort_icon_str;
        if (sort_order == 'asc')
            sort_icon_str = 'n'; // north
        else
            sort_icon_str = 's'; // south


        $('#sort-icon-' + column).attr('class', 'sort-icon ui-icon ui-icon-carat-1-' + sort_icon_str);  
    }

    if (sort_column)
    {
        var rows = [];
        $('#tbody_bins tr').each(function() { rows.push($(this).attr('bin-id')) } );
        
        rows.sort(customSortFunction);

        for (var i=0; i < rows.length; i++)
        {
            $('#tbody_bins').append($('#bin_row_' + rows[i]).detach());
        }
    }
}

var bases = {'K': 1000, 'M': 1000000, 'G': 1000000000};
function customSortFunction(a,b)
{
    var retval = 0;

    if (sort_column == 'contigs')
    {
        if (parseInt($('#contig_count_' + a).val()) > parseInt($('#contig_count_' + b).val()))
            retval = 1;
        else
            retval = -1;
    }
    else if (sort_column == 'name')
    {
        if ($('#bin_name_' + a).val() > $('#bin_name_' + b).val())
            retval = 1;
        else
            retval = -1;
    }
    else if (sort_column == 'length')
    {
        var vala = $('#contig_length_' + a).html();
        var valb = $('#contig_length_' + b).html();

        var basea = vala.replace(/[0-9\.]/g, "");
        if (basea == '')
            basea = 1;
        else
            basea = bases[basea];

        var baseb = valb.replace(/[0-9\.]/g, "");
        if (baseb == '')
            baseb = 1;
        else
            baseb = bases[baseb];

        vala = parseFloat('0' + vala) * basea;
        valb = parseFloat('0' + valb) * baseb;
        console.log(vala,valb);
        if (vala > valb)
            retval = 1;
        else
            retval = -1;
    }
    else
    {
        // comp/cont
        var obj_id;
        if (sort_column == 'cont')
        {
            obj_id = '#contamination_';
        }
        else
        {
            obj_id = '#completeness_'
        }

        if (parseFloat('0' + $(obj_id + a).val()) > parseFloat('0' + $(obj_id + b).val()))
            retval = 1;
        else
            retval = -1;
    }

    if (sort_order == 'desc')
    {
        retval = retval * -1;
    }
    return retval;
}
