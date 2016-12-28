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
var VERSION = '0.2.1';

var VIEWER_WIDTH;
var VIEWER_HEIGHT;
var dragging = false;

var zoomBox = {};
var drawing_zoom = false;

var LINE_COLOR='#888888';
var MONOSPACE_FONT_ASPECT_RATIO = 0.6;

var scale = 0;

var id_to_node_map = new Array();
var label_to_node_map = {};
var order_to_node_map = {};
var leaf_count;
var samples_id_to_node_map;

var angle_per_leaf;
var height_per_leaf;
var tree_type;
var margin;
var order_counter;

var total_radius = 0;
var layer_boundaries;

var SELECTED = new Array();
var clusteringData;
var hasTree;

var layerdata;
var contig_lengths;
var parameter_count;

var bin_counter = 0; // for id
var bin_count = 0;

var layer_types;

var categorical_data_colors = {};
var stack_bar_colors = {};

var context_menu_target_id = 0;
var context_menu_layer_id = 0;

var layerdata_title = {};
var layerdata_dict;
var empty_tooltip = "";

var last_settings;

var search_column;
var search_results = [];
var highlighted_splits = [];

var views = {};
var layers = {};
var current_view = '';
var layer_order;

var completeness_dict = {};
var PC_bins_summary_dict = {}

var sort_column;
var sort_order;

var bin_prefix;

var current_state_name = "";

var unique_session_id;
var ping_timer;
var autoload_state;
var autoload_collection;
var mode;
var samples_tree_hover = false;
var bbox;
//---------------------------------------------------------
//  Init
//---------------------------------------------------------

$(document).ready(function() {
    toastr.options = {
        "closeButton": true,
        "debug": false,
        "newestOnTop": true,
        "progressBar": false,
        "positionClass": "toast-top-right",
        "preventDuplicates": false,
        "onclick": null,
        "showDuration": "500",
        "hideDuration": "2000",
        "timeOut": "6000",
        "extendedTimeOut": "1000",
        "showEasing": "swing",
        "hideEasing": "linear",
        "showMethod": "fadeIn",
        "hideMethod": "fadeOut",
    }

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

    initData();

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

    $('#grid_color').colpick({
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

    document.body.addEventListener('click', function() {
        $('#control_contextmenu').hide();
    }, false);

    $("li[role='presentation']").click(function (e) {
        if ($(this).hasClass('disabled')) {
            e.preventDefault();
            return false;
        }  
    });

}); // document ready

function initData () {
    var timestamp = new Date().getTime(); 

    $.when(
    $.ajax({
            type: 'GET',
            cache: false,
            url: '/data/init?timestamp=' + timestamp,
        })
    )
    .then(
    function (response)
        {
        var titleResponse = [ response.title ];
        var clusteringsResponse = [ response.clusterings ];
        var viewsResponse = [ response.views ];
        var contigLengthsResponse = [ response.contigLengths ];
        var defaultViewResponse = [ response.defaultView ];
        var modeResponse = [ response.mode ];
        var readOnlyResponse = [ response.readOnly ];
        var prefixResponse = [ response.binPrefix ];
        var sessionIdResponse = [ response.sessionId ];
        var samplesOrderResponse = [ response.samplesOrder ];
        var sampleInformationResponse = [ response.sampleInformation ];
        var sampleInformationDefaultLayerOrderResponse = [ response.sampleInformationDefaultLayerOrder ];
        var stateAutoloadResponse = [ response.stateAutoload ];
        var collectionAutoloadResponse = [ response.collectionAutoload ];
        var inspectionAvailable = response.inspectionAvailable;
        var sequencesAvailable = response.sequencesAvailable;
            unique_session_id = sessionIdResponse[0];
            mode = modeResponse[0];

        if(!inspectionAvailable){
            toastr.info("Inspection of data items is not going to be available for this project.");
            $('.menuItemInspect').addClass('menu-disabled');
        }

        if(!sequencesAvailable && mode != "collection" && mode != "pan"){
            toastr.info("No sequence data is available. Some menu items will be disabled.");
            $('.menuItemSequence').addClass('menu-disabled');
        }

        if (! response.noPing) {
            ping_timer = setInterval(checkBackgroundProcess, 5000);
        }

            if (!$.browser.chrome)
            {
                toastr.warning("We tested anvi'o only on Google Chrome, and it seems you are using a different browser. For the best performance, and to avoid unexpected issues, please consider using anvi'o with the lastest version of Chrome.", "", { 'timeOut': '0', 'extendedTimeOut': '0' });
            }

            // hide all mode dependent divs:
            $('.full-mode').hide();
            $('.pan-mode').hide();
            $('.collection-mode').hide();
            $('.manual-mode').hide();
            $('.server-mode').hide();
            $('.refine-mode').hide();

            console.log("The running mode for the interface: " + mode);

            // mode switch:
            if (mode == 'refine')
            {
                $('.refine-mode').show();
                $('#fixed-navbar-div').css('background-image', 'url(images/refine-bg.png)');
            } else if (mode == 'server') {
                $('.server-mode').show();
                $('#fixed-navbar-div').css('background-image', 'url(images/server-bg.png)');
            } else if (mode == 'full') {
                $('.full-mode').show();
                $('#fixed-navbar-div').css('background-image', 'url(images/full-bg.png)');
            } else if (mode == 'pan') {
                $('.pan-mode').show();
                $('#fixed-navbar-div').css('background-image', 'url(images/pan-bg.png)');

                $('#completion_title').attr('title', 'PCs').html('PCs');
                $('#redundancy_title').attr('title', 'Gene Calls').html('Gene Calls');
                $('#splits_title').hide();
                $('#len_title').hide();

            } else if (mode == 'collection') {
                $('.collection-mode').show();
                $('#fixed-navbar-div').css('background-image', 'url(images/collection-bg.png)');
            } else if (mode == 'manual') {
                $('.manual-mode').show();
                $('#fixed-navbar-div').css('background-image', 'url(images/manual-bg.png)');
            }

            if (readOnlyResponse[0] == true)
            {
                toastr.info("It seems that this is a read-only instance, therefore the database-writing functions will be inaccessible.", "", { 'timeOut': '0', 'extendedTimeOut': '0' });
                $('[disabled-in-read-only=true]').addClass('disabled').prop('disabled', true);
            }

            bin_prefix = prefixResponse[0];

            document.title = titleResponse[0];
            contig_lengths = eval(contigLengthsResponse[0]);

            // if --state parameter given, autoload given state.
            autoload_state = stateAutoloadResponse[0];

            // if --collection parameter given, autoload given collection.
            autoload_collection = collectionAutoloadResponse[0];

            /* 
            //  Clusterings
            */
            var default_tree = clusteringsResponse[0][0];
            var available_trees = clusteringsResponse[0][1];
            var available_trees_combo = getComboBoxContent(default_tree, available_trees);

            $('#trees_container').append(available_trees_combo);
            $('#trees_container').change(function() {
                onTreeClusteringChange();
            });

            /* 
            //  Views
            */
            var default_view = viewsResponse[0][0];
            var available_views = viewsResponse[0][1];
            var available_views_combo = getComboBoxContent(default_view, available_views);

            $('#views_container').append(available_views_combo);
            $('#views_container').change(function() {
                onViewChange();
            });

            // make layers and samples table sortable
            var _notFirstSelector = ''
            if (mode != 'manual' && mode != 'pan' && mode != 'server') {
                _notFirstSelector = ':not(:first)';
            }
            $("#tbody_layers").sortable({helper: fixHelperModified, handle: '.drag-icon', items: "> tr" + _notFirstSelector}).disableSelection(); 
            $("#tbody_samples").sortable({helper: fixHelperModified, handle: '.drag-icon', items: "> tr"}).disableSelection(); 

            samples_order_dict = samplesOrderResponse[0];
            samples_information_dict = sampleInformationResponse[0];
            samples_information_default_layer_order = sampleInformationDefaultLayerOrderResponse[0];

            available_orders = Object.keys(samples_order_dict).sort();
            $('#samples_order').append(new Option('custom'));
            available_orders.forEach(function(order)
            {
                var order_name = order;
                if (samples_order_dict[order]['newick'] != null && samples_order_dict[order]['newick'] != '')
                    order_name += " (tree)";

                $('#samples_order').append(new Option(order_name, order));
            });
            buildSamplesTable(samples_information_default_layer_order);

            // load default data
            if (autoload_state !== null)
            {
                $.when({}).then(onViewChange)
                      .then(loadState)
                      .then(onViewChange)
                      .then(onTreeClusteringChange)
                      .then(function() {
                        drawTree();
                      });                
            }
            else
            {
                $.when({}).then(onViewChange)
                          .then(onTreeClusteringChange);          
            }

            /*
            //  Add bins
            */
            newBin();
        } // response callback
    ).fail(function() {
        toastr.error("One or more ajax request has failed, See console logs for details.", "", { 'timeOut': '0', 'extendedTimeOut': '0' });
        console.log(arguments);
    }); // promise
}

function onViewChange() {
    var defer = $.Deferred();
    console.log('View data ' + $('#views_container').val() + ' requested.');

    $('#views_container').prop('disabled', false);
    $('#btn_draw_tree').prop('disabled', true);

    waitingDialog.show('Requesting view data from the server ...', 
        {
            dialogSize: 'sm', 
            onHide: function() {
                defer.resolve(); 
            },
            onShow: function() {
                $.ajax({
                    type: 'GET',
                    cache: false,
                    url: '/data/view/' + $('#views_container').val() + '?timestamp=' + new Date().getTime(),
                    success: function(data) {
                        layerdata = eval(data);
                        parameter_count = layerdata[0].length;

                        // since we are painting parent layers odd-even, 
                        // we should remove single parents (single means no parent)
                        removeSingleParents(); // in utils.js

                        layer_order = Array.apply(null, Array(parameter_count-1)).map(function (_, i) {return i+1;}); // range(1, parameter_count)
                        layer_types = {};

                        // add layerdata columns to search window
                        $('#searchLayerList').empty();
                        for (var i=0; i < layerdata[0].length; i++)
                        {
                            $('#searchLayerList').append(new Option(layerdata[0][i],i));
                        }

                        $('#views_container').attr('disabled', false);
                        $('#btn_draw_tree').attr('disabled', false);

                        if (current_view != '') {
                            // backup current layer order and layers table to global views object
                            syncViews();
                        }
                        current_view = $('#views_container').val();

                        $("#tbody_layers").empty();

                        buildLayersTable(layer_order, views[current_view]);  

                        waitingDialog.hide();
                    }
                });
            },
        });

    return defer.promise();
}

function onTreeClusteringChange() {
    var defer = $.Deferred();
    console.log('Tree clustering data ' + $('#trees_container').val() + ' requested.');
    $('#trees_container').prop('disabled', true);
    $('#btn_draw_tree').prop('disabled', true);

    waitingDialog.show('Requesting the tree data ...', 
        {
            dialogSize: 'sm', 
            onHide: function() { 
                defer.resolve(); 
            },
            onShow: function() {    
                $.ajax({
                    type: 'GET',
                    cache: false,
                    url: '/tree/' + $('#trees_container').val() + '?timestamp=' + new Date().getTime(),
                    success: function(data) {
                        clusteringData = data;
                        $('#trees_container').attr('disabled', false);
                        $('#btn_draw_tree').attr('disabled', false); 
                        waitingDialog.hide();
                    }
                });
            },
        });

    return defer.promise();
}

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

            if (layers[layer_id]["type"] === 'text')
                layers[layer_id]["height"] = '0';
        }
    );    
}


function getComboBoxContent(default_item, available_items){
    var combo = '';
    var combo_item = '<option value="{val}"{sel}>{text}</option>';

    $.each(available_items, function(index, value) {
        if (index.indexOf(':') == -1) {
            text_val = getPrettyName(index);
        } else {
            text_val = getClusteringPrettyName(index);
        }

        if(index == default_item)
        {
            combo += combo_item
                        .replace('{val}', index)
                        .replace('{sel}', ' selected')
                        .replace('{text}', text_val);
        }
        else
        {
            combo += combo_item
                        .replace('{val}', index)
                        .replace('{sel}', '')
                        .replace('{text}', text_val);
        }
    });

    return combo;
}

function buildLayersTable(order, settings)
{
    for (var i = 0; i < order.length; i++) 
    {
        // common layer variables
        var layer_id = order[i];
        var layer_name = layerdata[0][layer_id];

        var short_name = (layer_name.indexOf('!') > -1) ? layer_name.split('!')[0] : layer_name;
        short_name = (short_name.length > 10) ? short_name.slice(0,10) + "..." : short_name;

        var hasViewSettings = false;
        if (typeof settings !== 'undefined' && typeof settings[layer_id] !== 'undefined') {
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
                '<td><input type="checkbox" class="layer_selectors"></input></td>' +
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
                var height = '300';
                var margin = '15';

                // pick random color for stack bar items
                if (!(layer_id in stack_bar_colors))
                {
                    stack_bar_colors[layer_id] = new Array();
                    var bars = (layer_name.indexOf('!') > -1) ? layer_name.split('!')[1].split(';') : layer_name.split(';');
                    for (var j=0; j < bars.length; j++)
                    {
                        stack_bar_colors[layer_id].push(randomColor({luminosity: 'dark'}));
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
                '        <option value="sqrt"{option-sqrt}>sqrt</option>' +
                '        <option value="log"{option-log}>log</option>' +
                '    </select>' +
                '</td>' +
                '<td><input class="input-height" type="text" size="3" id="height{id}" value="{height}"></input></td>' +
                '<td class="column-margin"><input class="input-margin" type="text" size="3" id="margin{id}" value="{margin}"></input></td>' +
                '<td>n/a</td>' +
                '<td>n/a</td>' +
                '<td><input type="checkbox" class="layer_selectors"></input></td>' +
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
        else if (layerdata[1][layer_id] === null || !isNumber(layerdata[1][layer_id]))
        { 
            layer_types[layer_id] = 2;

            if (hasLayerSettings)
            {
                var height = layer_settings['height'];
                var margin = layer_settings['margin'];
                var type = layer_settings['type'];
                var color = layer_settings['color'];
                var color_start = layer_settings['color-start'];
            }
            else
            {
                var color = "#000000";
                var height = '90';
                var margin = '15';
                var color_start = "#DDDDDD";

                if (mode == 'collection') {
                    var type = getNamedLayerDefaults(layer_name, 'type', 'color');
                    $('.max-font-size-input').show();
                } else {
                    var type = 'color';
                }

                // set default categorical layer type to 'text' 
                // if there are more than 11 unique values and leaf count is less than 300
                // 301 because layerdata has one extra row for the titles
                if (layerdata.length <= 301)
                {
                    var _unique_items = [];
                    for (var _pos = 1; _pos < layerdata.length; _pos++)
                    {
                        if (_unique_items.indexOf(layerdata[_pos][layer_id]) === -1)
                            _unique_items.push(layerdata[_pos][layer_id]);

                        if (_unique_items.length > 11) {
                            height = '0';
                            type = 'text';
                            // we have at least one text layer, we can show max font size input
                            $('.max-font-size-input').show();
                            break;
                        }
                    }
                }

                if (!(layer_id in categorical_data_colors))
                {
                    categorical_data_colors[layer_id] = {};
                }
            }
            
            var template = '<tr>' +
                '<td><img class="drag-icon" src="images/drag.gif" /></td>' +
                '<td title="{name}" class="titles" id="title{id}">{short-name}</td>' +
                '<td><div id="picker_start{id}" class="colorpicker picker_start" color="{color-start}" style="background-color: {color-start}; {color-start-hide}"></div><div id="picker{id}" class="colorpicker picker_end" color="{color}" style="background-color: {color}; {color-hide}"></div></td>' +
                '<td style="width: 50px;">' +
                '    <select id="type{id}" style="width: 50px;" class="type" onChange="togglePickerStart(this, true);">' +
                '        <option value="color"{option-type-color}>Color</option>' +
                '        <option value="text"{option-type-text}>Text</option>' +
                '    </select>' +
                '</td>' +
                '<td>n/a</td>' +
                '<td><input class="input-height" type="text" size="3" id="height{id}" value="{height}" style="{height-hide}"></input></td>' +
                '<td class="column-margin"><input class="input-margin" type="text" size="3" id="margin{id}" value="{margin}"></input></td>' +
                '<td>n/a</td>' +
                '<td>n/a</td>' +
                '<td><input type="checkbox" class="layer_selectors"></input></td>' +
                '</tr>';

            template = template.replace(new RegExp('{id}', 'g'), layer_id)
                               .replace(new RegExp('{name}', 'g'), layer_name)
                               .replace(new RegExp('{option-type-' + type + '}', 'g'), ' selected')
                               .replace(new RegExp('{option-type-([a-z]*)}', 'g'), '')
                               .replace(new RegExp('{short-name}', 'g'), short_name)
                               .replace(new RegExp('{color}', 'g'), color)
                               .replace(new RegExp('{color-start}', 'g'), color_start)
                               .replace(new RegExp('{color-hide}', 'g'), (type!='text') ? '; visibility: hidden;' : '')
                               .replace(new RegExp('{color-start-hide}', 'g'), (type!='text') ? '; visibility: hidden;' : '')
                               .replace(new RegExp('{height-hide}', 'g'), (type=='text') ? '; visibility: hidden;' : '')
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
                var min    = getNamedLayerDefaults(layer_name, 'min', 0);
                var max    = getNamedLayerDefaults(layer_name, 'max', 0);
                var min_disabled = getNamedLayerDefaults(layer_name, 'min_disabled', true);
                var max_disabled = getNamedLayerDefaults(layer_name, 'max_disabled', true);
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
                if (mode == 'collection') {
                    var type = getNamedLayerDefaults(layer_name, 'type', 'intensity');
                    var color_start = "#EEEEEE";
                } else {
                    var type = 'bar'
                    var color_start = "#FFFFFF";
                }
            }

            var template = '<tr>' +
                '<td><img class="drag-icon" src="images/drag.gif" /></td>' +
                '<td title="{name}" class="titles" id="title{id}">{short-name}</td>' +
                '<td><div id="picker_start{id}" class="colorpicker picker_start" color="{color-start}" style="background-color: {color-start}; {color-start-hide}"></div><div id="picker{id}" class="colorpicker" color="{color}" style="background-color: {color}"></div></td>' +
                '<td style="width: 50px;">' +
                '    <select id="type{id}" style="width: 50px;" class="type" onChange="togglePickerStart(this);">' +
                '        <option value="bar"{option-type-bar}>Bar</option>' +
                '        <option value="intensity"{option-type-intensity}>Intensity</option>' +
                '    </select>' +
                '</td>' +
                '<td>' +
                '    <select id="normalization{id}" onChange="clearMinMax(this);" class="normalization">' +
                '        <option value="none"{option-none}>none</option>' +
                '        <option value="sqrt"{option-sqrt}>sqrt</option>' +
                '        <option value="log"{option-log}>log</option>' +
                '    </select>' +
                '</td>' +
                '<td><input class="input-height" type="text" size="3" id="height{id}" value="{height}"></input></td>' +
                '<td class="column-margin"><input class="input-margin" type="text" size="3" id="margin{id}" value="{margin}"></input></td>' +
                '<td><input class="input-min" type="text" size="4" id="min{id}" value="{min}"{min-disabled}></input></td>' +
                '<td><input class="input-max" type="text" size="4" id="max{id}" value="{max}"{min-disabled}></input></td>' +
                '<td><input type="checkbox" class="layer_selectors"></input></td>' +
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
                               .replace(new RegExp('{margin}', 'g'), margin);


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

        $('#table_layers .drag-icon').on('mousedown', function() {
            $('#samples_order').val('custom').trigger('change');
        });
    }
}

function getLayerName(layer_id)
{
    return layerdata[0][layer_id];
}

function getLayerId(layer_name) 
{
    for (var i=0; i < parameter_count; i++)
    {
        if (layer_name == layerdata[0][i])
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
    state['outer-ring-margin'] = $('#outer-ring-margin').val();
    state['edge-normalization'] = $('#edge_length_normalization').is(':checked');
    state['custom-layer-margin'] = $('#custom_layer_margin').is(':checked');
    state['show-grid-for-bins'] = $('#show_grid_for_bins').is(':checked');
    state['grid-color'] = $('#grid_color').attr('color');
    state['grid-width'] = $('#grid_width').val();
    state['samples-order'] = $('#samples_order').val();
    state['max-font-size'] = $('#max_font_size').val();
    state['optimize-speed'] = $('#optimize_speed').is(':checked');
    state['show-bin-labels'] = $('#show_bin_labels').is(':checked');
    state['bin-labels-font-size'] = $('#bin_labels_font_size').val();
    state['autorotate-bin-labels'] = $('#autorotate_bin_labels').is(':checked');
    state['bin-labels-angle'] = $('#bin_labels_angle').val();

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

    state['samples-categorical-colors'] = samples_categorical_colors;
    state['samples-stack-bar-colors'] = samples_stack_bar_colors;
    state['samples-order'] = $('#samples_order').val();
    state['samples-edge-length-normalization'] = $('#samples_edge_length_normalization').is(':checked');
    state['samples-ignore-branch-length'] = $('#samples_ignore_branch_length').is(':checked');
    state['samples-tree-height'] = $('#samples_tree_height').val();

    state['samples-layer-order'] = [];
    state['samples-layers'] = {};
    $('#tbody_samples tr').each(
        function(index, tr) {
            var samples_layer_name = $(tr).attr('samples-layer-name');
            state['samples-layer-order'].push(samples_layer_name);
            state['samples-layers'][samples_layer_name] = {
                'data-type'     : $(tr).attr('data-type'),
                'height'        : parseFloat($(tr).find('.input-height').val()),
                'margin'        : parseFloat($(tr).find('.input-margin').val()),
                'normalization' : $(tr).find('.normalization').val(),
                'color'         : $(tr).find('.colorpicker:last').attr('color'),
                'min'           : {'value': parseFloat($(tr).find('.input-min').val()), 'disabled': $(tr).find('.input-min').is(':disabled') },
                'max'           : {'value': parseFloat($(tr).find('.input-max').val()), 'disabled': $(tr).find('.input-max').is(':disabled') },
                'type'          : $(tr).find('.type').val(),
                'color-start'   : $(tr).find('.colorpicker:first').attr('color'),
            };
        }
    );

    return state;
}

function drawTree() {
    // get current client size
    VIEWER_WIDTH = document.getElementById('svg').clientWidth;
    VIEWER_HEIGHT = document.getElementById('svg').clientHeight;

    var settings = serializeSettings();
    tree_type = settings['tree-type'];

    $('#draw_delta_time').html('');
    $('#btn_draw_tree').prop('disabled', true);
    $('#bin_settings_tab').removeClass("disabled"); // enable bins tab
    $('#sample_settings_tab').removeClass("disabled"); // enable bins tab
    $('#mouse_tooltips_tab').removeClass("disabled"); // enable bins tab
    $('#search_panel_tab').removeClass("disabled"); // enable bins tab


    // clear existing diagram, if any
    document.getElementById('svg').innerHTML = "";

    waitingDialog.show('Drawing ...', 
        {
            dialogSize: 'sm',
            onShow: function() {
                draw_tree(settings); // call treelib.js where the magic happens

                // last_settings used in export svg for layer information,
                // we didn't use "settings" sent to draw_tree because draw_tree updates layer's min&max
                // running serializeSettings() twice costs extra time but we can ignore it to keep code simple.
                last_settings = serializeSettings();

                redrawBins();

                waitingDialog.hide();
                $('#btn_draw_tree').prop('disabled', false);
                $('#btn_redraw_samples').prop('disabled', false);

                if (settings['tree-radius'] == 0)
                {
                    $('#tree-radius-container').show();
                    $('#tree-radius').val(Math.max(VIEWER_HEIGHT, VIEWER_WIDTH));
                }

                if (autoload_collection !== null)
                {
                    loadCollection(autoload_collection);
                    autoload_collection = null
                }
            },
        });
}


function getContigNames(bin_id) {
    var names = new Array();

    for (var j = 0; j < SELECTED[bin_id].length; j++) {
        if (label_to_node_map[SELECTED[bin_id][j]].IsLeaf()) {
            names.push(SELECTED[bin_id][j]);
        }
    }

    return names;
}


function showContigNames(bin_id, updateOnly) {
    if (typeof updateOnly === 'undefined')
        updateOnly = false;

    var title = 'Splits in "' + $('#bin_name_' + bin_id).val() + '"';

    if (updateOnly && !checkObjectExists('#modal' + title.hashCode()))
        return;

    var msg = '<table class="table table-striped">';
    var names = getContigNames(bin_id);

    for (var i in names)
        msg += "<tr><td><a href='#' class='no-link' onclick='highlightSplit(\"" + names[i] + "\");'>" + names[i] + "</a></td></tr>";

    msg = msg + '</table>';

    showDraggableDialog(title, msg, updateOnly);
}

function newBin(id, binState) {

    bin_count++;

    if (typeof id === 'undefined')
    {
        bin_counter++;
        var from_state = false;
        var id = bin_counter;
        var name = bin_prefix + id;
        var color = randomColor({luminosity: 'dark'});
        var contig_count = 0;
        var contig_length = "N/A";
        var completeness = '---';
        var redundancy = '---';

        SELECTED[bin_counter] = [];
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

    var template = '<tr bin-id="{id}" id="bin_row_{id}">' +
                   '    <td><input type="radio" name="active_bin" value="{id}" checked></td>' +
                   '    <td><div id="bin_color_{id}" class="colorpicker" color="{color}" style="background-color: {color}"></td>' +
                   '    <td data-value="{name}"><input type="text" onChange="redrawBins();" size="21" id="bin_name_{id}" value="{name}"></td>';

    if (mode != 'pan')
    {
        template +='    <td data-value="{count}"><input id="contig_count_{id}" type="button" value="{count}" title="Click for contig names" onClick="showContigNames({id});"></td> ' +
                   '    <td data-value="{length}"><span id="contig_length_{id}">{length}</span></td>';
    }

    template +=    '    <td data-value="{completeness}"><input id="completeness_{id}" type="button" value="{completeness}" title="Click for completeness table" onClick="showCompleteness({id});"></td> ' +
                   '    <td data-value="{redundancy}"><input id="redundancy_{id}" type="button" value="{redundancy}" title="Click for redundant hits" onClick="showRedundants({id});"></td> ' +
                   '    <td><center><span class="glyphicon glyphicon-trash" aria-hidden="true" alt="Delete this bin" title="Delete this bin" onClick="deleteBin({id});"></span></center></td>' +
                   '</tr>';

    template = template.replace(new RegExp('{id}', 'g'), id)
                       .replace(new RegExp('{name}', 'g'), name)
                       .replace(new RegExp('{color}', 'g'), color)
                       .replace(new RegExp('{count}', 'g'), contig_count)
                       .replace(new RegExp('{completeness}', 'g'), completeness)
                       .replace(new RegExp('{redundancy}', 'g'), redundancy)
                       .replace(new RegExp('{length}', 'g'), contig_length);

    $('#tbody_bins').append(template);

    if(!from_state){
        $('#completeness_' + id).attr("disabled", true);
        $('#redundancy_' + id).attr("disabled", true);
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

            showCollectionDetails('');
            $('#modGenerateSummary').modal('show');
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

        if (mode === 'pan'){
            updateProteinClustersBin(bin_id);
        } else {
            updateComplateness(bin_id);

            var contigs = 0;
            var length_sum = 0;

            for (var j = 0; j < SELECTED[bin_id].length; j++) {
                if (label_to_node_map[SELECTED[bin_id][j]].IsLeaf())
                {
                    contigs++;
                    length_sum += parseInt(contig_lengths[SELECTED[bin_id][j]]);
                }
            }

            $('#contig_count_' + bin_id).val(contigs).parent().attr('data-value', contigs);

            // it is likely in manual or server modes lenghts are not going to be available.
            if (isNaN(length_sum))
                $('#contig_length_' + bin_id).html('N/A').parent().attr('data-value', 0);
            else
                $('#contig_length_' + bin_id).html(readableNumber(length_sum)).parent().attr('data-value', length_sum);

        }

        showContigNames(bin_id, true);
    }

    $('#bin_settings_tab:not(.active) a').css('color', "#ff0000");
}


function updateProteinClustersBin(bin_id) {
    if (mode !== 'pan'){ 
        return;
    }

    $.ajax({
        type: "POST",
        url: "/data/proteinclusterssummary",
        cache: false,
        data: {split_names: JSON.stringify(getContigNames(bin_id)), bin_name: JSON.stringify($('#bin_name_' + bin_id).val())},
        success: function(data){
            data = JSON.parse(data);
            PC_bins_summary_dict[bin_id] = data;
            $('#redundancy_' + bin_id).val(data['num_gene_calls']).parent().attr('data-value', data['num_gene_calls']);
            $('#completeness_' + bin_id).val(data['num_PCs']).parent().attr('data-value', data['num_PCs']);

            $('#completeness_' + bin_id).attr("disabled", false);
            $('#redundancy_' + bin_id).attr("disabled", false);
        },
    });
}


function updateComplateness(bin_id) {
    if (mode === 'manual' || mode === 'pan' || mode === 'server'){ 
        // there is nothing to do here
        return;
    }

    $.ajax({
        type: "POST",
        url: "/data/completeness",
        cache: false,
        data: {split_names: JSON.stringify(getContigNames(bin_id)), bin_name: JSON.stringify($('#bin_name_' + bin_id).val())},
        success: function(data){
            completeness_info_dict = JSON.parse(data);

            stats = completeness_info_dict['stats'];
            refs = completeness_info_dict['refs'];
            averages = completeness_info_dict['averages'];

            completeness_dict[bin_id] = completeness_info_dict;

            average_completeness = averages['percent_completion']
            average_redundancy = averages['percent_redundancy']

            $('#completeness_' + bin_id).val(average_completeness.toFixed(1) + '%').parent().attr('data-value', average_completeness);
            $('#redundancy_' + bin_id).val(average_redundancy.toFixed(1) + '%').parent().attr('data-value', average_redundancy);

            $('#completeness_' + bin_id).attr("disabled", false);
            $('#redundancy_' + bin_id).attr("disabled", false);

            showCompleteness(bin_id, true);
            showRedundants(bin_id, true);
        },
    });
}

function showCompleteness(bin_id, updateOnly) {
    if (typeof updateOnly === 'undefined')
        updateOnly = false;

    if (!completeness_dict.hasOwnProperty(bin_id))
        return;

    var refs = completeness_dict[bin_id]['refs'];
    var stats = completeness_dict[bin_id]['stats'];
    var averages = completeness_dict[bin_id]['averages'];

    var title = 'Completeness of "' + $('#bin_name_' + bin_id).val() + '"';

    if (updateOnly && !checkObjectExists('#modal' + title.hashCode()))
        return;

    var msg = '<table class="table table-striped sortable">' +
        '<thead><tr><th data-sortcolumn="0" data-sortkey="0-0">Source</th><th data-sortcolumn="1" data-sortkey="1-0">SCG domain</th><th data-sortcolumn="2" data-sortkey="2-0">Percent complenetess</th></tr></thead><tbody>';

    for (var source in stats){
        if(stats[source]['domain'] != averages['domain'])
            // if the source is not matching the best domain recovered
            // don't show it in the interface
            continue;

        msg += "<tr><td data-value='" + source  + "'><a href='" + refs[source] + "' class='no-link' target='_blank'>" + source + "</a></td><td data-value='" + stats[source]['domain'] + "'>" + stats[source]['domain'] + "</td><td data-value='" + stats[source]['percent_complete'] + "'>" + stats[source]['percent_complete'].toFixed(2) + "%</td></tr>";
    }

    msg = msg + '</tbody></table>';

    showDraggableDialog(title, msg, updateOnly);
}

function showRedundants(bin_id, updateOnly) {
    if (typeof updateOnly === 'undefined')
        updateOnly = false;
    
    if (!completeness_dict.hasOwnProperty(bin_id))
        return;

    var stats = completeness_dict[bin_id]['stats'];
    var averages = completeness_dict[bin_id]['averages'];

    var output_title = 'Redundants of "' + $('#bin_name_' + bin_id).val() + '"';

    if (updateOnly && !checkObjectExists('#modal' + output_title.hashCode()))
        return;

    var output = '<div class="col-md-12">'
    var oddeven=0;

    for(var source in stats) {
        if(stats[source]['domain'] !== averages['domain'])
            continue;

        oddeven++;
        var tabletext = '<div class="table-responsive col-md-6"><table style="margin-bottom: 10px;"><tr><td>';
        tabletext += '<h5>' + source + ' (' + Object.keys(stats[source]['redundants']).length + ')</h5></td></tr>';

        var redundants_html = '';

        for (var redundant in stats[source]['redundants']) {
            var title = '';
            var split_array = '';
            for (var i = 0; i < stats[source]['redundants'][redundant].length; i++)
            {
                var contig = stats[source]['redundants'][redundant][i];

                for (var j = 0; j < contig.length; j++) {
                    // splits
                    title += contig[j] + '\n';
                    split_array += '\'' + contig[j] + '\', ';
                }
            }

            redundants_html += '<span style="cursor:pointer;" \
                                    data-toggle="tooltip" data-placement="top" title="' + title + '" \
                                    onclick="highlighted_splits = [' + split_array + ']; redrawBins();"> \
                                    ' + redundant + ' (' + stats[source]['redundants'][redundant].length + ') \
                                  </span><br />';
        }

        tabletext += '<tr><td valign="top">' + redundants_html + '</tr></td></table></div>';
        output += tabletext;

        if (oddeven%2==0)
        {
            output += '</div><div class="col-md-12">'
        }
    }

    output += '</div>';

    showDraggableDialog(output_title, output, updateOnly);
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
    drawLayerLegend(settings['layers'], settings['views'][current_view], settings['layer-order'], top, left);
    var detached_clones = $('#samples_tree path.clone').detach();

    svgCrowbar();

    $('#samples_tree').prepend(detached_clones);
    $('#bin_legend').remove();
    $('#layer_legend').remove();
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

            $('#modStoreCollection').modal('show');
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
        server_response = JSON.parse(server_response);
        if (server_response.status == -1){
            toastr.error(server_response.message, "You made the server upset :(");
        } else {
            toastr.info(server_response.message, "The server is on board");
        }
    });
}

function storeCollection() {
    var collection_name = $('#storeCollection_name').val();

    if (collection_name.length==0) {
        $('#storeCollection_name').focus();
        return;
    }
         
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
        toastr.info(server_response, "Server");
    });

    $('#modStoreCollection').modal('hide');    
}

function generateSummary() {
    var collection = $('#summaryCollection_list').val();

    if (collection === null)
        return;

    waitingDialog.show('Generating summary...', {dialogSize: 'sm'});

    $.ajax({
        type: 'GET',
        cache: false,
        url: '/summarize/' + collection + '?timestamp=' + new Date().getTime(),
        success: function(data) {
            if ('error' in data){
                $('#modGenerateSummary').modal('hide');
                waitingDialog.hide();
                toastr.error(data['error'], "", { 'timeOut': '0', 'extendedTimeOut': '0' });
            } else {
                $('#modGenerateSummary').modal('hide');
                waitingDialog.hide();

                // generate a full url using the window origin and collection path:
                var summary_url = window.location.origin + '/' + data['path'];

                $('#summary_link').html("Summary link: <a href='" + summary_url + "' target='_blank'>" + summary_url + "</a>");
                $('#modSummaryResult').modal('show');
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

            $('#loadCollection_list, #btn-load-collection').prop('disabled', false);
            showCollectionDetails('');
            $('#modLoadCollection').modal('show');
        }
    });
}

function showCollectionDetails(list) {

    var cname = $(list).val();

    if (cname=='' || typeof cname === 'undefined')
    {
        // clear details
        var tbl = '<div class="col-md-12">Collection Details</div><hr>' +
            '<div class="col-md-8">Number of Splits:</div><div class="col-md-4"><b>n/a</b></div>' +
            '<div class="col-md-8">Number of Bins:</div><div class="col-md-4"><b>n/a</b></div>';

        $('.collection-details').html(tbl);

        return;
    }

    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/collections?timestamp=' + new Date().getTime(),
        success: function(data) {
            var tbl = '<div class="col-md-12">Collection Details</div><hr>' +
                '<div class="col-md-8">Number of Splits:</div><div class="col-md-4"><b>' + data[cname]['num_splits'] + '</b></div>' +
                '<div class="col-md-8">Number of Bins:</div><div class="col-md-4"><b>' + data[cname]['num_bins'] + '</b></div>';

            $('.collection-details').html(tbl);
        }
    });
}

function loadCollection(default_collection) {
    if ($.isEmptyObject(label_to_node_map)) {
        toastr.warning('You should draw tree before load collection.');
        return;
    }

    var collection;
    if (default_collection) {
        collection = default_collection;
    } else {
        collection = $('#loadCollection_list').val();
    }

    if (collection === null) {
        toastr.warning('Please select a collection.');
        return;
    }

    $('#loadCollection_list, #btn-load-collection').prop('disabled', true);
    
    var bin_list = [];
    var total_selection = 0;
    
    $('#tbody_bins tr').each(
    function(index, bin) {
        bin_list.push(parseInt($(bin).attr('bin-id')));
    });
    for (var _i = 0; _i < bin_list.length; _i++) {
        var bin_id = bin_list[_i];
        for (var j = 0; j < SELECTED[bin_id].length; j++) {
            if (label_to_node_map[SELECTED[bin_id][j]].IsLeaf())
                total_selection++;
        }
    }

    if (total_selection > 0 && !confirm("You will lost current bins, please be sure you stored current bins. Do you want to continue?")) {
        showLoadCollectionWindow();
        return;
    }

    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/collection/' + collection + '?timestamp=' + new Date().getTime(),
        success: function(data) {
            $('#modLoadCollection').modal('hide');

            // clear bins tab
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
                    if (mode === 'manual' || mode === 'pan' || mode === 'server'){
                        contigs.push(data['data'][bin][index]);
                    } else if (typeof contig_lengths[data['data'][bin][index]] !== 'undefined') {
                        contigs.push(data['data'][bin][index]);
                        sum_contig_length += contig_lengths[data['data'][bin][index]];
                    }
                    
                }

                if (mode === 'manual' || mode === 'pan' || mode === 'server' || sum_contig_length >= threshold)
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

            $('#modSaveState').modal('show');
            if ($('#saveState_list').val() === null) {
                $('#saveState_name').val('default');
            } else {
                $('#saveState_list').trigger('change');
            }
        }
    });
}

function saveState() 
{
    var name = $('#saveState_name').val();

    if (name.length==0) {
        $('#saveState_name').focus();
        return;
    }

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
                toastr.error("Failed, Interface running in read only mode.");
            }
            else if (response['status_code']==1)
            {
                // successfull
                $('#modSaveState').modal('hide');

                current_state_name = name;
                $('#current_state').html('[current state: ' + current_state_name + ']');
                toastr.success("State '" + current_state_name + "' successfully saved.");
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

            $('#modLoadState').modal('show');
        }
    });
}

function loadState()
{
    $('#modLoadState').modal('hide');
    var defer = $.Deferred();
    var state_name;
    if (autoload_state !== null)
    {
        state_name = autoload_state;
        autoload_state = null // to prevent load again.
    }
    else
    {
        if ($('#loadState_list').val() == null) {
            defer.reject();
            return;
        }

        state_name = $('#loadState_list').val();
    }

    waitingDialog.show('Requesting state data from the server ...', 
        {
            dialogSize: 'sm', 
            onHide: function() {
                defer.resolve(); 
            },
            onShow: function() {
                $.ajax({
                        type: 'POST',
                        cache: false,
                        url: '/state/get?timestamp=' + new Date().getTime(),
                        data: {'name': state_name },
                        success: function(response) {
                            try{
                                var state = JSON.parse(response);
                            }catch(e){
                                toastr.error('Failed to parse state data, ' + e);
                                defer.reject();
                                return;
                            }
                            
                            if ((state['version'] !== VERSION) || !state.hasOwnProperty('version'))
                            {
                                toastr.error("Version of the given state file doesn't match with version of the interactive tree, ignoring state file.");
                                defer.reject();
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

                                // add layers that not exist in state and exist in layerdata
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
                                $('#tree_type').val(state['tree-type']).trigger('change');
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
                            if (state.hasOwnProperty('order-by') && $("#trees_container option[value='" + state['order-by'] + "']").length) {
                                $('#trees_container').val(state['order-by']);
                            }
                            if (state.hasOwnProperty('current-view') && $("#views_container option[value='" + state['current-view'] + "']").length) {
                                $('#views_container').val(state['current-view']);
                            }
                            if (state.hasOwnProperty('max-font-size')) {
                                $('#max_font_size').val(state['max-font-size']);
                            }
                            if (state.hasOwnProperty('layer-margin'))
                                $('#layer-margin').val(state['layer-margin']);
                            if (state.hasOwnProperty('outer-ring-height'))
                                $('#outer-ring-height').val(state['outer-ring-height']);
                            if (state.hasOwnProperty('outer-ring-margin'))
                                $('#outer-ring-margin').val(state['outer-ring-margin']);
                            if (state.hasOwnProperty('edge-normalization'))
                                $('#edge_length_normalization').prop('checked', state['edge-normalization']);
                            if (state.hasOwnProperty('optimize-speed'))
                                $('#optimize_speed').prop('checked', state['optimize-speed']);
                            if (state.hasOwnProperty('custom-layer-margin')) {
                                $('#custom_layer_margin').prop('checked', state['custom-layer-margin']).trigger('change');
                            }
                            if (state.hasOwnProperty('grid-color')) {
                                $('#grid_color').attr('color', state['grid-color']);
                                $('#grid_color').css('background-color', state['grid-color']);
                            }
                            if (state.hasOwnProperty('grid-width')) {
                                $('#grid_width').val(state['grid-width']);
                            }
                            if (state.hasOwnProperty('bin-labels-font-size')) {
                                $('#bin_labels_font_size').val(state['bin-labels-font-size']);
                            }
                            if (state.hasOwnProperty('bin-labels-angle')) {
                                $('#bin_labels_angle').val(state['bin-labels-angle']);
                            }
                            if (state.hasOwnProperty('show-bin-labels')) {
                                $('#show_bin_labels').prop('checked', state['show-bin-labels']).trigger('change');
                            }
                            if (state.hasOwnProperty('autorotate-bin-labels')) {
                                $('#autorotate_bin_labels').prop('checked', state['autorotate-bin-labels']).trigger('change');
                            }
                            if (state.hasOwnProperty('show-grid-for-bins')) {
                                $('#show_grid_for_bins').prop('checked', state['show-grid-for-bins']).trigger('change');
                            }
                            if (state.hasOwnProperty('samples-edge-length-normalization')) {
                                $('#samples_edge_length_normalization').prop('checked', state['samples-edge-length-normalization']);
                            }
                            if (state.hasOwnProperty('samples-ignore-branch-length')) {
                                $('#samples_ignore_branch_length').prop('checked', state['samples-ignore-branch-length']);
                            }
                            if (state.hasOwnProperty('samples-tree-height')) {
                                $('#samples_tree_height').val(state['samples-tree-height']);
                            }

                            // reload layers
                            var current_view = $('#views_container').val();
                            $("#tbody_layers").empty();

                            if (state.hasOwnProperty('samples-categorical-colors'))
                                samples_categorical_colors = state['samples-categorical-colors']; 
                            if (state.hasOwnProperty('samples-stack-bar-colors'))
                                samples_stack_bar_colors = state['samples-stack-bar-colors']; 

                            if (state.hasOwnProperty('samples-order'))
                                $('#samples_order').val(state['samples-order']);

                            buildLayersTable(layer_order, views[current_view]);
                            buildSamplesTable(state['samples-layer-order'], state['samples-layers']);

                            current_state_name = state_name;
                            $('#current_state').html('[current state: ' + current_state_name + ']');

                            toastr.success("State '" + current_state_name + "' successfully loaded.");
                            waitingDialog.hide();
                        }
                    });
            },
        }
    );
    

    return defer.promise();
}
