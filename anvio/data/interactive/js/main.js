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


var VERSION = '3';
var LINE_COLOR='#888888';
var MONOSPACE_FONT_ASPECT_RATIO = 0.6;

var VIEWER_WIDTH;
var VIEWER_HEIGHT;

var scale = 0;
var drawer;

var samples_id_to_node_map;
var total_radius = 0;

var bins;
var clusteringData;

var layerdata;
var item_lengths;
var parameter_count;

var tree_type;
var layer_types;

var categorical_data_colors = {};
var categorical_stats = {};
var stack_bar_colors = {};
var stack_bar_stats = {};
var legends = [];

var layerdata_title = {};
var empty_tooltip = "";

var last_settings;

var search_column;
var search_results = [];

var views = {};
var layers = {};
var current_view = '';
var layer_order;

var current_state_name = "";

var collapsedNodes = [];

var session_id;
var mode;
var server_mode = false;
var samples_tree_hover = false;
var bbox;

var request_prefix = getParameterByName('request_prefix');
//---------------------------------------------------------
//  Init
//---------------------------------------------------------

$(window).resize(function() {
     // get current client size
    VIEWER_WIDTH = document.getElementById('svg').clientWidth || document.getElementById('svg').width.baseVal.value;
    VIEWER_HEIGHT = document.getElementById('svg').clientHeight || document.getElementById('svg').height.baseVal.value;
});

$(document).ready(function() {

    $.ajaxPrefilter(function(options) {
        if (request_prefix) {
            options.url = request_prefix + options.url;
            if (options.type.toLowerCase() == 'post')
            {
                options.data += '&csrfmiddlewaretoken=' + getCookie('csrftoken');
            }
        }
        return options;
    });

    $(window).trigger('resize');
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

    // initialize colorpicker for search result highlight color.
    $('#picker_highlight').colpick({
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

    $('#grid_color').colpick({
        layout: 'hex',
        submit: 0,
        colorScheme: 'light',
        onChange: function(hsb, hex, rgb, el, bySetColor) {
            $(el).css('background-color', '#' + hex);
            $(el).attr('color', '#' + hex);

            if (!bySetColor) $(el).val(hex);
        },
        onHide: function() {
            emit('bin-settings-changed');
        }
    }).keyup(function() {
        $(this).colpickSetColor(this.value);
    });

    $("li[role='presentation']").click(function (e) {
        if ($(this).hasClass('disabled')) {
            e.preventDefault();
            return false;
        }  
    });

    if (!$.browser.chrome)
    {
        toastr.warning("We tested anvi'o only on Google Chrome, and it seems you are using a different browser.\
                        For the best performance, and to avoid unexpected issues, please consider using anvi'o with\
                        the lastest version of Chrome.", "", { 'timeOut': '0', 'extendedTimeOut': '0' });
    }

    initData();
});

function initData() {
    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/init',
        success: function(response) {
            mode = response.mode;
            server_mode = response.server_mode;
            switchUserInterfaceMode(response.project, response.title);
            setupDescriptionPanel(response.description);

            document.title = response.title;
            $('#title-panel-first-line').text(response.title);

            session_id = response.session_id;
            if (response.check_background_process) {
                setTimeout(checkBackgroundProcess, 5000);
            }

            if(!response.inspection_available){
                toastr.info("Inspection of data items is not going to be available for this project.");
                $('.menuItemInspect').addClass('menu-disabled');
            }

            if(!response.sequences_available && mode != "collection" && mode != "pan"){
                toastr.info("No sequence data is available. Some menu items will be disabled.");
                $('.menuItemSequence').addClass('menu-disabled');
            }

            if (response.read_only)
            {
                toastr.info("It seems that this is a read-only instance, therefore the database-writing \
                            functions will be inaccessible.", "", { 'timeOut': '0', 'extendedTimeOut': '0' });

                $('[disabled-in-read-only=true]').addClass('disabled').prop('disabled', true);
            }

            item_lengths = response.item_lengths;

            var default_tree  = response.item_orders[0];
            var available_trees = response.item_orders[2];
            $('#trees_container').append(getComboBoxContent(default_tree, available_trees));
            clusteringData = response.item_orders[1]['data'];
            loadOrderingAdditionalData(response.item_orders[1]);


            var default_view = response.views[0];
            var available_views = response.views[2];
            $('#views_container').append(getComboBoxContent(default_view, available_views));

            $("#tbody_layers").sortable({helper: fixHelperModified, handle: '.drag-icon', items: "> tr"}).disableSelection(); 
            $("#tbody_samples").sortable({helper: fixHelperModified, handle: '.drag-icon', items: "> tr"}).disableSelection(); 
                        
            let merged = samplesMergeStackbarLayers(response.layers_information, response.layers_information_default_order);
            
            samples_order_dict = response.layers_order;
            samples_information_dict = merged['dict'];
            let samples_information_default_layer_order = merged['default_order'];

            let samples_groups = Object.keys(samples_information_dict).sort();

            samples_groups.forEach(function (group_name) {
                $('#sample_groups_container').append(`
                    <div style="float: left; padding: 4px 4px;">
                        <input type="checkbox" onclick="toggleSampleGroups();" id="group_${group_name}" value="${group_name}" ${group_name == 'default' ? 'checked="checked"' : ''}>
                        <label style="margin-left: 2px;" onclick="toggleSampleGroups();" for="group_${group_name}">${group_name}</label>
                    </div>`);
            });

            let available_orders = Object.keys(samples_order_dict).sort();
            $('#samples_order').append(new Option('custom', 'custom'));
            available_orders.forEach(function(order)
            {
                var order_name = order;
                if (samples_order_dict[order]['newick'] != null && samples_order_dict[order]['newick'] != '')
                    order_name += " (tree)";

                $('#samples_order').append(new Option(order_name, order));
            });
            $('#samples_order').val('custom').trigger('change');

            // Populate function sources checkboxes in search functions panel.
            if (response.functions_sources.length > 0) {
                $('#functions_sources_list').empty();
            }
            
            response.functions_sources.forEach((source) => {
                $('#functions_sources_list').append(`<label style="margin: 5px;"><input type="checkbox" value="${source}" checked="checked" style="margin-right: 2px;"/>${source}</label>`);
            });
            
            buildSamplesTable(convert_samples_order_to_array(samples_information_default_layer_order));
            toggleSampleGroups();
            changeViewData(response.views[1]);

            if (response.state[0] && response.state[1]) {
                processState(response.state[0], response.state[1]);
            }

            $('.loading-screen').hide();

            bins = new Bins(response.bin_prefix, document.getElementById('tbody_bins'));
            bins.NewBin();

            if (response.autodraw)
            {
                $('#btn_draw_tree').removeClass('glowing-button');

                $.when()
                 .then(drawTree)
                 .then(function() {
                    if (response.collection !== null && mode !== 'refine' && mode !== 'gene')
                    {
                        bins.ImportCollection(response.collection);
                    }

                    if ($('#panel-left').is(':visible')) {
                        setTimeout(toggleLeftPanel, 500);
                    }
                 });
            } 
        }
    });
}

function switchUserInterfaceMode(project, title) {
    if (server_mode == false && (mode == 'pan' || mode == 'gene' || mode == 'full')) {
        $('#search_functions_button').attr('disabled', false);
        $('#searchFunctionsValue').attr('disabled', false);
        $('.functions-not-available-message').hide();
    }

    // hide all mode dependent divs:
    $('.full-mode, .pan-mode, .collection-mode, .manual-mode, .server-mode, .refine-mode').hide();

    console.log("The running mode for the interface: " + mode);

    $('.' + mode + '-mode').show();
    $('.nav-tabs').css('background-image', 'url(images/' + mode + '-bg.png)');

    if (mode == 'pan') {
        $('#completion_title').attr('title', 'Gene Clusters').html('Gene Clusters');
        $('#redundancy_title').attr('title', 'Gene Calls').html('Gene Calls');
        $('#splits_title').hide();
        $('#len_title').hide();
        $('.gene-filters-not-available-message').hide();
        $('.pan-filters button,input:checkbox').removeAttr('disabled')
    }

    if (server_mode) {
        $('.server-mode').show();
        $('.nav-tabs').css('background-image', 'url(images/server-bg.png)');
        $('#multiUser').show();
        $('#multiUser > span').html('<b>' + title + '</b><br /><i>(by <a href="/' + project.username + '" target="_blank">' + project.fullname + '</a>)</i>');
        $('#multiUser > img').attr('src', project.user_avatar);
        $('#multiUser > .download-button').attr('href', project.download_zip_url);
        $('#sidebar').css('margin-top', '81px');
        $('.upload-button').hide();
    }
}

function setupDescriptionPanel(description) {  
    $('#description-editor').val(description);
    $('#description-editor').markdown({
        'onShow': function (e) {
            $('[data-handler="bootstrap-markdown-cmdPreview"]').trigger('click');
        },
        'hiddenButtons': ['cmdUrl', 'cmdImage', 'cmdCode', 'cmdQuote'],
        'parser': function(content) {
            return renderMarkdown(content);
        },
        'additionalButtons': [
          [{
            data: [{
              name: 'cmdSave',
              title: 'Save',
              btnText: 'Save',
              btnClass: 'btn btn-success btn-sm',
              icon: {
                'glyph': 'glyphicon glyphicon-floppy-save',
              },
              callback: function(e) {
                $.ajax({
                        type: 'POST',
                        cache: false,
                        url: '/store_description',
                        data: {description: e.getContent()},
                        success: function(data) {
                            toastr.info("Description successfully saved to database.");
                        }
                    });
              }
            }]
          }]
        ],
        'fullscreen': {'enable': false},
    });

    if (description.length > 100) {
        toggleRightPanel('#description-panel');
    }
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
                    url: '/data/view/' + $('#views_container').val(),
                    success: function(data) {
                        changeViewData(data);
                        waitingDialog.hide();
                    }
                });
            },
        });

    return defer.promise();
}


function changeViewData(view_data) {
    layerdata = mergeStackbarLayers(view_data);
    parameter_count = layerdata[0].length;

    // since we are painting parent layers odd-even, 
    // we should remove single parents (single means no parent)
    removeSingleParents(); // in utils.js

    layer_order = Array.apply(null, Array(parameter_count-1)).map(function (_, i) {return i+1;}); // range(1, parameter_count)
    layer_types = {};

    // add layerdata columns to search window
    $('#searchLayerList').empty();
    for (var i=0; i < parameter_count; i++)
    {
        if (i == 0) {
            $('#searchLayerList').append(new Option("Item Name", i));
        } else {
            $('#searchLayerList').append(new Option(getPrettyName(layerdata[0][i]),i));
        }
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
    populateColorDicts();
    buildLegendTables();
}


function mergeStackbarLayers(view_data) {
    layerdata = []

    for (let i=0; i < view_data.length; i++) {
        layerdata.push([]);
    }

    for (let i=0; i < view_data[0].length; i++) {
        let pos = -1;
        let stack_group_name; // part before !
        let stack_layer_name; // after !

        if (view_data[0][i].indexOf('!') > -1) {
            stack_group_name = view_data[0][i].split('!')[0];
            stack_layer_name = view_data[0][i].split('!')[1];

            for (let j=0; j < layerdata[0].length; j++) {
                if (layerdata[0][j].startsWith(stack_group_name + '!')) {
                    pos = j;
                    break;
                }
            }
        }

        for (let j=0; j < view_data.length; j++) {
            if (pos == -1) {
                layerdata[j].push(view_data[j][i]);
            } else {
                if (j==0) {
                    layerdata[0][pos] += ';' + stack_layer_name;
                }
                else {
                    layerdata[j][pos] += ';' + view_data[j][i];
                }
            }
        }
    }

    return layerdata;
}


function populateColorDicts() {
    for (var layer_id=0; layer_id < parameter_count; layer_id++)
    {
        let layer_name = layerdata[0][layer_id];

        if (layer_types[layer_id] == 1) {
            if (!(layer_id in stack_bar_colors))
            {
                stack_bar_colors[layer_id] = {};
                stack_bar_stats[layer_id] = {};

                var bars = (layer_name.indexOf('!') > -1) ? layer_name.split('!')[1].split(';') : layer_name.split(';');
                for (var j=0; j < bars.length; j++)
                {
                    stack_bar_colors[layer_id][bars[j]] = getNamedCategoryColor(bars[j]);
                    
                    let sum = 0;
                    for (let i=1; i < layerdata.length; i++) {
                        sum += parseFloat(layerdata[i][layer_id].split(';')[j]);
                    }

                    stack_bar_stats[layer_id][bars[j]] = sum;
                } 
            }
        }

        if (layer_types[layer_id] == 2) {
            if (!(layer_id in categorical_data_colors))
            {
                categorical_stats[layer_id] = {};
                categorical_data_colors[layer_id] = {};
                for (var i=1; i < layerdata.length; i++)
                {
                    var _category_name = layerdata[i][layer_id];
                    if (_category_name == null || _category_name == '' || _category_name == 'null')
                        _category_name = 'None';
                    layerdata[i][layer_id] = _category_name;

                    if (typeof categorical_data_colors[layer_id][_category_name] === 'undefined'){
                        categorical_data_colors[layer_id][_category_name]  = getNamedCategoryColor(_category_name);
                        categorical_stats[layer_id][_category_name] = 0;
                    }

                    categorical_stats[layer_id][_category_name]++;
                }
            }
        }
    }

    for (let group in samples_information_dict) {
        var first_sample = Object.keys(samples_information_dict[group])[0];

        if (typeof first_sample !== 'undefined')
        {
            for (let sample_layer_name in samples_information_dict[group][first_sample])
            {
                if (isNumber(samples_information_dict[group][first_sample][sample_layer_name]))
                {
                    // no color table for numeric
                }
                else if (sample_layer_name.indexOf(';') > -1) // stack bar
                {
                    if (typeof samples_stack_bar_colors[group] === 'undefined') {
                        samples_stack_bar_colors[group] = {};
                        samples_stack_bar_stats[group] = {};
                    }

                    if (!(sample_layer_name in samples_stack_bar_colors[group]))
                    {
                        samples_stack_bar_colors[group][sample_layer_name] = {};
                        samples_stack_bar_stats[group][sample_layer_name] = {};

                        for (var j=0; j < sample_layer_name.split(";").length; j++)
                        {
                            let item_name = sample_layer_name.split('!')[1].split(';')[j];

                            samples_stack_bar_colors[group][sample_layer_name][item_name] = randomColor();

                            let sum = 0;
                            for (let _sample in samples_information_dict[group]) {
                                sum += parseFloat(samples_information_dict[group][_sample][sample_layer_name].split(';')[j]);
                            }

                            samples_stack_bar_stats[group][sample_layer_name][item_name] = sum;
                        } 
                    }
                }
                else // categorical
                {
                    if (typeof samples_categorical_colors[group] === 'undefined') {
                        samples_categorical_colors[group] = {};
                        samples_categorical_stats[group] = {};
                    }

                    if (typeof samples_categorical_colors[group][sample_layer_name] === 'undefined') {
                        samples_categorical_colors[group][sample_layer_name] = {};
                        samples_categorical_stats[group][sample_layer_name] = {};

                        for (let _sample in samples_information_dict[group])
                        {
                            var _category_name = samples_information_dict[group][_sample][sample_layer_name];
                            if (_category_name == null || _category_name == '' || _category_name == 'null')
                                _category_name = 'None';
                            samples_information_dict[group][_sample][sample_layer_name] = _category_name;

                            if (typeof samples_categorical_colors[group][sample_layer_name][_category_name] === 'undefined'){
                                samples_categorical_colors[group][sample_layer_name][_category_name] = getNamedCategoryColor(_category_name);
                                samples_categorical_stats[group][sample_layer_name][_category_name] = 0;
                            }

                            samples_categorical_stats[group][sample_layer_name][_category_name]++;
                        }
                    }
                }
            }
        }
    }
}

function buildLegendTables() {
    if(typeof $('#legend_settings').data("ui-accordion") != "undefined"){
        $('#legend_settings').accordion("destroy");
        $('#legend_settings').empty();
    }
    
    legends = [];

    for (let pindex in categorical_data_colors)
    {
        var names = Object.keys(categorical_stats[pindex]).sort(function(a,b){return categorical_stats[pindex][b]-categorical_stats[pindex][a]});

        names.push(names.splice(names.indexOf('None'), 1)[0]); // move null and empty categorical items to end

        legends.push({
            'name': getPrettyName(getLayerName(pindex)),
            'source': 'categorical_data_colors',
            'key': pindex,
            'item_names': names,
            'item_keys': names,
            'stats': categorical_stats[pindex]
        });
    }

    for (pindex in stack_bar_colors)
    {
        var layer_name = getLayerName(pindex);
        var names = (layer_name.indexOf('!') > -1) ? layer_name.split('!')[1].split(';') : layer_name.split(';');
        var pretty_name = getLayerName(pindex);
        pretty_name = (pretty_name.indexOf('!') > -1) ? pretty_name.split('!')[0] : pretty_name;

        legends.push({
            'name': getPrettyName(pretty_name),
            'source': 'stack_bar_colors',
            'key': pindex,
            'item_names': names,
            'item_keys': names,
            'stats': stack_bar_stats[pindex]
        });    
    }

    for (let group in samples_categorical_colors) {
        for (let sample in samples_categorical_colors[group])
        {
            var names = Object.keys(samples_categorical_colors[group][sample]);

            legends.push({
                'name': group + ' :: ' + getPrettyName(sample),
                'source': 'samples_categorical_colors',
                'group': group,
                'key': sample,
                'item_names': names,
                'item_keys': names,
                'stats': samples_categorical_stats[group][sample]
            });
        }
    }

    for (let group in samples_stack_bar_colors) {
        for (let sample in samples_stack_bar_colors[group])
        {
            var names = (sample.indexOf('!') > -1) ? sample.split('!')[1].split(';') : sample.split(';');
            var pretty_name = (sample.indexOf('!') > -1) ? sample.split('!')[0] : sample;

            legends.push({
                'name': group + ' :: ' + getPrettyName(pretty_name),
                'source': 'samples_stack_bar_colors',
                'group': group,
                'key': sample,
                'item_names': names,
                'item_keys': names,
                'stats': samples_stack_bar_stats[group][sample]
            });
        }
    }

    for (var i=0; i < legends.length; i++)
    {
        var legend = legends[i];
        var template = '<span>';

        if (legends[i]['source'].indexOf('samples') > -1) {
            template += '<span class="label label-default">Samples</span> '
        } else {
            template += '<span class="label label-default">Main</span> '
        }

        template += legend['name'] + '</span><div>';
        template += `Sort: <div class="btn-group" role="group">
                        <button type="button" class="btn btn-default" onClick="orderLegend(` + i + `, 'alphabetical');"><span class="glyphicon glyphicon-sort-by-alphabet"></span> Alphabetical</button>
                        <button type="button" class="btn btn-default" onClick="orderLegend(` + i + `, 'count');"><span class="glyphicon glyphicon-sort-by-order-alt"></span> Count</button>
                    </div>
                    <div class="btn-group" role="group">
                        <button type="button" class="btn btn-default" style="margin-left: 10px;" onClick="$('#batch_coloring_` + i + `').slideToggle();"><span class="glyphicon glyphicon-tint"></span> Batch coloring</button>
                    </div>
                    <div id="batch_coloring_` + i + `"  style="display: none; margin: 10px;">
                        <table class="col-md-12 table-spacing">
                            <tr>
                                <td class="col-md-2">Rule: </td>
                                <td class="col-md-10">
                                    <input type="radio" name="batch_rule_`+i+`" value="all" checked> All <br />
                                    <input type="radio" name="batch_rule_`+i+`" value="name"> Name contains <input type="text" id="name_rule_`+i+`" size="8"><br />
                                    <input type="radio" name="batch_rule_`+i+`" value="count"> Count 
                                        <select id="count_rule_`+i+`">
                                            <option selected>==</option>
                                            <option>&lt;</option>
                                            <option>&gt;</option>
                                        </select>
                                        <input type="text" id="count_rule_value_`+i+`" size="3">
                                    <br />
                                </td>
                            </tr>
                            <tr>
                                <td class="col-md-2">Color: </td>
                                <td class="col-md-10"><div id="batch_colorpicker_`+i+`" class="colorpicker" color="#FFFFFF" style="margin-right: 5px; background-color: #FFFFFF; float: none; "></div></td>
                            </tr>
                            <tr>
                                <td class="col-md-2"></td>
                                <td class="col-md-10"><input id="batch_randomcolor_`+i+`" type="checkbox" /> Assign random color</td>
                            </tr>
                            <tr>
                                <td class="col-md-2"></td>
                                <td class="col-md-10"><button type="button" class="btn btn-default" onclick="batchColor(`+i+`);">Apply</button></td>
                            </tr>
                        </table>
                    </div>
                    <div style="clear: both; display:block;"></div>
                    <hr style="margin-top: 4px; margin-bottom: 4px; "/>`;

        template += '<div id="legend_content_' + i + '"></div>';
        template = template + '<div style="clear: both; display:block;"></div>';
        $('#legend_settings').append(template + '</div>');

        createLegendColorPanel(i); // this fills legend_content_X
    }

    $('#legend_settings, #search_tab_content').accordion({heightStyle: "content", collapsible: true});

    $('.colorpicker').colpick({
        layout: 'hex',
        submit: 0,
        colorScheme: 'light',
        onChange: function(hsb, hex, rgb, el, bySetColor) {
            $(el).css('background-color', '#' + hex);
            $(el).attr('color', '#' + hex);
        }
    });
}

function batchColor(legend_id) {
    var rule = $('[name=batch_rule_'+legend_id+']:checked').val()
    var color = $('#batch_colorpicker_' + legend_id).attr('color');
    var randomize_color = $('#batch_randomcolor_' + legend_id).is(':checked');
    var legend = legends[legend_id];

    for (var i=0; i < legend['item_keys'].length; i++)
    {
        if (randomize_color) {
            color = randomColor();
        }

        if (rule == 'all') {
            if (typeof legend['group'] === 'undefined') {
                window[legend['source']][legend['key']][legend['item_keys'][i]] = color;
            } else {
                window[legend['source']][legend['group']][legend['key']][legend['item_keys'][i]] = color;
            }
        }
        else if (rule == 'name') {
            if (legend['item_names'][i].toLowerCase().indexOf($('#name_rule_' + legend_id).val().toLowerCase()) > -1) {
                if (typeof legend['group'] === 'undefined') {
                    window[legend['source']][legend['key']][legend['item_keys'][i]] = color;
                } else {
                    window[legend['source']][legend['group']][legend['key']][legend['item_keys'][i]] = color;
                }
            }
        } 
        else if (rule == 'count') {
            if (eval("legend['stats'][legend['item_keys'][i]] " + unescape($('#count_rule_'+legend_id).val()) + " " + parseFloat($('#count_rule_value_'+legend_id).val()))) {
                if (typeof legend['group'] === 'undefined') {
                    window[legend['source']][legend['key']][legend['item_keys'][i]] = color;
                } else {
                    window[legend['source']][legend['group']][legend['key']][legend['item_keys'][i]] = color;
                }
            }
        }
    }

    createLegendColorPanel(legend_id);
}

function createLegendColorPanel(legend_id) {
    var legend = legends[legend_id];
    var template = '';

    for (var j = 0; j < legend['item_names'].length; j++) {

        var _name = legend['item_names'][j];

        if (legend.hasOwnProperty('group')) {
            var _color = window[legend['source']][legend['group']][legend['key']][legend['item_keys'][j]];
        } else {
            var _color = window[legend['source']][legend['key']][legend['item_keys'][j]];
        }

        if (legend.hasOwnProperty('stats') && legend['stats'][_name] == 0) {
            continue;
        }

        if (legend['source'].indexOf('stack') > -1) {
            _name = _name.replace('Unknown_t_', '').replace('_', ' ') + ' <span title="' + legend['stats'][_name] + '">(Total: ' + readableNumber(legend['stats'][_name]) + ')</span>';
        } else {
            _name = _name + ' (' + legend['stats'][_name] + ')';
        }

        template = template + '<div style="float: left; width: 50%; display: inline-block; padding: 3px 5px;">' + 
                                '<div class="colorpicker legendcolorpicker" color="' + _color + '"' +
                                'style="margin-right: 5px; background-color: ' + _color + '"' +
                                'callback_source="' + legend['source'] + '"' +
                                'callback_group="' + ((typeof legend['group'] !== 'undefined') ? legend['group'] : '') + '"' +
                                'callback_pindex="' + legend['key'] + '"' +
                                'callback_name="' + legend['item_keys'][j] + '"' + 
                               '></div>' + _name + '</div>';
    }

    $('#legend_content_' + legend_id).empty();
    $('#legend_content_' + legend_id).html(template);

    $('.legendcolorpicker').colpick({
        layout: 'hex',
        submit: 0,
        colorScheme: 'light',
        onChange: function(hsb, hex, rgb, el, bySetColor) {
            $(el).css('background-color', '#' + hex);
            if (el.getAttribute('callback_group') !== '') {
                window[el.getAttribute('callback_source')][el.getAttribute('callback_group')][el.getAttribute('callback_pindex')][el.getAttribute('callback_name')] = '#' + hex;
            } else {
                window[el.getAttribute('callback_source')][el.getAttribute('callback_pindex')][el.getAttribute('callback_name')] = '#' + hex;
            }
        }
    });
}

function orderLegend(legend_id, type) {
    if (type == 'alphabetical') {
        legends[legend_id]['item_names'] = legends[legend_id]['item_names'].sort();
    }
    else if (type == 'count') {
        legends[legend_id]['item_names'] = legends[legend_id]['item_names'].sort(function(a,b){return legends[legend_id]['stats'][b]-legends[legend_id]['stats'][a]});
    }

    // in both cases we will move None to end.
    legends[legend_id]['item_names'].push(legends[legend_id]['item_names'].splice(legends[legend_id]['item_names'].indexOf('None'), 1)[0]);

    createLegendColorPanel(legend_id);
}

function loadOrderingAdditionalData(order) {
    collapsedNodes = [];
    
    if (order.hasOwnProperty('additional')) {
        let orders_additional = JSON.parse(order['additional']);

        if (orders_additional.hasOwnProperty('collapsedNodes')) {
            collapsedNodes = orders_additional['collapsedNodes'];
        }
    }
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
                    url: '/tree/' + $('#trees_container').val(),
                    success: function(order) {
                        clusteringData = order['data'];
                        loadOrderingAdditionalData(order);

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
    available_items = available_items.sort()
    var combo = '';
    var combo_item = '<option value="{val}"{sel}>{text}</option>';

    $.each(available_items, function(index, item) {
        let text_val;
        if (item.indexOf(':') == -1) {
            text_val = getPrettyName(item);
        } else {
            text_val = getClusteringPrettyName(item);
        }

        if(item == default_item)
        {
            combo += combo_item
                        .replace('{val}', item)
                        .replace('{sel}', ' selected')
                        .replace('{text}', text_val);
        }
        else
        {
            combo += combo_item
                        .replace('{val}', item)
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
            }

            if (hasViewSettings)
            {
                var norm = view_settings['normalization'];
            }
            else
            {
                var norm = (mode == 'full') ? 'log' : 'none';
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
                var norm   = getNamedLayerDefaults(layer_name, 'norm', (mode == 'full' || mode == 'refine') ? 'log' : 'none');
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
                '        <option value="line"{option-type-line}>Line</option>' +
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

        $('#tbody_layers .input-height:last').change(function (ev) {
            // setting height 0 changes samples order to custom, only if layer is in samples order
            if (ev.target.value == 0) {
                var layer_name = $(ev.target).parent().parent().find('td:nth(1)').attr('title');
                var layer_names_in_samples = null;

                if (samples_order_dict.hasOwnProperty($('#samples_order').val())) {
                    var samples_organization = samples_order_dict[$('#samples_order').val()];

                    if (samples_organization['basic'] != null && samples_organization['basic'] != "")
                    {
                        layer_names_in_samples = samples_organization['basic'].split(',');
                    }
                    else
                    {
                        layer_names_in_samples = get_newick_leaf_order(samples_organization['newick']);
                    }

                    if (layer_names_in_samples.indexOf(layer_name) > -1) {
                        samplesClusteringData = {'newick': '', 'basic': null};
                        $('#samples_order').val('custom').trigger('change');
                    }
                }
            }
        });

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
            colorScheme: 'light',
            onChange: function(hsb, hex, rgb, el, bySetColor) {
                $(el).css('background-color', '#' + hex);
                $(el).attr('color', '#' + hex);

                if (!bySetColor) $(el).val(hex);
            }
        }).keyup(function() {
            $(this).colpickSetColor(this.value);
        });

        $('#table_layers .drag-icon').on('mousedown', function() {
            samplesClusteringData = {'newick': '', 'basic': null};
            $('#samples_tree_modified_warning').hide();
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
    state['background-opacity'] = $('#background_opacity').val();
    state['max-font-size-label'] = $('#max_font_size_label').val();
    state['draw-guide-lines'] = $('#draw_guide_lines').val();
    
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
            var samples_group_name = $(tr).attr('samples-group-name');

            state['samples-layer-order'].push({
                'layer_name': samples_layer_name,
                'group': samples_group_name
            });

            if (!state['samples-layers'].hasOwnProperty(samples_group_name)) {
                state['samples-layers'][samples_group_name] = {};
            }

            state['samples-layers'][samples_group_name][samples_layer_name] = {
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

    state['samples-groups'] = {};
    $('#sample_groups_container input:checkbox').each((index, checkbox) => {
        state['samples-groups'][$(checkbox).val()] = $(checkbox).is(':checked');
    });

    return state;
}

function drawTree() {
    var defer = $.Deferred();
    var settings = serializeSettings();
    tree_type = settings['tree-type'];

    $('#btn_draw_tree').removeClass('glowing-button');
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
            onHide: function() {
                defer.resolve(); 
            },
            onShow: function() {
                drawer = new Drawer(settings);
                drawer.draw();

                // last_settings used in export svg for layer information,
                // we didn't use "settings" sent to draw_tree because draw_tree updates layer's min&max
                last_settings = serializeSettings();

                bins.RedrawBins();

                $('#btn_draw_tree').prop('disabled', false);
                $('#btn_redraw_samples').prop('disabled', false);

                if (settings['tree-radius'] == 0)
                {
                    $('#tree-radius-container').show();
                    $('#tree-radius').val(Math.max(VIEWER_HEIGHT, VIEWER_WIDTH));
                }
                waitingDialog.hide();
            },
        });

    return defer.promise();
}

function showContigNames(bin_id, updateOnly) {
    if (typeof updateOnly === 'undefined')
        updateOnly = false;

    var title = 'Splits in "' + $('#bin_name_' + bin_id).val() + '"';

    if (updateOnly && !checkObjectExists('#modal' + title.hashCode()))
        return;

    var msg = '<table class="table table-striped">';
    for (const label of bins.GetBinNodeLabels(bin_id)) {
        msg += `<tr><td><a href='#' class='no-link' onclick="bins.HighlightItems('${label}');">${label}</a></td></tr>`;
    }

    msg = msg + '</table>';

    showDraggableDialog(title, msg, updateOnly);
}

function showGenSummaryWindow() {
    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/collections',
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

function showCompleteness(bin_id, updateOnly) {
    if (typeof updateOnly === 'undefined')
        updateOnly = false;

    if (!bins.cache['completeness'].hasOwnProperty(bin_id))
        return;

    var refs = bins.cache['completeness'][bin_id]['refs'];
    var stats = bins.cache['completeness'][bin_id]['stats'];
    var averages = bins.cache['completeness'][bin_id]['averages'];

    var title = 'Completeness of "' + $('#bin_name_' + bin_id).val() + '"';

    if (updateOnly && !checkObjectExists('#modal' + title.hashCode()))
        return;

    var msg = '<table class="table table-striped sortable">' +
        '<thead><tr><th data-sortcolumn="0" data-sortkey="0-0">Source</th><th data-sortcolumn="1" data-sortkey="1-0">SCG domain</th><th data-sortcolumn="2" data-sortkey="2-0">Percent completion</th></tr></thead><tbody>';

    for (let source in stats){
        if(stats[source]['domain'] != averages['domain'])
            // if the source is not matching the best domain recovered
            // don't show it in the interface
            continue;

        msg += "<tr><td data-value='" + source  + "'><a href='" + refs[source] + "' class='no-link' target='_blank'>" + source + "</a></td><td data-value='" + stats[source]['domain'] + "'>" + stats[source]['domain'] + "</td><td data-value='" + stats[source]['percent_completion'] + "'>" + stats[source]['percent_completion'].toFixed(2) + "%</td></tr>";
    }

    msg = msg + '</tbody></table>';

    showDraggableDialog(title, msg, updateOnly);
}

function showRedundants(bin_id, updateOnly) {
    if (typeof updateOnly === 'undefined')
        updateOnly = false;
    
    if (!bins.cache['completeness'].hasOwnProperty(bin_id))
        return;

    var stats = bins.cache['completeness'][bin_id]['stats'];
    var averages = bins.cache['completeness'][bin_id]['averages'];

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
        var split_array_all = '';

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

                split_array_all += split_array;
            }

            redundants_html += '<span style="cursor:pointer;" \
                                    data-toggle="tooltip" data-placement="top" title="' + title + '" \
                                    onclick="bins.HighlightItems([' + split_array + ']);"> \
                                    ' + redundant + ' (' + stats[source]['redundants'][redundant].length + ') \
                                  </span><br />';
        }

        tabletext += '<tr><td valign="top">' + redundants_html + '<br /><br /><span style="cursor:pointer;" \
                                    onclick="bins.HighlightItems([' + split_array_all + ']);">(Highlight All)\
                                  </span></tr></td></table></div>';
        output += tabletext;

        if (oddeven%2==0)
        {
            output += '</div><div class="col-md-12">'
        }
    }

    output += '</div>';

    showDraggableDialog(output_title, output, updateOnly);
}

function exportSvg(dontDownload) {
    if (!drawer) 
        return;

    // draw bin and layer legend to output svg
    var settings = serializeSettings();

    var bins_to_draw = new Array();
    $('#tbody_bins tr').each(
        function(index, bin) {
            var bin_id = $(bin).attr('bin-id');

            var _bin_info = {
                'name': $('#bin_name_' + bin_id).val(),
                'color': $('#bin_color_' + bin_id).attr('color'),
            };

            if (mode == 'pan') {
                _bin_info['gene_clusters'] = $('#completeness_' + bin_id).val(); 
                _bin_info['gene-calls'] = $('#redundancy_' + bin_id).val(); 
            } else {
                _bin_info['contig-length'] = $('#contig_length_' + bin_id).html();
                _bin_info['contig-count'] = $('#contig_count_' + bin_id).val();
            }
            
            bins_to_draw.push(_bin_info);
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
    drawTitle(last_settings);
    drawLegend();

    var svg = document.getElementById('svg');
    var viewBox = svg.getBBox();
    svg.setAttribute('viewBox', viewBox['x'] + " " + viewBox['y'] + " " + viewBox['width'] + " " + viewBox['height'])

    if (dontDownload == true) {
        return;
    }

    svgCrowbar();

    svg.removeAttribute('viewBox');
    $('#samples_tree').prepend(detached_clones);
    $('#bin_legend').remove();
    $('#layer_legend').remove();
    $('#title_group').remove();
    $('#legend_group').remove();
}


function storeRefinedBins() {
    let collection_info = bins.ExportCollection();

    $.ajax({
        type: 'POST',
        cache: false,
        url: '/data/store_refined_bins',
        data: { 
            data: JSON.stringify(collection_info['data'], null, 4), 
            colors: JSON.stringify(collection_info['colors'], null, 4)
        },
        success: function(data) {
            if (data.status == -1){
                toastr.error(data.message, "You made the server upset :(");
            } else {
                toastr.info(data.message, "The server is on board");
            }
        }
    });
}


function generateSummary() {
    var collection = $('#summaryCollection_list').val();

    if (collection === null)
        return;

    waitingDialog.show('Generating summary...', {dialogSize: 'sm'});

    $.ajax({
        type: 'GET',
        cache: false,
        url: '/summarize/' + collection,
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

function showSaveStateWindow()
{
    $.ajax({
        type: 'GET',
        cache: false,
        url: '/state/all',
        success: function(state_list) {
            $('#saveState_list').empty();

            for (let state_name in state_list) {
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

function showGeneratePhylogeneticTreeWindow() {
    $('#phylogeny_gene_cluster').empty();
    $('#phylogeny_programs').empty();
    $('#modPhylogeneticTree :input').attr("disabled", false);
    $('.generating-tree').hide();
    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/phylogeny/programs',
        success: function(available_programs) {
              $.ajax({
                type: 'GET',
                cache: false,
                url: '/data/phylogeny/aligners',
                success: function(available_aligners) {
                    $('#phylogeny_programs').empty();
                    for (var i=0; i < available_programs.length; i++) {
                        if (available_programs[i] == 'default')
                            continue;

                        $('#phylogeny_programs').append(new Option(available_programs[i]))
                    }

                    $('#phylogeny_aligners').empty();
                    for (var i=0; i < available_aligners.length; i++) {
                        if (available_aligners[i] == 'default')
                            continue;

                        $('#phylogeny_aligners').append(new Option(available_aligners[i]))
                    }

                    $('#tbody_bins tr').each(
                        function(index, bin) {
                            var bin_id = $(bin).attr('bin-id');
                            var bin_name = $('#bin_name_' + bin_id).val();

                            $('#phylogeny_gene_cluster').append('<option value="' + bin_id + '">' + bin_name + '</option>');
                        }
                    );
                    $('#modPhylogeneticTree').modal('show');
                }});
        }});
}

function generatePhylogeneticTree() {
    var new_phylogeny_name = $('#phylogeny_name').val();
    var gene_cluster_list = [];
    var gene_clusters_id = $('#phylogeny_gene_cluster').val();
    
    for (const node of this.selections[gene_clusters_id].values()) {
        if (node.IsLeaf()) {
            gene_cluster_list.push(node.label);
        } 
    }

    if (gene_cluster_list.length == 0) {
        alert("The Bin you selected does not contain any gene_clusters.");
        return;
    }

    if (samples_order_dict.hasOwnProperty(new_phylogeny_name)) {
        alert("The name '" + new_phylogeny_name + "' already exists, please give another name. ");
        return;
    }

    $('#modPhylogeneticTree :input').attr("disabled", true);
    $('.generating-tree').show();
    $.ajax({
        type: 'POST',
        cache: false,
        url: '/data/phylogeny/generate_tree',
        data: {
            'name': $('#phylogeny_name').val(),
            'program': $('#phylogeny_programs').val(),
            'aligner': $('#phylogeny_aligners').val(),
            'gene_clusters': gene_cluster_list,
            'store_tree': $('#phylogeny_store_generated_tree').is(':checked'),
        },
        success: function(response) {
            if (response['status'] != 0) {
                alert(response['message']);
                showGeneratePhylogeneticTreeWindow();
                return;
            } else {
                samples_order_dict[$('#phylogeny_name').val()] = {'basic': '', 'newick': response['tree']};
                $('#samples_order').append('<option value="'+ new_phylogeny_name + '">' + new_phylogeny_name + '</option>');
                $('#samples_order').val(new_phylogeny_name);
                $('#samples_order').trigger('change');
                $('#modPhylogeneticTree').modal('hide');
                drawTree();
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
        url: '/state/all',
        success: function(state_list) {
            for (let state_name in state_list) {
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
        url: '/state/save/' + name,
        data: {
            'content': JSON.stringify(serializeSettings(true), null, 4)
        },
        success: function(response) {
            if (typeof response != 'object') {
                response = JSON.parse(response);
            }

            if (response['status_code']==0)
            {
                toastr.error("Failed, Interface running in read only mode.");
            }
            else if (response['status_code']==1)
            {
                // successfull
                $('#modSaveState').modal('hide');

                current_state_name = name;
                toastr.success("State '" + current_state_name + "' successfully saved.");
            }
        }
    });
}

function showUploadProject() {
    $('#upload_state').empty();
    $('#upload_collection').empty();
    $('#upload_view').empty();
    $('#upload_ordering').empty();

    $('#trees_container option').each(function(index, option) {
        $(option).clone().appendTo('#upload_ordering');
    });

    $('#views_container option').each(function(index, option) {
        $(option).clone().appendTo('#upload_view');
    });

    $('#upload_state').append('<option selected>Select State</option>');
    $('#upload_collection').append('<option>Select Collection</option>');

    $('#upload_project_name').val($('#title-panel-first-line').text());

    $.ajax({
        type: 'GET',
        cache: false,
        url: '/state/all',
        success: function(state_list) {
            for (let state_name in state_list) {
                $('#upload_state').append('<option>' + state_name + '</option>');
            }

            $.ajax({
                type: 'GET',
                cache: false,
                url: '/data/collections',
                success: function(collection_list) {
                    for (collection_name in collection_list) {
                        $('#upload_collection').append('<option>' + collection_name + '</option>');
                    }
                }
            });
        }
    });
    $('#modUpload').modal('show');
}

function uploadProject() {
    $.ajax({
        type: 'POST',
        cache: false,
        url: '/upload_project',
        data: {
            username: $('#username').val(),
            password: $('#password').val(),
            project_name: $('#upload_project_name').val(),
            ordering: $('#upload_ordering').val(),
            view: $('#upload_view').val(),
            state: $('#upload_state').val(),
            collection: $('#upload_collection').val(),
            delete_if_exists: $('#upload_delete_if_exists').is(':checked'),
            include_samples: $('#upload_include_samples').is(':checked'),
            include_description: $('#upload_include_description').is(':checked')
        },
        success: function(data) {
            if (data['status'] == 1) {
                $('.upload-message').removeClass('alert-success').addClass('alert-danger');
                $('.upload-message').show();
                $('.upload-message').html(data['message']);
            } else {
                $('.upload-message').removeClass('alert-danger').addClass('alert-success');
                $('.upload-message').show();
                $('.upload-message').html("Project successfully uploaded, to view your projects click <a href='https://anvi-server.org/projects' target='_blank'>here.</a>");
            }
        }
    });
}


function showLoadStateWindow()
{
    $.ajax({
        type: 'GET',
        cache: false,
        url: '/state/all',
        success: function(state_list) {
            $('#loadState_list').empty();

            for (let state_name in state_list) {
                $('#loadState_list').append('<option lastmodified="' + state_list[state_name]['last_modified'] + '">' + state_name + '</option>');
            }

            $('#modLoadState').modal('show');
        }
    });
}

function loadState()
{
    var defer = $.Deferred();
    $('#modLoadState').modal('hide');
    if ($('#loadState_list').val() == null) {
        defer.reject();
        return;
    }

    var state_name = $('#loadState_list').val();
    waitingDialog.show('Requesting state data from the server ...', 
        {
            dialogSize: 'sm', 
            onHide: function() {
                defer.resolve(); 
            },
            onShow: function() {
                $.ajax({
                        type: 'GET',
                        cache: false,
                        url: '/state/get/' + state_name,
                        success: function(response) {
                            try{
                                clusteringData = response[1]['data'];
                                loadOrderingAdditionalData(response[1]);
                                changeViewData(response[2]);
                                processState(state_name, response[0]);
                            }catch(e){
                                console.error("Exception thrown", e.stack);
                                toastr.error('Failed to parse state data, ' + e);
                                defer.reject();
                                return;
                            }
                            waitingDialog.hide();
                        }
                    });
            },
        }
    );

    return defer.promise();
}

function processState(state_name, state) {
    if (!state.hasOwnProperty('version'))
    {
        toastr.error("Interface received a state without version information, it will be not loaded.");
        throw "";
    }

    if (state['version'] == '0.2.1') {
        // switch to numerical versioning instead semantic one.
        state['version'] = '1';
    }

    if (state['version'] != VERSION) {
        toastr.info(`Interface received a state at version ${state['version']} but the current version of the
            interface is ${VERSION}. Anvi'o will try to upgrade it automatically.`);

        state = migrate_state(state);
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
        for (let view_key in state['views'])
        {
            views[view_key] = {};
            for (let key in state['views'][view_key])
            {
                let layer_id = getLayerId(key);
                if (layer_id != -1)
                {
                    views[view_key][layer_id] = state['views'][view_key][key];
                }
            }
        }
    }

    if (state.hasOwnProperty('layers')) {
        layers = {};
        for (let key in state['layers'])
        {
            let layer_id = getLayerId(key);
            if (layer_id != -1)
            {
                layers[layer_id] = state['layers'][key];
            }
        }
    }

    if (state.hasOwnProperty('categorical_data_colors')) {
        for (let key in state['categorical_data_colors'])
        {
            let layer_id = getLayerId(key);
            if (layer_id != -1)
            {
                categorical_data_colors[layer_id] = state['categorical_data_colors'][key];
            }
        }
    }

    if (state.hasOwnProperty('stack_bar_colors')) {
        for (let key in state['stack_bar_colors'])
        {
            let layer_id = getLayerId(key);
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
    if (state.hasOwnProperty('max-font-size-label')) {
        $('#max_font_size_label').val(state['max-font-size-label']);
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
    if (state.hasOwnProperty('background-opacity')) {
        $('#background_opacity').val(state['background-opacity']);
    }
    if (state.hasOwnProperty('draw-guide-lines')) {
        $('#draw_guide_lines').val(state['draw-guide-lines'])
    }

    // reload layers
    var current_view = $('#views_container').val();
    $("#tbody_layers").empty();

    if (state.hasOwnProperty('samples-categorical-colors')) {
        for (let key in state['samples-categorical-colors']) {
            if (key in samples_categorical_colors) {
                samples_categorical_colors[key] = state['samples-categorical-colors'][key];
            } 
        }
    }
    if (state.hasOwnProperty('samples-stack-bar-colors')) {
        for (let key in state['samples-stack-bar-colors']) {
            if (key in samples_stack_bar_colors) {
                samples_stack_bar_colors[key] = state['samples-stack-bar-colors'][key];
            } 
        }
    }

    if (state.hasOwnProperty('samples-order') && $(`#samples_order option[value='${state['samples-order']}']`).length > 0) {
        $('#samples_order').val(state['samples-order']).trigger('change');
    }

    buildLayersTable(layer_order, views[current_view]);
    buildSamplesTable(state['samples-layer-order'], state['samples-layers']);

    if (state.hasOwnProperty('samples-groups')) {
        for (let group_name in state['samples-groups']) {
            let checkbox = $('input:checkbox#group_' + group_name);

            if (checkbox) {
                $(checkbox).prop('checked', state['samples-groups'][group_name]);
            }
        }
    }

    toggleSampleGroups();
    buildLegendTables();

    current_state_name = state_name;

    toastr.success("State '" + current_state_name + "' successfully loaded.");
}


function restoreOriginalTree(type) {
    if (type == 'samples') {
        samplesClusteringData = {'newick': '', 'basic': null};
        $('#samples_order').val('custom').trigger('change');
        $('#samples_tree_modified_warning').hide();
        drawTree();
        return;
    }

    $.when()
     .then(onTreeClusteringChange)
     .then(
        function() {
            collapsedNodes = [];
            $('#tree_modified_warning').hide();
            drawTree();
        }
    );
}
