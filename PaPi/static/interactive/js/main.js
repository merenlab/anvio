//--------------------------------------------------------------------------------------------------
//  Globals
//--------------------------------------------------------------------------------------------------

var VIEWER_WIDTH = window.innerWidth || document.documentElement.clientWidth || document.getElementsByTagName('body')[0].clientWidth;
var VIEWER_HEIGHT = window.innerHeight || document.documentElement.clientHeight || document.getElementsByTagName('body')[0].clientHeight;

var ZOOM_IN = 1.33;
var ZOOM_OUT = 0.75;

var LINE_COLOR='#888888';
var project_title;

var SCALE_MATRIX = 0;
var id_to_node_map = new Array();
var angle_per_leaf;
var height_per_leaf;
var tree_type;

var total_radius = 0;

var SELECTED = new Array();

var newick;
var available_trees;
var selected_tree_id;
var default_tree;

var metadata;
var contig_lengths;
var parameter_count;

var group_counter = 0; // for id
var group_count = 0;

var categorical_data_ids = new Array();
var categorical_data_colors = {};

var stack_bar_ids = new Array();
var stack_bar_colors = {};

var has_parent_layer = false;

var context_menu_target_id = 0;

var metadata_title = {};
var metadata_dict;

var metadata_swap_log = new Array();
var metadata_swap_log_reverse = new Array();

//---------------------------------------------------------
//  Init
//---------------------------------------------------------

$(document).ready(function() {

    // create "Group_1"
    newGroup();

    // common settings for dialogs
    $('.dialogs').dialog({
        resizable: false,
        width: 'auto',
        collapseEnabled: true,
        closeOnEscape: false,
        beforeclose: function(event, ui) {
            return false;
        },
        dialogClass: "noclose"
    });

    // settings per dialog
    $("#zoomDialog").dialog("option", "title", "Zoom").dialog("option", "position", {
            my: "right center",
            at: "right center",
            of: window
        });

    $("#treeControls").dialog("option", "title", "Tree Settings").dialog("option", "position", {
            my: "left bottom",
            at: "left center",
            of: window
        });

    $("#groups").dialog("option", "title", "Groups").dialog("option", "position", {
            my: "left top",
            at: "left center",
            of: window
        });

    $("#contig_name_dialog").dialog("option", "title", "Contig Names").dialog("option", "position", {
            my: "center",
            at: "center",
            of: window
        }).dialog('close');

    // get project title
    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/title?timestamp=' + new Date().getTime(),
        success: function(data) {
            if (data.length > 0)
            {
                document.title = data;
                project_title = data;
            }
        }
    });

    // get available clusterings (different newick trees)
    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/clusterings?timestamp=' + new Date().getTime(),
        success: function(data) {
        	default_tree = data[0];
        	available_trees = data[1];
        }
    });

    // load contig lengths
    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/contig_lengths?timestamp=' + new Date().getTime(),
        success: function(data) {
            contig_lengths = data;
        }
    });

    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/state?timestamp=' + new Date().getTime(),
        success: function(state) {
    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/meta?timestamp=' + new Date().getTime(),
        success: function(data) {
            state = eval(state);
            metadata = eval(data);

            parameter_count = metadata[0].length;

            // remove all single parents from metadata
            for (var i = 1; i < parameter_count; i++) 
            {
                if (metadata[0][i] == '__parent__') 
                {
                    has_parent_layer = true;

                    var parent_count_dict = {};
                    for (var j=1; j < metadata.length; j++)
                    {
                        if (metadata[j][i]=='')
                            continue;

                        if (typeof parent_count_dict[metadata[j][i]] === 'undefined')
                        {
                            parent_count_dict[metadata[j][i]] = 1;
                        }
                        else
                        {
                            parent_count_dict[metadata[j][i]]++;
                        }
                    }

                    $.each(parent_count_dict, function(parent_name, count)
                    {
                        if (count==1)
                        {
                            for (var j=1; j < metadata.length; j++)
                            {
                                if (metadata[j][i]==parent_name)
                                {
                                    metadata[j][i]='';
                                }
                            }
                        }
                    });
                }
            }
            // all clear

            if (jQuery.isEmptyObject(state)) {
                // state is empty, build ui using metadata

                categorical_data_ids = [];
                separated_data_ids = [];
                
                for (var i = 1; i < parameter_count; i++) {
                    if (metadata[0][i] == '__parent__') // parent
                    {
                        var parent_row_str = '<tr style="display: none">' +
                            '<td></td>' +
                            '<td></td>' +
                            '<td></td>' +
                            '<td><input class="input-height" type="text" size="2" id="height{id}" value="30"></input></td>' +
                            '</tr>';

                        parent_row_str = parent_row_str.replace(new RegExp('{id}', 'g'), i);

                        $('#tbody_layers').prepend(parent_row_str);
                    }
                    else if (metadata[1][i].indexOf(';') > -1) // stack data
                    {
                        stack_bar_ids.push(i);
                        stack_bar_colors[i] = new Array();

                        for (var j=0; j < metadata[1][i].split(";").length; j++)
                        {
                            stack_bar_colors[i].push('#000000');
                        }

                        var stack_bar_row_str = '<tr>' +
                            '<td><img class="drag-icon" src="images/drag.gif" /></td>' +
                            '<td title="' + metadata[0][i] + '">' + ((metadata[0][i].length > 10) ? metadata[0][i].slice(0,10) + "..." : metadata[0][i]) + '</td>' +
                            '<td>n/a</td>' +
                            '<td>' +
                            '    <select id="normalization{id}">' +
                            '        <option value="none">none</option>' +
                            '        <option value="sqrt">Square root</option>' +
                            '        <option value="log" selected>Logarithm</option>' +
                            '    </select>' +
                            '</td>' +
                            '<td><input class="input-height" type="text" size="3" id="height{id}" value="150"></input></td>' +
                            '<td>n/a</td>' +
                            '<td>n/a</td>' +
                            '</tr>';

                        stack_bar_row_str = stack_bar_row_str.replace(new RegExp('{id}', 'g'), i);

                        $('#tbody_layers').append(stack_bar_row_str);
                    }
                    else if (metadata[1][i] === '' || !isNumber(metadata[1][i])) // categorical data
                    { 
                        categorical_data_ids.push(i);
                        categorical_data_colors[i] = {};

                        var categorical_data_row_str = '<tr>' +
                            '<td><img class="drag-icon" src="images/drag.gif" /></td>' +
                            '<td title="' + metadata[0][i] + '">' + ((metadata[0][i].length > 10) ? metadata[0][i].slice(0,10) + "..." : metadata[0][i]) + '</td>' +
                            '<td>n/a</td>' +
                            '<td>n/a</td>' +
                            '<td><input class="input-height" type="text" size="3" id="height{id}" value="30"></input></td>' +
                            '<td>n/a</td>' +
                            '<td>n/a</td>' +
                            '</tr>';

                        categorical_data_row_str = categorical_data_row_str.replace(new RegExp('{id}', 'g'), i);

                        $('#tbody_layers').append(categorical_data_row_str);
                    } 
                    else // numerical data
                    { 
                        var numerical_data_row_str = '<tr>' +
                            '<td><img class="drag-icon" src="images/drag.gif" /></td>' +
                            '<td title="' + metadata[0][i] + '">' + ((metadata[0][i].length > 10) ? metadata[0][i].slice(0,10) + "..." : metadata[0][i]) + '</td>' +
                            '<td><div id="picker{id}" class="colorpicker"></td>' +
                            '<td>' +
                            '    <select id="normalization{id}" onChange="clearMinMax(this)">' +
                            '        <option value="none">none</option>' +
                            '        <option value="sqrt">Square root</option>' +
                            '        <option value="log" selected>Logarithm</option>' +
                            '    </select>' +
                            '</td>' +
                            '<td><input class="input-height" type="text" size="3" id="height{id}" value="150"></input></td>' +
                            '<td><input class="input-min" type="text" size="4" id="min{id}" value="0" disabled></input></td>' +
                            '<td><input class="input-max" type="text" size="4" id="max{id}" value="0" disabled></input></td>' +
                            '</tr>';

                        numerical_data_row_str = numerical_data_row_str.replace(new RegExp('{id}', 'g'), i);

                        $('#tbody_layers').append(numerical_data_row_str);
                    }
                }

                // TREES COMBO
                var available_trees_combo = '';
                var available_trees_combo_item = '<option value="{val}"{sel}>{text}</option>';
                
                $.each(available_trees, function(index, value) {
                	if(index == default_tree)
                		available_trees_combo += available_trees_combo_item
                					.replace('{val}', index)
                					.replace('{sel}', ' selected')
                					.replace('{text}', index);
                	else
                		available_trees_combo += available_trees_combo_item
                					.replace('{val}', index)
                					.replace('{sel}', '')
                					.replace('{text}', index);
                }); 
                
                $('#trees_container').append(available_trees_combo);

                
                // COLOR PICKER
                $('.colorpicker').each(function(index, element) {
                    var color = '#000000';

                    $(element).css('background-color', color);
                    $(element).attr('color', color);
                });


                
            } else {
                // load state

                group_counter = state['group_counter'];
                group_count = state['group_count'];

                $('#treeControls').html(state['settings_html']);
                $('#groups').html(state['groups_html']);

                SELECTED = state['SELECTED'];

                categorical_data_ids = state['categorical_data_ids'];
                categorical_data_colors = state['categorical_data_colors'];

                stack_bar_colors = state['stack_bar_colors'];
                stack_bar_ids = state['stack_bar_ids'];

                metadata_swap_log = state['metadata_swap_log'];
                metadata_swap_log_reverse = state['metadata_swap_log_reverse'];

                var positions = {};
                for (var i=1; i < parameter_count; i++)
                    positions[i] = i;

                for (var i=0; i < state['metadata_swap_log_reverse'].length; i++)
                {
                    var new_positions = {};
                    for (var j=1; j < parameter_count; j++)
                    {
                        new_positions[j] = positions[state['metadata_swap_log_reverse'][i][j]];
                    }
                    positions = new_positions;
                }

                for (var i=0; i < metadata.length; i++)
                {
                    var new_line = new Array();
                    new_line.push(metadata[i][0]);

                    for (var pindex = 1; pindex < parameter_count; pindex++)
                    {
                        new_line.push(metadata[i][positions[pindex]]);    
                    }
                    metadata[i] = new_line.splice(0);
                }
            }

            // make table sortable
            var fixHelperModified = function(e, tr) {
                    var $originals = tr.children();
                    var $helper = tr.clone();
                    $helper.children().each(function(index) {
                        $(this).width($originals.eq(index).width());
                    });
                    return $helper;
                };

            $("#table_layers>tbody").sortable({helper: fixHelperModified}).disableSelection();

            $("#treeControls").dialog("option", "position", {
                my: "left center",
                at: "left center",
                of: window
            });
            $("#groups").dialog("option", "position", {
                my: "left top",
                at: "left bottom",
                of: "#treeControls"
            });

            $('.colorpicker').colpick({
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

            $('body').bind('click', function() {
                $('#control_contextmenu').hide();
            });

    }}); // meta
    }}); // state
});



//---------------------------------------------------------
//  ui callbacks
//---------------------------------------------------------

function saveCurrentState() {
    var state = {};

    state['group_counter'] = group_counter;
    state['group_count'] = group_count;

    // update dom 
    $(':text').each(function(){
         $(this).attr('value',$(this).val());
    });

    $('select').each(function(i, select) {
        var selected = $(select).val();
        $(select).find('option').removeAttr('selected').each(function(j, option) {
            if($(option).val()==selected)
            {
                $(option).attr('selected', true);
            }
        });
    });

    $(':checkbox').each(function() {
        if ($(this)[0].checked)
        {
            $(this).attr('checked', true);
        }
        else
        {
            $(this).removeAttr('checked');
        }
    });

    state['settings_html'] = $('#treeControls').html();
    state['groups_html'] = $('#groups').html();

    state['SELECTED'] = SELECTED;

    state['categorical_data_ids'] = categorical_data_ids;
    state['categorical_data_colors'] = categorical_data_colors;

    state['stack_bar_colors'] = stack_bar_colors;
    state['stack_bar_ids'] = stack_bar_ids;

    state['metadata_swap_log'] = metadata_swap_log;
    state['metadata_swap_log_reverse'] = metadata_swap_log_reverse;

    $.post("/save_state", {
        state: JSON.stringify(state),
    });
}

function menu_callback(action) {
    var contig_name = id_to_node_map[context_menu_target_id].label;

    switch (action) {

        case 'select':
            var fake_event = {'target': {'id': '#line' + context_menu_target_id}};
            lineClickHandler(fake_event);
            break;

        case 'remove':
            var fake_event = {'target': {'id': '#line' + context_menu_target_id}};
            lineContextMenuHandler(fake_event);
            break;

        case 'content':
            $.ajax({
                type: 'GET',
                cache: false,
                url: '/data/contig/' + contig_name + '?timestamp=' + new Date().getTime(),
                success: function(data) {
                    $("#contig_name_dialog").dialog("option", "title", contig_name);
                    $('#contig_names').val(data);
                    $('#contig_name_dialog').dialog('open');
                    $('#contig_names').click(); // focus & select all
                }
            });
            break;

        case 'metadata':
            $("#contig_name_dialog").dialog("option", "title", contig_name);
            $('#contig_names').val(strip(metadata_title[contig_name].join('\n')));
            $('#contig_name_dialog').dialog('open');
            $('#contig_names').click(); // focus & select all
            break;
        
        case 'inspect':
        	var layers = new Array();
        	metadata[0].forEach(function (e,i,a) {if(["contigs", "__parent__"].indexOf(e) < 0) layers.push(e);});
        	window.open('charts.html?contig=' + contig_name, '_blank');
        	break;
    }
}

function draw_tree_callback(){
	
    var trees_container_selection = $('#trees_container').val()
    if(selected_tree_id != trees_container_selection){
	    selected_tree_id = trees_container_selection;
	    
	    // empty the variable
        newick = '';
        window.newick = '';

	    // load data
	    $.ajax({
	        type: 'GET',
	        cache: false,
	        url: '/tree/' + selected_tree_id + '?timestamp=' + new Date().getTime(),
	        success: function(data) {
	            newick = data;
	        }
	    });
	
    }

    if (newick === '' || typeof metadata === 'undefined' || typeof contig_lengths === 'undefined') {
        setTimeout(draw_tree_callback, 200);
    } else {
        tree_type = $('#tree_type').val();
        draw_tree($('#tree_type').val());

        // enable export as svg button
        $('#btn_export_svg').attr('disabled', false);
    }
}

function showContigNames(gid) {
    var names = new Array();

    for (var j = 0; j < SELECTED[gid].length; j++) {
        if (id_to_node_map[SELECTED[gid][j]].IsLeaf()) {
            names.push(id_to_node_map[SELECTED[gid][j]].label);
        }
    }

    $("#contig_name_dialog").dialog("option", "title", "Contig Names");
    $('#contig_names').val(names.join("\n"));
    $('#contig_name_dialog').dialog('open');
    $('#contig_names').click(); // focus & select all
}

function newGroup() {
    group_counter++;
    group_count++;

    SELECTED[group_counter] = [];

    var clone = $('#group-template').clone();
    clone.appendTo('#groups-table');
    clone.attr('id', '');
    clone.css('display', '');

    var color = '#000000';

    clone.find('.colorpicker').css('background-color', color);
    clone.find('.colorpicker').attr('color', color);
    clone.find('.colorpicker').attr('id', 'group_color_' + group_counter);
    clone.find('.colorpicker').colpick({
        layout: 'hex',
        submit: 0,
        colorScheme: 'light',
        onChange: function(hsb, hex, rgb, el, bySetColor) {
            $(el).css('background-color', '#' + hex);
            $(el).attr('color', '#' + hex);

            if (!bySetColor) $(el).val(hex);
        },
        onHide: function() {
            redrawGroupColors();
        }
    }).keyup(function() {
        $(this).colpickSetColor(this.value);
    });

    clone.find('input[type=radio]').attr('value', group_counter).prop('checked', true);
    clone.find('input[type=text]').attr('value', "Group_" + group_counter).attr('id', 'group_name_' + group_counter);
    clone.find('.span_contig_length').attr('id', 'contig_length_' + group_counter);
    clone.find('input[type=button]').attr('id', 'contig_count_' + group_counter).click(function() {
        if (this.value == '0')
            return;
        showContigNames(group_counter);
    });
}

function deleteGroup(elm) {
    if (confirm('Are you sure you want to delete this group?')) {
        var id = $(elm).closest('tr').find('input[type=radio]').attr('value');
        $(elm).closest('tr').remove();
        $('input[type=radio]').last().prop('checked', true);
        group_count--;

        for (var i = 0; i < SELECTED[id].length; i++) {
            $("#line" + SELECTED[id][i]).css('stroke-width', '1');
            $("#arc" + SELECTED[id][i]).css('stroke-width', '1');
            $("#line" + SELECTED[id][i]).css('stroke', LINE_COLOR);
            $("#arc" + SELECTED[id][i]).css('stroke', LINE_COLOR);

            if (id_to_node_map[SELECTED[id][i]].IsLeaf())
            {
                $('.path_' + SELECTED[id][i] + "_background").css({'fill': '#FFFFFF', 'fill-opacity': '0.0'});
                $('.path_' + SELECTED[id][i] + "_outer_ring").css('fill', '#FFFFFF');
            }
        }

        SELECTED[id] = [];

        if (group_count==0)
        {
            newGroup();
        }
    }
}

function submitGroups(only_svg) {

    if (!only_svg) {
        var output = {};
        var msg_group_count = 0;
        var msg_contig_count = 0;

        for (var gid = 1; gid <= group_counter; gid++) {
            if (SELECTED[gid].length > 0) {
                msg_group_count++;
                var group_name = $('#group_name_' + gid).val();

                output[group_name] = new Array();
                for (var i = 0; i < SELECTED[gid].length; i++) {
                    if (id_to_node_map[SELECTED[gid][i]].IsLeaf()) {
                        output[group_name].push(id_to_node_map[SELECTED[gid][i]].label);
                        msg_contig_count++;
                    }
                }
            }
        }

        if (!confirm('You\'ve selected ' + msg_contig_count + ' contigs in ' + msg_group_count + ' group. You won\'t able to select more contigs after submit. Do you want to continue?')) {
            return;
        }
    }


    // draw group list to output svg
    drawGroupLegend();

    // move group highlights to new svg groups
    for (var gid = 1; gid <= group_counter; gid++) {

        createGroup('tree_group', 'selected_group_' + gid);

        for (var j = 0; j < SELECTED[gid].length; j++) {
            if (id_to_node_map[SELECTED[gid][j]].IsLeaf()) {
                $('.path_' + SELECTED[gid][j] + "_background").detach().appendTo('#selected_group_' + gid);
                $('.path_' + SELECTED[gid][j] + "_outer_ring").detach().appendTo('#selected_group_' + gid);
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

    if (!only_svg)
    {
        // disable group controls
        $('#group-template input[type=radio]').prop('checked', true);
        $('#groups-table input[type=radio]').attr("disabled", true);
        $('#submit-groups').attr("disabled", true);
        $('#btn_new_group').attr("disabled", true);

        $.post("/submit", {
            groups: JSON.stringify(output),
            svg: document.getElementById('svg').outerHTML
        });
    }
    else
    {
        $.post("/submit", {
            groups: "{}",
            svg: document.getElementById('svg').outerHTML
        });
        // add removed ungrouped backgrounds back
        $(detached_paths).appendTo('#tree_group');

        $('#group_legend').remove();
    }
}

function updateGroupWindow() { 
    // count contigs and update group labels
    for (var gid = 1; gid <= group_counter; gid++) {
        var contigs = 0;
        var length_sum = 0;

        for (var j = 0; j < SELECTED[gid].length; j++) {
            if (id_to_node_map[SELECTED[gid][j]].IsLeaf())
            {
                contigs++;
                length_sum += parseInt(contig_lengths[id_to_node_map[SELECTED[gid][j]].label]);
            }
        }

        $('#contig_count_' + gid).val(contigs);
        $('#contig_length_' + gid).html(readableNumber(length_sum));
    }
}

//---------------------------------------------------------
//  zoom and scale functions
//---------------------------------------------------------

function smooth_scale(x) {
    return (Math.pow(Math.abs(x), 2) * (3 - 2 * Math.abs(x)))
}

function getMatrix(viewport) {
    return viewport.getAttribute('transform').split('(')[1].split(')')[0].split(',').map(parseFloat);
}

function setMatrix(viewport, matrix) {
    viewport.setAttribute('transform', 'matrix(' + matrix.join(',') + ')');
}

function zoom(viewport, scale) {
    matrix = getMatrix(viewport);
    old_scale = matrix[0];
    if (scale * old_scale < 0.1) scale = 0.1 / old_scale;
    if (scale * old_scale > 100) scale = 100 / old_scale;

    for (var i = 0; i < matrix.length; i++) {
        matrix[i] *= scale;

    }

    bbox = viewport.getBBox();

    matrix[4] += (1 - scale) * (bbox.width - 50) / 2;
    matrix[5] += (1 - scale) * bbox.height / 2;

    setMatrix(viewport, matrix);
}

function pan(viewport, dx, dy) {
    matrix = getMatrix(viewport);

    matrix[4] += dx;
    matrix[5] += dy;

    setMatrix(viewport, matrix);
}

function zoom_in() {
    var viewport = document.getElementById('viewport');
    gradual_zoom(viewport, ZOOM_IN);
}

function zoom_out() {
    var viewport = document.getElementById('viewport');
    gradual_zoom(viewport, ZOOM_OUT);
}

function mouse_zoom_in(event) {
    pan_to_mouse(event, ZOOM_IN);
    zoom_in();
}

function mouse_zoom_out(event) {
    pan_to_mouse(event, ZOOM_OUT);
    zoom_out();
}

function pan_btn(dir) {
    if (dir == '11') {
        var svg = document.getElementById('svg');
        // Scale to fit window
        var bbox = svg.getBBox();

        // move drawing to centre of viewport
        var viewport = document.getElementById('viewport');
        baseMatrix = [1 * SCALE_MATRIX, 0, 0, 1 * SCALE_MATRIX, VIEWER_WIDTH / 2, VIEWER_HEIGHT / 2];
        setMatrix(viewport, baseMatrix);


        // centre
        bbox = svg.getBBox();
        if (bbox.x < 0) {
            pan(viewport, bbox.width / 2, bbox.height / 2);
        } else {
            pan(viewport, -25, 25);
        }
        return;
    }
    start_pan();
    var viewport = document.getElementById('viewport');
    matrix = getMatrix(viewport);
    scale = matrix[0];
    d = 100;
    dx = dir == 'l' ? d : (dir == 'r' ? -d : 0);
    dy = dir == 'u' ? d : (dir == 'd' ? -d : 0);
    gradual_pan(viewport, dx, dy);
}

function pan_to_mouse(e, scale) {
    var container = document.getElementById('treeContainer');
    var viewport = document.getElementById('viewport');
    offset = $(container).offset();

    dx = e.pageX - (offset.left + container.offsetWidth / 2);
    dy = e.pageY - (offset.top + container.offsetHeight / 2);

    gradual_pan(viewport, dx * (1 - scale), dy * (1 - scale));
}

var panning = 0;
var pan_pressed = false;

function start_pan() {
    pan_pressed = true;
}

function stop_pan() {
    pan_pressed = false;
}

function gradual_pan(viewport, dx, dy) {
    if (panning) return;
    steps = 25;
    duration = 0.5;
    panning = 0;
    last_time = new Date() / 1000;

    function frame() {
        this_time = new Date() / 1000;
        time_elapsed = this_time - last_time;
        last_time = this_time;
        panning += time_elapsed / duration;
        if (panning >= 1) panning = 1;
        d = smooth_scale(panning);
        pan(viewport, d * dx / steps, d * dy / steps)
        if (panning >= 1 && !(pan_pressed)) {
            panning = 0;
            clearInterval(id);
        }
    }
    var id = setInterval(frame, 1000 * duration / steps);
}

var zooming = 0;
var zoom_pressed = false;

function start_zoom() {
    zoom_pressed = true;
}

function stop_zoom() {
    zoom_pressed = false;
}

function gradual_zoom(viewport, scale) {
    if (zooming > 0 || panning || pan_pressed) return;

    var viewport = document.getElementById('viewport');
    matrix = getMatrix(viewport);

    initial_scale = matrix[0];
    new_scale = initial_scale * scale;
    old_scale = initial_scale;

    steps = 25;
    duration = 0.5;
    last_time = new Date() / 1000;

    function frame() {
        console.log(zoom_pressed);
        this_time = new Date() / 1000;
        time_elapsed = this_time - last_time;
        last_time = this_time;
        mult = scale < 1 ? zooming : Math.pow(zooming, 2.5);
        this_scale = initial_scale + (new_scale - initial_scale) * mult;
        this_step = this_scale / old_scale;
        old_scale = this_scale;
        zooming += time_elapsed / duration;
        zoom(viewport, this_step)
        if (zooming >= 1 && !(zoom_pressed)) {
            clearInterval(id);
            zooming = 0;
        }
    }
    var id = setInterval(frame, 1000 * duration / steps);
}
