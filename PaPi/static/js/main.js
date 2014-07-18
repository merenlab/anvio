//--------------------------------------------------------------------------------------------------
//  Globals
//--------------------------------------------------------------------------------------------------

var VIEWER_WIDTH = window.innerWidth || document.documentElement.clientWidth || document.getElementsByTagName('body')[0].clientWidth;
var VIEWER_HEIGHT = window.innerHeight || document.documentElement.clientHeight || document.getElementsByTagName('body')[0].clientHeight;

var ZOOM_IN = 1.33;
var ZOOM_OUT = 0.75;

var LINE_COLOR='#888888';

var SCALE_MATRIX = 0;
var id_to_node_map = new Array();
var angle_per_leaf;

var total_radius = 0;

var SELECTED = new Array();

var newick;
var metadata;
var parameter_count;

var group_counter = 0; // for id
var group_count = 0;

var categorical_data_ids = new Array();
var categorical_data_colors = new Array();

var stack_bar_ids = new Array();
var stack_bar_colors = new Array();

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

    // load data
    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/tree?timestamp=' + new Date().getTime(),
        success: function(data) {
            newick = data;
        }
    });

    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/meta?timestamp=' + new Date().getTime(),
        success: function(data) {
            metadata = eval(data);

            parameter_count = metadata[0].length;

            categorical_data_ids = [];
            separated_data_ids = [];

            for (var i = 1; i < parameter_count; i++) {
                if (metadata[1][i].indexOf(';') > -1) // stack data
                {
                    stack_bar_ids.push(i);
                    stack_bar_colors[i] = new Array();

                    for (var j=0; j < metadata[1][i].split(";").length; j++)
                    {
                        stack_bar_colors[i].push(randomColor());
                    }

                    var stack_bar_row_str = '<tr>' +
                        '<td><img class="drag-icon" src="images/drag.gif" /></td>' +
                        '<td title="' + metadata[0][i] + '">' + ((metadata[0][i].length > 10) ? metadata[0][i].slice(0,10) + "..." : metadata[0][i]) + '</td>' +
                        '<td>' +
                        '    <select id="normalization{id}">' +
                        '        <option value="none">none</option>' +
                        '        <option value="sqrt">Square root</option>' +
                        '        <option value="log">Logarithm</option>' +
                        '    </select>' +
                        '</td>' +
                        '<td><input type="text" size="2" id="height{id}" value="50"></td>' +
                        '</tr>';

                    stack_bar_row_str = stack_bar_row_str.replace(new RegExp('{id}', 'g'), i);

                    $('#tbody_layers').append(stack_bar_row_str);
                }
                else if (!isNumber(metadata[1][i])) // categorical data
                { 
                    categorical_data_ids.push(i);
                    categorical_data_colors[i] = new Array();

                    var categorical_data_row_str = '<tr>' +
                        '<td><img class="drag-icon" src="images/drag.gif" /></td>' +
                        '<td title="' + metadata[0][i] + '">' + ((metadata[0][i].length > 10) ? metadata[0][i].slice(0,10) + "..." : metadata[0][i]) + '</td>' +
                        '<td>n/a</td>' +
                        '<td>n/a</td>' +
                        '<td><input type="text" size="2" id="height{id}" value="30"></td>' +
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
                        '    <select id="normalization{id}">' +
                        '        <option value="none">none</option>' +
                        '        <option value="sqrt">Square root</option>' +
                        '        <option value="log">Logarithm</option>' +
                        '    </select>' +
                        '</td>' +
                        '<td><input type="text" size="2" id="height{id}" value="50"></td>' +
                        '</tr>';

                    numerical_data_row_str = numerical_data_row_str.replace(new RegExp('{id}', 'g'), i);

                    $('#tbody_layers').append(numerical_data_row_str);
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

            $('.colorpicker').each(function(index, element) {
                var color = randomColor();

                $(element).css('background-color', color);
                $(element).attr('color', color);
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
        }
    });
});



//---------------------------------------------------------
//  ui callbacks
//---------------------------------------------------------

function draw_tree_callback(){
    if (typeof newick === 'undefined' || typeof metadata === 'undefined') {
        setTimeout(draw_tree_callback, 200)
    } else {
        draw_tree($('#tree_type').val());

        // enable tooltips.
        $('path[title]').aToolTip();

        // remove title attributes after tooltips enabled
        $('path').each(function(index, elm) {
            $(elm).attr('title', '');
        });
    }
}

function showContigNames(gid) {
    var names = new Array();

    for (var j = 0; j < SELECTED[gid].length; j++) {
        if (id_to_node_map[SELECTED[gid][j]].IsLeaf()) {
            names.push(id_to_node_map[SELECTED[gid][j]].label);
        }
    }

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

    var color = randomColor();

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
            $("#line" + SELECTED[id][i]).css('stroke', LINE_COLOR);
            $("#arc" + SELECTED[id][i]).css('stroke', LINE_COLOR);
        }

        SELECTED[id] = [];
    }
}

function submitGroups() {
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

    if (confirm('You\'ve selected ' + msg_contig_count + ' contigs in ' + msg_group_count + ' group. You won\'t able to select more contigs after submit. Do you want to continue?')) {
        // draw group list to output svg
        drawGroupLegend();

        // move group highlights to new svg groups
        for (var gid = 1; gid <= group_counter; gid++) {

            createGroup('viewport', 'selected_group_' + gid);

            for (var j = 0; j < SELECTED[gid].length; j++) {
                if (id_to_node_map[SELECTED[gid][j]].IsLeaf()) {
                    $('.path_' + SELECTED[gid][j] + "_background").detach().appendTo('#selected_group_' + gid);
                    $('.path_' + SELECTED[gid][j] + "_outer_ring").detach().appendTo('#selected_group_' + gid);
                }
            }
        }

        // remove ungrouped backgrounds.
        $('#viewport > path').remove();

        // disable group controls
        $('#group-template input[type=radio]').prop('checked', true);
        $('#groups-table input[type=radio]').attr("disabled", true);
        $('#submit-groups').attr("disabled", true);

        $.post("/submit", {
            groups: JSON.stringify(output),
            svg: Base64.encode($('#svg')[0].outerHTML)
        });

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
        baseMatrix = [1 * SCALE_MATRIX, 0, 0, 1 * SCALE_MATRIX, 100, 100];
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