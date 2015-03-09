function getGroupId() {
    var radios = document.getElementsByName('active_group');
    for(var i=0; i < radios.length; i++)
    {
        if (radios[i].checked)
            return radios[i].value;
    }
}

function lineClickHandler(event) {
    var p = getNodeFromEvent(event);

    var group_id = getGroupId();

    if (group_id === 'undefined')
        return;

    if (p.child_nodes.length > 2500 && !confirm("You just made a very big selection. Please click OK if it was intentional."))
        return;

    var group_color = document.getElementById('group_color_' + group_id).getAttribute('color');

    var groups_to_update = [];
    for (var i = 0; i < p.child_nodes.length; i++) {
        var pos = SELECTED[group_id].indexOf(id_to_node_map[p.child_nodes[i]].label);
        if (pos == -1) {
            SELECTED[group_id].push(id_to_node_map[p.child_nodes[i]].label);

            if (groups_to_update.indexOf(group_id) == -1)
                groups_to_update.push(group_id);
        }

        // remove nodes from other groups
        for (var gid = 1; gid <= group_counter; gid++) {
            // don't remove nodes from current group
            if (gid == group_id)
                continue;

            var pos = SELECTED[gid].indexOf(id_to_node_map[p.child_nodes[i]].label);
            if (pos > -1) {
                SELECTED[gid].splice(pos, 1);

                if (groups_to_update.indexOf(gid) == -1)
                    groups_to_update.push(gid);
            }
        }
    }

    redrawGroups();
    updateGroupWindow(groups_to_update);
}

function lineContextMenuHandler(event) {
    if (event.preventDefault) event.preventDefault();
    var group_id = getGroupId();

    if (event.target.id.indexOf('path_') > -1) // if layer -> show popup
    {
        // hide tooltip
        var tooltip = document.getElementById('aToolTip');
        if (tooltip)
            tooltip.parentNode.removeChild(tooltip);

        context_menu_target_id = getNodeFromEvent(event).id;

        $('#control_contextmenu').show();

        if (group_id > 0)
        {
            var pos = SELECTED[group_id].indexOf(id_to_node_map[parseInt(context_menu_target_id)].label);

            if (pos == -1) {
                $('#control_contextmenu #select').show();
                $('#control_contextmenu #remove').hide();
            }
            else
            {
                $('#control_contextmenu #select').hide();
                $('#control_contextmenu #remove').show();
            }
        }
        else
        {
            $('#control_contextmenu #select').hide();
            $('#control_contextmenu #remove').hide();
        }

        $('#control_contextmenu').offset({left:event.pageX-2,top:event.pageY-2});
        return false;
    }

    var p = getNodeFromEvent(event);

    if (group_id === 'undefined')
        return;

    if (p.child_nodes.length > 2500 && !confirm("You are about to unselect a large number of contigs. Press OK to continue."))
        return;

    var groups_to_update = [];
    for (var i = 0; i < p.child_nodes.length; i++) {
        // remove nodes from all groups
        for (var gid = 1; gid <= group_counter; gid++) {
            var pos = SELECTED[gid].indexOf(id_to_node_map[p.child_nodes[i]].label);
            if (pos > -1) {
                SELECTED[gid].splice(pos, 1);

                if (groups_to_update.indexOf(gid) == -1)
                    groups_to_update.push(gid);
            }
        }
    }
    redrawGroups();
    updateGroupWindow(groups_to_update);
    lineMouseLeaveHandler(event);
    return false;
}

function lineMouseEnterHandler(event) {
    var p = getNodeFromEvent(event);

    $('#path_hover').remove();

    var group_id = getGroupId();

    if (group_id === 'undefined')
        return;

    var group_color = document.getElementById('group_color_' + group_id).getAttribute('color');

    var p1 = p;
    while (p1.child) {
        p1 = p1.child;
    }

    var p2 = p;

    while (p2.child) {
        p2 = p2.child.GetRightMostSibling();
    }

    if (tree_type == 'circlephylogram')
    {
        drawPie('tree_group',
            'hover',
            p1.angle - angle_per_leaf / 2,
            p2.angle + angle_per_leaf / 2,
            distance(p.backarc, {
                'x': 0,
                'y': 0
            }),
            total_radius,
            (p2.angle - p1.angle + angle_per_leaf > Math.PI) ? 1 : 0,
            group_color,
            0.3,
            false);
    }
    else
    {  
        drawPhylogramRectangle('tree_group',
            'hover',
            p.ancestor.xy.x,
            (p1.xy.y + p2.xy.y) / 2,
            p2.xy.y - p1.xy.y + height_per_leaf,
            total_radius - p.ancestor.xy.x,
            group_color,
            0.3,
            false);
   }

    for (var index = 0; index < p.child_nodes.length; index++) {
        var _line = document.getElementById('line' + p.child_nodes[index]);
        if (_line) {
            _line.style['stroke-width'] = '3';
            _line.style['stroke'] = group_color;       
        }

        var _arc = document.getElementById('arc' + p.child_nodes[index]);
        if (_arc) {
            _arc.style['stroke-width'] = '3';
            _arc.style['stroke'] = group_color;
        }
    }
}

function lineMouseLeaveHandler(event) {
    var p = getNodeFromEvent(event);

    $('#path_hover').remove();

    var group_id = getGroupId();

    if (group_id === 'undefined') {
        document.focus();
        return;
    }

    for (var index = 0; index < p.child_nodes.length; index++) {
        var _line = document.getElementById('line' + p.child_nodes[index]);
        if (_line) {
            _line.style['stroke-width'] = '1';       
        }

        var _arc = document.getElementById('arc' + p.child_nodes[index]);
        if (_arc) {
            _arc.style['stroke-width'] = '1';
        }
    }

    var node_stack = [];
    for (var gid = 1; gid <= group_counter; gid++) {
        var color_picker = document.getElementById('group_color_' + gid);

        if (!color_picker)
            continue;

        var group_color = color_picker.getAttribute('color');

        for (var i = 0; i < SELECTED[gid].length; i++) {
            node_stack.push(label_to_node_map[SELECTED[gid][i]].id);

            var _line = document.getElementById('line' + label_to_node_map[SELECTED[gid][i]].id);
            if (_line) {
                _line.style['stroke-width'] = '2';
                _line.style['stroke'] = group_color;       
            }

            var _arc = document.getElementById('arc' + label_to_node_map[SELECTED[gid][i]].id);
            if (_arc) {
                _arc.style['stroke-width'] = '2';
                _arc.style['stroke'] = group_color;
            }
        }
    }

    for (var i = 0; i < p.child_nodes.length; i++) {
        if (node_stack.indexOf(p.child_nodes[i]) > -1)
            continue;

        var _line = document.getElementById('line' + p.child_nodes[i]);
        if (_line) {
            _line.style['stroke'] = LINE_COLOR;       
        }

        var _arc = document.getElementById('arc' + p.child_nodes[i]);
        if (_arc) {
            _arc.style['stroke'] = LINE_COLOR;
        }
    }
}

function mouseMoveHandler(event) {
    var p = getNodeFromEvent(event);

    if (event.target.id && event.target.id == 'path_event')
        lineMouseEnterHandler(event);
    
    var tooltip = document.getElementById('aToolTip');
    if (tooltip)
        tooltip.parentNode.removeChild(tooltip);

    if (tree_type == 'circlephylogram')
    {    
        var _y = event.clientY - origin_y;
        var _x = event.clientX - origin_x;
        var distance = Math.sqrt(Math.pow(_x,2) + Math.pow(_y,2)) / zoom_factor;
    }
    else
    {
        var distance = (event.clientY - tree_top) / zoom_factor;
    }

    var layer_id;
    for (var i=0; i < layer_boundaries.length; i++)
    {
        if (distance > layer_boundaries[i][0] && distance < layer_boundaries[i][1])
        {
            layer_id = i;
            break;
        }
    }

    var tooltip_arr = metadata_title[id_to_node_map[p.id].label].slice(0);
    tooltip_arr[layer_id] = '<font color="lime">' + tooltip_arr[layer_id] + '</font>';
    var message = tooltip_arr.join('<br />\n');

    var tooltip = document.createElement('div');
    tooltip.setAttribute('id', 'aToolTip');
    tooltip.setAttribute('class', 'defaultTheme');
    tooltip.innerHTML = "<p class='aToolTipContent'>"+message+"</p>";

    tooltip.style['top'] = (event.y+10) + 'px';
    tooltip.style['left'] = (event.x+10) + 'px';

    document.body.appendChild(tooltip);
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
                    messagePopupShow(contig_name, data);
                }
            });
            break;

        case 'metadata':
            messagePopupShow(contig_name, strip(metadata_title[contig_name].join('\n')));
            break;
        
        case 'inspect':
            var layers = new Array();
            metadata[0].forEach(function (e,i,a) {if(["contigs", "__parent__"].indexOf(e) < 0) layers.push(e);});
            window.open('charts.html?contig=' + contig_name, '_blank');
            break;
    }
}

// globals related single background
var rect_left;
var rect_width;
var tree_top;

var origin_x;
var origin_y;

var zoom_factor;

function updateSingleBackgroundGlobals()
{
    zoom_factor = getMatrix()[0];

    if (tree_type == 'phylogram')
    {
        var path_event = document.getElementById('path_event');
        var rect = path_event.getBoundingClientRect();

        var tree = document.getElementById('tree');
        var tree_rect = tree.getBoundingClientRect();

        rect_left = rect.left;
        rect_width = rect.width;
        tree_top = tree_rect.top;
    }
    else // circlephylogram
    {
        var root = document.getElementById('line0');
        var rect = root.getBoundingClientRect();

        var angle = id_to_node_map[0].angle;

        var halfPI = Math.PI / 2;

        if (angle < halfPI)
        {
            origin_x = rect.left;
            origin_y = rect.top;
        }
        else if (angle < 2 * halfPI)
        {
            origin_x = rect.left + rect.width;
            origin_y = rect.top;
        }
        else if (angle < 3 * halfPI)
        {
            origin_x = rect.left + rect.width;
            origin_y = rect.top + rect.height;
        }
        else // 4 * halfPI
        {
            origin_x = rect.left;
            origin_y = rect.top + rect.height;
        }
    }
}

function getNodeFromEvent(event)
{
    if (event.target.id == 'path_event')
    {
        if (tree_type == 'phylogram')
        {
            return order_to_node_map[leaf_count - parseInt((event.clientX - rect_left) / (rect_width / leaf_count)) - 1];
        }
        else
        {
            var _y = event.clientY - origin_y;
            var _x = event.clientX - origin_x;

            var angle = Math.atan2(_y, _x) - angle_per_leaf / 2;
            if (angle < 0)
                angle = 2 * Math.PI + angle;

            var order = Math.ceil((angle - Math.toRadians(last_settings['angle-min'])) / angle_per_leaf);
            
            if (order < 1 || order > leaf_count)
                order = 0;

            return order_to_node_map[order]
        }
    }
    else
    {
        var id = event.target.id.match(/\d+/);

        if (id)
            return id_to_node_map[id[0]];
    }
}