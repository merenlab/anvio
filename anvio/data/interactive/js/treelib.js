 /**
 *
 * Javascript library to display phylogenetic trees
 *
 */

//--------------------------------------------------------------------------------------------------
function createGroup(parent, group_id) {
    var svgObject = document.getElementById(parent);
    var g = document.createElementNS('http://www.w3.org/2000/svg', 'g');
    g.setAttribute('id', group_id);
    svgObject.appendChild(g);
}
//--------------------------------------------------------------------------------------------------

function drawLegend(top, left) {
    var legend_counter=0;

    $.each(layer_types, function (i, _) {

        if (layer_types[i] != 2)
            return; //not categorical, return equal to continue in native loop

        // collect categorical
        var pindex = i;
        var categorical_data_title = metadata[0][pindex];

        var names = new Array();

        for (var name in categorical_data_colors[pindex]) {
            names.push(name);
        }

        names.sort();

        var group_id = 'legend_' + legend_counter;
        legend_counter++;

        createGroup('viewport', group_id);

        // draw border
        drawRectangle(group_id, left - 10, top - 20, (names.length + 2.5) * 20, 200, 'white', 1, 'black');

        drawText(group_id, {
            'x': left,
            'y': top
        }, categorical_data_title, '16px');

        for (var j = 0; j < names.length; j++) {
            var name = names[j];

            top = top + 20;
            var rect = drawRectangle(group_id, left, top, 16, 16, categorical_data_colors[pindex][name], 1, 'black',
                null,
                function() {
                    // mouseenter
                    $(this).css('stroke-width', '2');
                },
                function() {
                    // mouseleave
                    $(this).css('stroke-width', '1');
                });

            rect.setAttribute('callback_pindex', pindex);
            rect.setAttribute('callback_name', name);

            $(rect).colpick({
                layout: 'hex',
                submit: 0,
                colorScheme: 'dark',
                onChange: function(hsb, hex, rgb, el, bySetColor) {
                    $(el).css('fill', '#' + hex);
                    categorical_data_colors[$(el).attr('callback_pindex')][$(el).attr('callback_name')] = '#' + hex;
                }
            });

            drawText(group_id, {
                'x': left + 30,
                'y': top + 8
            }, names[j], '12px');
        }
        top = top + 70;
    });

    //draw stack bar
    $.each(layer_types, function (i, _) {

        if (layer_types[i] != 1)
            return; //not stack bar, return equal to continue in native loop

        var pindex = i;
        var stack_bar_title = metadata[0][pindex];

        var names = stack_bar_title.split(";");

        var group_id = 'legend_' + legend_counter;
        legend_counter++;

        createGroup('viewport', group_id);
        drawRectangle(group_id, left - 10, top - 20, (names.length + 2.5) * 20, 200, 'white', 1, 'black');
        drawText(group_id, {
            'x': left,
            'y': top
        }, stack_bar_title, '16px');

        for (var j = 0; j < names.length; j++) {
            var name = names[j];

            top = top + 20;
            var rect = drawRectangle(group_id, left, top, 16, 16, stack_bar_colors[pindex][j], 1, 'black',
                null,
                function() {
                    // mouseenter
                    $(this).css('stroke-width', '2');
                },
                function() {
                    // mouseleave
                    $(this).css('stroke-width', '1');
                });

            rect.setAttribute('callback_pindex', pindex);
            rect.setAttribute('callback_id', j);

            $(rect).colpick({
                layout: 'hex',
                submit: 0,
                colorScheme: 'dark',
                onChange: function(hsb, hex, rgb, el, bySetColor) {
                    $(el).css('fill', '#' + hex);
                    stack_bar_colors[$(el).attr('callback_pindex')][parseInt($(el).attr('callback_id'))] = '#' + hex;
                }
            });

            drawText(group_id, {
                'x': left + 30,
                'y': top + 8
            }, names[j], '12px');
        }
        top = top + 70;

    });

}

function drawGroupLegend(groups_to_draw, top, left) {
    createGroup('viewport', 'group_legend');

    drawRectangle('group_legend', left - 10, top - 20,20 + (groups_to_draw.length + 2.5) * 20, 300, 'white', 1, 'black');
    drawText('group_legend', {
        'x': left,
        'y': top
    }, "Groups", '16px');

    // table titles
    top = top + 28;
    drawText('group_legend', {'x': left, 'y': top }, 'Color', '10px');
    drawText('group_legend', {'x': left + 30, 'y': top}, 'Name', '10px');
    drawText('group_legend', {'x': left + 170, 'y': top}, 'Contigs', '10px');
    drawText('group_legend', {'x': left + 230, 'y': top}, 'Length', '10px');


    for (var gid=0; gid < groups_to_draw.length; gid++) {
        var group = groups_to_draw[gid];
        top = top + 20;

        drawRectangle('group_legend', left, top-8, 16, 16, group['color'], 1, 'black');
        drawText('group_legend', {'x': left + 30, 'y': top }, group['name'], '12px');
        drawText('group_legend', {'x': left + 170, 'y': top}, group['contig-count'], '12px');
        drawText('group_legend', {'x': left + 230, 'y': top}, group['contig-length'], '12px');
    }
}

function drawLayerLegend(layers, layer_order, top, left) {
    createGroup('viewport', 'layer_legend');

    // legend border
    drawRectangle('layer_legend', left - 10, top - 20,20 + (layer_order.length + 2.5) * 20, 300, 'white', 1, 'black');
    
    // legend title
    drawText('layer_legend', {
        'x': left,
        'y': top
    }, "Layers", '16px');

    // table titles
    top = top + 28;
    drawText('layer_legend', {'x': left, 'y': top }, 'Color', '10px');
    drawText('layer_legend', {'x': left + 30, 'y': top}, 'Name', '10px');
    drawText('layer_legend', {'x': left + 120, 'y': top}, 'Norm.', '10px');
    drawText('layer_legend', {'x': left + 160, 'y': top}, 'Height', '10px');
    drawText('layer_legend', {'x': left + 200, 'y': top}, 'Min', '10px');
    drawText('layer_legend', {'x': left + 245, 'y': top}, 'Max', '10px');

    // table items
    for (var i = 0; i < layer_order.length; i++) 
    {
        var pindex = layer_order[i];
        var layer = layers[pindex];
        var layer_name = metadata[0][pindex];
        var short_name = (layer_name.length > 10) ? layer_name.slice(0,10) + "..." : layer_name;

        top = top + 20;

        // color
        if (layer.hasOwnProperty('color') && typeof layer['color'] != 'undefined') 
            drawRectangle('layer_legend', left, top - 8, 16, 16, layer['color'], 1, 'black');

        // name
        drawText('layer_legend', {'x': left + 30, 'y': top}, short_name, '12px');

        // normalization
        if (layer.hasOwnProperty('normalization') && typeof layer['normalization'] != 'undefined') {
            var _norm = layer['normalization'];   
        } else {
            var _norm = "-";
        }
        drawText('layer_legend', {'x': left + 120, 'y': top}, _norm, '12px');

        // height
        drawText('layer_legend', {'x': left + 160, 'y': top}, layer['height'], '12px');

        // min & max
        if (layer['min'].hasOwnProperty('value') && typeof layer['min']['value'] != 'undefined') {
            var _min = layer['min']['value'];
            var _max = layer['max']['value'];

            // normalize floating numbers 
            if (_min % 1 !== 0)
                _min = parseFloat(_min).toFixed(4);
            if (_max % 1 !== 0)
                _max = parseFloat(_max).toFixed(4);
        } else {
            var _min = "-";
            var _max = "-"; 
        }
        drawText('layer_legend', {'x': left + 200, 'y': top}, _min, '12px');
        drawText('layer_legend', {'x': left + 245, 'y': top}, _max, '12px');
    }

}

function drawLine(svg_id, p, p0, p1, isArc) {
    var line = document.createElementNS('http://www.w3.org/2000/svg', 'path');

    if (isArc) 
    {
        line.setAttribute('id', 'arc' + p.id);
    }
    else
    {
        line.setAttribute('id', 'line' + p.id);
    }

    line.setAttribute('vector-effect', 'non-scaling-stroke');
    line.setAttribute('style', 'stroke:' + LINE_COLOR + ';stroke-width:1;');
    line.setAttribute('d', linePath(p0, p1));

    var svg = document.getElementById(svg_id);
    svg.appendChild(line);
}

//--------------------------------------------------------------------------------------------------
function drawText(svg_id, p, string, font_size, align, color) {

    var text = document.createElementNS('http://www.w3.org/2000/svg', 'text');
    //newLine.setAttribute('id','node' + p.id);

    if (typeof color !== 'undefined')
        text.setAttribute('fill', color);

    text.setAttribute('style', 'alignment-baseline:middle;');
    text.setAttribute('x', p['x']);
    text.setAttribute('y', p['y']);
    text.setAttribute('font-size', font_size);
    text.setAttribute('pointer-events', 'none');

    switch (align) {
        case 'left':
            text.setAttribute('text-anchor', 'start');
            break;
        case 'centre':
        case 'center':
            text.setAttribute('text-anchor', 'middle');
            break;
        case 'right':
            text.setAttribute('text-anchor', 'end');
            break;
        default:
            text.setAttribute('text-anchor', 'start');
            break;
    }

    var textNode = document.createTextNode(string)
    text.appendChild(textNode);

    var svg = document.getElementById(svg_id);
    svg.appendChild(text);
}

//--------------------------------------------------------------------------------------------------
function drawRotatedText(svg_id, p, string, angle, align, font_size) {
    var text = document.createElementNS('http://www.w3.org/2000/svg', 'text');
    //newLine.setAttribute('id','node' + p.id);
    text.setAttribute('style', 'alignment-baseline:middle');
    text.setAttribute('x', p['x']);
    text.setAttribute('y', p['y']);
    text.setAttribute('pointer-events', 'none');
    text.setAttribute('font-size', font_size);

    switch (align) {
        case 'left':
            text.setAttribute('text-anchor', 'start');
            break;
        case 'centre':
        case 'center':
            text.setAttribute('text-anchor', 'middle');
            break;
        case 'right':
            text.setAttribute('text-anchor', 'end');
            break;
        default:
            text.setAttribute('text-anchor', 'start');
            break;
    }

    if (angle != 0) {
        text.setAttribute('transform', 'rotate(' + angle + ' ' + p['x'] + ' ' + p['y'] + ')');
    }

    var textNode = document.createTextNode(string)
    text.appendChild(textNode);

    var svg = document.getElementById(svg_id);
    svg.appendChild(text);
}

function drawStraightGuideLine(svg_id, id, xy, max_x)  {
    var line = document.createElementNS('http://www.w3.org/2000/svg', 'line');
    line.setAttribute('id','guide' + id);

    line.setAttribute('x1', xy['x']);
    line.setAttribute('y1', xy['y']);
    line.setAttribute('x2', max_x);
    line.setAttribute('y2', xy['y']);

    line.setAttribute('stroke', LINE_COLOR);
    line.setAttribute('stroke-opacity', '0.2');
    line.setAttribute('vector-effect', 'non-scaling-stroke');
    line.setAttribute('stroke-width', '1');

    var svg = document.getElementById(svg_id);
    svg.appendChild(line);
}

function drawGuideLine(svg_id, id, angle, start_radius, end_radius) {
    var line = document.createElementNS('http://www.w3.org/2000/svg', 'line');
    line.setAttribute('id','guide' + id);

    var ax = Math.cos(angle) * start_radius;
    var ay = Math.sin(angle) * start_radius;

    var bx = Math.cos(angle) * end_radius;
    var by = Math.sin(angle) * end_radius;

    line.setAttribute('x1', ax);
    line.setAttribute('y1', ay);
    line.setAttribute('x2', bx);
    line.setAttribute('y2', by); 

    line.setAttribute('stroke', LINE_COLOR);
    line.setAttribute('stroke-opacity', '0.2');
    line.setAttribute('vector-effect', 'non-scaling-stroke');
    line.setAttribute('stroke-width', '1');

    var svg = document.getElementById(svg_id);
    svg.appendChild(line);
}

function drawPie(svg_id, id, start_angle, end_angle, inner_radius, outer_radius, large_arc_flag, color, fill_opacity, pointer_events) {
    var pie = document.createElementNS('http://www.w3.org/2000/svg', 'path');

    if (start_angle > end_angle) {
        // swap 
        var t = end_angle;
        end_angle = start_angle;
        start_angle = t;
    }

    // origin
    var ox = 0;
    var oy = 0;

    // calculate points
    var ax = ox + Math.cos(start_angle) * inner_radius;
    var ay = ox + Math.sin(start_angle) * inner_radius;

    var bx = ox + Math.cos(end_angle) * inner_radius;
    var by = ox + Math.sin(end_angle) * inner_radius;

    var cx = ox + Math.cos(end_angle) * outer_radius;
    var cy = ox + Math.sin(end_angle) * outer_radius;

    var dx = ox + Math.cos(start_angle) * outer_radius;
    var dy = ox + Math.sin(start_angle) * outer_radius;

    // generate path string

    var path = new Array("M", ax, ay, // start point
        "A", inner_radius, inner_radius, 0, large_arc_flag, 1, bx, by, // inner arc
        "L", cx, cy, // line 1
        "A", outer_radius, outer_radius, 0, large_arc_flag, 0, dx, dy, // outer arc
        "Z"); // close path line 2

    pie.setAttribute('id', 'path_' + id);
    pie.setAttribute('class', 'path_' + id);
    pie.setAttribute('fill', color);
    pie.setAttribute('stroke-width', '0');
    pie.setAttribute('shape-rendering', 'auto');
    //pie.setAttribute('stroke', 'black');
    pie.setAttribute('d', path.join(" "));
    pie.setAttribute('fill-opacity', fill_opacity);

    if (!pointer_events)
        pie.setAttribute('pointer-events', 'none');

    var svg = document.getElementById(svg_id);
    svg.appendChild(pie);
}

function drawPhylogramRectangle(svg_id, id, x, y, height, width, color, fill_opacity, pointer_events) {
    var rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');

    rect.setAttribute('id', 'path_' + id);
    rect.setAttribute('class', 'path_' + id);
    rect.setAttribute('fill', color);
    rect.setAttribute('stroke-width', 0);
    rect.setAttribute('fill-opacity', fill_opacity);

    rect.setAttribute('x', x);
    rect.setAttribute('y', y - height / 2); 
    rect.setAttribute('width', width);
    rect.setAttribute('height', height);

    if (!pointer_events)
        rect.setAttribute('pointer-events', 'none');

    var svg = document.getElementById(svg_id);
    svg.appendChild(rect);
}

function drawRectangle(svg_id, x, y, height, width, fill, stroke_width, stroke_color, f_click, f_mouseenter, f_mouseleave) {
    var rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
    rect.setAttribute('fill', fill);
    rect.setAttribute('stroke-width', stroke_width);
    rect.setAttribute('stroke', stroke_color);

    rect.setAttribute('x', x);
    rect.setAttribute('y', y);
    rect.setAttribute('width', width);
    rect.setAttribute('height', height);

    $(rect).click(f_click);
    $(rect).mouseenter(f_mouseenter);
    $(rect).mouseleave(f_mouseleave);

    var svg = document.getElementById(svg_id);
    svg.appendChild(rect);

    return rect;
}

//--------------------------------------------------------------------------------------------------
function circeArcPath(p0, p1, radius, large_arc_flag) {
    var path = 'M ' + p0['x'] + ' ' + p0['y'] + ' A ' + radius + ' ' + radius + ' 0 ';

    if (large_arc_flag) {
        path += ' 1 ';
    } else {
        path += ' 0 ';
    }

    path += ' 1 ' + p1['x'] + ' ' + p1['y'];

    return path;
}

function drawCircleArc(svg_id, p, p0, p1, radius, large_arc_flag) {
    var arc = document.createElementNS('http://www.w3.org/2000/svg', 'path');
    arc.setAttribute('id', 'arc' + p.id);

    arc.setAttribute('vector-effect', 'non-scaling-stroke');
    arc.setAttribute('style', 'stroke:' + LINE_COLOR + ';stroke-width:1;');
    arc.setAttribute('fill', 'none');

    var path = circeArcPath(p0, p1, radius, large_arc_flag);
    arc.setAttribute('d', path)

    var svg = document.getElementById(svg_id);
    svg.appendChild(arc);
}

//--------------------------------------------------------------------------------------------------
function drawPath(svg_id, pathString) {
    var path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
    //newLine.setAttribute('id','node' + p.id);
    path.setAttribute('vector-effect', 'non-scaling-stroke');
    path.setAttribute('style', 'stroke:blue;stroke-width:1;');
    path.setAttribute('d', pathString);
    var svg = document.getElementById(svg_id);
    svg.appendChild(path);
}

//--------------------------------------------------------------------------------------------------
// Remove NEXUS-style string formatting, e.g. underscores
function formatString(s) {
    s = s.replace(/_/g, ' ');
    return s;
}

//--------------------------------------------------------------------------------------------------
// http://stackoverflow.com/questions/894860/set-a-default-parameter-value-for-a-javascript-function
function Node(label) {
    this.ancestor = null;
    this.child = null;
    this.sibling = null;
    this.label = typeof label !== 'undefined' ? label : '';
    this.id = 0;
    this.weight = 0;
    this.xy = [];
    this.edge_length = 0.0;
    this.path_length = 0.0;
    this.depth = 0;
    this.order = null;
}

//--------------------------------------------------------------------------------------------------
Node.prototype.IsLeaf = function() {
    return (!this.child);
}

//--------------------------------------------------------------------------------------------------
Node.prototype.GetRightMostSibling = function() {
    var p = this;
    while (p.sibling) {
        p = p.sibling;
    }
    return p;
}

//--------------------------------------------------------------------------------------------------
function Tree() {
    this.root = null;
    this.num_leaves = 0;
    this.num_nodes = 0;
    this.nodes = [];
    this.rooted = true;
    this.has_edge_lengths = false;
    this.error = 0;
}

//--------------------------------------------------------------------------------------------------
Tree.prototype.NewNode = function(label) {
    var node = new Node(label);
    node.id = this.num_nodes++;
    this.nodes[node.id] = node;
    id_to_node_map[node.id] = node;
    return node;
}

//--------------------------------------------------------------------------------------------------
Tree.prototype.Parse = function(str) {
    str = str.replace('"', "");

    // Strip NEXUS-style comments
    str = str.replace(/\[[^\[]+\]/g, "");

    str = str.replace(/\(/g, "|(|");
    str = str.replace(/\)/g, "|)|");
    str = str.replace(/,/g, "|,|");
    str = str.replace(/:/g, "|:|");
    str = str.replace(/;/g, "|;|");
    str = str.replace(/\|\|/g, "|");
    str = str.replace(/^\|/, "");
    str = str.replace(/\|$/, "");

    //console.log(str);

    var token = str.split("|");
    var curnode = this.NewNode();
    this.root = curnode;

    var state = 0;
    var stack = [];
    var i = 0;
    var q = null;

    var edge_length_norm = $('#edge_length_normalization')[0].checked;

    this.error = 0;

    while ((state != 99) && (this.error == 0)) {
        switch (state) {
            case 0:
                if (ctype_alnum(token[i].charAt(0))) {
                    this.num_leaves++;
                    label = token[i];
                    curnode.label = label;
                    i++;
                    state = 1;
                } else {
                    if (token[i].charAt(0) == "'") {
                        label = token[i];
                        label = label.replace(/^'/, "");
                        label = label.replace(/'$/, "");
                        this.num_leaves++;
                        curnode.label = label;
                        i++;
                        state = 1;
                    } else {
                        switch (token[i]) {
                            case '(':
                                state = 2;
                                break;

                            default:
                                state = 99;
                                this.error = 1; // syntax
                                break;
                        }

                    }
                }
                break;


            case 1: // getinternode
                switch (token[i]) {
                    case ':':
                    case ',':
                    case ')':
                        state = 2;
                        break;
                    default:
                        state = 99;
                        this.error = 1; // syntax
                        break;
                }
                break;

            case 2: // nextmove
                switch (token[i]) {
                    case ':':
                        i++;
                        if (isNumber(token[i])) {
                            // nnormalization of edge lengths
                            if (edge_length_norm) {
                                curnode.edge_length = Math.sqrt(parseFloat(token[i]) * 1000000) / 1000000;
                            } else {
                                curnode.edge_length = parseFloat(token[i]);
                            }
                            this.has_edge_lengths = true;
                            i++;
                        }
                        break;

                    case ',':
                        q = this.NewNode();
                        curnode.sibling = q;
                        var c = stack.length;
                        if (c == 0) {
                            state = 99;
                            this.error = 2; // missing (
                        } else {
                            q.ancestor = stack[c - 1];
                            curnode = q;
                            state = 0;
                            i++;
                        }
                        break;

                    case '(':
                        stack.push(curnode);
                        q = this.NewNode();
                        curnode.child = q;
                        q.ancestor = curnode;
                        curnode = q;
                        state = 0;
                        i++;
                        break;

                    case ')':
                        if (stack.length == 0) {
                            state = 99;
                            this.error = 3; // unbalanced
                        } else {
                            curnode = stack.pop();
                            state = 3;
                            i++;
                        }
                        break;

                    case ';':
                        if (stack.length == 0) {
                            state = 99;
                        } else {
                            state = 99;
                            this.error = 4; // stack not empty
                        }
                        break;

                    default:
                        state = 99;
                        this.error = 1; // syntax
                        break;
                }
                break;

            case 3: // finishchildren
                if (ctype_alnum(token[i].charAt(0))) {
                    curnode.label = token[i];
                    i++;
                } else {
                    switch (token[i]) {
                        case ':':
                            i++;
                            if (isNumber(token[i])) {
                                // nnormalization of edge lengths
                                if (edge_length_norm) {
                                    curnode.edge_length = Math.sqrt(parseFloat(token[i]) * 1000000) / 1000000;
                                } else {
                                    curnode.edge_length = parseFloat(token[i]);
                                }
                                this.has_edge_lengths = true;
                                i++;
                            }
                            break;

                        case ')':
                            if (stack.length == 0) {
                                state = 99;
                                this.error = 3; // unbalanced
                            } else {
                                curnode = stack.pop();
                                i++;
                            }
                            break;

                        case ',':
                            q = this.NewNode();
                            curnode.sibling = q;

                            if (stack.length == 0) {
                                state = 99;
                                this.error = 2; // missing (
                            } else {
                                q.ancestor = stack[stack.length - 1];
                                curnode = q;
                                state = 0;
                                i++;
                            }
                            break;

                        case ';':
                            state = 2;
                            break;

                        default:
                            state = 99;
                            this.error = 1; // syntax
                            break;
                    }
                }
                break;
        }
    }
}

//--------------------------------------------------------------------------------------------------
Tree.prototype.ComputeWeights = function(p) {
    if (p) {
        p.weight = 0;

        this.ComputeWeights(p.child);
        this.ComputeWeights(p.sibling);

        if (p.IsLeaf()) {
            p.weight = 1;
        }
        if (p.ancestor) {
            p.ancestor.weight += p.weight;
        }
    }
}

//--------------------------------------------------------------------------------------------------
Tree.prototype.ComputeDepths = function() {
    for (var i in this.nodes) {
        if (this.nodes[i].IsLeaf()) {
            p = this.nodes[i].ancestor;
            var count = 1;
            while (p) {
                p.depth = Math.max(p.depth, count);
                count++;
                p = p.ancestor;
            }
        }
    }
}

//--------------------------------------------------------------------------------------------------

function TreeDrawer() {
    //this.t = tree;

    this.leaf_count = 0;
    this.leaf_gap = 0;
    this.node_gap = 0;
    this.last_y = 0;

    this.svg_id;

    this.draw_scale_bar = false;
}

//--------------------------------------------------------------------------------------------------
TreeDrawer.prototype.Init = function(tree, settings) {
    this.t = tree;

    // defaults
    this.settings = settings;

    this.left = 0;
    this.top = 0;
    /*
    if (this.settings.fontHeight)
    {
        this.top += this.settings.fontHeight/2.0;
        this.settings.height -= this.settings.fontHeight;
    }
    */

}

//--------------------------------------------------------------------------------------------------
TreeDrawer.prototype.CalcInternal = function(p) {
    var pt = [];
    pt['x'] = this.left + this.node_gap * (this.t.num_leaves - p.weight);
    pt['y'] = this.last_y - ((p.weight - 1) * this.leaf_gap) / 2;
    p.xy = pt;
}

//--------------------------------------------------------------------------------------------------
TreeDrawer.prototype.CalcLeaf = function(p) {
    var pt = [];

    pt['y'] = this.top + (this.leaf_count * this.leaf_gap);
    this.last_y = pt['y'];
    this.leaf_count++;

    // slanted cladogram
    pt['x'] = this.left + this.settings.width;
    p.xy = pt;
}

//--------------------------------------------------------------------------------------------------
TreeDrawer.prototype.CalcNodeGap = function() {
    if (this.t.rooted) {
        this.node_gap = this.settings.width / this.t.num_leaves;
        this.left += this.node_gap;
        this.settings.width -= this.node_gap;
    } else {
        this.node_gap = this.settings.width / (this.t.num_leaves - 1);
    }
}


//--------------------------------------------------------------------------------------------------
TreeDrawer.prototype.CalcCoordinates = function() {
    this.t.ComputeWeights(this.t.root);

    this.leaf_count = 0;
    this.leaf_gap = this.settings.height / (this.t.num_leaves - 1);

    this.CalcNodeGap();

    var n = new NodeIterator(this.t.root);
    var q = n.Begin();
    while (q != null) {
        if (q.IsLeaf()) {
            this.CalcLeaf(q);
        } else {
            this.CalcInternal(q);
        }
        q = n.Next();
    }
}

//--------------------------------------------------------------------------------------------------
TreeDrawer.prototype.DrawLeaf = function(p) {
    var p0 = p.xy
    var anc = p.ancestor;
    if (anc) {
        var p1 = anc.xy;

        drawLine(this.settings.svg_id, p, p0, p1);
    }
}

//--------------------------------------------------------------------------------------------------
TreeDrawer.prototype.DrawInternal = function(p) {
    var p0 = p.xy
    var anc = p.ancestor;
    if (anc) {
        var p1 = anc.xy;
        drawLine(this.settings.svg_id, p, p0, p1);
    }
}

//--------------------------------------------------------------------------------------------------
TreeDrawer.prototype.DrawRoot = function() {
    var p0 = this.t.root.xy
    var p1 = [];
    p1['x'] = p0['x'];
    p1['y'] = p0['y'];
    p1['x'] -= this.node_gap;

    drawLine(this.settings.svg_id, this.t.root, p0, p1);
}

//--------------------------------------------------------------------------------------------------
TreeDrawer.prototype.Draw = function() {
    var n = new NodeIterator(this.t.root);
    var q = n.Begin();

    while (q != null) {
        if (q.IsLeaf()) {
            this.DrawLeaf(q);
        } else {
            this.DrawInternal(q);
        }
        q = n.Next();
    }
    if (this.t.rooted) {
        this.DrawRoot();
    }
}

//--------------------------------------------------------------------------------------------------
TreeDrawer.prototype.DrawLabels = function(nexus) {
    var nxs = typeof nexus !== 'undefined' ? nexus : null;

    var n = new NodeIterator(this.t.root);
    var q = n.Begin();
    while (q != null) {
        if (q.IsLeaf()) {
            var label = q.label;

            if (nxs) {
                if (nxs.treesblock.translate) {
                    if (nxs.treesblock.translate[label]) {
                        label = nxs.treesblock.translate[label];
                    }
                }
            }
            // offset 
            label_xy = q.xy;
            label_xy['x'] += this.settings.fontHeight / 2.0;

            drawText('viewport', label_xy, formatString(label));
        }
        q = n.Next();
    }
}

//--------------------------------------------------------------------------------------------------
RectangleTreeDrawer.prototype = new TreeDrawer();

function RectangleTreeDrawer() {
    TreeDrawer.apply(this, arguments);

    this.max_depth = 0;
};

//--------------------------------------------------------------------------------------------------
RectangleTreeDrawer.prototype.CalcInternal = function(p) {
    var pt = [];
    pt['x'] = this.left + this.node_gap * (this.t.root.depth - p.depth);

    var pl = p.child.xy;
    var pr = p.child.GetRightMostSibling().xy;

    pt['y'] = pl['y'] + (pr['y'] - pl['y']) / 2;
    p.xy['x'] = pt['x'];
    p.xy['y'] = pt['y'];
}

//--------------------------------------------------------------------------------------------------
RectangleTreeDrawer.prototype.CalcNodeGap = function() {
    this.t.ComputeDepths();
    //console.log(this.t.root.depth);
    if (this.t.rooted) {
        this.node_gap = this.settings.width / (this.t.root.depth + 1);
        this.left += this.node_gap;
        this.settings.width -= this.node_gap;
    } else {
        this.node_gap = this.settings.width / this.t.root.depth;
    }
}

//--------------------------------------------------------------------------------------------------
RectangleTreeDrawer.prototype.DrawLeaf = function(p) {
    var p0 = p.xy
    var p1 = [];
    var anc = p.ancestor;
    if (anc) {
        p1['x'] = anc.xy['x'];
        p1['y'] = p0['y'];

        drawLine(this.settings.svg_id, p, p0, p1);
    }
}

//--------------------------------------------------------------------------------------------------
RectangleTreeDrawer.prototype.DrawInternal = function(p) {
    var p0 = [];
    var p1 = [];

    p0['x'] = p.xy['x'];
    p0['y'] = p.xy['y'];

    var anc = p.ancestor;
    if (anc) {
        p1['x'] = anc.xy['x'];
        p1['y'] = p0['y'];

        drawLine(this.settings.svg_id, p, p0, p1);
    }

    // vertical line
    var pl = p.child.xy;
    var pr = p.child.GetRightMostSibling().xy;

    p0['x'] = p0['x'];
    p0['y'] = pl['y'];
    p1['x'] = p0['x'];
    p1['y'] = pr['y'];

    drawLine(this.settings.svg_id, p, p0, p1,true);
}


//--------------------------------------------------------------------------------------------------
PhylogramTreeDrawer.prototype = new RectangleTreeDrawer();

function PhylogramTreeDrawer() {
    RectangleTreeDrawer.apply(this, arguments);

    this.max_path_length = 0;
    this.draw_scale_bar = true;
};


//--------------------------------------------------------------------------------------------------
PhylogramTreeDrawer.prototype.CalcInternal = function(p) {
    var pt = [];
    pt['x'] = this.left + (p.path_length / this.max_path_length) * this.settings.width;

    var pl = p.child.xy;
    var pr = p.child.GetRightMostSibling().xy;

    pt['y'] = pl['y'] + (pr['y'] - pl['y']) / 2;
    p.xy['x'] = pt['x'];
    p.xy['y'] = pt['y'];
}

//--------------------------------------------------------------------------------------------------
PhylogramTreeDrawer.prototype.CalcLeaf = function(p) {
    var pt = [];
    pt['x'] = this.left + (p.path_length / this.max_path_length) * this.settings.width;

    pt['y'] = this.top + (this.leaf_count * this.leaf_gap);
    this.last_y = pt['y'];
    this.leaf_count++;

    p.xy['x'] = pt['x'];
    p.xy['y'] = pt['y'];

}


//--------------------------------------------------------------------------------------------------
PhylogramTreeDrawer.prototype.CalcCoordinates = function() {
    this.max_path_length = 0;
    //console.log(this.max_path_length);    

    this.t.root.path_length = this.t.root.edge_length;

    // build path lengths
    var n = new PreorderIterator(this.t.root);
    var q = n.Begin();
    while (q != null) {
        var d = q.edge_length;
        if (d < 0.00001) {
            d = 0.0;
        }
        if (q != this.t.root) {
            q.path_length = q.ancestor.path_length + d;
        }

        //console.log(q.label + ' ' + q.path_length + ' ' + q.edge_length);

        this.max_path_length = Math.max(this.max_path_length, q.path_length);
        q = n.Next();
    }

    //console.log(this.max_path_length);    

    this.leaf_count = 0;
    this.leaf_gap = this.settings.height / (this.t.num_leaves - 1);

    n = new NodeIterator(this.t.root);
    var q = n.Begin();
    while (q != null) {
        if (q.IsLeaf()) {
            this.CalcLeaf(q);
        } else {
            this.CalcInternal(q);
        }
        q = n.Next();
    }
}

//--------------------------------------------------------------------------------------------------
PhylogramTreeDrawer.prototype.Draw = function() {
    // parent method
    RectangleTreeDrawer.prototype.Draw.call(this);

    // scale bar
    if (this.draw_scale_bar) {
        this.DrawScaleBar();
    }
}

//--------------------------------------------------------------------------------------------------
PhylogramTreeDrawer.prototype.DrawScaleBar = function() {
    var p0 = [];
    var p1 = [];

    var m = log10(this.max_path_length);
    var i = Math.floor(m);
    var bar = Math.pow(10, i);

    var scalebar = (bar / this.max_path_length) * this.settings.width;

    p0['x'] = this.left;
    p0['y'] = this.top + this.settings.height + this.leaf_gap;

    p1['x'] = p0['x'] + scalebar;
    p1['y'] = p0['y'];

    //drawLine(this.settings.svg_id, 0, p0, p1);    
}



//--------------------------------------------------------------------------------------------------
CircleTreeDrawer.prototype = new RectangleTreeDrawer();

function CircleTreeDrawer() {
    RectangleTreeDrawer.apply(this, arguments);

    this.leaf_angle = 0;
    this.leaf_radius = 0;

    this.max_path_length = 0;
    this.root_length = 0;
};

//--------------------------------------------------------------------------------------------------
CircleTreeDrawer.prototype.CalcInternalRadius = function(p) {
    p.radius = this.node_gap * (this.t.root.depth - p.depth);
}

//--------------------------------------------------------------------------------------------------
CircleTreeDrawer.prototype.CalcInternal = function(p) {
    var left_angle = p.child.angle;
    var right_angle = p.child.GetRightMostSibling().angle;

    p.angle = left_angle + (right_angle - left_angle) / 2;

    this.CalcInternalRadius(p);

    var pt = [];
    pt['x'] = p.radius * Math.cos(p.angle);
    pt['y'] = p.radius * Math.sin(p.angle);

    p.xy['x'] = pt['x'];
    p.xy['y'] = pt['y'];

    var q = p.child;
    while (q) {
        pt = [];

        pt['x'] = p.radius * Math.cos(q.angle);
        pt['y'] = p.radius * Math.sin(q.angle);

        q.backarc = [];
        q.backarc['x'] = pt['x'];
        q.backarc['y'] = pt['y'];

        q = q.sibling;
    }
}


//--------------------------------------------------------------------------------------------------
CircleTreeDrawer.prototype.CalcLeafRadius = function(p) {
    p.radius = this.leaf_radius;
}

//--------------------------------------------------------------------------------------------------
CircleTreeDrawer.prototype.CalcLeaf = function(p) {
    var angle_min = Math.toRadians(parseFloat($('#angle-min').val()));
    p.angle = angle_min + this.leaf_angle * this.leaf_count;
    this.leaf_count++;

    this.CalcLeafRadius(p);

    var pt = [];
    pt['x'] = p.radius * Math.cos(p.angle);
    pt['y'] = p.radius * Math.sin(p.angle);

    p.xy['x'] = pt['x'];
    p.xy['y'] = pt['y'];
}

//--------------------------------------------------------------------------------------------------
CircleTreeDrawer.prototype.DrawLeaf = function(p) {

    var p0 = p.xy
    var p1 = p.backarc;

    drawLine(this.settings.svg_id, p, p0, p1);
}

//--------------------------------------------------------------------------------------------------
CircleTreeDrawer.prototype.DrawInternal = function(p) {
    var p0 = [];
    var p1 = [];

    p0['x'] = p.xy['x'];
    p0['y'] = p.xy['y'];

    var anc = p.ancestor;
    if (anc) {
        p0 = p.xy;
        p1 = p.backarc;

        drawLine(this.settings.svg_id, p, p0, p1);
    }

    // draw arc

    p0 = p.child.backarc;
    p1 = p.child.GetRightMostSibling().backarc;


    var large_arc_flag = (Math.abs(p.child.GetRightMostSibling().angle - p.child.angle) > Math.PI) ? true : false;
    drawCircleArc(this.settings.svg_id, p, p0, p1, p.radius, large_arc_flag);
}

//--------------------------------------------------------------------------------------------------
CircleTreeDrawer.prototype.DrawRoot = function() {
    var p0 = this.t.root.xy
    var p1 = [];
    p1['x'] = 0;
    p1['y'] = 0;


    drawLine(this.settings.svg_id, this.t.root, p0, p1);
}

//--------------------------------------------------------------------------------------------------
CircleTreeDrawer.prototype.CalcCoordinates = function() {
    this.t.ComputeDepths();

    this.max_path_length = 0;
    //console.log(this.max_path_length);    

    this.t.root.path_length = this.t.root.edge_length;

    // build path lengths
    var n = new PreorderIterator(this.t.root);
    var q = n.Begin();
    while (q != null) {
        var d = q.edge_length;
        if (d < 0.00001) {
            d = 0.0;
        }
        if (q != this.t.root) {
            q.path_length = q.ancestor.path_length + d;
        }

        //console.log(q.label + ' ' + q.path_length + ' ' + q.edge_length);

        this.max_path_length = Math.max(this.max_path_length, q.path_length);
        q = n.Next();
    }


    this.leaf_count = 0;

    var angle_max = parseFloat($('#angle-max').val());
    var angle_min = parseFloat($('#angle-min').val());

    this.leaf_angle = Math.toRadians(Math.abs(angle_max - angle_min)) / this.t.num_leaves;
    //this.leaf_angle = 2 * Math.PI / this.t.num_leaves;

    this.leaf_radius = this.settings.width / 2;
    this.node_gap = this.leaf_radius / this.t.root.depth;


    n = new NodeIterator(this.t.root);
    var q = n.Begin();
    while (q != null) {
        if (q.IsLeaf()) {
            this.CalcLeaf(q);
        } else {
            this.CalcInternal(q);
        }
        q = n.Next();
    }
}

//--------------------------------------------------------------------------------------------------
CircleTreeDrawer.prototype.Draw = function() {
    // parent method
    TreeDrawer.prototype.Draw.call(this);

    // move drawing to centre of viewport
    var viewport = document.getElementById(this.settings.svg_id);
    viewport.setAttribute('transform', 'translate(' + (this.settings.width + this.root_length) / 2 + ' ' + this.settings.height / 2 + ')');
}

//--------------------------------------------------------------------------------------------------
CircleTreeDrawer.prototype.DrawLabels = function(nexus) {
    var nxs = typeof nexus !== 'undefined' ? nexus : null;

    var n = new NodeIterator(this.t.root);
    var q = n.Begin();
    while (q != null) {
        if (q.IsLeaf()) {
            var label = q.label;

            if (nxs) {
                if (nxs.treesblock.translate) {
                    if (nxs.treesblock.translate[label]) {
                        label = nxs.treesblock.translate[label];
                    }
                }
            }

            var align = 'left';
            var angle = q.angle * 180.0 / Math.PI;
            if ((q.angle > Math.PI / 2.0) && (q.angle < 1.5 * Math.PI)) {
                align = 'right';
                angle += 180.0;
            }

            // offset label 
            var r = q.radius + this.settings.fontHeight / 2.0;
            var label_xy = [];
            label_xy['x'] = Math.cos(q.angle) * r;
            label_xy['y'] = Math.sin(q.angle) * r;

            drawRotatedText('viewport', label_xy, formatString(label), angle, align);
        }
        q = n.Next();
    }
}


//--------------------------------------------------------------------------------------------------
CirclePhylogramDrawer.prototype = new CircleTreeDrawer();

function CirclePhylogramDrawer() {
    CircleTreeDrawer.apply(this, arguments)

    this.max_path_length = 0;
    this.draw_scale_bar = true;
};


//--------------------------------------------------------------------------------------------------
CirclePhylogramDrawer.prototype.CalcInternalRadius = function(p) {
    p.radius = this.root_length + (p.path_length / this.max_path_length) * (this.settings.width / 2)
}

//--------------------------------------------------------------------------------------------------
CirclePhylogramDrawer.prototype.CalcLeafRadius = function(p) {
    p.radius = this.root_length + (p.path_length / this.max_path_length) * (this.settings.width / 2)
}

//--------------------------------------------------------------------------------------------------
CirclePhylogramDrawer.prototype.CalcCoordinates = function() {
    this.max_path_length = 0;
    //console.log(this.max_path_length);    

    if (this.settings.root_length) {
        this.root_length = this.settings.root_length * (this.settings.width / 2);
        this.settings.width -= 2 * this.root_length;
    }

    this.t.root.path_length = this.t.root.edge_length;

    // build path lengths
    var n = new PreorderIterator(this.t.root);
    var q = n.Begin();
    while (q != null) {
        var d = q.edge_length;
        if (d < 0.00001) {
            d = 0.0;
        }
        if (q != this.t.root) {
            q.path_length = q.ancestor.path_length + d;
        }

        this.max_path_length = Math.max(this.max_path_length, q.path_length);
        q = n.Next();
    }

    this.leaf_count = 0;
    var angle_max = parseFloat($('#angle-max').val());
    var angle_min = parseFloat($('#angle-min').val());

    this.leaf_angle = Math.toRadians(Math.abs(angle_max - angle_min)) / this.t.num_leaves;
    //this.leaf_angle = 2 * Math.PI / this.t.num_leaves;

    n = new NodeIterator(this.t.root);
    var q = n.Begin();
    while (q != null) {
        if (q.IsLeaf()) {
            this.CalcLeaf(q);
        } else {
            this.CalcInternal(q);
        }
        q = n.Next();
    }
}

//--------------------------------------------------------------------------------------------------
CirclePhylogramDrawer.prototype.Draw = function() {
    // parent method
    TreeDrawer.prototype.Draw.call(this);

    // move drawing to centre of viewport
    //var viewport = document.getElementById(this.settings.svg_id);
    //viewport.setAttribute('transform', 'translate(' + (this.settings.width + this.root_length) / 2 + ' ' + this.settings.height / 2 + ')');
}

function draw_tree(settings) {
    var tree_draw_timer = new BasicTimer('tree_draw');

    id_to_node_map = new Array();
    label_to_node_map = {};
    order_to_node_map = {};

    var t = new Tree();

    newick = newick.trim(newick);
    t.Parse(newick);

    for (var index = 1; index < id_to_node_map.length; index++)
    {
        label_to_node_map[id_to_node_map[index].label] = id_to_node_map[index];
    }

    // generate tooltip text before normalization
    metadata_dict = new Array();

    empty_tooltip = '<table>';
    empty_tooltip += '<tr><td>Name: </td><td>n/a</td></tr>';
    empty_tooltip += '<tr><td>Parent: </td><td>n/a</td></tr>';

    for (var i = 1; i < settings['layer-order'].length; i++)
    {
        var pindex = settings['layer-order'][i];
        empty_tooltip += '<tr><td>' + metadata[0][pindex] + '</td><td>n/a</td></tr>';
    }
    empty_tooltip += '</table>';

    $('#tooltip_content').html(empty_tooltip);
    $('#tooltip').dialog('open');

    for (var index = 1; index < metadata.length; index++) 
    {
        var params = metadata[index];
        metadata_dict[params[0]] = params.slice(0);

        var title = [];
        title.push('<td>Name: </td><td>' + metadata[index][0] + '</td>');
        for (var i = 0; i < settings['layer-order'].length; i++) 
        {
            var pindex = settings['layer-order'][i];

            if (layer_types[pindex] == 0) // check if parent
            {   
                if (metadata[index][pindex] == '')
                {
                    title.push('<td>Parent: </td><td>n/a</td>');
                }
                else
                {
                    title.push('<td>Parent: </td><td>' + metadata[index][pindex] + '</td>');
                }
                
            }
            else
            {
                title.push('<td>' + metadata[0][pindex] + '</td><td>' + metadata[index][pindex] + '</td>');
            }
        }

        metadata_title[params[0]] = title;
    }

    $('#draw_delta_time').html('tooltips ready (took <b>' + tree_draw_timer.getDeltaSeconds('tooltips')['deltaSecondsPrev'] + '</b> seconds).');
 
    // normalization
    var param_max = new Array();

    for (var id in metadata_dict) 
    {
        for (var pindex = 1; pindex < parameter_count; pindex++) 
        {
            var layer = settings['views'][current_view][pindex];

            if (layer_types[pindex] == 0 || layer_types[pindex] == 2) 
            {
                // skip normalization for parent layer & categorical data type
                continue;
            }    
            if (layer_types[pindex] == 1) // stack bar
            {
                // convert ";" string to array after normalization
                var stack_bar_items = metadata_dict[id][pindex].split(";");

                if (layer['normalization'] == 'sqrt') {
                    for (var j=0; j < stack_bar_items.length; j++)
                    {
                        stack_bar_items[j] = Math.sqrt(parseFloat(stack_bar_items[j]));
                    }
                }
                if (layer['normalization'] == 'log') {
                    for (var j=0; j < stack_bar_items.length; j++)
                    {
                        stack_bar_items[j] = log10(parseFloat(stack_bar_items[j]) + 1);
                    }
                }

                metadata_dict[id][pindex] = stack_bar_items.slice(0);
                continue;
            }

            // numerical data
            if (layer['normalization'] == 'sqrt') 
            {
                metadata_dict[id][pindex] = Math.sqrt(parseFloat(metadata_dict[id][pindex]));
            }
            if (layer['normalization'] == 'log') 
            {
                metadata_dict[id][pindex] = log10(parseFloat(metadata_dict[id][pindex]) + 1);
            }
            if (typeof param_max[pindex] === 'undefined' || parseFloat(metadata_dict[id][pindex]) > parseFloat(param_max[pindex])) 
            {
                param_max[pindex] = parseFloat(metadata_dict[id][pindex]);
            }
        }
    }

    $('#draw_delta_time').html('normalizations done (took <b>' + tree_draw_timer.getDeltaSeconds('normalizations')['deltaSecondsPrev'] + '</b> seconds).');

    // calculate bar sizes according to given height
    for (var pindex = 1; pindex < parameter_count; pindex++) 
    {
        if (layer_types[pindex] == 0 || layer_types[pindex] == 2) 
        {
            // skip normalization for parent layer & categorical data type
            continue;
        }    
        else
        {
            var layer = settings['views'][current_view][pindex];

            var min_max_disabled = layer['min']['disabled'];

            var min = parseFloat(layer['min']['value']);
            var max = parseFloat(layer['max']['value']);

            var min_new = null;
            var max_new = null;

            for (var id in metadata_dict)
            {
                if (layer_types[pindex] == 1) // stack bar
                {
                    var total = 0;

                    for (var j=0; j < metadata_dict[id][pindex].length; j++)
                    {
                        total = total + parseFloat(metadata_dict[id][pindex][j]);
                    }

                    var multiplier = parseFloat(layer['height']) / total;

                    for (var j=0; j < metadata_dict[id][pindex].length; j++)
                    {
                        metadata_dict[id][pindex][j] = metadata_dict[id][pindex][j] * multiplier;
                    }
                }
                else // numerical data
                {
                    var bar_size = parseFloat(metadata_dict[id][pindex]);
                                        
                    if (!min_max_disabled)
                    {
                        if (bar_size > max) {
                            bar_size = max - min;
                        }
                        else if (bar_size < min) {
                            bar_size = 0;
                        }
                        else {
                            bar_size = bar_size - min;
                        }

                        if (bar_size == 0) {
                            metadata_dict[id][pindex] = 0;
                        } else {
                            metadata_dict[id][pindex] = bar_size *  parseFloat(layer['height']) / (max - min);
                        }                        
                    }
                    else
                    {  
                        if ((min_new == null) || bar_size < min_new) {
                            min_new = bar_size;
                        }
                        
                        if ((max_new == null) || bar_size > max_new) {
                            max_new = bar_size;
                        }

                        if (bar_size == 0) {
                            metadata_dict[id][pindex] = 0;
                        } else {
                            metadata_dict[id][pindex] = bar_size *  parseFloat(layer['height']) / param_max[pindex];
                        }

                        var min_max_str = "Min: " + min_new + " - Max: " + max_new;
                        $('#min' + pindex).attr('title', min_max_str);
                        $('#max' + pindex).attr('title', min_max_str);


                    }
                }
            }

            if (min_max_disabled)
            {
                $('#min' + pindex).prop('disabled', false);
                $('#max' + pindex).val(max_new).prop('disabled', false);        
            }
        }
    }

    $('#draw_delta_time').html('barsizes done (took <b>' + tree_draw_timer.getDeltaSeconds('barsizes')['deltaSecondsPrev'] + '</b> seconds).');

    if (t.error != 0) {
        alert('Error parsing tree');
    } else {
        t.ComputeWeights(t.root);

        var td = null;

        // clear existing diagram, if any
        var svg = document.getElementById('svg');
        while (svg.hasChildNodes()) {
            svg.removeChild(svg.lastChild);
        }

        // create new group
        createGroup('svg', 'viewport');
        createGroup('viewport', 'tree_group');
        createGroup('tree_group', 'tree');

        var width = settings['tree-width'];
        var height = settings['tree-height'];

        if (width == 0)
            width = VIEWER_WIDTH;
        if (height == 0)
            height = VIEWER_HEIGHT;

        switch (settings['tree-type']) {
            case 'phylogram':
                if (t.has_edge_lengths) {
                    td = new PhylogramTreeDrawer();
                } else {
                    td = new RectangleTreeDrawer();
                }

                td.Init(t, {
                    svg_id: 'tree',
                    width: height,
                    height: width,
                    fontHeight: 10,
                    root_length: 0.1
                });
                $('#tree_group').attr('transform', 'rotate(90)'); // height and width swapped because of this.
                break;

            case 'circlephylogram':
                if (t.has_edge_lengths) {
                    td = new CirclePhylogramDrawer();
                } else {

                    td = new CircleTreeDrawer();
                }
                td.Init(t, {
                    svg_id: 'tree',
                    width: VIEWER_WIDTH,
                    height: VIEWER_HEIGHT,
                    fontHeight: 10,
                    root_length: 0.1
                });
                break;
        }

        td.CalcCoordinates();
        td.Draw();

        // calculate max radius of tree
        layer_boundaries = new Array();
        var tree_radius = 0;
        var tree_max_x = 0;
        var tree_max_y = 0;

        var n = new NodeIterator(t.root);
        var q = n.Begin();

        order_counter = 0;
        switch (settings['tree-type']) {
            case 'phylogram':
                while (q != null)
                {
                    if (q.xy.x > tree_max_x)
                        tree_max_x = q.xy.x;
                    if (q.xy.y > tree_max_y)
                        tree_max_y = q.xy.y;

                    // childs
                    var _n = new NodeIterator(q);
                    var _q = _n.Begin();

                    q.child_nodes = [];
                    while (_q != null) {
                        q.child_nodes.push(_q.id);
                        _q = _n.Next();
                    }
                    // end of childs
                    if (q.IsLeaf()) {
                        q.order = order_counter++;
                        order_to_node_map[q.order] = q;
                    }

                    q = n.Next();

                }
                layer_boundaries.push( [0, tree_max_x] );
                
                // calculate height per leaf
                height_per_leaf = width / (t.num_leaves - 1);

                break;

            case 'circlephylogram':
                while (q != null) 
                {
                    if (q.radius > tree_radius)
                        tree_radius = q.radius;

                    // childs
                    var _n = new NodeIterator(q);
                    var _q = _n.Begin();

                    q.child_nodes = [];
                    while (_q != null) {
                        q.child_nodes.push(_q.id);
                        _q = _n.Next();
                    }
                    // end of childs
                    if (q.IsLeaf()) {
                        q.order = order_counter++;
                        order_to_node_map[q.order] = q;
                    }

                    q = n.Next();
                }
                layer_boundaries.push( [0, tree_radius] );

                // calculate angle per leaf
                var angle_min = parseFloat(settings['angle-min']);
                var angle_max = parseFloat(settings['angle-max']);
                angle_per_leaf = Math.toRadians(angle_max - angle_min) / t.num_leaves;
                break;
        }
        leaf_count = Object.keys(order_to_node_map).length;
                
        margin = parseFloat(settings['layer-margin']);

        // calculate layer boundries
        for (var i = 0; i < settings['layer-order'].length; i++) {
            var layer_index = i+1;
            var pindex = settings['layer-order'][i];
            var layer = settings['views'][current_view][pindex];

            var layer_margin = (parseFloat(layer['height'])==0) ? 0 : margin;
            
            layer_boundaries.push( [ layer_boundaries[i][1] + layer_margin, layer_boundaries[i][1] + layer_margin + parseFloat(layer['height']) ] );

            createGroup('tree_group', 'layer_background_' + layer_index);
            createGroup('tree_group', 'layer_' + layer_index);
            createGroup('tree_group', 'event_catcher_' + layer_index);

            // draw event catcher of the layer
            if (settings['tree-type']=='phylogram')
            {
                var color = layer['color'];

                drawPhylogramRectangle('event_catcher_' + layer_index,
                    'event',
                    layer_boundaries[layer_index][0],
                    tree_max_y / 2,
                    tree_max_y + height_per_leaf,
                    layer_boundaries[layer_index][1] - layer_boundaries[layer_index][0],
                    '#ffffff',
                    0,
                    true);  
            }
            if (settings['tree-type']=='circlephylogram')
            {
                var color = layer['color'];

                var _min = Math.toRadians(settings['angle-min']) - (angle_per_leaf / 2);
                var _max = Math.toRadians(settings['angle-max']) - (angle_per_leaf / 2);

                drawPie('event_catcher_' + layer_index,
                    'event',
                    _min,
                    _max,
                    layer_boundaries[layer_index][0],
                    layer_boundaries[layer_index][1],
                    (_max - _min > Math.PI) ? 1:0, // large arc flag
                    '#ffffff',
                    0,
                    true);
            }

            if (settings['tree-type']=='phylogram' && layer_types[pindex] == 3) // draw numerical bar backgroung for phylogram
            {
                var color = layer['color'];

                drawPhylogramRectangle('layer_background_' + layer_index,
                    'all',
                    layer_boundaries[layer_index][0],
                    tree_max_y / 2,
                    tree_max_y + height_per_leaf,
                    layer_boundaries[layer_index][1] - layer_boundaries[layer_index][0],
                    color,
                    0.3,
                    false);
            }
            
            if (settings['tree-type']=='circlephylogram' && layer_types[pindex] == 3)
            {
                var color = layer['color'];

                var _min = Math.toRadians(settings['angle-min']) - (angle_per_leaf / 2);
                var _max = Math.toRadians(settings['angle-max']) - (angle_per_leaf / 2);

                drawPie('layer_background_' + layer_index,
                    'all',
                    _min,
                    _max,
                    layer_boundaries[layer_index][0],
                    layer_boundaries[layer_index][1],
                    (_max - _min > Math.PI) ? 1:0, // large arc flag
                    color,
                    0.3,
                    false);
            }
        }

        total_radius = layer_boundaries[layer_boundaries.length - 1][1];
        beginning_of_layers = layer_boundaries[0][1];

        var n = new NodeIterator(t.root);
        var q = n.Begin();

        createGroup('tree_group', 'guide_lines');

        // parent things
        var parent_odd = '#888888';
        var parent_even = '#666666';
        var parent_residual = '#AAAAAA';
        var prev_parent_color = parent_odd;
        var prev_parent_name = '';
        var prev_parent_items = new Array();
        var parent_count = 0;

        odd_even_flag = -1

        var categorical_layers_ordered = {};

        for (var i = 0; i < settings['layer-order'].length; i++) {
            var layer_index = i+1;
            var pindex = settings['layer-order'][i];

            if (layer_types[pindex] == 2 || layer_types[pindex] == 0) // categorical or parent
            {
                categorical_layers_ordered[layer_index] = new Array();
            }
        }

        switch (settings['tree-type']) {
            case 'phylogram':
                while (q != null) {
                    if (q.IsLeaf()) {
                        odd_even_flag = odd_even_flag * -1;

                        if (odd_even_flag > 0)
                            drawStraightGuideLine('guide_lines', q.id, q.xy, tree_max_x);

                        for (var i = 0; i < settings['layer-order'].length; i++) {
                            var layer_index = i+1;
                            var pindex = settings['layer-order'][i];
                            var layer = settings['views'][current_view][pindex];

                            var isParent      = (layer_types[pindex] == 0) ? true : false;
                            var isStackBar    = (layer_types[pindex] == 1) ? true : false;
                            var isCategorical = (layer_types[pindex] == 2) ? true : false;

                            if(isStackBar)
                            {
                                var offset = 0;
                                for (var j=0; j < metadata_dict[q.label][pindex].length; j++)
                                {
                                    drawPhylogramRectangle('layer_' + layer_index,
                                        q.id,
                                        layer_boundaries[layer_index][1] - offset - metadata_dict[q.label][pindex][j],
                                        q.xy['y'],
                                        height_per_leaf,
                                        metadata_dict[q.label][pindex][j],
                                        stack_bar_colors[pindex][j],
                                        1,
                                        false);
                                    offset += metadata_dict[q.label][pindex][j];
                                } 
                        
                            }
                            else if(isCategorical)
                            {
                                if (typeof categorical_data_colors[pindex][metadata_dict[q.label][pindex]] === 'undefined'){
                                    if (typeof(metadata_dict[q.label][pindex]) == typeof(null))
                                        categorical_data_colors[pindex][metadata_dict[q.label][pindex]] = '#ffffff';
                                    else
                                        categorical_data_colors[pindex][metadata_dict[q.label][pindex]] = randomColor();
                                }

                                categorical_layers_ordered[layer_index].push(metadata_dict[q.label][pindex]);
                            }
                            else if (isParent)
                            {
                                categorical_layers_ordered[layer_index].push(metadata_dict[q.label][pindex]);
                            }
                            else // numerical
                            {
                                if (metadata_dict[q.label][pindex] > 0) {
                                    var color = layer['color'];

                                     drawPhylogramRectangle('layer_' + layer_index,
                                        q.id,
                                        layer_boundaries[layer_index][1] - metadata_dict[q.label][pindex],
                                        q.xy['y'],
                                        height_per_leaf,
                                        metadata_dict[q.label][pindex],
                                        color,
                                        1,
                                        false);
                                }
                            }
                        }
                    }
                    q = n.Next();
                }
                break;
            case 'circlephylogram':
                while (q != null) {
                    if (q.IsLeaf()) {
                        odd_even_flag = odd_even_flag * -1;

                        // draw guidelines for every other leaf.
                        if (odd_even_flag > 0)
                            drawGuideLine('guide_lines', q.id, q.angle, q.radius, beginning_of_layers);

                        for (var i = 0; i < settings['layer-order'].length; i++) {
                            var layer_index = i+1;
                            var pindex = settings['layer-order'][i];
                            var layer = settings['views'][current_view][pindex];

                            var isParent      = (layer_types[pindex] == 0) ? true : false;
                            var isStackBar    = (layer_types[pindex] == 1) ? true : false;
                            var isCategorical = (layer_types[pindex] == 2) ? true : false;

                            if(isStackBar)
                            {
                                var offset = 0;
                                for (var j=0; j < metadata_dict[q.label][pindex].length; j++)
                                {
                                    drawPie('layer_' + layer_index,
                                        q.id,
                                        q.angle - angle_per_leaf / 2,
                                        q.angle + angle_per_leaf / 2,
                                        layer_boundaries[layer_index][0] + offset,
                                        layer_boundaries[layer_index][0] + offset + metadata_dict[q.label][pindex][j],
                                        0,
                                        stack_bar_colors[pindex][j],
                                        1,
                                        false);
                                    offset += metadata_dict[q.label][pindex][j];
                                } 
                        
                            }
                            else if(isCategorical)
                            {
                                if (typeof categorical_data_colors[pindex][metadata_dict[q.label][pindex]] === 'undefined'){
                                    if (typeof(metadata_dict[q.label][pindex]) == typeof(null))
                                        categorical_data_colors[pindex][metadata_dict[q.label][pindex]] = '#ffffff';
                                    else
                                        categorical_data_colors[pindex][metadata_dict[q.label][pindex]] = randomColor();
                                }

                                categorical_layers_ordered[layer_index].push(metadata_dict[q.label][pindex]);
                            }
                            else if (isParent)
                            {
                                categorical_layers_ordered[layer_index].push(metadata_dict[q.label][pindex]);
                            }
                            else // numerical
                            {
                                var color = layer['color'];

                                if (metadata_dict[q.label][pindex] > 0) {
                                    drawPie('layer_' + layer_index,
                                        q.id,
                                        q.angle - angle_per_leaf / 2,
                                        q.angle + angle_per_leaf / 2,
                                        layer_boundaries[layer_index][0], 
                                        layer_boundaries[layer_index][0] + metadata_dict[q.label][pindex],
                                        0,
                                        color,
                                        1,
                                        false);
                                }
                            }

                        }
                    }
                    q = n.Next();
                }
                break;
        }

        // cluster categorical items and draw them
        
        for (var i = 0; i < settings['layer-order'].length; i++) {
            var layer_index = i+1;
            var pindex = settings['layer-order'][i];

            var isParent      = (layer_types[pindex] == 0) ? true : false;
            var isCategorical = (layer_types[pindex] == 2) ? true : false;

            if (isParent || isCategorical)
            {
                var layer_items = categorical_layers_ordered[layer_index];
                layer_items.push(null);

                var prev_value = layer_items[0];
                var prev_start = 0;

                var items_to_draw = new Array();

                for (var j=1; j < layer_items.length; j++)
                {
                    if (prev_value != layer_items[j])
                    {
                        if (prev_value != null && prev_value != '')
                            items_to_draw.push(new Array(prev_start, j - 1, prev_value)); // start, end, item;

                        prev_start = j;
                    }
                    prev_value = layer_items[j];
                }

                for (var j=0; j < items_to_draw.length; j++)
                {
                    var categorical_item = items_to_draw[j];

                    var color;

                    if (isCategorical)
                    {
                        color = categorical_data_colors[pindex][categorical_item[2]];
                    } 
                    else // parent
                    {
                        if (j % 2 == 1 && j == items_to_draw.length - 1)
                            color = '#AAAAAA';
                        else if (j % 2 == 1)
                            color = '#888888';
                        else
                            color = '#666666';
                    }

                    var start = order_to_node_map[categorical_item[0]];
                    var end = order_to_node_map[categorical_item[1]];

                    if (tree_type == 'circlephylogram')
                    {
                        drawPie('layer_' + layer_index,
                            'categorical_' + layer_index + '_' + j, // path_<layer>_<id>
                            start.angle - angle_per_leaf / 2,
                            end.angle + angle_per_leaf / 2,
                            layer_boundaries[layer_index][0], 
                            layer_boundaries[layer_index][1],
                            (end.angle - start.angle + angle_per_leaf > Math.PI) ? 1 : 0,
                            color,
                            1,
                            false);
                    }
                    else
                    {
                        var height = end.xy['y'] - start.xy['y'];

                        drawPhylogramRectangle('layer_' + layer_index,
                            'categorical_' + layer_index + '_' + j, // path_<layer>_<id>
                            layer_boundaries[layer_index][0],
                            start.xy['y'] + height / 2,
                            height + height_per_leaf,
                            layer_boundaries[layer_index][1] - layer_boundaries[layer_index][0],
                            color,
                            1,
                            false);
                    }
                }
            }
        }

        rebuildIntersections();
        createGroup('tree_group', 'group');
        redrawGroups();

        // observe for transform matrix change
        var observer = new MutationObserver(updateSingleBackgroundGlobals);
    
        observer.observe(document.getElementById('viewport'), {
            attributes:    true,
            attributeFilter: ["transform"]
        });

        // draw title
        switch (settings['tree-type']) {
            case 'phylogram':
                drawLegend(0, 100);
                drawText('viewport', {'x': -0.5 * td.settings.height, 'y': -150}, document.title, '72px', 'center');
                break;
            case 'circlephylogram':
                drawLegend(20 - total_radius, total_radius + 100);
                drawText('viewport', {'x': 0, 'y': -1 * total_radius - 150}, document.title, '72px', 'center');
                break;
        }

        // draw layer names (in circlephylogram with 0-270)
        if (settings['tree-type'] == 'circlephylogram' && settings['angle-min'] == 0 && settings['angle-max'] == 270)
        {
            for (var i = 0; i < settings['layer-order'].length; i++) {
                var layer_index = i+1;
                var pindex = settings['layer-order'][i];
                var layer = settings['views'][current_view][pindex];

                var layer_title = metadata[0][pindex];

                if (layer_title in named_layers && 'pretty_name' in named_layers[layer_title]) {
                    layer_title = named_layers[layer_title]['pretty_name'];
                } else {
                    layer_title = layer_title.replace(/_/g, " ");
                }

                drawText('tree_group', {
                    'x': 10,
                    'y': 0 - (layer_boundaries[layer_index][1] + layer_boundaries[layer_index][0]) / 2
                }, layer_title , layer['height'] + 'px', 'left', layer['color']);
            }
        }

        // Scale to fit window
        var bbox = svg.getBBox();
        scale = Math.min(td.settings.width / bbox.width, td.settings.height / bbox.height) * 0.80;

        zoom_reset();

        // pan and mouse zoom
        $('svg').svgPan('viewport');

        // tree group = tree + layers
        // bind events
        var tree_group = document.getElementById('tree_group');
        tree_group.addEventListener('click', lineClickHandler, false);
        tree_group.addEventListener('contextmenu', lineContextMenuHandler, false);
        tree_group.addEventListener('mouseover',lineMouseEnterHandler, false);
        tree_group.addEventListener('mouseout', lineMouseLeaveHandler, false);
        document.body.addEventListener('mousemove', mouseMoveHandler, false); // for tooltip
        document.body.addEventListener('click', function() { $('#control_contextmenu').hide(); }, false);

        // code below required to stop clicking on contings while panning.
        var viewport = document.getElementById('svg');
        viewport.addEventListener('mousedown', 
            function(event) { 
                dragging = false; 

                if (shiftPressed)
                {
                    drawing_zoom = true;

                    zoomBox['start_x'] = event.clientX;
                    zoomBox['start_y'] = event.clientY;

                    $('#divzoom').css({"top": 0, "left": 0, "width": 0, "height": 0 });
                    $('#divzoom').show();
                }
            });
        viewport.addEventListener('mousemove', 
            function(event) { 
                dragging = true; 

                if (shiftPressed && drawing_zoom)
                {
                    var _top = zoomBox['start_y'] > event.clientY ? event.clientY : zoomBox['start_y'];
                    var _left = zoomBox['start_x'] > event.clientX ? event.clientX : zoomBox['start_x'];
                    var _height = Math.abs(zoomBox['start_y'] - event.clientY);
                    var _width = Math.abs(zoomBox['start_x'] - event.clientX);

                    var divzoom = document.getElementById('divzoom');
                    var divzoom_inner = document.getElementById('divzoom_inner');

                    divzoom.style.top = _top + "px";
                    divzoom.style.left = _left + "px";
                    divzoom.style.width = _width + "px";
                    divzoom.style.height = _height + "px";

                    var w_ratio = _width / VIEWER_WIDTH;
                    var h_ratio = _height / VIEWER_HEIGHT;

                    if (w_ratio > h_ratio)
                    {
                        divzoom_inner.style.width = h_ratio * VIEWER_WIDTH + "px";
                        divzoom_inner.style.height = _height + "px";
                    }
                    else
                    {
                        divzoom_inner.style.width = _width + "px";
                        divzoom_inner.style.height = w_ratio * VIEWER_HEIGHT + "px";
                    }

                    // when you drawing rectangle, if you drag over text on the screen browser selects that text
                    // with this hack you can continue drawing.
                    clearTextSelection(); // in utils.js

                }
            });

        viewport.addEventListener('mouseup', 
            function() {
                if (drawing_zoom)
                {
                    var inner_rect = document.getElementById('divzoom_inner').getBoundingClientRect();
                    
                    if (inner_rect.width > 0 && inner_rect.height > 0)
                    {
                        var _dx = (VIEWER_WIDTH / 2) - (inner_rect.left + inner_rect.width / 2);
                        var _dy = (VIEWER_HEIGHT / 2) - (inner_rect.top + inner_rect.height / 2);
                        pan(_dx,_dy);
                        zoom(VIEWER_WIDTH / inner_rect.width);
                    }
                }

                drawing_zoom=false; 
                zoomBox = {}; 
                $('#divzoom').hide(); 

            });
    }

    $('#draw_delta_time').html('drawn in ' + tree_draw_timer.getDeltaSeconds('done')['deltaSecondsStart'] + ' seconds.');
}

function redrawGroups(search_results)
{
    // check if tree parsed, if not there is nothing to redraw.
    if ($.isEmptyObject(label_to_node_map)) 
        return;

    var leaf_list = Array.apply(null, new Array(order_counter+1)).map(Number.prototype.valueOf,0);

    // put group numbers of selected leaves to leaf list
    // maybe we should write directly into leaf_list in mouse events, instead of generate it everytime.
    for (var gid = 1; gid <= group_counter; gid++) {
        for (var j = 0; j < SELECTED[gid].length; j++) {
            if (label_to_node_map[SELECTED[gid][j]].IsLeaf()) {
                leaf_list[label_to_node_map[SELECTED[gid][j]].order] = gid;
            }
        }
    }

    // cluster groups and put them into groups_to_draw array with (start, end, gid);
    var prev_value = leaf_list[0];
    var prev_start = 0;

    var groups_to_draw = new Array();

    for (var i=1; i < leaf_list.length; i++)
    {
        if (prev_value != leaf_list[i])
        {
            if (prev_value != 0)
                groups_to_draw.push(new Array(prev_start, i - 1, prev_value)); // start, end, gid;

            prev_start = i;
        }
        prev_value = leaf_list[i];
    }

    // remove exist group drawings
    var group = document.getElementById('group');
    while (group.hasChildNodes()) {
        group.removeChild(group.lastChild);
    }

    // draw new groups
    for (var i=0; i < groups_to_draw.length; i++) {
        var start = order_to_node_map[groups_to_draw[i][0]];
        var end = order_to_node_map[groups_to_draw[i][1]];

        var color = document.getElementById('group_color_' + groups_to_draw[i][2]).getAttribute('color');
        var outer_ring_size = 4;

        if (tree_type == 'circlephylogram')
        {
            drawPie('group',
                'group_background_' + i,
                start.angle - angle_per_leaf / 2,
                end.angle + angle_per_leaf / 2,
                beginning_of_layers,
                total_radius,
                (end.angle - start.angle + angle_per_leaf > Math.PI) ? 1 : 0,
                color,
                0.1,
                false);

            drawPie('group',
                'group_outer_' + 1,
                start.angle - angle_per_leaf / 2,
                end.angle + angle_per_leaf / 2,
                total_radius + margin,
                total_radius + margin * outer_ring_size,
                (end.angle - start.angle + angle_per_leaf > Math.PI) ? 1 : 0,
                color,
                1,
                false);     
        }
        else
        {

            drawPhylogramRectangle('group',
                'group_background_' + i,
                beginning_of_layers,
                (start.xy.y + end.xy.y) / 2,
                end.xy.y - start.xy.y + height_per_leaf,
                total_radius + margin - beginning_of_layers,
                color,
                0.1,
                false);

            drawPhylogramRectangle('group',
                'group_outer_' + 1,
                total_radius + margin,
                (start.xy.y + end.xy.y) / 2,
                end.xy.y - start.xy.y + height_per_leaf,
                margin * outer_ring_size,
                color,
                1,
                false);
        }
    }


    // draw search results
    if (typeof search_results !== 'undefined')
    {
        for (var i=0; i < search_results.length; i++) {
            var start = order_to_node_map[search_results[i]];
            var end = start;

            var color = document.getElementById('picker_highlight').getAttribute('color');
            var outer_ring_size = 6;

            if (tree_type == 'circlephylogram')
            {
                drawPie('group',
                    'group_outer_' + 1,
                    start.angle - angle_per_leaf / 2,
                    end.angle + angle_per_leaf / 2,
                    total_radius + margin,
                    total_radius + margin * outer_ring_size,
                    (end.angle - start.angle + angle_per_leaf > Math.PI) ? 1 : 0,
                    color,
                    1,
                    false);     
            }
            else
            {
                drawPhylogramRectangle('group',
                    'group_outer_' + 1,
                    total_radius + margin,
                    (start.xy.y + end.xy.y) / 2,
                    end.xy.y - start.xy.y + height_per_leaf,
                    margin * outer_ring_size,
                    color,
                    1,
                    false);
            }
        }
    }
}


function rebuildIntersections()
{
    for (var gid = 1; gid <= group_counter; gid++) {

        // delete extra intersections
        var deleted;
        do {
            deleted = 0;
            var cursor = SELECTED[gid].length;
            while (cursor--)
            {
                var node = label_to_node_map[SELECTED[gid][cursor]];
                
                if (node.IsLeaf())
                    continue;

                if (node.child != null && SELECTED[gid].indexOf(node.child.label) > -1)
                    continue;

                if (node.sibling != null && SELECTED[gid].indexOf(node.sibling.label) > -1)
                    continue;

                SELECTED[gid].splice(cursor,1);
                deleted++;
            }

        } while (deleted > 0)

        // try to make new intersections
        var inserted;
        do {
            inserted = 0;
            var length = SELECTED[gid].length;
            for (var cursor = 0; cursor < length; cursor++)
            {
                var node = label_to_node_map[SELECTED[gid][cursor]];
                var parent = node.ancestor;

                if (parent.ancestor == null) 
                {
                    // skip root
                    continue;
                }

                if (SELECTED[gid].indexOf(parent.label) > -1)
                {
                    // parent already in selected list
                    continue;
                }

                if (node.sibling != null && SELECTED[gid].indexOf(node.sibling.label) > -1)
                {
                    SELECTED[gid].push(parent.label);
                    
                    inserted++;
                }
            }

        } while (inserted > 0)
    }
}
