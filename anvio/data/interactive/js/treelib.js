/**
 * Javascript library to display phylogenetic trees
 *
 *  Author: Özcan Esen <ozcanesen@gmail.com>
 *  Credits: A. Murat Eren
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
function createBin(parent, bin_id) {
    var svgObject = document.getElementById(parent);
    var g = document.createElementNS('http://www.w3.org/2000/svg', 'g');
    g.setAttribute('id', bin_id);
    svgObject.appendChild(g);
}

//--------------------------------------------------------------------------------------------------
function drawTitle(top, left, settings) {
    var _font_size = (total_radius / 25);
    top -= 3 * _font_size;
    drawText('viewport', {'x': left, 'y': top}, document.title, 2 * _font_size + 'px', 'center');
    top += 2 * _font_size;
    var _sub_title = "Tree order: " + getClusteringPrettyName(settings['order-by']) + " | ";
    _sub_title    += "Current view: " + settings['current-view'] + " | ";
    _sub_title    += "Sample order: " + settings['samples-order'];

    drawText('viewport', {'x': left, 'y': top}, _sub_title, _font_size + 'px', 'center');
}

//--------------------------------------------------------------------------------------------------
function drawLegend(top, left, line_end) {
    var _left = left;
    var line_height = (line_end - _left) / 80;
    var gap = line_height / 2;

    var top = top + line_height * 3;
    var _top = top;

    var legends = [];
   
    $.each(layer_types, function (i, _) {
        var pindex = i;
        if (layer_types[pindex] != 2)
            return; // skip if not categorical

        if (layers[pindex]['type'] == 'text')
            return; // skip if type is text

        var categorical_stats = {};

        categorical_stats['None'] = 0;
        for (var name in categorical_data_colors[pindex]) {
            categorical_stats[name] = 0;
        }
        for (var index = 1; index < layerdata.length; index++)
        {
            var taxonomy_name = layerdata[index][pindex];
            if (taxonomy_name == null || taxonomy_name == '' || taxonomy_name == 'null')
                taxonomy_name = 'None';
            categorical_stats[taxonomy_name] += 1;
        }
        var names = Object.keys(categorical_stats).sort(function(a,b){return categorical_stats[b]-categorical_stats[a]});

        names.push(names.splice(names.indexOf('None'), 1)[0]); // move null and empty categorical items to end

        legends.push({
            'name': getLayerName(pindex),
            'source': 'categorical_data_colors',
            'key': pindex,
            'item_names': names,
            'item_keys': names,
            'stats': categorical_stats
        });
    });

    $.each(layer_types, function (i, _) {
        if (layer_types[i] != 1)
            return; // skip if not stack bar

        var pindex = i;
        var layer_name = getLayerName(pindex);
        var names = (layer_name.indexOf('!') > -1) ? layer_name.split('!')[1].split(';') : layer_name.split(';');
        var keys = Array.apply(null, Array(names.length)).map(function (_, i) {return i;});

        var pretty_name = getLayerName(pindex);
        pretty_name = (pretty_name.indexOf('!') > -1) ? pretty_name.split('!')[0] : pretty_name;

        legends.push({
            'name': pretty_name,
            'source': 'stack_bar_colors',
            'key': pindex,
            'item_names': names,
            'item_keys': keys,
        });
    });

    for (sample in samples_categorical_colors)
    {
        var names = Object.keys(samples_categorical_colors[sample]);

        legends.push({
            'name': sample,
            'source': 'samples_categorical_colors',
            'key': sample,
            'item_names': names,
            'item_keys': names,
            //'stats': TO DO
        });
    }

    for (sample in samples_stack_bar_colors)
    {
        var names = (sample.indexOf('!') > -1) ? sample.split('!')[1].split(';') : sample.split(';');
        var keys = Array.apply(null, Array(names.length)).map(function (_, i) {return i;});
        var pretty_name = (sample.indexOf('!') > -1) ? sample.split('!')[0] : sample;

        legends.push({
            'name': pretty_name,
            'source': 'samples_stack_bar_colors',
            'key': sample,
            'item_names': names,
            'item_keys': keys,
            //'stats': TO DO
        });
    }

    for (var i=0; i < legends.length; i++)
    {
        var bin_id = 'legend_' + i;
        createBin('viewport', bin_id);
        var legend = legends[i];

        for (var j = 0; j < legend['item_names'].length; j++) {
            var _name = legend['item_names'][j];
            if (legend.hasOwnProperty('stats') && legend['stats'][_name] == 0)
                continue;

            if (_left > line_end)
            {
                _left = left;
                _top = _top + line_height + gap;
            }
            var rect = drawRectangle(bin_id, _left, _top, line_height, line_height, window[legend['source']][legend['key']][legend['item_keys'][j]], 1, 'black',
                null,
                function() {
                    // mouseenter
                    $(this).css('stroke-width', '2');
                },
                function() {
                    // mouseleave
                    $(this).css('stroke-width', '1');
                });

            rect.setAttribute('callback_source', legend['source']);
            rect.setAttribute('callback_pindex', legend['key']);
            rect.setAttribute('callback_name', legend['item_keys'][j]);

            $(rect).colpick({
                layout: 'hex',
                submit: 0,
                colorScheme: 'dark',
                onChange: function(hsb, hex, rgb, el, bySetColor) {
                    $(el).css('fill', '#' + hex);
                    window[el.getAttribute('callback_source')][el.getAttribute('callback_pindex')][el.getAttribute('callback_name')] = '#' + hex;
                }
            });

            _left += line_height + gap;

            if (legend.hasOwnProperty('stats'))
            {
                _name = ((_name == 'null' || _name == '') ? 'None' : _name) + ' (' + legend['stats'][_name] + ')';
            }

            var text = drawText(bin_id, {
                'x': _left,
                'y': _top + (line_height * 3/4)
            }, _name, line_height + 'px', 'left', '#000000', 'baseline');

            _left += text.getBBox().width + gap;
        }
        drawText(bin_id, {
            'x': left - line_height,
            'y': (_top + top + line_height) / 2 
        }, legend['name'], 2*line_height + 'px', 'right');

        _top = _top + 3 * line_height + gap;
        top = _top;
        _left = left;
    }
}

function drawBinLegend(bins_to_draw, top, left) {
    createBin('viewport', 'bin_legend');

    drawRectangle('bin_legend', left - 10, top - 20,20 + (bins_to_draw.length + 2.5) * 20, 300, 'white', 1, 'black');
    drawText('bin_legend', {
        'x': left,
        'y': top
    }, "Bins", '16px');

    // table titles
    top = top + 28;
    drawText('bin_legend', {'x': left, 'y': top }, 'Color', '10px');
    drawText('bin_legend', {'x': left + 30, 'y': top}, 'Name', '10px');
    drawText('bin_legend', {'x': left + 170, 'y': top}, 'Contigs', '10px');
    drawText('bin_legend', {'x': left + 230, 'y': top}, 'Length', '10px');


    for (var bin_id=0; bin_id < bins_to_draw.length; bin_id++) {
        var bin = bins_to_draw[bin_id];
        top = top + 20;

        drawRectangle('bin_legend', left, top-8, 16, 16, bin['color'], 1, 'black');
        drawText('bin_legend', {'x': left + 30, 'y': top }, bin['name'], '12px');
        drawText('bin_legend', {'x': left + 170, 'y': top}, bin['contig-count'], '12px');
        drawText('bin_legend', {'x': left + 230, 'y': top}, bin['contig-length'], '12px');
    }
}

function drawLayerLegend(_layers, _view, _layer_order, top, left) {
    createBin('viewport', 'layer_legend');

    // legend border
    drawRectangle('layer_legend', left - 10, top - 20,20 + (_layer_order.length + 2.5) * 20, 300, 'white', 1, 'black');
    
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
    for (var i = 0; i < _layer_order.length; i++) 
    {
        var pindex = _layer_order[i];

        var layer_settings = _layers[pindex];
        var layer_view = _view[pindex];

        var layer_name = getLayerName(pindex);
        var short_name = (layer_name.length > 10) ? layer_name.slice(0,10) + "..." : layer_name;

        top = top + 20;

        // color
        if (layer_settings.hasOwnProperty('color') && typeof layer_settings['color'] != 'undefined') 
            drawRectangle('layer_legend', left, top - 8, 16, 16, layer_settings['color'], 1, 'black');

        // name
        drawText('layer_legend', {'x': left + 30, 'y': top}, short_name, '12px');

        // normalization
        if (layer_view.hasOwnProperty('normalization') && typeof layer_view['normalization'] != 'undefined') {
            var _norm = layer_view['normalization'];   
        } else {
            var _norm = "-";
        }
        drawText('layer_legend', {'x': left + 120, 'y': top}, _norm, '12px');

        // height
        drawText('layer_legend', {'x': left + 160, 'y': top}, layer_settings['height'], '12px');

        // min & max
        if (layer_view['min'].hasOwnProperty('value') && typeof layer_view['min']['value'] != 'undefined') {
            var _min = layer_view['min']['value'];
            var _max = layer_view['max']['value'];

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

    return line;
}

//--------------------------------------------------------------------------------------------------
function drawText(svg_id, p, string, font_size, align, color, baseline) {

    var text = document.createElementNS('http://www.w3.org/2000/svg', 'text');
    //newLine.setAttribute('id','node' + p.id);

    if (typeof color !== 'undefined')
        text.setAttribute('fill', color);

    if (typeof baseline === 'undefined')
        baseline = 'middle';


    text.setAttribute('style', 'alignment-baseline:' + baseline + '; font-size:' + font_size);
    text.setAttribute('x', p['x']);
    text.setAttribute('y', p['y']);
    text.setAttribute('pointer-events', 'none');
    text.setAttribute('text-rendering', 'optimizeLegibility');
    text.setAttribute('font-family', 'Helvetica Neue, Helvetica, Arial, sans-serif;')

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

    return text;
}

//--------------------------------------------------------------------------------------------------
function drawFixedWidthText(svg_id, p, string, font_size, color, width, height) {
    var textObj = drawText(svg_id, p, string, font_size, 'left', color, 'baseline');

    //https://stackoverflow.com/questions/9241315/trimming-text-to-a-given-pixel-width-in-svg
    if (textObj.getSubStringLength(0,string.length)>=width){
        for (var x=string.length-3;x>0;x-=3){
            if (textObj.getSubStringLength(0,x)<=width){
                textObj.textContent=string.substring(0,x)+"...";
                return;
            }
        }
        textObj.textContent="..."; //can't place at all
    }

    return textObj
}

//--------------------------------------------------------------------------------------------------
function drawRotatedText(svg_id, p, string, angle, align, font_size, font_family, color, maxLength) {
    var text = document.createElementNS('http://www.w3.org/2000/svg', 'text');
    //newLine.setAttribute('id','node' + p.id);
    text.setAttribute('style', 'alignment-baseline:middle');
    text.setAttribute('x', p['x']);
    text.setAttribute('y', p['y']);
    text.setAttribute('pointer-events', 'none');
    text.setAttribute('font-size', font_size);
    text.setAttribute('font-family', font_family);
    text.setAttribute('fill', color);

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

    var textNode = document.createTextNode(string);
    text.appendChild(textNode);
    var svg = document.getElementById(svg_id);
    svg.appendChild(text);

    // trim long text

    if (maxLength > 0)
    {
        for (var x=0; x < string.length; x++){
            if (text.getSubStringLength(0, x + 1) >= maxLength){
                text.textContent = string.substring(0, x);
                break;
            }
        }
    }
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

/*
    start_angle = start_angle - (angle_per_leaf / 40);
    end_angle = end_angle + (angle_per_leaf / 40);
*/
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
    pie.setAttribute('shape-rendering', 'auto');
    pie.setAttribute('stroke-width', '0');
    pie.setAttribute('d', path.join(" "));
    pie.setAttribute('fill-opacity', fill_opacity);

    if (!pointer_events)
        pie.setAttribute('pointer-events', 'none');

    var svg = document.getElementById(svg_id);
    svg.appendChild(pie);

    return pie;
}

function drawPhylogramRectangle(svg_id, id, x, y, height, width, color, fill_opacity, pointer_events) {
    var rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
    //height = height + height_per_leaf / 20;
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

    return rect;
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
    return node;
}

//--------------------------------------------------------------------------------------------------
Tree.prototype.Parse = function(str, edge_length_norm) {
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
    p.angle = (angle_per_leaf / 2) + angle_min + this.leaf_angle * this.leaf_count;
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
    // var viewport = document.getElementById(this.settings.svg_id);
    // viewport.setAttribute('transform', 'translate(' + (this.settings.width + this.root_length) / 2 + ' ' + this.settings.height / 2 + ')');
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

    var width = parseFloat(settings['tree-width']);
    var height = parseFloat(settings['tree-height']);
    var radius = parseFloat(settings['tree-radius']);

    if (width == 0)
        width = VIEWER_WIDTH;
    if (height == 0)
        height = VIEWER_HEIGHT;

    if (radius == 0)
        radius = (height > width) ? height : width;

    id_to_node_map = new Array();
    label_to_node_map = {};
    order_to_node_map = {};

    if (clusteringData.constructor === Array)
    {
        leaf_count = clusteringData.length;
        hasTree = false;
        var t = new Tree();
        t.root = clusteringData;

        var angle_min = parseFloat(settings['angle-min']);
        var angle_max = parseFloat(settings['angle-max']);
        height_per_leaf = width / (leaf_count - 1);
        angle_per_leaf = Math.toRadians(angle_max - angle_min) / leaf_count;

        for (var i=0; i < leaf_count; i++)
        {
            var leaf_name = clusteringData[i];
            var leaf_node = new Node(leaf_name);
            
            leaf_node.id = i+1;
            leaf_node.order = i;
            leaf_node.child_nodes = [leaf_node.id];

            if (settings['tree-type']=='phylogram') {
                leaf_node.xy = {};
                leaf_node.xy.x = 0;
                leaf_node.xy.y = height_per_leaf * i;
            } 
            else 
            {
                leaf_node.angle = angle_per_leaf / 2 + Math.toRadians(angle_min) + angle_per_leaf * i;
                pt = [];
                pt['x'] = radius * Math.cos(leaf_node.angle);
                pt['y'] = radius * Math.sin(leaf_node.angle);
                leaf_node.backarc = [];
                leaf_node.backarc['x'] = pt['x'];
                leaf_node.backarc['y'] = pt['y'];
            }

            id_to_node_map[leaf_node.id] = leaf_node;
            label_to_node_map[leaf_name] = leaf_node;
            order_to_node_map[i] = leaf_node;

            t.nodes[i] = leaf_node;
            t.num_leaves++;
        }
    }
    else
    {
        hasTree = true;
        var t = new Tree();

        clusteringData = clusteringData.trim(clusteringData);
        t.Parse(clusteringData, settings['edge-normalization']);

        if (t.error != 0) {
            toastr.error('Error while parsing tree data.');
            return;
        }

        var n = new NodeIterator(t.root);
        var q = n.Begin();

        order_counter = 0;
        while (q != null)
        {
            label_to_node_map[q.label] = q;
            id_to_node_map[q.id] = q;
            
            if (q.IsLeaf()) {
                q.order = order_counter++;
                order_to_node_map[q.order] = q;
            }

            q=n.Next();
        }
        leaf_count = t.num_leaves;
    }


    // generate tooltip text before normalization
    layerdata_dict = new Array();

    empty_tooltip = '<tr><td>split_name</td><td>n/a</td></tr>';
    empty_tooltip += '<tr><td>parent</td><td>n/a</td></tr>';

    for (var i = 1; i < settings['layer-order'].length; i++)
    {
        var pindex = settings['layer-order'][i];
        var layer_title = layerdata[0][pindex];
        if (layer_title.indexOf('!') > -1)
            layer_title = layer_title.split('!')[0];
        empty_tooltip += '<tr><td>' + layer_title + '</td><td>n/a</td></tr>';
    }

    $('#tooltip_content').html(empty_tooltip);
    $('#mouse_hover_scroll').css('top', '0');

    for (var index = 1; index < layerdata.length; index++) 
    {
        var params = layerdata[index];
        layerdata_dict[params[0]] = params.slice(0);

        var title = [];
        title.push('<td>split_name</td><td>' + layerdata[index][0] + '</td>');
        for (var i = 0; i < settings['layer-order'].length; i++) 
        {
            var pindex = settings['layer-order'][i];

            if (layer_types[pindex] == 0) // check if parent
            {   
                if (layerdata[index][pindex] == '')
                {
                    title.push('<td>parent</td><td>n/a</td>');
                }
                else
                {
                    title.push('<td>parent</td><td>' + layerdata[index][pindex] + '</td>');
                }
                
            }
            else
            {
                var layer_title = layerdata[0][pindex];
                if (layer_title.indexOf('!') > -1)
                    layer_title = layer_title.split('!')[0];

                title.push('<td>' + layer_title + '</td><td>' + layerdata[index][pindex] + '</td>');
            }
        }

        layerdata_title[params[0]] = title;
    }

    $('#draw_delta_time').html('tooltips ready (took <b>' + tree_draw_timer.getDeltaSeconds('tooltips')['deltaSecondsPrev'] + '</b> seconds).');
 
    // normalization
    var param_max = new Array();

    for (var id in layerdata_dict) 
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
                var stack_bar_items = layerdata_dict[id][pindex].split(";");

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

                layerdata_dict[id][pindex] = stack_bar_items.slice(0);
                continue;
            }

            // numerical data
            if (layer['normalization'] == 'sqrt') 
            {
                layerdata_dict[id][pindex] = Math.sqrt(parseFloat(layerdata_dict[id][pindex]));
            }
            if (layer['normalization'] == 'log') 
            {
                layerdata_dict[id][pindex] = log10(parseFloat(layerdata_dict[id][pindex]) + 1);
            }
            if (typeof param_max[pindex] === 'undefined' || parseFloat(layerdata_dict[id][pindex]) > parseFloat(param_max[pindex])) 
            {
                param_max[pindex] = parseFloat(layerdata_dict[id][pindex]);
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

            for (var id in layerdata_dict)
            {
                if (layer_types[pindex] == 1) // stack bar
                {
                    var total = 0;

                    for (var j=0; j < layerdata_dict[id][pindex].length; j++)
                    {
                        total = total + parseFloat(layerdata_dict[id][pindex][j]);
                    }

                    var multiplier = parseFloat(layers[pindex]['height']) / total;

                    for (var j=0; j < layerdata_dict[id][pindex].length; j++)
                    {
                        layerdata_dict[id][pindex][j] = layerdata_dict[id][pindex][j] * multiplier;
                    }
                }
                else // numerical data
                {
                    var bar_size = parseFloat(layerdata_dict[id][pindex]);
                                        
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
                            layerdata_dict[id][pindex] = 0;
                        } else {
                            layerdata_dict[id][pindex] = bar_size *  parseFloat(layers[pindex]['height']) / (max - min);
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
                            layerdata_dict[id][pindex] = 0;
                        } else {
                            layerdata_dict[id][pindex] = bar_size *  parseFloat(layers[pindex]['height']) / param_max[pindex];
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
    createBin('svg', 'viewport');
    createBin('viewport', 'tree_bin');
    createBin('tree_bin', 'tree');
    drawLine('tree', {'id': '_origin'}, {'x': 0, 'y': 0}, {'x': 0, 'y': 0}, false);

    if (settings['tree-type'] == 'phylogram')
        $('#tree_bin').attr('transform', 'rotate(90)'); 

    if (hasTree) 
    {
        t.ComputeWeights(t.root);
        var td = null;

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

                // calculate height per leaf
                height_per_leaf = width / (t.num_leaves - 1);
                break;

            case 'circlephylogram':
                if (t.has_edge_lengths) {
                    td = new CirclePhylogramDrawer();
                } else {

                    td = new CircleTreeDrawer();
                }
                td.Init(t, {
                    svg_id: 'tree',
                    width: (radius == 0) ? VIEWER_WIDTH : radius,
                    height: (radius == 0) ? VIEWER_HEIGHT : radius,
                    fontHeight: 10,
                    root_length: 0.1
                });

                // calculate angle per leaf
                var angle_min = parseFloat(settings['angle-min']);
                var angle_max = parseFloat(settings['angle-max']);
                angle_per_leaf = Math.toRadians(angle_max - angle_min) / t.num_leaves;
                break;
        }

        td.CalcCoordinates();
        td.Draw();
    }

    // calculate max radius of tree
    layer_boundaries = new Array();
    var tree_radius = 0;
    var tree_max_x = 0;
    var tree_max_y = 0;

    if (hasTree)
    {
        var n = new NodeIterator(t.root);
        var q = n.Begin();

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
                    q = n.Next();

                }
                layer_boundaries.push( [0, tree_max_x] );

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
                    q = n.Next();
                }
                layer_boundaries.push( [0, tree_radius] );
                break;
        }
    }
    else
    {
        if (settings['tree-type'] == 'circlephylogram') {
            layer_boundaries.push([0, radius]);
        }
        else
        {
            tree_max_y = width;
            tree_max_x = 0;
            layer_boundaries.push([0, 0]);
        }
    }

    margin = parseFloat(settings['layer-margin']);
    var layer_fonts = {};

    // calculate layer boundries
    for (var i = 0; i < settings['layer-order'].length; i++) {
        var layer_index = i+1;
        var pindex = settings['layer-order'][i];
        var layer = settings['views'][current_view][pindex];

        if (settings['custom-layer-margin'])
        {
            var layer_margin = parseFloat(layers[pindex]['margin']);
        }
        else
        {
            var layer_margin = margin;
        }

        if (!(layer_types[pindex] == 2 && layers[pindex]['type'] == 'text') && parseFloat(layers[pindex]['height'])==0)
        {
            layer_margin = 0;
        }

        // calculate per layer font size
        if (settings['tree-type'] == 'circlephylogram')
        {
            var layer_perimeter = ((angle_max - angle_min) / 360) * (2 * Math.PI * (layer_boundaries[i][1] + layer_margin));
        }
        else
        {
            var layer_perimeter = tree_max_y;
        }

        var layer_font = Math.min((layer_perimeter / leaf_count), parseFloat(settings['max-font-size']));
        layer_fonts[layer_index] = layer_font;

        // calculate new layer height if text layer heigth is 0
        if (layer_types[pindex] == 2 && layers[pindex]['type'] == 'text' && parseFloat(layers[pindex]['height'])==0)
        {
            // find longest item.
            var longest_text_len = 0;
            for (var _pos = 1; _pos < layerdata.length; _pos++)
            {
                if (layerdata[_pos][pindex] != null && layerdata[_pos][pindex].length > longest_text_len)
                {
                    longest_text_len = layerdata[_pos][pindex].length;
                }
            }
            // make background bit longer than text
            longest_text_len += 2;

            layers[pindex]['height'] = Math.ceil(longest_text_len * MONOSPACE_FONT_ASPECT_RATIO * layer_font) + 1;
            $('#height' + pindex).val(layers[pindex]['height']);
        }

        layer_boundaries.push( [ layer_boundaries[i][1] + layer_margin, layer_boundaries[i][1] + layer_margin + parseFloat(layers[pindex]['height']) ] );

        createBin('tree_bin', 'layer_background_' + layer_index);
        createBin('tree_bin', 'layer_' + layer_index);
        createBin('tree_bin', 'event_catcher_' + layer_index);

        // draw event catcher of the layer
        if (settings['tree-type']=='phylogram')
        {
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
            var _min = Math.toRadians(settings['angle-min']);
            var _max = Math.toRadians(settings['angle-max']);

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

        if (layer_types[pindex] == 2 && layers[pindex]['type'] == 'text') {
            var _bgcolor = layers[pindex]['color-start'];
            var _opacity = 1;
        } else {
            var _bgcolor = layers[pindex]['color'];
            var _opacity = 0.3;
        }

        // draw backgrounds
        if (settings['tree-type']=='phylogram' && ((layer_types[pindex] == 3 && layers[pindex]['type'] == 'bar') || (layer_types[pindex] == 2 && layers[pindex]['type'] == 'text')))
        {
            // if phylogram and
            // (numerical (3) and bar type) or (categorical (2) and text type)
            drawPhylogramRectangle('layer_background_' + layer_index,
                'all',
                layer_boundaries[layer_index][0],
                tree_max_y / 2,
                tree_max_y + height_per_leaf,
                layer_boundaries[layer_index][1] - layer_boundaries[layer_index][0],
                _bgcolor,
                _opacity,
                false);
        }

        if (settings['tree-type']=='circlephylogram' && ((layer_types[pindex] == 3 && layers[pindex]['type'] == 'bar') || (layer_types[pindex] == 2 && layers[pindex]['type'] == 'text')))
        {
            // if circlephylogram and
            // (numerical (3) and bar type) or (categorical (2) and text type)
            var _min = Math.toRadians(settings['angle-min']);
            var _max = Math.toRadians(settings['angle-max']);

            drawPie('layer_background_' + layer_index,
                'all',
                _min,
                _max,
                layer_boundaries[layer_index][0],
                layer_boundaries[layer_index][1],
                (_max - _min > Math.PI) ? 1:0, // large arc flag
                _bgcolor,
                _opacity,
                false);
        }
    }

    total_radius = layer_boundaries[layer_boundaries.length - 1][1];
    beginning_of_layers = layer_boundaries[0][1];

    var n = new NodeIterator(t.root);
    var q = n.Begin();

    createBin('tree_bin', 'guide_lines');

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

    var numeric_cache = {};

    switch (settings['tree-type']) {
        case 'phylogram':
            while (q != null) {
                if (q.IsLeaf()) {
                    odd_even_flag = odd_even_flag * -1;

                    if (odd_even_flag > 0 && hasTree)
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
                            for (var j=0; j < layerdata_dict[q.label][pindex].length; j++)
                            {
                                drawPhylogramRectangle('layer_' + layer_index,
                                    q.id,
                                    layer_boundaries[layer_index][1] - offset - layerdata_dict[q.label][pindex][j],
                                    q.xy['y'],
                                    height_per_leaf,
                                    layerdata_dict[q.label][pindex][j],
                                    stack_bar_colors[pindex][j],
                                    1,
                                    false);
                                offset += layerdata_dict[q.label][pindex][j];
                            } 

                        }
                        else if(isCategorical)
                        {
                            if (typeof categorical_data_colors[pindex][layerdata_dict[q.label][pindex]] === 'undefined'){
                                categorical_data_colors[pindex][layerdata_dict[q.label][pindex]] = getNamedCategoryColor(layerdata_dict[q.label][pindex]);
                            }

                            if (layers[pindex]['type'] == 'color') 
                            {
                                categorical_layers_ordered[layer_index].push(layerdata_dict[q.label][pindex]);
                            }
                            else // text
                            {
                                var _offsetx = layer_boundaries[layer_index][0] + layer_fonts[layer_index] * MONOSPACE_FONT_ASPECT_RATIO;
                                var offset_xy = [];
                                offset_xy['x'] = _offsetx;
                                offset_xy['y'] = q.xy['y'];
                                var _label = (layerdata_dict[q.label][pindex] == null) ? '' : layerdata_dict[q.label][pindex];

                                drawRotatedText('layer_' + layer_index, offset_xy, _label, 0, 'left', layer_fonts[layer_index], "monospace", layers[pindex]['color'], layers[pindex]['height']);
                                
                            }
                        }
                        else if (isParent)
                        {
                            categorical_layers_ordered[layer_index].push(layerdata_dict[q.label][pindex]);
                        }
                        else // numerical
                        {
                            if (layers[pindex]['type'] == 'intensity')
                            {
                                var color = layers[pindex]['color'];

                                 drawPhylogramRectangle('layer_' + layer_index,
                                    q.id,
                                    layer_boundaries[layer_index][0] ,
                                    q.xy['y'],
                                    height_per_leaf,
                                    layer_boundaries[layer_index][1] - layer_boundaries[layer_index][0],
                                    getGradientColor(layers[pindex]['color-start'], color,  layerdata_dict[q.label][pindex] / layers[pindex]['height']),
                                    1,
                                    false);
                            }
                            else
                            {
                                if (settings['optimize-speed'])
                                {
                                    if (!numeric_cache.hasOwnProperty(layer_index)){
                                        numeric_cache[layer_index] = [];
                                    }

                                    if (numeric_cache[layer_index].length == 0)
                                    {
                                        numeric_cache[layer_index].push("M", layer_boundaries[layer_index][1], q.xy['y'] - height_per_leaf / 2);
                                    }

                                    numeric_cache[layer_index].push("L", layer_boundaries[layer_index][1] - layerdata_dict[q.label][pindex], q.xy['y'] - height_per_leaf / 2, 
                                        "L", layer_boundaries[layer_index][1] - layerdata_dict[q.label][pindex], q.xy['y'] + height_per_leaf / 2);

                                    if (q.order == leaf_count-1) {
                                        numeric_cache[layer_index].push("L", layer_boundaries[layer_index][1], q.xy['y'] + height_per_leaf / 2,
                                            "Z");
                                    }
                                }
                                else
                                {
                                    // x y height width
                                    if (layerdata_dict[q.label][pindex] > 0) {
                                        var color = layers[pindex]['color'];

                                         drawPhylogramRectangle('layer_' + layer_index,
                                            q.id,
                                            layer_boundaries[layer_index][1] - layerdata_dict[q.label][pindex],
                                            q.xy['y'],
                                            height_per_leaf,
                                            layerdata_dict[q.label][pindex],
                                            color,
                                            1,
                                            false);
                                    }
                                }
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
                    if (odd_even_flag > 0 && hasTree)
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
                            for (var j=0; j < layerdata_dict[q.label][pindex].length; j++)
                            {
                                drawPie('layer_' + layer_index,
                                    q.id,
                                    q.angle - angle_per_leaf / 2,
                                    q.angle + angle_per_leaf / 2,
                                    layer_boundaries[layer_index][0] + offset,
                                    layer_boundaries[layer_index][0] + offset + layerdata_dict[q.label][pindex][j],
                                    0,
                                    stack_bar_colors[pindex][j],
                                    1,
                                    false);
                                offset += layerdata_dict[q.label][pindex][j];
                            } 
                    
                        }
                        else if(isCategorical)
                        {

                            if (typeof categorical_data_colors[pindex][layerdata_dict[q.label][pindex]] === 'undefined'){
                                categorical_data_colors[pindex][layerdata_dict[q.label][pindex]] = getNamedCategoryColor(layerdata_dict[q.label][pindex]);
                            }

                            if (layers[pindex]['type'] == 'color') 
                            {
                                categorical_layers_ordered[layer_index].push(layerdata_dict[q.label][pindex]);
                            }
                            else // text
                            {
                                var align = 'left';
                                var new_angle = q.angle * 180.0 / Math.PI;
                                if ((q.angle > Math.PI / 2.0) && (q.angle < 1.5 * Math.PI)) {
                                    align = 'right';
                                    new_angle += 180.0;
                                }

                                var offset_xy = [];
                                var _radius = layer_boundaries[layer_index][0] + layer_fonts[layer_index] * MONOSPACE_FONT_ASPECT_RATIO;
                                offset_xy['x'] = Math.cos(q.angle) * _radius;
                                offset_xy['y'] = Math.sin(q.angle) * _radius;
                                var _label = (layerdata_dict[q.label][pindex] == null) ? '' : layerdata_dict[q.label][pindex];

                                drawRotatedText('layer_' + layer_index, offset_xy, _label, new_angle, align, layer_fonts[layer_index], "monospace", layers[pindex]['color'], layers[pindex]['height']);
                            }
                        }
                        else if (isParent)
                        {
                            categorical_layers_ordered[layer_index].push(layerdata_dict[q.label][pindex]);
                        }
                        else // numerical
                        {
                            var color = layers[pindex]['color'];

                            if (layers[pindex]['type'] == 'intensity')
                            {
                                    drawPie('layer_' + layer_index,
                                        q.id,
                                        q.angle - angle_per_leaf / 2,
                                        q.angle + angle_per_leaf / 2,
                                        layer_boundaries[layer_index][0], 
                                        layer_boundaries[layer_index][1],
                                        0,
                                        getGradientColor(layers[pindex]['color-start'], color,   layerdata_dict[q.label][pindex] / layers[pindex]['height']),
                                        1,
                                        false);
                            }
                            else
                            {
                                if (settings['optimize-speed'])
                                {
                                    if (!numeric_cache.hasOwnProperty(layer_index)){
                                        numeric_cache[layer_index] = [];
                                    }

                                    start_angle = q.angle - angle_per_leaf / 2;
                                    end_angle = q.angle + angle_per_leaf / 2;
                                    inner_radius = layer_boundaries[layer_index][0];
                                    outer_radius = layer_boundaries[layer_index][0] + layerdata_dict[q.label][pindex];

                                    if (numeric_cache[layer_index].length == 0)
                                    {
                                        var ax = Math.cos(start_angle) * inner_radius;
                                        var ay = Math.sin(start_angle) * inner_radius;

                                        numeric_cache[layer_index].push("M", ax, ay);
                                    }

                                    var cx = Math.cos(end_angle) * outer_radius;
                                    var cy = Math.sin(end_angle) * outer_radius;

                                    var dx = Math.cos(start_angle) * outer_radius;
                                    var dy = Math.sin(start_angle) * outer_radius;

                                    numeric_cache[layer_index].push("L", dx, dy, "A", outer_radius, outer_radius, 0, 0, 1, cx, cy);

                                    if (q.order == leaf_count-1) {
                                        var bx = Math.cos(end_angle) * inner_radius;
                                        var by = Math.sin(end_angle) * inner_radius;

                                        var _min = Math.toRadians(settings['angle-min']);
                                        var _max = Math.toRadians(settings['angle-max']);
                                        var large_arc_flag = (_max - _min > Math.PI) ? 1:0;

                                        numeric_cache[layer_index].push("L", bx, by, 
                                            "A", inner_radius, inner_radius, 0, large_arc_flag, 0, numeric_cache[layer_index][1], numeric_cache[layer_index][2], 
                                            "Z");
                                    }
                                }
                                else
                                {
                                    if (layerdata_dict[q.label][pindex] > 0) {
                                        drawPie('layer_' + layer_index,
                                            q.id,
                                            q.angle - angle_per_leaf / 2,
                                            q.angle + angle_per_leaf / 2,
                                            layer_boundaries[layer_index][0], 
                                            layer_boundaries[layer_index][0] + layerdata_dict[q.label][pindex],
                                            0,
                                            color,
                                            1,
                                            false);
                                    }
                                }
                            }
                        }

                    }
                }
                q = n.Next();
            }
            break;
    }

    //console.log(numeric_cache);

    // cluster categorical items and draw them
    
    for (var i = 0; i < settings['layer-order'].length; i++) {
        var layer_index = i+1;
        var pindex = settings['layer-order'][i];

        var isParent      = (layer_types[pindex] == 0) ? true : false;
        var isCategorical = (layer_types[pindex] == 2) ? true : false;
        var isNumerical = (layer_types[pindex] == 3) ? true : false;


        if (isParent || isCategorical)
        {
            var layer_items = categorical_layers_ordered[layer_index];

            // This algorithm finds changing points in the categorical layer
            // to draw sequential items with single object, so we have to add
            // a dummy item to the end in order to detect last group.
            // I have used 'null', but if the last group is already 'null', algorithm skips
            // that group since there is no change, so we have to check the last item
            // and if it is 'null' we have to add something different.

            if (layer_items[layer_items.length-1] == null) {
                layer_items.push(-1);
            } else {
                layer_items.push(null)
            }

            var prev_value = layer_items[0];
            var prev_start = 0;

            var items_to_draw = new Array();

            for (var j=1; j < layer_items.length; j++)
            {
                if (prev_value != layer_items[j])
                {
                    if (isCategorical || prev_value != '')
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
        else if (isNumerical && (settings['optimize-speed']) && layers[pindex]['type'] != 'intensity') {
                var path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
                path.setAttribute('stroke-width', '0');
                path.setAttribute('shape-rendering', 'auto');
                path.setAttribute('pointer-events', 'none');
                path.setAttribute('fill', layers[pindex]['color']);
                path.setAttribute('d', numeric_cache[layer_index].join(' '));
                var layer_group = document.getElementById('layer_' + layer_index);
                layer_group.appendChild(path);
        }
    }

    rebuildIntersections();
    createBin('tree_bin', 'layer_labels');
    createBin('tree_bin', 'bin');
    // redrawBins() call moved to draw_tree_callback because redrawBins needs last_settings to be set.

    // observe for transform matrix change
    var observer = new MutationObserver(updateSingleBackgroundGlobals);

    observer.observe(document.getElementById('viewport'), {
        attributes:    true,
        attributeFilter: ["transform"]
    });

    // draw layer names and samples in circlephylogram for special angles.
    if (settings['tree-type'] == 'circlephylogram')
    {
        if (settings['angle-max'] == 270)
        {
            for (var i = 0; i < settings['layer-order'].length; i++) {
                var layer_index = i+1;
                var pindex = settings['layer-order'][i];
                var layer = settings['views'][current_view][pindex];

                // !!IMPORTANT!! no label for hidden layers, font-size: 0px causes bug on inkscape 0.91
                if ((layers[pindex]['height']) == 0)
                    continue;

                var layer_title = layerdata[0][pindex];

                if (layer_title in named_layers && 'pretty_name' in named_layers[layer_title]) {
                    layer_title = named_layers[layer_title]['pretty_name'];
                } else if(layer_title.substring(0, 5) == "hmmx_") {
                    layer_title = layer_title.replace(/hmmx_/g, "").replace(/_/g, " ");
                } else if(layer_title.substring(0, 5) == "hmms_") {
                    layer_title = layer_title.replace(/hmms_/g, "").replace(/_/g, " ");
                } else {
                    layer_title = layer_title.replace(/_/g, " ");
                }

                if (layer_title.indexOf('!') > -1 )
                {
                    layer_title = layer_title.split('!')[0];
                }

                drawFixedWidthText('layer_labels', {
                        'x': 10,
                        'y': 0 - layer_boundaries[layer_index][1] + (layers[pindex]['height'] * 3/4) 
                    }, 
                    layer_title, 
                    layers[pindex]['height'] + 'px',
                    layers[pindex]['color'],
                    total_radius,
                    layers[pindex]['height']);
            }
        }

        if (settings['angle-min'] == 0 && settings['angle-max'] <= 270)
        {
            //draw samples layers
            createBin('viewport', 'samples');
            drawSamples();
        }
    }

    // draw title
    var legend_top = total_radius + parseFloat(settings['layer-margin']) + parseFloat(settings['outer-ring-height']) * 2;
    switch (settings['tree-type']) {
        case 'phylogram':
            drawLegend(legend_top, 0 - VIEWER_HEIGHT, 0);
            drawTitle(-150, -0.5 * VIEWER_HEIGHT, settings);
            break;
        case 'circlephylogram':
            drawLegend(legend_top, 0 - total_radius, total_radius - 40);
            drawTitle(-1 * total_radius - 150, 0, settings);
            break;
    }

    // Scale to fit window
    bbox = svg.getBBox();
    scale = Math.min(VIEWER_WIDTH / bbox.width, VIEWER_HEIGHT / bbox.height) * 0.80;

    zoom_reset();

    // pan and mouse zoom
    $('svg').svgPan('viewport');

    // tree bin = tree + layers
    // bind events
    var tree_bin = document.getElementById('tree_bin');
    tree_bin.addEventListener('click', lineClickHandler, false);
    tree_bin.addEventListener('contextmenu', lineContextMenuHandler, false);
    tree_bin.addEventListener('mouseover',lineMouseEnterHandler, false);
    tree_bin.addEventListener('mouseout', lineMouseLeaveHandler, false);
    document.body.addEventListener('mousemove', mouseMoveHandler, false); // for tooltip
    document.body.addEventListener('click', function() { $('#default_right_click_menu').hide(); }, false);
    document.body.addEventListener('click', function() { $('#collection_mode_right_click_menu').hide(); }, false);
    document.body.addEventListener('click', function() { $('#pan_mode_right_click_menu').hide(); }, false);

    // code below required to stop clicking on contings while panning.
    var viewport = document.getElementById('svg');
    viewport.addEventListener('mousedown', 
        function(event) { 
            dragging = false; 

            if (event.shiftKey)
            {
                drawing_zoom = true;

                zoomBox['start_x'] = event.clientX;
                zoomBox['start_y'] = event.clientY;

                $('#divzoom').css({"top": 0, "left": 0, "width": 0, "height": 0 });
                $('#divzoom_inner').css({"width": 0, "height": 0 });
                $('#divzoom').show();
            }
        });
    viewport.addEventListener('mousemove', 
        function(event) { 
            dragging = true; 

            if (event.shiftKey && drawing_zoom)
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
                
                if (inner_rect.width > 2 && inner_rect.height > 2)
                {
                    var _dx = (parseInt("0" + $('#svg').position().left) + (VIEWER_WIDTH / 2)) - (inner_rect.left + inner_rect.width / 2);
                    var _dy = (parseInt("0" + $('#svg').position().top)  + (VIEWER_HEIGHT / 2)) - (inner_rect.top + inner_rect.height / 2);
                    pan(_dx,_dy);
                    zoom(VIEWER_WIDTH / inner_rect.width);
                }
            }

            clearTextSelection();
            drawing_zoom=false; 
            zoomBox = {}; 
            $('#divzoom').hide(); 

        });

    var tree_object_count = document.getElementById('tree').getElementsByTagName('*').length + document.getElementById('guide_lines').getElementsByTagName('*').length;
    var total_object_count = document.getElementById('svg').getElementsByTagName('*').length;

    $('#draw_delta_time').html(leaf_count + ' splits and ' + total_object_count +' objects drawn in ' + tree_draw_timer.getDeltaSeconds('done')['deltaSecondsStart'] + ' seconds.');

    console.log('[info] Leaf count: ' + leaf_count);
    console.log('[info] Object count in tree (with guide lines): ' + tree_object_count);
    console.log('[info] Total objects in SVG: ' + total_object_count);

}

function redrawBins()
{
    // check if tree parsed, if not there is nothing to redraw.
    if ($.isEmptyObject(label_to_node_map)) 
        return;

    var leaf_list = Array.apply(null, new Array(leaf_count+1)).map(Number.prototype.valueOf,0);

    // put bin numbers of selected leaves to leaf list
    // maybe we should write directly into leaf_list in mouse events, instead of generate it everytime.
    for (var bin_id = 1; bin_id <= bin_counter; bin_id++) {
        for (var j = 0; j < SELECTED[bin_id].length; j++) {
            if (label_to_node_map[SELECTED[bin_id][j]].IsLeaf()) {
                leaf_list[label_to_node_map[SELECTED[bin_id][j]].order] = bin_id;
            }
        }
    }

    // cluster bins and put them into bins_to_draw array with (start, end, bin_id);
    var prev_value = leaf_list[0];
    var prev_start = 0;

    var bins_to_draw = new Array();

    for (var i=1; i < leaf_list.length; i++)
    {
        if (prev_value != leaf_list[i])
        {
            if (prev_value != 0) {
                bins_to_draw.push(new Array(prev_start, i - 1, prev_value)); // start, end, bin_id;
            }

            prev_start = i;
        }
        prev_value = leaf_list[i];
    }

    // remove exist bin drawings
    var bin = document.getElementById('bin');
    while (bin.hasChildNodes()) {
        bin.removeChild(bin.lastChild);
    }

    // draw new bins
    var show_grid = $('#show_grid_for_bins')[0].checked;
    var grid_color = document.getElementById('grid_color').getAttribute('color');
    var grid_width = $('#grid_width').val();
    var show_bin_labels = $('#show_bin_labels')[0].checked;
    var bin_labels_font_size = $('#bin_labels_font_size').val();
    var autorotate_bin_labels = $('#autorotate_bin_labels')[0].checked;
    var bin_labels_angle = $('#bin_labels_angle').val();
    
    var outer_ring_size = parseFloat($('#outer-ring-height').val());
    var outer_ring_margin = parseFloat($('#outer-ring-margin').val());

    for (var i=0; i < bins_to_draw.length; i++) {
        var start = order_to_node_map[bins_to_draw[i][0]];
        var end = order_to_node_map[bins_to_draw[i][1]];

        var color = document.getElementById('bin_color_' + bins_to_draw[i][2]).getAttribute('color');

        if (tree_type == 'circlephylogram')
        {

            drawPie('bin',
                'bin_outer_' + i,
                start.angle - angle_per_leaf / 2,
                end.angle + angle_per_leaf / 2,
                total_radius + outer_ring_margin,
                total_radius + outer_ring_margin + outer_ring_size,
                (end.angle - start.angle + angle_per_leaf > Math.PI) ? 1 : 0,
                color,
                1,
                false);

            var align = 'left';
            var angle = (end.angle + start.angle) / 2;
            var new_angle = angle * 180.0 / Math.PI;
            if ((angle > Math.PI / 2.0) && (angle < 1.5 * Math.PI)) {
                align = 'right';
                new_angle += 180.0;
            }

            if (show_bin_labels)
            {
                drawRotatedText(
                    'bin',
                    {
                        'x': (total_radius + outer_ring_margin * 1.5 + outer_ring_size * (highlighted_splits.length > 0 ? 2 : 1)) * Math.cos((end.angle + start.angle) / 2), 
                        'y': (total_radius + outer_ring_margin * 1.5 + outer_ring_size * (highlighted_splits.length > 0 ? 2 : 1)) * Math.sin((end.angle + start.angle) / 2), 
                    },
                    $('#bin_name_' + bins_to_draw[i][2]).val().replace("_", " "),
                    (autorotate_bin_labels) ? new_angle : bin_labels_angle,
                    align,
                    bin_labels_font_size + "px",
                    "HelveticaNeue-CondensedBold, Helvetica Neue, Helvetica, sans-serif",
                    color,
                    0
                    );

            }

            var pie = drawPie('bin',
                'bin_background_' + i,
                start.angle - angle_per_leaf / 2,
                end.angle + angle_per_leaf / 2,
                beginning_of_layers,
                (show_grid) ? total_radius + outer_ring_margin + outer_ring_size : total_radius,
                (end.angle - start.angle + angle_per_leaf > Math.PI) ? 1 : 0,
                color,
                (show_grid) ? 0 : 0.1,
                false);

            if (show_grid) {
                pie.setAttribute('vector-effect', 'non-scaling-stroke');
                pie.setAttribute('stroke-opacity', '1');
                pie.setAttribute('stroke-width', grid_width);
                pie.setAttribute('stroke', grid_color);
            }


        }
        else
        {

            drawPhylogramRectangle('bin',
                'bin_outer_' + i,
                total_radius + outer_ring_margin,
                (start.xy.y + end.xy.y) / 2,
                end.xy.y - start.xy.y + height_per_leaf,
                outer_ring_size,
                color,
                1,
                false);

            if (show_bin_labels)
            {
                drawRotatedText(
                    'bin',
                    {
                        'y': (start.xy.y + end.xy.y) / 2, 
                        'x': (total_radius + outer_ring_margin * 1.5 + outer_ring_size * (highlighted_splits.length > 0 ? 2 : 1)), 
                    },
                    $('#bin_name_' + bins_to_draw[i][2]).val().replace("_", " "),
                    (autorotate_bin_labels) ? 0 : bin_labels_angle,
                    'left',
                    bin_labels_font_size + "px",
                    "HelveticaNeue-CondensedBold, Helvetica Neue, Helvetica, sans-serif",
                    color,
                    0
                    );

            }

            var rect = drawPhylogramRectangle('bin',
                'bin_background_' + i,
                beginning_of_layers,
                (start.xy.y + end.xy.y) / 2,
                end.xy.y - start.xy.y + height_per_leaf,
                (show_grid) ? total_radius + outer_ring_margin + outer_ring_size - beginning_of_layers : total_radius + margin - beginning_of_layers,
                color,
                (show_grid) ? 0 : 0.1,
                false);

            if (show_grid) {
                rect.setAttribute('vector-effect', 'non-scaling-stroke');
                rect.setAttribute('stroke-opacity', '1');
                rect.setAttribute('stroke-width', grid_width);
                rect.setAttribute('stroke', grid_color);
            }
        }
    }


    // draw higlighted splits
    for (var i=0; i < highlighted_splits.length; i++) {
        // TO DO: more performance
        var start = label_to_node_map[highlighted_splits[i]];
        var end = start;

        var color = document.getElementById('picker_highlight').getAttribute('color');
        var outer_ring_size = parseInt(last_settings['outer-ring-height']);

        if (tree_type == 'circlephylogram')
        {
            drawPie('bin',
                'bin_outer_' + 1,
                start.angle - angle_per_leaf / 2,
                end.angle + angle_per_leaf / 2,
                total_radius + outer_ring_margin + outer_ring_size,
                total_radius + outer_ring_margin + outer_ring_size * 2,
                (end.angle - start.angle + angle_per_leaf > Math.PI) ? 1 : 0,
                color,
                1,
                false);     
        }
        else
        {
            drawPhylogramRectangle('bin',
                'bin_outer_' + 1,
                total_radius + outer_ring_margin + outer_ring_size,
                (start.xy.y + end.xy.y) / 2,
                end.xy.y - start.xy.y + height_per_leaf,
                outer_ring_size,
                color,
                1,
                false);
        }
    }
}


function rebuildIntersections()
{
    if (!hasTree)
        return;

    for (var bin_id = 1; bin_id <= bin_counter; bin_id++) {

        // delete extra intersections
        var deleted;
        do {
            deleted = 0;
            var cursor = SELECTED[bin_id].length;
            while (cursor--)
            {
                var node = label_to_node_map[SELECTED[bin_id][cursor]];
                
                if (node.IsLeaf())
                    continue;

                if (node.child != null && SELECTED[bin_id].indexOf(node.child.label) > -1)
                    continue;

                if (node.sibling != null && SELECTED[bin_id].indexOf(node.sibling.label) > -1)
                    continue;

                SELECTED[bin_id].splice(cursor,1);
                deleted++;
            }

        } while (deleted > 0)

        // try to make new intersections
        var inserted;
        do {
            inserted = 0;
            var length = SELECTED[bin_id].length;
            for (var cursor = 0; cursor < length; cursor++)
            {
                var node = label_to_node_map[SELECTED[bin_id][cursor]];
                var parent = node.ancestor;

                if (parent == null || parent.ancestor == null) 
                {
                    // skip root
                    continue;
                }

                if (SELECTED[bin_id].indexOf(parent.label) > -1)
                {
                    // parent already in selected list
                    continue;
                }

                if (node.sibling != null && SELECTED[bin_id].indexOf(node.sibling.label) > -1)
                {
                    SELECTED[bin_id].push(parent.label);
                    
                    inserted++;
                }
            }

        } while (inserted > 0)
    }
}
