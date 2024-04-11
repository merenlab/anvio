/**
 *  Helper functions to draw SVG objects.
 *
 *  Authors: Ozcan Esen
 *           Matthew Klein <mtt.l.kln@gmail.com>
 *           A. Murat Eren <a.murat.eren@gmail.com>
 *
 * Copyright 2015-2021, The anvi'o project (http://anvio.org)

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

function drawTitle(settings) {
    createBin('viewport', 'title_group');
    _bbox = document.getElementById('viewport').getBBox();

    var top = 0;
    var left = 0;
    var _font_size = 0;
    switch (settings['tree-type']) {
        case 'phylogram':
            top = _bbox.y - 160;
            left = _bbox.x + _bbox.width / 2;
            _font_size = _bbox.width / 40;
            break;
        case 'circlephylogram':
            top = _bbox.y - 160;
            left = 0;
            _font_size = _bbox.width / 80;
            break;
    }

    top -= 3 * _font_size;
    drawText('title_group', {'x': left, 'y': top}, document.title, 2 * _font_size + 'px', 'center');
    top += 2 * _font_size;
    var _sub_title = "Items order: " + getClusteringPrettyName(settings['order-by']) + " | ";
    _sub_title    += "Current view: " + settings['current-view'] + " | ";
    _sub_title    += "Samples order: " + settings['samples-order'];

    drawText('title_group', {'x': left, 'y': top}, _sub_title, _font_size + 'px', 'center');
}

//--------------------------------------------------------------------------------------------------
function drawLegend() {
    createBin('viewport', 'legend_group');
    _bbox = document.getElementById('viewport').getBBox();

    var top = _bbox.y + _bbox.height;
    var left = _bbox.x;
    var line_end = _bbox.x + _bbox.width;

    var _left = left;
    var line_height = (line_end - _left) / 80;
    var gap = line_height / 2;

    var top = top + line_height * 3;
    var _top = top;

    for (var i=0; i < legends.length; i++)
    {
        var bin_id = 'legend_' + i;
        createBin('legend_group', bin_id);
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

            let color = '#FFFFFF';
            if (typeof legend['group'] === 'undefined') {
                color = window[legend['source']][legend['key']][legend['item_keys'][j]];
            } else {
                color = window[legend['source']][legend['group']][legend['key']][legend['item_keys'][j]];
            }

            var rect = drawRectangle(bin_id, _left, _top, line_height, line_height, color, 1, 'black',
                null,
                function() {
                    // mouseenter
                    $(this).css('stroke-width', '2');
                },
                function() {
                    // mouseleave
                    $(this).css('stroke-width', '1');
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

    if (mode == 'pan') {
        drawText('bin_legend', {'x': left + 170, 'y': top}, 'PCs', '10px');
        drawText('bin_legend', {'x': left + 230, 'y': top}, 'Gene Calls', '10px');
    } else {
        drawText('bin_legend', {'x': left + 170, 'y': top}, 'Contigs', '10px');
        drawText('bin_legend', {'x': left + 230, 'y': top}, 'Length', '10px');
    }

    for (var bin_id=0; bin_id < bins_to_draw.length; bin_id++) {
        var bin = bins_to_draw[bin_id];
        top = top + 20;

        drawRectangle('bin_legend', left, top-8, 16, 16, bin['color'], 1, 'black');
        drawText('bin_legend', {'x': left + 30, 'y': top }, bin['name'], '12px');

        if (mode == 'pan') {
            drawText('bin_legend', {'x': left + 170, 'y': top}, bin['pcs'], '12px');
            drawText('bin_legend', {'x': left + 230, 'y': top}, bin['gene-calls'], '12px');
        } else {
            drawText('bin_legend', {'x': left + 170, 'y': top}, bin['contig-count'], '12px');
            drawText('bin_legend', {'x': left + 230, 'y': top}, bin['contig-length'], '12px');
        }
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

function drawSupportValue(svg_id, p, p0, p1, supportValueData) {
    function checkInRange(){
        /**
         * Check if the branch support values fall within the given number range.
         * @return {bool}    True if the branch support values are within the specified range, False otherwise.
        */
        if (typeof p.branch_support === 'string' && p.branch_support.includes('/')) {
            const [branch_support_value0, branch_support_value1] = p.branch_support.split('/').map(parseFloat);
            const min_support = Math.min(branch_support_value0, branch_support_value1);
            const max_support = Math.max(branch_support_value0, branch_support_value1);
            return min_support >= supportValueData.numberRange[0] && max_support <= supportValueData.numberRange[1];
        } else {
            return p.branch_support >= supportValueData.numberRange[0] && p.branch_support <= supportValueData.numberRange[1];
        }
    }

    if( supportValueData.showNumber && checkInRange()){ // only render text if in range AND selected by user
        if($('#tree_type').val() == 'circlephylogram'){
            if(supportValueData.textRotation == '0'){
                drawText(svg_id, p.xy, p.branch_support, supportValueData.fontSize, 'Roboto' ,'black', '' , 'baseline')
            } else {
                drawRotatedText(svg_id, p.xy, p.branch_support, parseInt(supportValueData.textRotation), 'right', supportValueData.fontSize, 'Roboto' ,'black', '' , 'baseline')
            }
        } else {
            if(supportValueData.textRotation == '0'){
                drawRotatedText(svg_id, p.xy, p.branch_support, -90, 'left', supportValueData.fontSize, 'Roboto' ,'black', '' , 'baseline')
            } else {
                drawRotatedText(svg_id, p.xy, p.branch_support, parseInt(supportValueData.textRotation), 'left', supportValueData.fontSize, 'Roboto' ,'black', '' , 'baseline')
            }
        }
    }
    if(supportValueData.showSymbol && checkInRange()){ // only render symbol if in range AND selected by user
        drawSymbol()
    }

    /**
     * Draw symbols on an SVG canvas based on branch support values.
     * @returns {Element} The first circle element drawn on the SVG canvas.
     */
    function drawSymbol() {
        let first_circle = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
        let second_circle = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
        let maxRadius = supportValueData.maxRadius;
        let minRadius = supportValueData.minRadius;
        let rangeLow = parseFloat(supportValueData.numberRange[0]);
        let rangeHigh = parseFloat(supportValueData.numberRange[1]);

        /**
         * Calculate the percentile of a support value within the range.
         * @param {number} support_value - The branch support value.
         * @returns {number} The calculated percentile.
         */
        function calculatePercentile(support_value) {
            return (support_value - rangeLow) / (rangeHigh - rangeLow);
        }

        /**
         * Set details for the circle based on the percentile.
         * @param {Element} circle - The SVG circle element.
         * @param {number} percentile - The percentile value.
         */
        function setDetails(circle, percentile) {
            let radius;
            if (percentile > 0.67) {
                supportValueData.invertSymbol ? radius = maxRadius * 0.4 : radius = maxRadius;
            } else if (percentile < 0.67 && percentile > 0.33) {
                radius = maxRadius * 0.7;
                if(radius < minRadius && !supportValueData.invertSymbol){
                    radius = minRadius;
                }
            } else {
                supportValueData.invertSymbol ? radius = maxRadius : radius = maxRadius * 0.4;
                if(radius < minRadius && !supportValueData.invertSymbol){
                    radius = minRadius;
                }
            }
            circle.setAttribute('r', radius);
            circle.setAttribute('fill', supportValueData.symbolColor);
            circle.setAttribute('opacity', 0.6);
        }

        /**
         * Create a circle element based on the support value.
         * @param {number} support_value - The branch support value.
         * @param {number} index - The index of the circle.
         * @returns {Element} The created circle element.
         */
        function makeCircle(support_value, index) {
            let circle = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
            setDetails(circle, calculatePercentile(support_value));
            circle.setAttribute('cx', p.xy.x);
            circle.setAttribute('cy', p.xy.y);
            circle.setAttribute('id', p.id + '_' + index);
            var svg = document.getElementById(svg_id);
            svg.appendChild(circle);
            return circle;
        }

        if (typeof p.branch_support === 'string' && p.branch_support.includes('/')) {
            let branch_support_values = p.branch_support.split('/').map(parseFloat);
            first_circle = makeCircle(branch_support_values[0], 0);
            second_circle = makeCircle(branch_support_values[1], 1);
            second_circle.setAttribute('fill', supportValueData.secondSymbolColor);

            if($('#tree_type').val() == 'circlephylogram'){
                let length_original = Math.sqrt(p.xy.x**2 + p.xy.y**2);
                let length_scaled = length_original + 2 * maxRadius;
                let scaling_factor = length_scaled/length_original;
                second_circle.setAttribute('cx', p.xy.x * scaling_factor);
                second_circle.setAttribute('cy', p.xy.y * scaling_factor);
            } else{
                let firstCircleCx = parseFloat(first_circle.getAttribute('cx'));
                let updatedCx = firstCircleCx + maxRadius * 2;
                second_circle.setAttributeNS(null, 'cx', updatedCx.toString());
            }

        } else {
            first_circle = makeCircle(p.branch_support, 0);
        }

        return first_circle;
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
function drawText(svg_id, p, string, font_size, align, color, baseline, supportClickable) {

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
    text.setAttribute('font-family', 'Roboto','Helvetica', 'Arial;')

    text.setAttribute('onclick', function(event){
        console.log(event)
    })



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
                return textObj;
            }
        }
        textObj.textContent="..."; //can't place at all
    }

    return textObj
}

//--------------------------------------------------------------------------------------------------
function drawRotatedText(svg_id, p, string, angle, align, font_size, font_family, color, maxLength, baseline) {
    var text = document.createElementNS('http://www.w3.org/2000/svg', 'text');
    //newLine.setAttribute('id','node' + p.id);
    text.setAttribute('style', 'alignment-baseline: ' + baseline);
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

    // responsible for mapping the circular arc at each branch bifurcation
    // change fill attribute to 'red' to visualize arc path

    var arc = document.createElementNS('http://www.w3.org/2000/svg', 'path');
    arc.setAttribute('id', 'arc' + p.id);

    arc.setAttribute('vector-effect', 'non-scaling-stroke');
    arc.setAttribute('style', 'stroke:' + LINE_COLOR + ';stroke-width:1;');
    arc.setAttribute('fill', 'none');

    var path = circeArcPath(p0, p1, radius, large_arc_flag);
    arc.setAttribute('d', path)

    var svg = document.getElementById(svg_id);
    svg.appendChild(arc);

    return arc;
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
// https://stackoverflow.com/questions/10894377/dynamically-adding-a-svg-gradient
function createGradient(svg,id,stops){
  var svgNS = svg.namespaceURI;
  var grad  = document.createElementNS(svgNS,'linearGradient');
  grad.setAttribute('id',id);
  for (var i=0;i<stops.length;i++){
    var attrs = stops[i];
    var stop = document.createElementNS(svgNS,'stop');
    for (var attr in attrs){
      if (attrs.hasOwnProperty(attr)) stop.setAttribute(attr,attrs[attr]);
    }
    grad.appendChild(stop);
  }

  var defs = svg.querySelector('defs') ||
      svg.insertBefore( document.createElementNS(svgNS,'defs'), svg.firstChild);
  return defs.appendChild(grad);
}


function getGradientColor(start_color, end_color, percent) {
   // strip the leading # if it's there
   start_color = start_color.replace(/^\s*#|\s*$/g, '');
   end_color = end_color.replace(/^\s*#|\s*$/g, '');

   // convert 3 char codes --> 6, e.g. `E0F` --> `EE00FF`
   if(start_color.length == 3){
     start_color = start_color.replace(/(.)/g, '$1$1');
   }

   if(end_color.length == 3){
     end_color = end_color.replace(/(.)/g, '$1$1');
   }

   // get colors
   var start_red = parseInt(start_color.substr(0, 2), 16),
       start_green = parseInt(start_color.substr(2, 2), 16),
       start_blue = parseInt(start_color.substr(4, 2), 16);

   var end_red = parseInt(end_color.substr(0, 2), 16),
       end_green = parseInt(end_color.substr(2, 2), 16),
       end_blue = parseInt(end_color.substr(4, 2), 16);

   // calculate new color
   var diff_red = end_red - start_red;
   var diff_green = end_green - start_green;
   var diff_blue = end_blue - start_blue;

   diff_red = Math.abs(( (diff_red * percent) + start_red )).toString(16).split('.')[0];
   diff_green = Math.abs(( (diff_green * percent) + start_green )).toString(16).split('.')[0];
   diff_blue = Math.abs(( (diff_blue * percent) + start_blue )).toString(16).split('.')[0];

   // ensure 2 digits by color
   if( diff_red.length == 1 )
     diff_red = '0' + diff_red

   if( diff_green.length == 1 )
     diff_green = '0' + diff_green

   if( diff_blue.length == 1 )
     diff_blue = '0' + diff_blue

     return '#' + diff_red + diff_green + diff_blue;
 };
