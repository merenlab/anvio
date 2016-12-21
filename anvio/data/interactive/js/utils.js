/**
 * Utility functions for anvi'o interactive interface
 *
 *  Author: Özcan Esen <ozcanesen@gmail.com>
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
// http://stackoverflow.com/questions/3019278/any-way-to-specify-the-base-of-math-log-in-javascript
function log10(val) {
  return Math.log(val) / Math.LN10;
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

   diff_red = ( (diff_red * percent) + start_red ).toString(16).split('.')[0];
   diff_green = ( (diff_green * percent) + start_green ).toString(16).split('.')[0];
   diff_blue = ( (diff_blue * percent) + start_blue ).toString(16).split('.')[0];

   // ensure 2 digits by color
   if( diff_red.length == 1 )
     diff_red = '0' + diff_red

   if( diff_green.length == 1 )
     diff_green = '0' + diff_green

   if( diff_blue.length == 1 )
     diff_blue = '0' + diff_blue

   return '#' + diff_red + diff_green + diff_blue;
 };

function fire_up_ncbi_blast(item_name, program, database, target)
{
    if (["gene", "contig"].indexOf(target) < 0){
        console.log("fire_up_ncbi_blast: Unrecognized target. Target must be either 'gene', or 'contig'.");
        return;
    }

    var post_variables = {
        'PROGRAM': 'blastn',
        'DATABASE': 'nr',
        'QUERY': '',
        'BLAST_PROGRAMS': 'megaBlast',
        'PAGE_TYPE': 'BlastSearch',
        'SHOW_DEFAULTS': 'on',
        'SHOW_OVERVIEW': 'on',
        'LINK_LOC': 'blasthome',
        'MAX_NUM_SEQ': '100',
        "FORMAT_NUM_ORG": "1",
        "CONFIG_DESCR": "2,3,4,5,6,7,8",
        "CLIENT": "web" ,
        "SERVICE": "plain",
        "CMD": "request",
        "PAGE": "MegaBlast",
        "MEGABLAST": "on" ,
        "WWW_BLAST_TYPE": "newblast",
        "DEFAULT_PROG": "megaBlast",
        "SELECTED_PROG_TYPE": "megaBlast",
        "SAVED_SEARCH": "true",
        "NUM_DIFFS": "0",
        "NUM_OPTS_DIFFS": "0",
        "PAGE_TYPE": "BlastSearch",
        "USER_DEFAULT_PROG_TYPE": "megaBlast"
    }

    if (typeof program !== 'undefined')
        post_variables['PROGRAM'] = program;
    
    if (typeof database !== 'undefined')
        post_variables['DATABASE'] = database;

    var blast_window = window.open('about:blank', '_blank');

    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/' + target + '/' + item_name + '?timestamp=' + new Date().getTime(),
        success: function(data) {
            if ('error' in data){
                toastr.error(data['error'], "", { 'timeOut': '0', 'extendedTimeOut': '0' });
            } else {
                post_variables['QUERY'] = '>' + data['header'] + '\n' + data['sequence'];

                var form = document.createElement('form');
                
                form.action = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi';
                form.method = 'POST';

                for (name in post_variables)
                {
                    $(form).append('<input type="hidden" name="' + name + '" value="' + post_variables[name] + '" />');
                }

                blast_window.document.body.appendChild(form);
                form.submit();
            }
        }
    });
}

//--------------------------------------------------------------------------------------------------

function showDraggableDialog(title, content, updateOnly)
{
    var randomID = title.hashCode();

    if (updateOnly)
    {
        if (checkObjectExists('#modal' + randomID))
            $('#modal' + randomID).find('.modal-body').html(content);
        
        return;
    }
    else
    {
        var template = '<div class="modal" id="modal' + randomID + '" data-backdrop="false" style="pointer-events: none;"> \
            <div class="modal-dialog" style="pointer-events: all;"> \
                <div class="modal-content no-shadow"> \
                    <div class="modal-header"> \
                        <button type="button" class="close" data-dismiss="modal" aria-hidden="true">&times;</button> \
                         <h4 class="modal-title">' + title + '</h4> \
                    </div> \
                    <div class="modal-body"> \
                        ' + content + ' \
                    </div> \
                    <div class="modal-footer"> \
                        <button type="button" class="btn btn-default" data-dismiss="modal">Close</button> \
                    </div> \
                </div> \
            </div> \
        </div>';
    }

    $('body').append(template);
    $('#modal' + randomID).modal({'show': true, 'backdrop': false, 'keyboard': false}).find('.modal-dialog').draggable({handle: '.modal-header'});
    $('#modal' + randomID).on('hidden.bs.modal', function () {
        $(this).remove();
    });
}

//--------------------------------------------------------------------------------------------------

function checkObjectExists(selector)
{
    if ($(selector).length > 0)
        return true;
    return false;
}

//--------------------------------------------------------------------------------------------------
// https://stackoverflow.com/questions/7616461/generate-a-hash-from-string-in-javascript-jquery
String.prototype.hashCode = function() {
  var hash = 0, i, chr, len;
  if (this.length == 0) return hash;
  for (i = 0, len = this.length; i < len; i++) {
    chr   = this.charCodeAt(i);
    hash  = ((hash << 5) - hash) + chr;
    hash |= 0; // Convert to 32bit integer
  }
  return Math.abs(hash);
};


//--------------------------------------------------------------------------------------------------

function checkBackgroundProcess()
{
    var message = "It seems background process is down or changed, you may lose your work.";

    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/session_id?timestamp=' + new Date().getTime(),
        success: function (data) {
            if (data != unique_session_id)
            {
                toastr.error(message, "", { 'timeOut': '0', 'extendedTimeOut': '0' });
                clearTimeout(ping_timer);
            }
        },
        error: function(data) {
            toastr.error(message, "", { 'timeOut': '0', 'extendedTimeOut': '0' });
            clearTimeout(ping_timer);
        }
    });
}

//--------------------------------------------------------------------------------------------------
// http://stackoverflow.com/questions/1303646/check-whether-variable-is-number-or-string-in-javascript
function isNumber (o) {
  return ! isNaN (o-0);
}

//--------------------------------------------------------------------------------------------------
// http://stackoverflow.com/questions/3169786/clear-text-selection-with-javascript
function clearTextSelection() {
    if (window.getSelection) {
      if (window.getSelection().empty) {  // Chrome
        window.getSelection().empty();
      } else if (window.getSelection().removeAllRanges) {  // Firefox
        window.getSelection().removeAllRanges();
      }
    } else if (document.selection) {  // IE?
      document.selection.empty();
    }
}

//--------------------------------------------------------------------------------------------------
function ctype_alnum (str)
{
    return (str.match(/^[a-z0-9]+$/i) != null);
}

//--------------------------------------------------------------------------------------------------
function strip(html)
{
   var tmp = document.createElement("DIV");
   tmp.innerHTML = html;
   return tmp.textContent || tmp.innerText || "";
}

//--------------------------------------------------------------------------------------------------
function clearMinMax(selectbox) 
{
    var tr = $(selectbox).parent().parent();

    $(tr).find('.input-min').val('0').prop('disabled', true);
    $(tr).find('.input-max').val('0').prop('disabled', true);     
}

function togglePickerStart(selectbox, togglePicker)
{
    var tr = $(selectbox).parent().parent();

    if(selectbox.value=='intensity' || selectbox.value=='text') {  
        $(tr).find('.picker_start').css('visibility', 'visible');
        if (togglePicker) {
            $(tr).find('.picker_end').css('visibility', 'visible');
            $(tr).find('.input-height').css('visibility', 'hidden');
            $('.max-font-size-input').show();
        }
    } else { 
        $(tr).find('.picker_start').css('visibility', 'hidden');
        if (togglePicker) {
            $(tr).find('.picker_end').css('visibility', 'hidden');
            $(tr).find('.input-height').css('visibility', 'visible');
            $(tr).find('.input-height').val('30');
        }
    }  
}

/* Poor man's timer.
 * 
 *     function ...(...) {
 *         var my_timer = new BasicTimer('my_func');
 *         (...)
 *         my_timer.getDeltaSeconds('X happened');
 *         (...)
 *         my_timer.getDeltaSeconds('Y happened');
 *         (...)
 *         my_timer.getDeltaSeconds('End');
 *     }
 *
 */ 
function BasicTimer(name) {
    this.name = name;
    this.start = new Date().getTime();
    this.previousDelta = this.start;

    this.getDeltaSeconds = function(event, consoleOutput) {
        this.now = new Date().getTime();
        deltaSecondsStart = (this.now - this.start) / 1000;
        deltaSecondsPrev = (this.now - this.previousDelta) / 1000;

        this.previousDelta = this.now;
        
        consoleOutput = typeof consoleOutput !== 'undefined' ? consoleOutput: true;

        prettyText = this.name + ' [' + event + ']: ' + readableNumber(deltaSecondsPrev) + ' seconds (' + readableNumber(deltaSecondsStart) + ' seconds since beginning)';

        if(consoleOutput)
            console.log(prettyText);

        return {'deltaSecondsStart': deltaSecondsStart, 'deltaSecondsPrev': deltaSecondsPrev, 'prettyText': prettyText};
    };
}

//--------------------------------------------------------------------------------------------------
// source: https://gist.github.com/cjthompson/9140248
function readableNumber(num) {
    if(num == 0)
        return 0;
    if(num < 1)
        return num;
    var s = ['', 'K', 'M', 'G'];
    var e = Math.floor(Math.log(num) / Math.log(1000));
    return (num / Math.pow(1000, e)).toPrecision(3) + s[e];
}

//--------------------------------------------------------------------------------------------------
function linePath(p0, p1)
{
    var path = 'M ' + p0['x'] + ' ' + p0['y'] + ' ' + p1['x'] + ' ' + p1['y'];
    return path;
}

//--------------------------------------------------------------------------------------------------
function distance(p0, p1)
{
    return Math.sqrt(Math.pow(p1['x'] - p0['x'],2) + Math.pow(p1['y'] + p0['y'],2));
}

// http://stackoverflow.com/questions/498970/how-do-i-trim-a-string-in-javascript
if (!String.prototype.trim)
{
    String.prototype.trim=function(){return this.replace(/^\s+|\s+$/g, '');};
}

Math.toRadians = function(degrees) {
    return degrees * Math.PI / 180;
};

Math.toDegrees = function(radians) {
    return radians * 180 / Math.PI;
};


//--------------------------------------------------------------------------------------------------
//  Iterator
//--------------------------------------------------------------------------------------------------
function NodeIterator(root)
{
    this.root = root;
    this.cur = null;
    this.stack = [];
}

//--------------------------------------------------------------------------------------------------
NodeIterator.prototype.Begin = function() 
{
    if (this.root.constructor === Array)
    {
        this.cur = 0;
        return label_to_node_map[this.root[0]];
    }
    this.cur = this.root;
    while (this.cur.child)
    {
        this.stack.push(this.cur);
        this.cur = this.cur.child;
    }
    return this.cur;
}

//--------------------------------------------------------------------------------------------------
NodeIterator.prototype.Next = function() 
{
    if (this.root.constructor === Array)
    {
        this.cur = this.cur + 1;
        if (this.cur >= this.root.length)
        {
            return null;
        }
        return label_to_node_map[this.root[this.cur]];
    }
    if (this.stack.length == 0)
    {
        this.cur = null;
    }
    else
    {
        if (this.cur.sibling)
        {
            var p = this.cur.sibling;
            while (p.child)
            {
                this.stack.push(p);
                p = p.child;
            }
            this.cur = p;
        }
        else
        {
            this.cur = this.stack.pop();
        }
    }
    return this.cur;
}

//--------------------------------------------------------------------------------------------------
PreorderIterator.prototype = new NodeIterator;

function PreorderIterator()
{
    NodeIterator.apply(this, arguments)
};

//--------------------------------------------------------------------------------------------------
PreorderIterator.prototype.Begin = function() 
{
    this.cur = this.root;
    return this.cur;
}

//--------------------------------------------------------------------------------------------------
PreorderIterator.prototype.Next = function() 
{
    if (this.cur.child)
    {
        this.stack.push(this.cur);
        this.cur = this.cur.child;
    }
    else
    {
        while (this.stack.length > 0 && this.cur.sibling == null)
        {
            this.cur = this.stack.pop();
        }
        if (this.stack.length == 0)
        {
            this.cur = null;
        }
        else
        {
            this.cur = this.cur.sibling;
        }
    }
    return this.cur;
}

//---------------------------------------------------------
// layerdata operations
//---------------------------------------------------------

function removeSingleParents()
{
    // layerdata and parameter count is global

    for (var i = 1; i < parameter_count; i++) 
    {
        if (layerdata[0][i] == '__parent__') 
        {
            var parent_count_dict = {};
            for (var j=1; j < layerdata.length; j++)
            {
                if (layerdata[j][i]=='')
                    continue;

                if (typeof parent_count_dict[layerdata[j][i]] === 'undefined')
                {
                    parent_count_dict[layerdata[j][i]] = 1;
                }
                else
                {
                    parent_count_dict[layerdata[j][i]]++;
                }
            }

            $.each(parent_count_dict, function(parent_name, count)
            {
                if (count==1)
                {
                    for (var j=1; j < layerdata.length; j++)
                    {
                        if (layerdata[j][i]==parent_name)
                        {
                            layerdata[j][i]='';
                        }
                    }
                }
            });
        }
    }
}

//---------------------------------------------------------
// jquery table sort helper
//---------------------------------------------------------
var fixHelperModified = function(e, tr) {
    var $originals = tr.children();
    var $helper = tr.clone();
    $helper.children().each(function(index) {
        $(this).width($originals.eq(index).width());
    });
    return $helper;
};

//---------------------------------------------------------
//  zoom and scale
//---------------------------------------------------------
function getMatrix() {
    var viewport = document.getElementById('viewport');
    return viewport.getAttribute('transform').split('(')[1].split(')')[0].split(',').map(parseFloat);
}

function setMatrix(matrix) {
    var viewport = document.getElementById('viewport');
    viewport.setAttribute('transform', 'matrix(' + matrix.join(',') + ')');
}

function zoom(scale) {
    matrix = getMatrix(viewport);

    for (var i = 0; i < matrix.length; i++) {
        matrix[i] *= scale;
    }

    bbox = viewport.getBBox();

    matrix[4] += (1 - scale) * VIEWER_WIDTH / 2;
    matrix[5] += (1 - scale) * VIEWER_HEIGHT / 2;

    setMatrix(matrix);
}

function pan(dx, dy) {
    matrix = getMatrix();

    matrix[4] += dx;
    matrix[5] += dy;

    setMatrix(matrix);
}

function zoom_reset() {
    baseMatrix = [1 * scale, 0, 0, 1 * scale, VIEWER_WIDTH / 2 - (bbox.x + bbox.width / 2) * scale, VIEWER_HEIGHT / 2 - (bbox.y + bbox.height / 2) * scale];
    setMatrix(baseMatrix);
}
