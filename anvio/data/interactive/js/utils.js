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

function fire_up_ncbi_blast(contig_name, program, database)
{
    var blast_window = window.open('about:blank', '_blank');

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

    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/contig/' + contig_name + '?timestamp=' + new Date().getTime(),
        success: function(data) {
            post_variables['QUERY'] = '>' + contig_name + '\n' + data;
            
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
    });
}

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
                alert(message);
                clearTimeout(ping_timer);
            }
        },
        error: function(data) {
            alert(message);
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
    var id = $(selectbox).attr('id').replace('normalization', '');

    $('#min' + id).val('0').prop('disabled', true);
    $('#max' + id).val('0').prop('disabled', true);     
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
// Metadata operations
//---------------------------------------------------------

function removeSingleParents()
{
    // metadata and parameter count is global

    for (var i = 1; i < parameter_count; i++) 
    {
        if (metadata[0][i] == '__parent__') 
        {
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
    baseMatrix = [1 * scale, 0, 0, 1 * scale, VIEWER_WIDTH / 2, VIEWER_HEIGHT / 2];
    setMatrix(baseMatrix);
}
