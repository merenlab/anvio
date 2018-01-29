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

function is_large_angle(a, b) {
    return (Math.abs(b - a) > Math.PI) ? 1 : 0;
}


function get_sequence_and_blast(item_name, program, database, target) {
    $.ajax({
        type: 'GET',
        cache: false,
        // async: false is important here. DO NOT REMOVE.
        // If direct call chain between event handler and the code that opens new window is broken
        // chrome popup blocker will not allow opening new window.
        // async: false does not use asynchronus callbacks so protects direct call chain.
        async: false,
        url: '/data/' + target + '/' + item_name + '?timestamp=' + new Date().getTime(),
        success: function(data) {
            if ('error' in data){
                toastr.error(data['error'], "", { 'timeOut': '0', 'extendedTimeOut': '0' });
            } else {
              var sequence = '>' + data['header'] + '\n' + data['sequence'];
              fire_up_ncbi_blast(sequence, program, database, target)
            }
        }
    });
}

function fire_up_ncbi_blast(sequence, program, database, target)
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
        "USER_DEFAULT_PROG_TYPE": "megaBlast",
        "QUERY": sequence,
    }

    if (typeof program !== 'undefined')
        post_variables['PROGRAM'] = program;
    
    if (typeof database !== 'undefined')
        post_variables['DATABASE'] = database;

    var blast_window = window.open('about:blank', '_blank');
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
//--------------------------------------------------------------------------------------------------

function generate_inspect_link(type, item_name) {
    if (self == top) {
        // local anvio
        var url = window.location.href.split('?')[0];

        if (url.endsWith('index.html')) {
            // on index page
            if (type == 'inspect') {
                return 'charts.html?id=' + item_name;
            } 
            else if (type == 'geneclusters') {
                return 'geneclusters.html?id=' + item_name;
            }
        }
        else
        {
            // on charts or gene cluster page, so changing the ?id= part enough
            return url + '?id=' + item_name;
        }
    }
    else
    {
        // anvi server
        var url = window.parent.location.href.split('?')[0];
        var new_url = "";

        if (url.endsWith('/inspect') || url.endsWith('/geneclusters')) {
            // on charts or gene cluster page
            new_url = url;
        }
        else
        {
            // on main page
            new_url = url + '/' + type;
        }

        new_url = new_url + '?id=' + item_name;

        var view_key = request_prefix.substr(request_prefix.lastIndexOf('/') + 1, request_prefix.length);
        if (view_key != 'no_view_key') {
            new_url = new_url + '&view_key=' + view_key;
        }

        return new_url;
    }
}

//--------------------------------------------------------------------------------------------------
function getParameterByName(name, url) {
    if (!url) url = window.location.href;
    name = name.replace(/[\[\]]/g, "\\$&");
    var regex = new RegExp("[?&]" + name + "(=([^&#]*)|&|#|$)"),
        results = regex.exec(url);
    if (!results) return null;
    if (!results[2]) return '';
    return decodeURIComponent(results[2].replace(/\+/g, " "));
}

//-------------------------------------------------------------------------------------------------
function renderMarkdown(content) {
    var renderer = new marked.Renderer();

    renderer.link = function( href, title, text ) {
        if (href.startsWith('item://')) {
            var item_name = href.split('//')[1];

            var html = '<a href="#" class="item-link">' + text + '<span class="tooltiptext"> \
                <span href="#" onclick="highlighted_splits = [\'' + item_name + '\']; redrawBins();">HIGHLIGHT</span>';

            if (mode == 'full' | mode == 'pan') {
                var target = (mode == 'pan') ? 'inspect_gene_cluster' : 'inspect_contig';
                html += ' | <span href="#" onclick="context_menu_target_id = label_to_node_map[\'' + item_name + '\'].id; \
                                                 menu_callback(\'' + target + '\');">INSPECT</span>';
            }

            return html + '</span></a>';
        }
        return '<a target="_blank" href="' + href + '" title="' + title + '">' + text + '</a>';
    }

    return marked(content, { renderer:renderer });
}

//--------------------------------------------------------------------------------------------------
function getCookie(name) {
    var cookieValue = null;
    if (document.cookie && document.cookie !== '') {
        var cookies = document.cookie.split(';');
        for (var i = 0; i < cookies.length; i++) {
            var cookie = jQuery.trim(cookies[i]);
            // Does this cookie string begin with the name we want?
            if (cookie.substring(0, name.length + 1) === (name + '=')) {
                cookieValue = decodeURIComponent(cookie.substring(name.length + 1));
                break;
            }
        }
    }
    return cookieValue;
}

//--------------------------------------------------------------------------------------------------
// https://stackoverflow.com/questions/14573223/set-cookie-and-get-cookie-with-javascript
function createCookie(name,value,days) {
    var expires = "";
    var date = new Date();
    if (days == -1) {
        date.setTime(date.getTime() + (10*365*24*60*60*1000));
    } else if (days) {
        date.setTime(date.getTime() + (days*24*60*60*1000));
    }
    expires = "; expires=" + date.toUTCString();
    document.cookie = name + "=" + value + expires + "; path=/";
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

    if(selectbox.value=='intensity' || selectbox.value=='line' || selectbox.value=='text') {  
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
function getReadableSeqSizeString(seqSizeInBases) {
    // function based on answer at http://stackoverflow.com/questions/10420352/converting-file-size-in-bytes-to-human-readable
    var i = -1;
    var baseUnits = [' K', ' M', ' G', ' T'];

    // if the number is less than a K, return as is
    if (seqSizeInBases < 1000)
        return seqSizeInBases;

    do {
        seqSizeInBases = seqSizeInBases / 1000;
        i++;
    } while (seqSizeInBases >= 1000);

    return Math.round(seqSizeInBases) + baseUnits[i];
};

//--------------------------------------------------------------------------------------------------
function getCommafiedNumberString(number, decimals, dec_point, thousands_sep) {
    // function modified from https://stackoverflow.com/a/2901136

    if(isNaN(parseInt(number)))
        return number;

    var n = !isFinite(+number) ? 0 : +number,
        prec = !isFinite(+decimals) ? 2 : Math.abs(decimals),
        sep = (typeof thousands_sep === 'undefined') ? ',' : thousands_sep,
        dec = (typeof dec_point === 'undefined') ? '.' : dec_point,
        toFixedFix = function (n, prec) {
            // Fix for IE parseFloat(0.55).toFixed(0) = 0;
            var k = Math.pow(10, prec);
            return Math.round(n * k) / k;
        },
        s = (prec ? toFixedFix(n, prec) : Math.round(n)).toString().split('.');
    if (s[0].length > 3) {
        s[0] = s[0].replace(/\B(?=(?:\d{3})+(?!\d))/g, sep);
    }
    if ((s[1] || '').length < prec) {
        s[1] = s[1] || '';
        s[1] += new Array(prec - s[1].length + 1).join('0');
    }

    if(s[1] > 0)
        return s.join(dec);
    else
        return s[0];
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
    scale = Math.min(VIEWER_WIDTH / bbox.width, VIEWER_HEIGHT / bbox.height) * 0.80;
    baseMatrix = [1 * scale, 0, 0, 1 * scale, VIEWER_WIDTH / 2 - (bbox.x + bbox.width / 2) * scale, VIEWER_HEIGHT / 2 - (bbox.y + bbox.height / 2) * scale];
    setMatrix(baseMatrix);
}
