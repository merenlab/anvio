/**
 *  Admin functions for anvi'o interactive interface
 *
 *  Author: Tobias Paczian <tobiaspaczian@googlemail.com>
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

function initContent() {
    window.table = { "limit": 25, "offset": 0, "order": "lastname", "dir": "ASC", "filter": {} };
    var content = document.getElementById('content');
    var html = [];
    if (user && user.clearance == 'admin') {
	html.push('<h3>Welcome back, master '+user.lastname+'</h3><p>I am awaiting your command.</p><h3>User Table</h3><div id="usertable"><div class="progress" style="width: 50%; margin-top: 100px; margin-left: 25%; margin-right: 25%;"><div class="progress-bar progress-bar-striped active" role="progressbar" aria-valuenow="100" aria-valuemin="0" aria-valuemax="100" style="width: 100%"><span class="sr-only">Loading...</span></div></div></div>');
	updateTable();
    } else {
	html.push('<div class="alert alert-danger col-sm-6" role="alert" style="margin-top: 100px;">You are not authorized to view this page.</div>');
    }
    content.innerHTML = html.join('\n');
}

function updateTable(key, value) {
    if (key) {
	if (key == 'offset') {
	    if (value == '-') {
		table.offset -= table.limit;
	    } else if (value == '--') {
		table.offset = 0;
	    } else if (value == '+') {
		table.offset += table.limit;
	    } else if (value == '++' && table.total) {
		table.offset = table.total - table.limit;
	    } else {
		table.offset = value;
	    }
	    if (table.total && table.offset >= table.total) {
		table.offset = table.total - 1;
	    }
	    if (table.offset < 0) {
		table.offset = 0;
	    }
	} else if (key == 'orderasc') {
	    table.order = value;
	    table.dir = 'ASC';
	} else if (key == 'orderdesc') {
	    table.order = value;
	    table.dir = 'DESC';
	} else if (key == 'limit') {
	    table.limit = value;
	}
    }
    
    filter = [];
    for (var i in table.filter) {
	if (table.filter.hasOwnProperty(i) && table.filter[i] !== null) {
	    filter.push(i+"="+table.filter[i]);
	}
    }

    if (filter.length > 0) {
	filter = "&" + filter.join("&");
    } else {
	filter = "";
    }

    $.ajax({
	url : '/adminData?limit='+table.limit+'&offset='+table.offset+'&order='+table.order+'&direction='+table.dir+filter,
	type : 'GET',
	processData: false,
	contentType: false,
	complete : function(jqXHR) {
	    var r = JSON.parse(jqXHR.responseText);
	    var data = r.data;
	    table.total = parseInt(r.total);
	    var html = [];

	    // table control structure
	    html.push("<div style='margin-bottom: 3px; text-align: center;'>"+(r.offset>0 ? "<button class='btn btn-xs btn-default pull-left' onclick='updateTable(\"offset\", \"--\");'><span class='glyphicon glyphicon-fast-backward'></span></button><button class='btn btn-xs btn-default pull-left' onclick='updateTable(\"offset\", \"-\");'><span class='glyphicon glyphicon-step-backward'></span></button>" : "")+"showing row "+(parseInt(r.offset) + 1)+" to "+(parseInt(r.limit) - parseInt(r.offset) < parseInt(r.total) ? parseInt(r.limit) : parseInt(r.total))+" of "+r.total+(parseInt(r.offset) + parseInt(r.limit) < parseInt(r.total) ? "<button class='btn btn-xs btn-default pull-right' onclick='updateTable(\"offset\", \"++\");'><span class='glyphicon glyphicon-fast-forward'></span></button><button class='btn btn-xs btn-default pull-right' onclick='updateTable(\"offset\", \"+\");'><span class='glyphicon glyphicon-step-forward'></span></button>" : "")+"</div>");

	    // generic columns
	    var columns = [ 'firstname', 'lastname', 'login', 'email', 'affiliation', 'projects', 'date', 'clearance' ];

	    // start table
	    html.push('<table class="table table-condensed table-bordered table-striped table-hover"><thead><tr>');

	    // header
	    for (var i=0; i<columns.length; i++) {
		html.push('<th><span onclick="this.style.display=\'none\';this.nextSibling.style.display=\'\';this.nextSibling.focus()" style="cursor: pointer;" title="click to filter">'+columns[i]+'</span><input type="text" style="display: none; font-weight: normal; padding: 0px;" class="col-md-10" value="'+(table.filter.hasOwnProperty(columns[i]) ? table.filter[columns[i]] : "")+'" onkeypress="checkTableFilter(event)"><span class="glyphicon glyphicon-triangle-top pull-right" onclick="updateTable(\'orderasc\', \''+columns[i]+'\')" style="cursor: pointer;"></span><span class="glyphicon glyphicon-triangle-bottom pull-right" onclick="updateTable(\'orderdesc\', \''+columns[i]+'\')" style="cursor: pointer;"></th>')
	    }
	    html.push('<th>impersonate</th></tr></thead><tbody>');

	    // rows
	    for (var i=0; i<data.length; i++) {
		html.push('<tr>');
		for (var h=0; h<columns.length; h++) {
		    html.push('<td>'+data[i][columns[h]]+'</td>');
		}
		html.push('<td><button class="btn btn-default btn-xs" onclick="impersonate(\''+data[i].login+'\');">impersonate</button></td></tr>');
	    }
	    html.push('</tbody></table>');

	    // update the DOM
	    document.getElementById('usertable').innerHTML = html.join('');
	}
    });
}

function checkTableFilter (event) {
    event = event || window.event;
    var input = event.target;
    var text = event.target.previousSibling;
    if (event.keyCode == 13) {
	if (input.value.length > 0) {
	    table.filter[text.innerHTML] = input.value;
	} else {
	    delete table.filter[text.innerHTML];
	}
	updateTable();
    } else if (event.keyCode == 27 || event.keyCode == 9) {
	input.style.display = 'none';
	text.style.display = '';
    }
}

function impersonate(login) {
    var formData = new FormData();
    formData.append('login', login);
    $.ajax({
	url : '/impersonate',
	type : 'POST',
	processData: false,
	contentType: false,
	data : formData,
	complete : function(jqXHR) {
	    window.location = 'user.html';
	}
    });
}
