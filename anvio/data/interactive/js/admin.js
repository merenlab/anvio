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
    var content = document.getElementById('content');
    var html = [];
    if (user && user.clearance == 'admin') {
	html.push('<h3>Welcome back, master '+user.lastname+'</h3><p>I am awaiting your command.</p><h3>User Table</h3><div id="usertable"><div class="progress" style="width: 50%; margin-top: 100px; margin-left: 25%; margin-right: 25%;"><div class="progress-bar progress-bar-striped active" role="progressbar" aria-valuenow="100" aria-valuemin="0" aria-valuemax="100" style="width: 100%"><span class="sr-only">Loading...</span></div></div></div>');
	$.ajax({
	    url : '/adminData',
	    type : 'GET',
	    processData: false,
	    contentType: false,
	    complete : function(jqXHR) {
		var data = JSON.parse(jqXHR.responseText);
		var html = [];
		html.push('<table class="table table-condensed table-bordered table-striped table-hover"><tr><th>firstname</th><th>lastname</th><th>login</th><th>email</th><th>affiliation</th><th>project</th><th>date</th><th>ip</th><th>clearance</th><th>accepted</th><th>impersonate</th></tr>');
		for (var i=0; i<data.length; i++) {
		    html.push('<tr><td>'+data[i].firstname+'</td><td>'+data[i].lastname+'</td><td>'+data[i].login+'</td><td>'+data[i].email+'</td><td>'+data[i].affiliation+'</td><td>'+(data[i].project || '- none -')+'</td><td>'+data[i].date+'</td><td>'+data[i].ip+'</td><td>'+data[i].clearance+'</td><td>'+data[i].accepted+'</td><td><button class="btn btn-default btn-xs" onclick="impersonate(\''+data[i].login+'\');">impersonate</button></td></tr>');
		}
		html.push('</table>');
		document.getElementById('usertable').innerHTML = html.join('\n');
	    }
	});
    } else {
	html.push('<div class="alert alert-danger col-sm-6" role="alert" style="margin-top: 100px;">You are not authorized to view this page.</div>');
    }
    content.innerHTML = html.join('\n');
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
