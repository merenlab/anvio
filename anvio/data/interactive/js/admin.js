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
    window.tables = {
                "user": { "limit": 50, "offset": 0, "order": "lastname", "dir": "ASC", "filter": {}, "url": "/adminData",
                "columns": [ 'firstname', 'lastname', 'login', 'email', 'affiliation', 'projects', 'date', 'clearance', 'visit', 'details' ],
                "special": { "column": "details", "code": '<button class="btn btn-default btn-xs" onclick="showUserDetails(\'$$\');">show details</button>', "data": 'login' },
                "response": null},

		            "project": { "limit": 100, "offset": 0, "order": "user", "dir": "ASC", "filter": {}, "url": "/adminProjectData",
                "columns": [ 'name', 'user', 'path', 'description', 'details' ],
                "special": { "column": "details", "code": '<button class="btn btn-default btn-xs" onclick="getProjectDetails(\'$$\');">show details</button>', "data": 'name' },
                "response": null }
            };
    var content = document.getElementById('content');
    var html = [];
    if (user && user.clearance == 'admin') {
	html.push('<h3>Welcome back ' + user.firstname + '</h3><h4 align="center">Users</h4><div id="usertable"><div class="progress" style="width: 50%; margin-top: 100px; margin-left: 25%; margin-right: 25%;"><div class="progress-bar progress-bar-striped active" role="progressbar" aria-valuenow="100" aria-valuemin="0" aria-valuemax="100" style="width: 100%"><span class="sr-only">Loading...</span></div></div></div>');
	html.push('<h4 align="center">Projects</h4><div id="projecttable"><div class="progress" style="width: 50%; margin-top: 100px; margin-left: 25%; margin-right: 25%;"><div class="progress-bar progress-bar-striped active" role="progressbar" aria-valuenow="100" aria-valuemin="0" aria-valuemax="100" style="width: 100%"><span class="sr-only">Loading...</span></div></div></div>');
	html.push('<div id="details"></div>');
	content.innerHTML = html.join('\n');
	updateTable('user').then(function(){updateTable('project');});
    } else {
	html.push('<div class="alert alert-danger col-sm-6" role="alert" style="margin-top: 100px;">You are not authorized to view this page.</div>');
	content.innerHTML = html.join('\n');
    }
}

function noUser() {
    window.location = 'home.html';
}

function updateTable(id, key, value) {
    var table = tables[id];
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

    return $.ajax({
	url : table.url+'?limit='+table.limit+'&offset='+table.offset+'&order='+table.order+'&direction='+table.dir+filter,
	type : 'GET',
	index: id,
	processData: false,
	contentType: false,
	complete : function(jqXHR) {
	    var table = tables[this.index];
	    var r = JSON.parse(jqXHR.responseText);
	    if (r.status == 'error') {
		alert(r.message);
		return;
	    }
	    r = r.data;
	    var data = r.data;
	    table.response = data;
	    table.total = parseInt(r.total);
	    var html = [];

	    // table control structure
	    html.push("<div style='margin-bottom: 3px; text-align: center;'>"+(r.offset>0 ? "<button class='btn btn-xs btn-default pull-left' onclick='updateTable(\""+this.index+"\", \"offset\", \"--\");'><span class='glyphicon glyphicon-fast-backward'></span></button><button class='btn btn-xs btn-default pull-left' onclick='updateTable(\""+this.index+"\", \"offset\", \"-\");'><span class='glyphicon glyphicon-step-backward'></span></button>" : "")+"showing row "+(parseInt(r.offset) + 1)+" to "+(parseInt(r.limit) + parseInt(r.offset) < parseInt(r.total) ? parseInt(r.limit) : parseInt(r.total))+" of "+r.total+(parseInt(r.offset) + parseInt(r.limit) < parseInt(r.total) ? "<button class='btn btn-xs btn-default pull-right' onclick='updateTable(\""+this.index+"\", \"offset\", \"++\");'><span class='glyphicon glyphicon-fast-forward'></span></button><button class='btn btn-xs btn-default pull-right' onclick='updateTable(\""+this.index+"\", \"offset\", \"+\");'><span class='glyphicon glyphicon-step-forward'></span></button>" : "")+"</div>");

	    // generic columns
	    var columns = table.columns;

	    // start table
	    html.push('<table class="table table-condensed table-bordered table-striped table-hover"><thead><tr>');

	    // header
	    for (var i=0; i<columns.length; i++) {
		if (table.special.column == columns[i]) {
		    html.push('<th>'+columns[i]+'</th>');
		} else {
		    html.push('<th><span onclick="this.style.display=\'none\';this.nextSibling.style.display=\'\';this.nextSibling.focus()" style="cursor: pointer;" title="click to filter">'+columns[i]+'</span><input type="text" style="display: none; font-weight: normal; padding: 0px;" class="col-md-10" value="'+(table.filter.hasOwnProperty(columns[i]) ? table.filter[columns[i]] : "")+'" onkeypress="checkTableFilter(event, \''+this.index+'\')"><span class="glyphicon glyphicon-triangle-top pull-right" onclick="updateTable(\''+this.index+'\', \'orderasc\', \''+columns[i]+'\')" style="cursor: pointer;"></span><span class="glyphicon glyphicon-triangle-bottom pull-right" onclick="updateTable(\''+this.index+'\', \'orderdesc\', \''+columns[i]+'\')" style="cursor: pointer;"></th>')
		}
	    }
	    html.push('</tr></thead><tbody>');

	    // rows
	    for (var i=0; i<data.length; i++) {
		html.push('<tr>');
		for (var h=0; h<columns.length; h++) {
		    if (table.special.column == columns[h]) {
			html.push('<td>'+table.special.code.replace(/\$\$/g, data[i][table.special.data])+'</td>');
		    } else {
			html.push('<td>'+data[i][columns[h]]+'</td>');
		    }
		}
		html.push('</tr>');
	    }
	    html.push('</tbody></table>');

	    // update the DOM
	    document.getElementById(this.index+'table').innerHTML = html.join('');
	}
    });
}

function checkTableFilter (event, index) {
    var table = tables[index];
    event = event || window.event;
    var input = event.target;
    var text = event.target.previousSibling;
    if (event.keyCode == 13) {
	if (input.value.length > 0) {
	    table.filter[text.innerHTML] = input.value;
	} else {
	    delete table.filter[text.innerHTML];
	}
	updateTable(index);
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

function getProjectDetails(projectname, user) {
    var project = null;
    var projects = tables['project'].response;
    if (user) {
	project = { 'name': projectname, 'user': user };
    } else {
	for (var i=0; i<projects.length; i++) {
	    if (projects[i].name == projectname) {
		project = projects[i];
		break;
	    }
	}
    }
    if (project) {
	$.ajax({
	url : '/adminProjectDetails?project='+project.name+'&user='+project.user,
	type : 'GET',
	processData: false,
	contentType: false,
	complete : function(jqXHR) {
	    var result = JSON.parse(jqXHR.responseText);
	    if (result.status == 'ok') {
		projectData = result.data;
		showProjectDetails();
	    } else {
		toastr.error(result.message);
	    }
	}
	});
    } else {
	toastr.error('invalid project link');
    }
}

function showProjectDetails() {
    var html = [];
    html.push('<div class="panel panel-primary"><div class="panel-heading"><h4 class="panel-title">Project Details</h4></div><div class="panel-body">');
    // project base data
    html.push('<button style="float: right;" class="btn btn-sm btn-danger" onclick="deleteProject(\''+projectData.name+'\', \''+projectData.user.login+'\')"><span class="glyphicon glyphicon-trash"></span></button>');
    html.push('<h4 style="margin-top: 0px;">Project '+projectData.name+'</h4>');
    html.push('<p>' + projectData.description + '</p>');

    // project user
    html.push('<h4><span class="glyphicon glyphicon-user" style="margin-right: 10px;"></span>Owner</h4><p><b>'+projectData.user.firstname+' '+projectData.user.lastname+'</b></p><table class="table table-condensed col-sm-8">');
    var fields = [ 'login', 'email', 'affiliation', 'visit', 'project' ];
    for (var i=0; i<fields.length; i++) {
	html.push('<tr><th>'+fields[i]+'</th><td>'+projectData.user[fields[i]]+'</td></tr>');
    }
    html.push('</table>');

    // project views
    html.push('<h4><span class="glyphicon glyphicon-share" style="margin-right: 10px;"></span>Views</h4>');
    if (projectData.views.length) {
	html.push('<ul class="list-group col-sm-8">');
	for (var i=0; i<projectData.views.length; i++) {
	    var li = window.location.origin+"/";
	    if (projectData.views[i]['public']) {
		li += 'public/'+projectData.user.login+'/'+projectData.views[i].name;
	    } else {
		li += 'private/'+projectData.user.login+'/'+projectData.views[i].name + '?code=' + projectData.views[i].token;
	    }
	    html.push('<li class="list-group-item"><a href="'+li+'" target=_blank>'+li+'</a></li>');
	}
	html.push('</ul><div style="clear: both;"></div>');
    } else {
	html.push('<p>This project is not shared</p>');
    }

    // project files
    html.push('<h4><span class="glyphicon glyphicon-list-alt" style="margin-right: 10px;"></span>Files</h4><ul class="list-group col-sm-8">');
    fields = Object.keys(projectData.files).sort();
    for (var i=0; i<fields.length; i++) {
	if (projectData.files[fields[i]]) {
	    html.push('<li class="list-group-item"><a href="#" onclick="saveAs(projectData.files.'+fields[i]+', \''+fields[i]+'.txt\')">'+fields[i]+'</a></li>');
	}
    }
    html.push('</ul>');
    
    html.push('</div></div>');

    document.getElementById('details').innerHTML = html.join("\n");
}

function deleteProject(project, user) {
    if (confirm("Really delete this project? This cannot be undone!")) {
	var formData = new FormData();
	formData.append('project', project);
	formData.append('user', user);
	$.ajax({
	    url : '/project',
	    type : 'DELETE',
	    data : formData,
	    processData: false,
	    contentType: false,
	    success : function(data) {
		if (data.status == 'ok') {
		    document.location.reload(true);
		} else {
		    toastr.error(data.message);
		}
	    }
	});
    }
}

function showUserDetails(login) {
    var user = null;
    for (var i=0; i<tables.user.response.length; i++) {
	if (tables.user.response[i]['login'] == login) {
	    user = tables.user.response[i];
	    break;
	}
    }
    if (! user) {
	toastr.error('data error selecting user');
    }

    var html = [];
    html.push('<div class="panel panel-primary"><div class="panel-heading"><h4 class="panel-title">User Details</h4></div><div class="panel-body">');

    html.push('<button class="btn btn-danger btn-sm" style="float: right; margin-left: 10px;" onclick="deleteUser(\''+user.login+'\')"><span class="glyphicon glyphicon-trash"></span></button>');
    html.push('<button class="btn btn-default btn-sm" style="float: right; margin-left: 10px;" onclick="impersonate(\''+user.login+'\')">impersonate</button>');
    html.push('<h4>'+user.firstname+' '+user.lastname+'</h4>');
    html.push('<table class="table col-sm-8">');
    html.push('<tr><th>login</th><td>'+user.login+'</td></tr>');
    html.push('<tr><th>email</th><td>'+user.email+'</td></tr>');
    html.push('<tr><th>affiliation</th><td>'+user.affiliation+'</td></tr>');
    html.push('<tr><th>entry date</th><td>'+user.date+'</td></tr>');
    html.push('<tr><th>last visit</th><td>'+user.visit+'</td></tr>');
    html.push('<tr><th>clearance</th><td><select id="clearanceSelect"><option'+(user.clearance=='admin' ? ' selected' : '')+'>admin</option><option'+(user.clearance=='user' ? ' selected' : '')+'>user</option></select><button class="btn btn-sm btn-default" style="margin-left: 10px;" onclick="changeClearance(\''+login+'\', this.previousSibling.options[this.previousSibling.selectedIndex].value);">change</button></td></tr>');
    html.push('<tr><th>current project</th><td>'+(user.project ? '<a href="#" onclick="getProjectDetails(\''+user.project+'\', \''+user.login+'\');">'+user.project+'</a>' : '- none -')+'</td></tr>');
    html.push('<tr><th># of projects</th><td>'+user.projects+'</td></tr>');
    html.push('<tr><th>registration ip</th><td>'+user.ip+'</td></tr>');
    html.push('</div></div>');

    document.getElementById('details').innerHTML = html.join("\n");
}

function deleteUser(user) {
    var retval = prompt("Really delete this user? Type 'CONFIRM' to proceed.")
    if (retval && retval == "CONFIRM") {
	var formData = new FormData();
	formData.append('user', user);
	$.ajax({
	    url : '/user',
	    type : 'DELETE',
	    data : formData,
	    processData: false,
	    contentType: false,
	    success : function(data) {
		if (data.status == 'ok') {
		    document.location.reload(true);
		} else {
		    toastr.error(data.message);
		}
	    }
	});
    }
}

function changeClearance(user, clearance) {
    var formData = new FormData();
    formData.append('user', user);
    formData.append('clearance', clearance);
    $.ajax({
	url : '/clearance',
	type : 'POST',
	data : formData,
	processData: false,
	contentType: false,
	success : function(data) {
	    if (data.status == 'ok') {
		document.location.reload(true);
	    } else {
		toastr.error(data.message);
	    }
	}
    });
}
