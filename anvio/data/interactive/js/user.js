/**
 *  User functions for anvi'o interactive interface
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

// user section on anvio user.html
function showUserData (user) {

    // username and logout button
    var html = [];
    html.push('<div style="padding-top: 2px; float: right;">');
    html.push('<img src="images/user.png" style="width: 32px; border-radius: 3px; margin-right: 5px;">');
    html.push('<span>logged in as '+user.firstname+' '+user.lastname+' ('+user.login+')</span>');
    html.push('<button type="button" style="margin-left: 5px;" class="btn btn-danger btn-sm" title="log out" onclick="performLogout();"><span class="glyphicon glyphicon-white glyphicon-off" aria-hidden="true"></span></button>');
    html.push('</div>');

    // divider
    html.push('<div style="clear: both; margin-top: 50px; margin-bottom: 50px;"></div>');

    // get all projects the user has access to
    html.push('<h3>Your Projects</h3>');
    if (user.project_names.length) {
	html.push('<ul class="list-group col-sm-6">');
	for (var i=0; i<user.project_names.length; i++) {
	    html.push('<li class="list-group-item">');
	    html.push('<a href="#" onclick="setActiveProject(\''+user.project_names[i]+'\');" title="view project">'+user.project_names[i]+'</a>');
	    html.push('<button type="button" style="margin-right: 5px; float: right; position: relative; bottom: 5px;" class="btn btn-default btn-sm" title="share project" onclick="$(\'#modShareProject\').modal(\'show\');"><span class="glyphicon glyphicon-share" aria-hidden="true"></span></button>');
	    html.push('<button type="button" style="margin-right: 5px; float: right; position: relative; bottom: 5px;" class="btn btn-default btn-sm" title="delete project" onclick="deleteProject(\''+user.project_names[i]+'\')"><span class="glyphicon glyphicon-trash" aria-hidden="true"></span></button>');
	    html.push('<button type="button" style="margin-right: 5px; float: right; position: relative; bottom: 5px;" class="btn btn-default btn-sm" title="add data" onclick="addDataToProject(\''+user.project_names[i]+'\')"><span class="glyphicon glyphicon-floppy-open" aria-hidden="true"></span></button>');

	    html.push('</button>');
	    html.push('</li>');
	}
	html.push('</ul>');
    } else {
	html.push('<p>You currently do not have any projects</p>');
    }

    // divider
    html.push('<div style="clear: both; margin-top: 15px;"></div>');
    
    // upload button
    html.push('<button type="button" style="margin-right: 5px;" class="btn btn-default btn-sm" title="upload data files" onclick="$(\'#modUploadData\').modal(\'show\');"><span class="glyphicon glyphicon-floppy-open" aria-hidden="true"></span> upload new project</button>');
    
    document.getElementById('content').innerHTML = html.join('');

}

/* Project Management */
function addDataToProject (project) {

};

function setActiveProject(p, view, token) {
    if (view) {
	window.location(window.location.origin+'/project?name=p'+(token ? '&code='+token : ''));
    }
    $.removeCookie('anvioView', { path: '/' });
    var formData = new FormData();
    formData.append('project', p);
    $.ajax({
	url : '/project',
	type : 'POST',
	data : formData,
	processData: false,
	contentType: false,
	success : function(data) {
	    window.location = 'index.html';
	}
    });
}

function deleteProject() {
    if (confirm("Really delete this project? This cannot be undone!")) {
	var formData = new FormData();
	formData.append('project', document.title);
	$.ajax({
	    url : '/project',
	    type : 'DELETE',
	    data : formData,
	    processData: false,
	    contentType: false,
	    success : function(data) {
		document.location.reload(true);
	    }
	});
    }
}

function shareProject() {
    var name = document.getElementById('projectName').value;
    var isPublic = document.getElementById('projectPublic').checked;
    if (! name.match(/^\w+$/)) {
	aler('The project name may only contain word characters');
    } else {
	$('#modShareProject').modal('hide');
	session.project.view = name;
	session.project.isPublic = isPublic;
	var formData = new FormData();
	formData.append('name', name);
	formData.append('public', isPublic ? 1 : 0);
	$.ajax({
	    url : '/share',
	    type : 'POST',
	    data : formData,
	    processData: false,
	    contentType: false,
	    complete : function(jqXHR) {
		var data = JSON.parse(jqXHR.responseText);
		if (data.hasOwnProperty('ERROR')) {
		    toastr.error(data.ERROR, 'share project failed');
		} else {
		    var baseURL = window.location.origin + '/project?name='+session.project.view;
		    var code = session.project.isPublic ? '' : '&code='+data.token;
		    var msg = 'Your project has been shared.<br>It is now available via the following link:<br><br><a href="'+baseURL+code+'" target=_blank>'+baseURL+code+'</a>';
		    document.getElementById('projectSettingsContent').innerHTML = msg;
		    $('#modProjectSettings').modal('show');
		}
	    }
	});
    }
}

/* Data upload */
function uploadFiles () {
    var formData = new FormData();
    if ($('#treeFileSelect')[0].files.length) {
	formData.append('treeFile', $('#treeFileSelect')[0].files[0]);
    } else {
	alert('You must provide a tree file');
	return;
    }
    if ($('#fastaFileSelect')[0].files.length) {
	formData.append('fastaFile', $('#fastaFileSelect')[0].files[0]);
    }
    if ($('#dataFileSelect')[0].files.length) {
	formData.append('dataFile', $('#dataFileSelect')[0].files[0]);
    }
    if ($('#uploadTitle')[0].value) {
	formData.append('title', $('#uploadTitle')[0].value);
    }
    $('#modUploadData').modal('hide');
    
    uploadProgress = document.createElement('div');
    uploadProgress.setAttribute('style', 'position: absolute; right: 0px; width: 400px; bottom: 0px; border: 1px solid lightgray; border-bottom: none; height: 25px; background-color: white;');
    uploadProgress.innerHTML = '<div class="progress" style="margin: 5px; margin-bottom: 0px;">\
  <div id="uploadProgressBar" class="progress-bar progress-bar-striped active" role="progressbar" aria-valuenow="0" aria-valuemin="0" aria-valuemax="100" style="min-width: 2em;">\
    0%\
  </div>\
</div>';
    document.body.appendChild(uploadProgress);
    
    $.ajax({
	url : '/upload',
	xhr: function() {
	    var xhr = new window.XMLHttpRequest();
	    xhr.upload.addEventListener("progress", dataFileUploadProgress, false);
	    return xhr;
	},
	type : 'POST',
	data : formData,
	processData: false,
	contentType: false,
	success : function(data) {
	    document.body.removeChild(uploadProgress);
	    uploadProgress = null;
            document.location.reload(true);
	}
    });
}

function dataFileUploadProgress (event) {
    var curr = parseInt(event.loaded / event.total * 100);
    var p = document.getElementById('uploadProgressBar');
    p.setAttribute('aria-valuenow', curr);
    p.style.width = curr + "%";
    p.innerHTML = curr + "%";
}

function uploadFileSelected (which) {
    $('#'+which+'FileName')[0].value = $('#'+which+'FileSelect')[0].files[0].name || "";
}
