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
function initContent () {

    // username and logout button
    var html = [];
    html.push('<div style="padding-top: 2px; float: right;">');
    html.push('<img src="images/user.png" style="width: 32px; border-radius: 3px; margin-right: 5px;">');
    html.push('<span>logged in as '+user.firstname+' '+user.lastname+' ('+user.login+')</span>');
    html.push('<button type="button" style="margin-left: 5px;" class="btn btn-danger btn-sm" title="log out" onclick="performLogout();"><span class="glyphicon glyphicon-white glyphicon-off" aria-hidden="true"></span></button><div style="font-size: 11px;"><a href="change.html">change password</a>'+(user.clearance=='admin'?'<a style="float: right;" href="admin.html">admin page</a>':'')+'</div>');
    html.push('</div>');

    // divider
    html.push('<div style="clear: both; margin-top: 50px; margin-bottom: 50px;"></div>');

    // get all projects the user has access to
    html.push('<h3>Your Projects</h3>');
    if (user.projects.length) {
	html.push('<ul class="list-group col-sm-8">');
	for (var i=0; i<user.projects.length; i++) {
	    html.push('<li class="list-group-item">');
	    html.push('<a href="#" onclick="setActiveProject(\''+i+'\');" title="view project">'+user.projects[i].name+'</a>');
	    html.push('<button type="button" style="margin-right: 5px; float: right; position: relative; bottom: 5px;" class="btn btn-danger btn-sm" title="delete project" onclick="deleteProject(\''+i+'\')"><span class="glyphicon glyphicon-trash" aria-hidden="true"></span></button>');
	    html.push('<button type="button" style="margin-right: 5px; float: right; position: relative; bottom: 5px;" class="btn btn-default btn-sm" title="add data" onclick="addDataToProject(\''+i+'\')"><span class="glyphicon glyphicon-floppy-open" aria-hidden="true"></span></button>');
	    html.push('<button type="button" style="margin-right: 5px; float: right; position: relative; bottom: 5px;" class="btn btn-default btn-sm" title="share project" onclick="selectedProject=\''+i+'\';$(\'#modShareProject\').modal(\'show\');"><span class="glyphicon glyphicon-share" aria-hidden="true"></span></button>');
	    html.push('<button type="button" style="margin-right: 5px; float: right; position: relative; bottom: 5px;" class="btn btn-default btn-sm" title="project settings" onclick="showProjectSettings(\''+i+'\');"><span class="glyphicon glyphicon-cog" aria-hidden="true"></span></button>');
	    
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
function showProjectSettings (index) {
    var project = user.projects[index];

    var html = "";
    if (project.views.length) {
	html = "<p>This project has been shared.</p>";
	for (var i=0; i<project.views.length; i++) {
	    var baseURL = window.location.origin + '/project?name='+project.views[i].name;
	    var code = project.views[i]['public'] ? '' : '&code='+project.views[i].token;
	    html += "<p style='margin-bottom: 20px;'>"+(project.views[i]['public'] ? 'public ' : '')+'share link: <a href="'+baseURL+code+'" target=_blank>'+baseURL+code+'</a><button type="button" style="margin-left: 5px; float: right; position: relative; bottom: 7px;" class="btn btn-danger btn-sm" title="remove view" onclick="removeProjectView(\''+index+'\', \''+i+'\');"><span class="glyphicon glyphicon-trash" aria-hidden="true"></span></button></p>';
	}
    } else {
	html = "<p>This project has not been shared</p>"
    }
    
    document.getElementById('projectSettingsContent').innerHTML = html;
    $('#modProjectSettings').modal('show');
    
};

function removeProjectView (pindex, vindex) {
    if (confirm("Really delete this share? This cannot be undone!")) {
	var formData = new FormData();
	formData.append('project', user.projects[pindex].name);
	formData.append('name', user.projects[pindex].views[vindex].name);
	$.ajax({
	    url : '/share',
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
};

function addDataToProject (index) {
    window.selectedProject = index;
    $('#modUploadAdditionalData').modal('show');
};

function setActiveProject(index) {
    $.removeCookie('anvioView', { path: '/' });
    var formData = new FormData();
    formData.append('project', user.projects[index].name);
    $.ajax({
	url : '/project',
	type : 'POST',
	data : formData,
	processData: false,
	contentType: false,
	success : function(data) {
	    if (data.status == 'ok') {
		window.location = 'index.html';
	    } else {
		toastr.error(data.message);
	    }
	}
    });
}

function deleteProject(index) {
    if (confirm("Really delete this project? This cannot be undone!")) {
	var formData = new FormData();
	formData.append('project', user.projects[index].name);
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

function shareProject() {
    var index = selectedProject;
    var name = document.getElementById('projectName').value;
    var isPublic = document.getElementById('projectPublic').checked;
    if (! name.match(/^\w+$/)) {
	alert('The project name may only contain word characters without spaces.');
    } else {
	$('#modShareProject').modal('hide');
	var formData = new FormData();
	formData.append('name', name);
	formData.append('project', user.projects[index].name);
	formData.append('public', isPublic ? 1 : 0);
	$.ajax({
	    url : '/share',
	    type : 'POST',
	    data : formData,
	    processData: false,
	    contentType: false,
	    complete : function(jqXHR) {
		var data = JSON.parse(jqXHR.responseText);
		if (data.status == 'error') {
		    toastr.error(data.message, 'share project failed');
		} else {
		    data = data.data;
		    var which = 0;
		    for (var i=0; i<user.projects.length; i++) {
			if (user.projects[i].name = data['project']) {
			    user.projects[i].views.push({ "name": data["name"], "token": data["token"], "public": data["public"] == "0" ? false : true });
			    which = i;
			    break;
			}
		    }
		    showProjectSettings(which);
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
    if ($('#samplesOrderFileSelect')[0].files.length) {
	formData.append('samplesOrderFile', $('#samplesOrderFileSelect')[0].files[0]);
    }
    if ($('#samplesInformationFileSelect')[0].files.length) {
	formData.append('samplesInformationFile', $('#samplesInformationFileSelect')[0].files[0]);
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

function uploadAdditional () {
    var index = window.selectedProject;
    var formData = new FormData();
    if ($('#additionalFileSelect')[0].files.length) {
	formData.append('additionalFile', $('#additionalFileSelect')[0].files[0]);
	formData.append('project', user.projects[index].name);
    } else {
	alert('You must provide a file');
	return;
    }

    if (! confirm("WARNING: If you already provided additional data it will be overwritten.\nDo you still want to proceed?")) {
	return false;
    }
    
    $('#modUploadAdditionalData').modal('hide');
    
    uploadProgress = document.createElement('div');
    uploadProgress.setAttribute('style', 'position: absolute; right: 0px; width: 400px; bottom: 0px; border: 1px solid lightgray; border-bottom: none; height: 25px; background-color: white;');
    uploadProgress.innerHTML = '<div class="progress" style="margin: 5px; margin-bottom: 0px;">\
  <div id="uploadProgressBar" class="progress-bar progress-bar-striped active" role="progressbar" aria-valuenow="0" aria-valuemin="0" aria-valuemax="100" style="min-width: 2em;">\
    0%\
  </div>\
</div>';
    document.body.appendChild(uploadProgress);
    
    $.ajax({
	url : '/uploadMore',
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
	    if (data.status == 'ok') {
		toastr.info(data.message);
	    } else {
		toastr.error(data.message);
	    }
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
