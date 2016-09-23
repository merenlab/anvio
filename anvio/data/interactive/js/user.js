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
    var formData = new FormData();
    formData.append('token', user.token);
    formData.append('verbose', true);
    $.ajax({
        url : '/token',
    processData: false,
    contentType: false,
    type : 'POST',
        data : formData,
        success : function(data) {
        if (data.status == 'ok') {
        window.user = data.data;
        showUserdata();
        }
        }
    });
}

function noUser() {
    window.location = 'home.html';
}

function showUserdata () {
    // username and logout button
    var html = [];
    html.push('<div style="padding-top: 10px; float: right;">');
    html.push('<img src="http://www.gravatar.com/avatar/'+SparkMD5.hash(user.email.toLowerCase())+'?s=30&d=mm" style="margin-right: 5px;" />');
    html.push('<span>logged in as '+user.firstname+' '+user.lastname+' ('+user.login+')</span>');
    html.push('<button type="button" style="margin-left: 5px;" class="btn btn-danger btn-sm" title="log out" onclick="performLogout();"><span class="glyphicon glyphicon-white glyphicon-off" aria-hidden="true"></span></button><div style="font-size: 11px;"><a href="change.html">change password</a>'+(user.clearance=='admin'?'<a style="float: right;" href="admin.html">admin page</a>':'')+'</div>');
    html.push('</div>');

    // divider
    html.push('<div style="clear: both; margin-top: 50px; margin-bottom: 50px;"></div>');

    // get all projects the user has access to
    html.push('<h3>Your Projects</h3>');
    if (user.projects && user.projects.length) {
    html.push('<ul class="list-group col-sm-8">');
    for (var i=0; i<user.projects.length; i++) {
        var isPublic = false;
        var isShared = false;
        var publicLink = "";
        for (var h=0; h<user.projects[i].views.length; h++) {
        isShared = true;
        if (user.projects[i].views[h]['public'] == 1) {
            isPublic = true;
            publicLink = window.location.origin + '/public/' +user.projects[i].views[h]['user'] + '/' + user.projects[i].views[h].name;
        }
        }
        
        html.push('<li class="list-group-item">');
        html.push('<a href="#" onclick="setActiveProject(\''+i+'\');" title="'+user.projects[i].description+'">'+user.projects[i].name+'</a>');
        
        html.push('<button type="button" style="margin-right: 5px; float: right; position: relative; bottom: 5px;" class="btn btn-danger btn-sm" title="delete project" onclick="deleteProject(\''+i+'\')" id="projectDelete'+i+'"><span class="glyphicon glyphicon-trash" aria-hidden="true"></span></button>');
        html.push('<button type="button" style="margin-right: 5px; float: right; position: relative; bottom: 5px;" class="btn btn-default btn-sm" title="add data" onclick="addDataToProject(\''+i+'\')"><span class="glyphicon glyphicon-floppy-open" aria-hidden="true"></span></button>');

        html.push('<button type="button" style="margin-right: 5px; float: right; position: relative; bottom: 5px;" class="btn btn-default btn-sm" title="project settings" onclick="showProjectSettings(\''+i+'\');"><span class="glyphicon glyphicon-cog" aria-hidden="true"></span></button>');
        
        html.push('<button type="button" style="margin-right: 5px; float: right; position: relative; bottom: 5px;" class="btn btn-default btn-sm'+(isPublic ? ' btn-success' : (isShared ? ' btn-info': ''))+'" title="'+(isPublic ? ' publicly shared, add share?' : (isShared ? ' privately shared, add share?': 'share'))+'" onclick="selectedProject=\''+i+'\';$(\'#modShareProject\').modal(\'show\');"><span class="glyphicon glyphicon-share" aria-hidden="true"></span></button>');
        
        if (isPublic) {
        html.push('<a href="http://twitter.com/home/?status=my anvi\'o project: '+publicLink+'" target=_blank type="button" style="margin-right: 5px; float: right; position: relative; bottom: 5px; padding: 0px;" class="btn btn-default btn-sm" title="tweet on twitter"><img src="images/twitter.png" style="width: 28px;"></a>');
        html.push('<a href="http://www.facebook.com/sharer/sharer.php?u='+publicLink+'" target=_blank type="button" style="margin-right: 5px; float: right; position: relative; bottom: 5px; padding: 0px;" class="btn btn-default btn-sm" title="share on facebook"><img src="images/facebook.png" style="width: 28px;"></a>');
        }
        
        html.push('</li>');
    }
    html.push('</ul>');
    } else {
    html.push('<p>You currently do not have any projects</p>');
    }

    // divider
    html.push('<div style="clear: both; margin-top: 15px;"></div>');
    
    // upload button
    html.push('<button type="button" style="margin-right: 5px;" class="btn btn-default btn-sm" title="upload data files" onclick="$(\'#modUploadData\').modal(\'show\');"><span class="glyphicon glyphicon-floppy-open" aria-hidden="true"></span> New Project</button>');
    
    document.getElementById('content').innerHTML = html.join('');

}

/* Project Management */
function showProjectSettings (index) {
    var project = user.projects[index];

    var html = "<h4 style='margin-top: 0px;'>"+project.name+"</h4><p><button class='btn btn-default btn-sm' title='edit description' style='float: right;' onclick='editProjectDescription("+index+");' id='projectDescriptionButton'><span class='glyphicon glyphicon-pencil' aria-hidden='true'></span></button><div id='projectDescription'>"+project.description+"</div></p>";
    html += "<h4 style='margin-top: 0px;'>Statistics</h4>";
    var k = Object.keys(project.metadata);
    if (k.length) {
    html += "<table class='table table-condensed'>";
    k = k.sort();
    for (var i=0; i<k.length; i++) {
        html += "<tr><th>"+k[i]+"</th><td>"+project.metadata[k[i]]+"</td></tr>";
    }
    html += "</table>";
    } else {
    html += "<p>No statistics available</p>";
    }
    html += "<h4 style='margin-top: 0px;'>Sharing</h4>";
    if (project.views.length) {
    html += "<p>This project has been shared.</p>";
    for (var i=0; i<project.views.length; i++) {
        var baseURL = window.location.origin + '/' + (project.views[i]['public'] ? 'public' : 'private') + '/'+project.views[i]['user'] + '/' + project.views[i].name;
        var code = project.views[i]['public'] ? '' : '?code='+project.views[i].token;
        html += "<p style='margin-bottom: 20px;'>"+(project.views[i]['public'] ? 'public ' : '')+'share link: <a href="'+baseURL+code+'" target=_blank>'+baseURL+code+'</a><button type="button" style="margin-left: 5px; float: right; position: relative; bottom: 7px;" class="btn btn-danger btn-sm" title="remove view" onclick="removeProjectView(\''+index+'\', \''+i+'\');"><span class="glyphicon glyphicon-trash" aria-hidden="true"></span></button></p>';
    }
    } else {
    html += "<p>This project has not been shared</p>"
    }
    html += "<h4 style='margin-top: 0px'>Files</h4><table class='table'><tr><th>type</th><th>size</th><th>created</th><th>view</th></tr>";
    var fnames = Object.keys(project.files).sort();
    for (var i=0; i<fnames.length; i++) {
    html += "<tr><td>"+fnames[i]+"</td><td>"+project.files[fnames[i]].size.byteSize()+"</td><td>"+project.files[fnames[i]].created+"</td><td><button class='btn btn-sm btn-default' onclick='$(\"#pfc"+i+"\").toggle();'><span class='glyphicon glyphicon-list-alt'></span></button></td></tr><tr style='display: none;' id='pfc"+i+"'><td colspan=4><pre style='white-space: pre-wrap'>"+project.files[fnames[i]].top+(project.files[fnames[i]].size > 1024 ? "..." : "")+"</pre></td></tr>";
    }
    html += "</table>";
    
    document.getElementById('projectSettingsContent').innerHTML = html;
    $('#modProjectSettings').modal('show');
    
};

function editProjectDescription (index) {
    var project = user.projects[index];
    var btn = document.getElementById('projectDescriptionButton');
    if (btn.getAttribute('title') == 'edit description') {
    btn.setAttribute('title', 'save changes');
    btn.innerHTML = "<span class='glyphicon glyphicon-ok' aria-hidden='true'></span>";
    document.getElementById('projectDescription').innerHTML = "<textarea id='newProjectDescription' style='height: 100px; width: 500px;'>"+project.description+"</textarea>";
    } else {
    btn.setAttribute('disabled', 'disabled');
    $('#modProjectSettings').modal('hide');
    $.ajax({
        url : '/project',
        type : 'PUT',
        data : '{ "project": "'+project.name.replace(/"/g, "\\\"")+'", "description": "'+document.getElementById('newProjectDescription').value.replace(/"/g, "\\\"").replace(/[\x00-\x1F\x7F-\x9F]/g, " ")+'"}',
        processData: false,
        contentType: false,
        success : function(data) {
        if (data.status == 'ok') {
            alert('description updated');
            document.location.reload(true);
        } else {
            toastr.error(data.message);
        }
        }
    });
    }
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
    var html = '';
    var fileNames = ['additionalFile', 'dataFile', 'treeFile', 'fastaFile', 'samplesOrderFile', 'samplesInformationFile'];
    for (var i=0; i<fileNames.length; i++) {
    html += '<option value="'+fileNames[i]+'">'+ fileNames[i];
    if (user.projects[index].files[fileNames[i]]) {
        html += ' (replace version from '+user.projects[index].files[fileNames[i]].created+')';
    }
    html += '</option>';
    }
    document.getElementById('additionalFileType').innerHTML = html;
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
    toastr.error('The view name may only contain word characters without spaces.');
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
            if (user.projects[i].name == data['project']) {
                user.projects[i].views.push({ "name": data["name"], "user": data["user"], "token": data["token"], "public": data["public"] == "0" ? false : true });
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
        toastr.error('You must provide a tree file');
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
    } else {
        toastr.error("You must provide a project name");
        return;
    }
    if ($('#uploadDescription')[0].value) {
        formData.append('description', $('#uploadDescription')[0].value);
    }
    $('#modUploadData').modal('hide');
    
    uploadProgress = document.createElement('div');
    uploadProgress.setAttribute('style', 'position: absolute; right: 0px; width: 400px; bottom: 0px; border: 1px solid lightgray; border-bottom: none; height: 25px; background-color: white;');
    uploadProgress.innerHTML = '<div class="progress" style="margin: 5px; margin-bottom: 0px;"><div id="uploadProgressBar" class="progress-bar progress-bar-striped active" role="progressbar" aria-valuenow="0" aria-valuemin="0" aria-valuemax="100" style="min-width: 2em;">0%</div></div>';
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
        if (data.hasOwnProperty('status') && data.status == 'ok') {
            toastr.info(data.message);
            document.location.reload(true);
        } else {
            toastr.error(data.message, null, {timeOut: 0, closeButton: true});
        }
    }
    });
}

function uploadAdditional () {
    var index = window.selectedProject;
    var formData = new FormData();
    if ($('#additionalFileSelect')[0].files.length) {
    formData.append('uploadFile', $('#additionalFileSelect')[0].files[0]);
    formData.append('project', user.projects[index].name);
    formData.append('type', $('#additionalFileType')[0].options[$('#additionalFileType')[0].selectedIndex].value);
    } else {
        alert('You must provide a file');
        return;
    }
    
    $('#modUploadAdditionalData').modal('hide');
    
    uploadProgress = document.createElement('div');
    uploadProgress.setAttribute('style', 'position: absolute; right: 0px; width: 400px; bottom: 0px; border: 1px solid lightgray; border-bottom: none; height: 25px; background-color: white;');
    uploadProgress.innerHTML = '<div class="progress" style="margin: 5px; margin-bottom: 0px;"><div id="uploadProgressBar" class="progress-bar progress-bar-striped active" role="progressbar" aria-valuenow="0" aria-valuemin="0" aria-valuemax="100" style="min-width: 2em;">0%</div></div>';
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
            document.location.reload(true);
        } else {
            toastr.error(data.message, null, {timeOut: 0, closeButton: true});
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
    // validate file contents
    var fileReader = new FileReader();
    var file = $('#'+which+'FileSelect')[0].files[0];

    if (file.size > 1500000)
    {
        toastr.error('Sorry, current version does not support files larger than 1.5 MB.');

        $('#'+which+'FileName').val("")
        $('#'+which+'FileSelect').val("");
        return;
    }

    fileReader.onerror = function (error) {
        toastr.error('the file could not be opened');
    }
    fileReader.onload = function(e) {
        var fileContent = e.target.result;
        if (which == 'tree') {
            if (isNewickTree(fileContent)) {
                toastr.info('valid newick tree file');
                $('#'+which+'FileName')[0].value = file.name;
            } else {
                toastr.error('file is not in newick tree format');
                $('#'+which+'FileName').val("")
                $('#'+which+'FileSelect').val("");
            }
        } else if (which == 'data') {
            if (isDataFile(fileContent)) {
                toastr.info('valid data file');
                $('#'+which+'FileName')[0].value = file.name;
            } else {
                toastr.error('invalid data file');
                $('#'+which+'FileName').val("")
                $('#'+which+'FileSelect').val("");
            }
        } else if (which == 'fasta') {
            if (isFastaFile(fileContent)) {
                toastr.info('valid FASTA file');
                $('#'+which+'FileName')[0].value = file.name;
            } else {
                toastr.error('invalid FASTA file');
                $('#'+which+'FileName').val("")
                $('#'+which+'FileSelect').val("");
            }
        } else if (which == 'samplesOrder') {
            if (isSamplesOrderFile(fileContent)) {
                toastr.info('valid samples order file');
                $('#'+which+'FileName')[0].value = file.name;
            } else {
                toastr.error('invalid samples order file');
                $('#'+which+'FileName').val("")
                $('#'+which+'FileSelect').val("");
            }
        } else if (which == 'samplesInformation') {
            if (isSamplesInfoFile(fileContent)) {
                toastr.info('valid samples information file');
                $('#'+which+'FileName')[0].value = file.name;
            } else {
                toastr.error('invalid samples information file');
                $('#'+which+'FileName').val("")
                $('#'+which+'FileSelect').val("");
            }
        } else {
            $('#'+which+'FileName')[0].value = file.name;
        }
    }
    fileReader.readAsText(file);
}

Number.prototype.byteSize = function() {
    var size = this;
    var magnitude = "B";

    if (size > 999) {
        size = size / 1024;
        magnitude = "KB";
    }
    if (size > 999) {
        size = size / 1024;
        magnitude = "MB";
    }
    if (size > 999) {
        size = size / 1024;
        magnitude = "GB";
    }
    if (size > 999) {
        size = size / 1024;
        magnitude = "TB";
    }
    if (size > 999) {
        size = size / 1024;
        magnitude = "PB";
    }

    size = size.toFixed(1);

    size += '';
    var x = size.split('.');
    var x1 = x[0];
    var x2 = x.length > 1 ? '.' + x[1] : '';
    var rgx = /(\d+)(\d{3})/;
    while (rgx.test(x1)) {
        x1 = x1.replace(rgx, '$1' + ',' + '$2');
    }
    size =  x1 + x2;

    return size + " " + magnitude;
};

function isSamplesInfoFile (string) {
    try {
    var lines = string.split(/\n/);
    if (! lines.length) {
        lines = string.split(/\r/);
    }
    if (! lines.length) {
        toastr.error('file contains no data', null, {timeOut: 0, closeButton: true});
        return false;
    }
    if (lines.length < 2) {
        toastr.error('file contains only one line', null, {timeOut: 0, closeButton: true});
        return false;
    }
    var header = lines[0].split(/\t/);
    if (header.length < 2) {
        toastr.error('sample info file must have at least two columns', null, {timeOut: 0, closeButton: true});
        return false;
    }
    if (header[0] !== 'samples') {
        toastr.error('the cell A1 in the file must contain the string "samples"', null, {timeOut: 0, closeButton: true});
        return false;
    }
    for (var i=1; i<lines.length; i++) {
        if (i + 1 == lines.length && lines[i].length == 0) {
        continue;
        }
        var cells = lines[i].split(/\t/);
        if (cells.length != header.length) {
        toastr.error('invalid number of columns in line '+(i+1)+' ('+cells.length+'), should be '+header.length, null, {timeOut: 0, closeButton: true});
        return false;
        }
    }
    return true;
    } catch (error) {
    toastr.error('could not read file', null, {timeOut: 0, closeButton: true});
    return false;
    }
}

function isSamplesOrderFile (string) {
    try {
    var lines = string.split(/\n/);
    if (! lines.length) {
        lines = string.split(/\r/);
    }
    if (! lines.length) {
        toastr.error('file contains no data', null, {timeOut: 0, closeButton: true});
        return false;
    }
    if (lines.length < 2) {
        toastr.error('file contains only one line', null, {timeOut: 0, closeButton: true});
        return false;
    }
    var validHeaders = { 'attributes': true, 'basic': true, 'newick': true };
    var headerPositions = {};
    var headers = lines[0].split(/\t/);
    if (headers.length > 3 || headers.length < 2) {
        toastr.error('header column must contain at least two and at most three columns', null, {timeOut: 0, closeButton: true});
        return false;
    }
    var hasAttributesColumn = false;
    for (var i=0; i<headers.length; i++) {
        if (! validHeaders[headers[i]]) {
        toastr.error('header column contains invalid term: '+headers[i], null, {timeOut: 0, closeButton: true});
        return false;
        } else {
        if (headers[i] == 'attributes') {
            hasAttributesColumn = true;
        }
        headerPositions[i] = headers[i];
        }
    }
    if (! hasAttributesColumn) {
        toastr.error('header must contain an attributes column', null, {timeOut: 0, closeButton: true});
        return false;
    }
    for (var i=1; i<lines.length; i++) {
        if (i + 1 == lines.length && lines[i].length == 0) {
        continue;
        }
        var cells = lines[i].split(/\t/);
        var containsInfo = false;
        for (var h=0; h<cells.length; h++) {
        if (headerPositions[h] == 'attributes') {
            if (! cells[h].length) {
            toastr.error('row '+(i+1)+' does not contain an attribute name', null, {timeOut: 0, closeButton: true});
            return false;
            }
        } else if (headerPositions[h] == 'basic') {
            if (cells[h].length) {
            var basics = cells[h].split(/,/);
            if (basics.length < 2) {
                toastr.error('basic cell in line '+(i+1)+' contains less than two items', null, {timeOut: 0, closeButton: true});
                return false;
            }
            containsInfo = true;
            }
        } else if (headerPositions[h] == 'newick') {
            if (cells[h].length) {
            if (! isNewickTree(cells[h])) {
                toastr.error('newick cell in line '+(i+1)+' is not valid newick format', null, {timeOut: 0, closeButton: true});
                return false;
            }
            containsInfo = true;
            }
        }
        }
        if (! containsInfo) {
        toastr.error('row '+(i+i)+' does not contain data', null, {timeOut: 0, closeButton: true});
        return false;
        }
    }
    return true;
    } catch (error) {
    toastr.error('could not read file', null, {timeOut: 0, closeButton: true});
    return false;
    }
}

function isDataFile (string) {
    try {
    var lines = string.split(/\n/);
    if (! lines.length) {
        lines = string.split(/\r/);
    }
    if (! lines.length) {
        toastr.error('file contains no data', null, {timeOut: 0, closeButton: true});
        return false;
    }
    var header = lines[0].split(/\t/);
    if (header.length < 2) {
        toastr.error('data files must have at least two columns', null, {timeOut: 0, closeButton: true});
        return false;
    }
    if (header[0] !== 'contig') {
        toastr.error('the cell A1 in the file must contain the string "contig"', null, {timeOut: 0, closeButton: true});
        return false;
    }
    for (var i=1; i<lines.length; i++) {
        if (i + 1 == lines.length && lines[i].length == 0) {
        continue;
        }
        var cells = lines[i].split(/\t/);
        if (cells.length != header.length) {
        toastr.error('invalid number of columns in line '+(i+1)+' ('+cells.length+'), should be '+header.length, null, {timeOut: 0, closeButton: true});
        return false;
        }
    }
    return true;
    } catch (error) {
    toastr.error('could not read file', null, {timeOut: 0, closeButton: true});
    return false;
    }
}

function isFastaFile (string) {
    try {
    var lines = string.split(/\n/);
    if (! lines.length) {
        lines = string.split(/\r/);
    }
    if (! lines.length) {
        toastr.error('file contains no data', null, {timeOut: 0, closeButton: true});
        return false;
    }
    for (var i=0; i<lines.length; i++) {
        if (i==0) {
        if (! lines[i].match(/^>[\w\d\s]+$/)) {
            toastr.error('no identifier in line 1', null, {timeOut: 0, closeButton: true});
            return false;
        }
        } else {
        if (! (lines[i].match(/^[acgturykmswbdhvnxACGTURYKMSWBDHVNX\*\-\n\r\s]*$/) || lines[i].match(/^>[\w\d\s]+$/))) {
            toastr.error('invalid FASTA in line '+(i+1)+':\n"'+lines[i]+'"', null, {timeOut: 0, closeButton: true});
            return false;
        }
        if (lines[i].match(/^>[\w\d\s]+$/) && lines[i - 1].match(/^>[\w\d\s]+$/)) {
            toastr.error('no sequence for ID '+lines[i-1]+' in line '+(i+1), null, {timeOut: 0, closeButton: true});
            return false;
        }
        }
    }
    return true;
    } catch (error) {
    toastr.error('could not read file', null, {timeOut: 0, closeButton: true});
    return false;
    }
}

function isNewickTree (string) {
    var parents = [];
    var tree = {};
    try {
    if (string.match(/^\(/)) {
        var items = string.split(/\s*(;|\(|\)|,|:)\s*/);
        for (var i=0; i<items.length; i++) {
        switch (items[i]) {
        case '(':
            var subtree = {};
            tree.children = [subtree];
            parents.push(tree);
            tree = subtree;
            break;
        case ',':
            var subtree = {};
            parents[parents.length-1].children.push(subtree);
            tree = subtree;
            break;
        case ')':
            tree = parents.pop();
            break;
        case ':':
            break;
        default:
            var x = items[i-1];
            if (x == ')' || x == '(' || x == ',') {
            tree.name = items[i];
            } else if (x == ':') {
            tree.value = parseFloat(items[i]);
            }
        }
        }
    }
    } catch (error) {
    toastr.error('could not read file', null, {timeOut: 0, closeButton: true});
    return false;
    }
    if (tree.hasOwnProperty('name') && tree.hasOwnProperty('children')) {
    return true;
    } else {
    return false;
    }
};

function isNumber(n) {
  return !isNaN(parseFloat(n)) && isFinite(n);
}

