/**
 *  functions for anvi'o interactive interface
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

function initContent () {
    var html = '<div style="height: 32px;">';
    if (user) {
	html += '<button type="button" style="margin-right: 5px; float: right;" class="btn btn-danger btn-sm" title="log out" onclick="performLogout();"><span class="glyphicon glyphicon-white glyphicon-off" aria-hidden="true"></span></button>';
	html += '<button type="button" style="margin-right: 5px; float: left;" class="btn btn-default btn-sm" title="project management" onclick="window.location=\'/\';"><span class="glyphicon glyphicon-home" aria-hidden="true"></span></button>';
	html += '<img src="http://www.gravatar.com/avatar/'+SparkMD5.hash(user.email.toLowerCase())+'?s=30&d=mm" style="margin-left: 10px; border: 1px solid #cccccc; border-radius: 3px; float: left; margin-right: 10px;" title="logged in as '+user.firstname+' '+user.lastname+' ('+user.login+')" />';
	html += '<div id="projectInfo" style="width: 330px; float: left;"><img src="images/loading.gif" style="width: 32px;"></div>';
    } else {
	html += '<div id="projectInfo" style="margin-right: 10px;"><img src="images/loading.gif" style="width: 32px;"></div>';
    }
    html += '</div>';
    document.getElementById('multiUser').innerHTML = html;

    // get the project base information
    var promise = $.Deferred();
    $.ajax({
	url : '/project',
	type : 'GET',
	processData: false,
	p: promise,
	contentType: false,
	success : function(data) {
	    if (data.status == 'error') {
		toastr.error(data.message);
		return;
	    }
	    window.pbasedata = data.data;

	    // add the content of the info section
	    document.getElementById('projectInfo').innerHTML = '<button style="float: right;" class="btn btn-info btn-sm" title="open project information" onclick="$(\'#modProjectInfo\').modal(\'show\');"><i class="glyphicon glyphicon-info-sign"></i></button><span style="font-weight: bold;">'+pbasedata.name + '</span><br/><i>by '+pbasedata.user+'</i>';
	    document.getElementById('sidebar').style.marginTop = "86px";
	    this.p.resolve();
	},
	error: function(jqXHR) {
	    document.getElementById('multiUser').style.display = 'none';
	    document.getElementById('sidebar').style.marginTop = "43px";
	}
    });

    promise.then( function () {
	// get the project files
	$.ajax({
	    url : '/projectfiles',
	    type : 'GET',
	    processData: false,
	    contentType: false,
	    success : function(data) {
		if (data.status == 'error') {
		    toastr.error(data.message);
		    return;
		}
		window.pdata = data.data;
		
		// create the HTML for the modal
		var html = [];
		
		// project base data
		html.push('<div style="padding: 10px;">');
		html.push('<p style="margin-top: -10px; margin-bottom: 0px; font-size: 22px; font-weight: bold;">' + pdata.name + '</p>');
		html.push('<hr style="margin: 0px; margin-top: 5px;" />');
		html.push('<p style="font-family: \'PT Serif\',serif; font-size: 14px;"><i>by '+ pdata.user + ' (hash: ' + pdata.path + ')</i></p>');
		html.push('<div class="desc">'+pdata.description+'</div>');
		// project files
		html.push('<p style="font-family: \'PT Serif\',serif; font-size: 16px; font-weight: bold;">Project Files<button class="btn btn-default" style="float: right; position: relative; bottom: 5px;" title="download all project files" onclick="downloadProjectZIP();"><i class="glyphicon glyphicon-floppy-save"></i></button></p>');
		var fields = Object.keys(pdata.files).sort();
		for (var i=0; i<fields.length; i++) {
		    if (pdata.files[fields[i]]) {
			html.push('<li style="cursor: pointer;font-family: \'PT Serif\',serif; font-size: 14px;" class="list-group-item" title="download this file" onclick="saveAs(pdata.files.'+fields[i]+', \''+fields[i]+'.txt\')">'+fields[i]+'<i style="float: right; margin-right: 10px;" class="glyphicon glyphicon-floppy-save"></i></li>');
		    }
		}
		html.push('</div>');
		
		document.getElementById('projectInfoContent').innerHTML = html.join("\n");
	    },
	    error: function(jqXHR) {
		document.getElementById('multiUser').style.display = 'none';
		document.getElementById('sidebar').style.marginTop = "43px";
	    }
	});
    });
};

function downloadProjectZIP() {
  window.location = '/app/downloadProjectFiles';
};
