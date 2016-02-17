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
	html += '<img src="http://www.gravatar.com/avatar/'+SparkMD5.hash(user.email.toLowerCase())+'?d=mm&s=30" style="margin-left: 10px; border: 1px solid #cccccc; border-radius: 3px; float: left; margin-right: 10px;" title="logged in as '+user.firstname+' '+user.lastname+' ('+user.login+')" />';
	html += '<div id="projectInfo" style="width: 320px; float: left;"></div>';
    } else {
	html += '<div id="projectInfo" style="margin-right: 10px;"><img src="images/loading.gif" style="width: 32px;"></div>';
    }
    html += '</div>';
    document.getElementById('multiUser').innerHTML = html;

    $.ajax({
	url : '/project',
	type : 'GET',
	processData: false,
	contentType: false,
	complete : function(jqXHR) {
	    var r = JSON.parse(jqXHR.responseText);
	    if (r.status == 'error') {
		toastr.error(r.message);
		return;
	    }
	    var data = r.data;
	    document.getElementById('projectInfo').innerHTML = '<a target="_blank" href="/data/project" style="float: right;" class="btn btn-default btn-sm" title="download project archive"><i class="glyphicon glyphicon-download"></i></a><span style="font-weight: bold; cursor: help;" id="projectDescription">'+data.name + '</span><br/><i>by '+data.user+'</i>';
	    $('#projectDescription').popover({"content":data.description, "title": data.name, "trigger": "hover"})
	}
    });
};
