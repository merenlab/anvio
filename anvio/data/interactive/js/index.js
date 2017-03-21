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

function initProjectInfoBar() {
    $.ajax({
        type: 'GET',
        cache: false,
        url: '/project?timestamp=' + new Date().getTime(),
        success: function(data) {
            var html = '<div style="height: 32px;">';
            html += '<button type="button" style="margin-right: 5px; float: left;" class="btn btn-default btn-sm" title="project management" onclick="top.location.href=\'/\';"><span class="glyphicon glyphicon-home" aria-hidden="true"></span></button>';
            html += ' <a href="/' + data['username'] + '" target="_top"><img src="http://www.gravatar.com/avatar/' + data['user_email_hash'] + '?s=30&d=mm" style="margin-left: 10px; border: 1px solid #cccccc; border-radius: 3px; float: left; margin-right: 10px;" title="" /></a>';
            html += '<div id="projectInfo" style="width: 330px; float: left;"><b>' + data['project_name'] + '</b><br><i>by <a href="/' + data['username'] + '" target="_top">' + data['username'] + '</a></i></div>';
            html += '</div>';

            $('#multiUser').show();
            $('#multiUser').html(html);
            $('#sidebar').addClass('sidebarPadding');
        }
    });
}