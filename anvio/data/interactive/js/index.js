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
    var html = '<div style="padding-top: 10px;">';
    html += '<button type="button" style="margin-right: 5px; float: right;" class="btn btn-danger btn-sm" title="log out" onclick="performLogout();"><span class="glyphicon glyphicon-white glyphicon-off" aria-hidden="true"></span></button>';
    html += '<button type="button" style="margin-right: 5px;" class="btn btn-default btn-sm" title="project management" onclick="window.location=\'user.html\';"><span class="glyphicon glyphicon-home" aria-hidden="true"></span></button>';
    html += '<img src="images/user.png" style="width: 20px; position:  relative; top: 1px; order-radius: 3px; margin-left: 10px;"><span style="margin-left: 10px; position: relative; top: 2px;">logged in as '+user.firstname+' '+user.lastname+' ('+user.login+')</span>';
    html += '</div>';
    document.getElementById('multiUser').innerHTML = html;
};
