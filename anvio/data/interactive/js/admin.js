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
    if (user.clearance == 'admin') {
	html.push('<h3>Welcome back, master '+user.lastname+'</h3><p>I am awaiting your command.</p>');
    } else {
	html.push('<div class="alert alert-danger col-sm-6" role="alert">You are not authorized to view this page.</div>');
    }
    content.innerHTML = html.join('\n');
}

function showUserTable() {
    alert('hello world');
}
