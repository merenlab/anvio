/**
 *  Session functions for anvi'o interactive interface
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

/* login */
function performLogin () {
    var formData = new FormData();
    if (! (document.getElementById('login').value && document.getElementById('password').value)) {
	alert('you must enter both password and login');
	return;
    }
    formData.append('login', document.getElementById('login').value);
    formData.append('password', document.getElementById('password').value);
    $.ajax({
    	url : '/login',
	processData: false,
	contentType: false,
	type : 'POST',
    	data : formData,
    	success : function(data) {
	    if (data[0]) {
		window.location = 'user.html';
	    } else {
		alert("login failed");
	    }
    	}
    });
}

function performLogout () {
    var formData = new FormData();
    formData.append('login', session.user.login);
    $.ajax({
    	url : '/logout',
	processData: false,
	contentType: false,
	type : 'POST',
    	data : formData,
    	success : function(data) {
	    $.removeCookie('anvioSession', { path: '/' });
	    window.location = 'home.html';
	}
    });
}

function checkCookie () {
    var cookie = $.cookie('anvioSession');
    if (cookie) {
	var formData = new FormData();
	formData.append('token', cookie);
	$.ajax({
    	    url : '/token',
	    processData: false,
	    contentType: false,
	    type : 'POST',
    	    data : formData,
    	    success : function(data) {
		if (data[0]) {
		    session = { user: data[1] };
		    if (document.title == 'anvio home') {
			window.location = 'user.html';
		    } else if (document.title == 'anvio user home') {
			showUserData(data[1]);
		    } else if (document.title == 'anvio change password') {
			window.user = data[1];
		    } else {
			showUserOptions(data[1]);
		    }
		}
    	    }
	});
    }
}

// user section on avio index.html
function showUserOptions (user) {
    var html = '<div style="padding-top: 2px;">';
    html += '<button type="button" style="margin-right: 5px; float: right;" class="btn btn-danger btn-sm" title="log out" onclick="performLogout();"><span class="glyphicon glyphicon-white glyphicon-off" aria-hidden="true"></span></button>';
    html += '<button type="button" style="margin-right: 5px;" class="btn btn-default btn-sm" title="project management" onclick="window.location=\'user.html\';"><span class="glyphicon glyphicon-home" aria-hidden="true"></span></button>';
    html += '<img src="images/user.png" style="width: 20px; position:  relative; top: 1px; order-radius: 3px; margin-left: 10px;"><span style="margin-left: 10px; position: relative; top: 2px;">logged in as '+user.firstname+' '+user.lastname+' ('+session.user.login+')</span>';
    html += '</div>';
    document.getElementById('multiUser').innerHTML = html;
};

!function(e){"function"==typeof define&&define.amd?define(["jquery"],e):e(jQuery)}(function(e){function n(e){return u.raw?e:encodeURIComponent(e)}function o(e){return u.raw?e:decodeURIComponent(e)}function i(e){return n(u.json?JSON.stringify(e):String(e))}function r(e){0===e.indexOf('"')&&(e=e.slice(1,-1).replace(/\\"/g,'"').replace(/\\\\/g,"\\"));try{return e=decodeURIComponent(e.replace(c," ")),u.json?JSON.parse(e):e}catch(n){}}function t(n,o){var i=u.raw?n:r(n);return e.isFunction(o)?o(i):i}var c=/\+/g,u=e.cookie=function(r,c,a){if(void 0!==c&&!e.isFunction(c)){if(a=e.extend({},u.defaults,a),"number"==typeof a.expires){var d=a.expires,f=a.expires=new Date;f.setTime(+f+864e5*d)}return document.cookie=[n(r),"=",i(c),a.expires?"; expires="+a.expires.toUTCString():"",a.path?"; path="+a.path:"",a.domain?"; domain="+a.domain:"",a.secure?"; secure":""].join("")}for(var s=r?void 0:{},p=document.cookie?document.cookie.split("; "):[],m=0,v=p.length;v>m;m++){var x=p[m].split("="),k=o(x.shift()),l=x.join("=");if(r&&r===k){s=t(l,c);break}r||void 0===(l=t(l))||(s[k]=l)}return s};u.defaults={},e.removeCookie=function(n,o){return void 0===e.cookie(n)?!1:(e.cookie(n,"",e.extend({},o,{expires:-1})),!e.cookie(n))}});

// pixelpushig
function carouselImages() {
    $('.carouselImage').each(function(index){
	var factor = $(this)[0].naturalWidth / $(this)[0].naturalHeight;
	var scaleX = 350 / $(this)[0].naturalWidth;
	var scaleY = 350 / $(this)[0].naturalHeight;
	if (factor < 1) {
	    $(this)[0].height = $(this)[0].naturalHeight * scaleY;
	    $(this)[0].width = $(this)[0].naturalWidth * scaleY;
	    if (scaleY < 1) {
		$(this).css('marginLeft', ((350 - $(this)[0].width) / 2) + 'px');
		$(this).css('marginRight', ((350 - $(this)[0].width) / 2) + 'px');
	    }
	} else {
	    $(this)[0].height = $(this)[0].naturalHeight * scaleX;
	    $(this)[0].width = $(this)[0].naturalWidth * scaleX;
	    if (scaleX < 1) {
		$(this).css('marginTop', ((350 - $(this)[0].height) / 2) + 'px');
		$(this).css('marginBottom', ((350 - $(this)[0].height) / 2) + 'px');
	    }
	}
    });
}
