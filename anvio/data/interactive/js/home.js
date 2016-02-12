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

$(window).load(function(){
    carouselImages();
});

function initContent () {
    $('#password').on('keypress', function(e) {e=e||window.event;if(e.keyCode==13){performLogin();}});
    $("#loginForm").submit(function(event){ event.preventDefault(); });
    if (window.hasOwnProperty('user') && user) {
	window.location = 'user.html';
    }
};

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

function loginForm() {
    return '\
<form class="navbar-form navbar-right" role="login">\
  <div class="form-group">\
    <input type="text" class="form-control input-sm" placeholder="Login" style="width: 100px;" id="login">\
    <input type="password" class="form-control input-sm" placeholder="Password" style="width: 100px;" id="password">\
    <button type="button" class="btn btn-default btn-sm" onclick="performLogin();" id="loginButton">login</button>\
    <div style="clear: both; text-align: left; font-size: 11px;"><a style="" href="forgot.html">forgot password?</a></div>\
  </div>\
</form>\
<ul class="nav navbar-nav navbar-right">\
  <li><a href="register.html">Register</a></li>\
  <li><a href="http://merenlab.org/projects/anvio/" target=_blank>Help</a></li>\
</ul>';
}
