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
    $("#resetPasswordForm").submit(function(event){ event.preventDefault(); });
    //	  grecaptcha.render('recap', { 'sitekey' : '6Lf1FL4SAAAAAO3ToArzXm_cu6qvzIvZF4zviX2z' });
}

function resetPassword() {
    document.getElementById('submit').setAttribute('disabled', 'disabled');
    jQuery.post("/resetPassword", {
	"email": document.getElementById('inputEmail').value,
	//	"response": grecaptcha.getResponse()
    }, function (result) {
	document.getElementById('submit').removeAttribute('disabled');
	if (result.status == 'error') {
	    alert("Resetting password failed: "+result.message);
	} else {
	    document.getElementById('main').innerHTML = "<h3>Password Reset Successful</h3><div class='alert alert-success col-sm-6'><p>Your password has successfully been reset. You will receiv your new credentials at the registered email address shortly.</p></div><div style='clear: both;'></div><br><button class='btn' type='button' onclick='window.location=\"home.html\";'>go to home page</button>";
	}
    }).fail(function(result){
	document.getElementById('submit').removeAttribute('disabled');	  
	if (result.status == 'error') {
	    alert("Resetting password failed: "+result.message);
	} else {
	    alert('An error occurred during password reset');
	}
    });
};
