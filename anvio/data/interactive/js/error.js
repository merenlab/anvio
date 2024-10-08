// * This file is part of anvi'o (<https://github.com/meren/anvio>).
// *
// * Anvi'o is a free software. You can redistribute this program
// * and/or modify it under the terms of the GNU General Public
// * License as published by the Free Software Foundation, either
// * version 3 of the License, or (at your option) any later version.
// *
// * You should have received a copy of the GNU General Public License
// * along with anvi'o. If not, see <http://opensource.org/licenses/GPL-3.0>.
// *
// * @license GPL-3.0+ <http://opensource.org/licenses/GPL-3.0>
// */

let ERROR_COUNT = 0

const issueCategories = [
    {
        'category' : 'Dependencies failed to load',
        'content' : `This can occur when npm packages that Anvi'o relies on fail to load. 
        Try running <b> npm install </b> in the <b>anvio/data/interactive</b> directory.`
    },
    {
        'category' : 'ReferenceError - ___ is not defined',
        'content' : `This can occur when Anvi'o wants to utilize some variable which it cannot resolve. Until we get some better troubleshooting advice here,
        try running your interactive session with the <b>--debug</b> flag`
    },
]

function alertDependencyError(dependencyError, isFinalDependency){
    if(dependencyError){
        ERROR_COUNT += 1
    }
    if(isFinalDependency && ERROR_COUNT){ // hacky way of 'iterating' all dependency calls before error messaging
        displayAlert('dependencies')
    }
}

function displayAlert(error){
    let reason;

    if(error == 'dependencies'){
        reason = 'loading dependencies'
    } else if (error.includes('is not defined')){
        reason = 'a variable reference error'
    }

    alert(`Anvi'o has encountered an error, possibly related to ${reason}. Anvi'o would like to offer some guidance in a new browser tab. Please make sure popups are enabled :)`)
    window.open('error-landing.html', '_blank')
}

function errorLandingContext(){ // onload function called by error-landing.html, generate help 'docs' from object above
    issueCategories.map((issue, idx) => {
        document.querySelector('#content-div').innerHTML +=
        `
            <h1 class='dropdown-category closed'>${issue.category} <span class='icon'>∆</span></h1>
            <div class='dropdown-content'>
                <p>${issue.content}</p>
            </div>
        `
    })
    document.addEventListener('click', (e) => {
        if(e.target.parentNode.classList.contains('closed')){
            e.target.parentNode.classList.remove('closed')
            e.target.parentNode.classList.add('open')
        } else {
            e.target.parentNode.classList.add('closed')
            e.target.parentNode.classList.remove('open')
        }
    })
}