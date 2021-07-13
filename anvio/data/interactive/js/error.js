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
let ERROR_REASONS = []

const issueCategories = [
    {
        'category' : 'Dependencies failed to load', 
        'content' : `Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. 
        Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. 
        Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.`
    },
    {
        'category' : 'There are bats in my terminal', 
        'content' : 'Did you hire a cut-rate terminal sweep?'
    }
]

function alertDependencyError(dependencyError, isFinalDependency){

    if(dependencyError){
        ERROR_COUNT += 1
    }
    
    if(isFinalDependency && ERROR_COUNT){ // hacky way of 'iterating' all dependency calls before error messaging  
        displayAlert('dependencies')
    }
}

function displayAlert(reason){
    alert(`Anvi'o has encountered an error, possibly related to ${reason}. Anvi'o would like to offer some guidance in a new browser tab. Please make sure popups are enabled :)`)
    window.open('error-landing.html', '_blank')
    ERROR_REASONS.push(reason)
}

function errorLandingContext(){
    issueCategories.map(issue => {
        document.querySelector('#content-div').innerHTML += 
        `
            <h1 class='dropdown-category closed'>${issue.category} <span class='icon'>∆</span></h1>
            <div class='dropdown-content'>
                <p>${issue.content}</p>
            </div>
        `
    })
    
    document.querySelector('.icon').addEventListener('click', (e) => {
        if(e.target.parentNode.classList.contains('closed')){
            e.target.parentNode.classList.remove('closed')
            e.target.parentNode.classList.add('open')
        } else {
            e.target.parentNode.classList.add('closed')
            e.target.parentNode.classList.remove('open')
        }
    })
}


