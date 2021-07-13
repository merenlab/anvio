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
        document.querySelector('#content-div').innerHTML += `
        <h1>${issue.category}</h1>
        <p>${issue.content}</p>
        `
    })

    // document.querySelector('#content-div').innerHTML += '<p> this can be programmatically changed based on current issue</p>'
}

const issueCategories = [
    {
        'category' : 'Dependencies failed to load', 
        'content' : 'Did you make sure to pay your dependency bill?'
    },
    {
        'category' : 'There are bats in my terminal', 
        'content' : 'Did you hire a cut-rate terminal sweep?'
    }
]