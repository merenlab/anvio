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
let FAILED_DEPENDENCIES = []

function alertDependencyError(dependencyError, isFinalDependency){

    if(dependencyError){
        ERROR_COUNT += 1
        FAILED_DEPENDENCIES.push(dependencyError)
    }
    
    if(isFinalDependency && ERROR_COUNT){ // hacky way of 'iterating' all dependency calls before error messaging  
        ERROR_COUNT === 1 ? 
        alert(`${ERROR_COUNT} dependency failed to load :(. The culprit is: ${FAILED_DEPENDENCIES}`)
        :
        alert(`${ERROR_COUNT} dependencies failed to load :(. These are the culprits: ${FAILED_DEPENDENCIES}`)
    }
}