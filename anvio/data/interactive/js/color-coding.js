/**
 * Javascript library to visualize anvi'o gene clusterss
 *
 *  Authors: Mahmoud Yousef
 *
 *  Copyright 2018-2021, The anvi'o project (http://anvio.org)
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
 **/

//Returns the length of the largest object in an array
function maxLength(arr){
    var length = 0;
    for (var i = 0; i < arr.length; i++)
    {
        if (arr[i].length > length){
            length = arr[i].length;
        }
    }
    return length;
}

//Use an array of dictionaries. The key is the amino acid, and the data is the color, or "null" if none
//Input: a 2D array of the amino acids
//Returns: a 2D array of dictionaries
function determineColor(sequences_array){
    var cols = sequences_array.length;
    var results = [];
    for (var n = 0; n < cols; n++){
        results.push([]);
    }

    var column = [];
    for(var c = 0; c < cols; c++){
        column = sequences_array[c];
        _positions = colorAlgorithm(column)
        results[c] = _positions
    }
    return results;
}


//compares all of the aa's and assigns colors or "null" to the respective amino acid
//input: the "column" array from determineColor
//output: the array with all of the color assignments\
function colorAlgorithm(positions){
    var _positions = []
    for (aa in positions){
        if (checked(positions[aa]) && aboveThreshold(positions, positions[aa])) {
            _positions[aa] = color(positions, positions[aa]);
        } else{
            var dict = {}
            dict[positions[aa]] = "BLACK";
            _positions[aa] = dict;
        }
    }
    return _positions;
}

function checked(letter){
    if (letter == undefined || letter == null){
        return false;
    }
    box =  document.getElementById(letter)
    if (box == null){
        return false;
    }
    return box.checked;
}


//does the actual comparisons
//This checks for amino acid conservation by common characteristics
function aboveThreshold(positions, aa) {
    var number = 0;
    for (acid in positions) {
        if (acid === '') {
            continue;
        }
        number++;
    }
    var count = 0.0;
    var count2 = 0.0;
    var count3 = 0.0;

    var letter = aa
    switch (letter) {
        case "A":
        case "I":
        case "L":
        case "M":
        case "F":
        case "W":
        case "V":
            for (amino in positions){
                var acid = positions[amino]
                if (acid == "W" || acid == "L" || acid == "V" || acid == "I" || acid == "M" || acid == "A" || acid == "F" || acid == "C" || acid == "H" || acid == "P") {
                    count = count + 1;
                }
                }
              if ( (100 * count) / number >= 60) {
                return true;
              }
              break;
        case "P":
        case "G":
              return true;
              break; //P and G have a 0% threshold
        case "R":
        case "K":
          for (amino in positions){
              var acid = positions[amino]
                if (acid == "K" || acid == "R") {
                    count = count + 1;
                    count2 = count2 + 1;
                }
                else if (acid == "Q") {
                    count2 = count2 + 1;
                }
              }
              if ((100 * count) / number >= 60 || (100 * count2) / number >= 80) {
                return true;
              }
              break;
        case "N":
          for (amino in positions){
              var acid = positions[amino]
                if (acid == "N"){
                    count = count + 1;
                    count2 = count2 + 1;
                }
                if (acid == "Y"){
                    count2 = count2 + 1;
                }
              }
              if ((count * 100) / number >= 50 || (count2 * 100) / number >= 85) {
                    return true;
              }
              break;
        case "C":
          for (amino in positions){
              var acid = positions[amino]
                if (acid == "W" || acid == "L" || acid == "V" || acid == "I" || acid == "M" || acid == "A" || acid == "F" || acid == "C" || acid == "H" || acid == "P") {
                    count = count + 1;
                }
                else if (acid == "C"){
                    count2 = count2 + 1;
                }
              }
              if ((count * 100) / number >= 60 || count2 == number) {
                return true;
              }
              break;
        case "Q":
          for (amino in positions){
              var acid = positions[amino]
                if (acid == "K" || acid == "R"){
                    count++;
                    count3++;
                }else if (acid == "Q" || acid == "E"){
                    count2++;
                    count3++;
                }
              }
              if ((count * 100) / number >= 60 || (count2 * 100) / number >= 50 || (count3 * 100) / number >= 85){
                return true;
              }
              break;
        case "E":
          for (amino in positions){
              var acid = positions[amino]
                if (acid == "K" || acid == "R"){
                    count++;
                }else if (acid == "Q" || acid == "E"){
                    count2++;
                    count3++;
                }else if (acid == "D"){
                    count3++;
                }
              }
              if ((count * 100) / number >= 60 || (count2 * 100) / number >= 50 || (count3 * 100) / number >= 85){
                return true;
              }
              break;
        case "D":
          for (amino in positions){
              var acid = positions[amino]
                if (acid == "K" || acid == "R"){
                    count++;
                    count2++;
                }else if (acid == "Q"){
                    count2++;
                }else if (acid == "E" || acid == "D"){
                    count3++;
                }
              }
              if ((count * 100) / number >= 60 || (count2 * 100) / number >= 85 || (count3 * 100) / number >= 50){
                return true;
              }
              break;
        case "H":
        case "Y":
          for (amino in positions){
              var acid = positions[amino]
                if (acid == "W" || acid == "L" || acid == "V" || acid == "I" || acid == "M" || acid == "A" || acid == "F" || acid == "C" || acid == "H" || acid == "P") {
                    count++;
                    count2++;
                }else if (acid == "Y" || acid == "Q"){
                    count2++;
                }
              }
              if ((count * 100) / number >= 60 || (count2 * 100) / number >= 85) {
                return true;
              }
              break;
        case "S":
        case "T":
          for (amino in positions){
              var acid = positions[amino]
                if (acid == "W" || acid == "L" || acid == "V" || acid == "I" || acid == "M" || acid == "A" || acid == "F" || acid == "C" || acid == "H" || acid == "P") {
                    count++;
                }else if (acid == "S" || acid == "T"){
                    count2++;
                }
              }
              if ((count * 100) / number >= 60 || (count2 * 100) / number >= 50) {
                return true;
              }
              break;
        default: break;
    }
    return false;

}

function color(positions, aa){
    var x = '';
    switch(aa){
        case "A":
        case "I":
        case "L":
        case "M":
        case "F":
        case "W":
        case "V":
              x = "BLUE";
              break;
        case "R":
        case "K":
              x = "RED";
              break;
        case "N":
        case "Q":
        case "S":
        case "T":
              x = "GREEN";
              break;
        case "E":
        case "D":
              x = "MAGENTA";
              break;
        case "G":
            x = "ORANGE";
              break;
        case "H":
        case "Y":
              x = "DARKTURQUOISE";
              break;
        case "P":
            x = "YELLOW";
              break;
        case "C":
                check: {
                for (acid in positions){
                    if (positions[acid] != "C"){
                        x = "BLUE";
                        break check;
                      }
                  }
                x = "HOTPINK";
              }
              break;
        default: x = null;
    }
    var dict = {}
    dict[aa] = x;
    return dict;
}


function initializeCheckBoxes(){
    var container = document.getElementById('display-conservation-controls');

    var labels = [A = "A", C = "C", D = "D", E = "E", F = "F", G = "G", H = "H",
         I = "I", K = "K", L = "L", M = "M", N = "N", P = "P", Q = "Q", R = "R",
          S = "S", T = "T", V = "V", W = "W", Y = "Y"]

    for (i = 0; i < labels.length; i++){
        var word = String(labels[i]);
        var box = document.createElement('input');
        box.type = "checkbox";
        box.name = word;
        box.value = word;
        box.id = word;
        box.style = "margin-left:6px;";
        box.onclick = ( function() {
                        return createDisplay();
                         } )

        var label = document.createElement('label')
        label.htmlFor = word;
        label.style = "margin-left:1px;";
        label.appendChild(document.createTextNode(word));

        container.appendChild(box);
        container.appendChild(label);
        box.checked = true
    }

    container.appendChild(document.createElement('br'));

    var all = document.createElement("button");
    all.innerHTML = "Check All";
    all.style = "margin-right:20px;"
    all.onclick  = (function() {
        for (i = 0; i < labels.length; i++){
            var letter = String(labels[i]);
            document.getElementById(letter).checked = true
        }
        createDisplay();
    } )
    container.appendChild(all);

    var none = document.createElement("button")
    none.innerHTML = "Uncheck All"
    none.onclick  = (function() {
        for (i = 0; i < labels.length; i++){
            var letter = String(labels[i]);
            document.getElementById(letter).checked = false;
        }
        createDisplay();
    } )
    container.appendChild(none);

}
