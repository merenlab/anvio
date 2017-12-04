/**
 * Javascript library to visualize anvi'o protein clusterss
 *
 *  Author: A. Murat Eren <a.murat.eren@gmail.com>
 *  Credits: Özcan Esen
 *  Copyright 2016, The anvio Project
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


/****************************************************************************
 * The following code is from "colorcode.js"
 ****************************************************************************/
/**
 * This function will be used to test if a dictionary is empty
 */
function isEmpty(obj){
    return Object.keys(obj).length === 0;
  }
  
  //Returns the length of the largest object (string, in our case) in an array
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
  //Input: an array of strings (characters)
  //Returns: a 2D array of dictionaries
  function determineColor(sequences_array){
      console.log(sequences_array);
      var cols = maxLength(sequences_array);
      var rows = sequences_array.length;
      var results = [];
      for (var n = 0; n < rows; n++){
          results.push([]);
          for (var c = 0; c < cols; c++){
              results[n].push("");
          }
      }
      console.log(results);
      //return;


      for(var i = 0; i < cols; i++) { 
          var column = [];
          for (var l = 0; l < rows; l++) {
              var val = ""
/*               if (sequences_array[l] === ""){
                  continue;
              } */
              if(i < sequences_array[l].length) {
                  var x = sequences_array[l][i];
                  column.push(x);
                  //console.log(positions[l]);
            }
          }
          console.log(column);
          var _positions = colorAlgorithm(column);
          console.log(_positions);
          for (var l = 0; l < rows; l++) {
              if (_positions[l] != undefined){
                  //console.log(_positions[l]);
                  results[l].push(_positions[l]); //why does this drop results[l].size to 0?
                  //console.log(results[l][i]);
                  //console.log(results);
              }
          }
      }

      console.log(results);
      return results;
  }
  
  
  //compares all of the aa's and assigns colors or "null" to the respective amino acid
  //input: the "positions" array from determineColor
  //output: the array with all of the color assignments
  //I can potentially combine all three of these into a single function to be more efficient.
  function colorAlgorithm(positions){
      for (aa in positions){
          if (aboveThreshold(positions,positions[aa])) {
              positions[aa] = color(positions, positions[aa]); 	
          }
      }
      return positions;
  }
  
  
  //does the actual comparisons
  function aboveThreshold(positions,aa) { 
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
      
      var letter = Object.keys(aa)[0];
      switch (letter) {
          case "A":
          case "I":
          case "L":
          case "M":
          case "F":
          case "W":
          case "V":
              for (amino in positions){
                  var acid = Object.keys(positions[amino])[0];
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
                var acid = Object.keys(positions[amino])[0];
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
                var acid = Object.keys(positions[amino])[0];
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
                var acid = Object.keys(positions[amino])[0];
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
                var acid = Object.keys(positions[amino])[0];
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
                var acid = Object.keys(positions[amino])[0];
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
                var acid = Object.keys(positions[amino])[0];
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
                var acid = Object.keys(positions[amino])[0];
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
                var acid = Object.keys(positions[amino])[0];
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
      /*Two problems: 
          2) Need to put in the color in the definition slot, not just return it
          */
      var x = '';
      var letter = Object.keys(aa)[0];
      switch(letter){
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
                x = "CYAN";
                break;
          case "P": 
              x = "YELLOW";
                break;
          case "C":
                  check: {
                  for (amino in positions){
                      var acid = Object.keys(positions[amino])[0];
                      if (acid != "C"){
                          x = "BLUE";
                          break check;
                        }
                    }
                  x = "PINK";
                }
                break;
          default: x = null;
      }
      aa[letter] = x;
      return aa;
  }
  
  //Returns the hexadecimal value of the colors
  //Precondition: the amino acid is of type {letter:color}
  function readColor(aa) {
      var letter = Object.keys(aa)[0];
      var color = aa[letter];
      var ans;
      switch(color){
          case "BLUE" : ans = 0x0000FF;
                                      break;
          case "RED" : ans = 0xFF0000;
                                      break;
          case "GREEN" : ans = 0x008000;
                                      break;
          case "MAGENTA" : ans = 0xFF00FF;
                                      break;
          case "ORANGE" : ans = 0xFFA500
                                      break;
          case "CYAN" : ans = 0x00FFFF;
                                      break;
          case "YELLOW" : ans = 0xFFFF00;
                                      break;
          case "PINK" : ans = 0xFFC0CB;
                                      break;
          default: ans = 0x000000;
      }
      return ans;
  }

  /********************************************************************
   * END COLORCODE.JS
   ********************************************************************/


var VIEWER_WIDTH = window.innerWidth || document.documentElement.clientWidth || document.getElementsByTagName('body')[0].clientWidth;

var genomes;
var gene_caller_ids;
var gene_caller_ids_in_genomes;
var aa_sequences_in_pc;
var previous_pc_name;
var next_pc_name;
var index;
var total;

var state;
var pc_data;


function loadAll() {
    pc_name = getUrlVars()["id"];
    document.title = pc_name + " detailed";

    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/proteinclusters/' + pc_name,
        success: function(_pc_data) {
            pc_data = _pc_data;
            genomes = pc_data.genomes;
            gene_caller_ids = pc_data.gene_caller_ids;
            gene_caller_ids_in_genomes = pc_data.gene_caller_ids_in_genomes;
            aa_sequences_in_pc = pc_data.aa_sequences_in_pc;
            previous_pc_name = pc_data.previous_pc_name;
            next_pc_name = pc_data.next_pc_name;
            index = pc_data.index;
            total = pc_data.total;

            if(genomes.length == 0){
                console.log('Warning: no genomes returned')
            }

            next_str = " | next &gt;&gt;&gt;";
            prev_str = "&lt;&lt;&lt; prev | ";
            position = index + " of " + total;

            if(next_pc_name)
                next_str = '<a onclick="sessionStorage.state = JSON.stringify(state, null, 4);" href="proteinclusters.html?id=' + next_pc_name + '"> | next &gt;&gt;&gt;</a>';

            if(previous_pc_name)
                prev_str = '<a onclick="sessionStorage.state = JSON.stringify(state, null, 4);" href="proteinclusters.html?id=' + previous_pc_name + '">&lt;&lt;&lt; prev | </a>';

            document.getElementById("header").innerHTML = "<strong>" + pc_name + "</strong> with " + gene_caller_ids.length + " genes detailed <br /><small><small>" + prev_str + position + next_str + "</small></small>";

            if (typeof sessionStorage.state === 'undefined')
            {
                alert("Something went wrong, couldn't access to sessionStorage");
            }
            else
            {
                // backup the state, if user changes the page (prev, next) we are going to overwrite it.
                state = JSON.parse(sessionStorage.state);
                createDisplay();
            }
        }
    });

}


function createDisplay(){
    var sequence_wrap_val = parseInt($('#wrap_length').val());
    var sequence_font_size_val = parseInt($('#font_size').val());

    var sequence_wrap = (isNumber(sequence_wrap_val) && sequence_wrap_val > 0) ? sequence_wrap_val : 140;
    var sequence_font_size = (isNumber(sequence_font_size_val) && sequence_font_size_val > 0) ? sequence_font_size_val : 12;
    
    var svg = document.getElementById('svg');

    // clear content
    while (svg.firstChild) {
        svg.removeChild(svg.firstChild);
    }
    
    var y_cord = 0;
    var offset = 0;

    var acid_sequences = [];

    for (var index = 0; index < state['layer-order'].length; index++){
        var seq = state['layer-order'][index];
        var id = pc_data.genomes.indexOf(seq);
        if (id === -1)
            continue;
       var call_id = pc_data.gene_caller_ids_in_genomes[seq]; //call_id is an array; empty, or one element with the id
       var str = (pc_data.aa_sequences_in_pc[seq])[call_id];
       if (str === undefined) {
           str = "";
       }
       acid_sequences.push(str);
    }
    var colored_sequences = determineColor(acid_sequences); //I have to implement caller ids and make a dictionary out of this
    console.log(colored_sequences);
    //colored_sequences fails
    for (var index = 0; index < state['layer-order'].length; index++){
        var seq = state['layer-order'][index];
        var id = pc_data.genomes.indexOf(seq);
        if (id === -1)
            continue;
       var call_id = pc_data.gene_caller_ids_in_genomes[seq];
       (pc_data.aa_sequences_in_pc[seq])[call_id] = colored_sequences[index];
    }
    console.log(pc_data.aa_sequences_in_pc);
    //deal with efficiency later

    while (true)
    {
        ignoreTable = true;
        var fragment = document.createDocumentFragment();
        for (var layer_id = 0; layer_id < state['layer-order'].length; layer_id++)
        {
            var layer = state['layer-order'][layer_id];

            if (pc_data.genomes.indexOf(layer) === -1)
                continue;
            
            var rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
            rect.setAttribute('x', 0);
            rect.setAttribute('y', y_cord);
            rect.setAttribute('fill', state['layers'][layer]['color']);
            rect.setAttribute('opacity', 0.2);
            rect.setAttribute('height', (Math.max(pc_data.gene_caller_ids_in_genomes[layer].length,1) * sequence_font_size * 1.5) + 10);
            rect.setAttribute('width', 400);
            rect.setAttribute('class', 'sequenceBackground');
            rect.setAttribute('rx', 10);
            fragment.appendChild(rect);

            var text = document.createElementNS('http://www.w3.org/2000/svg', 'text');
            text.setAttribute('x', 0);
            text.setAttribute('y', parseFloat(rect.getAttribute('y')) + parseFloat(rect.getAttribute('height')) / 2);
            text.setAttribute('font-size', "24px");
            text.setAttribute('font-family', "Lato, Arial");
            text.setAttribute('font-weight', '300')
            text.setAttribute('style', 'alignment-baseline:central');
            text.setAttribute('class', 'genomeTitle')
            text.appendChild(document.createTextNode(layer));
            fragment.appendChild(text);

            sub_y_cord = y_cord + 5;
	         pc_data.gene_caller_ids_in_genomes[layer].forEach(function (caller_id) {
                sequence = pc_data.aa_sequences_in_pc[layer][caller_id]; //the entire aa sequence of a specie
                var text = document.createElementNS('http://www.w3.org/2000/svg', 'text');
                text.setAttribute('x', 0);
                text.setAttribute('y', sub_y_cord);
                text.setAttribute('font-size', sequence_font_size);
                text.setAttribute('font-family', "Lato, Arial");
                text.setAttribute('font-weight', '300');
                text.setAttribute('style', 'alignment-baseline:text-before-edge');
                text.setAttribute('class', 'callerTitle')
                text.appendChild(document.createTextNode(caller_id));
                fragment.appendChild(text);

                var text = document.createElementNS('http://www.w3.org/2000/svg', 'text');
                text.setAttribute('x', 0);
                text.setAttribute('y', sub_y_cord);
                text.setAttribute('font-size', sequence_font_size);
                text.setAttribute('font-family', "monospace");
                text.setAttribute('font-weight', '100');
                text.setAttribute('style', 'alignment-baseline:text-before-edge');
                text.setAttribute('class', 'sequence');
                //text.appendChild(document.createTextNode(sequence.substr(offset, sequence_wrap)));
                //text.appendChild(document.createTextNode());

                _sequence = sequence.splice(offset, sequence_wrap);
                for (var _letter_index=0; _letter_index < _sequence.length; _letter_index++) {
                    var tspan = document.createElementNS('http://www.w3.org/2000/svg', 'tspan');
                    var num = Object.keys(_sequence)[_letter_index];
                    var acid = _sequence[num];
                    var txt = Object.keys(acid)[0];
                    tspan.setAttribute('fill', readColor(acid).toString()); //USE THIS to change color
                    tspan.appendChild(document.createTextNode(txt));
		            tspan.setAttribute('style', 'alignment-baseline:text-before-edge');
                    text.appendChild(tspan);
                } 
		

                fragment.appendChild(text);

                sub_y_cord = sub_y_cord + sequence_font_size * 1.5;

                if (offset < sequence.length) {
                    ignoreTable = false;
                }
            });

            y_cord = y_cord + parseFloat(rect.getAttribute('height')) + 5;

        }
        if (ignoreTable) {
            break;
        } else {
            svg.appendChild(fragment);
        }
        offset += sequence_wrap;
        y_cord = y_cord + 15;
    }

    calculateLayout();
}

function calculateLayout() {
    var HORIZONTAL_PADDING = 10;
    var svg = document.getElementById('svg');

    var max_genome_title_width = 0;
    $('.genomeTitle').each(function(i, d) {
        dwidth = d.getBBox().width;
        if (dwidth > max_genome_title_width) {
            max_genome_title_width = dwidth;
        };
    });

    $('.sequenceBackground').attr('x', max_genome_title_width + HORIZONTAL_PADDING);
    $('.callerTitle').attr('x', max_genome_title_width + (HORIZONTAL_PADDING * 2));

    var max_caller_width = 0;
    $('.callerTitle').each(function(i, d) {
        dwidth = d.getBBox().width;
        if (dwidth > max_caller_width) {
            max_caller_width = dwidth;
        };
    });

    $('.sequence').attr('x', max_genome_title_width + max_caller_width + (HORIZONTAL_PADDING * 3));

    var max_sequence_width = 0;
    $('.sequence').each(function(i, d) {
        dwidth = d.getBBox().width;
        if (dwidth > max_sequence_width) {
            max_sequence_width = dwidth;
        };
    });

    $('.sequenceBackground').attr('width', max_caller_width + max_sequence_width + (HORIZONTAL_PADDING * 3));

    bbox = svg.getBBox();
    svg.setAttribute('width', bbox.width);
    svg.setAttribute('height', bbox.height);
}


function removeGeneChart() {
  var node = document.getElementById("gene-arrow-chart");
  if (node && node.parentNode) {
    node.parentNode.removeChild(node);
  }
}
