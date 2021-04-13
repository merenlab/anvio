/**
 * Javascript library for anvi'o genome view
 *
 *  Authors: Isaac Fink <iafink@uchicago.edu>
 *           Matthew Klein <mtt.l.kln@gmail.com>
 *           A. Murat Eren <a.murat.eren@gmail.com>
 *
 * Copyright 2021, The anvi'o project (http://anvio.org)
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

 var VIEWER_WIDTH = window.innerWidth || document.documentElement.clientWidth || document.getElementsByTagName('body')[0].clientWidth;

 var canvas;
 var genomeMax = 0;
 var genomes;

 // Settings vars

 var spacing = 30; // genome spacing
 var showLabels = true; // show genome labels?
 var showScale = true; // show nt scale?
 var scale = 100; // nt scale intervals

 var alignToGC = null;

 var color_db;
 var cog_annotated = true, kegg_annotated = false;
 // doesn't make sense to iterate through every gene with dozens of genomes...
 // will need to find an efficient way to find these automatically

var genomeData;

$(document).ready(function() {
    initData();

    function loadAll() {
      canvas = new fabric.Canvas('myCanvas');
      $('#tooltip-body').hide() // set initual tooltip hide value
      $('#toggle_label_box').attr("checked", showLabels);

      // can either set it on the canvas to check for all arrows, or when arrow is created.
      canvas.on('mouse:down', function(options) {
        if(options.target && options.target.id == 'arrow') {
          options.target.set('fill', options.target.fill=="red" ? "blue" : 'red');
          var testid = options.target.gene.gene_callers_id;
          console.log(options.target.gene);
        }
      });

      canvas.on('mouse:over', (event) => {
        if(event.target && event.target.id === 'arrow'){
          showToolTip(event)
        }
      })

      canvas.on('mouse:out', (event) => {
        $('#tooltip-body').html('').hide()
      })

      function showToolTip(event){
        $('#tooltip-body').show().append(`
          <p></p>
          <style type="text/css">
            .tftable {font-size:12px;color:#333333;width:100%;border-width: 1px;border-color: #729ea5;border-collapse: collapse;}
            .tftable th {font-size:12px;background-color:#acc8cc;border-width: 1px;padding: 8px;border-style: solid;border-color: #729ea5;text-align:left;}
            .tftable tr {background-color:#d4e3e5;}
            .tftable td {font-size:12px;border-width: 1px;padding: 8px;border-style: solid;border-color: #729ea5;}
            .tftable tr:hover {background-color:#ffffff;}
          </style>

          <table class="tftable" border="1">
            <tr><th>Data</th><th>Value</th></tr>
            <tr><td>Split</td><td>${event.target.gene.split}</td></tr>
            <tr><td>Start in Contig</td><td>${event.target.gene.start_in_contig}</td></tr>
            <tr><td>Length</td><td>${event.target.gene.length}</td></tr>
            <tr><td>Gene Callers ID</td><td>${event.target.gene.gene_callers_id}</td></tr>
            <tr><td>Gene Cluster</td><td>${idToGC[event.target.gene.gene_callers_id]}</td></tr>
          </table>
          <button>some action</button>
          <button>some other action</button>
          `).css({'position' : 'absolute', 'left' : event.e.clientX, 'top' : event.e.clientY })
      }

      genomes = [contig437, contig1001, contig798];

      $('#gene_color_order').append($('<option>', {
        value: 'Source',
        text: 'Source'
      }));
      if(cog_annotated) {
        $('#gene_color_order').append($('<option>', {
          value: 'COG',
          text: 'COG'
        }));
      }
      if(kegg_annotated) {
        $('#gene_color_order').append($('<option>', {
          value: 'KEGG',
          text: 'KEGG'
        }));
      }

      draw();

      // zooming and panning
      // http://fabricjs.com/fabric-intro-part-5#pan_zoom
      canvas.on('mouse:down', function(opt) {
        var evt = opt.e;
        if (evt.shiftKey === true) {
          this.isDragging = true;
          this.selection = false;
          this.lastPosX = evt.clientX;
          this.lastPosY = evt.clientY;
        }
      });
      canvas.on('mouse:move', function(opt) {
        if (this.isDragging) {
          var e = opt.e;
          var vpt = this.viewportTransform;
          vpt[4] += e.clientX - this.lastPosX;
          vpt[5] += e.clientY - this.lastPosY;
          this.requestRenderAll();
          this.lastPosX = e.clientX;
          this.lastPosY = e.clientY;
        }
      });
      canvas.on('mouse:up', function(opt) {
        // on mouse up we want to recalculate new interaction
        // for all objects, so we call setViewportTransform
        this.setViewportTransform(this.viewportTransform);
        this.isDragging = false;
        this.selection = true;
      });
      canvas.on('mouse:wheel', function(opt) {
        var delta = opt.e.deltaY;
        var zoom = canvas.getZoom();
        zoom *= 0.999 ** delta;
        if (zoom > 20) zoom = 20;
        if (zoom < 0.01) zoom = 0.01;
        canvas.zoomToPoint({ x: opt.e.offsetX, y: opt.e.offsetY }, zoom);
        opt.e.preventDefault();
        opt.e.stopPropagation();
      });

      $('#geneClusterInput').on('keydown', function(e) {
        if(e.keyCode == 13) { // 13 = enter key
          alignToCluster($(this).val());
          $(this).blur();
        }
      });
      document.body.addEventListener("keydown", function(ev) {
          if(ev.which == 83 && !$('#geneClusterInput').is(':focus')) { // S = 83
            toggleSettingsPanel();
          }
      });
      $('#gene_color_order').on('change', function() {
          color_db = $(this).val();
          canvas.clear();
          draw();
          $(this).blur();
      });
      $('#toggle_label_box').on('change', function() {
        showLabels = !showLabels;
        alignToGC = null;
        canvas.clear();
        draw();
      });
    }

    function draw() {
      for(genes of genomes) {
        let g = genes[genes.length-1].stop_in_split;
        if(g > genomeMax) genomeMax = g;
      }

      var i;
      for(i = 0; i < genomes.length; i++) {
        var genome = genomes[i];
        var label = genome[0].split.substring(9,genome[0].split.indexOf('_', 9));
        addGenome(label, genome, i+1);
      }

      if(showScale) {
        for(var w = 0; w < genomeMax; w+=scale) {
          canvas.add(new fabric.Line([0,0,0,20], {left: w+(showLabels?120:0),
                top: (i+1)*(spacing)-24,
                stroke: 'black',
                strokeWidth: 1,
                fontSize: 10,
                fontFamily: 'sans-serif',
                selectable: false}));

          canvas.add(new fabric.Text(w/1000 + " kB", {left: w+5+(showLabels?120:0),
                top: (i+1)*(spacing)-24,
                stroke: 'black',
                strokeWidth: .25,
                fontSize: 15,
                fontFamily: 'sans-serif',
                selectable: false}));
        }

        canvas.add(new fabric.Line([0,0,100,0], {left: (showLabels?120:0),
              top: (i+1)*(1.25*spacing)-4,
              stroke: 'black',
              strokeWidth: 2,
              selectable: false}));
        canvas.add(new fabric.Text("100 nts", {left: 15+(showLabels?120:0),
              top: (i+1)*(1.25*spacing)-4,
              stroke: 'black',
              strokeWidth: 1,
              fontSize: 20,
              fontFamily: 'sans-serif',
              selectable: false}));
      }
    }

    function alignToCluster(gc) {
      if(!gc || gc in mock_gene_clusters) {
        alignToGC = gc;
        showLabels = !gc; // only show labels if changing to default view
        $('#toggle_label_box').attr("checked", showLabels);
      } else {
        console.log('Warning: ' + gc + ' is not a gene cluster in data structure');
      }
      canvas.clear();
      draw();
    }

    function addGenome(label, gene_list, y) {
      var offsetX = 0;
      if(alignToGC) {
        var targetGeneID = mock_gene_clusters[alignToGC][label];
        var targetGene = gene_list.find(gene => gene.gene_callers_id == targetGeneID);
        var genePos = targetGene.start_in_split + (targetGene.stop_in_split - targetGene.start_in_split) / 2;
        //var windowCenter = fabric.util.transformPoint({x:canvas.getWidth()/2,y:0}, canvas.viewportTransform)['x'];
        var windowCenter = canvas.getWidth()/2 - canvas.viewportTransform[4]; // canvas.getWidth()/2 is clientX of center of screen
        offsetX = windowCenter - genePos;
        console.log('offsetX: ' + offsetX + ', genePos: ' + genePos + ', windowCenter: ' + windowCenter);
      }

      // label
      if(showLabels) {
        canvas.add(new fabric.Text(label, {top: spacing*y-10, selectable: false, fontSize: 15, fontFamily: 'sans-serif', fontWeight: 'bold'}));
      }

      // line
      canvas.add(new fabric.Line([0,0,genomeMax,0], {left: offsetX+(showLabels?120:0),
            top: spacing*y - 4,
            stroke: 'black',
            strokeWidth: 2,
            selectable: false}));

      // here we can either draw genes individually, or collectively as a group
      // grouping allows you to transform them all at once, but makes selecting individual arrows difficult
      // so for now they are drawn individually

      //var geneGroup = new fabric.Group();
      for(gene of gene_list) {
        //addGene(gene, y);
        //geneGroup.addWithUpdate(geneArrow(gene,y));   // IMPORTANT: only way to select is to select the group or use indices. maybe don't group them but some alternative which lets me scale them all at once?
        var geneObj = geneArrow(gene,y);
        if(showLabels) {
          geneObj.left += 120;
        }
        if(alignToGC) {
          geneObj.left += offsetX;
        }
        canvas.add(geneObj);
      }
      //canvas.add(geneGroup.set('scaleX',canvas.getWidth()/genomeMax/3));
      //geneGroup.destroy();
    }

    function geneArrow(gene, y) {
      var cag = null;
      var color = 'gray';
      if(gene.functions) {
        switch(color_db) {
          case 'COG':
            if(gene.functions["COG14_CATEGORY"]) cag = gene.functions["COG14_CATEGORY"][0][0];
            if(gene.functions["COG20_CATEGORY"]) cag = gene.functions["COG20_CATEGORY"][0][0];
            color = cag in default_COG_colors ? default_COG_colors[cag] : 'gray';
            break;
          case 'KEGG':
            if(gene.functions.hasOwnProperty("KEGG_Class") && gene.functions.KEGG_Class != null) {
              cag = getCategoryForKEGGClass(gene.functions["KEGG_Class"][1]);
            }
            color = cag in default_KEGG_colors ? default_KEGG_colors[cag] : 'gray';
            break;
          default:
            if (gene.source.startsWith('Ribosomal_RNA')) {
              cag = 'rRNA';
            } else if (gene.source == 'Transfer_RNAs') {
              cag = 'tRNA';
            } else if (gene.functions !== null) {
              cag = 'Function';
            }
            color = cag in default_source_colors ? default_source_colors[cag] : 'gray';
        }
      }
      /* Issue here: each genome might be differentially annotated... how to make sure all have COG annotations for example? */

      var length = gene.stop_in_split-gene.start_in_split;
      var arrow = new fabric.Path('M 0 0 L ' + length + ' 0 L ' + length + ' 10 L 0 10 M ' + length + ' 0 L ' + length + ' 20 L ' + (25+length) + ' 5 L ' + length + ' -10 z');
      arrow.set({
        id: 'arrow',
        gene: gene,   // better not to store entire gene object, but a pointer/id to find it in the genomes dict?
        selectable: false,
        top: -11+spacing*y,
        left: 1.5+gene.start_in_split,
        scaleX: 0.5,
        scaleY: 0.5,
        fill: color,
        zoomX: 0.2,
        zoomY: 0.2
      });
      if(gene.direction == 'r') arrow.rotate(180);

      return arrow;
    }

    loadAll()
});

function initData() {
// initialize the bulk of the data.
    $.ajax({
        type: 'POST',
        cache: false,
        url: '/data/get_genome_view_data',
        success: function(data) {
            genomeData = data;
            console.log("Saved the following data:");
            console.log(data);
        }
    });
}

function toggleSettingsPanel() {
    $('#settings-panel').toggle();

    if ($('#settings-panel').is(':visible')) {
        $('#toggle-panel-settings').addClass('toggle-panel-settings-pos');
        $('#toggle-panel-settings-inner').html('&#9658;');
    } else {
        $('#toggle-panel-settings').removeClass('toggle-panel-settings-pos');
        $('#toggle-panel-settings-inner').html('&#9664;');
    }
}

function getCategoryForKEGGClass(class_str) {
  if(class_str == null) return null;

  var category_name = getClassFromKEGGAnnotation(class_str);
  return getKeyByValue(KEGG_categories, category_name);
}

function getClassFromKEGGAnnotation(class_str) {
  return class_str.substring(17, class_str.indexOf(';', 17));
}

// https://stackoverflow.com/questions/9907419/how-to-get-a-key-in-a-javascript-object-by-its-value/36705765
function getKeyByValue(object, value) {
  return Object.keys(object).find(key => object[key] === value);
}

var mock_gene_clusters = {'GC_X': {'contig437': 14902,
                                   'contig1001': 19391,
                                   'contig798': 18019},
                          'GC_Y': {'contig437': 14937,
                                   'contig1001': 19393,
                                   'contig798': 18011}
}

var idToGC = {
  14902: 'GC_X',
  19391: 'GC_X',
  18019: 'GC_X',
  14937: 'GC_Y',
  19393: 'GC_Y',
  18011: 'GC_Y'
}
