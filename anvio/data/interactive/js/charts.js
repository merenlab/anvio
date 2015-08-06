/**
 * Javascript library to visualize anvi'o charts
 *
 *  Author: A. Murat Eren <a.murat.eren@gmail.com>
 *  Credits: Özcan Esen, Gökmen Göksel
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

var VIEWER_WIDTH = window.innerWidth || document.documentElement.clientWidth || document.getElementsByTagName('body')[0].clientWidth;

var layers;
var coverage;
var variability;
var geneParser;
var contextSvg;

GeneParser = (function() {
  function GeneParser(data) {
    this.data = data;
  }

  GeneParser.prototype.filterData = function(start, stop) {
    var gene, _data, _i, _len, _ref;
    _data = [];
    _ref = this.data;
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      gene = _ref[_i];
      if (gene.start_in_split > start && gene.start_in_split < stop || gene.stop_in_split > start && gene.stop_in_split < stop || gene.start_in_split < start && gene.stop_in_split > stop) {
        _data.push(gene);
      }
    }
    return _data;
  };

  return GeneParser;

})();

function getUrlVars() {
    var map = {};
    var parts = window.location.href.replace(/[?&]+([^=&]+)=([^&]*)/gi, function(m,key,value) {
        map[key] = value;
    });
    return map;
}


function loadAll() {
    contig_id = getUrlVars()["contig"];
    document.title = contig_id + " detailed";

    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/charts/' + contig_id,
        success: function(data) {
            contig_data = JSON.parse(data);

            layers = contig_data.layers;
            coverage = contig_data.coverage;
            variability = contig_data.variability;
            competing_nucleotides = contig_data.competing_nucleotides;
            previous_contig_name = contig_data.previous_contig_name;
            next_contig_name = contig_data.next_contig_name;
            index = contig_data.index;
            total = contig_data.total;
            genes = contig_data.genes;

            if(layers.length == 0){
                console.log('Warning: no layers returned')
            }

            next_str = " | next &gt;&gt;&gt;";
            prev_str = "&lt;&lt;&lt; prev | ";
            position = index + " of " + total;

            if(next_contig_name)
                next_str = '<small><a href="charts.html?contig=' + next_contig_name + '"> | next &gt;&gt;&gt;</a>';

            if(previous_contig_name)
                prev_str = '<small><a href="charts.html?contig=' + previous_contig_name + '">&lt;&lt;&lt; prev | </a>';

            document.getElementById("header").innerHTML = "<strong>" + contig_id + "</strong> detailed <br /><small>" + prev_str + position + next_str + "</small>";

            $.ajax({
              type: 'GET',
              cache: false,
              url: '/data/charts/get_state?timestamp=' + new Date().getTime(),
              success: function(state) {
                createCharts(state);
              }
            });
        }
    });

}


function createCharts(state){
    /* Adapted from Tyler Craft's Multiple area charts with D3.js article:
    http://tympanus.net/codrops/2012/08/29/multiple-area-charts-with-d3-js/  */

    var layers_ordered = state['layer-order'];

    layers_ordered = layers_ordered.filter(function (value) { if (layers.indexOf(value)>-1) return true; return false; });

    var visible_layers = 0;
    for (i in layers_ordered)
    {
      var layer_id = layers_ordered[i];

      if (parseFloat(state['layers'][layer_id]['height']) > 0)
        visible_layers++;
    }

    geneParser = new GeneParser(genes);

    var margin = {top: 20, right: 50, bottom: 150, left: 50};
    var width = VIEWER_WIDTH * .80;
    var chartHeight = 200;
    var height = (chartHeight * visible_layers + 400);
    var contextHeight = 50;
    var contextWidth = width;

    var svg = d3.select("#chart-container").append("svg")
            .attr("width", width + margin.left + margin.right)
            .attr("height", (height + margin.top + margin.bottom));

    $('#chart-container').css("width", (width + 150) + "px");
    
    var charts = [];
    var maxVariability = 0;
    
    var layersCount = layers.length;
    
    coverage.forEach(function(d) {
        for (var prop in d) {
            if (d.hasOwnProperty(prop)) {
                d[prop] = parseFloat(d[prop]);
            }
        }
    });

    variability.forEach(function(d) {
        for (var prop in d) {
            if (d.hasOwnProperty(prop)) {
                d[prop] = parseFloat(d[prop]);
                
                if (d[prop] > maxVariability) {
                    maxVariability = d[prop];
                }
            }
        }
    });

    var j=0;
    for(var i = 0; i < layersCount; i++){
        var layer_index = layers.indexOf(layers_ordered[i]);

        if (parseFloat(state['layers'][layers_ordered[i]]['height']) == 0)
          continue;

        charts.push(new Chart({
                        name: layers[layer_index],
                        coverage: coverage[layer_index],
                        variability: variability[layer_index],
                        competing_nucleotides: competing_nucleotides[layer_index],
                        id: j++,
                        width: width,
                        height: chartHeight,
                        maxVariability: maxVariability,
                        svg: svg,
                        margin: margin,
                        showBottomAxis: (j == visible_layers - 1),
                        color: state['layers'][layers[layer_index]]['color']
                }));
        
    }


    contextSvg = d3.select("#context-container").append("svg")
            .attr("width", width + margin.left + margin.right)
            .attr("height", 150);

    var defs = contextSvg.append('svg:defs');

    contextSvg.append("rect")
       .attr("width", width)
       .attr("height", "60px")
       .attr("fill", "black")
       .attr("fill-opacity", "0.2")
       .attr('transform', 'translate(50, 10)');

    // Define arrow markers
    ['green', 'gray'].forEach(function(color){
      defs.append('svg:marker')
          .attr('id', 'arrow_' + color )
          .attr('markerHeight', 2)
          .attr('markerWidth', 2)
          .attr('orient', 'auto')
          .attr('refX', 0)
          .attr('refY', 0)
          .attr('viewBox', '-5 -5 10 10')
          .append('svg:path')
            .attr('d', 'M 0,0 m -5,-5 L 5,0 L -5,5 Z')
            .attr('fill', color);
    });

    $('#context-container').css("width", (width + 150) + "px");

    /* Context down below */
    var contextXScale = d3.scale.linear().range([0, contextWidth]).domain(charts[0].xScale.domain());

    var contextAxis = d3.svg.axis()
                .scale(contextXScale)
                .tickSize(contextHeight);

    var contextArea = d3.svg.area()
                .interpolate("monotone")
                .x(function(d) { return contextXScale(d); })
                .y0(contextHeight)
                .y1(0);

    var brush = d3.svg.brush()
                .x(contextXScale)
                .on("brushend", onBrush);

    var context = contextSvg.append("g")
                .attr("class","context")
                .attr("transform", "translate(" + (margin.left) + ", 80)");

    context.append("g")
                .attr("class", "x axis top")
                .attr("transform", "translate(0,0)")
                .call(contextAxis)

    context.append("g")
                .attr("class", "x brush")
                .call(brush)
                .selectAll("rect")
                .attr("y", 0)
                .attr("height", contextHeight);

    function onBrush(){
        /* this will return a date range to pass into the chart object */
        var b = brush.empty() ? contextXScale.domain() : brush.extent();
        b = [Math.floor(b[0]), Math.floor(b[1])];
        for(var i = 0; i < layersCount; i++){
            charts[i].showOnly(b);
        }
        drawArrows(b[0], b[1]);
    }

    drawArrows(0, charts[0].xScale.domain()[1]);
}

function removeGeneChart() {
  var node = document.getElementById("gene-arrow-chart");
  if (node && node.parentNode) {
    node.parentNode.removeChild(node);
  }
}

function drawArrows(_start, _stop) {

    width = VIEWER_WIDTH * 0.80;
    genes = geneParser.filterData(_start, _stop);

    console.log("Start/Stop:", _start, _stop);
    console.log("Filtered genes:", genes);

    removeGeneChart();

    paths = contextSvg.append('svg:g')
      .attr('id', 'gene-arrow-chart')
      .attr('transform', 'translate(50, -10)');

    paths.selectAll('path');

    // Find the max_split
    max_split = _stop - _start;

    // Find the ratio based on screen width and max_split
    ratio = width / max_split;

    // Draw arrows
    genes.forEach(function(gene) {

      diff = 0;

      if (gene.start_in_split < _start) {
        start = 0;
        diff  = _start - gene.start_in_split;
      }

      start = gene.start_in_split < _start ? 0 : (gene.start_in_split - _start) * ratio;
      // start = ratio * (start - gene.start_in_split);
      __stop= (gene.stop_in_split > _stop ? _stop : gene.stop_in_split);
      stop  = (__stop - gene.start_in_split - diff) * ratio;

      var y = 10 + (gene.level * 20);

      color = (gene.function !== null ? 'green' : 'gray');

      // M10 15 l20 0
      path = paths.append('svg:path')
           .attr('d', 'M' + start +' '+ y +' l'+ stop +' 0')
           .attr('stroke', color)
           .attr('stroke-width', 5)
           .attr('marker-end', function() {

             if ((gene.percentage_in_split == 100) &&
                 (gene.direction == 'r' && gene.start_in_split > _start) ||
                 (gene.direction == 'f' && gene.stop_in_split  < _stop)) {
                   return 'url(#arrow_' + color + ')';
                 }

              return '';
           })
           .attr('transform', function() {
               return gene.direction == 'r' ? "translate(" + (2*start+stop) + ", 0), scale(-1, 1)" : "";
             })
           .append('svg:title')
             .text(gene.function + '');
    });

}


var base_colors = ['#CCB48F', '#727EA3', '#65567A', '#CCC68F', '#648F7D', '#CC9B8F', '#A37297', '#708059'];

function get_comp_nt_color(nts){
    if(nts == "CT" || nts == "TC")
        return "red";
    if(nts == "GA" || nts == "AG")
        return "green";
    if(nts == "AT" || nts == "TA")
        return "blue";
    if(nts == "CA" || nts == "AC")
        return "purple";
    if(nts == "GT" || nts == "TG")
        return "orange";
    else
        return "black";
}

function Chart(options){
    this.coverage = options.coverage;
    this.variability = options.variability;
    this.competing_nucleotides = options.competing_nucleotides;
    this.width = options.width;
    this.height = options.height;
    this.maxVariability = options.maxVariability;
    this.svg = options.svg;
    this.id = options.id;
    this.name = options.name;
    this.margin = options.margin;
    this.showBottomAxis = options.showBottomAxis;
    this.color = options.color;
    
    var localName = this.name;
    var num_data_points = this.variability.length;
    
    this.xScale = d3.scale.linear()
                            .range([0, this.width])
                            .domain([0, this.coverage.length]);
   
    this.maxCoverage = Math.max.apply(null, this.coverage);
    if(this.maxCoverage < 20)
        this.maxCoverage = 20;
    this.yScale = d3.scale.linear()
                            .range([this.height,0])
                            .domain([0,this.maxCoverage]);

    this.yScaleLine = d3.scale.linear()
                            .range([this.height, 0])
                            .domain([0, this.maxVariability]);
    
    var xS = this.xScale;
    var yS = this.yScale;
    var ySL = this.yScaleLine;
    
    this.area = d3.svg.area()
                            .x(function(d, i) { return xS(i); })
                            .y0(this.height)
                            .y1(function(d) { return yS(d); });

    this.line = d3.svg.line()
                            .x(function(d, i) { return xS(i); })
                            .y(function(d, i) { if(i == 0) return ySL(0); if(i == num_data_points - 1) return ySL(0); return ySL(d); })
                            .interpolate('step-before');

    /*
        Assign it a class so we can assign a fill color
        And position it on the page
    */
    this.chartContainer = this.svg.append("g")
                        .attr('class',this.name.toLowerCase())
                        .attr("transform", "translate(" + this.margin.left + "," + (this.margin.top + (this.height * this.id) + (10 * this.id)) + ")");

    this.lineContainer = this.svg.append("g")
                        .attr('class',this.name.toLowerCase())
                        .attr("transform", "translate(" + this.margin.left + "," + (this.margin.top + (this.height * this.id) + (10 * this.id)) + ")");

    this.textContainer = this.svg.append("g")
                        .attr('class',this.name.toLowerCase())
                        .attr("transform", "translate(" + this.margin.left + "," + (this.margin.top + (this.height * this.id) + (10 * this.id)) + ")");

    /* Add both into the page */
    this.chartContainer.append("path")
                              .data([this.coverage])
                              .attr("class", "chart")
                              .style("fill", this.color)
                              .style("fill-opacity", "0.5")
                              .attr("d", this.area);
                                    
    this.lineContainer.append("path")
                              .data([this.variability])
                              .attr("class", "line")
                              .style("stroke", '#000000')
                              .style("stroke-width", "1")
                              .attr("d", this.line);

    this.textContainer.selectAll("text")
                            .data(d3.entries(this.competing_nucleotides))
                            .enter()
                            .append("text")
                            .attr("x", function (d) { return xS(d.key); })
                            .attr("y", function (d) { return 0; })
                            .attr("writing-mode", "tb")
                            .attr("font-size", "7px")
                            .attr("glyph-orientation-vertical", "0")
                            .attr("fill", function (d){return get_comp_nt_color(d.value);})
                            .text(function (d) {
                                return d.value;
                            });


    
    this.xAxisTop = d3.svg.axis().scale(this.xScale).orient("top");

    if(this.id == 0){
        this.chartContainer.append("g")
                    .attr("class", "x axis top")
                    .attr("transform", "translate(0,0)")
                    .call(this.xAxisTop);
    }
    
        
    this.yAxis = d3.svg.axis().scale(this.yScale).orient("left").ticks(5);
    this.yAxisLine = d3.svg.axis().scale(this.yScaleLine).orient("right").ticks(5);
        
    this.chartContainer.append("g")
                   .attr("class", "y axis")
                   .attr("transform", "translate(-15,0)")
                   .call(this.yAxis);

    this.lineContainer.append("g")
                   .attr("class", "y axis")
                   .attr("transform", "translate(" + (this.width + 15) + ",0)")
                   .call(this.yAxisLine);

    this.chartContainer.append("text")
                   .attr("class","country-title")
                   .attr("transform", "translate(0,20)")
                   .text(this.name);
    
}
    
Chart.prototype.showOnly = function(b){
        this.xScale.domain(b); var xS = this.xScale;
        this.chartContainer.selectAll("path").data([this.coverage]).attr("d", this.area);
        this.lineContainer.select("path").data([this.variability]).attr("d", this.line);
        this.textContainer.selectAll("text").data(d3.entries(this.competing_nucleotides)).attr("x", function (d) { return xS(d.key); });
        this.chartContainer.select(".x.axis.top").call(this.xAxisTop);
}

