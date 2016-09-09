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

var VIEWER_WIDTH = window.innerWidth || document.documentElement.clientWidth || document.getElementsByTagName('body')[0].clientWidth;

var genomes;
var gene_caller_ids_in_genomes;
var next_pc_name;
var index;
var total;


function getUrlVars() {
    var map = {};
    var parts = window.location.href.replace(/[?&]+([^=&]+)=([^&]*)/gi, function(m,key,value) {
        map[key] = value;
    });
    return map;
}


function loadAll() {
    pc_name = getUrlVars()["id"];
    document.title = pc_name + " detailed";

    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/proteinclusters/' + pc_name,
        success: function(data) {
            pc_data = JSON.parse(data);

            genomes = pc_data.genomes;
            gene_caller_ids = pc_data.gene_caller_ids;
            gene_caller_ids_in_genomes = pc_data.gene_caller_ids_in_genomes;
            gene_aa_sequences_for_gene_caller_ids = pc_data.gene_aa_sequences_for_gene_caller_ids;
            previous_pc_name = pc_data.previous_pc_name;
            next_pc_name = pc_data.next_pc_name;
            index = pc_data.index;
            total = pc_data.total;

            console.log(pc_data);

            if(genomes.length == 0){
                console.log('Warning: no genomes returned')
            }

            next_str = " | next &gt;&gt;&gt;";
            prev_str = "&lt;&lt;&lt; prev | ";
            position = index + " of " + total;

            if(next_pc_name)
                next_str = '<a href="proteinclusters.html?pc=' + next_pc_name + '"> | next &gt;&gt;&gt;</a>';

            if(previous_pc_name)
                prev_str = '<a href="proteinclusters.html?pc=' + previous_pc_name + '">&lt;&lt;&lt; prev | </a>';

            document.getElementById("header").innerHTML = "<strong>" + pc_name + "</strong> with " + gene_caller_ids.length + " genes detailed <br /><small><small>" + prev_str + position + next_str + "</small></small>";

            $.ajax({
              type: 'GET',
              cache: false,
              url: '/data/proteinclusters/get_state?timestamp=' + new Date().getTime(),
              success: function(state) {
                createDisplay(state);
              }
            });
        }
    });

}


function createDisplay(state){
    var genomes_ordered;

    genomes_ordered = state['layer-order'].filter(function (value) { if (genomes.indexOf(value)>-1) return true; return false; });

    var visible_genomes = 0;
    for (i in genomes_ordered)
    {
      var genome_id = genomes_ordered[i];

      if (parseFloat(state['layers'][genome_id]['height']) > 0)
        visible_genomes++;
    }

    var margin = {top: 20, right: 50, bottom: 150, left: 50};
    var width = VIEWER_WIDTH * .80;
    var chartHeight = 200;
    var height = (chartHeight * visible_genomes + 400);
    var contextHeight = 50;
    var contextWidth = width;

    var svg = d3.select("#chart-container").append("svg")
            .attr("width", width + margin.left + margin.right)
            .attr("height", (height + margin.top + margin.bottom));

    $('#chart-container').css("width", (width + 150) + "px");
    
    var proteinclusters = [];
    
    var genomesCount = genomes.length;
    
    var j=0;
    for(var i = 0; i < genomesCount; i++){
        var genome_index = genomes.indexOf(genomes_ordered[i]);

        if (parseFloat(state['layers'][genomes_ordered[i]]['height']) == 0)
          continue;

        proteinclusters.push(new Chart({
                        genome_name: genomes[genome_index],
                        gene_caller_ids_in_genome: gene_caller_ids_in_genomes[genome_name],
                        id: j++,
                        width: width,
                        height: chartHeight,
                        svg: svg,
                        margin: margin,
                        showBottomAxis: (j == visible_genomes - 1),
                        color: state['layers'][genomes[genome_index]]['color']
                }));
        
    }
}


function removeGeneChart() {
  var node = document.getElementById("gene-arrow-chart");
  if (node && node.parentNode) {
    node.parentNode.removeChild(node);
  }
}


function get_gene_functions_table_html(gene){
    functions_table_html = '<h2>Gene Call</h2>';
    functions_table_html += '<table class="table table-striped" style="width: 600px; text-align: center;">';
    functions_table_html += '<thead><th>ID</th><th>Source</th><th>Length</th><th>Direction</th><th>Start</th><th>Stop</th><th>Complete</th><th>% in split</th></thead>';
    functions_table_html += '<tbody>';
    functions_table_html += '<tr><td>' + gene.gene_callers_id
                          + '</td><td>' + gene.source
                          + '</td><td>' + gene.length
                          + '</td><td>' + gene.direction
                          + '</td><td>' + gene.start_in_pc
                          + '</td><td>' + gene.stop_in_pc
                          + '</td><td>' + gene.complete_gene_call
                          + '</td><td>' + gene.percentage_in_split.toFixed(2) + '%'
                          + '</td></tr></tbody></table>';

    functions_table_html += '<button type="button" class="btn btn-default btn-sm" onClick="show_sequence(' + gene.gene_callers_id + ');">Get sequence</button> ';
    functions_table_html += '<button type="button" class="btn btn-default btn-sm" onClick="fire_up_ncbi_blast(' + gene.gene_callers_id + ', \'blastn\', \'nr\', \'gene\');">blastn @ nr</button> ';
    functions_table_html += '<button type="button" class="btn btn-default btn-sm" onClick="fire_up_ncbi_blast(' + gene.gene_callers_id + ', \'blastn\', \'refseq_genomic\', \'gene\');">blastn @ refseq_genomic</button> ';
    functions_table_html += '<button type="button" class="btn btn-default btn-sm" onClick="fire_up_ncbi_blast(' + gene.gene_callers_id + ', \'blastx\', \'nr\', \'gene\');">blastx @ nr</button> ';
    functions_table_html += '<button type="button" class="btn btn-default btn-sm" onClick="fire_up_ncbi_blast(' + gene.gene_callers_id + ', \'blastn\', \'refseq_genomic\', \'gene\');">blastx @ refseq_genomic</button> ';

    if(!gene.functions)
        return functions_table_html;

    functions_table_html += '<h2>Annotation</h2>';
    functions_table_html += '<table class="table table-striped">';
    functions_table_html += '<thead><th>Source</th><th>Hit</th><th>Score</th></thead>';
    functions_table_html += '<tbody>';

    for (function_source in gene.functions){
        functions_table_html += '<tr>';

        functions_table_html += '<td><b>' + function_source + '</b></td>';
        if (gene.functions[function_source]) {
            functions_table_html += '<td>' + gene.functions[function_source][0] + '</td>';
            functions_table_html += '<td><em>' + gene.functions[function_source][1] + '</em></td>';
        } else {
            functions_table_html += '<td>&nbsp;</td>';
            functions_table_html += '<td>&nbsp;</td>';
        }

        functions_table_html += '</tr>';
    }

    functions_table_html += '</tbody></table>';

    return functions_table_html;
}


function show_sequence(gene_id) {
    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/gene/' + gene_id + '?timestamp=' + new Date().getTime(),
        success: function(data) {
            $('body').append('<div class="modal modal-sequence"> \
                <div class="modal-dialog"> \
                    <div class="modal-content"> \
                        <div class="modal-header"> \
                            <button class="close" data-dismiss="modal" type="button"><span>&times;</span></button> \
                            <h4 class="modal-title">Split Sequence</h4> \
                        </div> \
                        <div class="modal-body"> \
                            <div class="col-md-12"> \
                                <textarea class="form-control" rows="16" onclick="$(this).select();" readonly>&gt;' + data['header'] + '\n' + data['sequence'] + '</textarea> \
                            </div> \
                        </div> \
                        <div class="modal-footer"> \
                            <button class="btn btn-default" data-dismiss="modal" type="button">Close</button> \
                        </div> \
                    </div> \
                </div> \
            </div>');
            $('[data-toggle="popover"]').popover('hide');
            $('.modal-sequence').modal('show');
        }
    });
}

var base_colors = ['#CCB48F', '#727EA3', '#65567A', '#CCC68F', '#648F7D', '#CC9B8F', '#A37297', '#708059'];

function Chart(options){
    this.genome_name = options.genome_name;
    this.gene_caller_ids_in_genome = options.gene_caller_ids_in_genome
    this.id = options.id;
    this.width = options.width;
    this.height = options.height;
    this.svg = options.svg;
    this.margin = options.margin;
    this.showBottomAxis = options.showBottomAxis;
    this.color = options.color;

    console.log(this_genome_name)

    var num_data_points = this.variability_a.length;

    this.yScale = d3.scale.linear()
                            .range([this.height,0])
                            .domain([0, 20]);

    this.yScaleLine = d3.scale.linear()
                            .range([this.height, 0])
                            .domain([0, 20]);
    
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
                        .attr('class',this.genome_name.toLowerCase())
                        .attr("transform", "translate(" + this.margin.left + "," + (this.margin.top + (this.height * this.id) + (10 * this.id)) + ")");

    this.lineContainer = this.svg.append("g")
                        .attr('class',this.genome_name.toLowerCase())
                        .attr("transform", "translate(" + this.margin.left + "," + (this.margin.top + (this.height * this.id) + (10 * this.id)) + ")");

    this.textContainer = this.svg.append("g")
                        .attr('class',this.genome_name.toLowerCase())
                        .attr("transform", "translate(" + this.margin.left + "," + (this.margin.top + (this.height * this.id) + (10 * this.id)) + ")");

    /* Add both into the page */
    this.lineContainer.append("path")
        .data([this.variability_b])
        .attr("class", "line")
        .attr("name", "first_pos")
        .style("fill", '#990000')
        .attr("d", this.line);

    this.lineContainer.append("path")
        .data([this.variability_c])
        .attr("class", "line")
        .attr("name", "second_pos")
        .style("fill", '#990000')
        .attr("d", this.line);

    this.lineContainer.append("path")
        .data([this.variability_d])
        .attr("class", "line")
        .attr("name", "third_pos")
        .style("fill", '#004400')
        .attr("d", this.line);

    this.lineContainer.append("path")
        .data([this.variability_a])
        .attr("class", "line")
        .attr("name", "outside_gene")
        .style("stroke", '#666666')
        .style("stroke-width", "0.2")
        .attr("d", this.line);



    
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
                   .text(this.genome_name);
    
}
