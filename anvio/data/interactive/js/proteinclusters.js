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
var gene_caller_ids;
var gene_caller_ids_in_genomes;
var gene_aa_sequences_for_gene_caller_ids;
var previous_pc_name;
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
    // magic happens here.
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

