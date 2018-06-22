/**
 *  Handles right click menu
 *
 *  Author: Özcan Esen <ozcanesen@gmail.com>
 *  Credits: A. Murat Eren
 *  Copyright 2018, The anvio Project
 *
 * This file is part of anvi'o (<https://github.com/merenlab/anvio>).
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


ContextMenu = function(options) {
    this.container = options.container;
    this.event = event;
    this.node = options.node;
    this.layer = options.layer;
    this.isSample = options.isSample;

    this.menu_items = {
        'select': {
            'title': 'Add item to bin',
            'action': (node, layer, param) => {
                bins.AppendNode(node);
            }
        },
        'remove': {
            'title': 'Remove item from bin',
            'action': (node, layer, param) => {
                bins.RemoveNode(node);
            }
        },
        'select_layer': {
            'title': 'Select layer',
            'action': (node, layer, param) => {
                $(`#tbody_layers tr:nth-child(${layer}) input:checkbox`).prop('checked', true);
            }
        },
        'unselect_layer': {
            'title': 'Unselect layer',
            'action': (node, layer, param) => {
                $(`#tbody_layers tr:nth-child(${layer}) input:checkbox`).prop('checked', false);
            }
        },
        'inspect': {
            'title': 'Inspect',
            'action': (node, layer, param, show_snvs) => {
                if (typeof show_snvs === 'undefined') {
                    show_snvs = true;
                }

                localStorage.state = JSON.stringify(serializeSettings(true), null, 4);
                window.open(generate_inspect_link({'type': 'inspect_' + param, 'item_name': node.label, 'show_snvs': show_snvs}), '_blank'); 
            }
        },
        'inspect_split': {
            'title': 'Inspect split',
            'action': (node, layer, param) => {
                this.menu_items['inspect']['action'](node, layer, 'split');
            }
        },
        'inspect_split_quick': {
            'title': 'Inspect split (no SNVs)',
            'action': (node, layer, param) => {
                this.menu_items['inspect']['action'](node, layer, 'split', false);
            }
        },
        'inspect_geneclusters': {
            'title': 'Inspect gene cluster',
            'action': (node, layer, param) => {
                this.menu_items['inspect']['action'](node, layer, 'geneclusters');
            }
        },
        'inspect_context': {
            'title': 'Inspect context',
            'action': (node, layer, param) => {
                this.menu_items['inspect']['action'](node, layer, 'context');
            }
        },
        'inspect_gene': {
            'title': 'Inspect gene',
            'action': (node, layer, param) => {
                this.menu_items['inspect']['action'](node, layer, 'gene');
            }
        },
        'inspect_context_quick': {
            'title': 'Inspect context (no SNVs)',
            'action': (node, layer, param) => {
                this.menu_items['inspect']['action'](node, layer, 'context', false);
            }
        },
        'inspect_gene_quick': {
            'title': 'Inspect gene (no SNVs)',
            'action': (node, layer, param) => {
                this.menu_items['inspect']['action'](node, layer, 'gene', false);
            }
        },
        'get_hmm_sequence': {
            'title': 'Inspect gene',
            'action': (node, layer, param) => {
                $.ajax({
                    type: 'GET',
                    cache: false,
                    url: '/data/hmm/' + node.label + '/' + param,
                    success: function(data) {
                        if ('error' in data){
                            $('#modGenerateSummary').modal('hide');
                            waitingDialog.hide();
                            toastr.error(data['error'], "", { 'timeOut': '0', 'extendedTimeOut': '0' });
                        } else {
                            $('#splitSequence').val('>' + data['header'] + '\n' + data['sequence']);
                            $('#modSplitSequence').modal('show');
                        }
                    }
                });
            }
        },
        'get_gene_sequence': {
            'title': 'Get gene sequence',
            'action': (node, layer, param) => {
                $.ajax({
                    type: 'GET',
                    cache: false,
                    url: '/data/gene/' + node.label,
                    success: function(data) {
                        $('#modSplitSequence .modal-title').html('Gene Sequence');
                        $('#splitSequence').val('>' + data['header'] + '\n' + data['sequence']);
                        $('#modSplitSequence').modal('show');
                    }
                });
            }
        },
        'collapse': {
            'title': 'Collapse',
            'action': (node, layer, param) => {
                new CollapseNodeDialog(node);
            }
        },
        'expand': {
            'title': 'Expand',
            'action': (node, layer, param) => {
                for (let i=0; i < collapsedNodes.length; i++) {
                    if (collapsedNodes[i]['label'] == node.label) {
                        collapsedNodes.splice(i, 1);
                        break;
                    }
                }
                drawTree();
            }
        }, 
        'rotate': {
            'title': 'Rotate',
            'action': (node, layer, param) => {
                new_tree = new Tree();
                new_tree.Parse(clusteringData.trim(), false);
                new_tree.nodes[this.node.id].Rotate();
                clusteringData = new_tree.Serialize();
                $('#tree_modified_warning').show();
                drawTree();
            }
        },
        'reroot': {
            'title': 'Reroot',
            'action': (node, layer, param) => {
                let [left_most, right_most] = this.node.GetBorderNodes();

                $.ajax({
                    type: 'POST',
                    cache: false,
                    url: '/data/reroot_tree',
                    data: {
                        'newick': clusteringData,
                        'left_most': left_most.label,
                        'right_most': right_most.label  
                    },
                    success: function(data) {
                        collapsedNodes = [];
                        clusteringData = data['newick'];
                        $('#tree_modified_warning').show();
                        drawTree();
                    }
                });
            }
        },
        'get_split_sequence': {
            'title': 'Get split sequence',
            'action': (node, layer, param) => {
                $.ajax({
                    type: 'GET',
                    cache: false,
                    url: '/data/contig/' + node.label,
                    success: function(data) {
                        $('#modSplitSequence .modal-title').html('Split Sequence');
                        $('#splitSequence').val('>' + data['header'] + '\n' + data['sequence']);
                        $('#modSplitSequence').modal('show');
                    }
                });
            }
        },
        'blastn_nr': {
            'title': ' - blastn @ nr',
            'action': (node, layer, param) => { 
                get_sequence_and_blast(node.label, 'blastn', 'nr', (mode == 'gene') ? 'gene' : 'contig');
             }
        },
        'blastx_nr': {
            'title': ' - blastx @ nr',
            'action': (node, layer, param) => { 
                get_sequence_and_blast(node.label, 'blastx', 'nr', (mode == 'gene') ? 'gene' : 'contig');
             }
        },
        'blastn_refseq_genomic': {
            'title': ' - blastn @ refseq_genomic',
            'action': (node, layer, param) => { 
                get_sequence_and_blast(node.label, 'blastn', 'refseq_genomic', (mode == 'gene') ? 'gene' : 'contig');
             }
        },
        'blastx_refseq_protein': {
            'title': ' - blastx @ refseq_protein',
            'action': (node, layer, param) => { 
                get_sequence_and_blast(node.label, 'blastx', 'refseq_protein', (mode == 'gene') ? 'gene' : 'contig');
             }
        },
        'samples_rotate': {
            'title': 'Rotate',
            'action': (node, layer, param) => {
                new_tree = new Tree();
                new_tree.Parse(samplesClusteringData['newick'], false);
                new_tree.nodes[this.node.id].Rotate();
                samplesClusteringData['newick'] = new_tree.Serialize();
                $('#samples_order').val('custom').trigger('change');
                $('#samples_tree_modified_warning').show();
                drawTree();
            }
        },
        'samples_reroot': {
            'title': 'Reroot',
            'action': (node, layer, param) => { 
                let [left_most, right_most] = this.node.GetBorderNodes();

                $.ajax({
                    type: 'POST',
                    cache: false,
                    url: '/data/reroot_tree',
                    data: {
                        'newick': samplesClusteringData['newick'],
                        'left_most': left_most.label,
                        'right_most': right_most.label  
                    },
                    success: function(data) {
                        samplesClusteringData['newick'] = data['newick'];
                        $('#samples_order').val('custom').trigger('change');
                        $('#samples_tree_modified_warning').show();
                        drawTree();
                    }
                });
            }
        },
        'get_AA_sequences_for_gene_cluster': {
            'title': 'Get AA sequences',
            'action': (node, layer, param) => {
                $.ajax({
                    type: 'GET',
                    cache: false,
                    url: '/data/get_AA_sequences_for_gene_cluster/' + node.label,
                    success: function(data) {
                        var output = '';

                        for (var key in data)
                            output = output + ">" + key + "\n" + data[key] + "\n";

                        $('#splitSequence').val(output);
                        $('#modSplitSequence').modal('show');
                    }
                });
            }
        },
        'hmm_RecA': {
            'title': 'RecA',
            'action': (node, layer, param) => { this.menu_items['get_hmm_sequence']['action'](node, layer, 'RecA'); }
        },
        'hmm_Ribosomal_S3_C': {
            'title': 'Ribosomal S3',
            'action': (node, layer, param) => { this.menu_items['get_hmm_sequence']['action'](node, layer, 'Ribosomal_S3_C'); }
        },
        'hmm_Ribosomal_L10': {
            'title': 'Ribosomal L10',
            'action': (node, layer, param) => { this.menu_items['get_hmm_sequence']['action'](node, layer, 'Ribosomal_L10'); }
        }
    }
}

ContextMenu.prototype.BuildMenu = function() {
    var menu = [];

    if (this.isSample) {
        menu.push('samples_rotate');
        menu.push('samples_reroot');
    }
    else
    {
        if (this.node.IsLeaf() && !this.node.collapsed) {
            if (bins.IsNodeMemberOfBin(this.node)) {
                menu.push('remove');
            } else {
                menu.push('select');
            }
            menu.push('divider');
            menu.push('reroot');
            menu.push('divider');

            if (this.layer) {
                menu.push('select_layer');
                menu.push('unselect_layer');
                menu.push('divider');       
            }

            if (mode == 'gene') {
                menu.push('inspect_context');
                menu.push('inspect_gene');
                menu.push('inspect_context_quick');
                menu.push('inspect_gene_quick');
            }
            else if (mode == 'pan') {
                menu.push('inspect_geneclusters');
            }
            else {
                menu.push('inspect_split');
                menu.push('inspect_split_quick');
            }

            menu.push('divider');
            if (mode == 'pan') {
                menu.push('get_AA_sequences_for_gene_cluster');
            }
            else if (mode == 'collection') {
                menu.push('hmm_RecA');
                menu.push('hmm_Ribosomal_S3_C');
                menu.push('hmm_Ribosomal_L10');
            }
            else {
                menu.push('get_split_sequence');
                menu.push('blastn_nr');
                menu.push('blastx_nr');
                menu.push('blastn_refseq_genomic');
                menu.push('blastx_refseq_protein');
            }
        }
        else
        {
            if (this.node.collapsed) {
                menu.push('expand');
            }
            else {
                menu.push('select');
                menu.push('remove');

                menu.push('divider');
                
                menu.push('collapse');
                menu.push('rotate');
                menu.push('reroot');
            }
        }

    }
    
    return menu;
};


ContextMenu.prototype.Show = function() {
    this.container.querySelectorAll('.context-menu').forEach((menu) => {
        menu.remove();
    });

    var list = document.createElement('ul');
    list.setAttribute('class', 'dropdown-menu context-menu');
    list.setAttribute('role', 'menu');
    list.style.display = 'block';
    list.style.left = this.event.clientX + 'px';
    list.style.top = this.event.clientY + 'px';

    for (const item of this.BuildMenu()) {
        if (item == 'divider') {
            list.innerHTML += `<li class="divider"></li>`;
        } else {
            list.innerHTML += `<li><a href="#" item-name="${item}">${this.menu_items[item]['title']}</a></li>`;
        }
    }

    this.container.appendChild(list);

    list.addEventListener('click', (event) => {
        let item_name = event.target.getAttribute('item-name');

        this.menu_items[item_name]['action'](this.node, this.layer);
    });

    this.container.addEventListener('click', (event) => {
        list.remove();
    }, {once: true});
};
