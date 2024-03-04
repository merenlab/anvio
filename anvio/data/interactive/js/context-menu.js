/**
 *  Handles right click menu functions
 *
 *  Authors: Özcan Esen <ozcanesen@gmail.com>
 *           A. Murat Eren
 *           Matthew Klein <mtt.l.kln@gmail.com>
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

let outerLimit1;
let outerLimit2;

ContextMenu = function(options) {
    this.container = options.container;
    this.event = options.event;
    this.node = options.node;
    this.layer = options.layer;
    this.isSample = options.isSample;
    let all = options.all

    this.menu_items = {
        'select': {
            'title': 'Add item to active bin',
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
        'set_outer_limit_1' : {
            'title' : "Mark item as 'range start'",
            'action': (node, layer, param, all) => {
                bins.AppendNode(node); // TODO flag node in interface instead of appending
                outerLimit1 = node.order
            }
        },
        'remove_outer_limit_1' : {
            'title' : 'Remove marked range start',
            'action': (node, layer, param) => { removeOuterLimit1(all) }
        },
        'set_outer_limit_2_add' : {
            'title' : "Add items in 'range' to active bin",
            'action': (node, layer, param) => {
                bins.AppendNode(node);
                outerLimit2 = node.order
                setBinningRange(outerLimit1, outerLimit2, all, 'add')
            }
        },
        'set_outer_limit_2_remove' : {
            'title' : "Remove items in 'range' from any bin",
            'action': (node, layer, param) => {
                bins.AppendNode(node);
                outerLimit2 = node.order
                setBinningRange(outerLimit1, outerLimit2, all, 'remove')
            }
        },
        'select_layer': {
            'title': 'Select layer in the main panel',
            'action': (node, layer, param) => {
                $(`#tbody_layers tr:nth-child(${layer}) input:checkbox`).prop('checked', true);
            }
        },
        'unselect_layer': {
            'title': 'Unselect layer in the main panel',
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
            'title': 'Inspect gene and context',
            'action': (node, layer, param) => {
                this.menu_items['inspect']['action'](node, layer, 'context');
            }
        },
        'inspect_gene': {
            'title': 'Inspect gene alone',
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
        'open_in_legends': {
            'title': 'Change category color(s)',
            'action': (node, layer, param) => {
                if (!$('#panel-left').is(':visible')) {
                    toggleLeftPanel()
                }

                let splitHTML = layerdata_title[node.label][layer[0] -1].split('</td>')
                let legend = extractContent(splitHTML[0]).toLowerCase()
                let query = extractContent(splitHTML[1])
                let legendIndex = Number()

                var all = $(".ui-accordion-header").map(function() {
                    return this.innerHTML;
                }).get();
                all.map((i, idx) => {
                    if(i.toLowerCase().includes(legend.replaceAll('_', ' '))){
                        legendIndex = idx
                    }
                })
                $( ".ui-accordion" ).accordion( "option", "active", legendIndex );
                $( `#${legend.replaceAll('_','-').replaceAll(' ','-')}-query-input`).val(query)

                function extractContent(s) { // reference https://stackoverflow.com/questions/28899298/extract-the-text-out-of-html-string-using-javascript
                    var span = document.createElement('span');
                    span.innerHTML = s;
                    return span.textContent || span.innerText;
                };
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
            'title': 'Collapse these items',
            'action': (node, layer, param) => {
                new CollapseNodeDialog(node);
            }
        },
        'expand': {
            'title': 'Expand this collapsed node',
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
            'title': 'Rotate the tree/dendrogram here',
            'action': (node, layer, param) => {
                new_tree = new Tree();
                new_tree.Parse(clusteringData.trim(), false);
                new_tree.nodes[this.node.id].Rotate();
                clusteringData = new_tree.Serialize();
                $('#tree_modified_warning').show();
                drawTree();
            }
        },
        'reroot_disabled': {
            'title': 'Reroot the tree/dendogram here'
        },
        'reroot_with_internal_node': {
            'title': 'Tree has internal node names',
            'action': (node, layer, param) => {
                let [left_most, right_most] = this.node.GetBorderNodes();

                $.ajax({
                    type: 'POST',
                    cache: false,
                    url: '/data/reroot_tree',
                    data: {
                        'newick': clusteringData,
                        'left_most': left_most.label,
                        'right_most': right_most.label,
                        'internal_node': internal_node_value
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
        'reroot_with_support_values': {
            'title': 'Tree has branch support values',
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
                let target = (mode == 'gene') ? 'gene' : 'contig';
                $.ajax({
                    type: 'GET',
                    cache: false,
                    url: '/data/'+ target + '/' + node.label,
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
                search_gene_sequence_in_remote_dbs(node.label, 'blastn', 'nr', (mode == 'gene') ? 'gene' : 'contig');
             }
        },
        'blastx_nr': {
            'title': ' - blastx @ nr',
            'action': (node, layer, param) => {
                search_gene_sequence_in_remote_dbs(node.label, 'blastx', 'nr', (mode == 'gene') ? 'gene' : 'contig');
             }
        },
        'blastn_refseq_genomic': {
            'title': ' - blastn @ refseq_genomic',
            'action': (node, layer, param) => {
                search_gene_sequence_in_remote_dbs(node.label, 'blastn', 'refseq_genomic', (mode == 'gene') ? 'gene' : 'contig');
             }
        },
        'blastx_refseq_protein': {
            'title': ' - blastx @ refseq_protein',
            'action': (node, layer, param) => {
                search_gene_sequence_in_remote_dbs(node.label, 'blastx', 'refseq_protein', (mode == 'gene') ? 'gene' : 'contig');
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

function setBinningRange(limit1, limit2, all, action) {
    limit1 > limit2 ? countDown() : countUp()
    function countDown() {
        for(let i = limit1 - 1; i > limit2; i--){
            action === 'add' ? bins.AppendNode([all[i]]) : bins.RemoveNode([all[i]])
        }
    }
    function countUp() {
        for (let i = limit1 + 1; i < limit2; i++){
            action === 'add' ? bins.AppendNode([all[i]]) : bins.RemoveNode([all[i]])
        }
    }
    if(action === 'remove'){
        bins.RemoveNode([all[outerLimit1]])
        bins.RemoveNode([all[outerLimit2]])
    }
    outerLimit1 = null
    outerLimit2 = null
}

removeOuterLimit1 = (all) => {
    bins.RemoveNode(all[outerLimit1])
    outerLimit1 = null
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
            // start with menu items for inspection pages
            // as they are most frequently used
            if (mode == 'gene' && inspection_available) {
                menu.push('inspect_context');
                menu.push('inspect_gene');
                menu.push('inspect_context_quick');
                menu.push('inspect_gene_quick');
                menu.push('divider');
            }
            else if (mode == 'pan') {
                menu.push('inspect_geneclusters');
                menu.push('divider');
            }
            else if (inspection_available && mode != 'collection') {
                menu.push('inspect_split');
                menu.push('inspect_split_quick');
                menu.push('divider');
            }

            // next, we have adding/removing items to bins individually
            // or as a range
            if (bins.IsNodeMemberOfBin(this.node)) {
                menu.push('remove');
            } else {
                menu.push('select');
            }

            let limit1Exists;
            let currentBin = $('input[name="active_bin"]:checked').val();

            bins.selections[`${currentBin}`].forEach(node => { // iterate through bin selections to make sure limit1 hasnt been removed elsewhere
                if(node.order === outerLimit1){
                    limit1Exists = true
                    return
                }
            })
            if(limit1Exists){
                menu.push('set_outer_limit_2_add')
                menu.push('set_outer_limit_2_remove')
                menu.push('remove_outer_limit_1')
            } else {
                menu.push('set_outer_limit_1')
            }

            menu.push('divider');

            // menu items to select/unselect layers in the Settings tab
            if (this.layer) {
                menu.push('select_layer');
                menu.push('unselect_layer');

                // only show 'open in legends' option when target layer contains categorical data
                let splitHTML = layerdata_title[this.node.label][this.layer[0] -1].split('</td>')
                let legend = extractContent(splitHTML[0]).toLowerCase()
                let matchedLayerIndex = layerdata[0].map(l => l.toLowerCase()).indexOf(legend).toString()
                if(Object.keys(categorical_data_colors).includes(matchedLayerIndex)){
                    menu.push('divider');
                    menu.push('open_in_legends')
                }
                menu.push('divider');

                function extractContent(s) { // reference https://stackoverflow.com/questions/28899298/extract-the-text-out-of-html-string-using-javascript
                    var span = document.createElement('span');
                    span.innerHTML = s;
                    return span.textContent || span.innerText;
                };
            }

            // getting back sequences for things
            if (mode == 'pan') {
                menu.push('get_AA_sequences_for_gene_cluster');
            }
            else if (mode == 'collection') {
                menu.push('hmm_RecA');
                menu.push('hmm_Ribosomal_S3_C');
                menu.push('hmm_Ribosomal_L10');
            }
            else if (sequences_available) {
                menu.push('get_split_sequence');
                menu.push('blastn_nr');
                menu.push('blastx_nr');
                menu.push('blastn_refseq_genomic');
                menu.push('blastx_refseq_protein');
            }

            else if (mode == 'manual'){
                menu.push('reroot_disabled');
                menu.push('reroot_with_internal_node');
                menu.push('reroot_with_support_values');
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

                menu.push('divider');
                menu.push('reroot_disabled');
                menu.push('reroot_with_internal_node');
                menu.push('reroot_with_support_values');
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
    list.style.visibility = 'hidden';
    list.style.display = 'block';

    for (const item of this.BuildMenu()) {
        if (item == 'divider') {
            list.innerHTML += `<li class="dropdown-divider"></li>`;
        } else if (item.includes('disabled')){
            list.innerHTML += `<li><a class="dropdown-item disabled bold" href="#">${this.menu_items[item]['title']}</a></li>`;
        } else {
            list.innerHTML += `<li><a class="dropdown-item" href="#" item-name="${item}">${this.menu_items[item]['title']}</a></li>`;
        }
    }

    this.container.appendChild(list);
    list.style.left = Math.min(VIEWER_WIDTH - list.clientWidth, this.event.clientX) + 'px';
    list.style.top = Math.min(VIEWER_HEIGHT - list.clientHeight, this.event.clientY) + 'px';
    list.style.visibility = '';

    list.addEventListener('click', (event) => {
        let item_name = event.target.getAttribute('item-name');

        this.menu_items[item_name]['action'](this.node, this.layer);
    });

    this.container.addEventListener('click', (event) => {
        list.remove();
    }, {once: true});
};
