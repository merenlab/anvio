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
            'action': (node, layer, param) => {
                localStorage.state = JSON.stringify(serializeSettings(true), null, 4);
                window.open(generate_inspect_link('inspect_' + param, node.label), '_blank'); 
            }
        },
        'inspect_split': {
            'title': 'Inspect split',
            'action': (node, layer, param) => {
                this.menu_items['inspect']['action'](node, layer, 'split');
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
        'collapse': {
            'title': 'Collapse',
            'action': (node, layer, param) => {
                new_tree = new Tree();
                new_tree.Parse(clusteringData.trim(), false);
                new_tree.nodes[this.node.id].collapsed = true;
                clusteringData = new_tree.Serialize();
                $('#tree_modified_warning').show();
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

    }
}

ContextMenu.prototype.BuildMenu = function() {
    var menu = [];

    if (this.node.IsLeaf()) {
        if (bins.IsNodeMemberOfBin(this.node)) {
            menu.push('remove');
        } else {
            menu.push('select');
        }
        menu.push('divider');
        menu.push('select_layer');
        menu.push('unselect_layer');
        menu.push('divider');

        if (mode == 'gene') {
            menu.push('inspect_context');
            menu.push('inspect_gene');
        }
        else if (mode == 'pan') {
            menu.push('inspect_geneclusters');
        }
        else {
            menu.push('inspect_split');
        }

        menu.push('divider');
        menu.push('get_split_sequence');
        menu.push('blastn_nr');
        menu.push('blastx_nr');
        menu.push('blastn_refseq_genomic');
        menu.push('blastx_refseq_protein');
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
