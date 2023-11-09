/**
 *  Dialog for saving modified trees.
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


function SaveTreeDialog(tree_type) {
    if (!(mode == 'full' || mode == 'pan' || mode == 'manual')) {
        toastr.warning('Saving modified tree functionality is not available for this mode.');
        return;
    }

    if (typeof tree_type === 'undefined') {
        tree_type = 'items';
    }

    this.tree_type = tree_type

    this.current_tree_name = last_settings['order-by'];

    this.dialog = document.createElement('div');
    this.dialog.setAttribute('class', 'modal fade in');

    this.dialog.innerHTML = `<div class="modal-dialog">
            <div class="modal-content">
                <div class="modal-header">
                    <h4 class="modal-title">Save Tree ${this.tree_type=='samples' ? ` (Layers) ` : ``}</h4>
                    <button class="close" data-dismiss="modal" type="button"><span>&times;</span></button>
                </div>

                <div class="modal-body">
                    <div class="col-md-12">
                        <label class="col-md-4 settings-label"><input type="radio" name="overwrite[]" value="no" checked="checked">Create new tree with name</label>  
                        <div class="col-md-8">
                            <input type="text" id="tree_name" value="New tree">
                        </div>
                    </div>
                    ${this.tree_type!='samples' ? `
                    <div class="col-md-12">
                        <label class="col-md-4 settings-label"><input type="radio" name="overwrite[]" value="yes">Overwrite the current tree</label>  
                        <div class="col-md-8">
                            ${getClusteringPrettyName(this.current_tree_name)}
                        </div>
                    </div>
                    ` : ``}
                </div>

                <div class="modal-footer">
                    <button class="btn btn-primary" type="button">Save</button>
                    <button class="btn btn-default" data-dismiss="modal" type="button">Close</button>
                </div>
            </div>
        </div>`;

    this.dialog.querySelector('.btn-primary').addEventListener('click', (event) => { this.SaveTree(); });
    this.dialog.querySelectorAll('input[type="radio"]').forEach((input) => { 
        input.addEventListener('click', (event) => { 
            if (this.dialog.querySelectorAll('input[type="radio"]')[0].checked) {
                this.dialog.querySelector('#tree_name').disabled = false;
            } else {
                this.dialog.querySelector('#tree_name').disabled = true;
            }
        })
    });

    $(this.dialog).modal('show').on('hidden.bs.modal', () => this.dialog.remove());
};


SaveTreeDialog.prototype.SaveTree = function() {
    let new_tree_name;
    let overwrite = !this.dialog.querySelectorAll('input[type="radio"]')[0].checked;

    if (overwrite) 
    {
        new_tree_name = this.current_tree_name;
    } else {
        if (this.tree_type == 'items') {
            let parts = this.current_tree_name.split(':');
            new_tree_name = this.dialog.querySelector('#tree_name').value + ':' + parts[1] + ':' + parts[2];
        }
        else
        {
            new_tree_name = this.dialog.querySelector('#tree_name').value;   
        }
    }

    $.ajax({
        type: 'POST',
        cache: false,
        url: '/data/save_tree',
        data: {
            'tree_type': this.tree_type,
            'overwrite': overwrite,
            'name': new_tree_name,
            'data': (this.tree_type == 'samples') ? samplesClusteringData['newick'] : clusteringData,
            'additional': JSON.stringify({'collapsedNodes': collapsedNodes})
        },
        success: function(data) {
            if (data['status'] == 0) {
                toastr.success(data['message'], "Server");

                if (this.tree_type == 'samples') {
                    samples_order_dict[new_tree_name] = {'newick': samplesClusteringData['newick'], 'basic': ''};
                    $('#samples_order').append(`<option value="${new_tree_name}">${new_tree_name}</option>`);
                    $('#samples_order').val(new_tree_name);
                    $('#samples_tree_modified_warning').hide();
                } else {
                    $('#tree_modified_warning').hide();
                    $('#trees_container').append(`<option value="${new_tree_name}">${getClusteringPrettyName(new_tree_name)}</option>`);
                    $('#trees_container').val(new_tree_name);
                }
                
                $(this.dialog).modal('hide');
                this.dialog.remove();
            } else {
                toastr.warning(data['message'], "Server", { 'timeOut': '0', 'extendedTimeOut': '0' });
            }
        }.bind(this)
    });
};

