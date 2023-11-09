/**
 *  Dialog for loading collections.
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


function StoreCollectionDialog() {
    this.dialog = document.createElement('div');
    this.dialog.setAttribute('class', 'modal fade in');

    this.dialog.innerHTML = `<div class="modal-dialog">
            <div class="modal-content">
                <div class="modal-header">
                    <h4 class="modal-title">Store bin collection</h4>
                    <button class="close" data-dismiss="modal" type="button"><span>&times;</span></button>
                </div>

                <div class="modal-body">
                    <div class="form-group">
                        <select class="form-control collection-list" size="10"></select>
                    </div>
                    <div class="form-group">
                            <label class="col-md-2 settings-label label-big">Name: </label>  
                            <div class="col-md-10">
                                <input type="text" class="form-control collection-name" value="default" />
                            </div>
                    </div>
                </div>

                <div class="modal-footer">
                    <button type="button" class="btn btn-primary">Store</button>
                    <button class="btn btn-default" data-dismiss="modal" type="button">Close</button>
                </div>
            </div>
        </div>`;

    this.dialog.querySelector('.collection-list').addEventListener('change', () => { this.UpdateName(); });
    this.dialog.querySelector('.collection-list').addEventListener('dblclick', () => { this.StoreCollection(); });
    this.dialog.querySelector('.btn-primary').addEventListener('click', () => { this.StoreCollection(); });
};


StoreCollectionDialog.prototype.UpdateName = function() {
    this.dialog.querySelector('.collection-name').value = this.dialog.querySelector('.collection-list').value;
}


StoreCollectionDialog.prototype.Show = async function() {
    this.collections = await $.ajax({
        type: 'GET', 
        url: '/data/collections', 
        cache: false
    });

    for (let source in this.collections) {
        let read_only = this.collections[source]['read_only'];
        let name = read_only ? source + ' (read only)' : source;
        let disabled = read_only ? ' disabled="disabled"' : '';

        this.dialog.querySelector('.collection-list').innerHTML += `<option value="${source}" ${disabled}>${name}</option>`;
    }

    $(this.dialog).modal('show').on('hidden.bs.modal', () => this.dialog.remove());
};


StoreCollectionDialog.prototype.StoreCollection = async function() {
    let collection_name = this.dialog.querySelector('.collection-name').value.replace(/\W+/g, "_");
    
    if (collection_name == "") {
        toastr.warning('Please enter collection name or select one.');
        this.dialog.querySelector('.collection-name').focus();
    }

    if (this.collections.hasOwnProperty(collection_name) && this.collections[collection_name]['read_only']) {
        toastr.warning('This collection is read only.');
    }

    let collection_info = bins.ExportCollection();

    $.ajax({
        type: 'POST',
        url: '/store_collection',
        data: {
            source: collection_name,
            data: JSON.stringify(collection_info['data'], null, 4),
            colors: JSON.stringify(collection_info['colors'], null, 4),
        },
        success: (response) => {
            toastr.info(response, "Server");

            $(this.dialog).modal('hide');
            this.dialog.remove();
        }
    });
};