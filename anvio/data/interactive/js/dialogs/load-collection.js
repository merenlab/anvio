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


function LoadCollectionDialog() {
    this.dialog = document.createElement('div');
    this.dialog.setAttribute('class', 'modal fade in');

    this.dialog.innerHTML = `<div class="modal-dialog">
            <div class="modal-content">
                <div class="modal-header">
                    <button class="close" data-dismiss="modal" type="button"><span>&times;</span></button>
                    <h4 class="modal-title">Load bin collection</h4>
                </div>
                <div class="modal-body">
                    <div class="col-md-7">
                        <select class="form-control collection-list" size="10"></select>
                    </div>
                    <div class="col-md-5 collection-details">
                            <div class="col-md-12">Collection Details</div>
                            <hr>
                            <div class="col-md-8">Number of Items:</div><div class="col-md-4">n/a</div>
                            <div class="col-md-8">Number of Bins:</div><div class="col-md-4">n/a</div>
                    </div>
                    <div class="col-md-6 col-md-offset-3 form-inline" style="${mode != 'full' ? 'display:none' : ''}">
                        <br /> Minimum bin size: 
                        <input type="text" value="1" size="4" class="form-control input-xs threshold-value">
                        <select class="form-control input-xs threshold-base">
                            <option value="1000">K</option>
                            <option value="1000000">M</option>
                            <option value="1000000000">G</option>
                        </select> bases.
                    </div>
                </div>

                <div class="modal-footer">
                    <button type="button" class="btn btn-primary">Load</button>
                    <button class="btn btn-default" data-dismiss="modal" type="button">Close</button>
                </div>
            </div>
        </div>`;

    this.dialog.querySelector('.collection-list').addEventListener('change', (event) => { this.UpdateDetails(); });
    this.dialog.querySelector('.collection-list').addEventListener('dblclick', (event) => { this.LoadCollection(); });
    this.dialog.querySelector('.btn-primary').addEventListener('click', (event) => { this.LoadCollection(); });
};


LoadCollectionDialog.prototype.Show = async function() {
    this.collections = await $.ajax({
        type: 'GET', 
        url: '/data/collections', 
        cache: false
    });

    for (let source in this.collections) {
        let read_only = this.collections[source]["read_only"];
        let name = (read_only) ? source + ' (read only)' : source;

        this.dialog.querySelector('.collection-list').innerHTML += `<option value="${source}">${name}</option>`;
    }

    $(this.dialog).modal('show').on('hidden.bs.modal', () => this.dialog.remove());
};


LoadCollectionDialog.prototype.UpdateDetails = function() {
    let collection_name = this.dialog.querySelector('.collection-list').value;
    if (collection_name == "")
        return;

    let collection_info = this.collections[collection_name];
    this.dialog.querySelector('.collection-details').innerHTML = `<div class="col-md-12">Collection Details</div><hr>
        <div class="col-md-8">Number of Items:</div><div class="col-md-4"><b>${collection_info['num_splits']}</b></div>
        <div class="col-md-8">Number of Bins:</div><div class="col-md-4"><b>${collection_info['num_bins']}</b></div>`;
};


LoadCollectionDialog.prototype.LoadCollection = function() {
    let collection_name = this.dialog.querySelector('.collection-list').value;
    
    if (collection_name == "") {
        toastr.warning('Please select a collection.');
    }

    if (bins.GetTotalSelectionCount() > 0 && !confirm("You are trying to load a collection but you already have some selections,\
                                                       you will lose them if they are not stored. Do you want to continue?")) {
        return;
    }

    let threshold_value = this.dialog.querySelector('.threshold-value').value;
    let threshold_base  = this.dialog.querySelector('.threshold-base').value;
    let threshold = parseInt(threshold_value) * parseInt(threshold_base);

    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/collection/' + collection_name,
        success: function(data) {
            bins.ImportCollection(data, threshold);
            
            $(this.dialog).modal('hide');
            this.dialog.remove();
        }.bind(this)
    });
};

