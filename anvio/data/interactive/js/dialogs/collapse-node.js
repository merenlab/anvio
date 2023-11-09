/**
 *  Dialog for collapsing nodes
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


function CollapseNodeDialog(node) {
    this.node = node;

    this.dialog = document.createElement('div');
    this.dialog.setAttribute('class', 'modal fade in');

    this.dialog.innerHTML = `<div class="modal-dialog">
            <div class="modal-content">
                <div class="modal-header">
                    <h4 class="modal-title">Collapse Node</h4>
                    <button class="close" data-dismiss="modal" type="button"><span>&times;</span></button>
                </div>

                <div class="modal-body">
                    <div class="col-md-12">
                        <label class="col-md-4 settings-label">Label</label>  
                        <div class="col-md-8">
                            <input type="text" id="collapsed_node_label" value="Collapsed Node ${collapsedNodes.length+1}">
                        </div>
                    </div>
                    <div class="col-md-12">
                        <label class="col-md-4 settings-label">Font size</label>  
                        <div class="col-md-8">
                            <input type="text" id="collapsed_font_size" value="0"> (0: auto)
                        </div>
                    </div>
                    <div class="col-md-12">
                        <label class="col-md-4 settings-label">Font color</label>  
                        <div class="col-md-8">
                            <div class="colorpicker" color="#888888" style="background-color: #888888;"></div>
                        </div>
                    </div>
                    <div class="col-md-12">
                        <label class="col-md-4 settings-label">Reduce branch size to</label>  
                        <div class="col-md-8">
                            <input type="text" id="collapsed_node_size" value="0.25">
                        </div>
                    </div>
                </div>

                <div class="modal-footer">
                    <button class="btn btn-primary" type="button">Collapse</button>
                    <button class="btn btn-default" data-dismiss="modal" type="button">Close</button>
                </div>
            </div>
        </div>`;

    $(this.dialog.querySelector('.colorpicker')).colpick({
        layout: 'hex',
        submit: 0,
        colorScheme: 'light',
        onChange: function(hsb, hex, rgb, el, bySetColor) {
            $(el).css('background-color', '#' + hex);
            $(el).attr('color', '#' + hex);

            if (!bySetColor) $(el).val(hex);
        }
    }).keyup(function() {
        $(this).colpickSetColor(this.value);
    });

    this.dialog.querySelector('.btn-primary').addEventListener('click', (event) => { this.Collapse(); });
    $(this.dialog).modal('show').on('hidden.bs.modal', () => this.dialog.remove());
};


CollapseNodeDialog.prototype.IsNameExists = function(name) {
    if (drawer.tree.label_to_leaves.hasOwnProperty(name)) {
        return true;
    }

    collapsedNodes.forEach((entry) => {
        if (entry['name'] == name) {
            return true;
        }
    });
    return false;
}


CollapseNodeDialog.prototype.Collapse = function() {
    let label = this.dialog.querySelector('#collapsed_node_label').value;
    let size = this.dialog.querySelector('#collapsed_node_size').value;
    let font_size = this.dialog.querySelector('#collapsed_font_size').value;
    let color = this.dialog.querySelector('.colorpicker').getAttribute('color');

    if (this.IsNameExists(label)) {
        toastr.warning("This node label already exists in tree.", "Collapse", { 'timeOut': '0', 'extendedTimeOut': '0' });
        return;
    }

    let [left_most, right_most] = this.node.GetBorderNodes();

    collapsedNodes.push({
        'left_most': left_most.label,
        'right_most': right_most.label,
        'label': label,
        'size': size,
        'font_size': font_size,
        'color': color
    });

    $('#tree_modified_warning').show();
    $(this.dialog).modal('hide');
    this.dialog.remove();

    drawTree();
};

