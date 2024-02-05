/**
 * Draw bins, bin labels stuff.
 *
 *  Authors: Özcan Esen <ozcanesen@gmail.com>
 *           A. Murat Eren <a. murat.eren@gmail.com>
 *
 *  Copyright 2015-2021, The anvi'o project (http://anvio.org)
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

 const MAX_HISTORY_SIZE = 50;


function Bins(prefix, container) {
    this.selections = {}
    this.bin_counter = 0;
    this.prefix = prefix || "Bin_";
    this.higlighted_items = [];
    this.container = container || document.createElement("div");

    this.cache = {
        'completeness': {},
        'taxonomy': {},
    };

    this.keepHistory = false;
    this.allowRedraw = true;

    this.history = [];
    this.future = [];

    document.body.addEventListener('bin-settings-changed', (event) => this.RedrawBins());
    $(document).on('sorted', () => { this.BinsSorted(); });
};


Bins.prototype.NewBin = function(id, binState) {
    if (typeof id === 'undefined')
    {
        var from_state = false;
        var id = this.bin_counter;
        var name = this.prefix + (id + 1);
        var color = randomColor({luminosity: 'dark'});
        var contig_count = 0;
        var contig_length = "N/A";
        var num_gene_clusters = '---';
        var num_gene_calls = '---';
        var completeness = '---';
        var redundancy = '---';

        this.selections[id] = new Set();
        this.bin_counter++;
    }
    else
    {
        // we are adding bins from collection
        var from_state = true;
        var name = binState['name'];
        var color = binState['color'];
        var contig_count = 0;
        var contig_length = "N/A";
        var num_gene_clusters = "---";
        var num_gene_calls = "---";
        var completeness = "---";
        var redundancy = "---";
    }

    var template = `<tr bin-id="${id}" class="bin-row">
                       <td><input type="radio" name="active_bin" value="${id}"></td>
                       <td><div id="bin_color_${id}" class="colorpicker" color="${color}" style="background-color: ${color}"></td>
                       <td data-value="${name}">
                            <input type="text" class="bin-name" oninput="this.value = event.target.value.replaceAll(' ', '_');" onChange="emit('bin-settings-changed'); this.parentNode.setAttribute('data-value', this.value);" size="21" id="bin_name_${id}" value="${name}">
                        </td>
                       ${mode != 'pan' ? `
                           <td data-value="${contig_count}" class="num-items"><input type="button" value="${contig_count}" title="Click for contig names" onClick="showContigNames(${id});"></td>
                           <td data-value="${contig_length}" class="length-sum"><span>${contig_length}</span></td>
                       ` : ''}
                       ${mode == 'pan' ? `
                            <td data-value="${num_gene_clusters}" class="num-gene-clusters"><input type="button" value="${num_gene_clusters}" title="Click for quick gene cluster summaries" onClick="showGeneClusterDetails(${id});"></td>
                            <td data-value="${num_gene_calls}" class="num-gene-calls"><input type="button" value="${num_gene_calls}"></td>
                       ` : `
                            <td data-value="${completeness}" class="completeness"><input type="button" value="${completeness}" title="Click for completeness table" onClick="showCompleteness(${id});"></td>
                            <td data-value="${redundancy}" class="redundancy"><input type="button" value="${redundancy}" title="Click for redundant hits" onClick="showRedundants(${id}); "></td>
                       `}
                       <td><center><span class="default-bin-icon bi bi-trash-fill fa-lg" aria-hidden="true" alt="Delete this bin" title="Delete this bin" onClick="bins.DeleteBin(${id});"></span></center></td>
                    </tr>
                    <tr style="${ $('#estimate_taxonomy').is(':checked') ? `` : `display: none;`}" data-parent="${id}">
                            <td style="border-top: 0px;">&nbsp;</td>
                            <td style="border-top: 0px;">&nbsp;</td>
                            <td colspan="6" style="border-top: 0px; padding-top: 0px;">
                                <span bin-id="${id}" class="taxonomy-name-label">N/A</span>
                            </td>
                    </tr>`;

    this.container.insertAdjacentHTML('beforeend', template);
    this.SelectLastRadio();

    $('#bin_color_' + id).colpick({
        layout: 'hex',
        submit: 0,
        colorScheme: 'light',
        onChange: function(hsb, hex, rgb, el, bySetColor) {
            $(el).css('background-color', '#' + hex);
            $(el).attr('color', '#' + hex);

            if (!bySetColor) $(el).val(hex);
        },
        onShow: function() {
            $(this).attr('color-before', $(this).attr('color'));
        },
        onHide: function(picker) {
            let el = $(picker).data('colpick').el;
            let current_color = $(el).attr('color');
            let previous_color = $(el).attr('color-before');

            if (current_color != previous_color) {
                emit('bin-settings-changed');
                bins.PushHistory([{'type': 'ChangeColor',
                                   'bin_id': id,
                                   'color-before': previous_color,
                                   'color': current_color}]);
            }
        }
    }).keyup(function() {
        $(this).colpickSetColor(this.value);
    });

    this.PushHistory([{'type': 'NewBin',
                       'name': name,
                       'bin_id': id,
                       'color': color}]);
};


Bins.prototype.SelectLastRadio = function() {
    let radios = this.container.querySelectorAll('input[name=active_bin]');
    radios[radios.length - 1].checked = true;
};


Bins.prototype.GetSelectedBinId = function() {
    return this.container.querySelector('input[name=active_bin]:checked').value;
};


Bins.prototype.GetTotalSelectionCount = function() {
    let sum = 0;
    for (let bin_id in this.selections) {
        sum += this.selections[bin_id].size;
    }

    return sum;
};


Bins.prototype.GetSelectedBinColor = function() {
    return this.GetBinColor(this.GetSelectedBinId());
};


Bins.prototype.GetBinColor = function(bin_id) {
    return this.container.querySelector('#bin_color_' + bin_id.toString()).getAttribute('color');
}


Bins.prototype.DeleteBin = function(bin_id, show_confirm=true) {
    if (show_confirm && !confirm('Are you sure?')) {
        return;
    }

    let transaction = [];

    for (const node of this.selections[bin_id].values()) {
        node.ResetColor();

        if (this.keepHistory) {
            transaction.push({'type': 'RemoveNode',
                             'bin_id': bin_id,
                             'node': node});
        }
    }

    transaction.push({'type': 'DeleteBin',
                       'bin_id': bin_id,
                       'name': document.getElementById('bin_name_' + bin_id).value,
                       'color': document.getElementById('bin_color_' + bin_id).getAttribute('color')});

    this.selections[bin_id].clear();
    this.PushHistory(transaction);

    let bin_row = this.container.querySelector(`tr[bin-id='${bin_id}']`);

    bin_row.nextElementSibling.remove(); // remove taxonomy row too.
    bin_row.remove();

    if (!this.container.querySelectorAll('*').length) {
        // No bins left
        this.NewBin();
    }

    if (!this.container.querySelector('input[name=active_bin]:checked')) {
        this.SelectLastRadio();
    }

    this.RedrawBins();
};


Bins.prototype.PushHistory = function(transaction) {
    if (!this.keepHistory)
        return;

    this.history.push(transaction);

    if (this.history.length > MAX_HISTORY_SIZE) {
        this.history.shift();
    }

    // adding something to history always clears the future
    this.future = [];
}


Bins.prototype.Undo = function() {
    let transaction = this.history.pop();
    if (transaction) {
        this.ProcessTransaction(transaction, reversed=true);
        this.future.push(transaction);
    }
    else {
        toastr.warning('Can\'t do undo, history is empty.');
    }
}

Bins.prototype.Redo = function() {
    let transaction = this.future.pop();

    if (transaction) {
        this.ProcessTransaction(transaction);
        this.history.push(transaction);
    } else {
        toastr.warning('Can\'t do redo, future is empty.');
    }
}

Bins.prototype.ProcessTransaction = function(transaction, reversed=false) {
    this.keepHistory = false;
    this.allowRedraw = false;

    let updated_bins = new Set();
    let removed_bins = new Set();

    for (var i = 0; i < transaction.length; i++) {
        let operation = transaction[i];

        updated_bins.add(operation.bin_id);

        if (reversed) {
            switch (operation.type) {
                case 'AppendNode':
                    this.RemoveNode(operation.node, operation.bin_id);
                    break;
                case 'RemoveNode':
                    this.AppendNode(operation.node, operation.bin_id);
                    break;
                case 'ChangeColor':
                    $('#bin_color_' + operation.bin_id).attr('color', operation['color-before']);
                    $('#bin_color_' + operation.bin_id).css('background-color', operation['color-before']);
                    break;
                case 'DeleteBin':
                    this.NewBin(operation.bin_id, {'name': operation.name,
                                                  'color': operation.color});
                    removed_bins.delete(operation.bin_id);
                    break;
                case 'NewBin':
                    this.DeleteBin(operation.bin_id, show_confirm=false);
                    removed_bins.add(operation.bin_id);
                    break;
            }
        }
        else {
            switch (operation.type) {
                case 'AppendNode':
                    this.AppendNode(operation.node, operation.bin_id);
                    break;
                case 'RemoveNode':
                    this.RemoveNode(operation.node, operation.bin_id);
                    break;
                case 'ChangeColor':
                    $('#bin_color_' + operation.bin_id).attr('color', operation.color);
                    $('#bin_color_' + operation.bin_id).css('background-color', operation.color);
                    break;
                case 'DeleteBin':
                    this.DeleteBin(operation.bin_id, show_confirm=false);
                    removed_bins.add(operation.bin_id);
                    break;
                case 'NewBin':
                    this.NewBin(operation.bin_id, {'name': operation.name,
                                                  'color': operation.color});
                    removed_bins.delete(operation.bin_id);
                    break;
            }
        }
    }

    let bins_to_update = [];
    for (const bin_id of updated_bins) {
        if (!removed_bins.has(bin_id)) {
            bins_to_update.push(bin_id);
        }
    }

    this.keepHistory = true;
    this.allowRedraw = true;
    this.RebuildIntersections();
    this.RedrawBins();
    this.UpdateBinsWindow(bins_to_update);
};


Bins.prototype.DeleteAllBins = function() {
    if (!confirm('Are you sure you want to remove all bins?')) {
        return;
    }

    for (let tr of this.container.querySelectorAll('tr.bin-row')) {
        this.DeleteBin(tr.getAttribute('bin-id'), false);
    }
};


Bins.prototype.AppendNode = function(targets, bin_id) {
    if (typeof bin_id === 'undefined') {
        var bin_id = this.GetSelectedBinId();
    }
    var bin_color = this.GetBinColor(bin_id);
    var bins_to_update = new Set();

    if (!Array.isArray(targets)) {
        targets = [targets];
    }

    let transaction = [];

    for (const target of targets) {
        if (target.collapsed)
            continue;

        for (const node of target.IterateChildren()) {
            for (let other_bin_id in this.selections) {
                // remove node from other bins except the current one
                if (other_bin_id == bin_id) {
                    continue;
                }

                if (this.selections[other_bin_id].has(node)) {
                    this.selections[other_bin_id].delete(node);
                    bins_to_update.add(other_bin_id);

                    if (this.keepHistory && node.IsLeaf()) {
                        transaction.push({'type': 'RemoveNode',
                                          'bin_id': other_bin_id,
                                          'node': node});
                    }
                }
            }

            if (!this.selections[bin_id].has(node)) {
                this.selections[bin_id].add(node);
                bins_to_update.add(bin_id);

                if (this.keepHistory && node.IsLeaf()) {
                    transaction.push({'type': 'AppendNode',
                                      'bin_id': bin_id,
                                      'node': node});
                }
            }

            node.SetColor(bin_color);
        }
    }

    bins_to_update = Array.from(bins_to_update);
    this.PushHistory(transaction);
    this.RedrawBins();
    this.UpdateBinsWindow(bins_to_update);
};


Bins.prototype.RemoveNode = function(targets, bin_id) {
    if (typeof bin_id === 'undefined') {
        var bin_id = this.GetSelectedBinId();
    }
    var bins_to_update = new Set();
    let transaction = [];

    if (!Array.isArray(targets)) {
        targets = [targets];
    }

    for (const target of targets) {
        if (target.collapsed)
            continue;

        for (const node of target.IterateChildren()) {
            for (let bin_id in this.selections) {
                if (this.selections[bin_id].has(node)) {
                    this.selections[bin_id].delete(node);
                    bins_to_update.add(bin_id);

                    if (this.keepHistory && node.IsLeaf()) {
                        transaction.push({'type': 'RemoveNode',
                                          'bin_id': bin_id,
                                          'node': node});
                    }
                }
            }

            node.ResetColor();
        }
    }

    bins_to_update = Array.from(bins_to_update);
    this.PushHistory(transaction);
    this.RedrawBins();
    this.UpdateBinsWindow(bins_to_update);
};


Bins.prototype.IsNodeMemberOfBin = function(node) {
    for (let bin_id in this.selections) {
        if (this.selections[bin_id].has(node)) {
            return true;
        }
    }
    return false;
};


Bins.prototype.BinsSorted = function() {
    // move every taxonomy row after their original parent.
    this.container.querySelectorAll('[data-parent]').forEach((elem) => {
        let taxonomy_row = elem.parentNode.removeChild(elem);
        let parent_bin_id = elem.getAttribute('data-parent');

        let parent_row = this.container.querySelector(`tr[bin-id="${parent_bin_id}"]`);
        parent_row.insertAdjacentHTML('afterend', taxonomy_row.outerHTML);
    });
}


Bins.prototype.UpdateBinsWindow = function(bin_list) {
    if (!this.allowRedraw)
        return;

    bin_list = bin_list || Object.keys(this.selections);

    for (let i = 0; i < bin_list.length; i++) {
        let bin_id = bin_list[i];

        if (mode == 'pan') {
            let num_gene_clusters = 0;
            let num_gene_calls = 0;

            for (let node of this.selections[bin_id].values()) {
                if (node.IsLeaf()) {
                    num_gene_clusters++;
                    num_gene_calls += parseInt(item_lengths[node.label]);
                }
            }

            let bin_row = this.container.querySelector(`tr[bin-id="${bin_id}"]`);
            let bin_name = bin_row.querySelector('.bin-name').valie

            bin_row.querySelector('td.num-gene-clusters').setAttribute('data-value', num_gene_clusters);
            bin_row.querySelector('td.num-gene-clusters>input').value = num_gene_clusters;

            if (isNaN(num_gene_calls)) {
                bin_row.querySelector('td.num-gene-calls').setAttribute('data-value', 0);
                bin_row.querySelector('td.num-gene-calls>input').value = 'n/a';
            } else {
                bin_row.querySelector('td.num-gene-calls').setAttribute('data-value', num_gene_calls);
                bin_row.querySelector('td.num-gene-calls>input').value = num_gene_calls;
            }
        } else {
            let num_items = 0;
            let length_sum = 0;

            for (let node of this.selections[bin_id].values()) {
                if (node.IsLeaf()) {
                    num_items++;
                    length_sum += parseInt(item_lengths[node.label]);
                }
            }

            let bin_row = this.container.querySelector(`tr[bin-id="${bin_id}"]`);
            let bin_name = bin_row.querySelector('.bin-name').value;

            bin_row.querySelector('td.num-items').setAttribute('data-value', num_items);
            bin_row.querySelector('td.num-items>input').value = num_items;

            if (isNaN(length_sum)) {
                bin_row.querySelector('td.length-sum').setAttribute('data-value', 0);
                bin_row.querySelector('td.length-sum>span').innerHTML = 'n/a';
            } else {
                bin_row.querySelector('td.length-sum').setAttribute('data-value', length_sum);
                bin_row.querySelector('td.length-sum>span').innerHTML = readableNumber(length_sum);
            }

            if (mode == 'full' || mode == 'refine') {
                $.ajax({
                    type: "POST",
                    url: "/data/completeness",
                    cache: false,
                    data: {
                        'split_names': JSON.stringify(this.GetBinNodeLabels(bin_id)),
                        'bin_name': JSON.stringify($('#bin_name_' + bin_id).val())
                    },
                    success: (data) => {
                        this.cache['completeness'][bin_id] = data;
                        let average_completeness = data['averages']['percent_completion'];
                        let average_redundancy = data['averages']['percent_redundancy'];
                        let best_matching_domain = data['averages']['domain'];
                        let stats = data['stats'];
                        let domain_probabilities = data['averages']['domain_probabilities'];

                        bin_row.querySelector('td.completeness>input').style.background = "";
                        bin_row.querySelector('td.redundancy>input').style.background = "";

                        let _dc = null;
                        let _dr = null;
                        let _vc = null;
                        let _vr = null;
                        let _cc = null;
                        let _cr = null;
                        if (num_items == 0) {
                            _dc = -1; _dr = -1; _vc = '--'; _vr = '--';
                        } else if (best_matching_domain == '') {
                            _dc = -1; _dr = -1; _vc = '***'; _vr = '***'; _cc = '#DDAAAA'; _cr = '#DDAAAA';
                        } else if (best_matching_domain == "blank") {
                            _dc = -1; _dr = -1; _vc = '??'; _vr = '??'; _cc = '#DDDDDD'; _cr = '#DDDDDD';
                        } else if (best_matching_domain == "mixed") {
                            _dc = -1; _dr = -1; _vc = 'xx'; _vr = 'xx'; _cc = '#FF9999'; _cr = '#FF9999';
                        } else if (average_completeness != null && average_redundancy != null) {
                            _dc = average_completeness; _dr = average_redundancy;
                            _vc = average_completeness.toFixed(1);
                            _vr = average_redundancy.toFixed(1);

                            if (average_completeness > 80) {
                                _cc = '#AAFFAA';
                            } else if (average_completeness > 50) {
                                _cc = '#DDFFDD';
                            }

                            if (average_redundancy > 50) {
                                _cr = '#FF7777';
                            } else if (average_redundancy > 7) {
                                _cr = '#FFAAAA';
                            }
                        }

                        bin_row.querySelector('td.completeness').setAttribute('data-value', _dc);
                        bin_row.querySelector('td.completeness>input').style.background = _cc;
                        bin_row.querySelector('td.completeness>input').value = _vc;

                        bin_row.querySelector('td.redundancy').setAttribute('data-value', _dr);
                        bin_row.querySelector('td.redundancy>input').style.background = _cr;
                        bin_row.querySelector('td.redundancy>input').value = _vr;

                        showCompleteness(bin_id, true);
                        showRedundants(bin_id, true);
                    },
                });

                if (!$('#estimate_taxonomy').is(':checked')) {
                    continue;
                }

                let collection_data = {};
                collection_data[bin_name] = [];

                for (let node of this.selections[bin_id].values()) {
                    if (node.IsLeaf()) {
                        collection_data[bin_name].push(node.label);
                    }
                }

                $.ajax({
                    type: "POST",
                    url: "/data/get_taxonomy",
                    cache: false,
                    data: {
                        'collection': JSON.stringify(collection_data)
                    },
                    success: (data) => {
                        if (data.hasOwnProperty('status') && data.status != 0) {
                            if ($('#estimate_taxonomy').is(':checked')) {
                                toastr.error('"' + data.message + '", the server said.', "The anvi'o headquarters is upset");
                            }
                            $('#estimate_taxonomy').prop('checked', false);
                            toggleTaxonomyEstimation();
                            return;
                        }
                        if (data.hasOwnProperty(bin_name)) {
                            let order = ["t_domain", "t_phylum", "t_class",
                                         "t_order", "t_family", "t_genus", "t_species"];

                            this.container.querySelector(`span.taxonomy-name-label[bin-id="${bin_id}"]`).innerHTML = " (?) Unknown";

                            for (let i=order.length-1; i >= 0; i--) {
                                let level = order[i];

                                if (data[bin_name]['consensus_taxonomy'][level] !== null) {
                                    let label = level.split('_')[1][0];
                                    this.container.querySelector(`span.taxonomy-name-label[bin-id="${bin_id}"]`).innerHTML = " (" + label + ") " + data[bin_name]['consensus_taxonomy'][level].replace('_', ' ');
                                    return;
                                }
                            }
                        }
                    }
                });
            }
        }
    }

    $('#bin_settings_tab:not(.active) a').css('color', "#ff0000");
};


Bins.prototype.GetBinNodeLabels = function(bin_id) {
    let node_labels = [];

    for (const node of this.selections[bin_id].values()) {
        if (node.IsLeaf()) {
            node_labels.push(node.label);
        }
    }

    return node_labels;
};


Bins.prototype.MigrateCollection = function() {
    // When we redraw tree, tree gets parsed again and object ids change.
    // we need to migrate exist collection and bind them to newly created objects.

    let collection = this.ExportCollection(use_bin_id=true);
    if (collection){
        for (let bin_id of Object.keys(collection?.['data']).sort() || undefined) {
            if (typeof bin_id === undefined){
                return
            }
            this.selections[bin_id] = new Set(collection['data'][bin_id].map((node_label) => drawer.tree.GetLeafByName(node_label)));
        }
    }
}


Bins.prototype.ImportCollection = function(collection, threshold = 1000) {
    this.keepHistory = false;
    let bins_cleared = false;

    for (let bin_name of Object.keys(collection['data']).sort())
    {
        let nodes = [];
        let sum_length = 0;

        for (let i = 0; i < collection['data'][bin_name].length; i++)
        {
            if (mode != 'full')
            {
                nodes.push(collection['data'][bin_name][i]);
            }
            else if (typeof item_lengths[collection['data'][bin_name][i]] !== 'undefined')
            {
                nodes.push(collection['data'][bin_name][i]);
                sum_length += item_lengths[collection['data'][bin_name][i]];
            }
        }

        if ((mode != 'full' || threshold == 0 || sum_length >= threshold) && nodes.length > 0)
        {
            if (!bins_cleared)
            {
                this.bin_counter = 0;
                this.selections = {};
                this.container.innerHTML = '';
                this.history = [];
                bins_cleared = true;
            }

            this.selections[this.bin_counter] = new Set(nodes.map((node_label) => drawer.tree.GetLeafByName(node_label)));

            let bin_color = (collection['colors'][bin_name]) ? collection['colors'][bin_name] : randomColor();
            this.NewBin(this.bin_counter, {'name': bin_name, 'color': bin_color});
            this.bin_counter++;
        }
    }

    if (bins_cleared) {
        // means we have added new things.
        this.RebuildIntersections();
        this.UpdateBinsWindow();
        this.RedrawBins();
    }

    this.keepHistory = true;
};


Bins.prototype.ExportCollection = function(use_bin_id=false) {
    let data = {};
    let colors = {};

    for (let tr of this.container.querySelectorAll('tr.bin-row')) {
        let bin_id = tr.getAttribute('bin-id');
        let bin_name = tr.querySelector('.bin-name').value;
        let bin_color = tr.querySelector('.colorpicker').getAttribute('color');
        let items = [];

        for (let node of this.selections[bin_id].values()) {
            if (node === null) return;
            if (node.IsLeaf()) {
                items.push(node.label);
            }
        }

        if (items.length > 0) {
            if (use_bin_id) {
                data[bin_id] = items;
                colors[bin_id] = bin_color;
            }
            else
            {
                data[bin_name] = items;
                colors[bin_name] = bin_color;
            }
        }
    }

    return {'data': data, 'colors': colors};
};


Bins.prototype.ExportBin = function(bin_id_to_export, use_bin_id=false) {
    for (let tr of this.container.querySelectorAll('tr.bin-row')) {
        let bin_id = tr.getAttribute('bin-id');

        if (bin_id == bin_id_to_export) {
            let bin_name = tr.querySelector('.bin-name').value;
            let bin_color = tr.querySelector('.colorpicker').getAttribute('color');
            let items = [];

            for (let node of this.selections[bin_id].values()) {
                if (node.IsLeaf()) {
                    items.push(node.label);
                }
            }

            return {'items': items, 'color': bin_color, 'bin_name': bin_name};
        }
    }
};


Bins.prototype.HighlightItems = function(item_list) {
    if (!Array.isArray(item_list)) {
        item_list = [item_list];
    }

    this.higlighted_items = [];
    for (let name of item_list) {
        let node = drawer.tree.GetLeafByName(name);
        if (node) {
            this.higlighted_items.push(name);
        } else {
            console.error('Received highlight request for non-existed node: ' + name);
        }
    }

    this.RedrawBins();
};


Bins.prototype.ClearHighlightedItems = function() {
    this.higlighted_items = [];
    this.RedrawBins();
};


Bins.prototype.RedrawLineColors = function() {
    for (let bin_id in this.selections) {
        if (this.selections[bin_id].size) {
            let bin_color = this.GetBinColor(bin_id);
            for (let node of this.selections[bin_id].values()) {
                if (node === null) return;
                node.SetColor(bin_color);
            }
        }
    }
};

Bins.prototype.DrawInvertedNodes = function(leaf_list, rect_width){

    var inverse_fill_opacity = $('#inverse_fill_opacity').val();
    var inverse_color = document.getElementById('inverse_color').getAttribute('color');
    let nodes_for_inversion = []
    let calculatedRectX = 0

    for(let i = 0; i < leaf_list.length; i++){ //selecting all the nodes that aren't assigned to bins
        if(leaf_list[i] === -1){
            nodes_for_inversion.push(drawer.tree.leaves.filter(leaf => leaf.order === i))
        }
    }

    nodes_for_inversion.map((node, idx) => {
        // map through each node, use the highest 'x' value for all nodes in subsequent rect phylogram render
        // addresses bug of rect shades starting from their relative position in tree
        let p = node[0];
        if(p){
            p.xy.x > calculatedRectX ? calculatedRectX = p.xy.x : null
        }
    })

    nodes_for_inversion.map((node, idx) => {

        let p1;
        let p2;
        let p = node[0];

        if(idx == [nodes_for_inversion.length - 1]){
            return // last node does not have 2 adjacent border nodes, throws error on GetBorderNodes() call
        }
        else {
            [p1, p2] = p.GetBorderNodes();
        }

        if (tree_type == 'circlephylogram'){
            let pie = drawPie(
                'bin',
                'bin_outer_' + idx,
                p1.angle - p1.size / 2,
                p2.angle + p2.size / 2,
                distance(p.backarc, {
                    'x': 0,
                    'y': 0
                }),
                total_radius,
                (p2.angle - p1.angle + (p1.size / 2) + (p2.size / 2) > Math.PI) ? 1 : 0,
                inverse_color,
                inverse_fill_opacity,
                false
            );
            pie.setAttribute('vector-effect', 'non-scaling-stroke');
            pie.setAttribute('stroke-opacity', inverse_fill_opacity);
        } else {
           let rect =  drawPhylogramRectangle('bin',
                'bin_background_' + idx,
                calculatedRectX, // new variable calculated above
                p.xy.y,
                p.size,
                rect_width,
                inverse_color,
                inverse_fill_opacity,
                true
            );
                rect.setAttribute('vector-effect', 'non-scaling-stroke');
                rect.setAttribute('stroke-opacity', '1');
        }
    })
}

Bins.prototype.RedrawBins = function() {
    if (!drawer)
        return;

    if (!this.allowRedraw)
        return;

    var leaf_list = [];
    for (var i=0; i < drawer.tree.leaves.length + 1; i++) {
        leaf_list.push(-1);
    }

    for (let bin_id in this.selections) {
        for (let node of this.selections[bin_id].values()) {
            if (node === null) return;
            if (typeof node === 'undefined')
            {
                this.selections[bin_id].delete(node);
                continue;
            }

            if (node.IsLeaf() && !node.collapsed) {
                leaf_list[node.order] = bin_id;
            }
        }
    }

    // cluster bins and put them into bins_to_draw array with (start, end, bin_id);
    var prev_value = leaf_list[0];
    var prev_start = 0;

    var bins_to_draw = new Array();
    for (var i=1; i < leaf_list.length; i++)
    {
        if (prev_value != leaf_list[i])
        {
            if (prev_value != -1) {
                bins_to_draw.push(new Array(prev_start, i - 1, prev_value)); // start, end, bin_id;
            }
            prev_start = i;
        }
        prev_value = leaf_list[i];
    }

    var bin = document.getElementById('bin');
    while (bin.hasChildNodes()) {
        bin.removeChild(bin.lastChild);
    }

    // draw new bins
    var show_grid = $('#show_grid_for_bins')[0].checked;
    var show_shade = $('#show_shade_for_bins')[0].checked;
    var invert_shade = $('#invert_shade_for_bins')[0].checked;
    var shade_fill_opacity = $('#shade_fill_opacity').val();
    var grid_color = document.getElementById('grid_color').getAttribute('color');
    var grid_width = $('#grid_width').val();
    var show_bin_labels = $('#show_bin_labels')[0].checked;
    var bin_labels_font_size = (parseFloat($('#bin_labels_font_size').val())) ? parseFloat($('#bin_labels_font_size').val()) : 32;
    var autorotate_bin_labels = $('#autorotate_bin_labels')[0].checked;
    var bin_labels_angle = $('#bin_labels_angle').val();
    var background_starts_from_branch = $('#begins_from_branch').is(':checked');

    var outer_ring_size = parseFloat($('#outer-ring-height').val());
    var outer_ring_margin = parseFloat($('#outer-ring-margin').val());

    var width_with_grids;
    var width_no_grids;

    for (var i=0; i < bins_to_draw.length; i++) {
        var start = drawer.tree.leaves[bins_to_draw[i][0]];
        var end = drawer.tree.leaves[bins_to_draw[i][1]];

        if (background_starts_from_branch) {
            var startAncestors = new Set(start.GetAncestors());
            var endAncestors = new Set(end.GetAncestors());
            var intersection = new Set([...startAncestors].filter(x => endAncestors.has(x))).values().next().value;

            if (typeof intersection === 'undefined') {
                throw `It seems Node:${start.id} and Node:${end.id} does not have common ancestor.`;
            }
        }

        var color = document.getElementById('bin_color_' + bins_to_draw[i][2]).getAttribute('color');

        if (tree_type == 'circlephylogram')
        {
            drawPie('bin',
                'bin_outer_' + i,
                start.angle - start.size / 2,
                end.angle + end.size / 2,
                total_radius + outer_ring_margin,
                total_radius + outer_ring_margin + outer_ring_size,
                (end.angle - start.angle + (start.size / 2) + (end.size / 2) > Math.PI) ? 1 : 0,
                color,
                1,
                true);

            var align = 'left';
            var angle = (end.angle + end.size / 2 + start.angle - start.size / 2) / 2;
            var new_angle = angle * 180.0 / Math.PI;
            if ((angle > Math.PI / 2.0) && (angle < 1.5 * Math.PI)) {
                align = 'right';
                new_angle += 180.0;
            }

            if (show_bin_labels)
            {
                var bin_label_radius = total_radius + outer_ring_margin * 1.5 + outer_ring_size * (this.higlighted_items.length > 0 ? 2 : 1);
                var bin_label_angle = (end.angle + end.size / 2 + start.angle - start.size / 2) / 2;

                var bin_label_px = bin_label_radius * Math.cos(bin_label_angle);
                var bin_label_py = bin_label_radius * Math.sin(bin_label_angle);

                bin_label_px = bin_label_px - Math.cos(Math.PI / 2 + bin_label_angle) * (bin_labels_font_size / 3) * (align == 'right' ? 1 : -1);
                bin_label_py = bin_label_py - Math.sin(Math.PI / 2 + bin_label_angle) * (bin_labels_font_size / 3) * (align == 'right' ? 1 : -1);

                drawRotatedText(
                    'bin',
                    {
                        'x': bin_label_px,
                        'y': bin_label_py,
                    },
                    $('#bin_name_' + bins_to_draw[i][2]).val().replaceAll("_", " "),
                    (autorotate_bin_labels) ? new_angle : bin_labels_angle,
                    align,
                    bin_labels_font_size + "px",
                    "HelveticaNeue-CondensedBold, Helvetica Neue, Helvetica, sans-serif",
                    color,
                    0,
                    'baseline'
                    );

            }

            var pie = drawPie('bin',
                'bin_background_' + i,
                start.angle - start.size / 2,
                end.angle + end.size / 2,
                (background_starts_from_branch) ? intersection.radius : beginning_of_layers,
                (show_grid) ? total_radius + outer_ring_margin + outer_ring_size : total_radius,
                (Math.abs(end.angle - start.angle) + start.size / 2 + end.size / 2 > Math.PI) ? 1 : 0,
                color,
                (show_grid) ? 0 : shade_fill_opacity,
                false);

            if (show_grid && !show_shade) {
                pie.setAttribute('vector-effect', 'non-scaling-stroke');
                pie.setAttribute('stroke-opacity', '1');
                pie.setAttribute('stroke-width', grid_width);
                pie.setAttribute('stroke', grid_color);
            } else if(show_grid && show_shade){
                pie.setAttribute('vector-effect', 'non-scaling-stroke');
                pie.setAttribute('stroke-opacity', '1');
                pie.setAttribute('stroke-width', grid_width);
                pie.setAttribute('stroke', grid_color);
                pie.setAttribute('fill-opacity', shade_fill_opacity)
            } else if(!show_grid && !show_shade){
                pie.setAttribute('vector-effect', 'non-scaling-stroke');
                pie.setAttribute('stroke-opacity', '0');
                pie.setAttribute('stroke-width', 0);
                pie.setAttribute('stroke', grid_color);
                pie.setAttribute('fill-opacity', 0)
            }
        }
        else
        {
            var height = end.xy['y'] + end.size / 2 - start.xy['y'] + start.size / 2;

            drawPhylogramRectangle('bin',
                'bin_outer_' + i,
                total_radius + outer_ring_margin,
                start.xy['y'] - start.size / 2 + height / 2,
                height,
                outer_ring_size,
                color,
                1,
                true);

            if (show_bin_labels)
            {
                drawRotatedText(
                    'bin',
                    {
                        'y':  (start.xy.y - start.size / 2 + end.xy.y + end.size / 2) / 2 + (bin_labels_font_size / 3),
                        'x': (total_radius + outer_ring_margin * 1.5 + outer_ring_size * (this.higlighted_items.length > 0 ? 2 : 1)),
                    },
                    $('#bin_name_' + bins_to_draw[i][2]).val().replace("_", " "),
                    (autorotate_bin_labels) ? 0 : bin_labels_angle,
                    'left',
                    bin_labels_font_size + "px",
                    "HelveticaNeue-CondensedBold, Helvetica Neue, Helvetica, sans-serif",
                    color,
                    0,
                    'baseline'
                    );

            }

            let backgroundStart = (background_starts_from_branch) ? intersection.xy['x'] : beginning_of_layers;

            width_with_grids = total_radius + outer_ring_margin + outer_ring_size - backgroundStart
            width_no_grids = total_radius - backgroundStart
            // ^ these get passed to invert_shade method so we don't have to re-declare a bunch of other stuff.

            var rect = drawPhylogramRectangle('bin',
                'bin_background_' + i,
                backgroundStart,
                start.xy['y'] - start.size / 2 + height / 2,
                height,
                (show_grid) ? width_with_grids : width_no_grids,
                color,
                (show_grid) ? 0 : shade_fill_opacity,
                false);

            if (show_grid && !show_shade) {
                rect.setAttribute('vector-effect', 'non-scaling-stroke');
                rect.setAttribute('stroke-opacity', '1');
                rect.setAttribute('stroke-width', grid_width);
                rect.setAttribute('stroke', grid_color);
            } else if(show_grid && show_shade){
                rect.setAttribute('vector-effect', 'non-scaling-stroke');
                rect.setAttribute('stroke-opacity', '1');
                rect.setAttribute('stroke-width', grid_width);
                rect.setAttribute('stroke', grid_color);
                rect.setAttribute('fill-opacity', shade_fill_opacity)
            } else if(!show_grid && !show_shade){
                rect.setAttribute('vector-effect', 'non-scaling-stroke');
                rect.setAttribute('stroke-opacity', '1');
                rect.setAttribute('stroke-width', 0);
                rect.setAttribute('stroke', grid_color);
                rect.setAttribute('fill-opacity', 0)
            }
        }
    }

    for (const name of this.higlighted_items) {
        let node = drawer.tree.GetLeafByName(name);
        let color = document.getElementById('picker_highlight').getAttribute('color');

        if (tree_type == 'circlephylogram')
        {
            drawPie('bin',
                'bin_outer_' + 1,
                node.angle - node.size / 2,
                node.angle + node.size / 2,
                total_radius + outer_ring_margin + outer_ring_size,
                total_radius + outer_ring_margin + outer_ring_size * 2,
                (node.size / 2 > Math.PI) ? 1 : 0,
                color,
                1,
                true);
        }
        else
        {
            drawPhylogramRectangle('bin',
                'bin_outer_' + 1,
                total_radius + outer_ring_margin + outer_ring_size,
                node.xy['y'],
                node.size,
                outer_ring_size,
                color,
                1,
                true);
        }
    }
    invert_shade ? this.DrawInvertedNodes(leaf_list, show_grid ? width_with_grids : width_no_grids) : null

}

Bins.prototype.RebuildIntersections = function() {
    if (!this.allowRedraw)
        return;

    for (let bin_id in this.selections) {
        if (!this.selections[bin_id].size) {
            continue;
        }
        let bin_color = this.GetBinColor(bin_id);
        let inserted = true;
        let removed = true;

        while (removed) {
            removed = false;

            for (let node of this.selections[bin_id].values()) {
                if (node === null) return;
                if (node.IsLeaf()) {
                    continue;
                }

                let any_child_in_bin = false;
                let child = node.child;

                while (child.sibling) {
                    any_child_in_bin = any_child_in_bin || this.selections[bin_id].has(child);
                    child = child.sibling;
                }

                if (!any_child_in_bin) {
                    // node doesn't have any child in this bin so let's get rid of it.
                    this.selections[bin_id].delete(node);
                    node.ResetColor();
                }
            }
        }

        while (inserted) {
            inserted = false;
            for (let node of this.selections[bin_id].values()) {
                let parent = node.ancestor;

                if (!parent) {
                    // no parent to add
                    continue;
                }

                if (this.selections[bin_id].has(parent)) {
                    // parent already in bin
                    continue;
                }

                if (node.sibling && this.selections[bin_id].has(node.sibling)) {
                    // node and its sibling in same bin, so parent should too.
                    this.selections[bin_id].add(parent);
                    parent.SetColor(bin_color);
                    inserted = true;
                }
            }
        }
    }
}
