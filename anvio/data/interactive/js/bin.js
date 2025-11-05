/**
 * Draw bins, bin labels stuff.
 *
 *  Authors: Özcan Esen <ozcanesen@gmail.com>
 *           A. Murat Eren <a.murat.eren@gmail.com>
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

// ============================================================================
// Constants
// ============================================================================

const MAX_HISTORY_SIZE = 50;

const BIN_DEFAULTS = {
    PREFIX: 'Bin_',
    DEFAULT_VALUES: {
        contig_count: 0,
        contig_length: 'N/A',
        num_gene_clusters: '---',
        num_gene_calls: '---',
        completeness: '---',
        redundancy: '---'
    },
    COLORS: {
        COMPLETENESS: {
            HIGH: '#AAFFAA',      // > 80%
            MEDIUM: '#DDFFDD',    // > 50%
            ERROR: '#DDAAAA',     // Parse error
            BLANK: '#DDDDDD',     // Blank domain
            MIXED: '#FF9999'      // Mixed domain
        },
        REDUNDANCY: {
            HIGH: '#FF7777',      // > 50%
            MEDIUM: '#FFAAAA'     // > 7%
        }
    }
};

const TRANSACTION_TYPES = {
    APPEND_NODE: 'AppendNode',
    REMOVE_NODE: 'RemoveNode',
    CHANGE_COLOR: 'ChangeColor',
    DELETE_BIN: 'DeleteBin',
    NEW_BIN: 'NewBin'
};

// ============================================================================
// Bins Class Constructor
// ============================================================================

/**
 * Bins manager for handling bin selections and operations
 * @param {string} prefix - Bin name prefix
 * @param {HTMLElement} container - Container element for bins
 */
function Bins(prefix, container) {
    this.selections = {};
    this.bin_counter = 0;
    this.prefix = prefix || BIN_DEFAULTS.PREFIX;
    this.higlighted_items = [];
    this.container = container || document.createElement("div");

    this.cache = {
        'completeness': {},
        'taxonomy': {}
    };

    this.keepHistory = false;
    this.allowRedraw = true;

    this.history = [];
    this.future = [];

    // Event listeners
    document.body.addEventListener('bin-settings-changed', (event) => this.RedrawBins());
    $(document).on('sorted', () => { this.BinsSorted(); });
}

// ============================================================================
// Bin Creation and Management
// ============================================================================

/**
 * Create a new bin
 * @param {number} id - Bin ID (optional)
 * @param {Object} binState - Bin state object with name and color
 */
Bins.prototype.NewBin = function(id, binState) {
    const binData = this._prepareBinData(id, binState);
    const template = this._generateBinTemplate(binData);


    this.container.insertAdjacentHTML('beforeend', template);
    this.SelectLastRadio();

    this._initializeColorPicker(binData.id);
    this._recordBinCreation(binData);
};

/**
 * Prepare bin data for creation
 * @private
 */
Bins.prototype._prepareBinData = function(id, binState) {
    if (typeof id === 'undefined') {
        // Creating new bin
        const newId = this.bin_counter;
        this.selections[newId] = new Set();
        this.bin_counter++;

        return {
            id: newId,
            name: this.prefix + (newId + 1),
            color: randomColor({luminosity: 'dark'}),
            from_state: false,
            ...BIN_DEFAULTS.DEFAULT_VALUES
        };
    } else {
        // Loading bin from collection
        return {
            id: id,
            name: binState['name'],
            color: binState['color'],
            from_state: true,
            ...BIN_DEFAULTS.DEFAULT_VALUES
        };
    }
};

/**
 * Generate HTML template for bin
 * @private
 */
Bins.prototype._generateBinTemplate = function(binData) {
    const {id, name, color} = binData;

    let template = `
        <tr bin-id="${id}" class="bin-row">
            <td><input type="radio" name="active_bin" value="${id}"></td>
            <td><div id="bin_color_${id}" class="colorpicker" color="${color}"
                     style="background-color: ${color}"></div></td>
            <td data-value="${name}">
                <input type="text" class="bin-name"
                       oninput="this.value = event.target.value.replaceAll(' ', '_');"
                       onchange="emit('bin-settings-changed'); this.parentNode.setAttribute('data-value', this.value);"
                       size="21" id="bin_name_${id}" value="${name}">
            </td>`;

    // Mode-specific columns
    if (mode === 'pan' || mode === 'structure') {
        template += this._getPanModeColumns(id, binData);
    } else if (mode === 'codon-frequencies') {
        template += this._getCodonFrequencyColumns(id, binData);
    } else {
        template += this._getDefaultModeColumns(id, binData);
    }

    // Delete button
    template += `
            <td><center>
                <span class="default-bin-icon bi bi-trash-fill fa-lg"
                      aria-hidden="true" alt="Delete this bin" title="Delete this bin"
                      onclick="bins.DeleteBin(${id});">
                </span>
            </center></td>
        </tr>`;

    // Taxonomy row (if applicable)
    if (this._shouldShowTaxonomyRow()) {
        template += this._getTaxonomyRowTemplate(id);
    }

    return template;
};

/**
 * Get pan mode specific columns
 * @private
 */
Bins.prototype._getPanModeColumns = function(id, binData) {
    return `
        <td data-value="${binData.num_gene_clusters}" class="num-gene-clusters">
            <input type="button" value="${binData.num_gene_clusters}"
                   title="Click for quick gene cluster summaries"
                   onclick="showGeneClusterDetails(${id});">
        </td>
        <td data-value="${binData.num_gene_calls}" class="num-gene-calls">
            <input type="button" value="${binData.num_gene_calls}">
        </td>`;
};

/**
 * Get codon frequency mode columns
 * @private
 */
Bins.prototype._getCodonFrequencyColumns = function(id, binData) {
    return `
        <td data-value="${binData.contig_count}" class="num-items">
            <input type="button" value="${binData.contig_count}"
                   title="Click for contig names"
                   onclick="showGeneFunctions(${id});">
        </td>`;
};

/**
 * Get default mode columns
 * @private
 */
Bins.prototype._getDefaultModeColumns = function(id, binData) {
    return `
        <td data-value="${binData.contig_count}" class="num-items">
            <input type="button" value="${binData.contig_count}"
                   title="Click for contig names"
                   onclick="showContigNames(${id});">
        </td>
        <td data-value="${binData.contig_length}" class="length-sum">
            <span>${binData.contig_length}</span>
        </td>
        <td data-value="${binData.completeness}" class="completeness">
            <input type="button" value="${binData.completeness}"
                   title="Click for completeness table"
                   onclick="showCompleteness(${id});">
        </td>
        <td data-value="${binData.redundancy}" class="redundancy">
            <input type="button" value="${binData.redundancy}"
                   title="Click for redundant hits"
                   onclick="showRedundants(${id});">
        </td>`;
};

/**
 * Check if taxonomy row should be shown
 * @private
 */
Bins.prototype._shouldShowTaxonomyRow = function() {
    const validModes = ['full', 'refine', 'manual', 'pan'];
    return validModes.includes(mode);
};

/**
 * Get taxonomy row template
 * @private
 */
Bins.prototype._getTaxonomyRowTemplate = function(id) {
    const displayStyle = $('#estimate_taxonomy').is(':checked') ? '' : 'display: none;';
    return `
        <tr style="${displayStyle}" data-parent="${id}">
            <td style="border-top: 0px;">&nbsp;</td>
            <td style="border-top: 0px;">&nbsp;</td>
            <td colspan="6" style="border-top: 0px; padding-top: 0px;">
                <span bin-id="${id}" class="taxonomy-name-label">N/A</span>
            </td>
        </tr>`;
};

/**
 * Initialize color picker for bin
 * @private
 */
Bins.prototype._initializeColorPicker = function(id) {
    $(`#bin_color_${id}`).colpick({
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
        onHide: (picker) => {
            const el = $(picker).data('colpick').el;
            const current_color = $(el).attr('color');
            const previous_color = $(el).attr('color-before');

            if (current_color != previous_color) {
                emit('bin-settings-changed');
                this.PushHistory([{
                    'type': TRANSACTION_TYPES.CHANGE_COLOR,
                    'bin_id': id,
                    'color-before': previous_color,
                    'color': current_color
                }]);
            }
        }
    }).keyup(function() {
        $(this).colpickSetColor(this.value);
    });
};

/**
 * Record bin creation in history
 * @private
 */
Bins.prototype._recordBinCreation = function(binData) {
    this.PushHistory([{
        'type': TRANSACTION_TYPES.NEW_BIN,
        'name': binData.name,
        'bin_id': binData.id,
        'color': binData.color
    }]);
};

/**
 * Delete a bin
 * @param {number} bin_id - Bin ID to delete
 * @param {boolean} show_confirm - Show confirmation dialog
 */
Bins.prototype.DeleteBin = function(bin_id, show_confirm = true) {
    if (show_confirm && !confirm('Are you sure?')) {
        return;
    }

    const transaction = this._prepareDeletionTransaction(bin_id);
    this._clearBinSelections(bin_id);
    this.PushHistory(transaction);
    this._removeBinFromDOM(bin_id);
    this._ensureAtLeastOneBin();
    this.RedrawBins();
};

/**
 * Prepare deletion transaction for history
 * @private
 */
Bins.prototype._prepareDeletionTransaction = function(bin_id) {
    const transaction = [];

    // Record node removals
    for (const node of this.selections[bin_id].values()) {
        node.ResetColor();
        if (this.keepHistory) {
            transaction.push({
                'type': TRANSACTION_TYPES.REMOVE_NODE,
                'bin_id': bin_id,
                'node': node
            });
        }
    }

    // Record bin deletion
    const binElement = document.getElementById(`bin_name_${bin_id}`);
    const colorElement = document.getElementById(`bin_color_${bin_id}`);

    transaction.push({
        'type': TRANSACTION_TYPES.DELETE_BIN,
        'bin_id': bin_id,
        'name': binElement.value,
        'color': colorElement.getAttribute('color')
    });

    return transaction;
};

/**
 * Clear bin selections
 * @private
 */
Bins.prototype._clearBinSelections = function(bin_id) {
    this.selections[bin_id].clear();
};

/**
 * Remove bin from DOM
 * @private
 */
Bins.prototype._removeBinFromDOM = function(bin_id) {
    const bin_row = this.container.querySelector(`tr[bin-id='${bin_id}']`);

    // Remove taxonomy row if exists
    if (bin_row.nextElementSibling) {
        bin_row.nextElementSibling.remove();
    } else {
        console.error('Next sibling element not found for bin_row.');
    }

    bin_row.remove();
};

/**
 * Ensure at least one bin exists
 * @private
 */
Bins.prototype._ensureAtLeastOneBin = function() {
    if (!this.container.querySelectorAll('*').length) {
        this.NewBin();
    }

    if (!this.container.querySelector('input[name=active_bin]:checked')) {
        this.SelectLastRadio();
    }
};

/**
 * Delete all bins
 */
Bins.prototype.DeleteAllBins = function() {
    if (!confirm('Are you sure you want to remove all bins?')) {
        return;
    }

    const binRows = this.container.querySelectorAll('tr.bin-row');
    for (const tr of binRows) {
        this.DeleteBin(tr.getAttribute('bin-id'), false);
    }
};

// ============================================================================
// Selection Management
// ============================================================================

/**
 * Select the last radio button
 */
Bins.prototype.SelectLastRadio = function() {
    const radios = this.container.querySelectorAll('input[name=active_bin]');
    if (radios.length > 0) {
        radios[radios.length - 1].checked = true;
    }
};

/**
 * Get currently selected bin ID
 * @returns {string} Selected bin ID
 */
Bins.prototype.GetSelectedBinId = function() {
    const selected = this.container.querySelector('input[name=active_bin]:checked');
    return selected ? selected.value : null;
};

/**
 * Get total selection count across all bins
 * @returns {number} Total selection count
 */
Bins.prototype.GetTotalSelectionCount = function() {
    let sum = 0;
    for (const bin_id in this.selections) {
        sum += this.selections[bin_id].size;
    }
    return sum;
};

/**
 * Get selected bin color
 * @returns {string} Hex color
 */
Bins.prototype.GetSelectedBinColor = function() {
    return this.GetBinColor(this.GetSelectedBinId());
};

/**
 * Get bin color by ID
 * @param {number} bin_id - Bin ID
 * @returns {string} Hex color
 */
Bins.prototype.GetBinColor = function(bin_id) {
    const colorElement = this.container.querySelector(`#bin_color_${bin_id}`);
    return colorElement ? colorElement.getAttribute('color') : null;
};

// ============================================================================
// Node Operations
// ============================================================================

/**
 * Append nodes to a bin
 * @param {Array|Object} targets - Node(s) to append
 * @param {number} bin_id - Target bin ID (optional, uses selected)
 */
Bins.prototype.AppendNode = function(targets, bin_id) {
    if (typeof bin_id === 'undefined') {
        bin_id = this.GetSelectedBinId();
    }

    const bin_color = this.GetBinColor(bin_id);
    const bins_to_update = new Set();
    const transaction = [];

    // Ensure targets is an array
    if (!Array.isArray(targets)) {
        targets = [targets];
    }

    // Process each target
    for (const target of targets) {
        if (target.collapsed) continue;

        for (const node of target.IterateChildren()) {
            this._removeNodeFromOtherBins(node, bin_id, bins_to_update, transaction);
            this._addNodeToBin(node, bin_id, bin_color, bins_to_update, transaction);
        }
    }

    this._finalizeNodeOperation(bins_to_update, transaction);
};

/**
 * Remove node from other bins
 * @private
 */
Bins.prototype._removeNodeFromOtherBins = function(node, current_bin_id, bins_to_update, transaction) {
    for (const other_bin_id in this.selections) {
        if (other_bin_id == current_bin_id) continue;

        if (this.selections[other_bin_id].has(node)) {
            this.selections[other_bin_id].delete(node);
            bins_to_update.add(other_bin_id);

            if (this.keepHistory && node.IsLeaf()) {
                transaction.push({
                    'type': TRANSACTION_TYPES.REMOVE_NODE,
                    'bin_id': other_bin_id,
                    'node': node
                });
            }
        }
    }
};

/**
 * Add node to bin
 * @private
 */
Bins.prototype._addNodeToBin = function(node, bin_id, bin_color, bins_to_update, transaction) {
    if (!this.selections[bin_id].has(node)) {
        this.selections[bin_id].add(node);
        bins_to_update.add(bin_id);

        if (this.keepHistory && node.IsLeaf()) {
            transaction.push({
                'type': TRANSACTION_TYPES.APPEND_NODE,
                'bin_id': bin_id,
                'node': node
            });
        }
    }

    node.SetColor(bin_color);
};

/**
 * Finalize node operation
 * @private
 */
Bins.prototype._finalizeNodeOperation = function(bins_to_update, transaction) {
    const bins_array = Array.from(bins_to_update);
    this.PushHistory(transaction);
    this.RedrawBins();
    this.UpdateBinsWindow(bins_array);
};

/**
 * Remove nodes from bins
 * @param {Array|Object} targets - Node(s) to remove
 * @param {number} bin_id - Source bin ID (optional)
 */
Bins.prototype.RemoveNode = function(targets, bin_id) {
    if (typeof bin_id === 'undefined') {
        bin_id = this.GetSelectedBinId();
    }

    const bins_to_update = new Set();
    const transaction = [];

    // Ensure targets is an array
    if (!Array.isArray(targets)) {
        targets = [targets];
    }

    // Process each target
    for (const target of targets) {
        if (target.collapsed) continue;

        for (const node of target.IterateChildren()) {
            this._removeNodeFromAllBins(node, bins_to_update, transaction);
            node.ResetColor();
        }
    }

    this._finalizeNodeOperation(bins_to_update, transaction);
};

/**
 * Remove node from all bins
 * @private
 */
Bins.prototype._removeNodeFromAllBins = function(node, bins_to_update, transaction) {
    for (const bin_id in this.selections) {
        if (this.selections[bin_id].has(node)) {
            this.selections[bin_id].delete(node);
            bins_to_update.add(bin_id);

            if (this.keepHistory && node.IsLeaf()) {
                transaction.push({
                    'type': TRANSACTION_TYPES.REMOVE_NODE,
                    'bin_id': bin_id,
                    'node': node
                });
            }
        }
    }
};

/**
 * Check if node is member of any bin
 * @param {Object} node - Node to check
 * @returns {boolean} True if node is in any bin
 */
Bins.prototype.IsNodeMemberOfBin = function(node) {
    for (const bin_id in this.selections) {
        if (this.selections[bin_id].has(node)) {
            return true;
        }
    }
    return false;
};

// ============================================================================
// History Management
// ============================================================================

/**
 * Push transaction to history
 * @param {Array} transaction - Transaction to record
 */
Bins.prototype.PushHistory = function(transaction) {
    if (!this.keepHistory) return;

    this.history.push(transaction);

    if (this.history.length > MAX_HISTORY_SIZE) {
        this.history.shift();
    }

    // Adding to history clears the future
    this.future = [];
};

/**
 * Undo last operation
 */
Bins.prototype.Undo = function() {
    const transaction = this.history.pop();

    if (transaction) {
        this.ProcessTransaction(transaction, true);
        this.future.push(transaction);
    } else {
        toastr.warning('Can\'t do undo, history is empty.');
    }
};

/**
 * Redo previously undone operation
 */
Bins.prototype.Redo = function() {
    const transaction = this.future.pop();

    if (transaction) {
        this.ProcessTransaction(transaction);
        this.history.push(transaction);
    } else {
        toastr.warning('Can\'t do redo, future is empty.');
    }
};

/**
 * Process a transaction
 * @param {Array} transaction - Transaction to process
 * @param {boolean} reversed - Process in reverse
 */
Bins.prototype.ProcessTransaction = function(transaction, reversed = false) {
    // Disable history and redraw during processing
    this.keepHistory = false;
    this.allowRedraw = false;

    const updated_bins = new Set();
    const removed_bins = new Set();

    // Process each operation
    for (const operation of transaction) {
        updated_bins.add(operation.bin_id);

        if (reversed) {
            this._processReversedOperation(operation, removed_bins);
        } else {
            this._processForwardOperation(operation, removed_bins);
        }
    }

    // Update affected bins
    const bins_to_update = Array.from(updated_bins).filter(
        bin_id => !removed_bins.has(bin_id)
    );

    // Re-enable features and update
    this.keepHistory = true;
    this.allowRedraw = true;
    this.RebuildIntersections();
    this.RedrawBins();
    this.UpdateBinsWindow(bins_to_update);
};

/**
 * Process operation in reverse
 * @private
 */
Bins.prototype._processReversedOperation = function(operation, removed_bins) {
    const operationMap = {
        [TRANSACTION_TYPES.APPEND_NODE]: () => this.RemoveNode(operation.node, operation.bin_id),
        [TRANSACTION_TYPES.REMOVE_NODE]: () => this.AppendNode(operation.node, operation.bin_id),
        [TRANSACTION_TYPES.CHANGE_COLOR]: () => {
            $(`#bin_color_${operation.bin_id}`).attr('color', operation['color-before']);
            $(`#bin_color_${operation.bin_id}`).css('background-color', operation['color-before']);
        },
        [TRANSACTION_TYPES.DELETE_BIN]: () => {
            this.NewBin(operation.bin_id, {
                'name': operation.name,
                'color': operation.color
            });
            removed_bins.delete(operation.bin_id);
        },
        [TRANSACTION_TYPES.NEW_BIN]: () => {
            this.DeleteBin(operation.bin_id, false);
            removed_bins.add(operation.bin_id);
        }
    };

    const handler = operationMap[operation.type];
    if (handler) handler();
};

/**
 * Process operation forward
 * @private
 */
Bins.prototype._processForwardOperation = function(operation, removed_bins) {
    const operationMap = {
        [TRANSACTION_TYPES.APPEND_NODE]: () => this.AppendNode(operation.node, operation.bin_id),
        [TRANSACTION_TYPES.REMOVE_NODE]: () => this.RemoveNode(operation.node, operation.bin_id),
        [TRANSACTION_TYPES.CHANGE_COLOR]: () => {
            $(`#bin_color_${operation.bin_id}`).attr('color', operation.color);
            $(`#bin_color_${operation.bin_id}`).css('background-color', operation.color);
        },
        [TRANSACTION_TYPES.DELETE_BIN]: () => {
            this.DeleteBin(operation.bin_id, false);
            removed_bins.add(operation.bin_id);
        },
        [TRANSACTION_TYPES.NEW_BIN]: () => {
            this.NewBin(operation.bin_id, {
                'name': operation.name,
                'color': operation.color
            });
            removed_bins.delete(operation.bin_id);
        }
    };

    const handler = operationMap[operation.type];
    if (handler) handler();
};

// ============================================================================
// Sorting and DOM Updates
// ============================================================================

/**
 * Handle bins sorted event
 */
Bins.prototype.BinsSorted = function() {
    // Move taxonomy rows after their parent bins
    const taxonomyRows = this.container.querySelectorAll('[data-parent]');

    taxonomyRows.forEach((elem) => {
        const taxonomy_row = elem.parentNode.removeChild(elem);
        const parent_bin_id = elem.getAttribute('data-parent');
        const parent_row = this.container.querySelector(`tr[bin-id="${parent_bin_id}"]`);

        if (parent_row) {
            parent_row.insertAdjacentHTML('afterend', taxonomy_row.outerHTML);
        }
    });
};

// ============================================================================
// Window Updates
// ============================================================================

/**
 * Update bins window with current statistics
 * @param {Array} bin_list - List of bin IDs to update (optional)
 */
Bins.prototype.UpdateBinsWindow = function(bin_list) {
    if (!this.allowRedraw) return;

    bin_list = bin_list || Object.keys(this.selections);

    for (const bin_id of bin_list) {
        if (mode === 'pan' || mode === 'structure') {
            this._updatePanModeStatistics(bin_id);
        } else if (mode === 'codon-frequencies') {
            this._updateCodonFrequencyStatistics(bin_id);
        } else {
            this._updateDefaultModeStatistics(bin_id);
        }
    }

    $('#bin_settings_tab:not(.active) a').css('color', "#ff0000");
};

/**
 * Update pan mode statistics
 * @private
 */
Bins.prototype._updatePanModeStatistics = function(bin_id) {
    let num_gene_clusters = 0;
    let num_gene_calls = 0;

    for (const node of this.selections[bin_id].values()) {
        if (node.IsLeaf()) {
            num_gene_clusters++;
            num_gene_calls += parseInt(item_lengths[node.label]) || 0;
        }
    }

    const bin_row = this.container.querySelector(`tr[bin-id="${bin_id}"]`);

    this._updateTableCell(bin_row, 'td.num-gene-clusters', num_gene_clusters);
    this._updateTableCell(bin_row, 'td.num-gene-calls',
                          isNaN(num_gene_calls) ? 'n/a' : num_gene_calls,
                          isNaN(num_gene_calls) ? 0 : num_gene_calls);
};

/**
 * Update codon frequency statistics
 * @private
 */
Bins.prototype._updateCodonFrequencyStatistics = function(bin_id) {
    let num_items = 0;

    for (const node of this.selections[bin_id].values()) {
        if (node.IsLeaf()) {
            num_items++;
        }
    }

    const bin_row = this.container.querySelector(`tr[bin-id="${bin_id}"]`);
    this._updateTableCell(bin_row, 'td.num-items', num_items);
};

/**
 * Update default mode statistics
 * @private
 */
Bins.prototype._updateDefaultModeStatistics = function(bin_id) {
    let num_items = 0;
    let length_sum = 0;

    for (const node of this.selections[bin_id].values()) {
        if (node.IsLeaf()) {
            num_items++;
            length_sum += parseInt(item_lengths[node.label]) || 0;
        }
    }

    const bin_row = this.container.querySelector(`tr[bin-id="${bin_id}"]`);

    // Update item count
    this._updateTableCell(bin_row, 'td.num-items', num_items);

    // Update length sum
    const lengthCell = bin_row.querySelector('td.length-sum');
    if (lengthCell) {
        lengthCell.setAttribute('data-value', isNaN(length_sum) ? 0 : length_sum);
        lengthCell.querySelector('span').innerHTML =
            isNaN(length_sum) ? 'n/a' : readableNumber(length_sum);
    }

    // Update completeness and taxonomy if needed
    if (mode === 'full' || mode === 'refine') {
        this._updateCompletenessAndTaxonomy(bin_id, bin_row, num_items);
    }
};

/**
 * Update table cell helper
 * @private
 */
Bins.prototype._updateTableCell = function(row, selector, displayValue, dataValue) {
    const cell = row.querySelector(selector);
    if (cell) {
        cell.setAttribute('data-value', dataValue !== undefined ? dataValue : displayValue);
        const input = cell.querySelector('input');
        if (input) {
            input.value = displayValue;
        }
    }
};

/**
 * Update completeness and taxonomy data
 * @private
 */
Bins.prototype._updateCompletenessAndTaxonomy = function(bin_id, bin_row, num_items) {
    const bin_name = $(`#bin_name_${bin_id}`).val();

    // Fetch completeness data
    $.ajax({
        type: "POST",
        url: "/data/completeness",
        cache: false,
        data: {
            'split_names': JSON.stringify(this.GetBinNodeLabels(bin_id)),
            'bin_name': JSON.stringify(bin_name)
        },
        success: (data) => {
            this._processCompletenessData(bin_id, bin_row, data, num_items);
        }
    });

    // Fetch taxonomy if enabled
    if ($('#estimate_taxonomy').is(':checked')) {
        this._fetchTaxonomyData(bin_id, bin_name);
    }
};

/**
 * Process completeness data response
 * @private
 */
Bins.prototype._processCompletenessData = function(bin_id, bin_row, data, num_items) {
    this.cache['completeness'][bin_id] = data;

    const average_completeness = data['averages']['percent_completion'];
    const average_redundancy = data['averages']['percent_redundancy'];
    const best_matching_domain = data['averages']['domain'];

    const completenessData = this._calculateCompletenessMetrics(
        num_items,
        best_matching_domain,
        average_completeness,
        average_redundancy
    );

    this._updateCompletenessDisplay(bin_row, completenessData);

    // Trigger completeness display update
    showCompleteness(bin_id, true);
    showRedundants(bin_id, true);
};

/**
 * Calculate completeness metrics
 * @private
 */
Bins.prototype._calculateCompletenessMetrics = function(num_items, domain, completeness, redundancy) {
    if (num_items === 0) {
        return { dc: -1, dr: -1, vc: '--', vr: '--', cc: null, cr: null };
    }

    if (!domain) {
        return { dc: -1, dr: -1, vc: '***', vr: '***',
                cc: BIN_DEFAULTS.COLORS.COMPLETENESS.ERROR,
                cr: BIN_DEFAULTS.COLORS.COMPLETENESS.ERROR };
    }

    if (domain === "blank") {
        return { dc: -1, dr: -1, vc: '??', vr: '??',
                cc: BIN_DEFAULTS.COLORS.COMPLETENESS.BLANK,
                cr: BIN_DEFAULTS.COLORS.COMPLETENESS.BLANK };
    }

    if (domain === "mixed") {
        return { dc: -1, dr: -1, vc: 'xx', vr: 'xx',
                cc: BIN_DEFAULTS.COLORS.COMPLETENESS.MIXED,
                cr: BIN_DEFAULTS.COLORS.COMPLETENESS.MIXED };
    }

    if (completeness != null && redundancy != null) {
        const cc = completeness > 80 ? BIN_DEFAULTS.COLORS.COMPLETENESS.HIGH :
                   completeness > 50 ? BIN_DEFAULTS.COLORS.COMPLETENESS.MEDIUM : null;

        const cr = redundancy > 50 ? BIN_DEFAULTS.COLORS.REDUNDANCY.HIGH :
                   redundancy > 7 ? BIN_DEFAULTS.COLORS.REDUNDANCY.MEDIUM : null;

        return {
            dc: completeness,
            dr: redundancy,
            vc: completeness.toFixed(1),
            vr: redundancy.toFixed(1),
            cc: cc,
            cr: cr
        };
    }

    return { dc: -1, dr: -1, vc: '--', vr: '--', cc: null, cr: null };
};

/**
 * Update completeness display in UI
 * @private
 */
Bins.prototype._updateCompletenessDisplay = function(bin_row, metrics) {
    const completenessCell = bin_row.querySelector('td.completeness');
    const redundancyCell = bin_row.querySelector('td.redundancy');

    if (completenessCell) {
        completenessCell.setAttribute('data-value', metrics.dc);
        const input = completenessCell.querySelector('input');
        if (input) {
            input.style.background = metrics.cc || "";
            input.value = metrics.vc;
        }
    }

    if (redundancyCell) {
        redundancyCell.setAttribute('data-value', metrics.dr);
        const input = redundancyCell.querySelector('input');
        if (input) {
            input.style.background = metrics.cr || "";
            input.value = metrics.vr;
        }
    }
};

/**
 * Fetch taxonomy data
 * @private
 */
Bins.prototype._fetchTaxonomyData = function(bin_id, bin_name) {
    const collection_data = {};
    collection_data[bin_name] = [];

    for (const node of this.selections[bin_id].values()) {
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
            this._processTaxonomyData(bin_id, bin_name, data);
        }
    });
};

/**
 * Process taxonomy data response
 * @private
 */
Bins.prototype._processTaxonomyData = function(bin_id, bin_name, data) {
    if (data.status && data.status !== 0) {
        if ($('#estimate_taxonomy').is(':checked')) {
            toastr.error(`"${data.message}", the server said.`, "The anvi'o headquarters is upset");
        }
        $('#estimate_taxonomy').prop('checked', false);
        toggleTaxonomyEstimation();
        return;
    }

    if (!data.hasOwnProperty(bin_name)) return;

    const taxonomyLevels = ["t_domain", "t_phylum", "t_class", "t_order", "t_family", "t_genus", "t_species"];
    const taxonomyLabel = this.container.querySelector(`span.taxonomy-name-label[bin-id="${bin_id}"]`);

    if (!taxonomyLabel) return;

    // Default to unknown
    taxonomyLabel.innerHTML = " (?) Unknown";

    // Find the most specific taxonomy level
    for (let i = taxonomyLevels.length - 1; i >= 0; i--) {
        const level = taxonomyLevels[i];
        const taxonomyValue = data[bin_name]['consensus_taxonomy'][level];

        if (taxonomyValue !== null) {
            const levelLabel = level.split('_')[1][0];
            taxonomyLabel.innerHTML = ` (${levelLabel}) ${taxonomyValue.replace('_', ' ')}`;
            break;
        }
    }
};

// ============================================================================
// Data Management
// ============================================================================

/**
 * Get node labels for a bin
 * @param {number} bin_id - Bin ID
 * @returns {Array} Array of node labels
 */
Bins.prototype.GetBinNodeLabels = function(bin_id) {
    const node_labels = [];

    for (const node of this.selections[bin_id].values()) {
        if (node.IsLeaf()) {
            node_labels.push(node.label);
        }
    }

    return node_labels;
};

/**
 * Migrate collection after tree redraw
 */
Bins.prototype.MigrateCollection = function() {
    const collection = this.ExportCollection(true);

    if (!collection || !collection.data) return;

    const sortedBinIds = Object.keys(collection.data).sort();

    for (const bin_id of sortedBinIds) {
        const nodeLabels = collection.data[bin_id];
        this.selections[bin_id] = new Set(
            nodeLabels.map(label => drawer.tree.GetLeafByName(label))
        );
    }
};

/**
 * Import a collection
 * @param {Object} collection - Collection data
 * @param {number} threshold - Minimum threshold for bin size
 */
Bins.prototype.ImportCollection = function(collection, threshold = 1000) {
    this.keepHistory = false;
    let bins_cleared = false;

    const sortedBinNames = Object.keys(collection['data']).sort();

    for (const bin_name of sortedBinNames) {
        const binData = this._prepareImportedBinData(bin_name, collection, threshold);

        if (binData && binData.nodes.length > 0) {
            if (!bins_cleared) {
                this._clearAllBins();
                bins_cleared = true;
            }

            this._createImportedBin(bin_name, binData, collection);
        }
    }

    if (bins_cleared) {
        this.RebuildIntersections();
        this.UpdateBinsWindow();
        this.RedrawBins();
    }

    this.keepHistory = true;
};

/**
 * Prepare imported bin data
 * @private
 */
Bins.prototype._prepareImportedBinData = function(bin_name, collection, threshold) {
    const nodes = [];
    let sum_length = 0;

    for (const item of collection['data'][bin_name]) {
        if (mode !== 'full') {
            nodes.push(item);
        } else if (typeof item_lengths[item] !== 'undefined') {
            nodes.push(item);
            sum_length += item_lengths[item];
        }
    }

    const meetsThreshold = (mode !== 'full' || threshold === 0 || sum_length >= threshold);

    return meetsThreshold ? { nodes, sum_length } : null;
};

/**
 * Clear all bins for import
 * @private
 */
Bins.prototype._clearAllBins = function() {
    this.bin_counter = 0;
    this.selections = {};
    this.container.innerHTML = '';
    this.history = [];
};

/**
 * Create imported bin
 * @private
 */
Bins.prototype._createImportedBin = function(bin_name, binData, collection) {
    const bin_id = this.bin_counter;

    this.selections[bin_id] = new Set(
        binData.nodes.map(node_label => drawer.tree.GetLeafByName(node_label))
    );

    const bin_color = collection['colors'][bin_name] || randomColor();
    this.NewBin(bin_id, { 'name': bin_name, 'color': bin_color });
    this.bin_counter++;
};

/**
 * Export collection data
 * @param {boolean} use_bin_id - Use bin IDs instead of names
 * @returns {Object} Collection data
 */
Bins.prototype.ExportCollection = function(use_bin_id = false) {
    const data = {};
    const colors = {};

    const binRows = this.container.querySelectorAll('tr.bin-row');

    for (const tr of binRows) {
        const bin_id = tr.getAttribute('bin-id');
        const bin_name = tr.querySelector('.bin-name').value;
        const bin_color = tr.querySelector('.colorpicker').getAttribute('color');
        const items = [];

        for (const node of this.selections[bin_id].values()) {
            if (node === null) return null;
            if (node.IsLeaf()) {
                items.push(node.label);
            }
        }

        if (items.length > 0) {
            const key = use_bin_id ? bin_id : bin_name;
            data[key] = items;
            colors[key] = bin_color;
        }
    }

    return { 'data': data, 'colors': colors };
};

/**
 * Export single bin
 * @param {number} bin_id_to_export - Bin ID to export
 * @param {boolean} use_bin_id - Use bin IDs instead of names
 * @returns {Object} Bin data
 */
Bins.prototype.ExportBin = function(bin_id_to_export, use_bin_id = false) {
    const binRows = this.container.querySelectorAll('tr.bin-row');

    for (const tr of binRows) {
        const bin_id = tr.getAttribute('bin-id');

        if (bin_id == bin_id_to_export) {
            const bin_name = tr.querySelector('.bin-name').value;
            const bin_color = tr.querySelector('.colorpicker').getAttribute('color');
            const items = [];

            for (const node of this.selections[bin_id].values()) {
                if (node.IsLeaf()) {
                    items.push(node.label);
                }
            }

            return {
                'items': items,
                'color': bin_color,
                'bin_name': bin_name
            };
        }
    }

    return null;
};

// ============================================================================
// Highlighting
// ============================================================================

/**
 * Highlight specific items
 * @param {Array|string} item_list - Items to highlight
 */
Bins.prototype.HighlightItems = function(item_list) {
    if (!Array.isArray(item_list)) {
        item_list = [item_list];
    }

    this.higlighted_items = [];

    for (const name of item_list) {
        const node = drawer.tree.GetLeafByName(name);
        if (node) {
            this.higlighted_items.push(name);
        } else {
            console.error(`Received highlight request for non-existent node: ${name}`);
        }
    }

    this.RedrawBins();
};

/**
 * Clear highlighted items
 */
Bins.prototype.ClearHighlightedItems = function() {
    this.higlighted_items = [];
    this.RedrawBins();
};

// ============================================================================
// Drawing and Rendering
// ============================================================================

/**
 * Redraw line colors for all bins
 */
Bins.prototype.RedrawLineColors = function() {
    for (const bin_id in this.selections) {
        if (this.selections[bin_id].size === 0) continue;

        const bin_color = this.GetBinColor(bin_id);
        for (const node of this.selections[bin_id].values()) {
            if (node === null) return;
            node.SetColor(bin_color);
        }
    }
};

/**
 * Draw inverted nodes (nodes not in any bin)
 * @param {Array} leaf_list - List of leaf assignments
 * @param {number} rect_width - Rectangle width for phylogram
 */
Bins.prototype.DrawInvertedNodes = function(leaf_list, rect_width) {
    const inverse_fill_opacity = $('#inverse_fill_opacity').val();
    const inverse_color = document.getElementById('inverse_color').getAttribute('color');
    const nodes_for_inversion = [];
    let calculatedRectX = 0;

    // Find nodes not assigned to bins
    for (let i = 0; i < leaf_list.length; i++) {
        if (leaf_list[i] === -1) {
            const leafNode = drawer.tree.leaves.filter(leaf => leaf.order === i);
            nodes_for_inversion.push(leafNode);
        }
    }

    // Calculate starting X position for rectangles
    nodes_for_inversion.forEach(node => {
        const p = node[0];
        if (p && p.xy.x > calculatedRectX) {
            calculatedRectX = p.xy.x;
        }
    });

    // Draw inverted node backgrounds
    nodes_for_inversion.forEach((node, idx) => {
        if (idx === nodes_for_inversion.length - 1) {
            return; // Skip last node (no border nodes)
        }

        const p = node[0];
        const [p1, p2] = p.GetBorderNodes();

        if (tree_type === 'circlephylogram') {
            this._drawInvertedCircle(p, p1, p2, idx, inverse_color, inverse_fill_opacity);
        } else {
            this._drawInvertedRectangle(p, idx, calculatedRectX, rect_width,
                                       inverse_color, inverse_fill_opacity);
        }
    });
};

/**
 * Draw inverted circle for circlephylogram
 * @private
 */
Bins.prototype._drawInvertedCircle = function(p, p1, p2, idx, color, opacity) {
    const pie = drawPie(
        'bin',
        `bin_outer_${idx}`,
        p1.angle - p1.size / 2,
        p2.angle + p2.size / 2,
        distance(p.backarc, {'x': 0, 'y': 0}),
        total_radius,
        (p2.angle - p1.angle + (p1.size / 2) + (p2.size / 2) > Math.PI) ? 1 : 0,
        color,
        opacity,
        false
    );

    pie.setAttribute('vector-effect', 'non-scaling-stroke');
    pie.setAttribute('stroke-opacity', opacity);
};

/**
 * Draw inverted rectangle for phylogram
 * @private
 */
Bins.prototype._drawInvertedRectangle = function(p, idx, rectX, width, color, opacity) {
    const rect = drawPhylogramRectangle(
        'bin',
        `bin_background_${idx}`,
        rectX,
        p.xy.y,
        p.size,
        width,
        color,
        opacity,
        true
    );

    rect.setAttribute('vector-effect', 'non-scaling-stroke');
    rect.setAttribute('stroke-opacity', '1');
};

/**
 * Main redraw function for bins
 */
Bins.prototype.RedrawBins = function() {
    if (!drawer || !this.allowRedraw) return;

    // Build leaf assignment list
    const leaf_list = this._buildLeafAssignmentList();

    // Group consecutive bins
    const bins_to_draw = this._groupConsecutiveBins(leaf_list);

    // Clear existing bin drawings
    this._clearBinDrawings();

    // Get drawing settings
    const settings = this._getDrawingSettings();

    // Draw bins
    this._drawBins(bins_to_draw, settings);

    // Draw highlighted items
    this._drawHighlightedItems(settings);

    // Draw inverted nodes if enabled
    if (settings.invert_shade) {
        const rect_width = settings.show_grid ? settings.width_with_grids : settings.width_no_grids;
        this.DrawInvertedNodes(leaf_list, rect_width);
    }
};

/**
 * Build leaf assignment list
 * @private
 */
Bins.prototype._buildLeafAssignmentList = function() {
    const leaf_list = new Array(drawer.tree.leaves.length + 1).fill(-1);

    for (const bin_id in this.selections) {
        for (const node of this.selections[bin_id].values()) {
            if (node === null || typeof node === 'undefined') {
                this.selections[bin_id].delete(node);
                continue;
            }

            if (node.IsLeaf() && !node.collapsed) {
                leaf_list[node.order] = bin_id;
                node.SetColor(this.GetBinColor(bin_id));
            }
        }
    }

    return leaf_list;
};

/**
 * Group consecutive bins for drawing
 * @private
 */
Bins.prototype._groupConsecutiveBins = function(leaf_list) {
    const bins_to_draw = [];
    let prev_value = leaf_list[0];
    let prev_start = 0;

    for (let i = 1; i < leaf_list.length; i++) {
        if (prev_value !== leaf_list[i]) {
            if (prev_value !== -1) {
                bins_to_draw.push([prev_start, i - 1, prev_value]);
            }
            prev_start = i;
        }
        prev_value = leaf_list[i];
    }

    return bins_to_draw;
};

/**
 * Clear existing bin drawings
 * @private
 */
Bins.prototype._clearBinDrawings = function() {
    const binElement = document.getElementById('bin');
    if (binElement) {
        while (binElement.hasChildNodes()) {
            binElement.removeChild(binElement.lastChild);
        }
    }
};

/**
 * Get drawing settings from UI
 * @private
 */
Bins.prototype._getDrawingSettings = function() {
    return {
        show_grid: $('#show_grid_for_bins')[0].checked,
        show_shade: $('#show_shade_for_bins')[0].checked,
        invert_shade: $('#invert_shade_for_bins')[0].checked,
        shade_fill_opacity: $('#shade_fill_opacity').val(),
        grid_color: document.getElementById('grid_color').getAttribute('color'),
        grid_width: $('#grid_width').val(),
        show_bin_labels: $('#show_bin_labels')[0].checked,
        bin_labels_font_size: parseFloat($('#bin_labels_font_size').val()) || 32,
        autorotate_bin_labels: $('#autorotate_bin_labels')[0].checked,
        bin_labels_angle: $('#bin_labels_angle').val(),
        background_starts_from_branch: $('#begins_from_branch').is(':checked'),
        outer_ring_size: parseFloat($('#outer-ring-height').val()),
        outer_ring_margin: parseFloat($('#outer-ring-margin').val()),
        width_with_grids: 0,
        width_no_grids: 0
    };
};

/**
 * Draw all bins
 * @private
 */
Bins.prototype._drawBins = function(bins_to_draw, settings) {
    for (const [start_idx, end_idx, bin_id] of bins_to_draw) {
        const start = drawer.tree.leaves[start_idx];
        const end = drawer.tree.leaves[end_idx];

        const binColor = this.GetBinColor(bin_id);
        if (!binColor) {
            console.error(`Color not found for bin_id ${bin_id}`);
            continue;
        }

        const intersection = settings.background_starts_from_branch
            ? this._findCommonAncestor(start, end)
            : null;

        if (tree_type === 'circlephylogram') {
            this._drawCircleBin(start, end, bin_id, binColor, settings, intersection);
        } else {
            const widths = this._drawRectangleBin(start, end, bin_id, binColor, settings, intersection);
            settings.width_with_grids = widths.width_with_grids;
            settings.width_no_grids = widths.width_no_grids;
        }
    }
};

/**
 * Find common ancestor of two nodes
 * @private
 */
Bins.prototype._findCommonAncestor = function(start, end) {
    const startAncestors = new Set(start.GetAncestors());
    const endAncestors = new Set(end.GetAncestors());
    const intersection = [...startAncestors].filter(x => endAncestors.has(x))[0];

    if (!intersection) {
        throw new Error(`Nodes ${start.id} and ${end.id} do not have a common ancestor`);
    }

    return intersection;
};

/**
 * Draw circle-style bin (circlephylogram)
 * @private
 */
Bins.prototype._drawCircleBin = function(start, end, bin_id, color, settings, intersection) {
    const binName = $(`#bin_name_${bin_id}`).val();

    // Outer ring
    drawPie('bin', `bin_outer_${bin_id}`,
        start.angle - start.size / 2,
        end.angle + end.size / 2,
        total_radius + settings.outer_ring_margin,
        total_radius + settings.outer_ring_margin + settings.outer_ring_size,
        (end.angle - start.angle + (start.size / 2) + (end.size / 2) > Math.PI) ? 1 : 0,
        color, 1, true);

    // Label if enabled
    if (settings.show_bin_labels) {
        this._drawCircleLabel(start, end, binName, color, settings);
    }

    // Background
    const backgroundStart = intersection ? intersection.radius : beginning_of_layers;
    const pie = drawPie('bin', `bin_background_${bin_id}`,
        start.angle - start.size / 2,
        end.angle + end.size / 2,
        backgroundStart,
        settings.show_grid ? total_radius + settings.outer_ring_margin + settings.outer_ring_size : total_radius,
        (Math.abs(end.angle - start.angle) + start.size / 2 + end.size / 2 > Math.PI) ? 1 : 0,
        color,
        settings.show_grid ? 0 : settings.shade_fill_opacity,
        false);

    this._applyBinStyles(pie, settings);
};

/**
 * Draw rectangle-style bin (phylogram)
 * @private
 */
Bins.prototype._drawRectangleBin = function(start, end, bin_id, color, settings, intersection) {
    const binName = $(`#bin_name_${bin_id}`).val();
    const height = end.xy.y + end.size / 2 - start.xy.y + start.size / 2;

    // Outer rectangle
    drawPhylogramRectangle('bin', `bin_outer_${bin_id}`,
        total_radius + settings.outer_ring_margin,
        start.xy.y - start.size / 2 + height / 2,
        height, settings.outer_ring_size, color, 1, true);

    // Label if enabled
    if (settings.show_bin_labels) {
        this._drawRectangleLabel(start, end, binName, color, settings);
    }

    // Background
    const backgroundStart = intersection ? intersection.xy.x : beginning_of_layers;
    const width_with_grids = total_radius + settings.outer_ring_margin + settings.outer_ring_size - backgroundStart;
    const width_no_grids = total_radius - backgroundStart;

    const rect = drawPhylogramRectangle('bin', `bin_background_${bin_id}`,
        backgroundStart,
        start.xy.y - start.size / 2 + height / 2,
        height,
        settings.show_grid ? width_with_grids : width_no_grids,
        color,
        settings.show_grid ? 0 : settings.shade_fill_opacity,
        false);

    this._applyBinStyles(rect, settings);

    return { width_with_grids, width_no_grids };
};

/**
 * Apply styles to bin element
 * @private
 */
Bins.prototype._applyBinStyles = function(element, settings) {
    if (settings.show_grid && !settings.show_shade) {
        element.setAttribute('vector-effect', 'non-scaling-stroke');
        element.setAttribute('stroke-opacity', '1');
        element.setAttribute('stroke-width', settings.grid_width);
        element.setAttribute('stroke', settings.grid_color);
    } else if (settings.show_grid && settings.show_shade) {
        element.setAttribute('vector-effect', 'non-scaling-stroke');
        element.setAttribute('stroke-opacity', '1');
        element.setAttribute('stroke-width', settings.grid_width);
        element.setAttribute('stroke', settings.grid_color);
        element.setAttribute('fill-opacity', settings.shade_fill_opacity);
    } else if (!settings.show_grid && !settings.show_shade) {
        element.setAttribute('vector-effect', 'non-scaling-stroke');
        element.setAttribute('stroke-opacity', '0');
        element.setAttribute('stroke-width', 0);
        element.setAttribute('stroke', settings.grid_color);
        element.setAttribute('fill-opacity', 0);
    }
};

/**
 * Draw circle bin label
 * @private
 */
Bins.prototype._drawCircleLabel = function(start, end, binName, color, settings) {
    const labelRadius = total_radius + settings.outer_ring_margin * 1.5 +
                       settings.outer_ring_size * (this.higlighted_items.length > 0 ? 2 : 1);
    const labelAngle = (end.angle + end.size / 2 + start.angle - start.size / 2) / 2;

    let align = 'left';
    let displayAngle = labelAngle * 180.0 / Math.PI;

    if (labelAngle > Math.PI / 2.0 && labelAngle < 1.5 * Math.PI) {
        align = 'right';
        displayAngle += 180.0;
    }

    const px = labelRadius * Math.cos(labelAngle) -
              Math.cos(Math.PI / 2 + labelAngle) * (settings.bin_labels_font_size / 3) *
              (align === 'right' ? 1 : -1);
    const py = labelRadius * Math.sin(labelAngle) -
              Math.sin(Math.PI / 2 + labelAngle) * (settings.bin_labels_font_size / 3) *
              (align === 'right' ? 1 : -1);

    drawRotatedText('bin', {x: px, y: py},
        binName.replaceAll("_", " "),
        settings.autorotate_bin_labels ? displayAngle : settings.bin_labels_angle,
        align,
        settings.bin_labels_font_size + "px",
        "HelveticaNeue-CondensedBold, Helvetica Neue, Helvetica, sans-serif",
        color, 0, 'baseline');
};

/**
 * Draw rectangle bin label
 * @private
 */
Bins.prototype._drawRectangleLabel = function(start, end, binName, color, settings) {
    const labelX = total_radius + settings.outer_ring_margin * 1.5 +
                  settings.outer_ring_size * (this.higlighted_items.length > 0 ? 2 : 1);
    const labelY = (start.xy.y - start.size / 2 + end.xy.y + end.size / 2) / 2 +
                  (settings.bin_labels_font_size / 3);

    drawRotatedText('bin', {y: labelY, x: labelX},
        binName.replace("_", " "),
        settings.autorotate_bin_labels ? 0 : settings.bin_labels_angle,
        'left',
        settings.bin_labels_font_size + "px",
        "HelveticaNeue-CondensedBold, Helvetica Neue, Helvetica, sans-serif",
        color, 0, 'baseline');
};

/**
 * Draw highlighted items
 * @private
 */
Bins.prototype._drawHighlightedItems = function(settings) {
    const highlightColor = document.getElementById('picker_highlight').getAttribute('color');

    for (const name of this.higlighted_items) {
        const node = drawer.tree.GetLeafByName(name);

        if (tree_type === 'circlephylogram') {
            drawPie('bin', 'bin_outer_1',
                node.angle - node.size / 2,
                node.angle + node.size / 2,
                total_radius + settings.outer_ring_margin + settings.outer_ring_size,
                total_radius + settings.outer_ring_margin + settings.outer_ring_size * 2,
                (node.size / 2 > Math.PI) ? 1 : 0,
                highlightColor, 1, true);
        } else {
            drawPhylogramRectangle('bin', 'bin_outer_1',
                total_radius + settings.outer_ring_margin + settings.outer_ring_size,
                node.xy.y,
                node.size,
                settings.outer_ring_size,
                highlightColor, 1, true);
        }
    }
};

/**
 * Rebuild bin intersections
 */
Bins.prototype.RebuildIntersections = function() {
    if (!this.allowRedraw) return;

    for (const bin_id in this.selections) {
        if (this.selections[bin_id].size === 0) continue;

        const bin_color = this.GetBinColor(bin_id);

        // Remove non-leaf nodes without children in bin
        this._removeOrphanedNodes(bin_id);

        // Add parent nodes where both children are in bin
        this._addParentNodes(bin_id, bin_color);
    }
};

/**
 * Remove orphaned non-leaf nodes
 * @private
 */
Bins.prototype._removeOrphanedNodes = function(bin_id) {
    let removed = true;

    while (removed) {
        removed = false;

        for (const node of this.selections[bin_id].values()) {
            if (node === null) return;
            if (node.IsLeaf()) continue;

            // Check if any child is in this bin
            let any_child_in_bin = false;
            let child = node.child;

            while (child && child.sibling) {
                any_child_in_bin = any_child_in_bin || this.selections[bin_id].has(child);
                child = child.sibling;
            }

            if (!any_child_in_bin) {
                this.selections[bin_id].delete(node);
                node.ResetColor();
                removed = true;
            }
        }
    }
};

/**
 * Add parent nodes where appropriate
 * @private
 */
Bins.prototype._addParentNodes = function(bin_id, bin_color) {
    let inserted = true;

    while (inserted) {
        inserted = false;

        for (const node of this.selections[bin_id].values()) {
            const parent = node.ancestor;

            if (!parent || this.selections[bin_id].has(parent)) continue;

            if (node.sibling && this.selections[bin_id].has(node.sibling)) {
                this.selections[bin_id].add(parent);
                parent.SetColor(bin_color);
                inserted = true;
            }
        }
    }
};
