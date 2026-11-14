/**
 * Bin curation window for the anvi'o interactive interface (manual mode).
 *
 * When the user clicks on the number of items of a bin in the Bins panel, this window shows
 * the items in the bin with checkboxes, and (if the user provided a FASTA file and the relevant
 * option in the Bins panel is checked) the sequence of each item. The user can select items,
 * and stage their removal from the bin or their move to another bin. Nothing happens to the
 * actual bins until the user clicks 'Apply', and 'Discard' throws away all staged changes.
 *
 *  Copyright 2015-2026, The anvi'o project (http://anvio.org)
 *
 *  Anvi'o is a free software. You can redistribute this program
 *  and/or modify it under the terms of the GNU General Public
 *  License as published by the Free Software Foundation, either
 *  version 3 of the License, or (at your option) any later version.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with anvi'o. If not, see <http://opensource.org/licenses/GPL-3.0>.
 *
 *  @license GPL-3.0+ <http://opensource.org/licenses/GPL-3.0>
 */

const BIN_CURATION_MODAL_ID = 'bin-curation-modal';

// aligning sequences only makes sense if they are shown
$(document).ready(function() {
    $('#show_sequences_in_bin_items').on('change', function() {
        $('#align_sequences_in_bin_items').prop('disabled', !this.checked)
                                          .parent().css('opacity', this.checked ? 1 : 0.5);
    });
});

// the state of the currently open curation window (there is only one at a time)
let binCurationState = null;


function _binCurationEscape(text) {
    return String(text)
        .replace(/&/g, '&amp;')
        .replace(/</g, '&lt;')
        .replace(/>/g, '&gt;')
        .replace(/"/g, '&quot;')
        .replace(/'/g, '&#39;');
}


/**
 * Given an alignment for a superset of items, returns the alignment for a subset of them by
 * dropping columns that are made of gaps only once the other items are gone.
 */
function _binCurationSubAlignment(alignment, items) {
    const seqs = items.map(item => alignment[item]);
    const length = seqs.length ? seqs[0].length : 0;
    const keep = [];

    for (let i = 0; i < length; i++) {
        if (seqs.some(seq => seq[i] !== '-')) {
            keep.push(i);
        }
    }

    const result = {};
    items.forEach((item, j) => {
        result[item] = keep.map(i => seqs[j][i]).join('');
    });

    return result;
}


/**
 * Entry point: show the curation window for a bin.
 * @param {number} bin_id - The bin to curate
 * @param {Object} cached - Optional {sequences: {item: sequence}, aligned: bool, aligner: str} to
 *                          avoid re-fetching (and re-aligning) sequences after changes are applied
 */
function showBinCurationDialog(bin_id, cached = null) {
    const bin_info = bins.ExportBin(bin_id);

    if (!bin_info) {
        return;
    }

    if (!bin_info.items || bin_info.items.length === 0) {
        toastr.warning('There are no items in this bin yet, so there is nothing to show.', "The anvi'o headquarters has a note");
        return;
    }

    // order items the way they appear in the tree, so the order in this window makes sense
    const leaves = (typeof drawer !== 'undefined' && drawer && drawer.tree) ? drawer.tree.label_to_leaves : {};
    const items = bin_info.items.slice().sort((a, b) => {
        const order_a = leaves[a] ? leaves[a].order : Infinity;
        const order_b = leaves[b] ? leaves[b].order : Infinity;
        return order_a - order_b;
    });

    const showSequences = sequences_available && $('#show_sequences_in_bin_items').is(':checked');
    const alignSequences = showSequences && $('#align_sequences_in_bin_items').is(':checked');

    if (!showSequences) {
        _buildBinCurationDialog(bin_id, bin_info, items, null);
        return;
    }

    // if items were only removed from this bin since the last time we fetched its sequences, we
    // can reuse them (an alignment of a subset of sequences is the alignment of all of them
    // minus the columns that are gaps only in that subset)
    if (cached && cached.aligned === alignSequences && items.every(item => cached.sequences.hasOwnProperty(item))) {
        const sequences = cached.aligned ? _binCurationSubAlignment(cached.sequences, items) : cached.sequences;
        _buildBinCurationDialog(bin_id, bin_info, items, {sequences: sequences, aligned: cached.aligned, aligner: cached.aligner});
        return;
    }

    // we open the window right away with a placeholder for sequences, and fill it once the server
    // responds. NOTE: we do not use the 'waitingDialog' here on purpose: a second modal that is
    // hidden while it is still fading in can leave its backdrop behind, which blocks the entire
    // interface after the curation window is closed.
    _buildBinCurationDialog(bin_id, bin_info, items, {loading: true, aligning: alignSequences});

    const state = binCurationState;

    $.ajax({
        type: 'POST',
        url: '/data/get_sequences_for_splits',
        data: {'split_names': JSON.stringify(items), 'align': alignSequences ? 'true' : 'false'},
        success: function(response) {
            if (typeof response === 'string') {
                response = JSON.parse(response);
            }

            if (response.status != 0) {
                _setBinCurationSequencesError(state, `"${response.message}", the server said.`);
                return;
            }

            _setBinCurationSequences(state, {sequences: response.sequences,
                                             aligned: response.aligned,
                                             aligner: response.aligner});
        },
        error: function() {
            _setBinCurationSequencesError(state, 'Anvi\'o failed to fetch sequences for this bin :/');
        }
    });
}


/**
 * Fills the sequences section of the curation window once sequences arrive, as long as the
 * window they were requested for is still the one that is open.
 */
function _setBinCurationSequences(state, sequencesInfo) {
    if (binCurationState !== state) {
        return;
    }

    state.sequences = sequencesInfo.sequences;
    state.aligned = sequencesInfo.aligned;
    state.aligner = sequencesInfo.aligner;

    const $modal = _getBinCurationModal();
    $modal.find('.bin-curation-sequences-section').html(_buildBinCurationSequencesSection(state.items, state.sequences, state.aligned ? state.aligner : null));

    // carry over the current selections and staged changes to the new checkboxes
    state.items.forEach((item, idx) => {
        const checked = $modal.find(`.bin-curation-list .bin-curation-check[data-idx="${idx}"]`).is(':checked');
        $modal.find(`.bin-curation-seq-table .bin-curation-check[data-idx="${idx}"]`).prop('checked', checked);
    });

    for (const [idx, change] of state.staged.entries()) {
        _renderBinCurationStagedItem(idx, change);
    }
}


function _setBinCurationSequencesError(state, message) {
    if (binCurationState !== state) {
        return;
    }

    _getBinCurationModal().find('.bin-curation-sequences-section').html(`
        <div class="bin-modal-header" style="background: #e8f5e978; margin-top: 25px;">
            <span style="font-size: large;">Sequences</span>
        </div>
        <p style="color: #a94442;">${_binCurationEscape(message)} Everything else in this window works as usual, though.</p>`);
}


function _buildBinCurationDialog(bin_id, bin_info, items, sequencesInfo) {
    const loading = sequencesInfo && sequencesInfo.loading;
    const sequences = (sequencesInfo && !loading) ? sequencesInfo.sequences : null;

    // remove any curation window that may be lingering around (this should never happen since
    // a new window is only opened after the previous one is fully hidden, but if it does, we
    // also clean up after it so its backdrop does not block the interface)
    if ($(`#${BIN_CURATION_MODAL_ID}`).length) {
        $(`#${BIN_CURATION_MODAL_ID}`).remove();
        _removeStrayModalBackdrops();
    }

    binCurationState = {
        bin_id: bin_id,
        bin_name: bin_info.bin_name,
        items: items,
        sequences: sequences,
        aligned: sequences ? sequencesInfo.aligned : false,
        aligner: sequences ? sequencesInfo.aligner : null,
        // item index -> {type: 'remove'} or {type: 'move', bin_id: X, bin_name: Y, color: Z}
        staged: new Map(),
        allowClose: false,
        afterClose: null
    };

    const readOnlyAttr = interface_read_only ? 'disabled title="Not available in read-only mode"' : '';

    const controls = `
        <div class="bin-curation-controls">
            <div class="bin-curation-controls-row">
                <span class="bin-curation-controls-label">Selections:</span>
                <button type="button" class="btn btn-sm bin-curation-btn-selection" data-action="check-all">Check All</button>
                <button type="button" class="btn btn-sm bin-curation-btn-selection" data-action="check-none">Check None</button>
                <button type="button" class="btn btn-sm bin-curation-btn-selection" data-action="invert">Invert selections</button>
            </div>
            <div class="bin-curation-controls-row">
                <span class="bin-curation-controls-label">Operations:</span>
                <button type="button" class="btn btn-warning btn-sm" data-action="remove" ${readOnlyAttr}>Remove selected from bin</button>
                <button type="button" class="btn btn-warning btn-sm" data-action="move" ${readOnlyAttr}>Move selected to another bin</button>
            </div>
            <div class="bin-curation-controls-row">
                <span class="bin-curation-controls-label">Conclusion:</span>
                <button type="button" class="btn btn-sm bin-curation-btn-discard" data-action="discard">Discard all changes</button>
                <button type="button" class="btn btn-sm bin-curation-btn-apply" data-action="apply" ${readOnlyAttr}>Apply</button>
            </div>
            <div class="bin-curation-status"></div>
        </div>`;

    let sequences_section = '';
    if (sequences) {
        sequences_section = _buildBinCurationSequencesSection(items, sequences, binCurationState.aligned ? binCurationState.aligner : null);
    } else if (loading) {
        sequences_section = `
            <div class="bin-modal-header" style="background: #e8f5e978; margin-top: 25px;">
                <span style="font-size: large;">Sequences (${items.length})</span>
            </div>
            <p><span class="spinner-border spinner-border-sm" role="status"></span>
               ${sequencesInfo.aligning ? `Fetching and aligning ${items.length} sequences (this may take a while) ...` : 'Fetching sequences ...'}</p>`;
    }

    const content = _buildBinCurationSplitsSection(items) +
                    `<div class="bin-curation-sequences-section">${sequences_section}</div>`;

    const title = `Items in "${_binCurationEscape(bin_info.bin_name)}"`;

    const template = `
        <div class="modal fade" id="${BIN_CURATION_MODAL_ID}" role="dialog">
            <div class="modal-dialog modal-dialog-centered" style="pointer-events: all; max-width: 90vw; width: 90vw;">
                <div class="modal-content" style="max-height: 85vh; display: flex; flex-direction: column; position: relative;">
                    <div class="modal-header" style="flex-shrink: 0;">
                        <h4 class="modal-title">${title}</h4>
                        <button type="button" class="close" data-dismiss="modal" aria-hidden="true">${MODALS.CLOSE_BUTTON_HTML}</button>
                    </div>
                    <div style="flex-shrink: 0; padding: 10px 20px 0 20px;">
                        ${controls}
                    </div>
                    <div class="modal-body" style="overflow: auto; flex: 1; min-height: 0; padding: 10px 20px 20px 20px;">
                        ${content}
                    </div>
                    <div class="bin-curation-move-overlay" style="display: none;"></div>
                </div>
            </div>
        </div>`;

    $('body').append(template);

    const $modal = $(`#${BIN_CURATION_MODAL_ID}`);

    // keep checkboxes for the same item in sync across the 'Splits' and 'Sequences' sections
    $modal.on('change', '.bin-curation-check', function() {
        const idx = $(this).attr('data-idx');
        $modal.find(`.bin-curation-check[data-idx="${idx}"]`).prop('checked', this.checked);
        _updateBinCurationStatus();
    });

    $modal.on('click', '[data-action]', function() {
        const action = $(this).attr('data-action');

        switch (action) {
            case 'check-all':  _setBinCurationChecks(() => true); break;
            case 'check-none': _setBinCurationChecks(() => false); break;
            case 'invert':     _setBinCurationChecks((checked) => !checked); break;
            case 'remove':     _stageBinCurationRemoval(); break;
            case 'move':       _showBinCurationMoveOverlay(); break;
            case 'discard':    _discardBinCurationChanges(); break;
            case 'apply':      _applyBinCurationChanges(); break;
            case 'copy-fasta-all':     _copyBinCurationFASTA(this, false); break;
            case 'copy-fasta-checked': _copyBinCurationFASTA(this, true); break;
        }
    });

    // closing the window by any means other than 'Discard' or 'Apply' (the X button, the
    // escape key, or clicking outside) discards changes, so we ask the user first if there
    // are changes that would be lost
    $modal.on('hide.bs.modal', function(event) {
        if (!binCurationState || binCurationState.allowClose || binCurationState.staged.size === 0) {
            return;
        }

        const n = binCurationState.staged.size;
        if (!window.confirm(`You have ${n} pending change${n > 1 ? 's' : ''} that will be discarded if you close this window. Close anyway?`)) {
            event.preventDefault();
            event.stopImmediatePropagation();
        }
    });

    // bootstrap only closes the window with the escape key if the focus is within the window,
    // which is often not the case (e.g., after the move window is closed), so we listen to it
    // at the document level
    const escHandler = function(event) {
        if (event.key === 'Escape') {
            const $overlay = $modal.find('.bin-curation-move-overlay');

            if ($overlay.is(':visible')) {
                $overlay.hide().empty();
            } else {
                $modal.modal('hide');
            }
        }
    };

    $(document).on('keydown.binCuration', escHandler);

    $modal.on('hidden.bs.modal', function() {
        const afterClose = binCurationState ? binCurationState.afterClose : null;

        $(document).off('keydown.binCuration');
        $(this).remove();
        binCurationState = null;

        // make sure nothing is left behind that would block the interface
        _removeStrayModalBackdrops();

        // open whatever comes next only after bootstrap is done with this window
        if (afterClose) {
            setTimeout(afterClose, 0);
        }
    });

    $modal.modal({show: true, backdrop: true, keyboard: false})
          .find('.modal-dialog')
          .draggable({handle: '.modal-header'});

    _updateBinCurationStatus();
}


function _buildBinCurationSplitsSection(items) {
    const rows = items.map((item, idx) => `
        <li class="bin-curation-item" data-idx="${idx}">
            <label style="margin: 0; cursor: pointer;">
                <input type="checkbox" class="bin-curation-check" data-idx="${idx}">
                <span class="bin-curation-item-name">${_binCurationEscape(item)}</span>
            </label>
            <span class="bin-curation-tag"></span>
        </li>`).join('');

    return buildCopyableNamesSection('Splits', items, {
        background: '#ffe4c478',
        contentRenderer: () => `<ul class="bin-curation-list">${rows}</ul>`
    });
}


function _buildBinCurationSequencesSection(items, sequences, aligner) {
    const alignedNote = aligner
        ? `<span style="font-size: 0.85em; font-style: italic;">(aligned with ${_binCurationEscape(aligner)})</span>`
        : '';

    const rows = items.map((item, idx) => `
        <tr class="bin-curation-item" data-idx="${idx}">
            <td class="bin-curation-seq-name">
                <label style="margin: 0; cursor: pointer;">
                    <input type="checkbox" class="bin-curation-check" data-idx="${idx}">
                    <span class="bin-curation-item-name">${_binCurationEscape(item)}</span>
                </label>
                <span class="bin-curation-tag"></span>
            </td>
            <td class="bin-curation-seq">${_binCurationEscape(sequences[item])}</td>
        </tr>`).join('');

    return `
        <div class="bin-modal-header" style="background: #e8f5e978; margin-top: 25px;">
            <span style="font-size: large;">Sequences (${items.length})</span>
            ${alignedNote}
            <button type="button" class="btn btn-primary btn-sm" data-action="copy-fasta-all"
                    title="Copy sequences of all items that remain in this bin (i.e., excluding those staged to be removed or moved) as FASTA">Copy to Clipboard as FASTA</button>
            <span class="copy-feedback" style="min-width: 14px; font-size: 0.9em; line-height: 1;"></span>
            <button type="button" class="btn btn-outline-primary btn-sm" data-action="copy-fasta-checked"
                    title="Copy sequences of checked items only as FASTA">Copy checked as FASTA</button>
            <span class="copy-feedback" style="min-width: 14px; font-size: 0.9em; line-height: 1;"></span>
        </div>
        <div class="bin-curation-seq-container">
            <table class="bin-curation-seq-table">${rows}</table>
        </div>`;
}


function _getBinCurationModal() {
    return $(`#${BIN_CURATION_MODAL_ID}`);
}


/**
 * Returns the indices of checked items that are not already staged for a change.
 */
function _getBinCurationCheckedIndices() {
    const checked = new Set();

    _getBinCurationModal().find('.bin-curation-check:checked').each(function() {
        const idx = parseInt($(this).attr('data-idx'));
        if (!binCurationState.staged.has(idx)) {
            checked.add(idx);
        }
    });

    return Array.from(checked).sort((a, b) => a - b);
}


function _setBinCurationChecks(fn) {
    const $modal = _getBinCurationModal();

    binCurationState.items.forEach((item, idx) => {
        if (binCurationState.staged.has(idx)) {
            return;
        }

        const $checks = $modal.find(`.bin-curation-check[data-idx="${idx}"]`);
        $checks.prop('checked', fn($checks.first().is(':checked')));
    });

    _updateBinCurationStatus();
}


/**
 * Removes modal backdrops that are left behind when no modal is shown, as they would
 * block the entire interface.
 */
function _removeStrayModalBackdrops() {
    if ($('.modal.show').length === 0) {
        $('.modal-backdrop').remove();
        $('body').removeClass('modal-open').css('padding-right', '');
    }
}


function _renderBinCurationStagedItem(idx, change) {
    const $modal = _getBinCurationModal();

    const label = change.type === 'remove'
        ? '→ removed'
        : `→ <span class="bin-curation-tag-color" style="background: ${change.color};"></span>${_binCurationEscape(change.bin_name)}`;

    $modal.find(`.bin-curation-item[data-idx="${idx}"]`).addClass('bin-curation-staged');
    $modal.find(`.bin-curation-item[data-idx="${idx}"] .bin-curation-tag`).html(label);
    $modal.find(`.bin-curation-check[data-idx="${idx}"]`).prop('checked', false).prop('disabled', true);
}


function _stageBinCurationChange(indices, change) {
    for (const idx of indices) {
        binCurationState.staged.set(idx, change);
        _renderBinCurationStagedItem(idx, change);
    }

    _updateBinCurationStatus();
}


function _stageBinCurationRemoval() {
    const indices = _getBinCurationCheckedIndices();

    if (!indices.length) {
        toastr.warning('Please select some items first.');
        return;
    }

    _stageBinCurationChange(indices, {type: 'remove'});
}


function _getOtherBins() {
    const other_bins = [];

    for (const tr of bins.container.querySelectorAll('tr.bin-row')) {
        const bin_id = tr.getAttribute('bin-id');

        if (bin_id == binCurationState.bin_id) {
            continue;
        }

        other_bins.push({
            bin_id: bin_id,
            bin_name: tr.querySelector('.bin-name').value,
            color: tr.querySelector('.colorpicker').getAttribute('color')
        });
    }

    return other_bins;
}


/**
 * Returns new bins the user has staged moves to in this window (they will be created on 'Apply').
 * @returns {Array} [{bin_name, color}]
 */
function _getStagedNewBins() {
    const new_bins = new Map();

    for (const change of binCurationState.staged.values()) {
        if (change.type === 'move' && change.new_bin && !new_bins.has(change.bin_name)) {
            new_bins.set(change.bin_name, {bin_name: change.bin_name, color: change.color});
        }
    }

    return Array.from(new_bins.values());
}


/**
 * Suggests a name for a new bin that is not used by any existing bin, or by any new bin
 * that is already staged in this window.
 */
function _suggestNewBinName() {
    const taken = new Set();

    for (const tr of bins.container.querySelectorAll('tr.bin-row')) {
        taken.add(tr.querySelector('.bin-name').value);
    }

    for (const b of _getStagedNewBins()) {
        taken.add(b.bin_name);
    }

    let n = bins.bin_counter + 1;
    while (taken.has(bins.prefix + n)) {
        n++;
    }

    return bins.prefix + n;
}


function _showBinCurationMoveOverlay() {
    const indices = _getBinCurationCheckedIndices();

    if (!indices.length) {
        toastr.warning('Please select some items first.');
        return;
    }

    const other_bins = _getOtherBins();
    const staged_new_bins = _getStagedNewBins();

    const existing_buttons = other_bins.map(b => `
        <button type="button" class="btn btn-outline-dark btn-sm bin-curation-move-target" data-target-bin-id="${_binCurationEscape(b.bin_id)}">
            <span class="bin-curation-tag-color" style="background: ${b.color};"></span>${_binCurationEscape(b.bin_name)}
        </button>`).join('');

    const staged_new_buttons = staged_new_bins.map(b => `
        <button type="button" class="btn btn-outline-dark btn-sm bin-curation-move-target" data-target-new-bin="${_binCurationEscape(b.bin_name)}"
                title="This bin will be created when you click 'Apply'">
            <span class="bin-curation-tag-color" style="background: ${b.color};"></span>${_binCurationEscape(b.bin_name)} <i>(new)</i>
        </button>`).join('');

    const targets_section = (existing_buttons || staged_new_buttons)
        ? `<div class="bin-curation-move-targets">${existing_buttons}${staged_new_buttons}</div>
           <p style="margin: 15px 0 6px 0;">... or to a new bin:</p>`
        : `<p style="margin: 0 0 6px 0;">There are no other bins yet, but you can move them to a new one:</p>`;

    const $overlay = _getBinCurationModal().find('.bin-curation-move-overlay');

    $overlay.html(`
        <div class="bin-curation-move-window">
            <h5>Move ${indices.length} selected item${indices.length > 1 ? 's' : ''} to ...</h5>
            ${targets_section}
            <div class="bin-curation-move-targets">
                <input type="text" class="form-control form-control-sm bin-curation-new-bin-name" style="width: 220px;"
                       value="${_binCurationEscape(_suggestNewBinName())}">
                <button type="button" class="btn btn-primary btn-sm bin-curation-move-new">Move to this new bin</button>
            </div>
            <p style="font-size: 0.85em; font-style: italic; margin: 6px 0 0 0;">The new bin will be created when you click 'Apply'.</p>
            <div style="text-align: right; margin-top: 15px;">
                <button type="button" class="btn btn-outline-secondary btn-sm bin-curation-move-cancel">Cancel</button>
            </div>
        </div>`).show();

    const $name_input = $overlay.find('.bin-curation-new-bin-name');

    // same rule as the bin name inputs in the Bins panel
    $name_input.on('input', function() {
        this.value = this.value.replaceAll(' ', '_');
    });

    const stage_move_to_new_bin = function() {
        const name = $name_input.val().trim();

        if (!name.length) {
            toastr.warning('Please give the new bin a name.');
            return;
        }

        const existing = other_bins.find(b => b.bin_name === name);
        if (existing || name === binCurationState.bin_name) {
            toastr.warning(`There already is a bin called "${name}". You can pick it from the list, or give the new bin another name.`);
            return;
        }

        // if this name matches a new bin staged earlier, items go to that one
        const staged = staged_new_bins.find(b => b.bin_name === name);
        const color = staged ? staged.color : randomColor({luminosity: 'dark'});

        _stageBinCurationChange(indices, {type: 'move', new_bin: true, bin_id: null, bin_name: name, color: color});
        $overlay.hide().empty();
    };

    $name_input.on('keydown', function(event) {
        if (event.key === 'Enter') {
            event.preventDefault();
            stage_move_to_new_bin();
        }
    });

    $overlay.off('click').on('click', function(event) {
        const $target = $(event.target).closest('.bin-curation-move-target');

        if ($target.length) {
            if ($target.is('[data-target-new-bin]')) {
                const target = staged_new_bins.find(b => b.bin_name === $target.attr('data-target-new-bin'));
                _stageBinCurationChange(indices, {type: 'move', new_bin: true, bin_id: null, bin_name: target.bin_name, color: target.color});
            } else {
                const target = other_bins.find(b => b.bin_id == $target.attr('data-target-bin-id'));
                _stageBinCurationChange(indices, {type: 'move', bin_id: target.bin_id, bin_name: target.bin_name, color: target.color});
            }
            $overlay.hide().empty();
        } else if ($(event.target).closest('.bin-curation-move-new').length) {
            stage_move_to_new_bin();
        } else if ($(event.target).closest('.bin-curation-move-cancel').length || $(event.target).is($overlay)) {
            $overlay.hide().empty();
        }
    });

    $name_input.trigger('focus').trigger('select');
}


function _updateBinCurationStatus() {
    const $modal = _getBinCurationModal();
    const num_items = binCurationState.items.length;
    const num_checked = _getBinCurationCheckedIndices().length;

    let num_removed = 0, num_moved = 0;
    for (const change of binCurationState.staged.values()) {
        if (change.type === 'remove') num_removed++; else num_moved++;
    }

    const pending = num_removed + num_moved;
    const pending_text = pending
        ? `<b>${pending} pending change${pending > 1 ? 's' : ''}</b> (${num_removed} to be removed, ${num_moved} to be moved). Nothing will change until you click 'Apply'.`
        : 'No pending changes.';

    $modal.find('.bin-curation-status').html(`${num_checked} of ${num_items - pending} items selected. ${pending_text}`);
    $modal.find('[data-action="apply"]').text(pending ? `Apply (${pending})` : 'Apply');
}


function _copyBinCurationFASTA(btn, checkedOnly) {
    const indices = checkedOnly
        ? _getBinCurationCheckedIndices()
        : binCurationState.items.map((item, idx) => idx).filter(idx => !binCurationState.staged.has(idx));

    if (!indices.length) {
        toastr.warning(checkedOnly ? 'No items are checked :/' : 'There are no items left in this bin to copy :/');
        return;
    }

    const fasta = indices.map(idx => {
        const item = binCurationState.items[idx];
        return `>${item}\n${binCurationState.sequences[item]}`;
    }).join('\n') + '\n';

    copyTextWithFeedback(btn, fasta);
}


function _discardBinCurationChanges() {
    if (!binCurationState.staged.size) {
        toastr.info('There are no changes to discard.');
        return;
    }

    const $modal = _getBinCurationModal();

    for (const idx of binCurationState.staged.keys()) {
        $modal.find(`.bin-curation-item[data-idx="${idx}"]`).removeClass('bin-curation-staged');
        $modal.find(`.bin-curation-item[data-idx="${idx}"] .bin-curation-tag`).empty();
        $modal.find(`.bin-curation-check[data-idx="${idx}"]`).prop('disabled', false);
    }

    binCurationState.staged.clear();

    _updateBinCurationStatus();
}


function _applyBinCurationChanges() {
    if (interface_read_only) {
        return;
    }

    const state = binCurationState;

    if (!state.staged.size) {
        toastr.info('There are no changes to apply.');
        return;
    }

    const leaves = drawer.tree.label_to_leaves;
    const removals = [];
    const moves = {};
    const moves_to_new_bins = {};
    const missing = [];

    for (const [idx, change] of state.staged.entries()) {
        const item = state.items[idx];
        const node = leaves[item];

        if (!node) {
            missing.push(item);
            continue;
        }

        if (change.type === 'remove') {
            removals.push(node);
        } else if (change.new_bin) {
            if (!moves_to_new_bins.hasOwnProperty(change.bin_name)) {
                moves_to_new_bins[change.bin_name] = {color: change.color, nodes: []};
            }

            moves_to_new_bins[change.bin_name].nodes.push(node);
        } else {
            if (!bins.selections.hasOwnProperty(change.bin_id)) {
                toastr.error(`The bin "${change.bin_name}" does not seem to exist anymore :/ No changes were applied.`);
                return;
            }

            if (!moves.hasOwnProperty(change.bin_id)) {
                moves[change.bin_id] = [];
            }

            moves[change.bin_id].push(node);
        }
    }

    // the user may have created a bin with the same name in the meantime
    const existing_names = new Set(Array.from(bins.container.querySelectorAll('tr.bin-row')).map(tr => tr.querySelector('.bin-name').value));
    const clashing = Object.keys(moves_to_new_bins).filter(name => existing_names.has(name));
    if (clashing.length) {
        toastr.error(`There already is a bin called "${clashing[0]}", so anvi'o can't create a new one with that name :/ No changes were applied.`);
        return;
    }

    if (missing.length) {
        toastr.error(`Anvi'o could not find ${missing.length} of the items in the tree (e.g., "${missing[0]}"), so no changes were applied :/`);
        return;
    }

    // all changes go into a single history transaction, so a single 'undo' reverts them all
    bins.RunAsSingleTransaction(() => {
        if (removals.length) {
            bins.RemoveNode(removals, state.bin_id);
        }

        for (const target_bin_id in moves) {
            bins.AppendNode(moves[target_bin_id], target_bin_id);
        }

        // creating a new bin selects it in the Bins panel, so we remember the selected bin
        // and restore it afterwards
        const selected_bin_id = bins.GetSelectedBinId();

        for (const name in moves_to_new_bins) {
            const new_bin_id = bins.bin_counter;
            bins.selections[new_bin_id] = new Set();
            bins.bin_counter++;

            bins.NewBin(new_bin_id, {name: name, color: moves_to_new_bins[name].color});
            bins.AppendNode(moves_to_new_bins[name].nodes, new_bin_id);
        }

        if (Object.keys(moves_to_new_bins).length && selected_bin_id !== null) {
            $(bins.container).find(`input[name=active_bin][value="${selected_bin_id}"]`).prop('checked', true);
        }
    });

    const num_new_bins = Object.keys(moves_to_new_bins).length;

    toastr.success(`${removals.length} item${removals.length === 1 ? '' : 's'} removed, and ${state.staged.size - removals.length} moved` +
                   (num_new_bins ? ` (${num_new_bins} new bin${num_new_bins > 1 ? 's' : ''} created).` : '.'));

    // refresh the window with the current contents of the bin
    const bin_id = state.bin_id;
    const cached = state.sequences ? {sequences: state.sequences, aligned: state.aligned, aligner: state.aligner} : null;
    const $modal = _getBinCurationModal();

    state.allowClose = true;
    state.afterClose = function() {
        const bin_info = bins.ExportBin(bin_id);

        if (bin_info && bin_info.items.length) {
            showBinCurationDialog(bin_id, cached);
        } else {
            toastr.info('There are no items left in this bin.');
        }
    };

    $modal.modal('hide');
}
