/**
 * Utility functions for anvi'o interactive interface
 *
 *  Authors: Özcan Esen <ozcanesen@gmail.com>
 *           A. Murat Eren <a.murat.eren@gmail.com>
 *
 * Copyright 2015-2021, The anvi'o project (http://anvio.org)
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

const BLAST_DEFAULTS = {MAX_NUM_SEQ: '100',
                        FORMAT_NUM_ORG: '1',
                        CONFIG_DESCR: '2,3,4,5,6,7,8',
                        CLIENT: 'web',
                        SERVICE: 'plain',
                        CMD: 'request',
                        PAGE: 'MegaBlast',
                        MEGABLAST: 'on',
                        WWW_BLAST_TYPE: 'newblast',
                        DEFAULT_PROG: 'megaBlast',
                        SELECTED_PROG_TYPE: 'megaBlast',
                        SAVED_SEARCH: 'true',
                        NUM_DIFFS: '0',
                        NUM_OPTS_DIFFS: '0',
                        PAGE_TYPE: 'BlastSearch',
                        USER_DEFAULT_PROG_TYPE: 'megaBlast'
};

const UNITS = {
    READABLE_NUMBER: ['', 'K', 'M', 'G'],
    SEQ_SIZE: [' K', ' M', ' G', ' T']
};

const MODALS = {
    CLOSE_BUTTON_HTML: '&times;'
};

// ============================================================================
// Event Management
// ============================================================================

/**
 * Emit a custom event
 * @param {string} name - Event name
 * @param {HTMLElement} element - Target element (defaults to document.body)
 */
function emit(name, element = document.body) {
    const event = new Event(name);
    element.dispatchEvent(event);
}

// ============================================================================
// Math Utilities
// ============================================================================

/**
 * Calculate base-10 logarithm
 * @param {number} val - Input value
 * @returns {number} Base-10 logarithm
 */
function log10(val) {
    return Math.log(val) / Math.LN10;
}

/**
 * Check if angle difference is large (> π)
 * @param {number} a - First angle in radians
 * @param {number} b - Second angle in radians
 * @returns {number} 1 if large angle, 0 otherwise
 */
function is_large_angle(a, b) {
    return (Math.abs(b - a) > Math.PI) ? 1 : 0;
}

/**
 * Clamp a number between min and max values
 * @param {number} num - Number to clamp
 * @param {number} min - Minimum value
 * @param {number} max - Maximum value
 * @returns {number} Clamped value
 */
function clamp(num, min, max) {
    return Math.min(Math.max(num, min), max);
}

/**
 * Log progress information with timestamp
 * @param {string} step - Description of current step
 */
function info(step) {
    const timestamp = new Date(Date.now()).toLocaleString().substr(11, 7);
    console.log(`${step} (${timestamp}).`);
}

/**
 * Get object key by its value
 * @param {Object} object - Object to search
 * @param {*} value - Value to find
 * @returns {string|undefined} Key if found
 */
function getKeyByValue(object, value) {
    return Object.keys(object).find(key => object[key] === value);
}

// ============================================================================
// Gene Function Coloring
// ============================================================================

/**
 * Get category code for a functional annotation type
 * @param {Object} geneFunctions - Gene function data
 * @param {string} fn_type - Function type
 * @returns {string|null} Category code or null
 */
function getCagForType(geneFunctions, fn_type) {
    let cag = geneFunctions != null && geneFunctions[fn_type] != null
        ? geneFunctions[fn_type][1]
        : null;

    if (!cag) return null;

    // Handle multiple values - take first one
    const separators = [',', ';', '!!!'];
    for (const sep of separators) {
        if (cag.indexOf(sep) !== -1) {
            cag = cag.substr(0, cag.indexOf(sep));
        }
    }

    return cag;
}

/**
 * Get category code for a gene by genome and gene IDs
 * @param {string} genomeID - Genome identifier
 * @param {string} geneID - Gene identifier
 * @param {string} fn_type - Function type
 * @returns {string|null} Category code or null
 */
function getCagForID(genomeID, geneID, fn_type) {
    const genome = settings.genomeData.genomes.find(x => x[0] == genomeID);
    const functions = genome?.[1]?.genes?.functions?.[geneID];
    return getCagForType(functions, fn_type);
}

/**
 * Add a color row to the function colors table
 * @param {string} label - Row label
 * @param {string} cagCode - Category code
 * @param {string} color - Hex color
 * @param {boolean} prepend - Add to beginning instead of end
 */
function appendColorRow(label, cagCode, color, prepend = false) {
    const code = getCleanCagCode(cagCode);
    const tbody_content = `
        <tr id="picker_row_${code}">
            <td></td>
            <td>
                <div id="picker_${code}" class="colorpicker annotation_color"
                     color="${color}" background-color="${color}"
                     style="background-color: ${color}; margin-right:16px; margin-left:16px">
                </div>
            </td>
            <td>${label}</td>
        </tr>`;

    const method = prepend ? 'prepend' : 'append';
    $(`#tbody_function_colors`)[method](tbody_content);
}

/**
 * Clean category code for use as HTML ID
 * @param {string} code - Raw category code
 * @returns {string} Cleaned code safe for HTML IDs
 */
function getCleanCagCode(code) {
    if (!isNaN(code)) return code;

    const charsToReplace = [' ', '(', ')', ':', '/', '+', '.', "'", '"'];
    let cleanCode = code;

    for (const char of charsToReplace) {
        cleanCode = cleanCode.split(char).join('_');
    }

    return cleanCode;
}

// ============================================================================
// External Database Search Functions
// ============================================================================

/**
 * Search for protein structure in AlphaFold database
 * @param {string} gene_id - Gene identifier
 * @param {string} target - Target type (gene/contig)
 */
function search_protein_structure_in_alphafold(gene_id, target) {
    $.ajax({
        type: 'GET',
        cache: false,
        async: false,
        url: `/data/${target}/${gene_id}`,
        success: function(data) {
            if ('error' in data) {
                toastr.error(data['error'], "", {
                    'timeOut': '0',
                    'extendedTimeOut': '0'
                });
            } else {
                const sequence = data['aa_sequence'];
                const alphafoldUrl = `https://www.alphafold.ebi.ac.uk/search/text/${encodeURIComponent(sequence)}`;
                window.open(alphafoldUrl, '_blank');
            }
        }
    });
}

/**
 * Search gene sequence in remote databases
 * @param {string} gene_id - Gene identifier
 * @param {string} program - BLAST program type
 * @param {string} database - Target database
 * @param {string} target - Target type (gene/contig)
 */
function search_gene_sequence_in_remote_dbs(gene_id, program, database, target) {
    $.ajax({
        type: 'GET',
        cache: false,
        async: false, // Important: maintains direct call chain for popup blocker
        url: `/data/${target}/${gene_id}`,
        success: function(data) {
            if ('error' in data) {
                toastr.error(data['error'], "", {
                    'timeOut': '0',
                    'extendedTimeOut': '0'
                });
            } else {
                const isProtein = (program === 'blastp' || program === 'tblastn');
                const sequenceData = isProtein ? data['aa_sequence'] : data['sequence'];
                const sequence = `>${data['header']}\n${sequenceData}`;
                fire_up_ncbi_blast(sequence, program, database, target);
            }
        }
    });
}

/**
 * Open NCBI BLAST with given sequence
 * @param {string} sequence - FASTA formatted sequence
 * @param {string} program - BLAST program type
 * @param {string} database - Target database
 * @param {string} target - Target type (gene/contig)
 */
function fire_up_ncbi_blast(sequence, program, database, target) {
    const validTargets = ["gene", "contig"];
    if (!validTargets.includes(target)) {
        console.log("fire_up_ncbi_blast: Unrecognized target. Target must be either 'gene' or 'contig'.");
        return;
    }

    const post_variables = {
        ...BLAST_DEFAULTS,
        QUERY: sequence,
        PROGRAM: program || '',
        DATABASE: database || ''
    };

    const blast_window = window.open('about:blank', '_blank');
    const form = document.createElement('form');
    form.action = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi';
    form.method = 'POST';

    // Add form fields
    for (const [name, value] of Object.entries(post_variables)) {
        $(form).append(`<input type="hidden" name="${name}" value="${value}" />`);
    }

    // Add extra field for blastn
    if (program === 'blastn') {
        $(form).append('<input type="hidden" name="BLAST_PROGRAMS" value="megaBlast" />');
    }

    blast_window.document.body.appendChild(form);
    form.submit();
}

// ============================================================================
// URL and Navigation Utilities
// ============================================================================

/**
 * Generate inspection link for various view types
 * @param {Object} options - Link options
 * @param {boolean} options.show_snvs - Show SNVs
 * @param {string} options.type - Inspection type
 * @param {string} options.item_name - Item to inspect
 * @returns {string} Generated URL
 */
function generate_inspect_link(options) {
    const { show_snvs, type, item_name } = options;
    const base_url = window.location.href.split('?')[0];
    let new_url = "";

    if (self == top) {
        // Local anvi'o
        if (base_url.endsWith('index.html')) {
            // On index page
            const urlMap = {
                'inspect_split': `charts.html?id=${item_name}`,
                'inspect_context': `charts.html?id=${item_name}&highlight_gene=true`,
                'inspect_gene': `charts.html?id=${item_name}&highlight_gene=true&gene_mode=true`,
                'inspect_geneclusters': `geneclusters.html?id=${item_name}`
            };
            new_url = urlMap[type] || '';
        } else {
            // On charts or gene cluster page
            new_url = `${base_url}?id=${item_name}`;

            if (type === 'inspect_context') {
                new_url += '&highlight_gene=true';
            } else if (type === 'inspect_gene') {
                new_url += '&highlight_gene=true&gene_mode=true';
            }
        }
    } else {
        // Anvi server
        new_url = base_url.endsWith('/inspect') || base_url.endsWith('/geneclusters')
            ? base_url
            : `${base_url}/${type}`;

        new_url += `?id=${item_name}`;

        const view_key = request_prefix.substr(request_prefix.lastIndexOf('/') + 1);
        if (view_key !== 'no_view_key') {
            new_url += `&view_key=${view_key}`;
        }
    }

    if (type !== 'inspect_geneclusters') {
        new_url += `&show_snvs=${show_snvs}`;
    }

    return new_url;
}

/**
 * Get URL parameter by name
 * @param {string} name - Parameter name
 * @param {string} url - URL to parse (defaults to current)
 * @returns {string|null} Parameter value or null
 */
function getParameterByName(name, url) {
    if (!url) url = window.location.href;

    name = name.replace(/[\[\]]/g, "\\$&");
    const regex = new RegExp("[?&]" + name + "(=([^&#]*)|&|#|$)");
    const results = regex.exec(url);

    if (!results) return null;
    if (!results[2]) return '';

    return decodeURIComponent(results[2].replace(/\+/g, " "));
}

// ============================================================================
// Content Processing
// ============================================================================

/**
 * Render markdown content with custom link handling
 * @param {string} content - Markdown content
 * @returns {string} Rendered HTML
 */
function renderMarkdown(content) {
    const renderer = new marked.Renderer();

    renderer.link = function(hrefObj, title, text) {
        const href = typeof hrefObj === 'string' ? hrefObj : hrefObj.href;

        if (typeof href !== 'string') {
            console.error('Expected href to be a string, got:', hrefObj);
            return `<a href="#">${text}</a>`;
        }

        // Handle special item:// links
        if (href.startsWith('item://')) {
            const item_name = href.split('//')[1];
            let html = `<a href="#" class="item-link">${text}<span class="tooltiptext">
                <span href="#" onclick="bins.HighlightItems('${item_name}');">HIGHLIGHT</span>`;

            if (mode === 'full' || mode === 'pan') {
                const target = mode === 'pan' ? 'inspect_gene_cluster' : 'inspect_contig';
                html += ` | <span href="#" onclick="context_menu_target_id = label_to_node_map['${item_name}'].id;
                                                    menu_callback('${target}');">INSPECT</span>`;
            }

            return html + '</span></a>';
        }

        return `<a target="_blank" href="${href}" title="${title}">${text}</a>`;
    };

    return marked.parse(content, { renderer });
}

// ============================================================================
// Cookie Management
// ============================================================================

/**
 * Get cookie value by name
 * @param {string} name - Cookie name
 * @returns {string|null} Cookie value or null
 */
function getCookie(name) {
    if (!document.cookie) return null;

    const nameEQ = name + '=';
    const cookies = document.cookie.split(';');

    for (const cookie of cookies) {
        const trimmedCookie = cookie.trim();
        if (trimmedCookie.indexOf(nameEQ) === 0) {
            return decodeURIComponent(trimmedCookie.substring(nameEQ.length));
        }
    }

    return null;
}

/**
 * Create or update a cookie
 * @param {string} name - Cookie name
 * @param {string} value - Cookie value
 * @param {number} days - Expiration in days (-1 for 10 years)
 */
function createCookie(name, value, days) {
    let expires = "";
    const date = new Date();

    if (days === -1) {
        // 10 years
        date.setTime(date.getTime() + (10 * 365 * 24 * 60 * 60 * 1000));
    } else if (days) {
        date.setTime(date.getTime() + (days * 24 * 60 * 60 * 1000));
    }

    expires = `; expires=${date.toUTCString()}`;
    document.cookie = `${name}=${value}${expires}; path=/`;
}

// ============================================================================
// Dialog/Modal Functions
// ============================================================================

/**
 * Generic modal dialog creator
 * @param {Object} options - Modal configuration options
 * @private
 */
function _createModalDialog(options) {
    const {
        title,
        content,
        modalClass = 'genericDialog',
        dialogClass = 'modal-dialog',
        noteHTML = null
    } = options;

    const randomID = title.hashCode();

    const noteSection = noteHTML
        ? `<p style="margin: 20px; font-style: italic; flex-shrink: 0;">${noteHTML}</p>`
        : '';

    const template = `
        <div class="modal fade ${modalClass}" id="modal${randomID}" role="dialog">
            <div class="${dialogClass} modal-dialog modal-dialog-centered"
                 style="pointer-events: all; max-width: 90vw; width: 90vw;">
                <div class="modal-content" style="max-height: 80vh; display: flex; flex-direction: column;">
                    <div class="modal-header" style="flex-shrink: 0;">
                        <h4 class="modal-title">${title}</h4>
                        <button type="button" class="close" data-dismiss="modal"
                                aria-hidden="true">${MODALS.CLOSE_BUTTON_HTML}</button>
                    </div>
                    ${noteSection}
                    <div class="modal-body" style="overflow: auto; flex: 1; min-height: 0; padding: ${noteHTML ? '0 20px 20px 20px' : '20px'};">
                        <div class="table-responsive" style="max-height: none;">
                            ${content}
                        </div>
                    </div>

                    <div class="modal-footer" style="flex-shrink: 0;">
                        <button type="button" class="btn btn-outline-danger"
                                data-dismiss="modal">Close</button>
                    </div>
                </div>
            </div>
        </div>`;

    $('body').append(template);

    const $modal = $(`#modal${randomID}`);
    const escHandlerNamespace = `keydown.modalClose${randomID}`;

    // If a lingering waiting dialog is still around, clear it before/when this modal shows.
    if (typeof hideWaitingDialogSafely === 'function') {
        hideWaitingDialogSafely();
        $modal.on('shown.bs.modal', hideWaitingDialogSafely);
    }

    const escHandler = (event) => {
        if (event.key === 'Escape') {
            $modal.modal('hide');
        }
    };

    $modal
        .modal({ show: true, backdrop: true, keyboard: true })
        .on('click', function(event) {
            // Close when the click is on the backdrop (outside the dialog)
            if ($(event.target).is(this)) {
                $(this).modal('hide');
            }
        })
        .find('.modal-dialog')
        .draggable({ handle: '.modal-header' });

    $(document).on(escHandlerNamespace, escHandler);

    $modal.on('hidden.bs.modal', function() {
        $(document).off(escHandlerNamespace, escHandler);
        $(this).remove();
    });
}

/**
 * Show taxonomy table dialog
 * @param {string} title - Dialog title
 * @param {string} content - Dialog content HTML
 */
function showTaxonomyTableDialog(title, content) {
    _createModalDialog({
        title,
        content,
        modalClass: 'taxonomyTableDialog',
        dialogClass: 'taxonomy-modal-dialog'
    });
}

/**
 * Show gene functions summary table dialog
 * @param {string} title - Dialog title
 * @param {string} content - Dialog content HTML
 */
function showGeneFunctionsSummaryTableDialog(title, content) {
     const noteHTML = `Tables below show the functions encoded by the genes in this bin, and their involvement in
                       metabolic modules, if any...

                       Please note that this is just a quick view of the functions associated with your genes.
                       A much more appropriate way to summarize this information is to use the program
                       <a href="http://merenlab.org/software/anvio/help/programs/anvi-summarize/" target="_blank">
                       anvi-summarize</a> when applicable, and inspect the resulting TAB-delimited output file.`;
    _createModalDialog({
        title,
        content,
        modalClass: 'geneFunctionsSummaryDialog',
        dialogClass: 'gene-functions-modal-dialog',
        noteHTML
    });
}

/**
 * Show gene cluster functions summary table dialog
 * @param {string} title - Dialog title
 * @param {string} content - Dialog content HTML
 */
function showGeneClusterFunctionsSummaryTableDialog(title, content) {
    const noteHTML = `Tables below show the functions associated with the gene clusters in this bin, and their
                      involvement in metabolic modules, if any.

                      Please note that this is just a quick view of the functions associated with your gene
                      clusters. A much more appropriate way to summarize this information is to use the program
                      <a href="http://merenlab.org/software/anvio/help/programs/anvi-summarize/" target="_blank">
                      anvi-summarize</a> when applicable, and inspect the resulting TAB-delimited output file.`;

    _createModalDialog({
        title,
        content,
        modalClass: 'geneClusterFunctionsSummaryDialog',
        dialogClass: 'gene-cluster-functions-modal-dialog',
        noteHTML
    });
}

/**
 * Show gene functions in splits summary table dialog
 * @param {string} title - Dialog title
 * @param {string} content - Dialog content HTML
 */
function showGeneFunctionsInSplitsSummaryTableDialog(title, content) {
    const noteHTML = `Tables below show the functions associated with the genes in splits in this bin, and their
                      involvement in metabolic modules, if any.

                      Please note that this is just a quick view of the functions associated with your genes.
                      A much more appropriate way to summarize this information is to use the program
                      <a href="http://merenlab.org/software/anvio/help/programs/anvi-summarize/" target="_blank">
                      anvi-summarize</a> when applicable, and inspect the resulting TAB-delimited output file.`;

    _createModalDialog({
        title,
        content,
        modalClass: 'geneFunctionsInSplitsSummaryDialog',
        dialogClass: 'gene-function-in-splits-modal-dialog',
        noteHTML
    });
}

/**
 * Show draggable dialog
 * @param {string} title - Dialog title
 * @param {string} content - Dialog content HTML
 * @param {boolean} updateOnly - Only update existing dialog
 */
function showDraggableDialog(title, content, updateOnly) {
    const randomID = title.hashCode();

    if (updateOnly) {
        if (checkObjectExists(`#modal${randomID}`)) {
            $(`#modal${randomID}`).find('.modal-body').html(content);
        }
        return;
    }

    const template = `
        <div class="modal" id="modal${randomID}" data-backdrop="false" style="pointer-events: none;">
            <div class="modal-dialog" style="pointer-events: all;">
                <div class="modal-content no-shadow">
                    <div class="modal-header">
                        <h4 class="modal-title">${title}</h4>
                        <button type="button" class="close" data-dismiss="modal"
                                aria-hidden="true">${MODALS.CLOSE_BUTTON_HTML}</button>
                    </div>
                    <div class="modal-body">${content}</div>
                    <div class="modal-footer">
                        <button type="button" class="btn btn-outline-danger"
                                data-dismiss="modal">Close</button>
                    </div>
                </div>
            </div>
        </div>`;

    $('body').append(template);
    $(`#modal${randomID}`)
        .modal({ show: true, backdrop: false, keyboard: false })
        .find('.modal-dialog')
        .draggable({ handle: '.modal-header' });

    $(`#modal${randomID}`).on('hidden.bs.modal', function() {
        $(this).remove();
    });
}

// ============================================================================
// Utility Functions
// ============================================================================

/**
 * Check if Ctrl key is pressed (handles Mac vs PC)
 * @param {Event} event - Keyboard/mouse event
 * @returns {boolean} True if Ctrl/Cmd is pressed
 */
function IsCtrlPressed(event) {
    const isMac = navigator.platform.toUpperCase().indexOf('MAC') >= 0;
    return (isMac && event.metaKey) || event.ctrlKey;
}

/**
 * Check if jQuery selector finds any elements
 * @param {string} selector - jQuery selector
 * @returns {boolean} True if elements exist
 */
function checkObjectExists(selector) {
    return $(selector).length > 0;
}

/**
 * Generate hash code from string
 * @returns {number} Positive hash code
 */
String.prototype.hashCode = function() {
    let hash = 0;
    if (this.length === 0) return hash;

    for (let i = 0; i < this.length; i++) {
        const chr = this.charCodeAt(i);
        hash = ((hash << 5) - hash) + chr;
        hash |= 0; // Convert to 32bit integer
    }

    return Math.abs(hash);
};

/**
 * Check background process status
 */
function checkBackgroundProcess() {
    const errorMessage = "It seems the server that's been serving this page is no longer accessible. " +
                        "You may lose your unsaved changes in this window.";

    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/session_id',
        success: function(data) {
            if (data != session_id) {
                toastr.error(errorMessage, "", {
                    'timeOut': '0',
                    'extendedTimeOut': '0'
                });
            } else {
                setTimeout(checkBackgroundProcess, 5000);
            }
        },
        error: function() {
            toastr.error(errorMessage, "", {
                'timeOut': '0',
                'extendedTimeOut': '0'
            });
        }
    });
}

/**
 * Check if value is empty
 * @param {*} data - Value to check
 * @returns {boolean} True if empty
 */
function isEmpty(data) {
    if (typeof data === 'object') {
        if (!data || JSON.stringify(data) === '{}' || JSON.stringify(data) === '[]') {
            return true;
        }
        return false;
    } else if (typeof data === 'string') {
        return !data.trim();
    } else if (typeof data === 'undefined') {
        return true;
    }

    return false;
}

/**
 * Check if value is a number
 * @param {*} o - Value to check
 * @returns {boolean} True if number
 */
function isNumber(o) {
    return !isEmpty(o) && !isNaN(o - 0);
}

/**
 * Clear text selection in browser
 */
function clearTextSelection() {
    if (window.getSelection) {
        const selection = window.getSelection();
        if (selection.empty) {
            selection.empty(); // Chrome
        } else if (selection.removeAllRanges) {
            selection.removeAllRanges(); // Firefox
        }
    } else if (document.selection) {
        document.selection.empty(); // IE
    }
}

/**
 * Check if string is alphanumeric
 * @param {string} str - String to check
 * @returns {boolean} True if alphanumeric
 */
function ctype_alnum(str) {
    return /^[a-z0-9]+$/i.test(str);
}

/**
 * Strip HTML tags from string
 * @param {string} html - HTML string
 * @returns {string} Plain text
 */
function strip(html) {
    const tmp = document.createElement("DIV");
    tmp.innerHTML = html;
    return tmp.textContent || tmp.innerText || "";
}

/**
 * Clear min/max inputs for a selectbox row
 * @param {HTMLElement} selectbox - Select element
 */
function clearMinMax(selectbox) {
    const $tr = $(selectbox).parent().parent();
    $tr.find('.input-min').val('0').prop('disabled', true);
    $tr.find('.input-max').val('0').prop('disabled', true);
}

/**
 * Toggle picker visibility based on selection
 * @param {HTMLElement} selectbox - Select element
 * @param {boolean} togglePicker - Whether to toggle picker
 */
function togglePickerStart(selectbox, togglePicker) {
    const $tr = $(selectbox).parent().parent();
    const showPickers = ['intensity', 'line', 'text'].includes(selectbox.value);

    if (showPickers) {
        $tr.find('.picker_start').css('visibility', 'visible');
        if (togglePicker) {
            $tr.find('.picker_end').css('visibility', 'visible');
            $tr.find('.input-height').css('visibility', 'hidden');
            $('.max-font-size-input').show();
        }
    } else {
        $tr.find('.picker_start').css('visibility', 'hidden');
        if (togglePicker) {
            $tr.find('.picker_end').css('visibility', 'hidden');
            $tr.find('.input-height').css('visibility', 'visible').val('30');
        }
    }
}

/**
 * Simple performance timer
 * @param {string} name - Timer name
 */
function BasicTimer(name) {
    this.name = name;
    this.start = new Date().getTime();
    this.previousDelta = this.start;

    /**
     * Get time elapsed and log to console
     * @param {string} event - Event description
     * @param {boolean} consoleOutput - Whether to log to console
     * @returns {Object} Timing information
     */
    this.getDeltaSeconds = function(event, consoleOutput = true) {
        const now = new Date().getTime();
        const deltaSecondsStart = (now - this.start) / 1000;
        const deltaSecondsPrev = (now - this.previousDelta) / 1000;

        this.previousDelta = now;

        const prettyText = `${this.name} [${event}]: ${readableNumber(deltaSecondsPrev)} seconds ` +
                          `(${readableNumber(deltaSecondsStart)} seconds since beginning)`;

        if (consoleOutput) {
            console.log(prettyText);
        }

        return { deltaSecondsStart, deltaSecondsPrev, prettyText };
    };
}

// ============================================================================
// Number Formatting
// ============================================================================

/**
 * Format number in readable format with unit suffixes
 * @param {number} num - Number to format
 * @returns {string} Formatted number
 */
function readableNumber(num) {
    if (num === 0) return 0;
    if (num < 1) return num;

    const e = Math.floor(Math.log(num) / Math.log(1000));
    let result = (num / Math.pow(1000, e)).toPrecision(3);

    if (e < UNITS.READABLE_NUMBER.length) {
        result += UNITS.READABLE_NUMBER[e];
    } else {
        result += 'x10^' + Math.floor(Math.log(num));
    }

    return result;
}

/**
 * Get readable sequence size string
 * @param {number} seqSizeInBases - Sequence size in bases
 * @returns {string|number} Formatted size
 */
function getReadableSeqSizeString(seqSizeInBases) {
    if (seqSizeInBases < 1000) return seqSizeInBases;

    let i = -1;
    let size = seqSizeInBases;

    do {
        size = size / 1000;
        i++;
    } while (size >= 1000 && i < UNITS.SEQ_SIZE.length - 1);

    return Math.round(size) + UNITS.SEQ_SIZE[i];
}

/**
 * Format number with thousands separator
 * @param {number} number - Number to format
 * @param {number} decimals - Decimal places
 * @param {string} dec_point - Decimal separator
 * @param {string} thousands_sep - Thousands separator
 * @returns {string} Formatted number
 */
function getCommafiedNumberString(number, decimals = 2, dec_point = '.', thousands_sep = ',') {
    if (isNaN(parseInt(number))) return number;

    const n = !isFinite(+number) ? 0 : +number;
    const prec = !isFinite(+decimals) ? 2 : Math.abs(decimals);

    const toFixedFix = function(n, prec) {
        const k = Math.pow(10, prec);
        return Math.round(n * k) / k;
    };

    const s = (prec ? toFixedFix(n, prec) : Math.round(n)).toString().split('.');

    if (s[0].length > 3) {
        s[0] = s[0].replace(/\B(?=(?:\d{3})+(?!\d))/g, thousands_sep);
    }

    if ((s[1] || '').length < prec) {
        s[1] = s[1] || '';
        s[1] += new Array(prec - s[1].length + 1).join('0');
    }

    return s[1] > 0 ? s.join(dec_point) : s[0];
}

// ============================================================================
// Geometry and SVG Utilities
// ============================================================================

/**
 * Create SVG line path between two points
 * @param {Object} p0 - Start point {x, y}
 * @param {Object} p1 - End point {x, y}
 * @returns {string} SVG path string
 */
function linePath(p0, p1) {
    return `M ${p0.x} ${p0.y} ${p1.x} ${p1.y}`;
}

/**
 * Calculate distance between two points
 * @param {Object} p0 - First point {x, y}
 * @param {Object} p1 - Second point {x, y}
 * @returns {number} Distance
 */
function distance(p0, p1) {
    const dx = p1.x - p0.x;
    const dy = p1.y - p0.y;
    return Math.sqrt(dx * dx + dy * dy);
}

// ============================================================================
// String Polyfills and Math Extensions
// ============================================================================

// Add trim method if not available (old browsers)
if (!String.prototype.trim) {
    String.prototype.trim = function() {
        return this.replace(/^\s+|\s+$/g, '');
    };
}

/**
 * Convert degrees to radians
 * @param {number} degrees - Angle in degrees
 * @returns {number} Angle in radians
 */
Math.toRadians = function(degrees) {
    return degrees * Math.PI / 180;
};

/**
 * Convert radians to degrees
 * @param {number} radians - Angle in radians
 * @returns {number} Angle in degrees
 */
Math.toDegrees = function(radians) {
    return radians * 180 / Math.PI;
};

// ============================================================================
// Layer Data Operations
// ============================================================================

/**
 * Remove single parents from layer data
 */
function removeSingleParents() {
    // Note: layerdata and parameter_count are global
    for (let i = 1; i < parameter_count; i++) {
        if (layerdata[0][i] !== '__parent__') continue;

        const parent_count_dict = {};

        // Count occurrences of each parent
        for (let j = 1; j < layerdata.length; j++) {
            if (!layerdata[j][i]) continue;

            const parent = layerdata[j][i];
            parent_count_dict[parent] = (parent_count_dict[parent] || 0) + 1;
        }

        // Remove single parents
        for (const [parent_name, count] of Object.entries(parent_count_dict)) {
            if (count === 1) {
                for (let j = 1; j < layerdata.length; j++) {
                    if (layerdata[j][i] === parent_name) {
                        layerdata[j][i] = '';
                    }
                }
            }
        }
    }
}

// ============================================================================
// jQuery Table Sort Helper
// ============================================================================

/**
 * Helper for sortable table rows
 * @param {Event} e - Event
 * @param {jQuery} tr - Table row
 * @returns {jQuery} Helper element
 */
const fixHelperModified = function(e, tr) {
    const $originals = tr.children();
    const $helper = tr.clone();

    $helper.children().each(function(index) {
        $(this).width($originals.eq(index).width());
    });

    return $helper;
};

// ============================================================================
// Zoom and Transform Operations
// ============================================================================

/**
 * Get current transformation matrix
 * @returns {number[]} Matrix values
 */
function getMatrix() {
    const viewport = document.getElementById('viewport');
    const transform = viewport.getAttribute('transform');
    return transform
        .split('(')[1]
        .split(')')[0]
        .split(',')
        .map(parseFloat);
}

/**
 * Set transformation matrix
 * @param {number[]} matrix - Matrix values
 */
function setMatrix(matrix) {
    const viewport = document.getElementById('viewport');
    viewport.setAttribute('transform', `matrix(${matrix.join(',')})`);
}

/**
 * Apply zoom to viewport
 * @param {number} scale - Zoom scale factor
 */
function zoom(scale) {
    const matrix = getMatrix();

    // Scale all matrix values
    for (let i = 0; i < matrix.length; i++) {
        matrix[i] *= scale;
    }

    // Adjust position to zoom from center
    matrix[4] += (1 - scale) * VIEWER_WIDTH / 2;
    matrix[5] += (1 - scale) * VIEWER_HEIGHT / 2;

    setMatrix(matrix);
}

/**
 * Pan viewport
 * @param {number} dx - X offset
 * @param {number} dy - Y offset
 */
function pan(dx, dy) {
    const matrix = getMatrix();
    matrix[4] += dx;
    matrix[5] += dy;
    setMatrix(matrix);
}

/**
 * Reset zoom to fit content
 */
function zoom_reset() {
    const viewport = document.getElementById('viewport');
    const bbox = viewport.getBBox();

    const scale = Math.min(VIEWER_WIDTH / bbox.width, VIEWER_HEIGHT / bbox.height) * 0.80;
    const centerX = VIEWER_WIDTH / 2 - (bbox.x + bbox.width / 2) * scale;
    const centerY = VIEWER_HEIGHT / 2 - (bbox.y + bbox.height / 2) * scale;

    const baseMatrix = [scale, 0, 0, scale, centerX, centerY];
    setMatrix(baseMatrix);
}


// ─── Shared functions/metabolism display builders ────────────────────────────
// Used by both anvi-display-pan (main.js) and anvi-display-pan-graph (pangenome-graph.js).
// Kept here in utils.js so both interfaces can use them without loading main.js.

function copyTextWithFeedback(btn, text) {
    if (typeof text !== 'string' || !text.length) return;
    navigator.clipboard.writeText(text)
        .then(function () {
            const icon = btn.nextElementSibling;
            if (!icon) return;
            icon.textContent = "✔";
            icon.style.color = "green";
            icon.style.marginLeft = "4px";
            setTimeout(function () { icon.textContent = ""; }, 2000);
        })
        .catch(function (err) { console.error("Copy failed", err); });
}

function buildCopyableNamesSection(label, items, options = {}) {
    const safeItems = Array.isArray(items) ? items : [];
    const headerText = `${label} (${safeItems.length})`;
    const background = options.background || '#ffe4c478';
    const marginBottom = options.marginBottom || '20px';
    const copySeparator = options.copySeparator || '\n';
    const contentRenderer = options.contentRenderer;
    const copyReadyNames = safeItems.join(copySeparator);
    const header = `
        <div class="bin-modal-header" style="background: ${background};">
            <span style="font-size: large;">${headerText}</span>
            <button class="btn btn-primary btn-sm"
                onclick='copyTextWithFeedback(this, ${JSON.stringify(copyReadyNames)})'>Copy to Clipboard</button>
            <span class="copy-feedback" style="min-width: 14px; font-size: 0.9em; line-height: 1;"></span>
        </div>`;
    const body = typeof contentRenderer === 'function'
        ? contentRenderer(safeItems)
        : `<p style="margin-bottom: ${marginBottom};">${safeItems.join(', ') || 'No items found.'}</p>`;
    return header + body;
}

function buildCopyButton(label, items) {
    const safeItems = Array.isArray(items) ? items : [];
    const copyText = safeItems.join('\n');
    return `
        <button class="btn btn-primary btn-sm"
            style="margin-left: 12px;"
            title="Copy ${label}"
            onclick='copyTextWithFeedback(this, ${JSON.stringify(copyText)})'>
            Copy ${label}
        </button>
        <span class="copy-feedback" style="min-width: 14px; font-size: 0.9em; line-height: 1;"></span>`;
}

function buildFunctionsContent(response, config) {
    const fmtPct = (v) => (typeof v === 'number' && !isNaN(v)) ? (v * 100).toFixed(1) + '%' : 'NA';
    let content = '';
    if (config.dialogFunction === 'showGeneFunctionsInSplitsSummaryTableDialog') {
        content += buildContigAndSplitNamesTable(response, config);
    }
    content += buildMetabolismTable(response, config, fmtPct);
    content += '<hr style="margin: 30px !important;">';
    content += buildFunctionsTable(response, config);
    return content;
}

function buildContigAndSplitNamesTable(response, config) {
    const contigNames = response && response.contig_names;
    const splitNames  = response && response.split_names;
    let out = '';
    out += buildCopyableNamesSection('Contig Names', contigNames, {
        background: '#f5f5dc9c', marginBottom: '35px', copySeparator: '\n',
        contentRenderer: (names) => `<p style="margin-bottom: 35px;">${(names || []).join(', ')}</p>`
    });
    out += buildCopyableNamesSection('Split Names', splitNames, {
        background: '#f5f5dc9c', marginBottom: '35px', copySeparator: '\n',
        contentRenderer: (names) => `<p style="margin-bottom: 35px;">${(names || []).join(', ')}</p>`
    });
    return out;
}

function buildMetabolismTable(response, config, fmtPct) {
    const metabolism = response && response.metabolism;
    let out = `<p class="bin-modal-header" style="background: #f5f5dc9c">Metabolic module involvement</p>`;
    if (!metabolism || typeof metabolism !== 'object' || !Object.keys(metabolism).length) {
        out += '<p style="margin-bottom: 35px;">There are no metabolic insights to show here :/</p>';
        return out;
    }
    out += `<p style="margin-bottom: 35px;">${config.metabolismDescription}</p>
        <table class="table table-sm table-striped" style="width: 95%; margin-left: 10px;">
            <thead class="thead-light"><tr>
                <th>Metabolic module</th>
                <th style="text-align: center;">Contribution to pathway completeness</th>
                <th style="text-align: center;">Contribution to Stepwise completeness</th>
                <th style="text-align: center;">${config.itemLabel === 'gene clusters' ? 'Gene clusters' : 'Gene calls'} involved</th>
            </tr></thead><tbody>`;
    Object.keys(metabolism)
        .sort((a, b) => {
            const pa = metabolism[a]?.pathwise_percent_complete ?? 0;
            const pb = metabolism[b]?.pathwise_percent_complete ?? 0;
            if (pb !== pa) return pb - pa;
            const ga = Array.isArray(metabolism[a]?.gene_caller_ids) ? metabolism[a].gene_caller_ids.length : 0;
            const gb = Array.isArray(metabolism[b]?.gene_caller_ids) ? metabolism[b].gene_caller_ids.length : 0;
            return gb - ga;
        })
        .forEach((moduleId) => {
            const m = metabolism[moduleId] || {};
            const genes = Array.isArray(m.gene_caller_ids)
                ? [...m.gene_caller_ids].map(String)
                    .sort((a, b) => a.localeCompare(b, undefined, { numeric: true, sensitivity: "base" }))
                    .join(", ")
                : "NA";
            const pathwayPct  = fmtPct(m.pathwise_percent_complete);
            const stepwisePct = fmtPct(m.stepwise_completeness);
            const moduleName  = m.NAME ?? 'NA';
            const moduleClass = (typeof m.CLASS === 'string' ? m.CLASS : 'NA');
            const complete    = (m.pathwise_is_complete === true) || (m.stepwise_is_complete === true);
            const badge       = complete ? ` <span class="badge badge-success">complete</span>` : '';
            out += `<tr>
                <td style="padding-left: 20px;">
                   <b>${moduleId}</b>${badge}<br />
                   - Module function: ${moduleName}<br />
                   - Module class: ${moduleClass.split(";").map((p,i,a) => i===a.length-1 ? `<b>${p.trim()}</b>` : p).join("; ")}
                </td>
                <td style="text-align: center; vertical-align: middle;">${pathwayPct}</td>
                <td style="text-align: center; vertical-align: middle;">${stepwisePct}</td>
                <td style="text-align: center; vertical-align: middle; max-width: 420px; word-break: break-word;">${genes}</td>
            </tr>`;
        });
    out += `</tbody></table>`;
    return out;
}

function buildFilterControls(sources, functions, config, gene_clusters) {
    const totalItems = Object.keys(functions).length;
    const sourceCounts = {};
    Object.keys(sources).forEach(function(index) {
        let source = sources[index];
        sourceCounts[source] = 0;
        Object.keys(functions).forEach(function(item_id) {
            let d = functions[item_id] || {};
            if (d[source]) {
                const accession_string = config.getAccessionString(d, source);
                const function_string  = config.getFunctionString(d, source);
                if (accession_string !== 'N/A' || function_string !== 'N/A') {
                    sourceCounts[source]++;
                }
            }
        });
    });
    let out = `
        <div style="background-color: #f8f9fa; padding: 15px; margin: 20px 10px; border-radius: 5px;">
            <div style="margin-bottom: 10px;">
                <label><input type="checkbox" id="hideNA"> Hide entries with no annotation to simplify the display</label>
            </div>
            <hr>
            <div>
                <p>Annotation sources to display for a total of ${totalItems} ${config.itemLabel}:</p>
                <div style="margin-top: 5px; display: flex; flex-wrap: wrap; gap: 15px;">`;
    Object.keys(sources).forEach(function(index) {
        let source = sources[index];
        out += `<label style="margin-right: 15px;">
                    <input type="checkbox" class="source-filter" data-source="${source}" checked>
                    ${source} <span style="color: #666; font-size: 0.9em;">(${sourceCounts[source]})</span>
                </label>`;
    });
    out += `</div>
                <div style="margin-top: 10px;">
                    <button type="button" class="btn btn-sm btn-outline-secondary" id="selectAllSources">Select All</button>
                    <button type="button" class="btn btn-sm btn-outline-secondary" id="deselectAllSources">Deselect All</button>
                    <button type="button" class="btn btn-sm btn-primary" id="copyTableBelow" style="margin-left: 10px;">Copy table below</button><span class="copy-feedback" style="min-width:14px;font-size:0.9em;line-height:1;margin-right:6px;"></span>
                    <button type="button" class="btn btn-sm btn-primary" id="copyComprehensiveFunctions">Copy comprehensive function data</button><span class="copy-feedback" style="min-width:14px;font-size:0.9em;line-height:1;"></span>
                </div>
            </div>
        </div>`;
    return out;
}

function buildFunctionsTable(response, config) {
    if (config.dialogFunction === 'showGeneFunctionsInSplitsSummaryTableDialog') return "";
    const itemIds = Object.keys(response['functions'] || {});
    const copyButton = (config.itemLabel === 'gene clusters' || config.itemLabel === 'synteny gene clusters')
        ? buildCopyButton(`${config.itemLabel} names`, itemIds)
        : '';
    let content = `
        <div class="bin-modal-header" style="background: #ffe4c478;">
            <span style="font-size: large;">Functions per ${config.itemLabel}</span>
            ${copyButton}
        </div>
        <p>${config.functionsDescription}</p>`;
    content += buildFilterControls(response['sources'], response['functions'], config, response['gene_clusters']);
    content += `
        <table class="table" id="itemFunctionsTable" style="width: 95%; margin-left: 10px; table-layout: fixed;">
           <thead class="thead-light"><tr>
             <th style="width: 10%;">${config.itemIdLabel}</th>
             <th style="width: 10%;">Source</th>
             <th style="width: 10%;">Accession</th>
             <th style="width: 70%;">Function</th>
           </tr></thead><tbody>`;
    Object.keys(response['functions']).forEach(function(item_id) {
        let d = response['functions'][item_id] || {};
        Object.keys(response['sources']).forEach(function(index) {
            let function_source = response['sources'][index];
            let accession_string, function_string;
            let hasAnnotation = false;
            if (d[function_source]) {
                accession_string = config.getAccessionString(d, function_source);
                function_string  = config.getFunctionString(d, function_source);
                hasAnnotation    = (accession_string !== 'N/A' && function_string !== 'N/A');
            } else {
                accession_string = 'N/A';
                function_string  = 'N/A';
            }
            const dataAttrs = `data-source="${function_source}" data-has-annotation="${hasAnnotation}" data-item="${item_id}"`;
            if (index == 0) {
                content += `<tr style="border-top: 3px solid #d0d0d0;" ${dataAttrs}>
                    <td rowspan="${Object.keys(response['sources']).length}" style="vertical-align: middle; word-wrap: break-word;"><b>${item_id}</b></td>
                    <td style="word-wrap: break-word;">${function_source}</td>
                    <td style="word-wrap: break-word;">${accession_string}</td>
                    <td style="word-wrap: break-word;">${function_string}</td>
                    </tr>`;
            } else {
                content += `<tr ${dataAttrs}>
                    <td style="word-wrap: break-word;">${function_source}</td>
                    <td style="word-wrap: break-word;">${accession_string}</td>
                    <td style="word-wrap: break-word;">${function_string}</td>
                    </tr>`;
            }
        });
    });
    content += `</tbody></table>`;
    return content;
}

function setupItemTableFiltering(gene_clusters) {
    const table = document.getElementById('itemFunctionsTable');
    if (!table) return;
    const hideNACheckbox   = document.getElementById('hideNA');
    const sourceCheckboxes = document.querySelectorAll('.source-filter');
    const selectAllBtn     = document.getElementById('selectAllSources');
    const deselectAllBtn   = document.getElementById('deselectAllSources');
    const copyTableBtn     = document.getElementById('copyTableBelow');
    const copyCompBtn      = document.getElementById('copyComprehensiveFunctions');

    function applyFilters() {
        const hideNA = hideNACheckbox.checked;
        const activeSources = Array.from(sourceCheckboxes).filter(cb => cb.checked).map(cb => cb.dataset.source);
        table.querySelectorAll('tbody tr').forEach(row => {
            const hasAnnotation = row.dataset.hasAnnotation === 'true';
            const source = row.dataset.source;
            let shouldShow = true;
            if (hideNA && !hasAnnotation) shouldShow = false;
            if (!activeSources.includes(source)) shouldShow = false;
            row.style.display = shouldShow ? '' : 'none';
        });
        updateItemRowspans();
        updateTableStripes();
    }

    function updateTableStripes() {
        const visibleRows = Array.from(table.querySelectorAll('tbody tr')).filter(r => r.style.display !== 'none');
        const geneGroups = {};
        visibleRows.forEach(row => {
            const gene = row.dataset.item;
            if (!geneGroups[gene]) geneGroups[gene] = [];
            geneGroups[gene].push(row);
        });
        visibleRows.forEach((row, index) => {
            const cells = row.querySelectorAll('td');
            cells.forEach(cell => { cell.style.backgroundColor = ''; cell.classList.remove('table-striped-manual'); });
            row.style.borderTop = '';
            if (index % 2 === 1) {
                cells.forEach((cell, cellIndex) => {
                    if (cellIndex === 0 && (cell.getAttribute('rowspan') || cell.classList.contains('dynamic-gene-cell'))) return;
                    cell.style.backgroundColor = 'rgba(0,0,0,.05)';
                    cell.classList.add('table-striped-manual');
                });
            }
        });
        Object.keys(geneGroups).forEach(gene => {
            if (geneGroups[gene].length > 0) geneGroups[gene][0].style.borderTop = '3px solid #d0d0d0';
        });
    }

    function updateItemRowspans() {
        const allRows = Array.from(table.querySelectorAll('tbody tr'));
        const itemGroups = {};
        allRows.forEach(row => {
            const item = row.dataset.item;
            if (!itemGroups[item]) itemGroups[item] = [];
            itemGroups[item].push(row);
        });
        Object.keys(itemGroups).forEach(item => {
            const allRowsForItem     = itemGroups[item];
            const visibleRowsForItem = allRowsForItem.filter(r => r.style.display !== 'none');
            allRowsForItem.forEach(row => {
                const existingCells = row.querySelectorAll('td');
                if (existingCells.length > 3) {
                    existingCells[0].style.display = 'none';
                    existingCells[0].removeAttribute('rowspan');
                }
                const dynamicCell = row.querySelector('td.dynamic-gene-cell');
                if (dynamicCell) dynamicCell.remove();
            });
            if (visibleRowsForItem.length > 0) {
                const firstRow = visibleRowsForItem[0];
                let geneCell = firstRow.querySelector('td:first-child');
                if (geneCell && firstRow.querySelectorAll('td').length > 3) {
                    geneCell.style.display = '';
                    geneCell.style.verticalAlign = 'middle';
                    geneCell.style.wordWrap = 'break-word';
                    geneCell.setAttribute('rowspan', visibleRowsForItem.length);
                } else {
                    const newCell = document.createElement('td');
                    newCell.innerHTML = `<b>${item}</b>`;
                    newCell.setAttribute('rowspan', visibleRowsForItem.length);
                    newCell.className = 'dynamic-gene-cell';
                    newCell.style.borderTop = '3px solid #d0d0d0';
                    newCell.style.verticalAlign = 'middle';
                    newCell.style.wordWrap = 'break-word';
                    firstRow.insertBefore(newCell, firstRow.firstChild);
                }
            }
        });
    }

    hideNACheckbox.addEventListener('change', applyFilters);
    sourceCheckboxes.forEach(cb => cb.addEventListener('change', applyFilters));
    selectAllBtn.addEventListener('click', () => { sourceCheckboxes.forEach(cb => cb.checked = true); applyFilters(); });
    deselectAllBtn.addEventListener('click', () => { sourceCheckboxes.forEach(cb => cb.checked = false); applyFilters(); });

    if (copyTableBtn) {
        copyTableBtn.addEventListener('click', function() {
            const rows = Array.from(table.querySelectorAll('tbody tr')).filter(r => r.style.display !== 'none');
            let tsv = 'SynGC\tSource\tAccession\tFunction\n';
            rows.forEach(row => {
                const item = row.dataset.item;
                // get only visible cells that are not the item/rowspan cell
                const dataCells = Array.from(row.querySelectorAll('td')).filter(td =>
                    td.style.display !== 'none' && !td.classList.contains('dynamic-gene-cell') &&
                    !td.hasAttribute('rowspan')
                );
                if (dataCells.length >= 3) {
                    tsv += `${item}\t${dataCells[0].textContent.trim()}\t${dataCells[1].textContent.trim()}\t${dataCells[2].textContent.trim()}\n`;
                }
            });
            copyTextWithFeedback(copyTableBtn, tsv);
        });
    }

    if (copyCompBtn && !gene_clusters) {
        copyCompBtn.addEventListener('click', function() {
            if (typeof toastr !== 'undefined') toastr.warning('Per-gene data not available. Please restart the anvi\'o server.', 'No data');
        });
    }
    if (copyCompBtn && gene_clusters) {
        copyCompBtn.addEventListener('click', function() {
            const activeSources = Array.from(sourceCheckboxes).filter(cb => cb.checked).map(cb => cb.dataset.source);
            const hideNA = hideNACheckbox.checked;
            let tsv = 'SynGC\tGenome_Name\tGene_Callers_ID\tSource\tAccession\tFunction\n';
            Object.keys(gene_clusters).forEach(sgc => {
                const genomes = gene_clusters[sgc] || {};
                Object.keys(genomes).forEach(genome => {
                    const genes = genomes[genome] || {};
                    Object.keys(genes).forEach(gene_callers_id => {
                        const annots = genes[gene_callers_id] || {};
                        activeSources.forEach(source => {
                            const blob = annots[source] || '';
                            const parts = blob.split('|||');
                            const accession = parts[0] || 'N/A';
                            const func      = parts[1] || 'N/A';
                            if (hideNA && accession === 'N/A' && func === 'N/A') return;
                            tsv += `${sgc}\t${genome}\t${gene_callers_id}\t${source}\t${accession}\t${func}\n`;
                        });
                    });
                });
            });
            copyTextWithFeedback(copyCompBtn, tsv);
        });
    }

    applyFilters();
}
