/**
 * Javascript library to visualize gene structures
 *
 *  Authors: Evan Kiefl <kiefl.evan@gmail.com>
 *           Ozcan Esen
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

const mode = 'structure';
const MAX_NGL_WIDGETS = 16;

var stages = {};
var variability = {};
var histogram_data;
var residue_info;
var residue_info_types;
var column_info;
var sample_groups;
var pdb_content;
var cached_orientation_matrices = {};

var color_legend = {};
var size_legend = {};

var filter_backup = {};
var sample_groups_backup = {};

var current_state_name;


$(document).ready(function() {
    toastr.options = {
        "closeButton": true,
        "debug": false,
        "newestOnTop": true,
        "progressBar": false,
        "positionClass": "toast-top-right",
        "preventDuplicates": false,
        "onclick": null,
        "showDuration": "500",
        "hideDuration": "2000",
        "timeOut": "6000",
        "extendedTimeOut": "1000",
        "showEasing": "swing",
        "hideEasing": "linear",
        "showMethod": "fadeIn",
        "hideMethod": "fadeOut",
    }

    $('.colorpicker').colpick({
        layout: 'hex',
        submit: 0,
        colorScheme: 'light',
        onChange: function(hsb, hex, rgb, el, bySetColor) {
            $(el).css('background-color', '#' + hex);
            $(el).attr('color', '#' + hex);

            if (!bySetColor) $(el).val(hex);
        },
        onHide: function(cal) {
            let el = $(cal).data('colpick').el;
            $(el).trigger('change');
        }
    }).keyup(function() {
        $(this).colpickSetColor(this.value);
    });

    window.addEventListener( "resize", function( event ){
        for (let group in stages) {
            stages[group].handleResize();
        }
    }, false );

    $('#gene_callers_id_list').on('change', function(ev) {
        $.when({}).then(load_protein).then(load_gene_function_info).then(load_model_info).then(() => {
            create_ui();
            load_sample_group_widget($('#sample_groups_list').val());
        });
    });

    $('#sample_groups_list').on('change', function(ev) {
        backupGroupsWidget();
        $('.overlay').show();
        load_sample_group_widget($('#sample_groups_list').val());
    });

    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/get_initial_data?timestamp=' + new Date().getTime(),
        success: function(data) {
            let available_gene_callers_ids = data['available_gene_callers_ids'];
            let available_engines = data['available_engines'];
            sample_groups = data['sample_groups'];

            available_gene_callers_ids.forEach(function(gene_callers_id) {
                $('#gene_callers_id_list').append(`<option id=${gene_callers_id}>${gene_callers_id}</option>`);
            });

            $.when({}).then(load_protein).then(load_gene_function_info).then(load_model_info).then(() => {
                let default_engine = available_engines[0];
                available_engines.forEach(function(engine) {
                    $('#engine_list').append(`<input type="radio" name="engine" onclick="$.when({}).then(create_ui).then(() => { fetch_and_draw_variability(); });" value="${engine}" id="engine_${engine}" ${engine == default_engine ? 'checked="checked"' : ''}><label for="engine_${engine}">${engine}</label>`);
                });
                create_ui();

                for (let category in sample_groups) {
                    $('#sample_groups_list').append(`<option id=${category}>${category}</option>`);
                }
                $('#sample_groups_list').trigger('change');
            });
        }
    });
});


function load_sample_group_widget(category, trigger_create_ngl_views=true) {
    $('#sample_groups').empty();
    $('#sample_groups').attr('created-for-category', category);

    tableHtml = '<table class="table table-sm table-responsive"><tr><td><label class="col-md-4 settings-label">Groups</label></td><td><label class="col-md-4 settings-label">Samples</label></td></tr>';

    let counter=0;
    for (let group in sample_groups[category]) {
        let group_checked = true;

        if (sample_groups_backup.hasOwnProperty(category)) {
            group_checked = (sample_groups_backup[category]['groups'].indexOf(group) > -1);
        }

        if (group_checked) {
            counter++;
        }

        if (counter > 15) {
            group_checked = false;
        }

        tableHtml += `
            <tr>
                <td>
                    <input class="form-check-input"
                        checkbox-for="group"
                        id="${category}_${group}"
                        type="checkbox"
                        data-category="${category}"
                        data-group="${group}"
                        value="${group}"
                        ${ group_checked ? `checked="checked"` : `` }>
                    <label class="form-check-label" for="${category}_${group}">${group}</label>
                </td>
                <td>`;

        sample_groups[category][group].forEach((sample) => {
            let sample_checked = true;
            if (sample_groups_backup.hasOwnProperty(category) && sample_groups_backup[category]['samples'].hasOwnProperty(group)) {
                sample_checked = (sample_groups_backup[category]['samples'][group].indexOf(sample) > -1);
            }

            tableHtml += `
                <div class="table-group-checkbox" style="display: inline-block; float: left;">
                    <input class="form-check-input"
                            id="${category}_${group}_${sample}"
                            type="checkbox"
                            data-category="${category}"
                            data-group="${group}"
                            data-sample="${sample}"
                            value="${sample}"
                            ${sample_checked ? 'checked="checked"' : ''}>
                    <label class="form-check-label" for="${category}_${group}_${sample}">${sample}</label>
                </div>`;
        });

        tableHtml += '</td></tr>';
    }

    $('#sample_groups').append(tableHtml + '</table>');

    if (trigger_create_ngl_views) {
        create_ngl_views();
    }
}

function apply_orientation_matrix_to_all_stages(orientationMatrix) {
    for (let group in stages) {
        stages[group].viewerControls.orient(orientationMatrix);
    }
    cached_orientation_matrices[$('#gene_callers_id_list').val()] = orientationMatrix;
}

async function create_ngl_views(fetch_variability = true) {
    // if fetch_variability, fetch_and_draw_variability is called, otherwise draw_variability is called
    let selected_groups = $('[checkbox-for="group"]:checked');
    if (selected_groups.length > MAX_NGL_WIDGETS) {
        $('#maximum_ngl_widgets_error').show();
        return;
    } else {
        $('#maximum_ngl_widgets_error').hide();
    }

    for (let group in stages) {
        stages[group].dispose();
    }
    stages = {};

    $('#ngl-container').empty();

    let num_cells = selected_groups.length;
    let num_columns = Math.min(4, Math.ceil(Math.sqrt(num_cells)));
    let num_rows = Math.min(4, Math.ceil(num_cells / num_columns));

    for (let i=0; i < selected_groups.length; i++) {
        let group = selected_groups[i].getAttribute('data-group');
        await create_single_ngl_view(group, num_rows, num_columns);
    }

    if (fetch_variability) {
        fetch_and_draw_variability();
    } else {
        draw_variability();
    }
}

async function create_single_ngl_view(group, num_rows, num_columns) {
    var defer = $.Deferred();

    $('#ngl-container').append(`
        <div id="ngl_${group}_wrapper"
             class="col-md-${parseInt(12 / num_columns)} nopadding"
             style="height: ${parseFloat(100 / num_rows)}%; ">
             <div class="ngl-group-title">
                ${group}
             </div>
             <div class="ngl-group-fullscreen">
                <button type="button" class="btn btn-link btn-sm" onclick="stages['${group}'].toggleFullscreen();" title="Fullscreen">
                    <span class="bi bi-fullscreen"></span>
                 </button>
             </div>
             <div id="ngl_${group}" class="ngl-inner">

             </div>
        </div>`);

    var stage = new NGL.Stage(`ngl_${group}`);
    stage.setSize('100%', '1000px');
    var stringBlob = new Blob( [ pdb_content ], { type: 'text/plain'} );

    stage.setParameters({
        backgroundColor: "white"
    });

    stage.loadFile(stringBlob, { ext: "pdb" }).then((component) => {
        if ( component.type !== "structure" ) return;

        if ($('#show_surface').is(':checked')) {
            if ($('#surface_color_type').val() == 'Static') {
                var color_value = $('#color_static_surface').attr('color');
            } else if ($('#surface_color_type').val() == 'Dynamic') {
                // Show range error if min is greater than max
                if (parseFloat($('#surface_color_min').val()) >= parseFloat($('#surface_color_max').val())) {
                    $('#dynamic_surface_color_error').show();
                } else {
                    $('#dynamic_surface_color_error').hide();
                }
                var color_value = getSurfaceColorScheme(group);
            } else {
                var color_value = $('#surface_color_type').val()
            }

            surface_rep_params = {
                surfaceType: "av",
                smooth: 3,
                opaqueBack: false,
                side: 'double',
                probeRadius: parseFloat($('#surface_probe_radius').val()),
                scaleFactor: 3.0,
                opacity: parseFloat($('#surface_opacity').val()),
                lowResolution: false,
                color: color_value
            }
            component.addRepresentation("surface", surface_rep_params);
        }

        if ($('#show_residue_labels').is(':checked')) {
            component.addRepresentation("label", {
            sele: ".CA",
            color: "element",
            labelType: "format",
            labelFormat: "%(resname)s"
            });
        }

        if ($('#show_backbone').is(':checked')) {

            if ($('#backbone_color_type').val() == 'Static') {
                var color_value = $('#color_static_backbone').attr('color');
            } else if ($('#backbone_color_type').val() == 'Dynamic') {
                // Show range error if min is greater than max
                if (parseFloat($('#backbone_color_min').val()) >= parseFloat($('#backbone_color_max').val())) {
                    $('#dynamic_backbone_color_error').show();
                } else {
                    $('#dynamic_backbone_color_error').hide();
                }
                var color_value = getBackboneColorScheme(group);
            } else {
                var color_value = $('#backbone_color_type').val()
            }

            backbone_type = $('#backbone_type').val()

            if (backbone_type == 'rocket+loop') {
                // custom
                component.addRepresentation('rocket', {
                    color: color_value,
                    scale: 1.8,
                });
                component.addRepresentation('tube', {
                    color: color_value,
                    sele: 'not helix',
                    scale: 2.0,
                    aspectRatio: 1.0,
                });
            } else if (backbone_type == 'custom') {
                // custom
                null
            } else {
                component.addRepresentation(backbone_type, {
                    color: color_value,
                    scale: 1.5,
                    aspectRatio: 3.0,
                });
            }
        }

        if ($('#show_ballstick').is(':checked')) {
            if ($('#show_ballstick_when').val() == 'always') {
                component.addRepresentation("ball+stick", {
                    sele: "sidechainAttached"
                });
            }
        }

        if (cached_orientation_matrices.hasOwnProperty($('#gene_callers_id_list').val())) {
            stage.viewerControls.orient(cached_orientation_matrices[$('#gene_callers_id_list').val()]);
        } else {
            component.autoView();
        }

        // prevent default tooltip
        stage.mouseControls.remove("hoverPick");

        // add custom tooltip
        var previous_hovered_residue = null;
        stage.signals.hovered.add(function (pickingProxy) {
            let tooltip = document.getElementById('ngl-tooltip');

            if (pickingProxy && pickingProxy.atom) {
                if (pickingProxy.atom.resno != previous_hovered_residue && !(['always', 'variant residue'].includes($('#show_ballstick_when').val()))) {
                    // remove ball+stick if hovered residue changed or
                    if (pickingProxy.atom.resno != previous_hovered_residue) {
                        stage.compList[0].reprList.slice(0).forEach((rep) => {
                            if (rep.name == 'ball+stick') {
                                rep.dispose();
                            }
                        });
                    }
                }

                let residue = pickingProxy.atom.resno;
                let mp = pickingProxy.mouse.position;

                if ($('#show_ballstick').is(':checked') && !(['always', 'variant residue'].includes($('#show_ballstick_when').val()))) {
                    if ($('#show_ballstick_when').val() == 'hovered residue') {
                        var selection = "(" + residue + ")" + " and sidechainAttached";
                    }
                    else if ($('#show_ballstick_when').val() == 'hovered residue + contacts') {
                        var selection = "(" + residue_info[residue]['contact_numbers'].split(',').join(', ') + ")" + " and sidechainAttached";
                    }
                    if ($('#show_ballstick_when').val() == 'hovered residue + variant contacts') {
                        let contacts = residue_info[residue]['contact_numbers'].split(',');
                        let variant_contacts = [];
                        for (i in contacts) {
                            if (contacts[i] == String(residue)) {
                                variant_contacts.push(contacts[i]);
                            }
                            else if (variability[group].hasOwnProperty(parseInt(contacts[i]))) {
                                variant_contacts.push(contacts[i]);
                            }
                        }
                        var selection = "(" + variant_contacts.join(', ') + ")" + " and sidechainAttached";
                    }
                    stage.compList[0].addRepresentation("ball+stick", {
                        hydrogenBond: true,
                        sele: selection
                    });
                }

                // Reference data's availability depends on how anvi-gen-structure-database was
                // created. For example, if --skip-DSSP, there is no secondary structure
                var tooltip_HTML_title = `<h5>Reference info</h5>`
                var tooltip_HTML_body = `
                    <tr><td>Residue</td><td>${residue_info[residue]['amino_acid']} (${residue_info[residue]['codon']})</td></tr>
                    <tr><td>Residue No.</td><td>${residue_info[residue]['codon_number']}</td></tr>
                `
                if (residue_info[residue].hasOwnProperty('sec_struct')) {tooltip_HTML_body += `<tr><td>Secondary Structure</td><td class="d-block">${residue_info[residue]['sec_struct']}</td></tr>`}
                if (residue_info[residue].hasOwnProperty('rel_solvent_acc')) {tooltip_HTML_body += `<tr><td>Solvent Accessibility</td><td>${residue_info[residue]['rel_solvent_acc'].toFixed(2)}</td></tr>`}
                if (residue_info[residue].hasOwnProperty('phi')) {tooltip_HTML_body += `<tr><td>(Phi, Psi)</td><td>(${residue_info[residue]['phi'].toFixed(1)}, ${residue_info[residue]['psi'].toFixed(1)})</td></tr>`}
                if (residue_info[residue].hasOwnProperty('contact_numbers')) {tooltip_HTML_body += `<tr><td>Contacts With</td><td>${residue_info[residue]['contact_numbers']}</td></tr>`}

                if ($('#backbone_color_type').val() == 'Dynamic') {
                    let name = $('#backbone_color_variable').val();
                    if (!tooltip_HTML_body.includes(name)) {
                        let backbone_val;
                        if (residue_info[residue].hasOwnProperty(name)) {
                            // The selected dynamic variable is in residue_info
                            backbone_val = residue_info[residue][name];
                        } else if (variability[group].hasOwnProperty(residue)) {
                            // The selected dynamic variable is in the variability data
                            backbone_val = variability[group][residue][name];
                        } else {
                            // It's in neither. Not good
                            backbone_val = null;
                        }
                        tooltip_HTML_body += `<tr><td>${name}</td><td>${backbone_val}</td></tr>`
                    }
                }
                if ($('#surface_color_type').val() == 'Dynamic') {
                    let name = $('#surface_color_variable').val();
                    if (!tooltip_HTML_body.includes(name)) {
                        let surface_val;
                        if (residue_info[residue].hasOwnProperty(name)) {
                            // The selected dynamic variable is in residue_info
                            surface_val = residue_info[residue][name];
                        } else if (variability[group].hasOwnProperty(residue)) {
                            // The selected dynamic variable is in the variability data
                            surface_val = variability[group][residue][name];
                        } else {
                            // It's in neither. Not good
                            surface_val = null;
                        }
                        tooltip_HTML_body += `<tr><td>${name}</td><td>${surface_val}</td></tr>`
                    }
                }

                // Variant data is available if a variant exists at hovered residue
                if (variability[group].hasOwnProperty(residue)) {
                    var tooltip_HTML_variant_title = `<h5>Variant info</h5>`
                    var tooltip_HTML_variant_body = `
                        <tr><td>Mean Dfc</td><td>${variability[group][residue]['departure_from_consensus'].toFixed(2)}</td></tr>
                        <tr><td>Prevalence</td><td>${variability[group][residue]['occurrence']} of ${parseInt(Math.round(variability[group][residue]['occurrence'] / variability[group][residue]['prevalence']))} samples</td></tr>
                        <tr><td>Site Coverage</td><td>${variability[group][residue]['coverage'].toFixed(2)}</td></tr>
                        ${variability[group][residue].hasOwnProperty('mean_normalized_coverage') ? `<tr><td>Site Coverage / Gene Coverage</td><td>${variability[group][residue]['mean_normalized_coverage'].toFixed(2)}</td></tr>` : ``}
                        <tr><td>Mean Entropy</td><td>${variability[group][residue]['entropy'].toFixed(2)}</td></tr>
                    `
                    // add engine-specific data
                    if ($('[name=engine]:checked').val() == 'AA') {
                        // append to body
                        tooltip_HTML_variant_body += `<tr><td>Mean BLOSUM90</td><td>${variability[group][residue]['BLOSUM90'].toFixed(1)}</td></tr>`
                    } else {
                        // append to body
                        tooltip_HTML_variant_body += `<tr><td>log10(pN) [popular consensus])</td><td>${variability[group][residue]['log_pN_popular_consensus'].toFixed(4)}</td></tr>`
                        tooltip_HTML_variant_body += `<tr><td>log10(pS) [popular consensus])</td><td>${variability[group][residue]['log_pS_popular_consensus'].toFixed(4)}</td></tr>`
                    }

                    var tooltip_HTML_variant_freqs_title = `<h5>Variant frequencies</h5>`
                    if ($('[name=engine]:checked').val() == 'AA') {
                        var tooltip_HTML_variant_freqs_body = `
                            <tr><td>${variability[group][residue]['0_item']}</td><td>${variability[group][residue]['0_freq'].toFixed(3)}</td></tr>
                            <tr><td>${variability[group][residue]['1_item']}</td><td>${variability[group][residue]['1_freq'].toFixed(3)}</td></tr>
                            <tr><td>${variability[group][residue]['2_item']}</td><td>${variability[group][residue]['2_freq'].toFixed(3)}</td></tr>
                            <tr><td>${variability[group][residue]['3_item']}</td><td>${variability[group][residue]['3_freq'].toFixed(3)}</td></tr>
                        `
                    } else {
                        var tooltip_HTML_variant_freqs_body = `
                            <tr><td>${variability[group][residue]['0_item_AA']} (${variability[group][residue]['0_item']})</td><td>${variability[group][residue]['0_freq'].toFixed(3)}</td></tr>
                            <tr><td>${variability[group][residue]['1_item_AA']} (${variability[group][residue]['1_item']})</td><td>${variability[group][residue]['1_freq'].toFixed(3)}</td></tr>
                            <tr><td>${variability[group][residue]['2_item_AA']} (${variability[group][residue]['2_item']})</td><td>${variability[group][residue]['2_freq'].toFixed(3)}</td></tr>
                            <tr><td>${variability[group][residue]['3_item_AA']} (${variability[group][residue]['3_item']})</td><td>${variability[group][residue]['3_freq'].toFixed(3)}</td></tr>
                        `
                    }
                }

                if ($('#show_tooltip').is(':checked') && $('#show_tooltip_when').val() == 'all residues') {
                    tooltip_HTML_body = `<table class="tooltip-table">` + tooltip_HTML_body + `</table>`
                    tooltip_HTML = tooltip_HTML_title + tooltip_HTML_body

                    // Variant data is available if a variant exists at hovered residue
                    if (variability[group].hasOwnProperty(residue)) {
                        tooltip_HTML_variant_body = `<table class="tooltip-table">` + tooltip_HTML_variant_body + `</table>`
                        tooltip_HTML += tooltip_HTML_variant_title + tooltip_HTML_variant_body
                        tooltip_HTML_variant_freqs_body = `<table class="tooltip-table">` + tooltip_HTML_variant_freqs_body + `</table>`
                        tooltip_HTML += tooltip_HTML_variant_freqs_title + tooltip_HTML_variant_freqs_body
                    }

                    tooltip.innerHTML = tooltip_HTML;
                    tooltip.style.bottom = window.innerHeight - mp.y + 3 + "px";
                    tooltip.style.left = mp.x + 3 + "px";
                    tooltip.style.display = "block";
                }
                else if ($('#show_tooltip').is(':checked') && $('#show_tooltip_when').val() == 'variant residues') {
                    if (variability[group].hasOwnProperty(residue)) {
                        tooltip_HTML_body = `<table class="tooltip-table">` + tooltip_HTML_body + `</table>`
                        tooltip_HTML = tooltip_HTML_title + tooltip_HTML_body

                        tooltip_HTML_variant_body = `<table class="tooltip-table">` + tooltip_HTML_variant_body + `</table>`
                        tooltip_HTML += tooltip_HTML_variant_title + tooltip_HTML_variant_body
                        tooltip_HTML_variant_freqs_body = `<table class="tooltip-table">` + tooltip_HTML_variant_freqs_body + `</table>`
                        tooltip_HTML += tooltip_HTML_variant_freqs_title + tooltip_HTML_variant_freqs_body

                        tooltip.innerHTML = tooltip_HTML;
                        tooltip.style.bottom = window.innerHeight - mp.y + 3 + "px";
                        tooltip.style.left = mp.x + 3 + "px";
                        tooltip.style.display = "block";
                    }
                }

                previous_hovered_residue = residue;
            } else {
                tooltip.style.display = "none";
                if ($('#show_ballstick').is(':checked') && !(['always', 'variant residue'].includes($('#show_ballstick_when').val()))){
                    stage.compList[0].reprList.slice(0).forEach((rep) => {
                        if (rep.name == 'ball+stick') {
                            rep.dispose();
                        }
                    });
                }
            }
        });

        let func = () => { apply_orientation_matrix_to_all_stages( stage.viewerControls.getOrientation()); };
        stage.mouseObserver.signals.scrolled.add(() => { func(); });
        stage.mouseObserver.signals.dragged.add(() => { func(); });
        $(`#ngl_${group}`).mouseup(func);

        stages[group] = stage;

        defer.resolve();
    });

    return defer.promise();
}

function calcBackboneColor(residue, group=null) {
    let name = $('#backbone_color_variable').val();
    let val;

    if (residue_info[residue].hasOwnProperty(name)) {
        // The selected dynamic variable is in residue_info
        val = residue_info[residue][name];
    } else if (variability[group].hasOwnProperty(residue)) {
        // The selected dynamic variable is in the variability data
        val = variability[group][residue][name];
    } else {
        // It's in neither. Not good
        val = null;
    }

    if (val == null) {
        // This can be true either because the value was absent, or the name requested was
        // invalid. In this case we return the min value color
        return '0x' + $('#backbone_color_start').attr('color').substring(1, 7).toUpperCase()
    }

    let min_value = parseFloat($('#backbone_color_min').val());
    let max_value = parseFloat($('#backbone_color_max').val());
    let val_normalized = (parseFloat(val) - min_value) / (max_value - min_value);
    val_normalized = Math.max(0, Math.min(1, val_normalized));

    var hex = getGradientColor($('#backbone_color_start').attr('color'), $('#backbone_color_end').attr('color'), val_normalized);
    return '0x' + hex.substring(1, 7).toUpperCase()
}

function getBackboneColorScheme(group=null) {
    // group is only needed if the selected color variable has a group-specific value, i.e. entropy
    // is a group specific parameter but predicted ligand binding frequency is not

    var schemeId_backbone = NGL.ColormakerRegistry.addScheme(function (params) {
      this.atomColor = function (atom) {
        return calcBackboneColor(atom.resno, group);
      }
    })

    return schemeId_backbone;
}

function calcSurfaceColor(residue, group=null) {
    let name = $('#surface_color_variable').val();
    let val;

    if (residue_info[residue].hasOwnProperty(name)) {
        // The selected dynamic variable is in residue_info
        val = residue_info[residue][name];
    } else if (variability[group].hasOwnProperty(residue)) {
        // The selected dynamic variable is in the variability data
        val = variability[group][residue][name];
    } else {
        // It's in neither. Not good
        val = null;
    }

    if (val == null) {
        // This can be true either because the value was absent, or the name requested was
        // invalid. In this case we return the min value color
        return '0x' + $('#surface_color_start').attr('color').substring(1, 7).toUpperCase()
    }

    let min_value = parseFloat($('#surface_color_min').val());
    let max_value = parseFloat($('#surface_color_max').val());
    let val_normalized = (parseFloat(val) - min_value) / (max_value - min_value);
    val_normalized = Math.max(0, Math.min(1, val_normalized));

    var hex = getGradientColor($('#surface_color_start').attr('color'), $('#surface_color_end').attr('color'), val_normalized);
    return '0x' + hex.substring(1, 7).toUpperCase()
}

function getSurfaceColorScheme(group=null) {
    // group is only needed if the selected color variable has a group-specific value, i.e. entropy
    // is a group specific parameter but predicted ligand binding frequency is not

    var schemeId_surface = NGL.ColormakerRegistry.addScheme(function (params) {
      this.atomColor = function (atom) {
        return calcSurfaceColor(atom.resno, group);
      }
    })

    return schemeId_surface;
}

function backbone_rule(element) {
    if ($(element).val() == 'Static') {
        $('.static-color-backbone').show();
        $('.dynamic-color-backbone').hide();
    } else if ($(element).val() == 'Dynamic') {
        $('.static-color-backbone').hide();
        $('.dynamic-color-backbone').show();
    } else {
        $('.static-color-backbone').hide();
        $('.dynamic-color-backbone').hide();
    }

    create_ngl_views(fetch_variability=false);
}


function surface_rule(element) {
    if ($(element).val() == 'Static') {
        $('.static-color-surface').show();
        $('.dynamic-color-surface').hide();
    } else if ($(element).val() == 'Dynamic') {
        $('.static-color-surface').hide();
        $('.dynamic-color-surface').show();
    } else {
        $('.static-color-surface').hide();
        $('.dynamic-color-surface').hide();
    }

    create_ngl_views(fetch_variability=false);
}


function load_protein() {
    let gene_callers_id = $('#gene_callers_id_list').val();
    var defer = $.Deferred();
    $('.overlay').show();

    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/get_structure/' + gene_callers_id,
        success: function(data) {
            $('.overlay').hide();
            histogram_data = data['histograms'];
            pdb_content = data['pdb_content'];
            residue_info = move_codon_number_to_index(JSON.parse(data['residue_info']));
            residue_info_types = JSON.parse(data['residue_info_types']);
            defer.resolve();
        }
    });

    return defer.promise();
}

function load_gene_function_info() {
    let gene_callers_id = $('#gene_callers_id_list').val();
    var defer = $.Deferred();

    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/get_gene_function_info/' + gene_callers_id,
        success: function(gene_data) {
            var geneFunctionHtml = get_gene_functions_table_html_for_structure(gene_data);
            $("#gene_function_info").html(geneFunctionHtml);
            defer.resolve();
        }
    });

    return defer.promise();
}

function load_model_info() {
    let gene_callers_id = $('#gene_callers_id_list').val();
    var defer = $.Deferred();

    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/get_model_info/' + gene_callers_id,
        success: function(model_data) {
            model_data = model_data;
            var geneModelHtml = get_model_info_table_html(model_data);
            $("#model_info").html(geneModelHtml);
            defer.resolve();
        },
        error: function(xhr, status, error) {
            console.error("Error loading model info:", status, error);
            defer.reject();
        }
    });

    return defer.promise();
}

function get_model_info_table_html(model_data) {
    var templates = JSON.parse(model_data['templates']);
    var models = JSON.parse(model_data['models'])[0];

    var geneModelHtml = '';

    /* TEMPLATES */
    geneModelHtml += '<div class="widget">'
    geneModelHtml += '<span class="settings-header"><h4>Templates Used</h4></span>'
    geneModelHtml += '<table class="table table-sm table-responsive" id="model_info_table"><tbody>';

    var header = '<tr>';
    for (const col_name of Object.keys(templates[0])) {
        header += '<td><label class="col-md-4 settings-label">' + col_name + '</label></td>';
    }
    header += '</tr>';
    geneModelHtml += header;

    for (var i = 0; i < Object.keys(templates).length; i++) {
        var row = templates[i];
        var tr = '<tr>';

        for (const [col_name, value] of Object.entries(row)) {
            if (col_name == 'PDB') {
                formatted_value = '<a href="https://www.rcsb.org/structure/' + value + '" target=_"blank">' + value.toUpperCase() + '</a>'
            } else if (col_name == '%Identity') {
                formatted_value = Number(value).toFixed(2);
            } else if (col_name == 'Align fraction') {
                formatted_value = Number(value).toFixed(3);
            } else {
                formatted_value = value
            }
            tr += '<td>' + formatted_value + '</td>'
        }

        tr += '</tr>'
        geneModelHtml += tr;
    }

    geneModelHtml += "</tbody></table></div>";

    /* MODELS */
    geneModelHtml += '<div class="widget">'
    geneModelHtml += '<span class="settings-header"><h4>Model Scores</h4></span>'
    geneModelHtml += '<table class="table table-sm table-responsive" id="model_info_table"><tbody>';

    var header = '<tr>';
    var row = '<tr>';
    for (const [col_name, value] of Object.entries(models)) {
        header += '<td><label class="col-md-4 settings-label">' + col_name + '</label></td>';
        row += '<td>' + Number(value).toFixed(2) + '</td>';
    }
    header += '</tr>';
    row += '</tr>';

    geneModelHtml += header;
    geneModelHtml += row;

    geneModelHtml += "</tbody></table></div>";
    return geneModelHtml;
}



function serialize_checked_groups() {
    let output = {};

    let category = $('#sample_groups_list').val();
    for (let group in sample_groups[category]) {
        if ($(`input[type=checkbox]#${category}_${group}`).is(':checked')) {
            output[group] = [];

            sample_groups[category][group].forEach((sample) => {
                if ($(`input[type=checkbox]#${category}_${group}_${sample}`).is(':checked')) {
                    output[group].push(sample);
                }
            });
        }
    }

    return output;
}

function serialize_filtering_widgets() {
    let output = {};

    $('#controls .widget').each((index, widget) => {
        let column = $(widget).attr('data-column');
        let controller = $(widget).attr('data-controller');

        if (controller == 'slider') {
            output[column] = {}
            output[column]["min_" + column] = $(widget).find('input').val().split(',')[0];
            output[column]["max_" + column] = $(widget).find('input').val().split(',')[1];
        }
        else if (controller == 'checkbox') {
            output[column] = {}
            output[column][column + "s_of_interest"] = $(widget).find('input:checkbox:checked').toArray().map((checkbox) => { return $(checkbox).val(); });
        }
    });

    return output;
}

function fetch_and_draw_variability() {
    $('.overlay').show();
    let gene_callers_id = $('#gene_callers_id_list').val();
    let engine = $('[name=engine]:checked').val();

    // serialize options programatically
    let options = {
        'gene_callers_id': gene_callers_id,
        'engine': engine,
        'groups': serialize_checked_groups(),
        'filter_params': serialize_filtering_widgets()
    };

    $.ajax({
        type: 'POST',
        cache: false,
        data: {'options': JSON.stringify(options)},
        url: '/data/get_variability',
        success: function(response_all) {
            $('.overlay').hide();

            variability = {};
            for (let group in response_all) {
                let response = response_all[group];

                let data = JSON.parse(response['data']);
                let total_entries = response['total_entries'];
                let entries_after_filtering = response['entries_after_filtering'];
                variability[group] = move_codon_number_to_index(data);
            }

            draw_variability();
        },
        error: function(request, status, error) {
            console.log(request, status, error);
            $('.overlay').hide();
        }
    });
}

function draw_variability() {
    let gene_callers_id = $('#gene_callers_id_list').val();
    let engine = $('[name=engine]:checked').val();

    if (Object.keys(stages).length == 0)
        return;

    for (let group in variability) {
        if (!stages.hasOwnProperty(group))
            return;

        let data = variability[group];
        let compList = stages[group].compList;

        if (compList.length == 0)
            return;

        let component = compList[0];

        if (!component)
            return;

        component.reprList.slice(0).forEach((rep) => {
            if (rep.name == 'spacefill') {
                rep.dispose();
            }
            if (rep.name == 'ball+stick' && $('#show_ballstick_when').val() == 'variant residue') {
                rep.dispose();
            }
        });

        if (Object.keys(data).length > 0) {
            for (let index in data) {
                let spacefill_options = {
                    sele: data[index]['codon_number'] + " and .CA",
                    scale: 1
                }

                if ($('#color_type').val() == 'Dynamic') {
                    let column = $('#color_target_column').val();
                    let column_value = data[index][column];
                    let selected_column_info = column_info.find(function(el) {if (el['name'] == column) {return el}})

                    if (selected_column_info['data_type'] == 'float' || selected_column_info['data_type'] == 'integer') {
                        let min_value = parseFloat($('#color_min').val());
                        let max_value = parseFloat($('#color_max').val());
                        let middle_value = (min_value + max_value) / 2;

                        if (min_value >= max_value) {
                            $('#dynamic_color_error').show();
                        } else {
                            $('#dynamic_color_error').hide();
                        }

                        let val = (parseFloat(column_value) - min_value) / (max_value - min_value);
                        val = Math.max(0, Math.min(1, val));

                        if ($('#show_3_color').is(':checked')) {
                            if (val <= 0.5) {
                                spacefill_options['color'] = getGradientColor(
                                    $('#color_start').attr('color'),
                                    $('#color_middle').attr('color'),
                                    val);
                            } else {
                                spacefill_options['color'] = getGradientColor(
                                    $('#color_middle').attr('color'),
                                    $('#color_end').attr('color'),
                                    val);
                            }
                        } else {
                            spacefill_options['color'] = getGradientColor(
                                $('#color_start').attr('color'),
                                $('#color_end').attr('color'),
                                val);
                        }
                    }
                    else
                    {   
                        debugger;
                        spacefill_options['color'] = color_legend[engine][column][column_value];
                    }
                } else {
                    spacefill_options['color'] = $('#color_static').attr('color');
                }

                if ($('#size_type').val() == 'Dynamic') {
                    let column = $('#size_target_column').val();
                    let column_value = data[index][column];
                    let selected_column_info = column_info.find(function(el) {if (el['name'] == column) {return el}})

                    if (selected_column_info['data_type'] == 'float' || selected_column_info['data_type'] == 'integer') {
                        let min_value = parseFloat($('#size_min').val());
                        let max_value = parseFloat($('#size_max').val());
                        let start_value = parseFloat($('#size_start').val());
                        let end_value = parseFloat($('#size_end').val());

                        if (min_value >= max_value) {
                            $('#dynamic_size_error').show();
                        } else {
                            $('#dynamic_size_error').hide();
                        }

                        let val = (parseFloat(column_value) - min_value) / (max_value - min_value);
                        val = Math.max(0, Math.min(1, val));

                        spacefill_options['scale'] = start_value + (val * (end_value - start_value));
                    }
                    else {
                        spacefill_options['scale'] = size_legend[engine][column][column_value];
                    }
                } else {
                    spacefill_options['scale'] = parseFloat($('#size_static').val());
                }

                let representation = component.addRepresentation("spacefill", spacefill_options);
                representation['variability'] = data[index];

                if ($('#show_ballstick').is(':checked')) {
                    if ($('#show_ballstick_when').val() == 'variant residue') {
                        var selection = data[index]['codon_number'] + " and sidechainAttached";
                        component.addRepresentation("ball+stick", {
                            sele: selection
                        });
                    }
                }
            }
        }
    }
};


function draw_histogram() {
    let engine = $('[name=engine]:checked').val();

    for (let column in histogram_data[engine]) {
        let svg = d3.select('#histogram_' + column);

        if (svg.empty()) {
            continue;
        }

        let width = 200;
        let height = 30;

        // clear existing drawing
        svg.selectAll('*').remove();

        let bins = histogram_data[engine][column]['bins'];
        let counts = histogram_data[engine][column]['counts'];

        var min_count = Math.min(...counts);
        var max_count = Math.max(...counts);
        var max_slider = parseFloat(document.getElementById(column).dataset.sliderMax);
        var min_slider = parseFloat(document.getElementById(column).dataset.sliderMin);

        let normalized_counts = counts.map(v => (v / max_count) * height);
        let normalized_bins = bins.map(v => ((v - min_slider) / (max_slider - min_slider)) * width);

        let data_points = [];

        for (let i=0; i < normalized_counts.length; i++) {
            data_points.push({'x': normalized_bins[i], 'y': height - normalized_counts[i]});
            data_points.push({'x': normalized_bins[i+1], 'y': height - normalized_counts[i]});
        }
        data_points.push({'x': width, 'y': height})

        var make_bar_chart = d3.line()
                            .x(function(d) { return d.x; })
                            .y(function(d) { return d.y; });

        var bar_chart = make_bar_chart(data_points);
        bar_chart += `L ${data_points[normalized_bins.length - 1]['x']} ${height} L ${data_points[0]['x']} ${height}`;

        svg.append("path")
            .style("fill","#337ab7")
            .style("stroke","#182943")
            .attr("d",function(d,i){ return bar_chart; });
    }
};


function create_ui() {
    var defer = $.Deferred();
    let gene_callers_id = $('#gene_callers_id_list').val();
    let engine = $('[name=engine]:checked').val();

    backupFilters();

    $.ajax({
        type: 'POST',
        cache: false,
        url: '/data/get_column_info',
        data: {
            'gene_callers_id': gene_callers_id,
            'engine': engine
        },
        success: function(data) {
            column_info = data
            let container = $('#controls');
            $('#controls').attr('created-for-engine', engine);
            $('#controls').attr('created-for-gene-caller-id', gene_callers_id);

            // remove widgets
            container.empty();

            $('#color_target_column').empty();
            $('#size_target_column').empty();
            $('#backbone_color_variable').empty();
            $('#surface_color_variable').empty();

            if (!color_legend.hasOwnProperty(engine)) {
                color_legend[engine] = {};
            }

            if (!size_legend.hasOwnProperty(engine)) {
                size_legend[engine] = {};
            }

            for (const [info_name, types] of Object.entries(residue_info_types)) {
                if (types["dtype"] != "object") {
                    // Ok this is some type of numerical data form. We can include it
                    $('#backbone_color_variable').append(`<option value="${info_name}">${info_name}</item>`);
                    $('#surface_color_variable').append(`<option value="${info_name}">${info_name}</item>`);
                }
            }
            column_info.forEach((item) => {
                if ((item['as_view'] | item['as_filter']) & item['data_type'] != 'text') {
                    $('#backbone_color_variable').append(`<option value="${item['name']}">${item['title']}</item>`);
                    $('#surface_color_variable').append(`<option value="${item['name']}">${item['title']}</item>`);
                }
            })

            column_info.forEach((item) => {
                if (item['as_view']) {
                    $('#color_target_column').append(`<option value="${item['name']}">${item['title']}</item>`);
                    $('#size_target_column').append(`<option value="${item['name']}">${item['title']}</item>`);
                }

                if (item['as_filter'] == 'slider') {
                    let min_val = item['min'];
                    let max_val = item['max'];

                    if (filter_backup.hasOwnProperty(gene_callers_id) && filter_backup[gene_callers_id].hasOwnProperty(engine) && filter_backup[gene_callers_id][engine].hasOwnProperty(item['name'])) {
                        min_val = filter_backup[gene_callers_id][engine][item['name']]['min_' + item['name']];
                        max_val = filter_backup[gene_callers_id][engine][item['name']]['max_' + item['name']];
                    }

                    $(container).append(`
                        <div class="widget" data-column="${item['name']}" data-controller="${item['as_filter']}">
                            <span class="settings-header"><h5>${item['title']}</h5></span><br />
                            <svg id="histogram_${item['name']}" width="100%" height="30" style="position: relative; top: 6;" viewBox="0 0 200 30" preserveAspectRatio="none"></svg>
                            <input id="${item['name']}"
                                    type="${item['data_type']}"
                                    data-provide="slider"
                                    data-slider-min="${item['min']}"
                                    data-slider-max="${item['max']}"
                                    data-slider-step="${item['step']}"
                                    data-slider-value="[${min_val},${max_val}]"
                                    >
                        </div>
                    `);
                    $(`#${item['name']}`).slider({});
                }
                if (item['as_filter'] == 'checkbox') {
                    let checked_choices = item['choices'];

                    if (filter_backup.hasOwnProperty(gene_callers_id) && filter_backup[gene_callers_id].hasOwnProperty(engine) && filter_backup[gene_callers_id][engine].hasOwnProperty(item['name'])) {
                        checked_choices = filter_backup[gene_callers_id][engine][item['name']][item['name'] + 's_of_interest'];
                    }

                    $(container).append(`
                        <div class="widget" data-column="${item['name']}" data-controller="${item['as_filter']}">
                            <span class="settings-header"><h5>${item['title']}</h5></span><br />
                            <div class="ml-3 d-flex">
                            ${item['choices'].map((choice) => { return `
                                <div>
                                    <input class="form-check-input" type="checkbox" id="${item['name']}_${choice}" value="${choice}" ${ checked_choices.indexOf(choice) > -1 ? 'checked="checked"' : ''}>
                                    <label class="form-check-label" for="${item['name']}_${choice}">${choice}</label>`; }).join('')}
                                </div>    
                                <br />
                                <div>
                                    <button class="btn btn-xs" onclick="$(this).closest('.widget').find('input:checkbox').prop('checked', true);">Check All</button>
                                    <button class="btn btn-xs" onclick="$(this).closest('.widget').find('input:checkbox').prop('checked', false);">Uncheck All</button>
                                </div>
                            </div>
                        </div>
                    `);

                    if (!color_legend[engine].hasOwnProperty(item['name'])) {
                        color_legend[engine][item['name']] = {};

                        item['choices'].forEach((choice) => {
                            color_legend[engine][item['name']][choice] = randomColor();
                        });
                    }

                    if (!size_legend[engine].hasOwnProperty(item['name'])) {
                        size_legend[engine][item['name']] = {};

                        item['choices'].forEach((choice) => {
                            size_legend[engine][item['name']][choice] = 1;
                        });
                    }
                }
            });

            draw_histogram();
            defer.resolve();
        }
    });

   return defer.promise();
}


function onTargetResidueInfoChange(element) {
    let name = $(element).val();

    $(`#backbone_numerical_panel`).show();

    if (residue_info[1].hasOwnProperty(name)) {
      // The selected dynamic variable is in residue_info's elements
      $(`#backbone_color_min`).val(residue_info_types[name]['amin']);
      $(`#backbone_color_max`).val(residue_info_types[name]['amax']);
    } else {
      for (i in column_info) {
        let item = column_info[i];
        if (item['name'] == name) {
          $(`#backbone_color_min`).val(item['min']);
          $(`#backbone_color_max`).val(item['max']);
          break;
        }
      }
    }

    $('.colorpicker-legend').colpick({
        layout: 'hex',
        submit: 0,
        colorScheme: 'light',
        onChange: function(hsb, hex, rgb, el, bySetColor) {
            $(el).css('background-color', '#' + hex);
            $(el).attr('color', '#' + hex);

            if (!bySetColor) $(el).val(hex);
        },
        onHide: function(cal) {
            let el = $(cal).data('colpick').el;
            color_legend[$(el).attr('data-engine')][$(el).attr('data-column')][$(el).attr('data-key')] = $(el).attr('color');
            draw_variability();
        }
    }).keyup(function() {
        $(this).colpickSetColor(this.value);
    });

    create_ngl_views(fetch_variability=false);
}


function onTargetColumnChange(element) {
    // this on change event is shared between color_target_column and size_target_column.
    let engine = $('[name=engine]:checked').val();
    let column = $(element).val();
    let selected_column_info = column_info.find(function(el) {if (el['name'] == column) {return el}})

    // color or size
    let prefix = element.getAttribute('id').split('_')[0];

    // show column related panel.
    // for linear values, show numerical panel
    // for discreete values, show legend panel
    // read the min/max from numerical and put into prefixed input in perspective
    if (selected_column_info['data_type'] == 'float' || selected_column_info['data_type'] == 'integer') {
        $(`#${prefix}_numerical_panel`).show();
        $(`#${prefix}_legend_panel`).hide();

        $(`#${prefix}_min`).val(selected_column_info['min']);
        $(`#${prefix}_max`).val(selected_column_info['max']);
    }
    else
    {
        $(`#${prefix}_numerical_panel`).hide();
        $(`#${prefix}_legend_panel`).show();

        // populate color legend.
        $(`#${prefix}_legend_panel`).empty();

        let legend_items = window[`${prefix}_legend`][engine][column];

        for (let key in legend_items) {
            let value = legend_items[key];

            if (prefix == 'color') {
                $(`#color_legend_panel`).append(`<div class="col-md-4">
                                                 <div class="colorpicker colorpicker-legend"
                                                      color="${value}"
                                                      style="background-color:${value}"
                                                      data-engine="${engine}"
                                                      data-column="${column}"
                                                      data-key=${key}>
                                                 </div><span style="position: relative; top: -0.4em;">${key}</span>
                                                 </div>`);
            } else {
                $(`#size_legend_panel`).append(`<div class="col-md-6 row no-gutter" style="padding: 3px;">
                                                    <div class="col-md-5">
                                                        <input type="text" value="${value}" class="form-control input-sm"
                                                        onblur="size_legend['${engine}']['${column}']['${key}'] = parseFloat(this.value); draw_variability();">
                                                    </div>
                                                    <div class="col-md-7" style="text-align: left; padding-left: 4px;">
                                                        ${key}
                                                    </div>
                                                </div>`);
            }
        }
    }

    $('.colorpicker-legend').colpick({
        layout: 'hex',
        submit: 0,
        colorScheme: 'light',
        onChange: function(hsb, hex, rgb, el, bySetColor) {
            $(el).css('background-color', '#' + hex);
            $(el).attr('color', '#' + hex);

            if (!bySetColor) $(el).val(hex);
        },
        onHide: function(cal) {
            let el = $(cal).data('colpick').el;
            color_legend[$(el).attr('data-engine')][$(el).attr('data-column')][$(el).attr('data-key')] = $(el).attr('color');
            draw_variability();
        }
    }).keyup(function() {
        $(this).colpickSetColor(this.value);
    });

    draw_variability();
}


function getGradientColor(start_color, end_color, percent) {
    // strip the leading # if it's there
    start_color = start_color.replace(/^\s*#|\s*$/g, '');
    end_color = end_color.replace(/^\s*#|\s*$/g, '');

    // convert 3 char codes --> 6, e.g. `E0F` --> `EE00FF`
    if(start_color.length == 3){
      start_color = start_color.replace(/(.)/g, '$1$1');
    }

    if(end_color.length == 3){
      end_color = end_color.replace(/(.)/g, '$1$1');
    }

    // get colors
    var start_red = parseInt(start_color.substr(0, 2), 16),
        start_green = parseInt(start_color.substr(2, 2), 16),
        start_blue = parseInt(start_color.substr(4, 2), 16);

    var end_red = parseInt(end_color.substr(0, 2), 16),
        end_green = parseInt(end_color.substr(2, 2), 16),
        end_blue = parseInt(end_color.substr(4, 2), 16);

    // calculate new color
    var diff_red = end_red - start_red;
    var diff_green = end_green - start_green;
    var diff_blue = end_blue - start_blue;

    diff_red = Math.abs(( (diff_red * percent) + start_red )).toString(16).split('.')[0];
    diff_green = Math.abs(( (diff_green * percent) + start_green )).toString(16).split('.')[0];
    diff_blue = Math.abs(( (diff_blue * percent) + start_blue )).toString(16).split('.')[0];

     // ensure 2 digits by color
    if( diff_red.length == 1 )
      diff_red = '0' + diff_red

    if( diff_green.length == 1 )
      diff_green = '0' + diff_green

    if( diff_blue.length == 1 )
      diff_blue = '0' + diff_blue

    return '#' + diff_red + diff_green + diff_blue;
 };


function move_codon_number_to_index(data) {
    // dataframes converted to JSON in python use index as the key. this converts key to the
    // 'codon_number' entry.
    let codon_indexed = {};
    for (let index in data) {
        codon_indexed[data[index]['codon_number']] = data[index];
    }
    return codon_indexed
}


async function make_image(group, sample) {
    let blob;
    let image_options = {
        'trim': $('#trim').is(':checked'),
        'factor': $('#scale').val(),
        'transparent': $('#transparent').is(':checked'),
        'antialias': $('#antialias').is(':checked')
    }

    if (typeof sample === 'undefined') {
        // no sample requested, we generate image for merged.
        blob = await stages[group].viewer.makeImage(image_options);
    }
    else {
        // sample requested.
        let listRepresentations = stages[group].compList[0].reprList.slice(0);

        // hide all representations besides what we want.
        // variability_information dictionary linked during creation of representation
        // in draw_variability.
        listRepresentations.forEach((rep) => {
            if (rep.name == 'spacefill') {
                if (rep.variability.sample_ids.split(', ').indexOf(sample) === -1) {
                    rep.setVisibility(false);
                }
            }
        });

        // take the image
        blob = await stages[group].viewer.makeImage(image_options);

        // restore representations
        listRepresentations.forEach((rep) => {
            if (rep.name == 'spacefill') {
                if (rep.variability.sample_ids.split(', ').indexOf(sample) === -1) {
                    rep.setVisibility(true);
                }
            }
        });
    }

    return blob;
}

function get_gene_functions_table_html_for_structure(gene){
    if (gene.functions == null) {
        functions_table_html = '<div class="alert alert-info" role="alert" style="margin: 10px;" id="no_function_annotations_warning">Functional annotations: This gene doesn\'t have any.</div>';
        return functions_table_html
    }

    functions_table_html  = '<table class="table table-striped table-responsive">';
    functions_table_html += '<thead><th>Source</th>';
    functions_table_html += '<th>Accession</th>';
    functions_table_html += '<th>Annotation</th></thead>';
    functions_table_html += '<tbody>';

    for (function_source in gene.functions){
        functions_table_html += '<tr>';

        functions_table_html += '<td><b>' + function_source + '</b></td>';
        if (gene.functions[function_source]) {
            functions_table_html += '<td>' + decorateAccession(function_source, gene.functions[function_source][0]) + '</td>';
            functions_table_html += '<td><em>' + decorateAnnotation(function_source, gene.functions[function_source][1]) + '</em></td>';
        } else {
            functions_table_html += '<td>&nbsp;</td>';
            functions_table_html += '<td>&nbsp;</td>';
        }

        functions_table_html += '</tr>';
    }

    functions_table_html += '</tbody></table>';

    return functions_table_html;
}

function store_variability() {
    $('.overlay').show();
    let gene_callers_id = $('#gene_callers_id_list').val();
    let engine = $('[name=engine]:checked').val();
    let output_path = $('#var_output_path').val();

    // serialize options programatically
    let options = {
        'gene_callers_id': gene_callers_id,
        'engine': engine,
        'groups': serialize_checked_groups(),
        'filter_params': serialize_filtering_widgets(),
        'path': output_path,
    };

    $.ajax({
        type: 'POST',
        cache: false,
        data: {'options': JSON.stringify(options)},
        url: '/data/store_variability',
        success: function(msg) {
            $('.overlay').hide();
            if (typeof(msg['success']) != 'undefined') {
                $('#store_var_failure').hide();
                $('#store_var_success').html(msg['success']).show().fadeOut(3000);
            } else {
                $('#store_var_failure').html(msg['failure']).show();
            }
        },
        error: function(request, status, error) {
            console.log(request, status, error);
            $('.overlay').hide();
        }
    });
}

function store_structure_as_pdb(path_id, success_id, failure_id) {
    $('.overlay').show();
    let gene_callers_id = $('#gene_callers_id_list').val();
    let output_path = $(path_id).val();

    // serialize options programatically
    let options = {
        'gene_callers_id': gene_callers_id,
        'path': output_path,
    };

    $.ajax({
        type: 'POST',
        cache: false,
        data: {'options': JSON.stringify(options)},
        url: '/data/store_structure_as_pdb',
        success: function(msg) {
            $('.overlay').hide();
            if (typeof(msg['success']) != 'undefined') {
                $(failure_id).hide();
                $(success_id).html(msg['success']).show().fadeOut(3000);
            } else {
                $(failure_id).html(msg['failure']).show();
            }
        },
        error: function(request, status, error) {
            console.log(request, status, error);
            $('.overlay').hide();
        }
    });
}

function showPymolWindow() {
    gen_pymol_script_html(gen_pymol_script());
    $('#pymol_export_page').modal('show');
}

function gen_pymol_script_html(script) {
    var pymol_script_html = `
    <div class="modal-body">
        <textarea class="form-control" style="width: 100%; height: 100%; font-family: "Roboto", Helvetica, Arial;" rows="16" onclick="$(this).select();" readonly>${script}</textarea>
    </div>
    `

    $("#pymol_script_area").html(pymol_script_html);
}

function gen_pymol_script() {
    // from https://stackoverflow.com/questions/5623838/rgb-to-hex-and-hex-to-rgb
    function componentToHex(c) {
        var hex = c.toString(16);
        if (hex.length < 6) {
            hex = "0".repeat(6 - hex.length) + hex
        }
        return hex;
    }

    // from https://stackoverflow.com/questions/5623838/rgb-to-hex-and-hex-to-rgb
    // NOTE: RGB values are normalized to 1 for PyMOL convention
    function hexToRgb(hex) {
        var result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
        return result ? {
            r: parseInt(result[1], 16) / 256,
            g: parseInt(result[2], 16) / 256,
            b: parseInt(result[3], 16) / 256
        } : null;
    }

    var prefix = $('#gene_callers_id_list').val();
    var main_object = prefix;
    var surface_transparency = 1 - parseFloat($('#surface_opacity').val());
    var surface_probe_radius = parseFloat($('#surface_probe_radius').val());
    var s = `bg_color white\n` +
            `main_obj = cmd.get_object_list(selection='(${main_object})')[0]\n` +
            `cmd.hide('everything', main_obj)\n` +
            `cmd.show('cartoon', main_obj)\n` +
            `set transparency, ${surface_transparency}\n` +
            `set solvent_radius, ${surface_probe_radius}\n`;

    // s is the full PyMOL script string
    //var bb_color = hexToRgb($('#color_static_backbone').attr('color'));
    //var s += `cmd.set_color('bb_color', [${bb_color.r},${bb_color.g},${bb_color.b}])\n` +
    //        `color bb_color, structure_${$('#gene_callers_id_list').val()}\n`

    // list of the PyMOL objects--one for each group
    var group_object_list = [];

    for (let group in variability) {
        let compList = stages[group].compList;
        let component = compList[0];

        var res_attrs = '';
        var res_list = [];
        var group_object = `${prefix}_${group}`
        var group_var_object = `${prefix}_${group}_var`
        var group_selection = `${prefix}_${group}_sele`
        group_object_list.push(group_object)

        // Initialize a selection with 0 atoms
        s += `sele ${group_selection}, name CA and name CB\n`;

        // Initialize a dictionary to hold each variant's attributes
        s += `res_attrs = {}\n`;

        var counter = 1;
        component.reprList.slice(0).forEach((rep) => {
            if (rep.name == 'spacefill') {
                var res = rep.variability.codon_number;
                var color = hexToRgb(componentToHex(rep.repr.colorValue));
                var scale = rep.repr.scale;

                res_attrs += `${res}:{'color':[${color.r},${color.g},${color.b}],'scale':${scale}},`;
                res_list.push(res.toString());

                if (counter % 20 == 0) {
                    s += `select ${group_selection}, ${group_selection} | (${main_object} and name CA and resi ${res_list.join('+')})\n`;
                    s += `res_attrs.update({${res_attrs}})\n`;

                    res_list = [];
                    res_attrs = '';
                }
                counter += 1;
            }
        });

        s += `res_attrs.update({${res_attrs}})\n` +
             `select ${group_selection}, ${group_selection} | (${main_object} and name CA and resi ${res_list.join('+')})\n` +
             `create ${group_var_object}, ${group_selection}\n` +
             `for res in res_attrs: cmd.set_color('${group_var_object}' + str(res), res_attrs[res]['color'])\n` +
             `alter ${group_var_object}, color = cmd.get_color_index('${group_var_object}' + resi)\n` +
             `alter ${group_var_object}, s.sphere_scale = res_attrs[int(resi)]['scale']\n` +
             `hide everything, ${group_var_object}\n` +
             `show spheres, ${group_var_object}\n` +
             `create ${group_object}, ${group_var_object} or ${main_object}\n`;

        // group's backbone
        if ($('#show_backbone').is(':checked')) {
            if ($('#backbone_color_type').val() == 'Dynamic') {
                for (let residue in residue_info) {
                    var backbone_color = hexToRgb("#".concat(calcBackboneColor(residue, group).substring(2,8)));
                    s += `set_color ${prefix}_${group}_backbone_${residue}, [${backbone_color.r},${backbone_color.g},${backbone_color.b}]\n`;
                    s += `set cartoon_color, ${prefix}_${group}_backbone_${residue}, ${group_object} and resi ${residue}\n`;
                }
            } else {
                var static_backbone_color = hexToRgb($('#color_static_backbone').attr('color'));
                s += `set_color ${prefix}_${group}_backbone_static, [${static_backbone_color.r},${static_backbone_color.g},${static_backbone_color.b}]\n`;
                s += `set cartoon_color, ${prefix}_${group}_backbone_static, ${group_object}\n`;
            }
        } else {
            s += `hide cartoon, ${group_object}\n`
        }

        // group's surface
        if ($('#show_surface').is(':checked')) {
            if ($('#surface_color_type').val() == 'Dynamic') {
                for (let residue in residue_info) {
                    var surface_color = hexToRgb("#".concat(calcSurfaceColor(residue, group).substring(2,8)));
                    s += `set_color ${prefix}_${group}_surf_${residue}, [${surface_color.r},${surface_color.g},${surface_color.b}]\n`;
                    s += `set surface_color, ${prefix}_${group}_surf_${residue}, ${group_object} and resi ${residue}\n`;
                }
            } else {
                var static_surface_color = hexToRgb($('#color_static_surface').attr('color'));
                s += `set_color ${prefix}_${group}_surf_static, [${static_surface_color.r},${static_surface_color.g},${static_surface_color.b}]\n`;
                s += `set surface_color, ${prefix}_${group}_surf_static, ${group_object}\n`;
            }
            s += `show surface, ${group_object}\n`;
        }
    }
    s += `cmd.disable('all')\n` +
         `cmd.enable('${group_object_list[0]}')\n` +
         `orient\n` +
         `rebuild`;

    return s
}

async function generate_summary() {
    $('.overlay').show();

    let serialized_groups = serialize_checked_groups();
    var zip = new JSZip();

    for (let group in stages) {
        zip.file($('#zip_name').val() + `/${group}/merged.png`, await make_image(group));

        if (!$('#merged_view_only').is(':checked')) {
            // generate per sample.
            for (let i=0; i < serialized_groups[group].length; i++) {
                let sample_id = serialized_groups[group][i];
                zip.file($('#zip_name').val() + `/${group}/${sample_id}.png`, await make_image(group, sample_id));
            }
        }
    }

    let zip_options = {
        type: "blob",
        platform: 'UNIX',
        compression: "DEFLATE",
        compressionOptions: {
           level: 6
        }
    }

    zip.generateAsync(zip_options).then(function(content) {
        saveAs(content, $('#zip_name').val() + '.zip');
    });

    $('.overlay').hide();
}

function serializeAuxiliaryInputs() {
    let backup = {};

    ['tab_views', 'tab_output'].forEach((tab) => {
        backup[tab] = {};

        $(`#${tab} :input`).each((index, elem) => {
            let tag = elem.tagName;
            let id = elem.getAttribute('id');

            if (tag && id) {
                if (tag == 'SELECT') {
                    backup[tab][id] = $(elem).val();
                }
                else if (tag == 'INPUT') {
                    let type = elem.getAttribute('type');

                    if (type == 'checkbox') {
                        backup[tab][id] = $(elem).is(':checked');
                    }
                    else if (type == 'input' || type == 'text') {
                        backup[tab][id] = $(elem).val();
                    }
                }
            }
        });
         $(`#${tab} .colorpicker`).each((index, elem) => {
            let id = elem.getAttribute('id');
            if (id) {
                backup[tab][id] = $(elem).attr('color');
            }
        });
    });

    return backup;
}


function backupGroupsWidget() {
    let backup = {
        'groups': [],
        'samples': {}
    };

    $('[checkbox-for="group"]:checked').each((index, element) => {
        backup['groups'].push($(element).attr('data-group'));
    });

    $('[data-sample]:checked').each((index, element) => {
        if (!backup['samples'].hasOwnProperty($(element).attr('data-group'))) {
            backup['samples'][$(element).attr('data-group')] = [];
        }
        backup['samples'][$(element).attr('data-group')].push($(element).attr('data-sample'));
    });

    if (backup['groups'].length > 0 || Object.keys(backup['samples']).length > 0) {
        sample_groups_backup[$('#sample_groups').attr('created-for-category')] = backup;
    }
}

function backupFilters() {
    let backup = serialize_filtering_widgets();
    if (Object.keys(backup).length > 0) {
        if (!filter_backup.hasOwnProperty($('#controls').attr('created-for-gene-caller-id'))) {
            filter_backup[$('#controls').attr('created-for-gene-caller-id')] = {};
        }

        filter_backup[$('#controls').attr('created-for-gene-caller-id')][$('#controls').attr('created-for-engine')] = backup;
    }
}


function serializeState() {
    backupGroupsWidget();
    backupFilters();

    let state = {
        'version': '1',
        'gene_callers_id': $('#gene_callers_id_list').val(),
        'engine': $('[name=engine]:checked').val(),
        'category': $('#sample_groups_list').val(),
        'sample_groups_backup': sample_groups_backup,
        'filter_backup': filter_backup,
        'color_legend': color_legend,
        'size_legend': size_legend,
        'auxiliary': serializeAuxiliaryInputs()
    };

    return state;
}

function showLoadStateWindow()
{
    $.ajax({
        type: 'GET',
        cache: false,
        url: '/state/all',
        success: function(state_list) {
            $('#loadState_list').empty();

            for (let state_name in state_list) {
                $('#loadState_list').append('<option lastmodified="' + state_list[state_name]['last_modified'] + '">' + state_name + '</option>');
            }

            $('#modLoadState').modal('show');
        }
    });
}

function showSaveStateWindow()
{
    $.ajax({
        type: 'GET',
        cache: false,
        url: '/state/all',
        success: function(state_list) {
            $('#saveState_list').empty();

            for (let state_name in state_list) {
                var _select = "";
                if (state_name == current_state_name)
                {
                    _select = ' selected="selected"';
                }
                $('#saveState_list').append('<option ' + _select + '>' + state_name + '</option>');
            }

            $('#modSaveState').modal('show');
            if ($('#saveState_list').val() === null) {
                $('#saveState_name').val('default');
            } else {
                $('#saveState_list').trigger('change');
            }
        }
    });
}

function saveState()
{
    var name = $('#saveState_name').val();

    if (name.length==0) {
        $('#saveState_name').focus();
        return;
    }

    $.ajax({
        type: 'POST',
        cache: false,
        url: '/state/save/' + name,
        data: {
            'content': JSON.stringify(serializeState(), null, 4)
        },
        success: function(response) {
            if (typeof response != 'object') {
                response = JSON.parse(response);
            }

            if (response['status_code']==0)
            {
                toastr.error("Failed, Interface running in read only mode.");
            }
            else if (response['status_code']==1)
            {
                // successfull
                $('#modSaveState').modal('hide');
                toastr.success("State '" + name + "' successfully saved.");
            }
        }
    });
}

function loadState()
{
    $('#modLoadState').modal('hide');
    if ($('#loadState_list').val() == null) {
        return;
    }

    var state_name = $('#loadState_list').val();

    $.ajax({
        type: 'GET',
        cache: false,
        url: '/state/get/' + state_name,
        success: function(response) {
            cached_orientation_matrices = {};
            $('#controls').empty();
            $('#controls').removeAttr('created-for-gene-caller-id');
            $('#controls').removeAttr('created-for-engine')
            $('#sample_groups').empty();
            $('#sample_groups').removeAttr('created-for-category');

            state = JSON.parse(response['content']);
            current_state_name = state_name;

            sample_groups_backup = state['sample_groups_backup'];
            filter_backup = state['filter_backup'];
            color_legend = state['color_legend'];
            size_legend = state['size_legend'];

            if($(`#sample_groups_list option[id='${state['category']}']`).length > 0) {
                $('#sample_groups_list').val(state['category']);
            }

            if($(`[name=engine][value='${state['engine']}']`).length > 0) {
                $(`[name=engine][value='${state['engine']}']`).prop('checked', true);
            }

            if($(`#gene_callers_id_list option[id='${state['gene_callers_id']}']`).length > 0) {
                $('#gene_callers_id_list').val(state['gene_callers_id']);
            }

            for (let tab_name in state['auxiliary']) {
                for (let object_id in state['auxiliary'][tab_name]) {
                    let selector = `#${tab_name} #${object_id}`;

                    if ($(selector).length > 0) {
                        let elem = $(selector)[0];

                        if (elem.tagName == 'SELECT') {
                            $(elem).val(state['auxiliary'][tab_name][object_id]);
                            $(elem).trigger('change');
                        }
                        else if (elem.tagName == 'INPUT') {
                            let type = elem.getAttribute('type');

                            if (type == 'checkbox') {
                                $(elem).prop('checked', state['auxiliary'][tab_name][object_id]);
                            }
                            else if (type == 'input' || type == 'text') {
                                $(elem).val(state['auxiliary'][tab_name][object_id]);
                            }
                        }
                        else if (elem.tagName == 'DIV' && elem.className.indexOf('colorpicker') > -1) {
                            $(elem).attr('color', state['auxiliary'][tab_name][object_id]);
                            $(elem).css('background-color', state['auxiliary'][tab_name][object_id]);
                        }
                    }
                }
            }

            // start the chain reaction...
            $('#gene_callers_id_list').trigger('change');
        }
    });
}

