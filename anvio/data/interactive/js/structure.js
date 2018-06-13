var stages = {};
var variability = {};
var histogram_data;
var sample_groups;
var pdb_content;

var color_legend = {};
var size_legend = {};

$(document).ready(function() {
    $('.colorpicker').colpick({
        layout: 'hex',
        submit: 0,
        colorScheme: 'light',
        onChange: function(hsb, hex, rgb, el, bySetColor) {
            $(el).css('background-color', '#' + hex);
            $(el).attr('color', '#' + hex);

            if (!bySetColor) $(el).val(hex);
        },
        onHide: function() {
            draw_variability();
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
        $.when({}).then(load_protein).then(() => {
            create_ui();
            load_sample_group_widget($('#sample_groups_list').val());
        });
    });

    $('#sample_groups_list').on('change', function(ev) {
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

            $.when({}).then(load_protein).then(() => {
                let default_engine = available_engines[0];
                available_engines.forEach(function(engine) {
                    $('#engine_list').append(`<input type="radio" name="engine" onclick="create_ui();" value="${engine}" id="engine_${engine}" ${engine == default_engine ? 'checked="checked"' : ''}><label for="engine_${engine}">${engine}</label>`);
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


function load_sample_group_widget(category) {
    $('#sample_groups').empty();
    tableHtml = '<table class="table table-condensed"><tr><td>Groups</td><td>Samples</td></tr>';

    for (let group in sample_groups[category]) {
        tableHtml += `
            <tr>
                <td>
                    <input class="form-check-input" 
                        onclick="create_ngl_views();"
                        checkbox-for="group"
                        id="${category}_${group}"
                        type="checkbox" 
                        data-category="${category}"
                        data-group="${group}"
                        value="${group}" 
                        checked="checked">
                    <label class="form-check-label" for="${category}_${group}">${group}</label>
                </td>
                <td>`;

        sample_groups[category][group].forEach((sample) => {
            tableHtml += `
                <div class="table-group-checkbox" style="display: inline-block; float: left;">
                    <input class="form-check-input" 
                            id="${category}_${group}_${sample}"
                            onclick="draw_variability();"
                            type="checkbox" 
                            data-category="${category}"
                            data-group="${group}"
                            data-sample="${sample}"
                            value="${sample}" 
                            checked="checked">
                    <label class="form-check-label" for="${category}_${group}_${sample}">${sample}</label>
                </div>`;
        });

        tableHtml += '</td></tr>';
    }

    $('#sample_groups').append(tableHtml + '</table>');
    create_ngl_views();
}

function apply_orientation_matrix_to_all_stages(orientationMatrix) {
    for (let group in stages) {
        stages[group].viewerControls.orient(orientationMatrix); 
    }
}

function create_ngl_views() {
    for (let group in stages) {
        stages[group].dispose();
    }
    stages = {};

    $('#ngl-container').empty();

    let selected_groups = $('[checkbox-for="group"]:checked');

    $(selected_groups).each((index, element) => {
        let group = $(element).attr('data-group');
        let num_cells = selected_groups.length;
        let num_columns = Math.min(4, Math.ceil(Math.sqrt(num_cells)));
        let num_rows = Math.min(4, Math.ceil(num_cells / num_columns));

        $('#ngl-container').append(`
            <div id="ngl_${group}_wrapper" 
                 class="col-md-${parseInt(12 / num_columns)} nopadding" 
                 style="height: ${parseFloat(100 / num_rows)}%; ">
                 <div class="ngl-group-title">
                    ${group}
                 </div>
                 <div class="ngl-group-fullscreen">
                    <button type="button" class="btn btn-link btn-sm" onclick="stages['${group}'].toggleFullscreen();" title="Fullscreen">
                        <span class="glyphicon glyphicon-fullscreen"></span>
                     </button>
                 </div>
                 <div id="ngl_${group}" class="ngl-inner">

                 </div> 
            </div>`);

        var stage = new NGL.Stage(`ngl_${group}`);
        var stringBlob = new Blob( [ pdb_content ], { type: 'text/plain'} );
        
        stage.loadFile(stringBlob, { ext: "pdb" })
             .then((component) => {
                if( component.type !== "structure" ) return;

                // component.addRepresentation("surface", {
                //     surfaceType: "ms",
                //     smooth: 2,
                //     probeRadius: 1.4,
                //     scaleFactor: 2.0,
                //     flatShaded: false,
                //     opacity: 0.8,
                //     lowResolution: false,
                //     colorScheme: "element"
                // });

                component.addRepresentation("cartoon", {
                    //colorScheme: 'residueindex',
                    metalness: 0.1,
                    aspectRatio: 3.0,
                    scale: 1.5
                });

                var pa = component.structure.getPrincipalAxes();     
                component.setRotation(pa.getRotationQuaternion());
                stage.autoView();
             });
        
        stage.setParameters({
            backgroundColor: "white"
        });

        // prevent default tooltip
        stage.mouseControls.remove("hoverPick");

        // add custom tooltip
        stage.signals.hovered.add(function (pickingProxy) {
            let tooltip = document.getElementById('ngl-tooltip');

            if (pickingProxy && pickingProxy.atom) {
                let residue = pickingProxy.atom.resno;
                let mp = pickingProxy.mouse.position;

                if (variability[group].hasOwnProperty(residue)) {
                    let HTML_reference_title = `<h4>Reference info</h4>`
                    let HTML_reference_body = `
                        <tr><td>Secondary Structure</td><td>${variability[group][residue]['sec_struct']}</td></tr>
                        <tr><td>Solvent Accessibility</td><td>${variability[group][residue]['rel_solvent_acc'].toFixed(2)}</td></tr>
                        <tr><td>(Phi, Psi)</td><td>(${variability[group][residue]['phi'].toFixed(1)}, ${variability[group][residue]['psi'].toFixed(1)})</td></tr>
                        `

                    let HTML_variant_title = `<h4>Variant info</h4>`
                    let HTML_variant_body = `
                        <tr><td>Consensus</td><td>${variability[group][residue]['consensus']}</td></tr>
                        <tr><td>Mean Dfc</td><td>${variability[group][residue]['departure_from_consensus'].toFixed(2)}</td></tr>
                        <tr><td>Prevalence</td><td>${variability[group][residue]['occurrence']} of ${parseInt(Math.round(variability[group][residue]['occurrence'] / variability[group][residue]['prevalence']))} samples</td></tr>
                        <tr><td>Mean Coverage</td><td>${variability[group][residue]['coverage'].toFixed(2)}</td></tr>
                        <tr><td>Mean Entropy</td><td>${variability[group][residue]['entropy'].toFixed(2)}</td></tr>
                        `

                    // add CDN specific columns
                    if ($('[name=engine]:checked').val() == 'CDN') {
                        // prepend to body
                        HTML_reference_body = `<tr><td>Codon</td><td>${variability[group][residue]['reference']}</td></tr>` + HTML_reference_body
                        HTML_reference_body = `<tr><td>Codon No.</td><td>${variability[group][residue]['codon_number']}</td></tr>` + HTML_reference_body
                        // append to body
                        HTML_variant_body += `<tr><td>Synonymity</td><td>${variability[group][residue]['synonymity'].toFixed(2)}</td></tr>`
                    }

                    // add AA specific columns
                    if ($('[name=engine]:checked').val() == 'AA') {
                        // prepend to body
                        HTML_reference_body = `<tr><td>Residue</td><td>${variability[group][residue]['reference']}</td></tr>` + HTML_reference_body
                        HTML_reference_body = `<tr><td>Residue No.</td><td>${variability[group][residue]['codon_number']}</td></tr>` + HTML_reference_body
                        // append to body
                        HTML_variant_body += `<tr><td>Mean BLOSUM90</td><td>${variability[group][residue]['BLOSUM90'].toFixed(1)}</td></tr>`
                    }

                    HTML_reference_body = `<table class="table table-striped">` + HTML_reference_body + `</table>`
                    HTML_variant_body = `<table class="table table-striped">` + HTML_variant_body + `</table>`

                    tooltip.innerHTML = HTML_reference_title + HTML_reference_body + HTML_variant_title + HTML_variant_body;
                    tooltip.style.bottom = window.innerHeight - mp.y + 3 + "px";
                    tooltip.style.left = mp.x + 3 + "px";
                    tooltip.style.display = "block";
                }
            } else {
                tooltip.style.display = "none";
            }
        });


        let func = () => { apply_orientation_matrix_to_all_stages( stage.viewerControls.getOrientation()); };
        stage.mouseObserver.signals.scrolled.add(() => { func(); });
        stage.mouseObserver.signals.dragged.add(() => { func(); });
        $(`#ngl_${group}`).mouseup(func);

        stages[group] = stage;
    });

    draw_variability();
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
            defer.resolve();
        }
    });

    return defer.promise();
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

function draw_variability() {
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
            variability = {};
            for (let group in response_all) {
                let response = response_all[group];
                
                let data = JSON.parse(response['data']);
                let total_entries = response['total_entries'];
                let entries_after_filtering = response['entries_after_filtering'];

                let component = stages[group].compList[0];

                component.reprList.slice(0).forEach((rep) => {
                    if (rep.name == 'spacefill') {
                        rep.dispose();
                    }
                });

                let codon_to_variability = {};
                for (let index in data) {
                    codon_to_variability[data[index]['codon_number']] = data[index];
                }

                variability[group] = codon_to_variability;

                if (Object.keys(data).length > 0) {
                    for (let index in data) {
                        let spacefill_options = {
                            sele: data[index]['codon_number'] + " and .CA",
                            scale: 1
                        }

                        if ($('#color_type').val() == 'Dynamic') {
                            let column = $('#color_target_column').val();
                            let widget = $('.widget[data-column="' + column + '"]');
                            let controller = $(widget).attr('data-controller');
                            let column_value = data[index][column];

                            if (controller == 'slider') {
                                let min_value = parseFloat($('#color_min').val());
                                let max_value = parseFloat($('#color_max').val());

                                if (min_value >= max_value) {
                                    $('#dynamic_color_error').show();
                                } else {
                                    $('#dynamic_color_error').hide();
                                }

                                let val = (parseFloat(column_value) - min_value) / (max_value - min_value);
                                
                                val = Math.max(0, Math.min(1, val));

                                spacefill_options['color'] = getGradientColor(
                                    $('#color_start').attr('color'),
                                    $('#color_end').attr('color'),
                                    val);
                            }
                            else
                            {
                                spacefill_options['color'] = color_legend[engine][column][column_value];   
                            }
                        } else {
                            spacefill_options['color'] = $('#color_static').attr('color');
                        }

                        if ($('#size_type').val() == 'Dynamic') {
                            let column = $('#size_target_column').val();
                            let widget = $('.widget[data-column="' + column + '"]');
                            let controller = $(widget).attr('data-controller');
                            let column_value = data[index][column];

                            if (controller == 'slider') {
                                let min_value = parseFloat($('#size_min').val());
                                let max_value = parseFloat($('#size_max').val());
                                let start_value = parseFloat($('#size_start').val());
                                let end_value = parseFloat($('#size_end').val());

                                if (min_value >= max_value) {
                                    $('#dynamic_size_error').show();
                                } else {
                                    $('#dynamic_size_error').hide();
                                }

                                let val = (parseFloat(column_value) - min_value) / max_value - min_value;
                                
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
                    }
                }
            }

            $('.overlay').hide();
        },
        error: function(request, status, error) {
            console.log(request, status, error);
            $('.overlay').hide();
        }
    });
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

        for (let i=0; i < normalized_bins.length; i++) {
            data_points.push({'x': normalized_bins[i], 'y': height - normalized_counts[i]});
        }

        var interpolate = d3.line()
                            .x(function(d) { return d.x; })
                            .y(function(d) { return d.y; })
                            .curve(d3.curveCardinal);

        var interpolated_line = interpolate(data_points);
        interpolated_line += `L ${data_points[normalized_bins.length - 1]['x']} ${height} L ${data_points[0]['x']} ${height}`;

        svg.append("path")
            .style("fill","#337ab7")
            .attr("d",function(d,i){ return interpolated_line; });
    }
};

function create_ui() {
    let gene_callers_id = $('#gene_callers_id_list').val();
    let engine = $('[name=engine]:checked').val();

   $.ajax({
        type: 'POST',
        cache: false,
        url: '/data/get_column_info',
        data: {
            'gene_callers_id': gene_callers_id,
            'engine': engine
        },
        success: function(data) {
            let container = $('#controls');

            // remove widgets
            container.empty();

            $('#color_target_column').empty();
            $('#size_target_column').empty();

            if (!color_legend.hasOwnProperty(engine)) {
                color_legend[engine] = {};
            }

            if (!size_legend.hasOwnProperty(engine)) {
                size_legend[engine] = {};
            }

            data.forEach((item) => {
                $('#color_target_column').append(`<option value="${item['name']}">${item['title']}</item>`);
                $('#size_target_column').append(`<option value="${item['name']}">${item['title']}</item>`);

                if (item['controller'] == 'slider') {
                    $(container).append(`
                        <div class="widget" data-column="${item['name']}" data-controller="${item['controller']}">
                            ${item['title']}<br />
                            <svg id="histogram_${item['name']}" width="100%" height="30" style="position: relative; top: 6;" viewBox="0 0 200 30" preserveAspectRatio="none"></svg>   
                            <input id="${item['name']}" 
                                    type="${item['data_type']}" 
                                    data-provide="slider"
                                    data-slider-min="${item['min']}" 
                                    data-slider-max="${item['max']}" 
                                    data-slider-step="${item['step']}" 
                                    data-slider-value="[${item['min']},${item['max']}]"
                                    >
                        </div>
                    `);
                    $(`#${item['name']}`).slider({}).on('slideStop', () => { draw_variability(); });
                }
                if (item['controller'] == 'checkbox') {
                    $(container).append(`
                        <div class="widget" data-column="${item['name']}" data-controller="${item['controller']}">
                            ${item['title']}<br />
                            ${item['choices'].map((choice) => { return `
                                <input class="form-check-input" type="checkbox" id="${item['name']}_${choice}" value="${choice}" onclick="draw_variability();" checked="checked">
                                <label class="form-check-label" for="${item['name']}_${choice}">${choice}</label>`; }).join('')}
                            <br />
                            <button class="btn btn-xs" onclick="$(this).closest('.widget').find('input:checkbox').prop('checked', true); draw_variability();">Check All</button>
                            <button class="btn btn-xs" onclick="$(this).closest('.widget').find('input:checkbox').prop('checked', false); draw_variability();">Uncheck All</button>
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
        }
    });   
}


function onTargetColumnChange(element) {
    // this on change event shared between color_target_column, size_target_column.
    let engine = $('[name=engine]:checked').val();
    let column = $(element).val();
    let widget = $('.widget[data-column="' + column + '"]');
    let controller = $(widget).attr('data-controller');

    // color or size
    let prefix = element.getAttribute('id').split('_')[0];

    // show column related panel.
    // for linear values, show slider panel
    // for discreete values, show legend panel
    // read the min/max from slider and put into prefixed input in perspective
    if (controller == 'slider') {
        $(`#${prefix}_slider_panel`).show();
        $(`#${prefix}_legend_panel`).hide();
        let slider = $(widget).find('input');

        $(`#${prefix}_min`).val($(slider).attr('data-slider-min'));
        $(`#${prefix}_max`).val($(slider).attr('data-slider-max'));
    } 
    else 
    {
        $(`#${prefix}_slider_panel`).hide();
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

        // hide all representations besides we want.
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


async function generate_summary() {
    let serialized_groups = serialize_checked_groups();
    var zip = new JSZip();

    for (let group in stages) {
        zip.file(`images/${group}/merged.png`, await make_image(group));

        // generate per sample.
        for (let i=0; i < serialized_groups[group].length; i++) {
            let sample_id = serialized_groups[group][i];
            zip.file(`images/${group}/${sample_id}.png`, await make_image(group, sample_id));
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
        saveAs(content, $('#zip_name').val());
    });
}
