var stages = {};
var variability = {};
var histogram_data;
var sample_groups;
var pdb_content;

$(document).ready(function() {
    $('.colorpicker').colpick({
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

                component.addRepresentation( "cartoon", {
                    aspectRatio: 3.0,
                    scale: 1.5
                } );

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
                    tooltip.innerHTML = `<table>
                            <tr><td>Coverage</td><td>${variability[group][residue]['coverage']}</td></tr>
                            <tr><td>Cove</td><td>${variability[group][residue]['coverage']}</td></tr>
                            <tr><td>Coverage</td><td>${variability[group][residue]['coverage']}</td></tr>
                            </table>`;
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
    
    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/get_structure/' + gene_callers_id,
        success: function(data) {
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

                variability[group] = data;

                if (Object.keys(data).length > 0) {
                    for (let index in data) {
                        let spacefill_options = {
                            sele: data[index]['codon_order_in_gene'] + " and .CA",
                            scale: 1
                        }

                        if ($('#color_type').val() == 'Dynamic') {
                            let column = $('#color_target_column').val();
                            let min_value = parseFloat($('#color_min').val());
                            let max_value = parseFloat($('#color_max').val());
                            let val = Math.abs(parseFloat(data[index][column]) - min_value) / Math.abs(max_value - min_value);
                            
                            val = Math.max(0, Math.min(1, val));

                            spacefill_options['color'] = getGradientColor(
                                $('#color_start').attr('color'),
                                $('#color_end').attr('color'),
                                val);
                        } else {
                            spacefill_options['color'] = $('#color_static').attr('color');
                        }

                        if ($('#size_type').val() == 'Dynamic') {
                            let column = $('#size_target_column').val();
                            let min_value = parseFloat($('#size_min').val());
                            let max_value = parseFloat($('#size_max').val());
                            let val = Math.abs(parseFloat(data[index][column]) - min_value) / Math.abs(max_value - min_value);
                            
                            val = Math.max(0, Math.min(1, val));

                            spacefill_options['scale'] = parseFloat($('#size_start').val()) + (val * Math.abs(parseFloat($('#size_end').val()) - parseFloat($('#size_start').val())));
                        } else {
                            spacefill_options['scale'] = parseFloat($('#size_static').val());
                        }

                        component.addRepresentation("spacefill", spacefill_options);
                    }
                }
            }
        },
        error: function(request, status, error) {
            console.log(request, status, error);
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
                }
            });

            draw_histogram();
            draw_variability();
        }
    });   
}


function onTargetColumnChange(element) {
    // this on change event shared between color_target_column, size_target_column.

    let column = $(element).val();
    let widget = $('.widget[data-column="' + column + '"]');
    let controller = $(widget).attr('data-controller');

    // color or size
    let prefix = element.getAttribute('id').split('_')[0];

    // read the min/max from slider and put into prefixed input in perspective
    if (controller == 'slider') {
        let slider = $(widget).find('input');

        $(`#${prefix}_min`).val($(slider).attr('data-slider-min'));
        $(`#${prefix}_max`).val($(slider).attr('data-slider-max'));
    }
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
