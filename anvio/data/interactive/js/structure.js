var stages = {};
var histogram_data;
var sample_groups;
var pdb_content;

$(document).ready(function() {
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
            output[column] = $(widget).find('input').val();
        }
        else if (controller == 'checkbox') {
            output[column] = $(widget).find('input:checkbox:checked').toArray().map((checkbox) => { return $(checkbox).val(); });
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
                let variant_residues = [];

                for (let index in data) {
                    variant_residues.push(data[index]['codon_order_in_gene']);
                }

                component.reprList.forEach((rep) => {
                    if (rep.name == 'spacefill') {
                        rep.dispose();
                    }
                });

                if (variant_residues.length > 0) {
                    component.addRepresentation("spacefill", {
                        sele: "(" + variant_residues.join(', ') + ") and .CA",
                        scale: 1
                    });

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

            container.empty();

            data.forEach((item) => {
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
                        </div>
                    `);
                }
            });

            draw_histogram();
            draw_variability();
        }
    });   
}
