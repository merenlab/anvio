var stage;
var histogram_data;

$(document).ready(function() {
    stage = new NGL.Stage("viewport");
    stage.setParameters({
        backgroundColor:"white"
    });

    window.addEventListener( "resize", function( event ){
        stage.handleResize();
    }, false );

    $('#gene_callers_id_list').on('change', function(ev) {
        load_protein($('#gene_callers_id_list').val());
    });

    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/get_available_genes_and_samples?timestamp=' + new Date().getTime(),
        success: function(data) {
            let available_gene_callers_ids = data['available_gene_callers_ids'];
            let available_sample_ids = data['available_sample_ids'];
            let available_engines = data['available_engines']; 

            available_gene_callers_ids.forEach(function(gene_callers_id) {
                $('#gene_callers_id_list').append(`<option id=${gene_callers_id}>${gene_callers_id}</option>`);
            });

            $('#gene_callers_id_list').trigger('change');

            let default_engine = available_engines[0];
            available_engines.forEach(function(engine) {
                $('#engine_list').append(`<input type="radio" name="engine" onclick="draw_histogram();" value="${engine}" id="engine_${engine}" ${engine == default_engine ? 'checked="checked"' : ''}><label for="engine_${engine}">${engine}</label>`);
            });
                        
            available_sample_ids.forEach(function(sample_id) {

                $('#sample_id_list').append(`<input class="form-check-input" type="checkbox" id="sample_${sample_id}" value="${sample_id}" checked="checked"><label class="form-check-label" for="sample_${sample_id}">${sample_id}</label><br />`);
            });
            $('#sample_id_list').trigger('change');
        }
    });
});


function defaultStructureRepresentation( component ){
    // bail out if the component does not contain a structure
    if( component.type !== "structure" ) return;

    component.addRepresentation( "cartoon", {
        aspectRatio: 3.0,
        scale: 1.5
    } );

    // add annotation to a protein chain
    var chainText = {
    "A": "how do i put gene_caller_id here?",
    }
    var ap = component.structure.getAtomProxy()
    component.structure.eachChain(function (cp) {
        // annotation is anchored to the residue index equal to half the total number of residues
        ap.index = cp.atomOffset + Math.floor(cp.atomCount / 2)
        component.addAnnotation(ap.positionToVector3(), chainText[ cp.chainname ])
    }, new NGL.Selection("polymer"));

    var pa = component.structure.getPrincipalAxes();
    stage.animationControls.rotate(pa.getRotationQuaternion(), 1000);
    stage.autoView();        
}


function load_protein(gene_callers_id) {
    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/get_structure/' + gene_callers_id,
        success: function(data) {
            histogram_data = data['histograms'];
            draw_histogram();

            // create tooltip element and add to document body
            var tooltip = document.createElement("div");
            Object.assign(tooltip.style, {
                display: "none",
                position: "fixed",
                zIndex: 10,
                pointerEvents: "none",
                backgroundColor: "rgba( 0, 0, 0, 0.6 )",
                color: "lightgrey",
                padding: "8px",
                fontFamily: "sans-serif"
            })
            document.body.appendChild(tooltip)

            stage.removeAllComponents();

            var stringBlob = new Blob( [ data['pdb_content'] ], { type: 'text/plain'} );
            stage.loadFile(stringBlob, { ext: "pdb", name: gene_callers_id }).then(defaultStructureRepresentation).then();
            // remove default hoverPick mouse action
            stage.mouseControls.remove("hoverPick")
            // listen to `hovered` signal to move tooltip around and change its text
            stage.signals.hovered.add(function (pickingProxy) {
                if (pickingProxy && (pickingProxy.atom || pickingProxy.bond)) {
                    var atom = pickingProxy.atom || pickingProxy.closestBondAtom
                    var mp = pickingProxy.mouse.position
                    tooltip.innerText = "SAAV: \n" + atom.qualifiedName() + "\nAla:\t19" + "\nIle:\t5"
                    tooltip.style.bottom = window.innerHeight - mp.y + 3 + "px"
                    tooltip.style.left = mp.x + 3 + "px"
                    tooltip.style.display = "block"
                } else {
                tooltip.style.display = "none"
              }
            });
        }
    });
}

function draw_variability() {
    let gene_callers_id = $('#gene_callers_id_list').val()
    $.ajax({
        type: 'POST',
        cache: false,
        data: {
            'gene_callers_id': gene_callers_id,
            'engine': $('[name=engine]:checked').val(),
            'selected_samples': $('#sample_id_list input:checkbox:checked').toArray().map((checkbox) => { return $(checkbox).val(); }),
            'departure_from_consensus': $('#departure_from_consensus').val(),
        },
        url: '/data/get_variability',
        success: function(data) {
            let component = stage.compList[0];
            let variant_residues = [];

            console.log(data);

            for (let index in data) {
                variant_residues.push(data[index]['codon_order_in_gene']);
            }

            component.reprList.forEach((rep) => {
                if (rep.name == 'spacefill') {
                    rep.dispose();
                }
            });

            component.addRepresentation("spacefill", {
                sele: "(" + variant_residues.join(', ') + ") and .CA",
                scale: 1
            });

        },
        error: function(request, status, error) {
            console.log(request, status, error);
        }
    });
};


function draw_histogram() {
    let engine = $('[name=engine]:checked').val();

    for (let column in histogram_data[engine]) {
        let canvas = document.getElementById('histogram_' + column);
        var ctx = canvas.getContext("2d");

        let width = parseFloat(canvas.width);
        let height = parseFloat(canvas.height);

        // clear canvas
        ctx.clearRect(0, 0, width, height);

        let bins = histogram_data[engine][column]['bins'];
        let counts = histogram_data[engine][column]['counts'];
        let curve_x = histogram_data[engine][column]['curve_x'];
        let curve_y = histogram_data[engine][column]['curve_y'];

        // normalize counts and curve_y to pixels so they fit in canvas frame
        var max_count = Math.max(...counts);
        var max_curve_y = Math.max(...curve_y);
        counts = counts.map(v => (v / max_count) * height);
        curve_y = curve_y.map(v => (v / max_curve_y) * height);

        // normalize bins and curve_x to pixels so they fit in canvas frame
        var max_slider = parseFloat(document.getElementById(column).dataset.sliderMax);
        var min_slider = parseFloat(document.getElementById(column).dataset.sliderMin);
        bins = bins.map(v => (v / max_slider) * width);
        curve_x = curve_x.map(v => (v / max_slider) * width);

        // draw histogram
        for (let i=0; i < bins.length - 1; i++) {
            ctx.fillRect(bins[i],              // x
                         height - counts[i],   // y
                         bins[i + 1] - bins[i],// width
                         counts[i]);           // height
            ctx.stroke();
        }

        // draw curve
        ctx.lineWidth=1.0;
        ctx.beginPath();
        for (let j=0; j < curve_x.length; j++) {
            ctx.moveTo(curve_x[j],   height-curve_y[j]);
            ctx.lineTo(curve_x[j+1], height-curve_y[j+1]);
            ctx.stroke();
        }

        // draw points outside range over which curve is defined
        var min_curve_x = Math.min(...curve_x)
        if (Math.min(...curve_x) > 0) {
            ctx.moveTo(min_curve_x, height - curve_y[0]);
            ctx.lineTo(min_curve_x, height)
            ctx.moveTo(min_curve_x, height);
            ctx.lineTo(width, height)
            ctx.stroke();
        }

        var max_curve_x = Math.max(...curve_x)
        if (Math.max(...curve_x) < width) {
            ctx.moveTo(max_curve_x, height - curve_y[curve_y.length - 1]);
            ctx.lineTo(max_curve_x, height)
            ctx.moveTo(max_curve_x, height);
            ctx.lineTo(width, height)
            ctx.stroke();
        }
    }
};
