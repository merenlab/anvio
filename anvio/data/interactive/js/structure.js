var stage;

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
                $('#engine_list').append(`<input type="radio" name="engine" value="${engine}" id="engine_${engine}" ${engine == default_engine ? 'checked="checked"' : ''}><label for="engine_${engine}">${engine}</label>`);
            });
            
            available_sample_ids.forEach(function(sample_id) {

                $('#sample_id_list').append(`<input class="form-check-input" type="checkbox" id="sample_${sample_id}" value="${sample_id}"><label class="form-check-label" for="sample_${sample_id}">${sample_id}</label><br />`);
            });

            $('#sample_id_list').trigger('change');
        }
    });
});


/*        var SAAVs_array = [5, 43, 22, 30];
        var atomPair = [[ "22.CA", "5.CA" ], [ "22.CA", "43.CA" ],  [ "22.CA", "30.CA" ]];

        var SAAVs_string = SAAVs_array.join(", ");

        var schemeId = NGL.ColormakerRegistry.addScheme(function (params) {
          this.atomColor = function (atom) {
            if (atom.resno == SAAVs_array[0]) {
              return 0x0000FF  // blue
            } else if (atom.resno == SAAVs_array[1]) {
              return 0xFF0000  // red
            } else {
              return 0x00FF00  // green
            }
          }
        })*/

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
        url: '/data/get_queried_structure/' + gene_callers_id,
        success: function(data) {
            // create tooltip element and add to document body
            var tooltip = document.createElement("div")
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
            stage.loadFile(stringBlob, { ext: "pdb" }).then(defaultStructureRepresentation).then();
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
    $.ajax({
        type: 'POST',
        cache: false,
        data: {
            engine: $('[name=engine]:checked').val(),
            samples_of_interest: $('#sample_id_list input:checkbox:checked').toArray().map((checkbox) => { return $(checkbox).val(); }),
            departure_from_consensus_min_max: $('#departure_from_consensus').val(),
        },
        url: '/data/get_variability',
        success: function(data) {

        }
    });
};