//  Edit Attributes For Multiple Layers

// get numeric part from id
function getNumericPart(id){
    var $num = id.replace(/[^\d]+/, '');

    return $num;
}

$(document).ready(function() {
    $('#select_layer').on('change', function() {
        var layer_name = $('#select_layer').val();

        // clean prior selections
        $('.layer_selectors').prop('checked', false);

        if(layer_name){ // if layer_name is empty, there is nothing to select, move on.
            $('.titles').each(
                function(){
                    if (this.title.indexOf(layer_name) > -1)
                    {
                        $('#select_this_' + getNumericPart(this.id)).prop('checked','checked');
                    }
                }
            );
        }
    });

    $('#picker_multiple').colpick({
        layout: 'hex',
        submit: 0,
        colorScheme: 'light',
        onChange: function(hsb, hex, rgb, el, bySetColor) {
            $(el).css('background-color', '#' + hex);
            $(el).attr('color', '#' + hex);
            if (!bySetColor) $(el).val(hex);

            $('.layer_selectors:checked').each(
                function(){
                    $('#' + 'picker' + getNumericPart(this.id)).attr('color', '#' + hex);
                    $('#' + 'picker' + getNumericPart(this.id)).css('background-color', '#' + hex);
                }
            );
        }
    }).keyup(function() {
            $(this).colpickSetColor(this.value);
    });

    $('#picker_start_multiple').colpick({
        layout: 'hex',
        submit: 0,
        colorScheme: 'light',
        onChange: function(hsb, hex, rgb, el, bySetColor) {
            $(el).css('background-color', '#' + hex);
            $(el).attr('color', '#' + hex);
            if (!bySetColor) $(el).val(hex);

            $('.layer_selectors:checked').each(
                function(){
                    $('#' + 'picker_start' + getNumericPart(this.id)).attr('color', '#' + hex);
                    $('#' + 'picker_start' + getNumericPart(this.id)).css('background-color', '#' + hex);
                }
            );
        }
    }).keyup(function() {
            $(this).colpickSetColor(this.value);
    });

    $('#min_multiple').on('change', function(){
        var intend_value = $('#min_multiple').val();
        $('.layer_selectors:checked').each(
            function(){
                $('#' + 'min' + getNumericPart(this.id)).attr('value', intend_value);
            }
        );
    });

    $('#max_multiple').on('change', function(){
        var intend_value = $('#max_multiple').val();
        $('.layer_selectors:checked').each(
            function(){
                $('#' + 'max' + getNumericPart(this.id)).attr('value', intend_value);
            }
        );
    });

    $('#height_multiple').on('change', function(){
        var intend_value = $('#height_multiple').val();
        $('.layer_selectors:checked').each(
            function(){
                $('#' + 'height' + getNumericPart(this.id)).attr('value', intend_value);
            }
        );
    });

    $('#margin_multiple').on('change', function(){
        var intend_value = $('#margin_multiple').val();
        $('.layer_selectors:checked').each(
            function(){
                $('#' + 'margin' + getNumericPart(this.id)).attr('value', intend_value);
            }
        );
    });

    $('#type_multiple').on('change', function(){
        var intend_value = $('#type_multiple').val();
        $('.layer_selectors:checked').each(
            function(){
                $('#type' + getNumericPart(this.id)).val(intend_value).trigger('change');
            }
        );
    });

    $('#normalization_multiple').on('change', function(){
        var intend_value = $('#normalization_multiple option:selected').val();

        $('.layer_selectors:checked').each(
            function(){
                $('#normalization' + getNumericPart(this.id)).val(intend_value).trigger('change');
            }
        );
    });

    $('#select_all').on('click', function() {
        if(this.checked) {
            $('.layer_selectors').each(
                function() {
                    this.checked = true;
                }
            );

        }else{
            $('.layer_selectors').each(
                function() {
                    this.checked = false;
                }
            );
        }
    });
});