//  Edit Attributes For Multiple Layers

$(document).ready(function() {
    $('.select_layer').on('change', function() {
        var table = $(this).closest('table');
        var layer_name = this.value;
        // clean prior selections
        $(table).find('.layer_selectors').prop('checked', false);

        if(layer_name){ // if layer_name is empty, there is nothing to select, move on.
            $(table).find('.titles').each(
                function(){
                    if (this.title.toLowerCase().indexOf(layer_name.toLowerCase()) > -1)
                    {
                        $(this).parent().find('.layer_selectors').prop('checked','checked');
                    }
                }
            );
        }
    });

    $('.picker_multiple, .picker_start_multiple').colpick({
        layout: 'hex',
        submit: 0,
        colorScheme: 'light',
        onChange: function(hsb, hex, rgb, el, bySetColor) {
            $(el).css('background-color', '#' + hex);
            $(el).attr('color', '#' + hex);
            if (!bySetColor) $(el).val(hex);

            var table = $(el).closest('table');
            $(table).find('.layer_selectors:checked').each(
                function(){
                    var selector = '.colorpicker';

                    if ($(el).hasClass("picker_start_multiple"))
                    {
                        selector += ':first';
                    }
                    else
                    {
                        selector += ':last';
                    }
                    var picker = $(this).parent().parent().find(selector);

                    $(picker).attr('color', '#' + hex);
                    $(picker).css('background-color', '#' + hex);
                }
            );
        }
    }).keyup(function() {
            $(this).colpickSetColor(this.value);
    });


    $('.input-height-multiple, .input-margin-multiple, .input-min-multiple, .input-max-multiple').on('change', function() {
        var new_val = this.value;
        var target_selector = '.' + this.getAttribute('class').replace('-multiple', '') + ':enabled';
        var table = $(this).closest('table');

        $(table).find('.layer_selectors:checked').each(
            function(){
                var row = $(this).parent().parent();
                $(row).find(target_selector).val(new_val);
            }
        );
    });

    $('.normalization_multiple, .type_multiple').on('change', function() {
        var new_val = this.value;
        var target_selector = '.' + this.getAttribute('class').replace('_multiple', '');
        var table = $(this).closest('table');

        $(table).find('.layer_selectors:checked').each(
            function(){
                var row = $(this).parent().parent();
                $(row).find(target_selector).val(new_val).trigger('change');
            }
        );
    });

    $('.select_all').on('click', function() {
        var new_val = this.checked;
        var table = $(this).closest('table');

        $(table).find('.layer_selectors').each(
            function() {
                this.checked = new_val;
            }
        );
    });
});