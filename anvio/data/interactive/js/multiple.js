/**
 * Edit Attributes For Multiple Layers
 *
 *  Authors: Ozcan Esen
 *           Dogan Can Kilment
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

$(document).ready(function() {
    $('.select_layer').on('change', function() {
        var table = $(this).closest('table');
        var layer_name = this.value;
        // clean prior selections
        $(table).find('.layer_selectors:visible').prop('checked', false);

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

    $('.select_layer, .input-height-multiple, .input-min-multiple, .input-max-multiple').on('keydown', function(e) {
        var keyCode = e.keyCode || e.which;
        if (keyCode == '13') {
            $(e.target).trigger('change');
        }
    });

    $('#picker_multiple, #picker_start_multiple, #picker_multiple_samples, #picker_start_multiple_samples').colpick({
        layout: 'hex',
        submit: 0,
        colorScheme: 'light',
        onChange: function(hsb, hex, rgb, el, bySetColor) {
            $(el).css('background-color', '#' + hex);
            $(el).attr('color', '#' + hex);
            if (!bySetColor) $(el).val(hex);

            var table = $(el).closest('table');
            $(table).find('.layer_selectors:checked:visible').each(
                function(){
                    var selector = '.colorpicker';

                    if ($(el).attr('id').startsWith('picker_start'))
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
        var target_selector = '.' + this.getAttribute('class').split(' ').pop().replace('-multiple', '') + ':enabled';
        var table = $(this).closest('table');

        $(table).find('.layer_selectors:checked:visible').each(
            function(){
                var row = $(this).parent().parent();
                $(row).find(target_selector).val(new_val).trigger('change');
            }
        );
    });

    $('.normalization_multiple, .type_multiple').on('change', function() {
        var new_val = this.value;
        if (new_val == "")
            return;
        // get last item of class name, which is the target selector
        // Make sure to remove _multiple from the class name
        var target_selector = '.' + this.getAttribute('class').split(' ').pop().replace('_multiple', '');
        var table = $(this).closest('table');

        $(table).find('.layer_selectors:checked:visible').each(
            function(){
                var row = $(this).parent().parent();
                var combo = $(row).find(target_selector);

                if (combo.find('option[value="' + new_val + '"]').length > 0)
                {
                    combo.val(new_val).trigger('change');
                    return;
                }
            }
        );
    });

    $('.select_all').on('click', function() {
        var new_val = this.checked;
        var table = $(this).closest('table');

        $(table).find('.layer_selectors:visible').each(
            function() {
                this.checked = new_val;
            }
        );
    });
});

$(document).ready(function() {
    $('.select_bins').on('input', function() {
        var input_value = this.value.toLowerCase();
        var table = $('#bins-table');
        var colorPicker = $('#picker_multiple-bins');

        $('#picker_multiple-bins').css('background-color', '#FFFFFF').attr('color', '#FFFFFF');
        let matchCount = 0;

        table.find('.bin-row').each(function() {
            var bin_name = $(this).find('.bin-name').val().toLowerCase();
            var bin_color = $(this).find('.colorpicker').attr('color');

            if (bin_name.includes(input_value)) {
                $(this).show();
                matchCount++;
                colorPicker.css('background-color', bin_color).attr('color', bin_color);
            } else {
                $(this).hide();
            }
        });

        updateMatchCount(matchCount);

    });

    $('.select_bins').on('keypress', function(e) {
        if (e.which === 13) {
            applyColorToBins();
        }
    });

    $('#apply_color_button').on('click', function() {
        applyColorToBins();
    });

    function applyColorToBins() {
        var input_value = $('.select_bins').val().toLowerCase();
        var selectedColor = $('#picker_multiple-bins').attr('color');
        let matchCount = 0;

        $('#bins-table').find('.bin-row').each(function() {
            var bin_name = $(this).find('.bin-name').val().toLowerCase();
            if (bin_name.includes(input_value)) {
                $(this).find('.colorpicker').css('background-color', selectedColor).attr('color', selectedColor);
                matchCount++;
            }
        });

        updateMatchCount(matchCount);

    }

    function updateMatchCount(count) {
        $('#apply_color_button').text('Apply (' + count + ' matches)');
    }
});