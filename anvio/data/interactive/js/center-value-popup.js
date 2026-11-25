/**
 * A small popup that appears when the user clicks on a 'Center' input box in the layers table,
 * to help them choose a center value for diverging layers. It offers an input box to enter a
 * value manually, and buttons to set the center to the mean or the median of the layer's values
 * in the current view (after the normalization selected for the layer, if any, is applied).
 *
 * The popup works for the center box of a single layer, as well as for the center box in the
 * 'edit attributes for multiple layers' row, in which case every selected numerical layer gets
 * its own mean or median.
 *
 *  Copyright 2026, The anvi'o project (http://anvio.org)
 *
 *  Anvi'o is a free software. You can redistribute this program
 *  and/or modify it under the terms of the GNU General Public
 *  License as published by the Free Software Foundation, either
 *  version 3 of the License, or (at your option) any later version.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with anvi'o. If not, see <http://opensource.org/licenses/GPL-3.0>.
 *
 *  @license GPL-3.0+ <http://opensource.org/licenses/GPL-3.0>
 */


var CenterValuePopup = (function() {
    var $popup = null;
    var $anchor = null;

    /**
     * Returns the rows of the layers table the popup should act on: the row of the anchor if it is
     * the center box of a single layer, or the rows of all selected numerical layers if it is the
     * center box of the 'edit attributes for multiple layers' row
     */
    function getTargetRows() {
        if ($anchor.hasClass('input-center-multiple')) {
            return $('#tbody_layers .layer_selectors:checked:visible').closest('tr').filter(function() {
                return $(this).find('.input-center').length > 0;
            });
        }

        return $anchor.closest('tr');
    }

    /**
     * Returns the values of a layer in the current view, after applying the normalization that is
     * selected for it, in the same way the drawer does. Values that are missing, or that can't be
     * normalized (e.g., negative values with sqrt), are left out.
     */
    function getNormalizedLayerValues($tr) {
        var layer_id = $tr.find('.input-center').attr('id').replace('center', '');
        var normalization = $tr.find('.normalization').val();
        var values = [];

        for (var i = 1; i < layerdata.length; i++) {
            if (!isNumber(layerdata[i][layer_id]))
                continue;

            var value = parseFloat(layerdata[i][layer_id]);

            if (normalization == 'sqrt') {
                value = Math.sqrt(value);
            } else if (normalization == 'log') {
                value = log10(value + 1);
            }

            if (isFinite(value))
                values.push(value);
        }

        return values;
    }

    function getMean(values) {
        return values.reduce(function(a, b) { return a + b; }, 0) / values.length;
    }

    function getMedian(values) {
        var sorted = values.slice().sort(function(a, b) { return a - b; });
        var middle = Math.floor(sorted.length / 2);

        return (sorted.length % 2) ? sorted[middle] : (sorted[middle - 1] + sorted[middle]) / 2;
    }

    /**
     * Sets the center value of a layer through its own input box, so everything that happens
     * when a user types a center value happens here, too
     */
    function setCenterValue($tr, value) {
        $tr.find('.input-center').val(value).trigger('change');
    }

    function applyStatistic(statistic_name, statistic_function) {
        var $rows = getTargetRows();
        var layers_without_values = [];

        $rows.each(function() {
            var values = getNormalizedLayerValues($(this));

            if (values.length === 0) {
                layers_without_values.push($(this).find('td.titles').attr('title'));
                return;
            }

            // six significant digits are more than enough to center a figure, and they are much easier to read
            setCenterValue($(this), Number(statistic_function(values).toPrecision(6)));
        });

        if (layers_without_values.length) {
            toastr.warning(`Anvi'o could not compute the ${statistic_name} for ${layers_without_values.length == 1 ? 'this layer' : 'these layers'}, ` +
                           `since ${layers_without_values.length == 1 ? 'it has' : 'they have'} no numerical values in the current view: ` +
                           `${layers_without_values.map(name => `'${name}'`).join(', ')}.`);
        }

        // the multiple layers box can't show one value for many layers with their own means or medians
        if ($anchor.hasClass('input-center-multiple')) {
            $anchor.val('');
        }

        close();
    }

    function applyManualValue() {
        var $input = $popup.find('.center-popup-value');
        var value = $input.val().trim();

        if (!isValidNumericalInputValue(value, true)) {
            $input.addClass('is-invalid');
            return;
        }

        getTargetRows().each(function() {
            setCenterValue($(this), value);
        });

        if ($anchor.hasClass('input-center-multiple')) {
            $anchor.val(value);
        }

        close();
    }

    function clear() {
        getTargetRows().each(function() {
            setCenterValue($(this), '');
        });

        if ($anchor.hasClass('input-center-multiple')) {
            $anchor.val('');
        }

        close();
    }

    function close() {
        if ($popup) {
            $popup.remove();
            $popup = null;
            $anchor = null;
        }
    }

    function open(anchor) {
        // clicking on the same box again should not rebuild the popup
        if ($popup && $anchor && $anchor[0] === anchor)
            return;

        close();
        $anchor = $(anchor);

        var $rows = getTargetRows();
        var is_multiple = $anchor.hasClass('input-center-multiple');
        var title, note;

        if (is_multiple) {
            title = `Center for ${$rows.length} selected layer${$rows.length == 1 ? '' : 's'}`;
            note = `Mean and median are computed for each selected layer separately. A value entered above applies to all of them.`;
        } else {
            title = `Center for '${$rows.find('td.titles').attr('title')}'`;
            note = `Mean and median are computed from the values of this layer in the current view, after its normalization (if any) is applied.`;
        }

        $popup = $(`<div class="center-popup card shadow-sm">
                        <div class="card-body p-2">
                            <div class="small font-weight-bold mb-2 text-truncate"></div>
                            <div class="input-group input-group-sm mb-2">
                                <input type="text" class="form-control center-popup-value" placeholder="enter a value">
                                <div class="input-group-append">
                                    <button type="button" class="btn btn-outline-secondary center-popup-apply">Set</button>
                                </div>
                            </div>
                            <div class="btn-group btn-group-sm d-flex mb-2" role="group">
                                <button type="button" class="btn btn-outline-primary w-100 center-popup-mean">Set Mean</button>
                                <button type="button" class="btn btn-outline-primary w-100 center-popup-median">Set Median</button>
                            </div>
                            <button type="button" class="btn btn-sm btn-outline-danger btn-block center-popup-clear">Clear</button>
                            <div class="small text-muted mt-2 center-popup-note"></div>
                        </div>
                    </div>`);

        // text is set separately, so layer names can't break the markup
        $popup.find('.font-weight-bold').text(title).attr('title', title);
        $popup.find('.center-popup-note').text(note);
        $popup.find('.center-popup-value').val(is_multiple ? '' : $anchor.val());

        if ($rows.length === 0) {
            $popup.find('button').prop('disabled', true);
            $popup.find('.center-popup-note').text('Please select one or more numerical layers first.');
        }

        $popup.css({'position': 'absolute', 'z-index': 2000, 'width': '240px'});
        $('body').append($popup);

        // place the popup right below the box that was clicked, without going off the screen
        var offset = $anchor.offset();
        var left = Math.min(offset.left, $(window).scrollLeft() + $(window).width() - $popup.outerWidth() - 10);
        $popup.css({'top': offset.top + $anchor.outerHeight() + 4, 'left': Math.max(left, 10)});

        $popup.find('.center-popup-mean').on('click', function() { applyStatistic('mean', getMean); });
        $popup.find('.center-popup-median').on('click', function() { applyStatistic('median', getMedian); });
        $popup.find('.center-popup-apply').on('click', applyManualValue);
        $popup.find('.center-popup-clear').on('click', clear);
        $popup.find('.center-popup-value').on('keydown', function(e) {
            $(this).removeClass('is-invalid');

            if (e.key === 'Enter') {
                e.preventDefault();
                applyManualValue();
            } else if (e.key === 'Escape') {
                close();
            }
        });

        $popup.find('.center-popup-value').trigger('focus').trigger('select');
    }

    $(document).ready(function() {
        $(document).on('focus click', '#tbody_layers .input-center, .input-center-multiple', function() {
            open(this);
        });

        // close the popup when the user clicks anywhere else, or presses escape
        $(document).on('mousedown', function(e) {
            if ($popup && !$(e.target).closest('.center-popup').length && !$(e.target).is($anchor)) {
                close();
            }
        });

        $(document).on('keydown', function(e) {
            if (e.key === 'Escape') {
                close();
            }
        });
    });

    return {'open': open, 'close': close};
})();
