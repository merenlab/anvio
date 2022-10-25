/**
 * Zooming in an out.
 *
 *  Authors: Özcan Esen <ozcanesen@gmail.com>
 *
 *  Copyright 2015-2021, The anvi'o project (http://anvio.org)
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


var mouse_event_origin_x = 0;
var mouse_event_origin_y = 0;

var dragging = false;
var zoomBox = {};
var drawing_zoom = false;

function initialize_area_zoom() {
    var viewport = document.getElementById('svg');

    viewport.addEventListener('mousedown',
        function(event) {
            dragging = false;
            document.activeElement.blur();

            mouse_event_origin_x = event.clientX;
            mouse_event_origin_y = event.clientY;

            if (event.shiftKey)
            {
                drawing_zoom = true;

                zoomBox['start_x'] = event.clientX;
                zoomBox['start_y'] = event.clientY;

                $('#divzoom').css({"top": 0, "left": 0, "width": 0, "height": 0 });
                $('#divzoom').show();
            }
        });

    viewport.addEventListener('mousemove',
        function(event) {
            if (Math.abs(mouse_event_origin_x - event.clientX) + Math.abs(mouse_event_origin_y - event.clientY) > 2)
            {
                dragging = true;
            }

            if (event.shiftKey && drawing_zoom)
            {
                var _top = zoomBox['start_y'] > event.clientY ? event.clientY : zoomBox['start_y'];
                var _left = zoomBox['start_x'] > event.clientX ? event.clientX : zoomBox['start_x'];
                var _height = Math.abs(zoomBox['start_y'] - event.clientY);
                var _width = Math.abs(zoomBox['start_x'] - event.clientX);

                var divzoom = document.getElementById('divzoom');

                divzoom.style.top = _top + "px";
                divzoom.style.left = _left + "px";
                divzoom.style.width = _width + "px";
                divzoom.style.height = _height + "px";

                // when you drawing rectangle, if you drag over text on the screen browser selects that text
                // with this hack you can continue drawing.
                clearTextSelection(); // in utils.js

            }
        });

    viewport.addEventListener('mouseup',
        function() {
            if (drawing_zoom)
            {
                var zoom_rect = document.getElementById('divzoom').getBoundingClientRect();

                if (zoom_rect.width > 2 && zoom_rect.height > 2)
                {
                    var _dx = (parseInt("0" + $('#svg').position().left) + (VIEWER_WIDTH / 2)) - (zoom_rect.left + zoom_rect.width / 2);
                    var _dy = (parseInt("0" + $('#svg').position().top)  + (VIEWER_HEIGHT / 2)) - (zoom_rect.top + zoom_rect.height / 2);
                    pan(_dx,_dy);
                    zoom(Math.min(VIEWER_WIDTH / zoom_rect.width, VIEWER_HEIGHT / zoom_rect.height));
                }
            }

            clearTextSelection();
            drawing_zoom=false;
            zoomBox = {};
            $('#divzoom').hide();

        });
}
