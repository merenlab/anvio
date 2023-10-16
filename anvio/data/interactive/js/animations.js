/**
 * Functions for panel animations.
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

var ANIMATIONS_ENABLED = true;
var SLIDE_INTERVAL = 4;
var SLIDE_STEP_SIZE = 15;

var is_left_panel_sliding = false;

function toggleLeftPanel() {
    if (is_left_panel_sliding)
        return;

    is_left_panel_sliding = true;

    if ($('#panel-left').is(':visible')) {
        var animation_frame = function(){
            if (ANIMATIONS_ENABLED && $('#panel-left')[0].getBoundingClientRect().right > 0) {
                $('#panel-left').css('left', parseInt($('#panel-left').css('left')) - SLIDE_STEP_SIZE);
                $('#toggle-panel-left').css('left', $('#sidebar')[0].getBoundingClientRect().right + 'px');
                setTimeout(animation_frame, SLIDE_INTERVAL);
            }
            else {
                $('#panel-left').hide();
                $('#toggle-panel-left').css('left', '');
                $('#toggle-panel-left').removeClass('toggle-panel-left-pos');
                $('#toggle-panel-left-inner').html('&#9658;');
                is_left_panel_sliding = false;
            }
        };
        animation_frame();
    } else {
        $('#panel-left').show();
        var animation_frame = function(){
            if (ANIMATIONS_ENABLED && $('#panel-left')[0].getBoundingClientRect().left < 0) {
                $('#panel-left').css('left', parseInt($('#panel-left').css('left')) + SLIDE_STEP_SIZE);
                $('#toggle-panel-left').css('left', $('#sidebar')[0].getBoundingClientRect().right + 'px');
                setTimeout(animation_frame, SLIDE_INTERVAL);
            }
            else {
                $('#panel-left').css('left', '0px');
                $('#toggle-panel-left').css('left', '');
                $('#toggle-panel-left').addClass('toggle-panel-left-pos');
                $('#toggle-panel-left-inner').html("&#9664;");
                is_left_panel_sliding = false;
            }
        };
        animation_frame();
    }
}

function switchNavigationTabs() {
    if ($('#panel-left').is(':visible')) {
        $('a').each(function(){
            if ($(this).prop('href') == window.location.href) {
                $(this).addClass('active'); 
                $(this).parents('li').addClass('active');
            } else 
            { 
              $(this).removeClass('active');
              $(this).parents('li').removeClass('active');
            }
        });
    }
}

// Outdated with new UI but lets keep it for future development
function toggleRightPanel(name) {
    ['#mouse_hover_panel', '#description-panel', '#news-panel'].forEach(function(right_panel) {
        if (right_panel == name)
            return;

        $(right_panel).hide();
    });

    $(name).toggle();

    if ($('#mouse_hover_panel').is(':visible')) {
        $('#toggle-panel-right').addClass('toggle-panel-right-pos');
        $('#toggle-panel-right-inner').html('&#9658;');
    } else {
        $('#toggle-panel-right').removeClass('toggle-panel-right-pos');
        $('#toggle-panel-right-inner').html('&#9664;');
    }

    if ($('#description-panel').is(':visible')) {
        $('#toggle-panel-right-2').addClass('toggle-panel-right-pos-2');
        $('#toggle-panel-right-2-inner').html('&#9658;');
    } else {
        $('#toggle-panel-right-2').removeClass('toggle-panel-right-pos-2');
        $('#toggle-panel-right-2-inner').html('&#9664;');
    }

    if ($('#news-panel').is(':visible')) {
        checkNews();
        $('#toggle-panel-right-3').addClass('toggle-panel-right-pos-3');
        $('#toggle-panel-right-3-inner').html('&#9658;');
    } else {
        $('#toggle-panel-right-3').removeClass('toggle-panel-right-pos-3');
        $('#toggle-panel-right-3-inner').html('&#9664;');
    }
}

