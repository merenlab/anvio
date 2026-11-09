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
var is_panel_open = true;

function toggleLeftPanel() {
    if (is_left_panel_sliding)
        return;

    is_left_panel_sliding = true;

    const panel = document.getElementById('panel-left');
    const svgContainer = document.getElementById('svgbox');
    if (svgContainer) svgContainer.style.pointerEvents = 'none';

    let cleaned_up = false;
    let fallback_timer = null;

    const cleanup = () => {
        if (cleaned_up) return;   // run exactly once (transitionend vs. fallback race)
        cleaned_up = true;
        if (fallback_timer !== null) clearTimeout(fallback_timer);
        panel.removeEventListener('transitionend', onDone);
        if (svgContainer) svgContainer.style.pointerEvents = '';
        is_left_panel_sliding = false;
    };

    const onDone = (ev) => {
        if (ev && ev.target !== panel) return;   // ignore transitions bubbling from descendants
        cleanup();
    };

    panel.addEventListener('transitionend', onDone);
    // Fallback: if no CSS 'transition' is defined on #panel-left, transitionend never
    // fires. Slightly longer than the 0.3s transition so the event wins when it exists.
    fallback_timer = setTimeout(cleanup, 350);

    if (is_panel_open) {
        is_panel_open = false;
        panel.classList.add('panel-closed');
        $('#toggle-panel-left').css('left', '');
        $('#toggle-panel-left').removeClass('toggle-panel-left-pos');
        $('#toggle-panel-left-inner').html('&#9658;');
    } else {
        is_panel_open = true;
        panel.classList.remove('panel-closed');
        $('#toggle-panel-left').css('left', '');
        $('#toggle-panel-left').addClass('toggle-panel-left-pos');
        $('#toggle-panel-left-inner').html('&#9664;');
    }
}

function switchNavigationTabs(tab_number) {
    if (is_panel_open) {
        $('a').each(function(){
            if (tab_number){
                $(tab_number).addClass('active');
                $(tab_number).tab('show');
            }
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

